# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import json
import logging
import subprocess
from functools import cached_property
from pathlib import Path

from vasp_manager.utils import LoggerAdapter, change_directory

logger = logging.getLogger(__name__)


class JobManager:
    """
    Handles job submission and status monitoring
    """

    def __init__(
        self,
        calc_dir,
        exe_name="vasp.q",
        jobid_name="jobid",
        config_dir=None,
        manager_name=None,
        ignore_personal_errors=True,
    ):
        """
        Args:
            calc_dir (str | Path): base directory of job
            exe_name (str): filename for slurm job submission
            jobid_name (str): filename for storing slurm jobid
            config_dir (str | Path): path to directory containing configuration files
            manager_name (str): name for logging purposes
            ignore_personal_errors (bool): if True, ignore job submission errors
                if on personal computer
        """
        self.calc_dir = Path(calc_dir)
        self.config_dir = Path(config_dir) if config_dir else self.calc_dir.parents[1]
        self.ignore_personal_errors = ignore_personal_errors
        self.exe_name = exe_name
        self.jobid_name = jobid_name
        self.manager_name = manager_name if manager_name else str(self.calc_dir)
        self._jobid = None
        self._job_complete = None
        self.logger = LoggerAdapter(logging.getLogger(__name__), self.manager_name)

    @property
    def computing_config_dict(self):
        fname = "computing_config.json"
        fpath = self.config_dir / fname
        if fpath.exists():
            with open(fpath) as fr:
                computing_config = json.load(fr)
        else:
            raise Exception(f"No {fname} found in {self.config_dir.absolute()}")
        return computing_config

    @cached_property
    def computer(self):
        return self.computing_config_dict["computer"]

    @cached_property
    def user_id(self):
        return self.computing_config_dict[self.computer]["user_id"]

    @cached_property
    def mode(self):
        return self.calc_dir.name

    @property
    def job_exists(self):
        jobid_path = self.calc_dir / self.jobid_name
        if not jobid_path.exists():
            return False
        if jobid_path.stat().st_size == 0:
            return False

        with open(jobid_path) as fr:
            jobid = fr.read().strip()
        self.jobid = jobid
        return True

    @property
    def jobid(self):
        if not self.job_exists:
            raise Exception(f"jobid has not been set in {self.calc_dir}")
        return self._jobid

    @jobid.setter
    def jobid(self, job_value):
        # some criteria to make sure jobid actually looks reasonable?
        # but I think SLURM will throw and error if sbatch fails
        try:
            jobid_int = int(job_value)
        except Exception as e:
            raise RuntimeError(f"Tried to set jobid={job_value}\n{e}")
        self._jobid = jobid_int

    def submit_job(self):
        """
        Submits job, making sure to not make duplicate jobs
        """
        if self.job_exists:
            self.logger.info(f"{self.mode.upper()} job already exists")
            return True

        if "personal" in self.computer:
            msg = (
                f"Cannot submit {self.mode.upper()} job for on personal computer\n"
                "\tIgnoring job submission..."
            )
            self.logger.debug(msg)
            return True

        qpath = self.calc_dir / self.exe_name
        if not qpath.exists():
            self.logger.info(f"No {self.exe_name} file in {self.calc_dir}")
            # return False here instead of catching an exception
            # This enables job resubmission by letting the calling function
            # know that the calculation needs to be restarted
            return False

        submission_call = f"sbatch {self.exe_name}"
        with change_directory(self.calc_dir):
            submission_call_output = (
                subprocess.check_output(submission_call, shell=True)
                .decode("utf-8")
                .strip()
            )
            jobid = submission_call_output.split(" ")[3]
            self.jobid = jobid
            with open(self.jobid_name, "w+") as fw:
                fw.write(f"{jobid}\n")
        self.logger.info(f"Submitted job {jobid}")
        return True

    @property
    def job_complete(self):
        if self._job_complete is None and self.job_exists:
            self._job_complete = self._check_job_complete()
        return self._job_complete

    def _check_job_complete(self):
        """Returns True if job done"""
        if self.computer == "personal":
            msg = "Cannot check job on personal computer\n\tIgnoring job status check..."
            self.logger.debug(msg)
            # This enables job resubmission by letting the calling function
            # continue anyways
            return True

        check_queue_call = f"squeue -u {self.user_id}"
        queue_call = (
            subprocess.check_output(check_queue_call, shell=True)
            .decode("utf-8")
            .splitlines()
        )
        for line in queue_call:
            line = line.strip().split()
            if str(self.jobid) in line:
                return False
        return True
