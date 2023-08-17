# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import json
import logging
import subprocess
from functools import cached_property

from vasp_manager.utils import change_directory

logger = logging.getLogger(__name__)


class JobManager:
    """
    Handles job submission and status monitoring
    """

    def __init__(self, calc_path, ignore_personal_errors=True):
        """
        Args:
            calc_path (str): base directory of job
            ignore_personal_errors (bool): if True, ignore job submission errors
                if on personal computer
        """
        self.calc_path = calc_path
        self.ignore_personal_errors = ignore_personal_errors
        self._jobid = None
        self._job_complete = None

    @property
    def computing_config_dict(self):
        all_calcs_dir = self.calc_path.parent.parent
        fname = "computing_config.json"
        fpath = all_calcs_dir / "computing_config.json"
        if fpath.exists():
            with open(fpath) as fr:
                computing_config = json.load(fr)
        else:
            raise Exception(f"No {fname} found in path {all_calcs_dir.absolute()}")
        return computing_config

    @cached_property
    def computer(self):
        return self.computing_config_dict["computer"]

    @cached_property
    def user_id(self):
        return self.computing_config_dict[self.computer]["user_id"]

    @cached_property
    def mode(self):
        return self.calc_path.name

    @property
    def job_exists(self):
        jobid_path = self.calc_path / "jobid"
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
            raise Exception(f"jobid has not been set in {self.calc_path}")
        return self._jobid

    @jobid.setter
    def jobid(self, job_value):
        # some criteria to make sure jobid actually looks reasonable?
        # but I think SLURM will throw and error if sbatch fails
        try:
            jobid_int = int(job_value)
        except Exception as e:
            raise Exception(f"Tried to set jobid={job_value}\n{e}")
        self._jobid = jobid_int

    def submit_job(self):
        """
        Submits job, making sure to not make duplicate jobs
        """
        if self.job_exists:
            logger.info(f"{self.mode.upper()} Job already exists")
            return True

        if "personal" in self.computer:
            error_msg = f"Cannot submit {self.mode.upper()} job for on personal computer"
            error_msg += "\n\tIgnoring job submission..."
            logger.debug(error_msg)
            return True

        vaspq_path = self.calc_path / "vasp.q"
        if not vaspq_path.exists():
            logger.info(f"No vasp.q file in {self.calc_path}")
            # return False here instead of catching an exception
            # This enables job resubmission by letting the calling function
            # know that the calculation needs to be restarted
            return False

        submission_call = "sbatch vasp.q | awk '{ print $4 }'"
        with change_directory(self.calc_path):
            jobid = (
                subprocess.check_output(submission_call, shell=True)
                .decode("utf-8")
                .strip()
            )
            self.jobid = jobid
            with open("jobid", "w+") as fw:
                fw.write(jobid)
        logger.info(f"Submitted job {jobid}")
        return True

    @property
    def job_complete(self):
        if self._job_complete is None and self.job_exists:
            self._job_complete = self._check_job_complete()
        return self._job_complete

    def _check_job_complete(self):
        """Returns True if job done"""
        if self.computer == "personal":
            error_msg = "Cannot check job on personal computer"
            error_msg += "\n\tIgnoring job status check..."
            logger.debug(error_msg)
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
