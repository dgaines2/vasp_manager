# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

from __future__ import annotations

import json
import logging
import subprocess
from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING

from vasp_manager.utils import LoggerAdapter, change_directory

if TYPE_CHECKING:
    from vasp_manager.types import SourceDirectory, WorkingDirectory

logger = logging.getLogger(__name__)


class JobManager:
    """
    Handles job submission and status monitoring
    """

    def __init__(
        self,
        calc_dir: WorkingDirectory,
        exe_name: str = "vasp.q",
        jobid_name: str = "jobid",
        config_dir: SourceDirectory | None = None,
        manager_name: str = None,
        ignore_personal_errors: bool = True,
    ) -> None:
        """
        Args:
            calc_dir: base directory of job
            exe_name: filename for slurm job submission
            jobid_name: filename for storing slurm jobid
            config_dir: path to directory containing configuration files. If
                None, use the parent directory of calc_dir
            manager_name: name for logging purposes
            ignore_personal_errors: if True, ignore job submission errors
                if on personal computer
        """
        self.calc_dir = Path(calc_dir)
        self.config_dir = Path(config_dir) if config_dir else self.calc_dir.parents[1]
        self.ignore_personal_errors = ignore_personal_errors
        self.exe_name = exe_name
        self.jobid_name = jobid_name
        self.manager_name = manager_name if manager_name else str(self.calc_dir)
        self._jobid: int
        self._job_complete: bool
        self.logger = LoggerAdapter(logging.getLogger(__name__), self.manager_name)

    @property
    def computing_config_dict(self) -> dict:
        fname = "computing_config.json"
        fpath = self.config_dir / fname
        if fpath.exists():
            with open(fpath) as fr:
                computing_config = json.load(fr)
        else:
            raise Exception(f"No {fname} found in {self.config_dir.absolute()}")
        return computing_config

    @cached_property
    def computer(self) -> str:
        return self.computing_config_dict["computer"]

    @cached_property
    def user_id(self) -> str:
        return self.computing_config_dict[self.computer]["user_id"]

    @cached_property
    def mode(self) -> str:
        return self.calc_dir.name

    @property
    def job_exists(self) -> bool:
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
    def jobid(self) -> int:
        if not self.job_exists:
            raise Exception(f"jobid has not been set in {self.calc_dir}")
        return self._jobid

    @jobid.setter
    def jobid(self, job_value: str | int) -> None:
        try:
            jobid_int = int(job_value)
        except Exception as e:
            raise RuntimeError(f"Tried to set jobid={job_value}\n{e}")
        self._jobid = jobid_int

    def submit_job(self) -> bool:
        """
        Submits job, making sure to not make duplicate jobs

        Returns:
            job_submitted_successfully
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
    def job_complete(self) -> bool:
        if getattr(self, "_job_complete", None) is None and self.job_exists:
            self._job_complete = self._check_job_complete()
        else:
            self._job_complete = False
        return self._job_complete

    def _check_job_complete(self) -> bool:
        """Returns True if job complete"""
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
            line_as_str = line.strip().split()
            if str(self.jobid) in line_as_str:
                return False
        return True
