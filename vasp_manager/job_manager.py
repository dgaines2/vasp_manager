# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import json
import logging
import os
import pkgutil
import subprocess

from vasp_manager.utils import change_directory

logger = logging.getLogger(__name__)


class JobManager:
    """
    A JobManager -- handles job submission and status
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

    @property
    def computing_config_dict(self):
        fpath = "calculations/computing_config.json"
        if os.path.exists(fpath):
            with open(fpath) as fr:
                computing_config = json.load(fr)
        else:
            raise Exception(f"No {fpath} found in current path")
        return computing_config

    @property
    def computer(self):
        return self.computing_config_dict["computer"]

    @property
    def user_id(self):
        return self.computing_config_dict[self.computer]["user_id"]

    @property
    def mode(self):
        return self.calc_path.split("/")[-1]

    @property
    def job_exists(self):
        jobid_path = os.path.join(self.calc_path, "jobid")
        if os.path.exists(jobid_path):
            with open(jobid_path) as fr:
                jobid = fr.read().strip()
            self.jobid = jobid
            return True
        else:
            return False

    @property
    def jobid(self):
        if self._jobid is None:
            raise Exception("jobid has not been set")
        else:
            return self._jobid

    @jobid.setter
    def jobid(self, job_value):
        # some criteria to make sure jobid actually looks reasonable?
        # but I think SLURM will throw and error if sbatch fails
        try:
            jobid_int = int(job_value)
        except Exception as e:
            raise Exception(f"{e}")
        self._jobid = jobid_int
        return self._jobid

    def submit_job(self):
        if self.job_exists:
            logger.info(f"{self.mode.upper()} Job already exists")
            return True

        if "personal" in self.computer:
            error_msg = f"Cannot submit {self.mode.upper()} job for on personal computer"
            error_msg += "\n\tIgnoring job submission..."
            logger.debug(error_msg)
            return True

        vaspq_path = os.path.join(self.calc_path, "vasp.q")
        if not os.path.exists(vaspq_path):
            logger.info(f"No vasp.q file in {self.calc_path}")
            # return False here instead of catching an exception
            # This enables job resubmission by letting the calling function
            # know that the calculation needs to be restarted
            return False

        submission_call = "sbatch vasp.q | awk '{ print $4 }' | tee jobid"
        with change_directory(self.calc_path):
            jobid = subprocess.check_output(submission_call, shell=True).decode("utf-8")
        logger.info(f"Submitted job {jobid}")
        self.jobid = jobid
        return True

    def _check_job_complete(self):
        """Return True if job done"""
        if self.computer == "personal":
            error_msg = "Cannot check job on personal computer"
            error_msg += "\n\tIgnoring job status check..."
            logger.debug(error_msg)
            # This enables job resubmission by letting the calling function
            # continue anyways
            return True
        else:
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

    @property
    def job_complete(self):
        return self._check_job_complete()
