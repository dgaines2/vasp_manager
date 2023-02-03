# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import os
import shutil
from abc import ABC, abstractmethod
from functools import cached_property

import numpy as np

from vasp_manager.job_manager import JobManager


class BaseCalculationManager(ABC):
    """
    Runs vasp job workflow for a single material
    """

    def __init__(
        self,
        material_path,
        to_rerun,
        to_submit,
        primitive=True,
        ignore_personal_errors=True,
        from_scratch=False,  # DANGEROUS, WILL DELETE PREVIOUS CALCULATION
    ):
        """
        Args:
            material_path (str): path for a single material
                ex. calculations/AlAs
            to_rerun (bool): if True, rerun failed calculations
            to_submit (bool): if True, submit calculations to job manager
            primitive (bool): if True, find primitive cell, else find conventional cell
            ignore_personal_errors (bool): if True, ignore job submission errors
                if on personal computer
            from_scratch (bool): if True, remove the calculation's folder and
                restart
                note: DANGEROUS
        """
        self.material_path = material_path
        self.to_rerun = to_rerun
        self.to_submit = to_submit
        self.primitive = primitive
        self.job_manager = JobManager(
            calc_path=self.calc_path, ignore_personal_errors=ignore_personal_errors
        )

        if from_scratch:
            self._from_scratch()

    @property
    @abstractmethod
    def mode(self):
        pass

    @property
    @abstractmethod
    def is_done(self):
        pass

    @property
    @abstractmethod
    def vasp_input_creator(self):
        pass

    @cached_property
    def calc_path(self):
        return os.path.join(self.material_path, self.mode)

    @cached_property
    def material_name(self):
        return os.path.basename(self.material_path)

    @abstractmethod
    def setup_calc(self):
        pass

    @abstractmethod
    def check_calc(self):
        pass

    @property
    def job_exists(self):
        return self.job_manager.job_exists

    @property
    def job_complete(self):
        return self.job_manager.job_complete

    @property
    def stopped(self):
        stopped_path = os.path.join(self.material_path, "STOP")
        stopped = os.path.exists(stopped_path)
        return stopped

    def submit_job(self):
        return self.job_manager.submit_job()

    def _cancel_previous_job(self):
        jobid_path = os.path.join(self.calc_path, "jobid")
        if os.path.exists(jobid_path):
            with open(jobid_path) as fr:
                jobid = fr.read().strip()
            cancel_job_call = f"scancel {jobid}"
            os.system(cancel_job_call)

    def _from_scratch(self):
        self._cancel_previous_job()
        shutil.rmtree(self.calc_path)

    def _check_vasp_errors(self, stdout_path=None, stderr_path=None):
        """
        Find VASP errors in stdout and stderr
        """
        if stdout_path is None:
            stdout_path = os.path.join(self.calc_path, "stdout.txt")
        if stderr_path is None:
            stderr_path = os.path.join(self.calc_path, "stderr.txt")
        errors = set()

        with open(stdout_path) as fr:
            stdout = fr.readlines()
        for line in stdout:
            if "Sub-Space-Matrix" in line:
                errors.add("Sub-Space-Matrix")

        with open(stderr_path) as fr:
            stderr = fr.readlines()
        for line in stderr:
            if "oom-kill" in line:
                errors.add("OOM")

        return errors

    def _address_vasp_errors(self, errors):
        """
        Args:
            errors (set): set of errors found in stdout or stderr
        Returns:
            all_errors_addressed (bool): if True, all errors could
                be fixed automatically. If False, some errors
                could not be handled
        """
        vic = self.vasp_input_creator

        errors_addressed = {e: False for e in errors}
        for error in errors:
            match error:
                case "Sub-Space-Matrix":
                    vic.calc_config_dict[self.mode]["algo"] = "Fast"
                    errors_addressed[error] = True
                case "OOM":
                    if vic.computer == "quest":
                        ncore_per_node_for_memory = 8
                    else:
                        ncore_per_node_for_memory = 32
                    vic.ncore_per_node_for_memory = ncore_per_node_for_memory
                    errors_addressed[error] = True
                case _:
                    errors_addressed[error] = False
        all_errors_addressed = np.all([v for v in errors_addressed.values()])
        return all_errors_addressed
