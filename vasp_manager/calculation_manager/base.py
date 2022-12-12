# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import os
import shutil
from abc import ABC, abstractmethod
from functools import cached_property

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
