# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import os
import shutil
from abc import ABC, abstractmethod

from vasp_manager.job_manager import JobManager


class BaseCalculationManager(ABC):
    """
    Run vasp job workflow for a single material
    """

    def __init__(
        self,
        base_path,
        to_rerun,
        to_submit,
        ignore_personal_errors=True,
        from_scratch=False,  ## DANGEROUS, WILL DELETE PREVIOUS CALCULATION
    ):
        self.base_path = base_path
        self.to_rerun = to_rerun
        self.to_submit = to_submit
        self.job_manager = JobManager(
            calc_path=self.calc_path, ignore_personal_errors=ignore_personal_errors
        )

        if from_scratch:
            self._from_scratch()

    @property
    @abstractmethod
    def mode(self):
        pass

    @abstractmethod
    def setup_calc(self):
        pass

    @abstractmethod
    def check_calc(self):
        pass

    @property
    @abstractmethod
    def is_done(self):
        pass

    @property
    def calc_path(self):
        return os.path.join(self.base_path, self.mode)

    @property
    def material_name(self):
        return self.base_path.split("/")[-1]

    @property
    def job_exists(self):
        return self.job_manager.job_exists

    @property
    def job_complete(self):
        return self.job_manager.job_complete

    def submit_job(self):
        return self.job_manager.submit_job()

    def _from_scratch(self):
        shutil.rmtree(self.calc_path)
