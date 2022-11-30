# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import glob
import logging
import os
from functools import cached_property

from vasp_manager.calculation_manager.base import BaseCalculationManager
from vasp_manager.utils import ptail
from vasp_manager.vasp_input_creator import VaspInputCreator

logger = logging.getLogger(__name__)


class RlxCoarseCalculationManager(BaseCalculationManager):
    """
    Runs coarse relaxation job workflow for a single material
    """

    def __init__(
        self,
        material_path,
        to_rerun,
        to_submit,
        primitive=True,
        ignore_personal_errors=True,
        from_scratch=False,
        tail=5,
    ):
        """
        For material_path, to_rerun, to_submit, ignore_personal_errors, and from_scratch,
        see BaseCalculationManager

        Args:
            tail (int): number of last lines to log in debugging if job failed
        """
        self.tail = tail
        super().__init__(
            material_path=material_path,
            to_rerun=to_rerun,
            to_submit=to_submit,
            primitive=primitive,
            ignore_personal_errors=ignore_personal_errors,
            from_scratch=from_scratch,
        )
        self._is_done = None
        self._results = None

    @cached_property
    def mode(self):
        return "rlx-coarse"

    @cached_property
    def poscar_source_path(self):
        poscar_source_path = os.path.join(self.material_path, "POSCAR")
        return poscar_source_path

    def setup_calc(self):
        """
        Sets up a coarse relaxation
        """
        vasp_input_creator = VaspInputCreator(
            self.calc_path,
            mode=self.mode,
            poscar_source_path=self.poscar_source_path,
            primitive=self.primitive,
            name=self.material_name,
        )
        if self.to_rerun:
            archive_made = vasp_input_creator.make_archive_and_repopulate()
            if not archive_made:
                # set rerun to False to not make an achive and instead
                # continue to make the input files
                self.to_rerun = False
                self.setup_calc()
                return
        else:
            vasp_input_creator.create()

        if self.to_submit:
            job_submitted = self.submit_job()
            # job status returns True if sucessfully submitted, else False
            if not job_submitted:
                # if job didn't submit, try rerunning setup
                self.to_rerun = False
                self.setup_calc()

    def check_calc(self):
        """
        Checks if calculation has finished and reached required accuracy

        Returns:
            relaxation_successful (bool): if True, relaxation completed successfully
        """
        stdout_path = os.path.join(self.calc_path, "stdout.txt")
        if not os.path.exists(stdout_path):
            logger.info(f"{self.mode.upper()} not started")
            return False

        if not self.job_complete:
            logger.info(f"{self.mode.upper()} not finished")
            return False

        tail_output = ptail(stdout_path, n_tail=self.tail, as_string=True)
        if "reached required accuracy" not in tail_output:
            archive_dirs = glob.glob(os.path.join(self.calc_path, "archive*"))
            if len(archive_dirs) >= 3:
                logger.warning("Many archives exist, suggest force based relaxation")
                if self.to_rerun:
                    self.setup_calc()
                return True

            logger.warning(f"{self.mode.upper()} FAILED")
            logger.debug(tail_output)
            if self.to_rerun:
                logger.info(f"Rerunning {self.calc_path}")
                self.setup_calc()
            return False

        logger.info(f"{self.mode.upper()} Calculation: reached required accuracy")
        logger.debug(tail_output)
        return True

    @property
    def is_done(self):
        if self._is_done is None:
            self._is_done = self.check_calc()
        return self._is_done

    @property
    def results(self):
        if not self.is_done:
            self._results = "not finished"
        else:
            self._results = "done"
        return self._results
