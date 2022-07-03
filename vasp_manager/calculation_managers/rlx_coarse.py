# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import glob
import logging
import os

from vasp_manager.calculation_managers.base import BaseCalculationManager
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
            ignore_personal_errors=ignore_personal_errors,
            from_scratch=from_scratch,
        )
        self._results = "not complete"

    @property
    def mode(self):
        return "rlx-coarse"

    @property
    def poscar_source_path(self):
        poscar_source_path = os.path.join(self.calc_path, "POSCAR")
        return poscar_source_path

    def setup_calc(self):
        """
        Sets up a coarse relaxation
        """
        vasp_input_creator = VaspInputCreator(
            self.calc_path,
            mode=self.mode,
            poscar_source_path=self.poscar_source_path,
            name=self.material_name,
        )
        if self.to_rerun:
            archive_made = vasp_input_creator.make_archive()
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
        if not self.job_complete:
            logger.info(f"{self.mode.upper()} job not finished")
            return False

        stdout_path = os.path.join(self.calc_path, "stdout.txt")
        if os.path.exists(stdout_path):
            if not self.job_complete:
                logger.info(f"{self.mode.upper()} not finished")
                return False

            tail_output = ptail(stdout_path, n_tail=self.tail, as_string=True)
            if "reached required accuracy" in tail_output:
                logger.info(
                    f"{self.mode.upper()} Calculation: reached required accuracy"
                )
                logger.debug(tail_output)
                self.results = "done"
                return True
            else:
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
        else:
            logger.info(f"{self.mode.upper()} not started")
            return False

    @property
    def is_done(self):
        return self.check_calc()

    @property
    def results(self):
        return self._results

    @results.setter
    def results(self, value):
        if not "done" in value:
            raise Exception
        self._results = value
        return self._results
