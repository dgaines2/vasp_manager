# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import logging
from functools import cached_property

from vasp_manager.calculation_manager.base import BaseCalculationManager
from vasp_manager.utils import pgrep, ptail
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
        max_reruns=3,
    ):
        """
        For material_path, to_rerun, to_submit, ignore_personal_errors, and from_scratch,
        see BaseCalculationManager

        Args:
            tail (int): number of last lines to log in debugging if job failed
        """
        self.tail = tail
        self.max_reruns = max_reruns
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
        return self.material_path / "POSCAR"

    @cached_property
    def vasp_input_creator(self):
        return VaspInputCreator(
            self.calc_path,
            mode=self.mode,
            poscar_source_path=self.poscar_source_path,
            primitive=self.primitive,
            name=self.material_name,
        )

    def setup_calc(
        self,
        increase_nodes_by_factor=1,
        increase_walltime_by_factor=1,
        make_archive=False,
    ):
        """
        Sets up a coarse relaxation
        """
        self.vasp_input_creator.increase_nodes_by_factor = increase_nodes_by_factor
        self.vasp_input_creator.increase_walltime_by_factor = increase_walltime_by_factor

        if make_archive:
            self.vasp_input_creator.make_archive_and_repopulate()
        else:
            self.vasp_input_creator.create()

        if self.to_submit:
            job_submitted = self.submit_job()
            # job status returns True if sucessfully submitted, else False
            if not job_submitted:
                # if job didn't submit, try rerunning setup
                self.setup_calc()

    def check_calc(self):
        """
        Checks if calculation has finished and reached required accuracy

        Returns:
            relaxation_successful (bool): if True, relaxation completed successfully
        """
        if not self.job_complete:
            logger.info(f"{self.mode.upper()} not finished")
            return False

        stdout_path = self.calc_path / "stdout.txt"
        if not stdout_path.exists():
            # calculation never actually ran
            # shouldn't get here unless function was called with submit=False
            # or job was manually cancelled
            logger.info(f"{self.mode.upper()} Calculation: No stdout.txt available")
            if self.to_rerun:
                self._cancel_previous_job()
                self.setup_calc()
            return False

        vasp_errors = self._check_vasp_errors()
        if len(vasp_errors) > 0:
            all_errors_addressed = self._address_vasp_errors(vasp_errors)
            if all_errors_addressed:
                if self.to_rerun:
                    logger.info(f"Rerunning {self.calc_path}")
                    self.setup_calc(make_archive=True)
            else:
                msg = (
                    f"{self.mode.upper()} Calculation: ",
                    "Couldn't address all VASP Errors\n",
                    "\tRefusing to continue...\n",
                    f"\tVasp Errors: {vasp_errors}\n",
                )
                logger.error(msg)
                self.stop()
            return False

        tail_output = ptail(stdout_path, n_tail=self.tail, as_string=True)
        grep_output = pgrep(
            stdout_path, "reached required accuracy", stop_after_first_match=True
        )
        if len(grep_output) == 0:
            archive_dirs = list(self.calc_path.glob("archive*"))
            if len(archive_dirs) >= self.max_reruns - 1:
                logger.warning(
                    "Many archives exist, continuing to force based relaxation..."
                )
                return True

            logger.warning(f"{self.mode.upper()} FAILED")
            logger.debug(tail_output)
            if self.to_rerun:
                logger.info(f"Rerunning {self.calc_path}")
                # increase nodes as its likely the calculation failed
                self.setup_calc(increase_walltime_by_factor=2, make_archive=True)
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
            if self.stopped:
                return "STOPPED"
            else:
                return "not finished"
        else:
            self._results = "done"
        return self._results
