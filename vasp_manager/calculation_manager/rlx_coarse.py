# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

from __future__ import annotations

import logging
from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING

from vasp_manager.calculation_manager.base import BaseCalculationManager
from vasp_manager.utils import LoggerAdapter, pgrep, ptail
from vasp_manager.vasp_input_creator import VaspInputCreator

if TYPE_CHECKING:
    from vasp_manager.types import CalculationType, WorkingDirectory

logger = logging.getLogger(__name__)


class RlxCoarseCalculationManager(BaseCalculationManager):
    """
    Runs coarse relaxation job workflow for a single material
    """

    def __init__(
        self,
        material_dir: WorkingDirectory,
        to_rerun: bool,
        to_submit: bool,
        primitive: bool = True,
        ignore_personal_errors: bool = True,
        from_scratch: bool = False,
        tail: int = 5,
        max_reruns: int = 3,
    ) -> None:
        """
        Args:
            material_dir: path to a directory for a single material
            to_rerun: if True, rerun failed calculations
            to_submit: if True, submit calculations to job manager
            primitive: if True, find primitive cell, else find conventional cell
            ignore_personal_errors: if True, ignore job submission errors
                if on personal computer
            from_scratch: if True, remove the calculation's directory and
                restart
            tail: number of last lines to log in debugging if job failed
            max_reruns: maximum number of times to rerun rlx-coarse before
                continuing to rlx
        """
        self.tail = tail
        self.max_reruns = max_reruns
        super().__init__(
            material_dir=material_dir,
            to_rerun=to_rerun,
            to_submit=to_submit,
            primitive=primitive,
            ignore_personal_errors=ignore_personal_errors,
            from_scratch=from_scratch,
        )
        self._is_done: bool
        self._results: str
        self.logger = LoggerAdapter(logging.getLogger(__name__), self.material_name)

    @cached_property
    def mode(self) -> CalculationType:
        return "rlx-coarse"

    @cached_property
    def poscar_source_path(self) -> Path:
        return self.material_dir / "POSCAR"

    @cached_property
    def vasp_input_creator(self) -> VaspInputCreator:
        return VaspInputCreator(
            self.calc_dir,
            mode=self.mode,
            poscar_source_path=self.poscar_source_path,
            primitive=self.primitive,
            name=self.material_name,
        )

    def setup_calc(
        self,
        increase_nodes_by_factor: int = 1,
        increase_walltime_by_factor: int = 1,
        make_archive: bool = False,
    ) -> None:
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

    def check_calc(self) -> bool:
        """
        Checks if calculation has finished and reached required accuracy

        Returns:
            relaxation_successful (bool): if True, relaxation completed successfully
        """
        if not self.job_complete:
            self.logger.info(f"{self.mode.upper()} not finished")
            return False

        stdout_path = self.calc_dir / "stdout.txt"
        if not stdout_path.exists():
            # calculation never actually ran
            # shouldn't get here unless function was called with submit=False
            # or job was manually cancelled
            self.logger.info(f"{self.mode.upper()} Calculation: No stdout.txt available")
            if self.to_rerun:
                self._cancel_previous_job()
                self.setup_calc()
            return False

        vasp_errors = self._check_vasp_errors()
        if len(vasp_errors) > 0:
            all_errors_addressed = self._address_vasp_errors(vasp_errors)
            if all_errors_addressed:
                if self.to_rerun:
                    self.logger.info(f"Rerunning {self.calc_dir}")
                    self.setup_calc(make_archive=True)
            else:
                msg = (
                    f"{self.mode.upper()} Calculation: "
                    "Couldn't address all VASP Errors\n"
                    f"\tVASP Errors: {vasp_errors}\n"
                    "\tRefusing to continue...\n"
                )
                self.logger.error(msg)
                self.stop()
            return False

        tail_output = ptail(stdout_path, n_tail=self.tail, as_string=True)
        grep_output = pgrep(
            stdout_path, "reached required accuracy", stop_after_first_match=True
        )
        if len(grep_output) == 0:
            archive_dirs = list(self.calc_dir.glob("archive*"))
            if len(archive_dirs) >= self.max_reruns - 1:
                self.logger.warning(
                    "Many archives exist, continuing to force based relaxation..."
                )
                return True

            self.logger.warning(f"{self.mode.upper()} FAILED")
            self.logger.debug(tail_output)
            if self.to_rerun:
                self.logger.info(f"Rerunning {self.calc_dir}")
                # increase nodes as its likely the calculation failed
                self.setup_calc(increase_walltime_by_factor=2, make_archive=True)
            return False

        self.logger.info(f"{self.mode.upper()} Calculation: reached required accuracy")
        self.logger.debug(tail_output)
        return True

    @property
    def is_done(self) -> bool:
        if getattr(self, "_is_done", None) is None:
            self._is_done = self.check_calc()
        return self._is_done

    @property
    def results(self) -> str:
        if not self.is_done:
            if self.stopped:
                self._results = "STOPPED"
            else:
                self._results = "not finished"
        else:
            self._results = "done"
        return self._results
