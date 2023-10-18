# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

from __future__ import annotations

import logging
from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING

from pymatgen.core import Structure

from vasp_manager.calculation_manager.base import BaseCalculationManager
from vasp_manager.utils import LoggerAdapter, pgrep, ptail
from vasp_manager.vasp_input_creator import VaspInputCreator

if TYPE_CHECKING:
    from vasp_manager.types import CalculationType, WorkingDirectory

logger = logging.getLogger(__name__)


class StaticCalculationManager(BaseCalculationManager):
    """
    Runs static job workflow for a single material
    """

    def __init__(
        self,
        material_dir: WorkingDirectory,
        to_rerun: bool,
        to_submit: bool,
        primitive: bool = True,
        ignore_personal_errors: bool = True,
        from_scratch: bool = False,
        from_relax: bool = True,
        tail: int = 5,
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
            from_relax: if True, use CONTCAR from relax
            tail: number of last lines to log in debugging if job failed
        """
        self.from_relax = from_relax
        self.tail = tail
        super().__init__(
            material_dir=material_dir,
            to_rerun=to_rerun,
            to_submit=to_submit,
            primitive=primitive,
            ignore_personal_errors=ignore_personal_errors,
            from_scratch=from_scratch,
        )
        self._is_done: bool
        self._results: None | str | dict
        self.logger = LoggerAdapter(logging.getLogger(__name__), self.material_name)

    @cached_property
    def mode(self) -> CalculationType:
        return "static"

    @cached_property
    def poscar_source_path(self) -> Path:
        if self.from_relax:
            poscar_source_path = self.material_dir / "rlx" / "CONTCAR"
        else:
            poscar_source_path = self.material_dir / "POSCAR"
        return poscar_source_path

    @cached_property
    def vasp_input_creator(self) -> VaspInputCreator:
        return VaspInputCreator(
            self.calc_dir,
            mode=self.mode,
            poscar_source_path=self.poscar_source_path,
            primitive=self.primitive,
            name=self.material_name,
        )

    def _check_use_spin(self) -> bool:
        if self.from_relax:
            rlx_stdout = self.material_dir / "rlx" / "stdout.txt"
            rlx_mags = pgrep(rlx_stdout, "mag=", stop_after_first_match=True)
            use_spin = len(rlx_mags) != 0
        else:
            use_spin = True
        return use_spin

    def setup_calc(
        self,
        increase_nodes_by_factor: int = 1,
        increase_walltime_by_factor: int = 1,
    ) -> None:
        """
        Runs a static SCF calculation through VASP

        By default, requires previous relaxation run
        """
        self.vasp_input_creator.increase_nodes_by_factor = increase_nodes_by_factor
        self.vasp_input_creator.increase_walltime_by_factor = increase_walltime_by_factor
        self.vasp_input_creator.use_spin = self._check_use_spin()
        self.vasp_input_creator.create()

        if self.to_submit:
            job_submitted = self.submit_job()
            # job status returns True if sucessfully submitted, else False
            if not job_submitted:
                self.setup_calc()

    def check_calc(self) -> bool:
        """
        Checks result of static calculation

        Returns
            static_successful: if True, static calculation completed successfully
        """
        if not self.job_complete:
            self.logger.info(f"{self.mode.upper()} not finished")
            return False

        stdout_path = self.calc_dir / "stdout.txt"
        if not stdout_path.exists():
            # shouldn't get here unless function was called with submit=False
            self.logger.info(f"{self.mode.upper()} Calculation: No stdout.txt available")
            if self.to_rerun:
                self._from_scratch()
                self.setup_calc()
            return False

        vasp_errors = self._check_vasp_errors(extra_errors=["NELM"])
        if len(vasp_errors) > 0:
            all_errors_addressed = self._address_vasp_errors(vasp_errors)
            if all_errors_addressed:
                if self.to_rerun:
                    self.logger.info(f"Rerunning {self.calc_dir}")
                    self._from_scratch()
                    self.setup_calc()
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
        grep_output = pgrep(stdout_path, "1 F=", stop_after_first_match=True)
        if len(grep_output) == 0:
            self.logger.warning(f"{self.mode.upper()} FAILED")
            self.logger.debug(tail_output)
            if self.to_rerun:
                self.logger.info(f"Rerunning {self.calc_dir}")
                self._from_scratch()
                # increase nodes as its likely the calculation failed
                self.setup_calc(increase_nodes_by_factor=2)
            return False

        self._results = {}
        final_energy = float(grep_output[0].split()[2])
        num_atoms = len(Structure.from_file(self.calc_dir / "POSCAR"))
        magmom_per_atom = self._parse_magmom_per_atom()
        self._results["final_energy"] = final_energy
        self._results["final_energy_pa"] = final_energy / num_atoms
        self._results["magmom_pa"] = magmom_per_atom
        self.logger.info(f"{self.mode.upper()} Calculation: SCF converged")
        self.logger.debug(tail_output)
        return True

    @property
    def is_done(self) -> bool:
        if getattr(self, "_is_done", None) is None:
            self._is_done = self.check_calc()
        return self._is_done

    @property
    def results(self) -> None | str | dict:
        if not self.is_done:
            if self.stopped:
                return "STOPPED"
            else:
                return None
        return self._results
