# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

from __future__ import annotations

import logging
import os
import shutil
from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from vasp_manager.analyzer.bulkmod_analyzer import BulkmodAnalyzer
from vasp_manager.calculation_manager.base import BaseCalculationManager
from vasp_manager.utils import LoggerAdapter, change_directory, pgrep
from vasp_manager.vasp_input_creator import VaspInputCreator

if TYPE_CHECKING:
    from vasp_manager.types import CalculationType, WorkingDirectory

logger = logging.getLogger(__name__)


class BulkmodCalculationManager(BaseCalculationManager):
    """
    Runs bulk modulus job workflow for a single material
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
        strains: None | NDArray = None,
    ):
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
            strains: fractional strain along each axis for each deformation
                if None, use np.linspace(start=0.925, stop=1.075, number=11)**(1/3)
                len(strains) must be odd and strains must be centered around 0
        """
        self.from_relax = from_relax
        self.strains = (
            strains
            if strains is not None
            else np.power(np.linspace(0.925, 1.075, 11), 1 / 3)
        )
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
        return "bulkmod"

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

    @property
    def strains(self) -> NDArray:
        return self._strains

    @strains.setter
    def strains(self, values: NDArray) -> None:
        if np.any(values < 0.8) or np.any(values > 1.2):
            raise ValueError("Strains not in expected bounds")
        if (middle := values[int(len(values) / 2)]) != 1.0:
            raise ValueError(f"Strains not centered around 1.0: middle is {middle}")
        self._strains = values

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
        Sets up an EOS bulkmod calculation
        """
        if not self.from_relax:
            msg = (
                "Running bulk modulus calculation without previous relaxation\n"
                "\tStarting structure must be fairly close to equilibrium volume!"
            )
            self.logger.warning(msg)

        self.vasp_input_creator.increase_nodes_by_factor = increase_nodes_by_factor
        self.vasp_input_creator.increase_walltime_by_factor = increase_walltime_by_factor
        self.vasp_input_creator.use_spin = self._check_use_spin()
        self.vasp_input_creator.create()
        self._make_bulkmod_strains()

        if self.to_submit:
            job_submitted = self.submit_job()
            # job status returns True if sucessfully submitted, else False
            if not job_submitted:
                self.setup_calc()

    def check_calc(self) -> bool:
        """
        Checks result of bulk modulus calculation

        Returns:
            bulkmod_sucessful: if True, bulkmod calculation completed successfully
        """
        if not self.job_complete:
            self.logger.info(f"{self.mode.upper()} job not finished")
            return False

        for i, strain in enumerate(self.strains):
            middle = int(len(self.strains) / 2)
            strain_index = i - middle
            strain_name = f"strain_{strain_index}"
            strain_dir = self.calc_dir / strain_name
            stdout_path = strain_dir / "stdout.txt"
            stderr_path = strain_dir / "stderr.txt"
            if not stdout_path.exists():
                return False

            vasp_errors = self._check_vasp_errors(
                stdout_path=stdout_path, stderr_path=stderr_path, extra_errors=["NELM"]
            )
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

            grep_output = pgrep(stdout_path, "1 F=", stop_after_first_match=True)
            if len(grep_output) == 0:
                if self.to_rerun:
                    self.logger.info(f"Rerunning {self.calc_dir}")
                    # increase nodes as its likely the calculation failed
                    self._from_scratch()
                    self.setup_calc(increase_nodes_by_factor=2)
                return False
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
        try:
            ba = BulkmodAnalyzer(calc_dir=self.calc_dir)
            self._results = ba.results
            self.logger.info(f"{self.mode.upper()} Calculation: Success")
            self.logger.info(f"BULK MODULUS: {ba.results.get('B')}")
        except Exception as e:
            self.logger.warning(e)
            self._results = None
        return self._results

    def _make_bulkmod_strains(self) -> None:
        """
        Creates a set of strain directories for fitting the E-V info
        """
        self.logger.info("Making strain directories")
        for i, strain in enumerate(self.strains):
            middle = int(len(self.strains) / 2)
            strain_index = i - middle
            strain_name = f"strain_{strain_index}"
            strain_dir = self.calc_dir / strain_name
            self.logger.info(strain_dir)

            if not strain_dir.exists():
                strain_dir.mkdir()
            orig_poscar_path = self.calc_dir / "POSCAR"
            strain_poscar_path = strain_dir / "POSCAR"
            shutil.copy(orig_poscar_path, strain_poscar_path)

            # change second line to be {strain} rather than 1.0
            with open(strain_poscar_path, "r") as fr:
                strain_poscar = fr.read().splitlines()
            strain_poscar[1] = f"{strain}"
            self.logger.debug(f"{strain}")
            with open(strain_poscar_path, "w+") as fw:
                fw.write("\n".join(strain_poscar))

            with change_directory(strain_dir):
                for f in ["POTCAR", "INCAR"]:
                    if Path(f).exists():
                        os.remove(f)
                    orig_path = Path("..") / f
                    os.symlink(orig_path, f, target_is_directory=False)
