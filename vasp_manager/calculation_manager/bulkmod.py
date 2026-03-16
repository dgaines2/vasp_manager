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
from vasp_manager.job_manager import JobManager
from vasp_manager.utils import LoggerAdapter, pgrep
from vasp_manager.vasp_input_creator import VaspInputCreator
from vasp_manager.vasp_run import VaspRun

if TYPE_CHECKING:
    from vasp_manager.types import CalculationType, WorkingDirectory

logger = logging.getLogger(__name__)


class BulkmodCalculationManager(BaseCalculationManager):
    """
    Runs bulk modulus job workflow for a single material.

    Each strain is an independent static calculation with its own
    VaspInputCreator, JobManager, and SLURM job.
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
        self._vasp_runs: dict[str, VaspRun] | None = None
        self.logger = LoggerAdapter(logging.getLogger(__name__), self.material_name)

    @cached_property
    def mode(self) -> CalculationType:
        return "bulkmod"

    @property
    def job_prefix(self) -> str:
        return "b"

    @cached_property
    def poscar_source_path(self) -> Path:
        if self.from_relax:
            poscar_source_path = self.material_dir / "rlx" / "CONTCAR"
        else:
            poscar_source_path = self.material_dir / "POSCAR"
        return poscar_source_path

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

    @property
    def strain_names(self) -> list[str]:
        middle = int(len(self.strains) / 2)
        return [f"strain_{i - middle}" for i in range(len(self.strains))]

    @property
    def vasp_runs(self) -> dict[str, VaspRun]:
        """One VaspRun per strain directory."""
        if self._vasp_runs is None:
            self._vasp_runs = self._build_vasp_runs()
        return self._vasp_runs

    def _build_vasp_runs(self) -> dict[str, VaspRun]:
        """Build VaspRun dict from existing strain directories."""
        runs = {}
        base_structure = self._load_structure()
        for i, strain in enumerate(self.strains):
            strain_name = self.strain_names[i]
            strain_dir = self.calc_dir / strain_name
            if not strain_dir.exists():
                continue
            strained = base_structure.copy()
            strained.scale_lattice(base_structure.volume * strain**3)
            vic = VaspInputCreator(
                calc_dir=strain_dir,
                mode="static",
                structure=strained,
                config_dir=self.config_dir,
                name=f"{self.material_name}_s{strain_name.split('_')[1]}",
                job_prefix=self.job_prefix,
            )
            jm = JobManager(
                calc_dir=strain_dir,
                manager_name=f"{self.material_name} BULKMOD {strain_name}",
                config_dir=self.config_dir,
                ignore_personal_errors=self._ignore_personal_errors,
            )
            runs[strain_name] = VaspRun(
                run_dir=strain_dir,
                structure=strained,
                vasp_input_creator=vic,
                job_manager=jm,
            )
        return runs

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
        Sets up an EOS bulkmod calculation with independent strain jobs.
        """
        if not self.from_relax:
            msg = (
                "Running bulk modulus calculation without previous relaxation\n"
                "\tStarting structure must be fairly close to equilibrium volume!"
            )
            self.logger.warning(msg)

        if not self.calc_dir.exists():
            self.calc_dir.mkdir()

        base_structure = self._load_structure()
        use_spin = self._check_use_spin()
        runs: dict[str, VaspRun] = {}

        self.logger.info("Making strain directories")
        for i, strain in enumerate(self.strains):
            strain_name = self.strain_names[i]
            strain_dir = self.calc_dir / strain_name
            self.logger.info(strain_dir)

            if not strain_dir.exists():
                strain_dir.mkdir()

            # Scale structure for this strain
            strained = base_structure.copy()
            strained.scale_lattice(base_structure.volume * strain**3)

            vic = VaspInputCreator(
                calc_dir=strain_dir,
                mode="static",
                structure=strained,
                config_dir=self.config_dir,
                name=f"{self.material_name}_s{strain_name.split('_')[1]}",
                job_prefix=self.job_prefix,
                increase_nodes_by_factor=increase_nodes_by_factor,
                increase_walltime_by_factor=increase_walltime_by_factor,
                use_spin=use_spin,
            )
            vic.create()

            jm = JobManager(
                calc_dir=strain_dir,
                manager_name=f"{self.material_name} BULKMOD {strain_name}",
                config_dir=self.config_dir,
                ignore_personal_errors=self._ignore_personal_errors,
            )

            runs[strain_name] = VaspRun(
                run_dir=strain_dir,
                structure=strained,
                vasp_input_creator=vic,
                job_manager=jm,
            )

        self._vasp_runs = runs

        if self.to_submit:
            all_submitted = self.submit_job()
            if not all_submitted:
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

        for strain_name, run in self.vasp_runs.items():
            stdout_path = run.run_dir / "stdout.txt"
            if not stdout_path.exists():
                return False

            vasp_errors = run.check_vasp_errors(extra_errors=["NELM"])
            if len(vasp_errors) > 0:
                all_errors_addressed = run.address_vasp_errors(vasp_errors)
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
                    self.setup_calc(increase_walltime_by_factor=2)
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

    def _cancel_previous_job(self) -> None:
        """Cancel all strain jobs."""
        for strain_name in self.strain_names:
            strain_dir = self.calc_dir / strain_name
            jobid_path = strain_dir / "jobid"
            if jobid_path.exists():
                with open(jobid_path) as fr:
                    jobid = fr.read().strip()
                cancel_job_call = f"scancel {jobid}"
                os.system(cancel_job_call)
                os.remove(jobid_path)

    def _from_scratch(self) -> None:
        self._cancel_previous_job()
        shutil.rmtree(self.calc_dir)
        self._vasp_runs = None
