# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

from __future__ import annotations

import os
import shutil
from abc import ABC, abstractmethod
from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING

from pymatgen.core import Structure

from vasp_manager.config import load_computing_config
from vasp_manager.job_manager import JobManager
from vasp_manager.utils import get_pmg_structure_from_poscar
from vasp_manager.vasp_input_creator import VaspInputCreator
from vasp_manager.vasp_run import VaspRun

if TYPE_CHECKING:
    from vasp_manager.types import CalculationType, WorkingDirectory


class BaseCalculationManager(ABC):
    """Runs vasp job workflow for a single material"""

    def __init__(
        self,
        material_dir: WorkingDirectory,
        to_rerun: bool,
        to_submit: bool,
        primitive: bool = True,
        ignore_personal_errors: bool = True,
        from_scratch: bool = False,  # DANGEROUS, WILL DELETE PREVIOUS CALCULATION
    ):
        """Args:
        material_dir: path to a directory for a single material
            ex. calculations/AlAs
        to_rerun: if True, rerun failed calculations
        to_submit: if True, submit calculations to job manager
        primitive: if True, find primitive cell, else find conventional cell
        ignore_personal_errors: if True, ignore job submission errors
            if on personal computer
        from_scratch: if True, remove the calculation's directory and
            restart
            note: DANGEROUS
        """
        self.material_dir = Path(material_dir)
        self.to_rerun = to_rerun
        self.to_submit = to_submit
        self.primitive = primitive
        self._ignore_personal_errors = ignore_personal_errors

        self.from_scratch = from_scratch
        if from_scratch:
            self._from_scratch()

    @property
    @abstractmethod
    def mode(self) -> CalculationType:
        pass

    @property
    @abstractmethod
    def job_prefix(self) -> str:
        """Short prefix for SLURM job names (e.g. 'r', 's', 'b')."""

    @property
    @abstractmethod
    def is_done(self) -> bool:
        pass

    @cached_property
    @abstractmethod
    def poscar_source_path(self) -> Path:
        pass

    @cached_property
    def calc_dir(self) -> Path:
        return self.material_dir / self.mode

    @cached_property
    def config_dir(self) -> Path:
        return self.material_dir.parent

    @cached_property
    def material_name(self) -> str:
        return self.material_dir.name

    @abstractmethod
    def setup_calc(self) -> None:
        pass

    @abstractmethod
    def check_calc(self) -> bool:
        pass

    def _rerun_resource_kwargs(self) -> dict[str, int]:
        """Build increase_*_by_factor kwargs for a timeout rerun.

        Reads the rerun strategy from computing_config to determine whether
        to increase nodes or walltime, and by what factor.
        """
        cc = load_computing_config(self.config_dir)
        if cc.rerun_increase == "nodes":
            return {"increase_nodes_by_factor": cc.rerun_increase_factor}
        return {"increase_walltime_by_factor": cc.rerun_increase_factor}

    def _load_structure(self, poscar_path: Path | None = None) -> Structure:
        """Load a pymatgen Structure from a POSCAR file.

        Args:
            poscar_path: path to POSCAR/CONTCAR. Defaults to
                self.poscar_source_path.

        Returns:
            structure: pymatgen Structure
        """
        if poscar_path is None:
            poscar_path = self.poscar_source_path
        try:
            structure = get_pmg_structure_from_poscar(
                poscar_path, primitive=self.primitive
            )
        except Exception as e:
            raise Exception(f"Cannot load POSCAR in {poscar_path}: {e}")
        assert isinstance(structure, Structure)
        return structure

    def _load_structure_for_rerun(self) -> Structure:
        """Load structure for a rerun, checking for archives first.

        If archives exist in calc_dir, loads the CONTCAR from the latest
        archive. Otherwise falls back to poscar_source_path.
        """
        num_archives = len(list(self.calc_dir.glob("archive*")))
        if num_archives > 0:
            archive_name = f"archive_{num_archives - 1}"
            contcar_path = self.calc_dir / archive_name / "CONTCAR"
            return self._load_structure(contcar_path)
        return self._load_structure()

    @cached_property
    def vasp_runs(self) -> dict[str, VaspRun]:
        """All runs for this calculation. Override for multi-run managers."""
        structure = self._load_structure_for_rerun()
        vic = VaspInputCreator(
            calc_dir=self.calc_dir,
            mode=self.mode,
            structure=structure,
            config_dir=self.config_dir,
            name=self.material_name,
            job_prefix=self.job_prefix,
        )
        jm = JobManager(
            calc_dir=self.calc_dir,
            manager_name=f"{self.material_name} {self.mode.upper()}",
            config_dir=self.config_dir,
            ignore_personal_errors=self._ignore_personal_errors,
        )
        run = VaspRun(
            run_dir=self.calc_dir,
            structure=structure,
            vasp_input_creator=vic,
            job_manager=jm,
        )
        return {self.mode: run}

    @property
    def vasp_run(self) -> VaspRun:
        """Convenience accessor for single-run managers."""
        runs = self.vasp_runs
        if len(runs) != 1:
            raise AttributeError(
                f"{type(self).__name__} has {len(runs)} runs. "
                "Use vasp_runs to access individual runs."
            )
        return next(iter(runs.values()))

    def _invalidate_vasp_runs(self) -> None:
        """Clear the cached vasp_runs so it will be recreated on next access."""
        if "vasp_runs" in self.__dict__:
            del self.__dict__["vasp_runs"]

    @property
    def vasp_input_creator(self) -> VaspInputCreator:
        return self.vasp_run.vasp_input_creator

    @property
    def job_manager(self) -> JobManager:
        return self.vasp_run.job_manager

    @property
    def job_exists(self) -> bool:
        """True if all runs have a jobid file."""
        if not self.vasp_runs:
            return False
        return all(r.job_manager.job_exists for r in self.vasp_runs.values())

    @property
    def job_complete(self) -> bool:
        """True if all runs have finished (no longer in the SLURM queue)."""
        if not self.vasp_runs:
            return False
        return all(r.job_manager.job_complete for r in self.vasp_runs.values())

    def submit_job(self) -> bool:
        """Submit all runs to the job manager.

        Returns:
            True if all jobs were submitted successfully
        """
        if not self.vasp_runs:
            return False
        return all(r.job_manager.submit_job() for r in self.vasp_runs.values())

    @property
    def stopped(self) -> bool:
        """True if a STOP file exists in the material directory."""
        return (self.material_dir / "STOP").exists()

    def stop(self) -> None:
        """Create a STOP file to prevent further calculation attempts."""
        with open(self.material_dir / "STOP", "w+"):
            pass

    def _cancel_previous_job(self) -> None:
        """Cancel the SLURM job for this calculation and remove its jobid file."""
        jobid_path = self.calc_dir / "jobid"
        if jobid_path.exists():
            with open(jobid_path) as fr:
                jobid = fr.read().strip()
            cancel_job_call = f"scancel {jobid}"
            os.system(cancel_job_call)
            os.remove(jobid_path)

    def _from_scratch(self) -> None:
        """Cancel the existing job and delete the calculation directory."""
        self._cancel_previous_job()
        shutil.rmtree(self.calc_dir)
