# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

from __future__ import annotations

import os
import shutil
from abc import ABC, abstractmethod
from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING

from vasp_manager.job_manager import JobManager
from vasp_manager.utils import get_pmg_structure_from_poscar, pgrep
from vasp_manager.vasp_input_creator import VaspInputCreator

if TYPE_CHECKING:
    from vasp_manager.types import CalculationType, WorkingDirectory


class BaseCalculationManager(ABC):
    """
    Runs vasp job workflow for a single material
    """

    def __init__(
        self,
        material_dir: WorkingDirectory,
        to_rerun: bool,
        to_submit: bool,
        primitive: bool = True,
        ignore_personal_errors: bool = True,
        from_scratch: bool = False,  # DANGEROUS, WILL DELETE PREVIOUS CALCULATION
    ):
        """
        Args:
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
        self.job_manager = JobManager(
            calc_dir=self.calc_dir,
            manager_name=f"{self.material_name} {self.mode.upper()}",
            ignore_personal_errors=ignore_personal_errors,
        )

        self.from_scratch = from_scratch
        if from_scratch:
            self._from_scratch()

    @property
    @abstractmethod
    def mode(self) -> CalculationType:
        pass

    @property
    @abstractmethod
    def is_done(self) -> bool:
        pass

    @property
    @abstractmethod
    def vasp_input_creator(self) -> VaspInputCreator:
        pass

    @cached_property
    def calc_dir(self) -> Path:
        return self.material_dir / self.mode

    @cached_property
    def material_name(self) -> str:
        return self.material_dir.name

    @cached_property
    @abstractmethod
    def poscar_source_path(self) -> Path:
        pass

    @abstractmethod
    def setup_calc(self) -> None:
        pass

    @abstractmethod
    def check_calc(self) -> bool:
        pass

    @property
    def job_exists(self) -> bool:
        return self.job_manager.job_exists

    @property
    def job_complete(self) -> bool:
        return self.job_manager.job_complete

    @property
    def stopped(self) -> bool:
        return (self.material_dir / "STOP").exists()

    def stop(self) -> None:
        with open(self.material_dir / "STOP", "w+"):
            pass

    def submit_job(self) -> bool:
        return self.job_manager.submit_job()

    def _cancel_previous_job(self) -> None:
        jobid_path = self.calc_dir / "jobid"
        if jobid_path.exists():
            with open(jobid_path) as fr:
                jobid = fr.read().strip()
            cancel_job_call = f"scancel {jobid}"
            os.system(cancel_job_call)
            os.remove(jobid_path)

    def _from_scratch(self) -> None:
        self._cancel_previous_job()
        shutil.rmtree(self.calc_dir)

    def _parse_magmom(self) -> float | None:
        stdout_path = self.calc_dir / "stdout.txt"
        mag_lines = pgrep(stdout_path, "mag=")
        # if "mag=" not found in stdout, set magmom=None
        if len(mag_lines) == 0:
            magmom = None
        else:
            total_mag = mag_lines[-1].split()[-1]
            magmom = float(total_mag)
        return magmom

    def _parse_magmom_per_atom(self) -> float | None:
        total_magmom = self._parse_magmom()
        if total_magmom is None:
            return None
        structure = get_pmg_structure_from_poscar(self.poscar_source_path)
        magmom_per_atom = total_magmom / len(structure)
        return magmom_per_atom

    def _parse_incar_tag(self, tag) -> str:
        incar_path = self.calc_dir / "INCAR"
        with open(incar_path) as fr:
            incar = fr.readlines()
        for line in incar:
            if tag in line:
                tag_value = line.split("=")[1].strip()
        return tag_value

    def _check_vasp_errors(
        self,
        stdout_path: Path | None = None,
        stderr_path: Path | None = None,
        extra_errors: list[str] | None = None,
    ) -> set[str]:
        """
        Find VASP errors in stdout and stderr

        Args:
            stdout_path
            stderr_path
            extra_errors: names of other errors to include in search

        Returns:
            errors_found: True if any errors found
        """
        if stdout_path is None:
            stdout_path = self.calc_dir / "stdout.txt"
        if stderr_path is None:
            stderr_path = self.calc_dir / "stderr.txt"
        if extra_errors is None:
            extra_errors = []
        errors_found = set()

        with open(stdout_path) as fr:
            for line in fr:
                errors = [
                    "Sub-Space-Matrix",
                    "Inconsistent Bravais",
                    "num prob",
                    "BRMIX",
                    "SICK JOB",
                    "VERY BAD NEWS",
                    "Fatal error",
                ]
                errors.extend(extra_errors)
                for error in errors:
                    if error in line:
                        errors_found.add(error)

        with open(stderr_path) as fr:
            for line in fr:
                for error in [
                    "oom-kill",
                    "SETYLM",
                    "Segmentation",
                    "command not found",
                ]:
                    if error in line:
                        errors_found.add(error)

        return errors_found

    def _address_vasp_errors(self, errors: set[str]) -> bool:
        """
        Args:
            errors: set of errors found in stdout or stderr

        Returns:
            all_errors_addressed: if True, all errors could be fixed
                automatically. If False, some errors could not be handled
        """
        vic = self.vasp_input_creator

        errors_addressed = {e: False for e in errors}
        for error in errors:
            match error:
                case "Sub-Space-Matrix":
                    new_algo = "Fast"
                    previous_algo = self._parse_incar_tag("ALGO")
                    if previous_algo == new_algo:
                        errors_addressed[error] = False
                    else:
                        vic.calc_config["algo"] = new_algo
                        errors_addressed[error] = True
                case "Inconsistent Bravais":
                    new_symprec = "1e-08"
                    previous_symprec = self._parse_incar_tag("SYMPREC")
                    if previous_symprec == new_symprec:
                        errors_addressed[error] = False
                    else:
                        vic.calc_config["symprec"] = new_symprec
                        errors_addressed[error] = True
                case "oom-kill":
                    if vic.computer == "quest":
                        ncore_per_node_for_memory = 28
                    else:
                        ncore_per_node_for_memory = 64
                    vic.ncore_per_node_for_memory = ncore_per_node_for_memory
                    errors_addressed[error] = True
                case _:
                    errors_addressed[error] = False
        all_errors_addressed = all([v for v in errors_addressed.values()])
        return all_errors_addressed
