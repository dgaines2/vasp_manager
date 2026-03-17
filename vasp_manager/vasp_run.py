# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from pymatgen.core import Structure

from vasp_manager.utils import pgrep

if TYPE_CHECKING:
    from vasp_manager.job_manager import JobManager
    from vasp_manager.types import WorkingDirectory
    from vasp_manager.vasp_input_creator import VaspInputCreator

STDOUT_ERRORS = (
    "Sub-Space-Matrix",
    "Inconsistent Bravais",
    "num prob",
    "BRMIX",
    "SICK JOB",
    "VERY BAD NEWS",
    "Fatal error",
)

STDERR_ERRORS = (
    "oom-kill",
    "SETYLM",
    "Segmentation",
    "command not found",
)


class VaspRun:
    """A single VASP execution: one directory, one structure, one input set, one job."""

    def __init__(
        self,
        run_dir: WorkingDirectory,
        structure: Structure,
        vasp_input_creator: VaspInputCreator,
        job_manager: JobManager,
    ) -> None:
        self.run_dir = Path(run_dir)
        self.structure = structure
        self.vasp_input_creator = vasp_input_creator
        self.job_manager = job_manager

    def check_vasp_errors(
        self,
        extra_errors: list[str] | None = None,
    ) -> set[str]:
        """
        Find VASP errors in stdout and stderr

        Args:
            extra_errors: names of other errors to include in search

        Returns:
            errors_found: set of error strings found
        """
        stdout_path = self.run_dir / "stdout.txt"
        stderr_path = self.run_dir / "stderr.txt"
        if extra_errors is None:
            extra_errors = []
        errors_found = set()

        with open(stdout_path) as fr:
            for line in fr:
                errors = list(STDOUT_ERRORS) + extra_errors
                for error in errors:
                    if error in line:
                        errors_found.add(error)

        with open(stderr_path) as fr:
            for line in fr:
                for error in STDERR_ERRORS:
                    if error in line:
                        errors_found.add(error)

        return errors_found

    def address_vasp_errors(self, errors: set[str]) -> bool:
        """
        Try to automatically fix known VASP errors by modifying VIC settings.

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
                        vic.calc_config = vic.calc_config.model_copy(
                            update={"algo": new_algo}
                        )
                        errors_addressed[error] = True
                case "Inconsistent Bravais":
                    new_symprec = "1e-08"
                    previous_symprec = self._parse_incar_tag("SYMPREC")
                    if previous_symprec == new_symprec:
                        errors_addressed[error] = False
                    else:
                        vic.calc_config = vic.calc_config.model_copy(
                            update={"symprec": new_symprec}
                        )
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

    def parse_magmom(self) -> float | None:
        """Parse total magnetization from stdout.txt"""
        stdout_path = self.run_dir / "stdout.txt"
        mag_lines = pgrep(stdout_path, "mag=")
        if len(mag_lines) == 0:
            return None
        total_mag = mag_lines[-1].split()[-1]
        return float(total_mag)

    def parse_magmom_per_atom(self) -> float | None:
        """Parse magnetization per atom from stdout.txt"""
        total_magmom = self.parse_magmom()
        if total_magmom is None:
            return None
        return total_magmom / len(self.structure)

    def _parse_incar_tag(self, tag: str) -> str:
        """Read a tag value from the INCAR file in calc_dir"""
        incar_path = self.run_dir / "INCAR"
        with open(incar_path) as fr:
            incar = fr.readlines()
        for line in incar:
            if tag in line:
                tag_value = line.split("=")[1].strip()
        return tag_value
