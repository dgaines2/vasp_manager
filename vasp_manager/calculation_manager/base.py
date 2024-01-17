# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import os
import shutil
from abc import ABC, abstractmethod
from functools import cached_property
from pathlib import Path

import numpy as np

from vasp_manager.job_manager import JobManager
from vasp_manager.utils import get_pmg_structure_from_poscar, pgrep


class BaseCalculationManager(ABC):
    """
    Runs vasp job workflow for a single material
    """

    def __init__(
        self,
        material_path,
        to_rerun,
        to_submit,
        primitive=True,
        ignore_personal_errors=True,
        from_scratch=False,  # DANGEROUS, WILL DELETE PREVIOUS CALCULATION
    ):
        """
        Args:
            material_path (str | Path): path for a single material
                ex. calculations/AlAs
            to_rerun (bool): if True, rerun failed calculations
            to_submit (bool): if True, submit calculations to job manager
            primitive (bool): if True, find primitive cell, else find conventional cell
            ignore_personal_errors (bool): if True, ignore job submission errors
                if on personal computer
            from_scratch (bool): if True, remove the calculation's folder and
                restart
                note: DANGEROUS
        """
        self.material_path = Path(material_path)
        self.to_rerun = to_rerun
        self.to_submit = to_submit
        self.primitive = primitive
        self.job_manager = JobManager(
            calc_path=self.calc_path, ignore_personal_errors=ignore_personal_errors
        )

        self.from_scratch = from_scratch
        if from_scratch:
            self._from_scratch()

    @property
    @abstractmethod
    def mode(self):
        pass

    @property
    @abstractmethod
    def is_done(self):
        pass

    @property
    @abstractmethod
    def vasp_input_creator(self):
        pass

    @cached_property
    def calc_path(self):
        return self.material_path / self.mode

    @cached_property
    def material_name(self):
        return self.material_path.name

    @abstractmethod
    def poscar_source_path(self):
        pass

    @abstractmethod
    def setup_calc(self):
        pass

    @abstractmethod
    def check_calc(self):
        pass

    @property
    def job_exists(self):
        return self.job_manager.job_exists

    @property
    def job_complete(self):
        return self.job_manager.job_complete

    @property
    def stopped(self):
        return (self.material_path / "STOP").exists()

    def stop(self):
        with open(self.material_path / "STOP", "w+"):
            pass

    def submit_job(self):
        return self.job_manager.submit_job()

    def _cancel_previous_job(self):
        jobid_path = self.calc_path / "jobid"
        if jobid_path.exists():
            with open(jobid_path) as fr:
                jobid = fr.read().strip()
            cancel_job_call = f"scancel {jobid}"
            os.system(cancel_job_call)
            os.remove(jobid_path)

    def _from_scratch(self):
        self._cancel_previous_job()
        shutil.rmtree(self.calc_path)

    def _parse_magmom(self):
        stdout_path = self.calc_path / "stdout.txt"
        mag_lines = pgrep(stdout_path, "mag=")
        # if "mag=" not found in stdout, set magmom=None
        if len(mag_lines) == 0:
            magmom = None
        else:
            total_mag = mag_lines[-1].split()[-1]
            magmom = float(total_mag)
        return magmom

    def _parse_magmom_per_atom(self):
        total_magmom = self._parse_magmom()
        if total_magmom is None:
            return None
        structure = get_pmg_structure_from_poscar(self.poscar_source_path)
        magmom_per_atom = total_magmom / len(structure)
        return magmom_per_atom

    def _parse_incar_tag(self, tag):
        incar_path = self.calc_path / "INCAR"
        with open(incar_path) as fr:
            incar = fr.readlines()
        tag_value = None
        for line in incar:
            if tag in line:
                tag_value = line.split("=")[1].strip()
        return tag_value

    def _check_vasp_errors(self, stdout_path=None, stderr_path=None, extra_errors=None):
        """
        Find VASP errors in stdout and stderr
        """
        if stdout_path is None:
            stdout_path = self.calc_path / "stdout.txt"
        if stderr_path is None:
            stderr_path = self.calc_path / "stderr.txt"
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

    def _address_vasp_errors(self, errors):
        """
        Args:
            errors (set): set of errors found in stdout or stderr
        Returns:
            all_errors_addressed (bool): if True, all errors could
                be fixed automatically. If False, some errors
                could not be handled
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
                        # total of 16 (as vic adds 4 if on quest)
                        ncore_per_node_for_memory = 12
                    else:
                        ncore_per_node_for_memory = 64
                    vic.ncore_per_node_for_memory = ncore_per_node_for_memory
                    errors_addressed[error] = True
                case _:
                    errors_addressed[error] = False
        all_errors_addressed = np.all([v for v in errors_addressed.values()])
        return all_errors_addressed
