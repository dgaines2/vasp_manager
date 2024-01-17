# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import logging
import os
import shutil
from functools import cached_property
from pathlib import Path

import numpy as np

from vasp_manager.analyzer.bulkmod_analyzer import BulkmodAnalyzer
from vasp_manager.calculation_manager.base import BaseCalculationManager
from vasp_manager.utils import change_directory, pgrep
from vasp_manager.vasp_input_creator import VaspInputCreator

logger = logging.getLogger(__name__)


class BulkmodCalculationManager(BaseCalculationManager):
    """
    Runs bulk modulus job workflow for a single material
    """

    def __init__(
        self,
        material_path,
        to_rerun,
        to_submit,
        primitive=True,
        ignore_personal_errors=True,
        from_relax=True,
        from_scratch=False,
        strains=None,
        tail=5,
    ):
        """
        For material_path, to_rerun, to_submit, ignore_personal_errors, and from_scratch,
        see BaseCalculationManager

        Args:
            from_relax (bool): if True, use CONTCAR from relax
            strains (array-like): optional, fractional strain along each axis for
                each deformation
                if None, use default
                default is np.linspace(start=0.925, stop=1.075, number=11)**(1/3)
                len(strains) must be odd and strains must be centered around 0
            tail (int): number of last lines to log in debugging if job failed
        """
        self.from_relax = from_relax
        self.strains = (
            strains
            if strains is not None
            else np.power(np.linspace(0.925, 1.075, 11), 1 / 3)
        )
        self.tail = tail
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
        return "bulkmod"

    @cached_property
    def poscar_source_path(self):
        if self.from_relax:
            poscar_source_path = self.material_path / "rlx" / "CONTCAR"
        else:
            poscar_source_path = self.material_path / "POSCAR"
        return poscar_source_path

    @cached_property
    def vasp_input_creator(self):
        return VaspInputCreator(
            self.calc_path,
            mode=self.mode,
            poscar_source_path=self.poscar_source_path,
            primitive=self.primitive,
            name=self.material_name,
        )

    @property
    def strains(self):
        return self._strains

    @strains.setter
    def strains(self, values):
        if np.any(values < 0.8) or np.any(values > 1.2):
            raise ValueError("Strains not in expected bounds")
        if (middle := values[int(len(values) / 2)]) != 1.0:
            raise ValueError(f"Strains not centered around 1.0: middle is {middle}")
        self._strains = values

    def _check_use_spin(self):
        if self.from_relax:
            rlx_stdout = self.material_path / "rlx" / "stdout.txt"
            rlx_mags = pgrep(rlx_stdout, "mag=", stop_after_first_match=True)
            use_spin = len(rlx_mags) != 0
        else:
            use_spin = True
        return use_spin

    def setup_calc(
        self,
        increase_nodes_by_factor=1,
        increase_walltime_by_factor=1,
    ):
        """
        Sets up an EOS bulkmod calculation
        """
        if not self.from_relax:
            msg = (
                "Running bulk modulus calculation without previous relaxation"
                "\n\t starting structure must be fairly close to equilibrium volume!"
            )
            logger.warning(msg)

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

    def check_calc(self):
        """
        Checks result of bulk modulus calculation

        Returns:
            bulkmod_sucessful (bool): if True, bulkmod calculation completed successfully
        """
        if not self.job_complete:
            logger.info(f"{self.mode.upper()} job not finished")
            return False

        for i, strain in enumerate(self.strains):
            middle = int(len(self.strains) / 2)
            strain_index = i - middle
            strain_name = f"strain_{strain_index}"
            strain_path = self.calc_path / strain_name
            stdout_path = strain_path / "stdout.txt"
            stderr_path = strain_path / "stderr.txt"
            if not stdout_path.exists():
                return False

            vasp_errors = self._check_vasp_errors(
                stdout_path=stdout_path, stderr_path=stderr_path, extra_errors=["NELM"]
            )
            if len(vasp_errors) > 0:
                all_errors_addressed = self._address_vasp_errors(vasp_errors)
                if all_errors_addressed:
                    if self.to_rerun:
                        logger.info(f"Rerunning {self.calc_path}")
                        self._from_scratch()
                        self.setup_calc()
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

            grep_output = pgrep(stdout_path, "1 F=", stop_after_first_match=True)
            if len(grep_output) == 0:
                if self.to_rerun:
                    logger.info(f"Rerunning {self.calc_path}")
                    # increase nodes as its likely the calculation failed
                    self._from_scratch()
                    self.setup_calc(increase_walltime_by_factor=2)
                return False
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
                return None
        try:
            ba = BulkmodAnalyzer(calc_path=self.calc_path)
            self._results = ba.results
        except Exception as e:
            logger.warning(e)
            self._results = None
        return self._results

    def _make_bulkmod_strains(self):
        """
        Creates a set of strain directory for fitting the E-V info

        Args:
            strains (iterable of floats)
        """
        logger.info("Making strain directories")
        for i, strain in enumerate(self.strains):
            middle = int(len(self.strains) / 2)
            strain_index = i - middle
            strain_name = f"strain_{strain_index}"
            strain_path = self.calc_path / strain_name
            logger.info(strain_path)

            if not strain_path.exists():
                strain_path.mkdir()
            orig_poscar_path = self.calc_path / "POSCAR"
            strain_poscar_path = strain_path / "POSCAR"
            shutil.copy(orig_poscar_path, strain_poscar_path)

            # change second line to be {strain} rather than 1.0
            with open(strain_poscar_path, "r") as fr:
                strain_poscar = fr.read().splitlines()
            strain_poscar[1] = f"{strain}"
            logger.debug(f"{strain}")
            with open(strain_poscar_path, "w+") as fw:
                fw.write("\n".join(strain_poscar))

            with change_directory(strain_path):
                for f in ["POTCAR", "INCAR"]:
                    if Path(f).exists():
                        os.remove(f)
                    orig_path = Path("..") / f
                    os.symlink(orig_path, f, target_is_directory=False)
