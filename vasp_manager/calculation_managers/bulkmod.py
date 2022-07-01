# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import glob
import logging
import os
import shutil

import numpy as np
import pymatgen as pmg

from vasp_manager.calculation_managers.base import BaseCalculationManager
from vasp_manager.utils import change_directory
from vasp_manager.vasp_input_creator import VaspInputCreator

logger = logging.getLogger(__name__)


class BulkmodCalculationManager(BaseCalculationManager):
    def __init__(
        self,
        base_path,
        to_rerun,
        to_submit,
        ignore_personal_errors=True,
        from_relax=True,
        from_scratch=False,
    ):
        super().__init__(
            base_path=base_path,
            to_rerun=to_rerun,
            to_submit=to_submit,
            ignore_personal_errors=ignore_personal_errors,
            from_scratch=from_scratch,
        )
        self.from_relax = from_relax

    @property
    def mode(self):
        mode_str = "bulkmod"
        if not self.from_relax:
            mode_str += "_standalone"
        return mode_str

    @property
    def poscar_source_path(self):
        if self.from_relax:
            poscar_source_path = os.path.join(self.base_path, "rlx", "CONTCAR")
        else:
            poscar_source_path = os.path.join(self.base_path, "POSCAR")
        return poscar_source_path

    def setup_calc(self):
        """
        Set up EOS bulkmod calculation
        """
        vasp_input_creator = VaspInputCreator(
            self.calc_path,
            mode=self.mode,
            poscar_source_path=self.poscar_source_path,
            name=self.material_name,
        )
        vasp_input_creator.create()

        strains = np.power(np.linspace(0.8, 1.2, 11), 1 / 3)
        self._make_bulkmod_strains(strains)

        if self.to_submit:
            job_submitted = self.submit_job()
            # job status returns True if sucessfully submitted, else False
            if not job_submitted:
                self.setup_calc()

    def check_calc(self):
        # //TODO: make sure results look reasonable here
        # perhaps call self.bulk_modulus to make sure
        # everything parses correctly
        if not self.job_complete:
            logger.info(f"{self.mode.upper()} job not finished")
            return False
        else:
            return True

    @property
    def is_done(self):
        return self.check_calc()

    def _make_bulkmod_strains(self, strains):
        """
        Create a set of strain directory for fitting the E-V info

        Args:
            strains (iterable of floats)
        """
        logger.info("Making strain directories")
        for i, strain in enumerate(strains):
            middle = int(len(strains) / 2)
            strain_index = i - middle
            strain_name = f"strain_{strain_index}"
            strain_path = os.path.join(self.calc_path, strain_name)
            logger.info(strain_path)

            if not os.path.exists(strain_path):
                os.mkdir(strain_path)
            orig_poscar_path = os.path.join(self.calc_path, "POSCAR")
            strain_poscar_path = os.path.join(strain_path, "POSCAR")
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
                    orig_path = "../" + f
                    f_path = f
                    os.symlink(orig_path, f_path, target_is_directory=False)


def _analyze_bulkmod(self):
    """
    Fit an EOS to calculate the bulk modulus from a finished bulkmod calculation
    """
    strain_paths = [
        path for path in glob.glob(self.calc_path + "/strain*") if os.path.isdir(path)
    ]
    strain_paths = sorted(strain_paths, key=lambda d: int(d.split("_")[-1]))
    volumes = []
    final_energies = []
    for i, strain_path in enumerate(strain_paths):
        poscar_path = os.path.join(strain_path, "POSCAR")
        vasprun_path = os.path.join(strain_path, "vasprun.xml")
        volume = pmg.core.Structure.from_file(poscar_path).volume
        vasprun = pmg.io.vasp.outputs.Vasprun(
            filename=vasprun_path,
            parse_dos=False,
            parse_eigen=False,
            parse_potcar_file=False,
        )
        final_energy = vasprun.final_energy
        volumes.append(volume)
        final_energies.append(final_energy)
    logger.debug("Volumes")
    logger.debug(f"{volumes}")
    logger.debug("Final Energies")
    logger.debug(f"{final_energies}")
    eos_analyzer = pmg.analysis.eos.BirchMurnaghan(volumes, final_energies)
    eos_analyzer.fit()
    bulk_modulus = np.round(eos_analyzer.b0_GPa, 3)
    logger.info(f"{self.mode.upper()} Calculation: Successful")
    logger.info(f"BULK MODULUS: {bulk_modulus}")
    return bulk_modulus

    @property
    def bulk_modulus(self):
        return self._analyze_bulkmod(self)
