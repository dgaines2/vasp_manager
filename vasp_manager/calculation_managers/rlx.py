# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import glob
import logging
import os
import shutil
import subprocess

import numpy as np
import pymatgen as pmg
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from vasp_manager.calculation_managers.base import BaseCalculationManager
from vasp_manager.utils import get_pmg_structure_from_poscar
from vasp_manager.vasp_utils import make_archive, make_incar, make_potcar, make_vaspq

logger = logging.getLogger(__name__)


class RlxCalculationManager(BaseCalculationManager):
    def __init__(
        self,
        base_path,
        to_rerun,
        to_submit,
        ignore_personal_errors=True,
        from_coarse_relax=True,
        tail=5,
        from_scratch=False,
    ):
        super().__init__(
            base_path=base_path,
            to_rerun=to_rerun,
            to_submit=to_submit,
            ignore_personal_errors=ignore_personal_errors,
            from_scratch=from_scratch,
        )
        self.from_coarse_relax = from_coarse_relax
        self.tail = tail

    @property
    def mode(self):
        return "rlx"

    def setup_calc(self):
        """
        Set up a fine relaxation
        """
        if not os.path.exists(self.calc_path):
            os.mkdir(self.calc_path)

        if self.to_rerun:
            archive_made = make_archive(self.base_path, mode=self.mode)
            if not archive_made:
                # set rerun to not make an achive and instead
                # continue to make the input files
                self.to_rerun = False
                self.setup_calc()
        else:
            # POSCAR, POTCAR, INCAR, vasp.q
            if self.from_coarse_relax:
                # POSCAR
                crelax_path = os.path.join(self.base_path, "rlx-coarse")
                orig_poscar_path = os.path.join(crelax_path, "CONTCAR")
                final_poscar_path = os.path.join(self.calc_path, "POSCAR")
                shutil.copy(orig_poscar_path, final_poscar_path)
            else:
                # POSCAR
                orig_poscar_path = os.path.join(self.base_path, "POSCAR")
                final_poscar_path = os.path.join(self.calc_path, "POSCAR")
                shutil.copy(orig_poscar_path, final_poscar_path)
            structure = get_pmg_structure_from_poscar(final_poscar_path)

            # POTCAR
            potcar_path = os.path.join(self.calc_path, "POTCAR")
            make_potcar(structure, potcar_path)

            # INCAR
            incar_path = os.path.join(self.calc_path, "INCAR")
            make_incar(incar_path, mode=self.mode)

            # vasp.q
            vaspq_path = os.path.join(self.calc_path, "vasp.q")
            make_vaspq(vaspq_path, mode=self.mode, jobname=self.material_name)

        if self.to_submit:
            job_submitted = self.submit_job()
            if not job_submitted:
                self.to_rerun = False
                self.setup_calc()

    def check_calc(self):
        """
        Check if calculation has finished and reached required accuracy
        (No real automatic logging or fixing of VASP errors)

        Returns:
            relaxation_successful (bool): if True, relaxation completed successfully
        """
        stdout_path = os.path.join(self.calc_path, "stdout.txt")

        if os.path.exists(stdout_path):
            if not self.job_complete:
                logger.info(f"{self.mode.upper()} not finished")
                return False

            grep_call = f"tail -n{self.tail} {stdout_path}"
            grep_output = (
                subprocess.check_output(grep_call, shell=True).decode("utf-8").strip()
            )
            if "reached required accuracy" in grep_output:
                logger.info(
                    f"{self.mode.upper()} Calculation: reached required accuracy"
                )
                logger.debug(grep_output)
                return True
            else:
                archive_dirs = glob.glob(f"{self.calc_path}/archive*")
                if len(archive_dirs) >= 3:
                    logger.warning("Many archives exist, suggest force based relaxation")
                    if self.to_rerun:
                        self.setup_calc()
                    return True

                logger.warning(f"{self.mode.upper()} FAILED")
                logger.debug(grep_output)
                if self.to_rerun:
                    logger.info(f"Rerunning {self.calc_path}")
                    self.setup_calc()
                return False
        else:
            logger.info(f"{self.mode.upper()} not started")
            return False

    def check_volume_difference(self):
        """
        check relaxation runs for volume difference

        if abs(volume difference) is >= 5%, rerun relaxation
        only check for mode='rlx' as that's the structure for elastic analysis

        Returns:
            volume_converged (bool): if True, relaxation completed successfully
        """
        poscar_path = os.path.join(self.calc_path, "POSCAR")
        contcar_path = os.path.join(self.calc_path, "CONTCAR")
        try:
            p_structure, p_spacegroup = get_pmg_structure_from_poscar(
                poscar_path, return_sg=True
            )
            c_structure, c_spacegroup = get_pmg_structure_from_poscar(
                contcar_path, return_sg=True
            )
        except Exception as e:
            logger.info(f"  RLX CONTCAR doesn't exist or is empty: {e}")
            return False

        if p_spacegroup == c_spacegroup:
            logger.info(f"  Spacegroups match {p_spacegroup}=={c_spacegroup}")
        else:
            logger.warning(
                f"   Warning: spacegroups do not match {p_spacegroup} != {c_spacegroup}"
            )

        volume_diff = (c_structure.volume - p_structure.volume) / p_structure.volume
        if np.abs(volume_diff) >= 0.05:
            logger.warning(f"  NEED TO RE-RELAX: dV = {volume_diff:.4f}")
            volume_converged = False
            if self.to_rerun:
                self.setup_calc()
        else:
            logger.info(f"  RLX volume converged")
            logger.info(f"  dV = {volume_diff:.4f}")
            volume_converged = True
        return volume_converged

    @property
    def is_done(self):
        calc_done = self.check_calc()
        if calc_done:
            volume_converged = self.check_volume_difference()
            if volume_converged:
                return True
        else:
            return False
