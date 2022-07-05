# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import glob
import logging
import os
from functools import cached_property

import numpy as np

from vasp_manager.calculation_managers.base import BaseCalculationManager
from vasp_manager.utils import get_pmg_structure_from_poscar, ptail
from vasp_manager.vasp_input_creator import VaspInputCreator

logger = logging.getLogger(__name__)


class RlxCalculationManager(BaseCalculationManager):
    """
    Runs relaxation job workflow for a single material
    """

    def __init__(
        self,
        material_path,
        to_rerun,
        to_submit,
        ignore_personal_errors=True,
        from_coarse_relax=True,
        from_scratch=False,
        tail=5,
    ):
        """
        For material_path, to_rerun, to_submit, ignore_personal_errors, and from_scratch,
        see BaseCalculationManager

        Args:
            from_coarse_relax (bool): if True, use CONTCAR from coarse relax
            tail (int): number of last lines to log in debugging if job failed
        """
        self.from_coarse_relax = from_coarse_relax
        self.tail = tail
        super().__init__(
            material_path=material_path,
            to_rerun=to_rerun,
            to_submit=to_submit,
            ignore_personal_errors=ignore_personal_errors,
            from_scratch=from_scratch,
        )
        self._is_done = None
        self._results = None

    @cached_property
    def mode(self):
        return "rlx"

    @cached_property
    def poscar_source_path(self):
        if self.from_coarse_relax:
            poscar_source_path = os.path.join(
                self.material_path, "rlx-coarse", "CONTCAR"
            )
        else:
            poscar_source_path = os.path.join(self.material_path, "POSCAR")
        return poscar_source_path

    def setup_calc(self):
        """
        Sets up a fine relaxation
        """
        vasp_input_creator = VaspInputCreator(
            self.calc_path,
            mode=self.mode,
            poscar_source_path=self.poscar_source_path,
            name=self.material_name,
        )
        if self.to_rerun:
            archive_made = vasp_input_creator.make_archive()
            if not archive_made:
                # set rerun to not make an achive and instead
                # continue to make the input files
                self.to_rerun = False
                self.setup_calc()
                return
        else:
            vasp_input_creator.create()

        if self.to_submit:
            job_submitted = self.submit_job()
            if not job_submitted:
                self.to_rerun = False
                self.setup_calc()

    def check_calc(self):
        """
        Checks if calculation has finished and reached required accuracy

        Returns:
            relaxation_successful (bool): if True, relaxation completed successfully
        """
        if not self.job_complete:
            logger.info(f"{self.mode.upper()} job not finished")
            return False

        stdout_path = os.path.join(self.calc_path, "stdout.txt")
        if os.path.exists(stdout_path):
            if not self.job_complete:
                logger.info(f"{self.mode.upper()} not finished")
                return False

            tail_output = ptail(stdout_path, n_tail=self.tail, as_string=True)
            if "reached required accuracy" in tail_output:
                logger.info(
                    f"{self.mode.upper()} Calculation: reached required accuracy"
                )
                logger.debug(tail_output)
                return True
            else:
                archive_dirs = glob.glob(os.path.join(self.calc_path, "archive*"))
                if len(archive_dirs) >= 3:
                    logger.warning("Many archives exist, suggest force based relaxation")
                    if self.to_rerun:
                        self.setup_calc()
                    return True

                logger.warning(f"{self.mode.upper()} FAILED")
                logger.debug(tail_output)
                if self.to_rerun:
                    logger.info(f"Rerunning {self.calc_path}")
                    self.setup_calc()
                return False
        else:
            logger.info(f"{self.mode.upper()} not started")
            return False

    def check_volume_difference(self):
        """
        Checks relaxation runs for volume difference

        if abs(volume difference) is >= 5%, reruns relaxation
        only checks for mode='rlx' as that's the structure for further calculation

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
            logger.error(f"  RLX CONTCAR doesn't exist or is empty: {e}")
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
        if self._is_done is None:
            self._is_done = False
            calc_done = self.check_calc()
            if calc_done:
                volume_converged = self.check_volume_difference()
                if volume_converged:
                    self._is_done = True
        return self._is_done

    @property
    def results(self):
        if self._results is None:
            if self.is_done:
                self._results = "done"
            else:
                self._results = "not finished"
        return self._results
