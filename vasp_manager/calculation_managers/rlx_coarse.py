# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import glob
import logging
import os
import shutil
import subprocess

import pymatgen as pmg

from ..vasp_utils import make_archive, make_incar, make_potcar, make_vaspq
from .base import BaseCalculationManager

logger = logging.getLogger(__name__)


class RlxCoarseCalculationManager(BaseCalculationManager):
    def __init__(
        self,
        base_path,
        to_rerun,
        to_submit,
        ignore_personal_errors=True,
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
        self.tail = tail

    @property
    def mode(self):
        return "rlx-coarse"

    def setup_calc(self):
        """
        Set up a coarse relaxation
        """
        # POSCAR, POTCAR, INCAR, vasp.q
        if not os.path.exists(self.calc_path):
            os.mkdir(self.calc_path)

        # POSCAR
        if self.to_rerun:
            archive_made = make_archive(self.calc_path, mode=self.mode)
            if not archive_made:
                # set rerun to False to not make an achive and instead
                # continue to make the input files
                self.to_rerun = False
                self.setup_calc()
        else:
            orig_poscar_path = os.path.join(self.base_path, "POSCAR")
            final_poscar_path = os.path.join(self.calc_path, "POSCAR")
            shutil.copy(orig_poscar_path, final_poscar_path)

            # POTCAR
            structure = pmg.core.Structure.from_file(final_poscar_path)
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
            # job status returns True if sucessfully submitted, else False
            if not job_submitted:
                # if job didn't submit, try rerunning setup
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

    @property
    def is_done(self):
        return self.check_calc()
