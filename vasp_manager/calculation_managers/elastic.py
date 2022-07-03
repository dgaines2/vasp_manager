# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import logging
import os
import subprocess

from vasp_manager.calculation_managers.base import BaseCalculationManager
from vasp_manager.elastic_analysis import analyze_elastic_file, make_elastic_constants
from vasp_manager.vasp_input_creator import VaspInputCreator

logger = logging.getLogger(__name__)


class ElasticCalculationManager(BaseCalculationManager):
    """
    Runs elastic deformation job workflow for a single material
    """

    def __init__(
        self,
        base_path,
        to_rerun,
        to_submit,
        ignore_personal_errors=True,
        from_scratch=False,
        tail=5,
    ):
        """
        For base_path, to_rerun, to_submit, ignore_personal_errors, and from_scratch,
        see BaseCalculationManager

        Args:
            tail (int): number of last lines to log in debugging if job failed
        """
        self.tail = tail
        super().__init__(
            base_path=base_path,
            to_rerun=to_rerun,
            to_submit=to_submit,
            ignore_personal_errors=ignore_personal_errors,
            from_scratch=from_scratch,
        )
        self._results = None

    @property
    def mode(self):
        return "elastic"

    @property
    def poscar_source_path(self):
        poscar_source_path = os.path.join(self.base_path, "rlx", "CONTCAR")
        return poscar_source_path

    def setup_calc(self, increase_nodes=False):
        """
        Runs elastic constants routine through VASP

        By default, requires relaxation (as the elastic constants routine needs
            the cell to be nearly at equilibrium)
        """
        vasp_input_creator = VaspInputCreator(
            self.calc_path,
            mode=self.mode,
            poscar_source_path=self.poscar_source_path,
            name=self.material_name,
            increase_nodes=increase_nodes,
        )
        vasp_input_creator.create()

        if self.to_submit:
            job_submitted = self.submit_job()
            # job status returns True if sucessfully submitted, else False
            if not job_submitted:
                self.setup_calc()

    def check_calc(self):
        """
        Checks result of elastic calculation

        Returns:
            elastic_successful (bool): if True, elastic calculation completed successfully
        """
        if not self.job_complete:
            logger.info(f"{self.mode.upper()} job not finished")
            return False

        stdout_path = os.path.join(self.calc_path, "stdout.txt")
        if os.path.exists(stdout_path):
            grep_call = f"grep 'Total' {stdout_path}"
            grep_output = (
                subprocess.check_output(grep_call, shell=True)
                .decode("utf-8")
                .splitlines()
            )
            last_grep_line = grep_output[-1].strip().split()
            # last grep line looks something like 'Total: 36/ 36'
            finished_deformations = int(last_grep_line[-2].replace("/", ""))
            total_deformations = int(last_grep_line[-1])
            logger.debug(last_grep_line)
            if finished_deformations == total_deformations:
                logger.info(f"{self.mode.upper()} Calculation: Success")
                return True
            else:
                tail_call = f"tail -n{self.tail} {stdout_path}"
                tail_output = (
                    subprocess.check_output(tail_call, shell=True)
                    .decode("utf-8")
                    .strip()
                )
                logger.info(tail_output)
                logger.info(f"{self.mode.upper()} Calculation: FAILED")
                if self.to_rerun:
                    # increase nodes as its likely the calculation failed
                    # setup_elastic(elastic_path, submit=submit, increase_nodes=True)
                    self.setup_calc(increase_nodes=True)
                return False
        else:
            # shouldn't get here unless function was called with submit=False
            logger.info(f"{self.mode.upper()} Calculation: No stdout.txt available")
            if self.to_rerun:
                # setup_elastic(elastic_path, submit=submit, increase_nodes=False)
                self.setup_calc(increase_nodes=False)
            return False

    @property
    def is_done(self):
        return self.check_calc()

    def _analyze_elastic(self):
        """
        Gets results from elastic calculation
        """
        elastic_file = os.path.join(self.calc_path, "elastic_constants.txt")
        if not os.path.exists(elastic_file):
            outcar_file = os.path.join(self.calc_path, "OUTCAR")
            make_elastic_constants(outcar_file)
        results = analyze_elastic_file(elastic_file)
        return results

    @property
    def results(self):
        if self._results is None:
            self._results = self._analyze_elastic()
        return self._results
