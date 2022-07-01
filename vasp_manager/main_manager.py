# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import glob
import logging
import os

from vasp_manager.calculation_managers import (
    BulkmodCalculationManager,
    ElasticCalculationManager,
    RlxCalculationManager,
    RlxCoarseCalculationManager,
)

logger = logging.getLogger(__name__)


class VaspManager:
    """
    Handles set up and execution of each CalculationManager (rlx-coarse, rlx, bulkmod, elastic)
    """

    def __init__(
        self,
        calculation_types,
        material_paths=None,
        to_rerun=True,
        to_submit=True,
        ignore_personal_errrors=True,
        tail=5,
        from_scratch=False,
    ):
        """
        Args:
            calculation_types (list):
            material_paths (list):
            to_rerun (bool):
            to_submit (bool):
            ignore_personal_errors (bool):
            tail (int):
            from_scratch (bool): if True, remove the first calculation's folder
                DANGEROUS
        """
        self.calculation_types = calculation_types
        self.to_rerun = to_rerun
        self.to_submit = to_submit
        self.tail = tail
        self.ignore_personal_errors = ignore_personal_errrors
        self.from_scratch = from_scratch
        self.material_paths = (
            self._get_material_paths() if material_paths is None else material_paths
        )
        self.calculation_managers = self._get_all_calculation_managers()

    def _get_material_paths(self):
        """
        Defaults to all paths under calculations/*
        """
        material_paths = [d for d in glob.glob("calculations/*") if os.path.isdir(d)]
        # Sort the paths by name
        material_paths = sorted(material_paths, key=lambda d: int(d.split("/")[1]))
        return material_paths

    def _get_material_name_from_path(self, material_path):
        material_name = material_path.split("/")[-1]
        return material_name

    def _get_calculation_managers(self, material_path):
        """
        Get calculation managers for a single material
        """
        calc_managers = []
        for calc_type in self.calculation_types:
            match calc_type:
                case "rlx-coarse":
                    manager = RlxCoarseCalculationManager(
                        base_path=material_path,
                        to_rerun=self.to_rerun,
                        to_submit=self.to_submit,
                        ignore_personal_errors=self.ignore_personal_errors,
                        tail=self.tail,
                        from_scratch=self.from_scratch,
                    )
                case "rlx-fine":
                    if "rlx-coarse" in self.calculation_types:
                        from_coarse_relax = True
                    else:
                        from_coarse_relax = False
                    manager = RlxCalculationManager(
                        base_path=material_path,
                        to_rerun=self.to_rerun,
                        to_submit=self.to_submit,
                        ignore_personal_errors=self.ignore_personal_errors,
                        from_coarse_relax=from_coarse_relax,
                        tail=self.tail,
                        from_scratch=self.from_scratch,
                    )
                case "bulkmod":
                    if "rlx-fine" in self.calculation_types:
                        from_relax = True
                    else:
                        from_relax = False
                        msg = (
                            "Running bulk modulus calculation without previous"
                            " relaxation"
                        )
                        msg += (
                            "\n\t starting structure must be fairly close to equilibrium"
                            " volume!"
                        )
                        logger.warning(msg)
                    manager = BulkmodCalculationManager(
                        base_path=material_path,
                        to_rerun=self.to_rerun,
                        to_submit=self.to_submit,
                        ignore_personal_errors=self.ignore_personal_errors,
                        from_relax=from_relax,
                        from_scratch=self.from_scratch,
                    )
                case "elastic":
                    if "rlx-fine" not in self.calculation_types:
                        msg = (
                            "Cannot perform elastic calculation without mode='rlx-fine'"
                            " first"
                        )
                        raise Exception(msg)
                    manager = ElasticCalculationManager(
                        base_path=material_path,
                        to_rerun=self.to_rerun,
                        to_submit=self.to_submit,
                        ignore_personal_errors=self.ignore_personal_errors,
                        tail=self.tail,
                        from_scratch=self.from_scratch,
                    )
                case _:
                    raise Exception(f"Calc type {calc_type} not supported")

            calc_managers.append(manager)
        return calc_managers

    def _get_all_calculation_managers(self):
        """
        Get calculation managers for all materials
        """
        calc_managers = {}
        for material_path in self.material_paths:
            material_name = self._get_material_name_from_path(material_path)
            calc_managers[material_name] = self._get_calculation_managers(material_path)
        return calc_managers

    def _manage_calculations(self, material_name):
        """
        Run vasp job workflow for a single material

        Faulty logic here..

        """
        for calc_manager in self.calculation_managers[material_name]:
            if not calc_manager.job_exists:
                logger.info(
                    f"{calc_manager.material_name} Setting up"
                    f" {calc_manager.mode.upper()}"
                )
                calc_manager.setup_calc()
                break

            if not calc_manager.is_done:
                break

            logger.info(f"Calculation {calc_manager.mode.upper()} successful")

    def run_calculations(self):
        """
        Run vasp job workflow for all materials
        """
        for material_path in self.material_paths:
            material_name = self._get_material_name_from_path(material_path)
            print(material_name)
            self._manage_calculations(material_name)
            print("\n")
