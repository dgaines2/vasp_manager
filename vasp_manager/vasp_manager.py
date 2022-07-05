# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import glob
import json
import logging
import os
from multiprocessing import Pool, cpu_count

import numpy as np
from tqdm import tqdm

from vasp_manager.calculation_manager import (
    BulkmodCalculationManager,
    ElasticCalculationManager,
    RlxCalculationManager,
    RlxCoarseCalculationManager,
    StaticCalculationManager,
)
from vasp_manager.utils import NumpyEncoder

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
        from_scratch=False,
        tail=5,
        write_results=True,
        ncore=None,
    ):
        """
        Args:
            calculation_types (list[str]): list of calculation types
            material_paths (list[str] | str): list of material paths OR name of calculations dir
            to_rerun (bool): if True, rerun failed calculations
            to_submit (bool): if True, submit calculations
            ignore_personal_errors (bool): if True, ignore job submission errors
                if on personal computer
            from_scratch (bool): if True, remove the first calculation's folder
                DANGEROUS
            tail (int): number of last lines to log in debugging if job failed
            write_results (bool): if True, dump results to results.json
            ncore (int): if ncore, use {ncore} for multiprocessing
                defaults to cpu_count()
        """
        self.calculation_types = calculation_types
        self.to_rerun = to_rerun
        self.to_submit = to_submit
        self.ignore_personal_errors = ignore_personal_errrors
        self.from_scratch = from_scratch
        self.tail = tail
        self.write_results = write_results
        self._ncore = ncore

        self.material_paths = self._get_material_paths(material_paths)
        self.calculation_managers = self._get_all_calculation_managers()

    @property
    def ncore(self):
        if self._ncore is None:
            self.ncore = cpu_count()
        return self._ncore

    @ncore.setter
    def ncore(self, value):
        if not isinstance(value, int):
            raise Exception
        self._ncore = value

    def _get_material_paths(self, _material_paths):
        """
        Gets paths for all materials

        Args:
            _material_paths (list[str] | str): list of material paths OR name of calculations dir
                if is list, use that list directly
                if is string, find folders inside of that directory named {_material_paths}
        """
        match _material_paths:
            case str():
                self.base_path = _material_paths
                material_paths = [
                    d
                    for d in glob.glob(os.path.join(_material_paths, "*"))
                    if os.path.isdir(d)
                ]
            case list() | np.array():
                mat_path = _material_paths[0]
                self.base_path = os.path.dirname(mat_path)
                material_paths = _material_paths
            case _:
                raise TypeError(
                    "material_paths must be a directory name or a list of paths"
                )
        # Sort the paths by name
        material_paths = sorted(material_paths)
        return material_paths

    def _get_material_name_from_path(self, material_path):
        material_name = os.path.basename(material_path)
        return material_name

    def _get_calculation_managers(self, material_path):
        """
        Gets calculation managers for a single material
        """
        calc_managers = []
        for calc_type in self.calculation_types:
            match calc_type:
                case "rlx-coarse":
                    manager = RlxCoarseCalculationManager(
                        material_path=material_path,
                        to_rerun=self.to_rerun,
                        to_submit=self.to_submit,
                        ignore_personal_errors=self.ignore_personal_errors,
                        from_scratch=self.from_scratch,
                        tail=self.tail,
                    )
                case "rlx":
                    if "rlx-coarse" in self.calculation_types:
                        from_coarse_relax = True
                    else:
                        from_coarse_relax = False
                    manager = RlxCalculationManager(
                        material_path=material_path,
                        to_rerun=self.to_rerun,
                        to_submit=self.to_submit,
                        ignore_personal_errors=self.ignore_personal_errors,
                        from_coarse_relax=from_coarse_relax,
                        from_scratch=self.from_scratch,
                        tail=self.tail,
                    )
                case "static":
                    if "rlx" not in self.calculation_types:
                        msg = (
                            "Cannot perform static calculation without mode='rlx' first"
                        )
                        raise Exception(msg)
                    manager = StaticCalculationManager(
                        material_path=material_path,
                        to_rerun=self.to_rerun,
                        to_submit=self.to_submit,
                        ignore_personal_errors=self.ignore_personal_errors,
                        from_scratch=self.from_scratch,
                        tail=self.tail,
                    )
                case "bulkmod" | "bulkmod_standalone":
                    if calc_type == "bulkmod":
                        if "rlx" in self.calculation_types:
                            from_relax = True
                        else:
                            from_relax = False
                    else:
                        from_relax = False
                        if len(self.calculation_types) > 1:
                            msg = (
                                "bulkmod_standalone must be run alone -- remove other"
                                " calculation types to fix"
                            )
                            raise Exception()
                    manager = BulkmodCalculationManager(
                        material_path=material_path,
                        to_rerun=self.to_rerun,
                        to_submit=self.to_submit,
                        ignore_personal_errors=self.ignore_personal_errors,
                        from_relax=from_relax,
                        from_scratch=self.from_scratch,
                    )
                case "elastic":
                    if "rlx" not in self.calculation_types:
                        msg = (
                            "Cannot perform elastic calculation without mode='rlx-fine'"
                            " first"
                        )
                        raise Exception(msg)
                    manager = ElasticCalculationManager(
                        material_path=material_path,
                        to_rerun=self.to_rerun,
                        to_submit=self.to_submit,
                        ignore_personal_errors=self.ignore_personal_errors,
                        from_scratch=self.from_scratch,
                        tail=self.tail,
                    )
                case _:
                    raise Exception(f"Calc type {calc_type} not supported")

            calc_managers.append(manager)
        return calc_managers

    def _get_all_calculation_managers(self):
        """
        Gets calculation managers for all materials
        """
        calc_managers = {}
        for material_path in self.material_paths:
            material_name = self._get_material_name_from_path(material_path)
            calc_managers[material_name] = self._get_calculation_managers(material_path)
        return calc_managers

    def _manage_calculations(self, material_name):
        """
        Runs vasp job workflow for a single material
        """
        results = {}
        for calc_manager in self.calculation_managers[material_name]:
            if not calc_manager.job_exists:
                logger.info(
                    f"{calc_manager.material_name} Setting up"
                    f" {calc_manager.mode.upper()}"
                )
                calc_manager.setup_calc()
                match calc_manager.mode:
                    case "rlx-coarse" | "rlx":
                        break
                    case _:
                        pass

            if not calc_manager.is_done:
                match calc_manager.mode:
                    case "rlx-coarse" | "rlx":
                        # don't check further modes as they rely on rlx-coarse
                        # or rlx to be done
                        break
                    case _:
                        # go ahead and check the other modes as they are
                        # independent of each other
                        continue

            results[calc_manager.mode] = calc_manager.results
        return results

    def _manage_calculations_wrapper(self):
        material_names = [
            self._get_material_name_from_path(p) for p in self.material_paths
        ]
        with Pool(self.ncore) as pool:
            results = pool.map(self._manage_calculations, tqdm(material_names))

        results_dict = {}
        for material_name, result in zip(material_names, results):
            results_dict[material_name] = result
        return results_dict

    def run_calculations(self):
        """
        Runs vasp job workflow for all materials
        """
        all_results = self._manage_calculations_wrapper()

        json_str = json.dumps(all_results, indent=2, cls=NumpyEncoder)
        logger.info(json_str)
        if self.write_results:
            results_path = os.path.join(self.base_path, "results.json")
            with open(results_path, "w+") as fw:
                fw.write(json_str)
            logger.info("Dumping to results.json")

        self.results = all_results
        return all_results

    def summary(self):
        """
        Create a string summary of all calculations
        """
        n_materials = len(self.material_paths)
        summary_dict = {}
        for calc_type in self.calculation_types:
            summary_dict[calc_type] = {}
            summary_dict[calc_type]["n_finished"] = 0
            summary_dict[calc_type]["finished"] = []
            summary_dict[calc_type]["unfinished"] = []

            for material, mat_results in self.results.items():
                # need to account for case key doesn't yet exist
                if calc_type in mat_results:
                    match calc_type:
                        case "rlx-coarse" | "rlx" | "static":
                            case_condition = mat_results[calc_type] == "done"
                        case "bulkmod" | "bulkmod_standalone" | "elastic":
                            case_condition = mat_results[calc_type] is not None
                        case _:
                            raise ValueError(f"Unexpected calc_type: {calc_type}")
                    if case_condition:
                        summary_dict[calc_type]["n_finished"] += 1
                        summary_dict[calc_type]["finished"].append(material)
                    else:
                        summary_dict[calc_type]["unfinished"].append(material)
                else:
                    summary_dict[calc_type]["unfinished"].append(material)

        summary_str = ""
        summary_str += f"Total Materials = {n_materials}\n"
        summary_str += "-" * 30 + "\n"
        for calc_type in summary_dict:
            name = calc_type.upper()
            n_finished = summary_dict[calc_type]["n_finished"]
            summary_str += f"{name: <12} {n_finished}/{n_materials} completed\n"
        return summary_str
