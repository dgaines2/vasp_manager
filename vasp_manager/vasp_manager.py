# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import glob
import json
import logging
import os
from functools import cached_property
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
    Handles set up and execution of each CalculationManager
        (rlx-coarse, rlx, bulkmod, elastic)
    """

    def __init__(
        self,
        calculation_types,
        material_paths=None,
        to_rerun=True,
        to_submit=True,
        ignore_personal_errrors=True,
        tail=5,
        write_results=True,
        ncore=None,
        use_multiprocessing=False,
        calculation_manager_kwargs={},
    ):
        """
        Args:
            calculation_types (list[str]): list of calculation types
            material_paths (list[str] | str): list of material paths OR
                name of calculations dir
            to_rerun (bool): if True, rerun failed calculations
            to_submit (bool): if True, submit calculations
            ignore_personal_errors (bool): if True, ignore job submission errors
                if on personal computer
            tail (int): number of last lines to log in debugging if job failed
            write_results (bool): if True, dump results to results.json
            ncore (int): if ncore, use {ncore} for multiprocessing
                if None, defaults to minimum(number of materials, cpu_count)
            use_multiprocessing (bool): if True, use pool.map()
                Can be useful to set to false for debugging
                WARNING
                I've had issues with job submission using multiprocessing, so consider
                it an experimental feature for now
                \\TODO: use a multiprocessing queue manager to handle the jobs
                to ensure concurrency isn't an issue
            calculation_manager_kwargs (dict): contains subdictionaries for each
                calculation type. Eeach subdictorary can be filled with extra kwargs
                to pass to its associated CalculationManager during instantiation
        """
        self.calculation_types = calculation_types
        self.material_paths = material_paths
        self.to_rerun = to_rerun
        self.to_submit = to_submit
        self.ignore_personal_errors = ignore_personal_errrors
        self.tail = tail
        self.write_results = write_results
        self.ncore = (
            ncore
            if ncore is not None
            else int(np.min([len(self.material_paths), cpu_count()]))
        )
        self.use_multiprocessing = use_multiprocessing
        self.calculation_manager_kwargs = calculation_manager_kwargs

        self.calculation_managers = self._get_all_calculation_managers()
        self.results_path = os.path.join(self.base_path, "results.json")
        self.results = None

    @property
    def calculation_types(self):
        return self._calculation_types

    @calculation_types.setter
    def calculation_types(self, values):
        if not isinstance(values, list):
            raise TypeError("calculation_types must be a list")
        self._calculation_types = values

    @property
    def ncore(self):
        return self._ncore

    @ncore.setter
    def ncore(self, value):
        if not isinstance(value, int):
            raise Exception
        self._ncore = value

    @property
    def calculation_manager_kwargs(self):
        return self._calculation_manager_kwargs

    @calculation_manager_kwargs.setter
    def calculation_manager_kwargs(self, values):
        if not isinstance(values, dict):
            raise TypeError("calculation_manager_kwargs must be a dictionary")

        supported_kwargs = ["from_scratch", "strains"]
        for calc_type in self.calculation_types:
            if calc_type not in values:
                values[calc_type] = {}
            else:
                for kwarg in values[calc_type]:
                    if kwarg not in supported_kwargs:
                        raise ValueError(
                            f"kwarg={kwarg} is not supported for mode={calc_type}"
                        )
        self._calculation_manager_kwargs = values

    @property
    def material_paths(self):
        return self._material_paths

    @material_paths.setter
    def material_paths(self, values):
        """
        Sets paths for all materials

        Args:
            values (list[str] | str): list of material paths OR name
                of calculations dir

                if is list, use that list directly
                if is string, find folders inside of that directory named
                    {_material_paths}
        """
        match values:
            case str():
                self.base_path = values
                material_paths = [
                    d for d in glob.glob(os.path.join(values, "*")) if os.path.isdir(d)
                ]
            case list() | np.array():
                mat_path = values[0]
                self.base_path = os.path.dirname(mat_path)
                material_paths = values
            case _:
                raise TypeError(
                    "material_paths must be a directory name or a list of paths"
                )
        # Sort the paths by name
        self._material_paths = sorted(material_paths)

    def _get_material_name_from_path(self, material_path):
        material_name = os.path.basename(material_path)
        return material_name

    @cached_property
    def material_names(self):
        return [self._get_material_name_from_path(mpath) for mpath in self.material_paths]

    @property
    def results(self):
        return self._results

    @results.setter
    def results(self, value):
        if value is None:
            if os.path.exists(self.results_path):
                with open(self.results_path, "r") as fr:
                    self._results = json.load(fr)
            else:
                self._results = {mat_name: {} for mat_name in self.material_names}

    def _check_calc_by_result(self, material_name, calc_type):
        """
        Checks if job has been completed and analyzed

        Args:
            material_name (str): name of material to check
            calc_type (str): calculation type to check
        Returns:
            is_done (bool)
        """
        match calc_type:
            case "rlx-coarse" | "rlx":
                is_done = self.results[material_name][calc_type] == "done"
            case "static" | "bulkmod" | "elastic":
                is_done = self.results[material_name][calc_type] is not None
            case _:
                raise ValueError("Can't find mode {mode} in result")
        return is_done

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
                        tail=self.tail,
                        **self.calculation_manager_kwargs[calc_type],
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
                        tail=self.tail,
                        **self.calculation_manager_kwargs[calc_type],
                    )
                case "static":
                    if "rlx" not in self.calculation_types:
                        msg = "Cannot perform static calculation without mode='rlx' first"
                        raise Exception(msg)
                    manager = StaticCalculationManager(
                        material_path=material_path,
                        to_rerun=self.to_rerun,
                        to_submit=self.to_submit,
                        ignore_personal_errors=self.ignore_personal_errors,
                        tail=self.tail,
                        **self.calculation_manager_kwargs[calc_type],
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
                        **self.calculation_manager_kwargs[calc_type],
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
                        tail=self.tail,
                        **self.calculation_manager_kwargs[calc_type],
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
        for material_name, material_path in zip(self.material_names, self.material_paths):
            calc_managers[material_name] = self._get_calculation_managers(material_path)
        return calc_managers

    def _manage_calculations(self, material_name):
        """
        Runs vasp job workflow for a single material
        """
        for calc_manager in self.calculation_managers[material_name]:
            if calc_manager.mode in self.results[material_name].keys():
                calc_is_done = self._check_calc_by_result(
                    material_name, calc_manager.mode
                )
                if calc_is_done:
                    logger.info(
                        f"{material_name} -- {calc_manager.mode.upper()} Successful"
                    )
                    continue

            if not calc_manager.job_exists:
                logger.info(f"{material_name} Setting up" f" {calc_manager.mode.upper()}")
                calc_manager.setup_calc()
                match calc_manager.mode:
                    case "rlx-coarse" | "rlx":
                        break
                    case _:
                        pass

            if not calc_manager.is_done:
                match calc_manager.mode:
                    case "rlx-coarse" | "rlx" | "elastic":
                        # don't check further modes as they rely on rlx-coarse
                        # or rlx to be done
                        # break if elastic not done to avoid analysis
                        break
                    case _:
                        # go ahead and check the other modes as they are
                        # independent of each other
                        pass

            self.results[material_name][calc_manager.mode] = calc_manager.results

    def _manage_calculations_wrapper(self):
        if self.use_multiprocessing:
            with Pool(self.ncore) as pool:
                pool.map(self._manage_calculations, tqdm(self.material_names))
        else:
            for i, material_name in enumerate(self.material_names):
                print(f"{i+1}/{len(self.material_names)} -- {material_name}")
                self._manage_calculations(material_name)
                logger.info("")
                logger.info("")
                logger.info("")

    def run_calculations(self):
        """
        Runs vasp job workflow for all materials
        """
        self._manage_calculations_wrapper()

        json_str = json.dumps(self.results, indent=2, cls=NumpyEncoder)
        logger.debug(json_str)
        with open(self.results_path, "w+") as fw:
            fw.write(json_str)
        print(f"Dumped to {self.results_path}")
        return self.results

    def summary(self, as_string=True):
        """
        Create a string summary of all calculations

        Args:
            as_string (bool):
                if as_string, return string summary
                else, return dict summary
        """
        if not os.path.exists(self.results_path):
            raise ValueError(f"Can't find results in {self.results_path}")
        with open(self.results_path) as fr:
            results = json.load(fr)

        summary_dict = {}
        summary_dict["n_total"] = len(self.material_paths)
        for calc_type in self.calculation_types:
            summary_dict[calc_type] = {}
            summary_dict[calc_type]["n_finished"] = 0
            summary_dict[calc_type]["finished"] = []
            summary_dict[calc_type]["unfinished"] = []

            for material, mat_results in results.items():
                # need to account for case key doesn't yet exist
                if calc_type not in mat_results:
                    summary_dict[calc_type]["unfinished"].append(material)
                else:
                    is_done = self._check_calc_by_result(material, calc_type)
                    if is_done:
                        summary_dict[calc_type]["n_finished"] += 1
                        summary_dict[calc_type]["finished"].append(material)
                    else:
                        summary_dict[calc_type]["unfinished"].append(material)

        if as_string:
            n_materials = summary_dict["n_total"]
            summary_str = ""
            summary_str += f"Total Materials = {n_materials}\n"
            summary_str += "-" * 30 + "\n"
            for calc_type in self.calculation_types:
                name = calc_type.upper()
                n_finished = summary_dict[calc_type]["n_finished"]
                summary_str += f"{name: <12} {n_finished}/{n_materials} completed\n"
            return summary_str
        else:
            return summary_dict
