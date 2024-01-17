# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import json
import logging
from functools import cached_property
from multiprocessing import Pool
from pathlib import Path

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
        use_multiprocessing=False,
        ncore=None,
        calculation_manager_kwargs={},
        max_reruns=3,
        magmom_per_atom_cutoff=0.0,
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
            tail (int): number of last lines from stdout.txt to log in debugging
                if job failed
            use_multiprocessing (bool): if True, use pool.map()
            ncore (int): if ncore, use {ncore} for multiprocessing
                if None, defaults to minimum(number of materials, 4)
            calculation_manager_kwargs (dict): contains subdictionaries for each
                calculation type. Each subdictorary can be filled with extra kwargs
                to pass to its associated CalculationManager during instantiation
            max_reruns (int): the maximum number of times a rlx-coarse or rlx
                calculation can run before refusing to continue
                Note: other modes don't make archives, so they are not affected
                by this
            magmom_per_atom_cutoff (float): calculations that result in
                magmom_per_atom less than this parameter will be automatically
                rerun without spin-polarization
        """
        self.calculation_types = calculation_types
        self.material_paths = material_paths
        self.to_rerun = to_rerun
        self.to_submit = to_submit
        self.ignore_personal_errors = ignore_personal_errrors
        self.tail = tail
        self.use_multiprocessing = use_multiprocessing
        self.ncore = ncore
        self.calculation_manager_kwargs = calculation_manager_kwargs
        self.max_reruns = max_reruns
        self.magmom_per_atom_cutoff = magmom_per_atom_cutoff

        self.calculation_managers = self._get_all_calculation_managers()
        # self.base_path is set in material_paths.setter
        self.results_path = self.base_path / "results.json"
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
        if value is None:
            value = int(np.min([len(self.material_paths), 4]))
            if self.use_multiprocessing:
                print(
                    "WARNING: setting default ncore for multiprocessing to "
                    + f"{value}\n"
                    + "    We strongly recommend you set ncore to the number of "
                    + "available cores."
                )
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
                material_paths = [d for d in Path(values).glob("*") if d.is_dir()]
            case list() | np.array():
                values = [Path(p) for p in values]
                base_path = values[0].parent
                for mat_path in values:
                    if not mat_path.parent == base_path:
                        raise Exception(
                            "All material paths must be in the same directory"
                        )
                self.base_path = base_path
                material_paths = values
            case _:
                raise TypeError(
                    "material_paths must be a directory name or a list of paths"
                )
        # Sort the paths by name
        self._material_paths = sorted(material_paths)

    def _get_material_name_from_path(self, material_path):
        return material_path.name

    @cached_property
    def material_names(self):
        return [self._get_material_name_from_path(mpath) for mpath in self.material_paths]

    @property
    def results(self):
        return self._results

    @results.setter
    def results(self, value):
        if value is None:
            if self.results_path.exists():
                with open(self.results_path, "r") as fr:
                    self._results = json.load(fr)
                for mat_name in self.material_names:
                    if mat_name not in self._results:
                        self._results[mat_name] = {}
                self._results = dict(sorted(self._results.items()))
            else:
                self._results = {mat_name: {} for mat_name in self.material_names}

    def _check_calc_by_result(self, material_name, calc_type):
        """
        Checks if job has been completed and analyzed

        Args:
            material_name (str): name of material to check
            calc_type (str): calculation type to check
        Returns:
            (is_done (bool), is_stopped (bool))
        """
        match calc_type:
            case "rlx-coarse" | "rlx" | "static" | "bulkmod" | "elastic":
                is_done = self.results[material_name][calc_type] not in [
                    None,
                    "STOPPED",
                    "not finished",
                ]
                is_stopped = self.results[material_name][calc_type] == "STOPPED"
            case _:
                raise ValueError("Can't find mode {mode} in result")
        return (is_done, is_stopped)

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
                        max_reruns=self.max_reruns,
                        **self.calculation_manager_kwargs[calc_type],
                    )
                case "rlx":
                    from_coarse_relax = "rlx-coarse" in self.calculation_types
                    manager = RlxCalculationManager(
                        material_path=material_path,
                        to_rerun=self.to_rerun,
                        to_submit=self.to_submit,
                        ignore_personal_errors=self.ignore_personal_errors,
                        from_coarse_relax=from_coarse_relax,
                        tail=self.tail,
                        max_reruns=self.max_reruns,
                        magmom_per_atom_cutoff=self.magmom_per_atom_cutoff,
                        **self.calculation_manager_kwargs[calc_type],
                    )
                case "static":
                    from_relax = "rlx" in self.calculation_types
                    manager = StaticCalculationManager(
                        material_path=material_path,
                        to_rerun=self.to_rerun,
                        to_submit=self.to_submit,
                        ignore_personal_errors=self.ignore_personal_errors,
                        from_relax=from_relax,
                        tail=self.tail,
                        **self.calculation_manager_kwargs[calc_type],
                    )
                case "bulkmod":
                    from_relax = "rlx" in self.calculation_types
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
                            "Cannot perform elastic calculation without mode='rlx'"
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
        material_results = {}
        for calc_manager in self.calculation_managers[material_name]:
            if calc_manager.stopped:
                logger.info(f"{material_name} -- STOPPED")
                break

            if calc_manager.mode in self.results[material_name].keys():
                calc_is_done, calc_is_stopped = self._check_calc_by_result(
                    material_name, calc_manager.mode
                )
                if calc_is_done and not calc_manager.from_scratch:
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

            calc_is_done = calc_manager.is_done
            calc_is_stopped = calc_manager.stopped
            if calc_is_stopped:
                logger.info(f"{material_name} -- STOPPED")
            else:
                if not calc_is_done:
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

            material_results[calc_manager.mode] = calc_manager.results
        return (material_name, material_results)

    def _manage_calculations_wrapper(self):
        if self.use_multiprocessing:
            with Pool(self.ncore) as pool:
                results = pool.map(
                    self._manage_calculations, tqdm(self.material_names), 1
                )
        else:
            results = []
            for i, material_name in enumerate(self.material_names):
                print(f"{i+1}/{len(self.material_names)} -- {material_name}", flush=True)
                results.append(self._manage_calculations(material_name))
                print()
        return results

    def run_calculations(self):
        """
        Runs vasp job workflow for all materials
        """
        results = self._manage_calculations_wrapper()
        for material_name, material_result in results:
            self.results[material_name].update(material_result)

        json_str = json.dumps(self.results, indent=2, cls=NumpyEncoder)
        logger.debug(json_str)
        with open(self.results_path, "w+") as fw:
            fw.write(json_str)
        print(f"Dumped to {self.results_path}")
        return self.results

    def summary(self, as_string=True, print_unfinished=False, print_stopped=True):
        """
        Create a string summary of all calculations

        Args:
            as_string (bool):
                if as_string, return string summary
                else, return dict summary
        """
        if not self.results_path.exists():
            raise ValueError(f"Can't find results in {self.results_path}")
        with open(self.results_path) as fr:
            results = json.load(fr)

        summary_dict = {}
        summary_dict["n_total"] = len(results)
        for calc_type in self.calculation_types:
            summary_dict[calc_type] = {}
            summary_dict[calc_type]["n_finished"] = 0
            summary_dict[calc_type]["finished"] = []
            summary_dict[calc_type]["unfinished"] = []
            summary_dict[calc_type]["stopped"] = []

            for material, mat_results in results.items():
                # need to account for case key doesn't yet exist
                if calc_type not in mat_results:
                    summary_dict[calc_type]["unfinished"].append(material)
                else:
                    is_done, is_stopped = self._check_calc_by_result(material, calc_type)
                    if is_done:
                        summary_dict[calc_type]["n_finished"] += 1
                        summary_dict[calc_type]["finished"].append(material)
                    else:
                        if is_stopped:
                            summary_dict[calc_type]["stopped"].append(material)
                        summary_dict[calc_type]["unfinished"].append(material)

        if as_string:
            n_materials = summary_dict["n_total"]
            summary_str = ""
            summary_str += f"Total Materials = {n_materials}\n"
            summary_str += "-" * 30 + "\n"
            for calc_type in self.calculation_types:
                name = calc_type.upper()
                n_finished = summary_dict[calc_type]["n_finished"]
                summary_str += f"{name: <12}{n_finished}/{n_materials} completed\n"
                if print_unfinished:
                    unfinished = summary_dict[calc_type]["unfinished"]
                    if len(unfinished) != 0:
                        summary_str += (
                            " " * 12 + f"{n_materials - n_finished} not completed\n"
                        )
                        summary_str += (
                            f"Unfinished {calc_type.upper()}: " + f"{unfinished}\n"
                        )
                if print_stopped:
                    stopped = summary_dict[calc_type]["stopped"]
                    if len(stopped) != 0:
                        summary_str += f"Stopped {calc_type.upper()}: " + f"{stopped}\n"
            return summary_str
        else:
            return summary_dict
