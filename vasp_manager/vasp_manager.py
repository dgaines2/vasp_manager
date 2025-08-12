# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

from __future__ import annotations

import json
import logging
from collections.abc import Callable
from functools import cached_property
from importlib.metadata import version
from multiprocessing import Pool
from pathlib import Path
from typing import TYPE_CHECKING, Any

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

if TYPE_CHECKING:
    from vasp_manager.types import (
        CalculationManager,
        CalculationType,
        Filepaths,
        WorkingDirectory,
    )

logger = logging.getLogger(__name__)

ASCII_LOGO = r"""
 __      __             __  __
 \ \    / /            |  \/  |
  \ \  / /_ _ ___ _ __ | \  / | __ _ _ __   __ _  __ _  ___ _ __
   \ \/ / _` / __| '_ \| |\/| |/ _` | '_ \ / _` |/ _` |/ _ \ '__|
    \  / (_| \__ \ |_) | |  | | (_| | | | | (_| | (_| |  __/ |
     \/ \__,_|___/ .__/|_|  |_|\__,_|_| |_|\__,_|\__, |\___|_|
                 | |                              __/ |
                 |_|                             |___/ v{}
""".format(
    version("vasp_manager")
)


class VaspManager:
    """
    Handles set up and execution of each CalculationManager
        (rlx-coarse, rlx, static, bulkmod, elastic)
    """

    def __init__(
        self,
        calculation_types: list[CalculationType],
        material_dirs: Filepaths | WorkingDirectory,
        to_rerun: bool = True,
        to_submit: bool = True,
        ignore_personal_errors: bool = True,
        tail: int = 5,
        use_multiprocessing: bool = False,
        ncore: None | int = None,
        calculation_manager_kwargs: dict | None = None,
        max_reruns: int = 3,
        magmom_per_atom_cutoff: float = 0.0,
        sort_by: Callable[[str], str] = str,
    ) -> None:
        """
        Args:
            calculation_types: list of calculation types
            material_dirs: list of material directory paths or name of
                calculations directory
            to_rerun: if True, rerun failed calculations
            to_submit: if True, submit calculations
            ignore_personal_errors: if True, ignore job submission errors if on
                personal computer
            tail: number of last lines from stdout.txt to log in debugging
                if job failed
            use_multiprocessing: if True, use pool.map()
            ncore: if ncore, use {ncore} processes for multiprocessing
                if None, defaults to minimum(number of materials, 4)
            calculation_manager_kwargs: contains subdictionaries for each
                calculation type. Each subdictorary can be filled with extra kwargs
                to pass to its associated CalculationManager during instantiation
            max_reruns: the maximum number of times a rlx-coarse or rlx
                calculation can run before refusing to continue
                Note: other modes don't make archives, so they are not affected
                by this
            magmom_per_atom_cutoff: calculations that result in
                magmom_per_atom less than this parameter will be automatically
                rerun without spin-polarization
            sort_by: function to sort the keys of the result dictionary
        """
        print(ASCII_LOGO)
        self.sort_by = sort_by
        self.calculation_types = calculation_types
        self.material_dirs = material_dirs
        self.to_rerun = to_rerun
        self.to_submit = to_submit
        self.ignore_personal_errors = ignore_personal_errors
        self.tail = tail
        self.use_multiprocessing = use_multiprocessing
        self.ncore = ncore
        self.calculation_manager_kwargs = (
            calculation_manager_kwargs if calculation_manager_kwargs else {}
        )
        self.max_reruns = max_reruns
        self.magmom_per_atom_cutoff = magmom_per_atom_cutoff

        self.calculation_managers = self._get_all_calculation_managers()
        # self.base_dir is set in material_dirs.setter
        self.results_path = self.base_dir / "results.json"
        self.results = None

    @property
    def calculation_types(self) -> list[CalculationType]:
        return self._calculation_types

    @calculation_types.setter
    def calculation_types(self, values: list[CalculationType]) -> None:
        if not isinstance(values, list):
            raise TypeError("calculation_types must be a list")
        proper_order: list[CalculationType] = [
            "rlx-coarse",
            "rlx",
            "static",
            "bulkmod",
            "elastic",
        ]
        for calc_type in values:
            if calc_type not in proper_order:
                raise ValueError(f"Calculation type {calc_type} not supported")
        sorted_values: list[CalculationType] = [
            calc_type for calc_type in proper_order if calc_type in values
        ]
        self._calculation_types = sorted_values

    @property
    def ncore(self) -> int:
        return self._ncore

    @ncore.setter
    def ncore(self, value: None | int) -> None:
        if value is None:
            value = int(min(len(self.material_dirs), 4))
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
    def calculation_manager_kwargs(self) -> dict:
        return self._calculation_manager_kwargs

    @calculation_manager_kwargs.setter
    def calculation_manager_kwargs(self, values: dict) -> None:
        if not isinstance(values, dict):
            raise TypeError("calculation_manager_kwargs must be a dictionary")

        supported_kwargs = ["primitive", "from_scratch", "strains"]
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
    def material_dirs(self) -> Filepaths:
        return self._material_dirs

    @material_dirs.setter
    def material_dirs(self, values: Filepaths | WorkingDirectory) -> None:
        """
        Sets paths for all materials

        Args:
            values: list of material directory
                paths OR name of calculations dir

                if is list, use that list directly
                if is string, find folders inside of that directory
        """
        match values:
            case str() | Path():
                self.base_dir = Path(values)
                material_dirs = [d for d in Path(values).glob("*") if d.is_dir()]
            case list() | np.ndarray():
                values_as_paths = [Path(p) for p in values]
                base_dir = values_as_paths[0].parent
                for mat_dir in values_as_paths:
                    if not mat_dir.parent == base_dir:
                        raise Exception("All material_dirs must be in the same directory")
                self.base_dir = base_dir
                material_dirs = values_as_paths
            case _:
                raise TypeError(
                    "material_dirs must be a directory name or a list of directories"
                )
        # Sort the paths by name
        self._material_dirs = sorted(material_dirs, key=lambda x: self.sort_by(x.name))

    def _get_material_name_from_path(self, material_dir: Path) -> str:
        return material_dir.name

    @cached_property
    def material_names(self) -> list[str]:
        return [self._get_material_name_from_path(mpath) for mpath in self.material_dirs]

    @property
    def results(self) -> dict:
        return self._results

    @results.setter
    def results(self, value) -> None:
        if self.results_path.exists():
            with open(self.results_path, "r") as fr:
                self._results = json.load(fr)
            for mat_name in self.material_names:
                if mat_name not in self._results:
                    self._results[mat_name] = {}
            self._results = dict(
                sorted(self._results.items(), key=lambda kv: self.sort_by(kv[0]))
            )
        else:
            self._results = {mat_name: {} for mat_name in self.material_names}

    def _check_calc_by_result(
        self,
        material_name: str,
        calc_type: str,
    ) -> tuple[bool, bool]:
        """
        Checks if job has been completed and analyzed

        Args:
            material_name: name of material to check
            calc_type: calculation type to check

        Returns:
            (is_done, is_stopped)
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

    def _get_calculation_managers(self, material_dir: Path) -> list[CalculationManager]:
        """
        Gets calculation managers for a single material
        """
        calc_managers: list[CalculationManager] = []
        for calc_type in self.calculation_types:
            manager: CalculationManager
            match calc_type:
                case "rlx-coarse":
                    manager = RlxCoarseCalculationManager(
                        material_dir=material_dir,
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
                        material_dir=material_dir,
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
                        material_dir=material_dir,
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
                        material_dir=material_dir,
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
                        material_dir=material_dir,
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

    def _get_all_calculation_managers(self) -> dict[str, list[CalculationManager]]:
        """
        Gets calculation managers for all materials
        """
        calc_managers = {}
        for material_name, material_dir in zip(self.material_names, self.material_dirs):
            calc_managers[material_name] = self._get_calculation_managers(material_dir)
        return calc_managers

    def _manage_calculations(self, material_name: str) -> tuple[str, None | str | dict]:
        """
        Runs vasp job workflow for a single material
        """
        material_results: dict[str, Any] = {}
        for calc_manager in self.calculation_managers[material_name]:
            if calc_manager.mode in self.results[material_name].keys():
                calc_is_done, calc_is_stopped = self._check_calc_by_result(
                    material_name, calc_manager.mode
                )
                if calc_is_done and not calc_manager.from_scratch:
                    logger.info(
                        f"{material_name} -- {calc_manager.mode.upper()} Successful"
                    )
                    continue

            if calc_manager.stopped:
                logger.info(f"{material_name} -- {calc_manager.mode.upper()} STOPPED")
                material_results[calc_manager.mode] = "STOPPED"
                break

            if not calc_manager.job_exists:
                logger.info(f"{material_name} -- Setting up {calc_manager.mode.upper()}")
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

            material_results[calc_manager.mode] = calc_manager.results
        return (material_name, material_results)

    def _manage_calculations_wrapper(self) -> list[tuple]:
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

    def run_calculations(self) -> dict:
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

    def summary(
        self,
        as_string: bool = True,
        print_unfinished: bool = False,
        print_stopped: bool = True,
    ) -> str | dict:
        """
        Create a string summary of all calculations

        Args:
            as_string: if True, return string summary. Else, return dict summary
            print_unfinished: if True, include a list of unfinished
                materials for each calculation type in the summary
            print_stopped: if True, include a list of stopped materials for
                each calculation type in the summary
        """
        if not self.results_path.exists():
            raise ValueError(f"Can't find results in {self.results_path}")
        with open(self.results_path) as fr:
            results = json.load(fr)

        summary_dict: dict[str, Any] = {}
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
                        summary_str += f"Unfinished {calc_type.upper()}: {unfinished}\n"
                if print_stopped:
                    stopped = summary_dict[calc_type]["stopped"]
                    if len(stopped) != 0:
                        summary_str += f"\tStopped {calc_type.upper()}: {stopped}\n"
            return summary_str
        else:
            return summary_dict
