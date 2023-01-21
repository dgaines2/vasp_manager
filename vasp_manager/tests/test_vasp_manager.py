import glob
import os
import shutil
import subprocess
import sys

from vasp_manager import VaspManager
from vasp_manager.utils import change_directory

sys.path.append("..")
from run_vasp_calculations import make_calculations_folder


def test_vmg_in_order(tmpdir):
    original_calculations_folder = "calculations"
    temp_calculation_folder = tmpdir
    if os.path.exists(temp_calculation_folder):
        shutil.rmtree(temp_calculation_folder)
    make_calculations_folder(calcs_path=temp_calculation_folder)

    for f in ["computing_config.json", "calc_config.json", "unzip_outputs.sh"]:
        shutil.copy(
            os.path.join(original_calculations_folder, f), temp_calculation_folder
        )

    calculation_types = [
        "rlx-coarse",
        "rlx",
        "static",
        "bulkmod",
        "elastic",
    ]
    calc_paths = os.path.join(temp_calculation_folder, "*")
    material_paths = [p for p in sorted(glob.glob(calc_paths)) if os.path.isdir(p)]

    for i, calculation_type in enumerate(calculation_types):
        calculation_type_subset = calculation_types[: i + 1]
        vmg = VaspManager(
            calculation_types=calculation_type_subset,
            material_paths=material_paths,
            use_multiprocessing=True,
            to_rerun=True,
            to_submit=True,
        )
        results = vmg.run_calculations()
        if i != 0:
            previous_calc_types = calculation_type_subset[:i]
            # check results
            for material in results:
                for pc in previous_calc_types:
                    match pc:
                        case "rlx-coarse":
                            assert results[material][pc] == "done"
                        case "rlx":
                            assert isinstance(results[material][pc]["total_dV"], float)
                        case "static":
                            assert isinstance(
                                results[material][pc]["final_energy"], float
                            )
                        case "bulkmod":
                            assert isinstance(results[material][pc]["B"], float)
                        case _:
                            raise Exception
            # check summary
            summary = vmg.summary(as_string=False)
            for pc in previous_calc_types:
                assert summary["n_total"] == summary[pc]["n_finished"]

        for p in material_paths:
            mat_name = os.path.basename(p)
            orig_mode_path = os.path.join("calculations", mat_name, calculation_type)
            temp_mode_path = os.path.join(p, calculation_type)
            shutil.rmtree(temp_mode_path)
            shutil.copytree(orig_mode_path, temp_mode_path)
        unzip_call = f"bash unzip_outputs.sh"
        with change_directory(temp_calculation_folder):
            subprocess.call(unzip_call, shell=True)

    vmg = VaspManager(
        calculation_types=calculation_types,
        material_paths=material_paths,
        to_rerun=True,
        to_submit=True,
    )
    results = vmg.run_calculations()
    for material in results:
        assert results[material]["rlx-coarse"] == "done"
        assert isinstance(results[material]["rlx"]["initial_spacegroup"], int)
        assert isinstance(results[material]["rlx"]["relaxed_spacegroup"], int)
        assert isinstance(results[material]["rlx"]["total_dV"], float)
        assert isinstance(results[material]["static"]["final_energy"], float)
        assert isinstance(results[material]["bulkmod"]["B"], float)
        assert isinstance(results[material]["elastic"]["B_VRH"], float)

    summary = vmg.summary(as_string=False)
    for calculation_type in calculation_types:
        assert summary["n_total"] == summary[calculation_type]["n_finished"]


def test_vmg_with_skipping(tmpdir):
    original_calculations_folder = "calculations"
    temp_calculation_folder = tmpdir
    if os.path.exists(temp_calculation_folder):
        shutil.rmtree(temp_calculation_folder)
    make_calculations_folder(calcs_path=temp_calculation_folder)

    for f in ["computing_config.json", "calc_config.json", "unzip_outputs.sh"]:
        shutil.copy(
            os.path.join(original_calculations_folder, f), temp_calculation_folder
        )

    calculation_types = [
        "rlx-coarse",
        "rlx",
        "static",
        "bulkmod",
        "elastic",
    ]
    calc_paths = os.path.join(temp_calculation_folder, "*")
    material_paths = [p for p in sorted(glob.glob(calc_paths)) if os.path.isdir(p)]

    # skip static and bulkmod when copying
    completed_calc_types = ["rlx-coarse", "rlx", "elastic"]
    for calculation_type in completed_calc_types:
        for p in material_paths:
            mat_name = os.path.basename(p)
            orig_mode_path = os.path.join("calculations", mat_name, calculation_type)
            temp_mode_path = os.path.join(p, calculation_type)
            shutil.copytree(orig_mode_path, temp_mode_path)

    unzip_call = f"bash unzip_outputs.sh"
    with change_directory(temp_calculation_folder):
        subprocess.call(unzip_call, shell=True)

    vmg = VaspManager(
        calculation_types=calculation_types,
        material_paths=material_paths,
        use_multiprocessing=True,
        to_rerun=True,
        to_submit=True,
    )
    results = vmg.run_calculations()
    for material in results:
        for calc_type in calculation_types:
            match calc_type:
                case "rlx-coarse":
                    assert results[material][calc_type] == "done"
                case "rlx":
                    assert isinstance(results[material][calc_type]["total_dV"], float)
                case "static" | "bulkmod":
                    assert results[material][calc_type] is None
                case "elastic":
                    assert isinstance(results[material][calc_type]["B_VRH"], float)
                case _:
                    raise Exception

    summary = vmg.summary(as_string=False)
    for pc in calculation_types:
        match pc:
            case "rlx-coarse" | "rlx" | "elastic":
                assert summary["n_total"] == summary[pc]["n_finished"]
            case "bulkmod" | "static":
                assert summary[pc]["n_finished"] == 0
            case _:
                raise Exception
