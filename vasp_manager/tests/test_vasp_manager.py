import shutil
from textwrap import dedent

import importlib_resources

from vasp_manager import VaspManager


def test_vmg_in_order(tmp_path):
    original_calculation_folder = (
        importlib_resources.files("vasp_manager") / "tests" / "calculations"
    )
    temp_calculation_folder = tmp_path / "calculations"
    shutil.copytree(
        original_calculation_folder,
        temp_calculation_folder,
        dirs_exist_ok=True,
        symlinks=True,
        ignore=shutil.ignore_patterns(
            "rlx*",
            "static",
            "bulkmod",
            "elastic",
            "material_hit_errors",
            "material_needs_archive",
            "material_needs_rerun",
        ),
    )

    calculation_types = [
        "rlx-coarse",
        "rlx",
        "static",
        "bulkmod",
        "elastic",
    ]
    material_paths = [
        p for p in sorted(list(temp_calculation_folder.glob("*"))) if p.is_dir()
    ]

    for i, calculation_type in enumerate(calculation_types):
        calculation_type_subset = calculation_types[: i + 1]
        vmg = VaspManager(
            calculation_types=calculation_type_subset,
            material_paths=material_paths,
            use_multiprocessing=False,
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
            mat_name = p.name
            orig_mode_path = original_calculation_folder / mat_name / calculation_type
            temp_mode_path = p / calculation_type
            shutil.copytree(
                orig_mode_path,
                temp_mode_path,
                dirs_exist_ok=True,
                symlinks=True,
                ignore=shutil.ignore_patterns("elastic_constants.txt"),
            )

        # test string summary when unfinished
        if i == (len(calculation_types) - 1):
            summary = vmg.summary(as_string=True, print_unfinished=True)
            correct_summary = dedent(
                """\
            Total Materials = 2
            ------------------------------
            RLX-COARSE  2/2 completed
            RLX         2/2 completed
            STATIC      2/2 completed
            BULKMOD     2/2 completed
            ELASTIC     0/2 completed
                        2 not completed
            Unfinished ELASTIC: ['material', 'material_spinu']
            """
            )
            assert summary == correct_summary

    vmg = VaspManager(
        calculation_types=calculation_types,
        material_paths=material_paths,
        use_multiprocessing=False,
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


def test_vmg_with_skipping(tmp_path):
    original_calculation_folder = (
        importlib_resources.files("vasp_manager") / "tests" / "calculations"
    )
    temp_calculation_folder = tmp_path / "calculations"
    # skip static and bulkmod when copying
    shutil.copytree(
        original_calculation_folder,
        temp_calculation_folder,
        dirs_exist_ok=True,
        symlinks=True,
        ignore=shutil.ignore_patterns(
            "static",
            "bulkmod",
            "material_hit_errors",
            "material_needs_archive",
            "material_needs_rerun",
        ),
    )

    calculation_types = [
        "rlx-coarse",
        "rlx",
        "static",
        "bulkmod",
        "elastic",
    ]
    material_paths = [
        p for p in sorted(list(temp_calculation_folder.glob("*"))) if p.is_dir()
    ]

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
