import shutil

import importlib_resources
import numpy as np
import pytest

from vasp_manager.analyzer import ElasticAnalyzer

SUPPORTED_ELASTIC_PROPERTIES = [
    "B_Reuss",
    "B_Voigt",
    "B_VRH",
    "G_Reuss",
    "G_Voigt",
    "G_VRH",
    "unstable",
    "elastic_tensor",
    "vl",
    "vt",
    "vs",
]


@pytest.fixture(scope="session")
def stable_material_dir(tmp_path_factory):
    material_dir = (
        importlib_resources.files("vasp_manager") / "tests" / "calculations" / "material"
    )
    new_material_dir = tmp_path_factory.mktemp("material")
    shutil.copytree(
        material_dir,
        new_material_dir,
        dirs_exist_ok=True,
        symlinks=True,
        ignore=shutil.ignore_patterns("elastic_constants.txt"),
    )
    return new_material_dir


@pytest.fixture(scope="session")
def unstable_material_dir(tmp_path_factory):
    material_dir = (
        importlib_resources.files("vasp_manager")
        / "tests"
        / "calculations"
        / "material_spinu"
    )
    new_material_dir = tmp_path_factory.mktemp("material_spinu")
    shutil.copytree(
        material_dir,
        new_material_dir,
        dirs_exist_ok=True,
        symlinks=True,
        ignore=shutil.ignore_patterns("elastic_constants.txt"),
    )
    return new_material_dir


def test_elastic_analyzer_stable(stable_material_dir):
    """
    Make sure elastic is parsed properly for an elastically stable material
    """
    elastic_dir = stable_material_dir / "elastic"
    assert elastic_dir.exists()
    ea = ElasticAnalyzer.from_calc_dir(elastic_dir)
    results = ea.results
    for property in results:
        assert property in SUPPORTED_ELASTIC_PROPERTIES

    assert results["B_Reuss"] == 24.007
    assert results["B_Voigt"] == 24.007
    assert results["B_VRH"] == 24.007
    assert results["G_Reuss"] == 14.286
    assert results["G_Voigt"] == 14.734
    assert results["G_VRH"] == 14.51
    assert results["unstable"] == False
    assert np.array_equal(
        results["elastic_tensor"],
        [
            [47.9984, 12.0109, 12.0109, -0.0, 0.0, 0.0],
            [12.0109, 47.9984, 12.0109, -0.0, 0.0, 0.0],
            [12.0109, 12.0109, 47.9984, -0.0, -0.0, 0.0],
            [-0.0, -0.0, -0.0, 12.5609, 0.0, -0.0],
            [0.0, 0.0, -0.0, 0.0, 12.5609, -0.0],
            [0.0, 0.0, 0.0, -0.0, -0.0, 12.5609],
        ],
    )
    assert results["vl"] == 4.545
    assert results["vt"] == 2.629
    assert results["vs"] == 2.919


def test_elastic_analyzer_unstable(unstable_material_dir):
    """
    Make sure elastic is parsed properly for an elastically unstable material
    """
    elastic_dir = unstable_material_dir / "elastic"
    assert elastic_dir.exists()
    ea = ElasticAnalyzer.from_calc_dir(elastic_dir)
    results = ea.results
    for property in results:
        assert property in SUPPORTED_ELASTIC_PROPERTIES

    assert results["B_Reuss"] == 176.275
    assert results["B_Voigt"] == 176.275
    assert results["B_VRH"] == 176.275
    assert results["G_Reuss"] == -199.401
    assert results["G_Voigt"] == 25.595
    assert results["G_VRH"] == -86.903
    assert results["unstable"] == True
    assert np.array_equal(
        results["elastic_tensor"],
        [
            [140.5505, 194.1377, 194.1377, -0.0, -0.0, 0.0],
            [194.1377, 140.5505, 194.1377, -0.0, -0.0, 0.0],
            [194.1377, 194.1377, 140.5505, -0.0, 0.0, 0.0],
            [-0.0, -0.0, -0.0, 60.5209, 0.0, 0.0],
            [-0.0, -0.0, 0.0, 0.0, 60.5209, -0.0],
            [0.0, 0.0, 0.0, -0.0, -0.0, 60.5209],
        ],
    )
    assert results["vl"] == 2.933
    assert np.isnan(results["vt"])
    assert np.isnan(results["vs"])
