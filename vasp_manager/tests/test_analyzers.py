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

# stable_material_dir and unstable_material_dir fixtures are in conftest.py


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


# ---------------------------------------------------------------------------
# Failure modes: missing or malformed output files
# ---------------------------------------------------------------------------


def test_elastic_analyzer_missing_outcar_raises(tmp_path):
    """
    from_calc_dir raises FileNotFoundError when no OUTCAR* is present
    """
    import shutil

    import importlib_resources

    material_dir = (
        importlib_resources.files("vasp_manager") / "tests" / "calculations" / "material"
    )
    elastic_src = material_dir / "elastic"
    elastic_dir = tmp_path / "elastic"
    shutil.copytree(
        str(elastic_src),
        str(elastic_dir),
        symlinks=True,
    )
    # Remove the OUTCAR
    for outcar in elastic_dir.glob("OUTCAR*"):
        outcar.unlink()

    with pytest.raises(FileNotFoundError, match="No OUTCAR"):
        ElasticAnalyzer.from_calc_dir(elastic_dir)


def test_elastic_analyzer_missing_poscar_raises(tmp_path):
    """
    from_calc_dir raises an error when POSCAR is absent
    """
    import shutil

    import importlib_resources

    material_dir = (
        importlib_resources.files("vasp_manager") / "tests" / "calculations" / "material"
    )
    elastic_src = material_dir / "elastic"
    elastic_dir = tmp_path / "elastic"
    shutil.copytree(
        str(elastic_src),
        str(elastic_dir),
        symlinks=True,
    )
    poscar_path = elastic_dir / "POSCAR"
    poscar_path.unlink()

    with pytest.raises(Exception):
        ElasticAnalyzer.from_calc_dir(elastic_dir)


def test_elastic_analyzer_nonexistent_calc_dir_raises():
    """
    from_calc_dir raises ValueError for a directory that doesn't exist
    """
    with pytest.raises(ValueError, match="does not exist"):
        ElasticAnalyzer.from_calc_dir("/does/not/exist/elastic")


def test_elastic_analyzer_truncated_outcar_raises(tmp_path):
    """
    from_calc_dir raises RuntimeError when OUTCAR contains no elastic data
    """
    import shutil

    import importlib_resources

    material_dir = (
        importlib_resources.files("vasp_manager") / "tests" / "calculations" / "material"
    )
    elastic_src = material_dir / "elastic"
    elastic_dir = tmp_path / "elastic"
    shutil.copytree(
        str(elastic_src),
        str(elastic_dir),
        symlinks=True,
    )
    # Replace the OUTCAR with a truncated/empty one (no elastic tensor block)
    for outcar in elastic_dir.glob("OUTCAR*"):
        outcar.unlink()
    outcar_path = elastic_dir / "OUTCAR"
    outcar_path.write_text("truncated OUTCAR with no elastic constants\n")
    # Also remove any cached elastic_constants.txt
    ec_file = elastic_dir / "elastic_constants.txt"
    if ec_file.exists():
        ec_file.unlink()

    with pytest.raises((RuntimeError, AssertionError)):
        ElasticAnalyzer.from_calc_dir(elastic_dir)


def test_elastic_analyzer_cij_wrong_shape_raises(stable_material_dir):
    """
    ElasticAnalyzer raises ValueError if cij is not 6x6
    """
    import numpy as np
    from pymatgen.core import Structure

    poscar_path = stable_material_dir / "elastic" / "POSCAR"
    structure = Structure.from_file(poscar_path)
    bad_cij = np.ones((3, 3), dtype=float)
    with pytest.raises(ValueError, match="6x6"):
        ElasticAnalyzer(cij=bad_cij, structure=structure)
