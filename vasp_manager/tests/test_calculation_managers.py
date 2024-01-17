import filecmp
import shutil

import importlib_resources
import pytest
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Structure

from vasp_manager.calculation_manager import (
    BulkmodCalculationManager,
    ElasticCalculationManager,
    RlxCalculationManager,
    RlxCoarseCalculationManager,
    StaticCalculationManager,
)
from vasp_manager.utils import pgrep

"""
Four material paths exist in calculations/
1) material: NaCl completed successfully
2) material_needs_rerun: NaCl with failed runs to rerun from scratch
3) material_needs_archive: NaCl with failed runs to make archive folders
4) material_hit_erros: NaCl with errors in stdout.txt and stderr.txt
5) material_spinu: CuO with spin and DFT+U completed sucessfully
"""

INPUT_FILES = ["POSCAR", "POTCAR", "INCAR", "vasp.q"]


@pytest.fixture(scope="session")
def calc_dir(tmp_path_factory):
    """
    Complete calculations folder
    """
    calculation_folder = (
        importlib_resources.files("vasp_manager") / "tests" / "calculations"
    )
    new_calculation_folder = tmp_path_factory.mktemp("calculations")
    shutil.copytree(
        calculation_folder,
        new_calculation_folder,
        dirs_exist_ok=True,
        symlinks=True,
    )
    return new_calculation_folder


@pytest.fixture(scope="session")
def base_calc_dir(tmp_path_factory):
    """
    Calculations folder that needs each run type to be set up
    """
    calculation_folder = (
        importlib_resources.files("vasp_manager") / "tests" / "calculations"
    )
    new_calculation_folder = tmp_path_factory.mktemp("calculations")
    shutil.copytree(
        calculation_folder,
        new_calculation_folder,
        dirs_exist_ok=True,
        symlinks=True,
        ignore=shutil.ignore_patterns("rlx*", "static", "bulkmod", "elastic"),
    )
    return new_calculation_folder


"""
Do testing of the results first as they don't modify the folders
"""


def test_rlx_coarse_results(calc_dir):
    """
    Assert rlx-coarse results are parsed properly
    """
    material_path = calc_dir / "material"
    rlx_coarse_path = material_path / "rlx-coarse"
    assert rlx_coarse_path.exists()
    rlx_coarse_manager = RlxCoarseCalculationManager(
        material_path=material_path,
        to_rerun=True,
        to_submit=False,
    )
    assert rlx_coarse_manager.is_done
    assert rlx_coarse_manager.results == "done"


def test_rlx_results(calc_dir):
    """
    Assert rlx results are parsed properly
    """
    material_path = calc_dir / "material"
    rlx_dir = material_path / "rlx"
    assert rlx_dir.exists()
    rlx_manager = RlxCalculationManager(
        material_path=material_path,
        to_rerun=True,
        to_submit=False,
    )
    assert rlx_manager.is_done
    assert rlx_manager.results["initial_spacegroup"] == 225
    assert rlx_manager.results["relaxed_spacegroup"] == 225
    assert rlx_manager.results["total_dV"] == 0.0379


def test_static_results(calc_dir):
    """
    Assert static results are parsed properly
    """
    material_path = calc_dir / "material"
    static_dir = material_path / "static"
    assert static_dir.exists()
    static_manager = StaticCalculationManager(
        material_path=material_path,
        to_rerun=True,
        to_submit=False,
    )
    assert static_manager.is_done
    assert static_manager.results["final_energy"] == -6.7778924
    assert static_manager.results["final_energy_pa"] == -3.3889462
    assert static_manager.results["magmom_pa"] == None


def test_static_spin_results(calc_dir):
    """
    Assert static results w/ spin are parsed properly
    """
    material_path = calc_dir / "material_spinu"
    static_dir = material_path / "static"
    assert static_dir.exists()
    static_manager = StaticCalculationManager(
        material_path=material_path,
        to_rerun=True,
        to_submit=False,
    )
    assert static_manager.is_done
    assert static_manager.results["final_energy"] == -8.1574727
    assert static_manager.results["final_energy_pa"] == -4.07873635
    assert static_manager.results["magmom_pa"] == 0.0001


def test_bulkmod_results(calc_dir):
    """
    Assert bulkmod results are parsed properly
    """
    material_path = calc_dir / "material"
    bulkmod_dir = material_path / "bulkmod"
    assert bulkmod_dir.exists()
    bulkmod_manager = BulkmodCalculationManager(
        material_path=material_path,
        to_rerun=True,
        to_submit=False,
    )
    assert bulkmod_manager.is_done


def test_elastic_results(calc_dir):
    """
    Assert elastic results are parsed properly
    """
    material_path = calc_dir / "material"
    elastic_dir = material_path / "elastic"
    assert elastic_dir.exists()
    elastic_manager = ElasticCalculationManager(
        material_path=material_path,
        to_rerun=True,
        to_submit=False,
    )
    assert elastic_manager.is_done


"""
Do testing of the errors next as they are independent
"""


def test_hit_errors_and_restart(calc_dir):
    """
    Test parsing and handling of VASP errors
    """
    material_path = calc_dir / "material_hit_errors"
    rlx_coarse_dir = material_path / "rlx-coarse"
    stdout_path = rlx_coarse_dir / "stdout.txt"
    stderr_path = rlx_coarse_dir / "stderr.txt"
    assert rlx_coarse_dir.exists()
    rlx_coarse_manager = RlxCoarseCalculationManager(
        material_path=material_path,
        to_rerun=False,
        to_submit=False,
    )
    errors = rlx_coarse_manager._check_vasp_errors(stdout_path, stderr_path)
    all_errors_addressed = rlx_coarse_manager._address_vasp_errors(errors)
    for error in ["Sub-Space-Matrix", "Inconsistent Bravais", "oom-kill"]:
        assert error in errors
    assert all_errors_addressed
    assert not rlx_coarse_manager.is_done
    assert not rlx_coarse_manager.stopped
    assert rlx_coarse_manager.vasp_input_creator.calc_config["algo"] == "Fast"
    assert rlx_coarse_manager.vasp_input_creator.calc_config["symprec"] == "1e-08"
    assert rlx_coarse_manager.vasp_input_creator.ncore_per_node_for_memory == 64
    assert rlx_coarse_manager.results == "not finished"


def test_hit_errors_and_stop(calc_dir):
    """
    Test parsing and handling of VASP errors
    """
    material_path = calc_dir / "material_hit_errors"
    rlx_dir = material_path / "rlx"
    stdout_path = rlx_dir / "stdout.txt"
    stderr_path = rlx_dir / "stderr.txt"
    assert rlx_dir.exists()
    rlx_manager = RlxCalculationManager(
        material_path=material_path,
        to_rerun=False,
        to_submit=False,
    )
    errors = rlx_manager._check_vasp_errors(stdout_path, stderr_path)
    all_errors_addressed = rlx_manager._address_vasp_errors(errors)
    for error in [
        "num prob",
        "BRMIX",
        "SICK JOB",
        "VERY BAD NEWS",
        "Fatal error",
        "SETYLM",
        "Segmentation",
        "command not found",
    ]:
        assert error in errors
    assert not all_errors_addressed
    assert not rlx_manager.is_done
    assert rlx_manager.stopped
    assert rlx_manager.results == "STOPPED"


def test_rlx_coarse_hit_errors(calc_dir):
    """
    Test parsing and handling of VASP errors for rlx-coarse
    """
    material_path = calc_dir / "material_hit_errors"
    rlx_coarse_dir = material_path / "rlx-coarse"
    assert rlx_coarse_dir.exists()
    rlx_coarse_manager = RlxCoarseCalculationManager(
        material_path=material_path,
        to_rerun=False,
        to_submit=False,
    )
    assert not rlx_coarse_manager.is_done


def test_rlx_hit_errors(calc_dir):
    """
    Test parsing and handling of VASP errors for rlx
    """
    material_path = calc_dir / "material_hit_errors"
    rlx_dir = material_path / "rlx"
    assert rlx_dir.exists()
    rlx_manager = RlxCalculationManager(
        material_path=material_path,
        to_rerun=False,
        to_submit=False,
    )
    assert not rlx_manager.is_done


def test_static_hit_errors(calc_dir):
    """
    Test parsing and handling of VASP errors for static
    """
    material_path = calc_dir / "material_hit_errors"
    static_dir = material_path / "static"
    assert static_dir.exists()
    static_manager = StaticCalculationManager(
        material_path=material_path,
        to_rerun=False,
        to_submit=False,
    )
    assert not static_manager.is_done


def test_bulkmod_hit_errors(calc_dir):
    """
    Test parsing and handling of VASP errors for bulkmod
    """
    material_path = calc_dir / "material_hit_errors"
    bulkmod_dir = material_path / "bulkmod"
    assert bulkmod_dir.exists()
    bulkmod_manager = BulkmodCalculationManager(
        material_path=material_path,
        to_rerun=False,
        to_submit=False,
    )
    assert not bulkmod_manager.is_done


def test_elastic_hit_errors(calc_dir):
    """
    Test parsing and handling of VASP errors for elastic
    """
    material_path = calc_dir / "material_hit_errors"
    elastic_dir = material_path / "elastic"
    assert elastic_dir.exists()
    elastic_manager = ElasticCalculationManager(
        material_path=material_path,
        to_rerun=False,
        to_submit=False,
    )
    assert not elastic_manager.is_done


def test_parse_magmom(calc_dir):
    """
    Test parsing of magmom
    """
    for material_name, expected_magmom_pa in zip(
        ["material", "material_spinu"], [None, -0.0]
    ):
        material_path = calc_dir / material_name
        rlx_dir = material_path / "rlx"
        assert rlx_dir.exists()
        rlx_manager = RlxCalculationManager(
            material_path=material_path,
            to_rerun=True,
            to_submit=False,
        )
        assert rlx_manager._parse_magmom_per_atom() == expected_magmom_pa


"""
Do testing of the reruns/archives last as they modify their folders
"""


def test_too_many_rlx_coarse_archives(calc_dir):
    """
    Assert rlx-coarse handles too many archives correctly
        (Continues to rlx)
    """
    material_path = calc_dir / "material_needs_archive"
    rlx_coarse_dir = material_path / "rlx-coarse"
    assert rlx_coarse_dir.exists()
    rlx_coarse_manager = RlxCoarseCalculationManager(
        material_path=material_path, to_rerun=True, to_submit=False, max_reruns=2
    )
    # archive_0 already exists but allow material to continue to rlx instead
    assert rlx_coarse_manager.is_done


def test_too_many_rlx_archives(calc_dir):
    """
    Assert rlx handles too many archives correctly
        (Errors out and refuses to continue)
    """
    material_path = calc_dir / "material_needs_archive"
    rlx_dir = material_path / "rlx"
    assert rlx_dir.exists()
    rlx_manager = RlxCalculationManager(
        material_path=material_path, to_rerun=True, to_submit=False, max_reruns=2
    )
    # archive_0 already exists
    assert not rlx_manager.is_done


def test_rlx_coarse_archive(calc_dir):
    """
    Assert unfinished rlx-coarse makes archive and that the symmetrized
        structure for the new POSCAR matches the CONTCAR in the archive
    """
    material_path = calc_dir / "material_needs_archive"
    rlx_coarse_dir = material_path / "rlx-coarse"
    assert rlx_coarse_dir.exists()
    rlx_coarse_manager = RlxCoarseCalculationManager(
        material_path=material_path,
        to_rerun=True,
        to_submit=False,
    )
    assert not rlx_coarse_manager.is_done
    assert rlx_coarse_manager.results == "not finished"
    # check archive made correctly
    archive_dir = rlx_coarse_dir / "archive_1"
    assert archive_dir.exists()
    for dir in [rlx_coarse_dir, archive_dir]:
        for f_name in INPUT_FILES:
            assert (dir / f_name).exists()
    archive_dir_contcar = Structure.from_file(archive_dir / "CONTCAR")
    rlx_coarse_poscar = Structure.from_file(rlx_coarse_dir / "POSCAR")
    # compare structures to make sure they're the same
    structure_matcher = StructureMatcher(ltol=0.01, stol=0.1, angle_tol=1.0)
    structures_are_equivalent = structure_matcher.fit(
        archive_dir_contcar, rlx_coarse_poscar
    )
    assert structures_are_equivalent


def test_rlx_archive(calc_dir):
    """
    Assert unfinished rlx makes archive
    Also ensure rlx that finishes with magmom_pa below magmom_per_atom_cutoff
        restarts without spin
    """
    material_path = calc_dir / "material_spinu"
    rlx_dir = material_path / "rlx"
    assert rlx_dir.exists()
    rlx_manager = RlxCalculationManager(
        material_path=material_path,
        to_rerun=True,
        to_submit=False,
        magmom_per_atom_cutoff=1.0,
    )
    assert not rlx_manager.is_done
    assert rlx_manager.results == None
    assert not rlx_manager.vasp_input_creator.use_spin
    # check archive made correctly
    archive_dir = rlx_dir / "archive_0"
    assert archive_dir.exists()
    for dir in [rlx_dir, archive_dir]:
        for f_name in INPUT_FILES:
            assert (dir / f_name).exists()
    incar_path = rlx_dir / "INCAR"
    ispin_line = pgrep(incar_path, "ISPIN", as_string=True)
    ispin_setting = int(ispin_line.split("=")[1].strip())
    assert ispin_setting == 1


def test_static_rerun(calc_dir):
    """
    Assert static reruns from scratch if failed
    """
    material_path = calc_dir / "material_needs_rerun"
    static_dir = material_path / "static"
    assert static_dir.exists()
    static_manager = StaticCalculationManager(
        material_path=material_path,
        to_rerun=True,
        to_submit=False,
    )
    assert not static_manager.is_done
    assert static_manager.results == None
    static_dir_files = list(static_dir.glob("*"))
    assert len(static_dir_files) == 4
    for file in static_dir_files:
        assert file.name in INPUT_FILES


def test_bulkmod_rerun(calc_dir):
    """
    Assert bulkmod reruns from scratch if failed
    """
    material_path = calc_dir / "material_needs_rerun"
    bulkmod_dir = material_path / "bulkmod"
    assert bulkmod_dir.exists()
    bulkmod_manager = BulkmodCalculationManager(
        material_path=material_path,
        to_rerun=True,
        to_submit=False,
    )
    assert not bulkmod_manager.is_done
    assert bulkmod_manager.results == None

    bulkmod_dir_contents = list(bulkmod_dir.glob("*"))
    bulkmod_dir_files = [f for f in bulkmod_dir_contents if f.is_file()]
    strain_folders = [d for d in bulkmod_dir_contents if d.is_dir()]
    assert len(bulkmod_dir_files) == 4
    for file in bulkmod_dir_files:
        assert file.name in INPUT_FILES
    assert len(strain_folders) == len(bulkmod_manager.strains)
    for folder in strain_folders:
        strain_folder_files = list(folder.glob("*"))
        # no vasp.q
        assert len(strain_folder_files) == 3
        for file in strain_folder_files:
            assert file.name in INPUT_FILES


def test_elastic_rerun(calc_dir):
    """
    Assert elastic reruns from scratch if failed
    """
    material_path = calc_dir / "material_needs_rerun"
    elastic_dir = material_path / "elastic"
    assert elastic_dir.exists()
    elastic_manager = ElasticCalculationManager(
        material_path=material_path,
        to_rerun=True,
        to_submit=False,
    )
    assert not elastic_manager.is_done
    assert elastic_manager.results == None
