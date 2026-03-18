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

# Material paths in calculations/:
# 1) NaCl: completed successfully (no spin, no DFT+U, stable)
# 2) Fe: completed successfully (spin, stable)
# 3) NiO: completed successfully (spin + DFT+U, stable)
# 4) BCC_Ti: completed successfully (spin -> no spin re-relax, unstable elastic)
# 5) NaCl_needs_rerun: NaCl with failed runs to rerun from scratch
# 6) NaCl_needs_archive: NaCl with failed runs to make archive directories
# 7) NaCl_hit_errors: NaCl with errors in stdout.txt and stderr.txt
# 8) BCC_Ti_needs_archive: BCC_Ti with spin-polarized rlx to trigger re-relax

INPUT_FILES = ["POSCAR", "POTCAR", "INCAR", "vasp.q"]


# ---------------------------------------------------------------------------
# Do testing of the results first as they don't modify the directories
# ---------------------------------------------------------------------------


def test_rlx_coarse_results(calcs_dir):
    """
    Assert rlx-coarse results are parsed properly
    """
    material_dir = calcs_dir / "NaCl"
    rlx_coarse_dir = material_dir / "rlx-coarse"
    assert rlx_coarse_dir.exists()
    rlx_coarse_manager = RlxCoarseCalculationManager(
        material_dir=material_dir,
        to_rerun=True,
        to_submit=False,
    )
    assert rlx_coarse_manager.is_done
    assert rlx_coarse_manager.results == "done"


def test_rlx_results(calcs_dir):
    """
    Assert rlx results are parsed properly
    """
    material_dir = calcs_dir / "NaCl"
    rlx_dir = material_dir / "rlx"
    assert rlx_dir.exists()
    rlx_manager = RlxCalculationManager(
        material_dir=material_dir,
        to_rerun=True,
        to_submit=False,
    )
    assert rlx_manager.is_done
    assert rlx_manager.results["initial_spacegroup"] == 225
    assert rlx_manager.results["relaxed_spacegroup"] == 225
    assert rlx_manager.results["total_dV"] == 0.0372


def test_static_results(calcs_dir):
    """
    Assert static results are parsed properly
    """
    material_dir = calcs_dir / "NaCl"
    static_dir = material_dir / "static"
    assert static_dir.exists()
    static_manager = StaticCalculationManager(
        material_dir=material_dir,
        to_rerun=True,
        to_submit=False,
    )
    assert static_manager.is_done
    assert static_manager.results["final_energy"] == -6.777893
    assert static_manager.results["final_energy_pa"] == -3.3889465
    assert static_manager.results["magmom_pa"] is None


def test_static_spin_results(calcs_dir):
    """
    Assert static results w/ spin are parsed properly for Fe
    """
    material_dir = calcs_dir / "Fe"
    static_dir = material_dir / "static"
    assert static_dir.exists()
    static_manager = StaticCalculationManager(
        material_dir=material_dir,
        to_rerun=True,
        to_submit=False,
    )
    assert static_manager.is_done
    assert static_manager.results["final_energy"] == -8.2423959
    assert static_manager.results["final_energy_pa"] == -8.2423959
    assert static_manager.results["magmom_pa"] == 2.1782


def test_static_dftu_results(calcs_dir):
    """
    Assert static results w/ spin and DFT+U are parsed properly for NiO
    """
    material_dir = calcs_dir / "NiO"
    static_dir = material_dir / "static"
    assert static_dir.exists()
    static_manager = StaticCalculationManager(
        material_dir=material_dir,
        to_rerun=True,
        to_submit=False,
    )
    assert static_manager.is_done
    assert static_manager.results["final_energy"] == -10.139662
    assert static_manager.results["final_energy_pa"] == -5.069831
    assert static_manager.results["magmom_pa"] == 1.0


def test_static_rerelax_results(calcs_dir):
    """
    Assert static results for BCC_Ti (re-relaxed without spin)
    """
    material_dir = calcs_dir / "BCC_Ti"
    static_dir = material_dir / "static"
    assert static_dir.exists()
    static_manager = StaticCalculationManager(
        material_dir=material_dir,
        to_rerun=True,
        to_submit=False,
    )
    assert static_manager.is_done
    assert static_manager.results["final_energy"] == -7.7249028
    assert static_manager.results["final_energy_pa"] == -7.7249028
    assert static_manager.results["magmom_pa"] is None


def test_bulkmod_results(calcs_dir):
    """
    Assert bulkmod results are parsed properly
    """
    material_dir = calcs_dir / "NaCl"
    bulkmod_dir = material_dir / "bulkmod"
    assert bulkmod_dir.exists()
    bulkmod_manager = BulkmodCalculationManager(
        material_dir=material_dir,
        to_rerun=True,
        to_submit=False,
    )
    assert bulkmod_manager.is_done


def test_elastic_results(calcs_dir):
    """
    Assert elastic results are parsed properly
    """
    material_dir = calcs_dir / "NaCl"
    elastic_dir = material_dir / "elastic"
    assert elastic_dir.exists()
    elastic_manager = ElasticCalculationManager(
        material_dir=material_dir,
        to_rerun=True,
        to_submit=False,
    )
    assert elastic_manager.is_done


# ---------------------------------------------------------------------------
# Do testing of the errors next as they are independent
# ---------------------------------------------------------------------------


def test_hit_errors_and_restart(calcs_dir):
    """
    Test parsing and handling of VASP errors
    """
    material_dir = calcs_dir / "NaCl_hit_errors"
    rlx_coarse_dir = material_dir / "rlx-coarse"
    assert rlx_coarse_dir.exists()
    rlx_coarse_manager = RlxCoarseCalculationManager(
        material_dir=material_dir,
        to_rerun=False,
        to_submit=False,
    )
    errors = rlx_coarse_manager.vasp_run.check_vasp_errors()
    all_errors_addressed = rlx_coarse_manager.vasp_run.address_vasp_errors(errors)
    for error in ["Sub-Space-Matrix", "Inconsistent Bravais", "oom-kill"]:
        assert error in errors
    assert all_errors_addressed
    assert not rlx_coarse_manager.is_done
    assert not rlx_coarse_manager.stopped
    assert rlx_coarse_manager.vasp_input_creator.calc_config.algo == "Fast"
    assert rlx_coarse_manager.vasp_input_creator.calc_config.symprec == "1e-08"
    assert rlx_coarse_manager.vasp_input_creator.ncore_per_node_for_memory == 64
    assert rlx_coarse_manager.results == "not finished"


def test_hit_errors_and_stop(calcs_dir):
    """
    Test parsing and handling of VASP errors
    """
    material_dir = calcs_dir / "NaCl_hit_errors"
    rlx_dir = material_dir / "rlx"
    assert rlx_dir.exists()
    rlx_manager = RlxCalculationManager(
        material_dir=material_dir,
        to_rerun=False,
        to_submit=False,
    )
    errors = rlx_manager.vasp_run.check_vasp_errors()
    all_errors_addressed = rlx_manager.vasp_run.address_vasp_errors(errors)
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


def test_rlx_coarse_hit_errors(calcs_dir):
    """
    Test parsing and handling of VASP errors for rlx-coarse
    """
    material_dir = calcs_dir / "NaCl_hit_errors"
    rlx_coarse_dir = material_dir / "rlx-coarse"
    assert rlx_coarse_dir.exists()
    rlx_coarse_manager = RlxCoarseCalculationManager(
        material_dir=material_dir,
        to_rerun=False,
        to_submit=False,
    )
    assert not rlx_coarse_manager.is_done


def test_rlx_hit_errors(calcs_dir):
    """
    Test parsing and handling of VASP errors for rlx
    """
    material_dir = calcs_dir / "NaCl_hit_errors"
    rlx_dir = material_dir / "rlx"
    assert rlx_dir.exists()
    rlx_manager = RlxCalculationManager(
        material_dir=material_dir,
        to_rerun=False,
        to_submit=False,
    )
    assert not rlx_manager.is_done


def test_static_hit_errors(calcs_dir):
    """
    Test parsing and handling of VASP errors for static
    """
    material_dir = calcs_dir / "NaCl_hit_errors"
    static_dir = material_dir / "static"
    assert static_dir.exists()
    static_manager = StaticCalculationManager(
        material_dir=material_dir,
        to_rerun=False,
        to_submit=False,
    )
    assert not static_manager.is_done


def test_bulkmod_hit_errors(calcs_dir):
    """
    Test parsing and handling of VASP errors for bulkmod
    """
    material_dir = calcs_dir / "NaCl_hit_errors"
    bulkmod_dir = material_dir / "bulkmod"
    assert bulkmod_dir.exists()
    bulkmod_manager = BulkmodCalculationManager(
        material_dir=material_dir,
        to_rerun=False,
        to_submit=False,
    )
    assert not bulkmod_manager.is_done


def test_elastic_hit_errors(calcs_dir):
    """
    Test parsing and handling of VASP errors for elastic
    """
    material_dir = calcs_dir / "NaCl_hit_errors"
    elastic_dir = material_dir / "elastic"
    assert elastic_dir.exists()
    elastic_manager = ElasticCalculationManager(
        material_dir=material_dir,
        to_rerun=False,
        to_submit=False,
    )
    assert not elastic_manager.is_done


def test_parse_magmom(calcs_dir):
    """
    Test parsing of magmom
    """
    for material_name, expected_magmom_pa in zip(
        ["NaCl", "Fe", "NiO", "BCC_Ti"], [None, 2.1486, 1.0, None]
    ):
        material_dir = calcs_dir / material_name
        rlx_dir = material_dir / "rlx"
        assert rlx_dir.exists()
        rlx_manager = RlxCalculationManager(
            material_dir=material_dir,
            to_rerun=True,
            to_submit=False,
        )
        assert rlx_manager.vasp_run.parse_magmom_per_atom() == expected_magmom_pa


# ---------------------------------------------------------------------------
# Do testing of the reruns/archives last as they modify their directories
# ---------------------------------------------------------------------------


def test_too_many_rlx_coarse_archives(calcs_dir):
    """
    Assert rlx-coarse handles too many archives correctly
        (Continues to rlx)
    """
    material_dir = calcs_dir / "NaCl_needs_archive"
    rlx_coarse_dir = material_dir / "rlx-coarse"
    assert rlx_coarse_dir.exists()
    rlx_coarse_manager = RlxCoarseCalculationManager(
        material_dir=material_dir, to_rerun=True, to_submit=False, max_reruns=2
    )
    # archive_0 already exists but allow material to continue to rlx instead
    assert rlx_coarse_manager.is_done


def test_too_many_rlx_archives(calcs_dir):
    """
    Assert rlx handles too many archives correctly
        (Errors out and refuses to continue)
    """
    material_dir = calcs_dir / "NaCl_needs_archive"
    rlx_dir = material_dir / "rlx"
    assert rlx_dir.exists()
    rlx_manager = RlxCalculationManager(
        material_dir=material_dir, to_rerun=True, to_submit=False, max_reruns=2
    )
    # archive_0 already exists
    assert not rlx_manager.is_done


def test_rlx_coarse_archive(calcs_dir):
    """
    Assert unfinished rlx-coarse makes archive and that the symmetrized
        structure for the new POSCAR matches the CONTCAR in the archive
    """
    material_dir = calcs_dir / "NaCl_needs_archive"
    rlx_coarse_dir = material_dir / "rlx-coarse"
    assert rlx_coarse_dir.exists()
    rlx_coarse_manager = RlxCoarseCalculationManager(
        material_dir=material_dir,
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


def test_rlx_archive(calcs_dir):
    """
    Assert unfinished rlx makes archive
    Also ensure rlx that finishes with magmom_pa below magmom_per_atom_cutoff
        restarts without spin
    """
    material_dir = calcs_dir / "BCC_Ti_needs_archive"
    rlx_dir = material_dir / "rlx"
    assert rlx_dir.exists()
    rlx_manager = RlxCalculationManager(
        material_dir=material_dir,
        to_rerun=True,
        to_submit=False,
        magmom_per_atom_cutoff=1.0,
    )
    assert not rlx_manager.is_done
    assert rlx_manager.results is None
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


def test_static_rerun(calcs_dir):
    """
    Assert static reruns from scratch if failed
    """
    material_dir = calcs_dir / "NaCl_needs_rerun"
    static_dir = material_dir / "static"
    assert static_dir.exists()
    static_manager = StaticCalculationManager(
        material_dir=material_dir,
        to_rerun=True,
        to_submit=False,
    )
    assert not static_manager.is_done
    assert static_manager.results is None
    static_dir_files = list(static_dir.glob("*"))
    assert len(static_dir_files) == 4
    for file in static_dir_files:
        assert file.name in INPUT_FILES


def test_bulkmod_rerun(calcs_dir):
    """
    Assert bulkmod reruns only the failed strain, not all strains
    """
    material_dir = calcs_dir / "NaCl_needs_rerun"
    bulkmod_dir = material_dir / "bulkmod"
    assert bulkmod_dir.exists()
    bulkmod_manager = BulkmodCalculationManager(
        material_dir=material_dir,
        to_rerun=True,
        to_submit=False,
    )
    assert not bulkmod_manager.is_done
    assert bulkmod_manager.results is None

    # The failed strain (strain_-1) should have been cleaned and recreated
    failed_strain_dir = bulkmod_dir / "strain_-1"
    assert failed_strain_dir.exists()
    failed_files = [f.name for f in failed_strain_dir.glob("*") if f.is_file()]
    for f_name in INPUT_FILES:
        assert f_name in failed_files


def test_elastic_rerun(calcs_dir):
    """
    Assert elastic reruns from scratch if failed
    """
    material_dir = calcs_dir / "NaCl_needs_rerun"
    elastic_dir = material_dir / "elastic"
    assert elastic_dir.exists()
    elastic_manager = ElasticCalculationManager(
        material_dir=material_dir,
        to_rerun=True,
        to_submit=False,
    )
    assert not elastic_manager.is_done
    assert elastic_manager.results is None
