# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import shutil

import importlib_resources
import pytest

from vasp_manager.utils import get_pmg_structure_from_poscar
from vasp_manager.vasp_input_creator import VaspInputCreator


@pytest.fixture(scope="session")
def nacl_poscar_path():
    """
    Path to NaCl POSCAR from test data
    """
    return (
        importlib_resources.files("vasp_manager")
        / "tests"
        / "calculations"
        / "NaCl"
        / "POSCAR"
    )


@pytest.fixture(scope="session")
def nio_poscar_path():
    """
    Path to NiO POSCAR from test data (has d-block elements)
    """
    return (
        importlib_resources.files("vasp_manager")
        / "tests"
        / "calculations"
        / "NiO"
        / "POSCAR"
    )


@pytest.fixture(scope="session")
def config_dir():
    """
    Path to the shared test config directory (has calc_config.json, computing_config.json)
    """
    return importlib_resources.files("vasp_manager") / "tests" / "calculations"


@pytest.fixture
def nacl_vic(tmp_path, nacl_poscar_path, config_dir):
    """
    VaspInputCreator for NaCl, rlx-coarse mode, writing to a fresh tmp dir
    """
    calc_dir = tmp_path / "NaCl" / "rlx-coarse"
    calc_dir.mkdir(parents=True)
    # copy POTCAR so make_incar can read it
    potcar_src = (
        importlib_resources.files("vasp_manager")
        / "tests"
        / "calculations"
        / "NaCl"
        / "rlx-coarse"
        / "POTCAR"
    )
    potcar_dst = calc_dir / "POTCAR"
    shutil.copy(str(potcar_src), str(potcar_dst))
    structure = get_pmg_structure_from_poscar(nacl_poscar_path, primitive=True)
    return VaspInputCreator(
        calc_dir=calc_dir,
        mode="rlx-coarse",
        structure=structure,
        config_dir=config_dir,
        job_prefix="rc",
    )


@pytest.fixture
def nio_vic(tmp_path, nio_poscar_path, config_dir):
    """
    VaspInputCreator for NiO, rlx mode
    """
    calc_dir = tmp_path / "NiO" / "rlx"
    calc_dir.mkdir(parents=True)
    potcar_src = (
        importlib_resources.files("vasp_manager")
        / "tests"
        / "calculations"
        / "NiO"
        / "rlx"
        / "POTCAR"
    )
    potcar_dst = calc_dir / "POTCAR"
    shutil.copy(str(potcar_src), str(potcar_dst))
    structure = get_pmg_structure_from_poscar(nio_poscar_path, primitive=True)
    return VaspInputCreator(
        calc_dir=calc_dir,
        mode="rlx",
        structure=structure,
        config_dir=config_dir,
        job_prefix="r",
    )


# ---------------------------------------------------------------------------
# Pure-logic: spin polarization detection
# ---------------------------------------------------------------------------


def test_check_needs_spin_polarization_no_d_f_elements(nacl_vic):
    """
    NaCl (Na, Cl) should not need spin polarization
    """
    composition_dict = nacl_vic.structure.composition.as_dict()
    result = nacl_vic._check_needs_spin_polarization(composition_dict)
    assert result is False


def test_check_needs_spin_polarization_with_d_block(nio_vic):
    """
    NiO (Ni is d-block) should need spin polarization
    """
    composition_dict = nio_vic.structure.composition.as_dict()
    result = nio_vic._check_needs_spin_polarization(composition_dict)
    assert result is True


# ---------------------------------------------------------------------------
# Pure-logic: DFT+U detection
# ---------------------------------------------------------------------------


def test_check_needs_dftu_no_hubbards(nacl_vic):
    """
    No hubbards -- no DFT+U regardless of composition
    """
    composition_dict = nacl_vic.structure.composition.as_dict()
    assert nacl_vic._check_needs_dftu(None, composition_dict) is False
    assert nacl_vic._check_needs_dftu(False, composition_dict) is False


def test_check_needs_dftu_no_oxygen(nacl_vic):
    """
    Wang DFT+U requires oxygen; NaCl has no oxygen -- no DFT+U
    """
    composition_dict = nacl_vic.structure.composition.as_dict()
    assert nacl_vic._check_needs_dftu("wang", composition_dict) is False


def test_check_needs_dftu_with_oxygen_and_d_block(nio_vic):
    """
    NiO has Ni (d-block) + O -- needs DFT+U
    """
    composition_dict = nio_vic.structure.composition.as_dict()
    assert nio_vic._check_needs_dftu("wang", composition_dict) is True


# ---------------------------------------------------------------------------
# Pure-logic: LMAXMIX
# ---------------------------------------------------------------------------


def test_get_lmaxmix_sp_elements_only(nacl_vic):
    """
    NaCl has only s/p elements -- LMAXMIX = 2
    """
    composition_dict = nacl_vic.structure.composition.as_dict()
    assert nacl_vic._get_lmaxmix(composition_dict) == 2


def test_get_lmaxmix_d_block(nio_vic):
    """
    NiO has Ni (d-block) -- LMAXMIX = 4
    """
    composition_dict = nio_vic.structure.composition.as_dict()
    assert nio_vic._get_lmaxmix(composition_dict) == 4


# ---------------------------------------------------------------------------
# Pure-logic: MAGMOM string generation
# ---------------------------------------------------------------------------


def test_get_auto_magmom_no_d_f_elements(nacl_vic):
    """
    NaCl: all elements should get 0.0 initial moment
    """
    composition_dict = nacl_vic.structure.composition.as_dict()
    magmom = nacl_vic._get_auto_magmom(composition_dict)
    assert magmom.startswith("MAGMOM =")
    assert "5.0" not in magmom
    assert "7.0" not in magmom
    # All sites should have 0.0
    assert "0.0" in magmom


def test_get_auto_magmom_d_block_element(nio_vic):
    """
    NiO: Ni (d-block) should get 5.0 initial moment
    """
    composition_dict = nio_vic.structure.composition.as_dict()
    magmom = nio_vic._get_auto_magmom(composition_dict)
    assert "5.0" in magmom


# ---------------------------------------------------------------------------
# File writing: make_poscar
# ---------------------------------------------------------------------------


def test_make_poscar_creates_file(nacl_vic):
    """
    make_poscar should write a POSCAR file
    """
    nacl_vic.make_poscar()
    poscar_path = nacl_vic.calc_dir / "POSCAR"
    assert poscar_path.exists()
    assert poscar_path.stat().st_size > 0


def test_make_poscar_valid_structure(nacl_vic):
    """
    POSCAR written by make_poscar should be readable as a pymatgen Structure
    """
    from pymatgen.core import Structure

    nacl_vic.make_poscar()
    poscar_path = nacl_vic.calc_dir / "POSCAR"
    structure = Structure.from_file(poscar_path)
    assert structure.composition.reduced_formula == "NaCl"


# ---------------------------------------------------------------------------
# File writing: make_kpoints
# ---------------------------------------------------------------------------


def test_make_kpoints_creates_file(nacl_vic):
    """
    make_kpoints should write a KPOINTS file
    """
    nacl_vic.make_kpoints()
    kpoints_path = nacl_vic.calc_dir / "KPOINTS"
    assert kpoints_path.exists()


def test_make_kpoints_format(nacl_vic):
    """
    KPOINTS file should have the correct structure (Gamma-centered grid)
    """
    nacl_vic.make_kpoints()
    kpoints_path = nacl_vic.calc_dir / "KPOINTS"
    text = kpoints_path.read_text()
    lines = text.strip().splitlines()
    assert "kpoints generated from kspacing" in lines[0]
    assert lines[1] == "0"
    assert lines[2] == "Gamma"
    # Fourth line should be three positive integers
    kpts = list(map(int, lines[3].split()))
    assert len(kpts) == 3
    assert all(k > 0 for k in kpts)


# ---------------------------------------------------------------------------
# File writing: make_vaspq
# ---------------------------------------------------------------------------


def test_make_vaspq_creates_file(nacl_vic):
    """
    make_vaspq should write a vasp.q file
    """
    nacl_vic.make_vaspq()
    vaspq_path = nacl_vic.calc_dir / "vasp.q"
    assert vaspq_path.exists()
    assert vaspq_path.stat().st_size > 0


def test_make_vaspq_is_executable(nacl_vic):
    """
    vasp.q written by make_vaspq should have execute permission
    """
    import stat

    nacl_vic.make_vaspq()
    vaspq_path = nacl_vic.calc_dir / "vasp.q"
    mode = vaspq_path.stat().st_mode
    assert mode & stat.S_IEXEC


def test_make_vaspq_jobname_contains_formula(nacl_vic):
    """
    vasp.q jobname should contain the material formula
    """
    nacl_vic.make_vaspq()
    vaspq_path = nacl_vic.calc_dir / "vasp.q"
    text = vaspq_path.read_text()
    # name is "NaCl", mode is "rlx-coarse" -- pad_string is "rc" -- jobname "rcNaCl"
    assert "rcNaCl" in text


# ---------------------------------------------------------------------------
# File writing: make_incar
# ---------------------------------------------------------------------------


def test_make_incar_creates_file(nacl_vic):
    """
    make_incar should write an INCAR file
    """
    nacl_vic.make_incar()
    incar_path = nacl_vic.calc_dir / "INCAR"
    assert incar_path.exists()
    assert incar_path.stat().st_size > 0


def test_make_incar_contains_required_tags(nacl_vic):
    """
    INCAR should contain key VASP tags
    """
    nacl_vic.make_incar()
    incar_path = nacl_vic.calc_dir / "INCAR"
    text = incar_path.read_text()
    for tag in ["ENCUT", "ISPIN", "NSW", "IBRION", "KSPACING", "LMAXMIX"]:
        assert tag in text, f"Expected tag {tag} not found in INCAR"


def test_make_incar_no_spin_for_nacl(nacl_vic):
    """
    NaCl has no d/f-block elements, so ISPIN should be 1
    """
    from vasp_manager.utils import pgrep

    nacl_vic.make_incar()
    incar_path = nacl_vic.calc_dir / "INCAR"
    ispin_line = pgrep(incar_path, "ISPIN", as_string=True)
    ispin_value = int(ispin_line.split("=")[1].strip())
    assert ispin_value == 1


def test_make_incar_spin_for_nio(nio_vic):
    """
    NiO has Ni (d-block), so ISPIN should be 2 and MAGMOM present
    """
    from vasp_manager.utils import pgrep

    nio_vic.make_incar()
    incar_path = nio_vic.calc_dir / "INCAR"
    ispin_line = pgrep(incar_path, "ISPIN", as_string=True)
    ispin_value = int(ispin_line.split("=")[1].strip())
    assert ispin_value == 2
    text = incar_path.read_text()
    assert "MAGMOM" in text


def test_make_incar_lmaxmix_always_present(nacl_vic, nio_vic):
    """
    LMAXMIX should appear in all INCARs regardless of DFT+U:
    NaCl (s/p only) -> LMAXMIX = 2, NiO (d-block) -> LMAXMIX = 4
    """
    nacl_vic.make_incar()
    assert "LMAXMIX = 2" in (nacl_vic.calc_dir / "INCAR").read_text()

    nio_vic.make_incar()
    assert "LMAXMIX = 4" in (nio_vic.calc_dir / "INCAR").read_text()
