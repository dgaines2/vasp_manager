# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import shutil

import importlib_resources
import pytest

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
        / "material"
        / "POSCAR"
    )


@pytest.fixture(scope="session")
def cuo_poscar_path():
    """
    Path to CuO POSCAR from test data (has d-block elements)
    """
    return (
        importlib_resources.files("vasp_manager")
        / "tests"
        / "calculations"
        / "material_spinu"
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
    calc_dir = tmp_path / "material" / "rlx-coarse"
    calc_dir.mkdir(parents=True)
    # copy POTCAR so make_incar can read it
    potcar_src = (
        importlib_resources.files("vasp_manager")
        / "tests"
        / "calculations"
        / "material"
        / "rlx-coarse"
        / "POTCAR"
    )
    potcar_dst = calc_dir / "POTCAR"
    shutil.copy(str(potcar_src), str(potcar_dst))
    return VaspInputCreator(
        calc_dir=calc_dir,
        mode="rlx-coarse",
        poscar_source_path=nacl_poscar_path,
        config_dir=config_dir,
    )


@pytest.fixture
def cuo_vic(tmp_path, cuo_poscar_path, config_dir):
    """
    VaspInputCreator for CuO, rlx mode
    """
    calc_dir = tmp_path / "material_spinu" / "rlx"
    calc_dir.mkdir(parents=True)
    potcar_src = (
        importlib_resources.files("vasp_manager")
        / "tests"
        / "calculations"
        / "material_spinu"
        / "rlx"
        / "POTCAR"
    )
    potcar_dst = calc_dir / "POTCAR"
    shutil.copy(str(potcar_src), str(potcar_dst))
    return VaspInputCreator(
        calc_dir=calc_dir,
        mode="rlx",
        poscar_source_path=cuo_poscar_path,
        config_dir=config_dir,
    )


# ---------------------------------------------------------------------------
# Pure-logic: spin polarization detection
# ---------------------------------------------------------------------------


def test_check_needs_spin_polarization_no_d_f_elements(nacl_vic):
    """
    NaCl (Na, Cl) should not need spin polarization
    """
    composition_dict = nacl_vic.source_structure.composition.as_dict()
    result = nacl_vic._check_needs_spin_polarization(composition_dict)
    assert result is False


def test_check_needs_spin_polarization_with_d_block(cuo_vic):
    """
    CuO (Cu is d-block) should need spin polarization
    """
    composition_dict = cuo_vic.source_structure.composition.as_dict()
    result = cuo_vic._check_needs_spin_polarization(composition_dict)
    assert result is True


# ---------------------------------------------------------------------------
# Pure-logic: DFT+U detection
# ---------------------------------------------------------------------------


def test_check_needs_dftu_no_hubbards(nacl_vic):
    """
    No hubbards → no DFT+U regardless of composition
    """
    composition_dict = nacl_vic.source_structure.composition.as_dict()
    assert nacl_vic._check_needs_dftu(None, composition_dict) is False
    assert nacl_vic._check_needs_dftu(False, composition_dict) is False


def test_check_needs_dftu_no_oxygen(nacl_vic):
    """
    Wang DFT+U requires oxygen; NaCl has no oxygen → no DFT+U
    """
    composition_dict = nacl_vic.source_structure.composition.as_dict()
    assert nacl_vic._check_needs_dftu("wang", composition_dict) is False


def test_check_needs_dftu_with_oxygen_and_d_block(cuo_vic):
    """
    CuO has Cu (d-block) + O → needs DFT+U
    """
    composition_dict = cuo_vic.source_structure.composition.as_dict()
    assert cuo_vic._check_needs_dftu("wang", composition_dict) is True


# ---------------------------------------------------------------------------
# Pure-logic: LMAXMIX
# ---------------------------------------------------------------------------


def test_get_lmaxmix_sp_elements_only(nacl_vic):
    """
    NaCl has only s/p elements → LMAXMIX = 2
    """
    composition_dict = nacl_vic.source_structure.composition.as_dict()
    assert nacl_vic._get_lmaxmix(composition_dict) == 2


def test_get_lmaxmix_d_block(cuo_vic):
    """
    CuO has Cu (d-block) → LMAXMIX = 4
    """
    composition_dict = cuo_vic.source_structure.composition.as_dict()
    assert cuo_vic._get_lmaxmix(composition_dict) == 4


# ---------------------------------------------------------------------------
# Pure-logic: MAGMOM string generation
# ---------------------------------------------------------------------------


def test_get_auto_magmom_no_d_f_elements(nacl_vic):
    """
    NaCl: all elements should get 0.0 initial moment
    """
    composition_dict = nacl_vic.source_structure.composition.as_dict()
    magmom = nacl_vic._get_auto_magmom(composition_dict)
    assert magmom.startswith("MAGMOM =")
    assert "5.0" not in magmom
    assert "7.0" not in magmom
    # All sites should have 0.0
    assert "0.0" in magmom


def test_get_auto_magmom_d_block_element(cuo_vic):
    """
    CuO: Cu (d-block) should get 5.0 initial moment
    """
    composition_dict = cuo_vic.source_structure.composition.as_dict()
    magmom = cuo_vic._get_auto_magmom(composition_dict)
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
    # name is "NaCl", mode is "rlx-coarse" → pad_string is "rc" → jobname "rcNaCl"
    assert "rcNaCl" in text


# ---------------------------------------------------------------------------
# File writing: make_incar
# ---------------------------------------------------------------------------


def _make_fake_potcar(element_nelectrons: list[tuple[str, float]]):
    """
    Return a fake Potcar-like list of objects with .element and .nelectrons.
    The POTCAR files in the test fixtures are minimal stubs (header lines only),
    so we monkeypatch Potcar.from_file to avoid pymatgen parsing failures.
    """
    from types import SimpleNamespace

    return [SimpleNamespace(element=el, nelectrons=ne) for el, ne in element_nelectrons]


def test_make_incar_creates_file(nacl_vic, monkeypatch):
    """
    make_incar should write an INCAR file
    """
    monkeypatch.setattr(
        "vasp_manager.vasp_input_creator.Potcar.from_file",
        lambda path: _make_fake_potcar([("Na", 7.0), ("Cl", 7.0)]),
    )
    nacl_vic.make_incar()
    incar_path = nacl_vic.calc_dir / "INCAR"
    assert incar_path.exists()
    assert incar_path.stat().st_size > 0


def test_make_incar_contains_required_tags(nacl_vic, monkeypatch):
    """
    INCAR should contain key VASP tags
    """
    monkeypatch.setattr(
        "vasp_manager.vasp_input_creator.Potcar.from_file",
        lambda path: _make_fake_potcar([("Na", 7.0), ("Cl", 7.0)]),
    )
    nacl_vic.make_incar()
    incar_path = nacl_vic.calc_dir / "INCAR"
    text = incar_path.read_text()
    for tag in ["ENCUT", "ISPIN", "NSW", "IBRION", "KSPACING"]:
        assert tag in text, f"Expected tag {tag} not found in INCAR"


def test_make_incar_no_spin_for_nacl(nacl_vic, monkeypatch):
    """
    NaCl has no d/f-block elements, so ISPIN should be 1
    """
    from vasp_manager.utils import pgrep

    monkeypatch.setattr(
        "vasp_manager.vasp_input_creator.Potcar.from_file",
        lambda path: _make_fake_potcar([("Na", 7.0), ("Cl", 7.0)]),
    )
    nacl_vic.make_incar()
    incar_path = nacl_vic.calc_dir / "INCAR"
    ispin_line = pgrep(incar_path, "ISPIN", as_string=True)
    ispin_value = int(ispin_line.split("=")[1].strip())
    assert ispin_value == 1


def test_make_incar_spin_for_cuo(cuo_vic, monkeypatch):
    """
    CuO has Cu (d-block), so ISPIN should be 2 and MAGMOM present
    """
    from vasp_manager.utils import pgrep

    monkeypatch.setattr(
        "vasp_manager.vasp_input_creator.Potcar.from_file",
        lambda path: _make_fake_potcar([("Cu", 11.0), ("O", 6.0)]),
    )
    cuo_vic.make_incar()
    incar_path = cuo_vic.calc_dir / "INCAR"
    ispin_line = pgrep(incar_path, "ISPIN", as_string=True)
    ispin_value = int(ispin_line.split("=")[1].strip())
    assert ispin_value == 2
    text = incar_path.read_text()
    assert "MAGMOM" in text
