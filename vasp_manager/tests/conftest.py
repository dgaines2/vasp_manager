# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import shutil

import importlib_resources
import pytest


def pytest_configure(config):
    """Suppress UnknownPotcarWarning from pymatgen for minimal test POTCARs."""
    config.addinivalue_line(
        "filterwarnings",
        "ignore::pymatgen.io.vasp.inputs.UnknownPotcarWarning",
    )


@pytest.fixture(scope="session")
def calcs_dir(tmp_path_factory):
    """
    Complete calculations directory (all calc types finished)
    """
    calculations_dir = (
        importlib_resources.files("vasp_manager") / "tests" / "calculations"
    )
    new_calculations_dir = tmp_path_factory.mktemp("calculations")
    shutil.copytree(
        calculations_dir,
        new_calculations_dir,
        dirs_exist_ok=True,
        symlinks=True,
    )
    return new_calculations_dir


@pytest.fixture(scope="session")
def stable_material_dir(tmp_path_factory):
    """
    NaCl material directory with a complete, elastically stable calculation
    """
    material_dir = (
        importlib_resources.files("vasp_manager") / "tests" / "calculations" / "NaCl"
    )
    new_material_dir = tmp_path_factory.mktemp("NaCl")
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
    """
    BCC_Ti material directory with a complete, elastically unstable calculation
    """
    material_dir = (
        importlib_resources.files("vasp_manager") / "tests" / "calculations" / "BCC_Ti"
    )
    new_material_dir = tmp_path_factory.mktemp("BCC_Ti")
    shutil.copytree(
        material_dir,
        new_material_dir,
        dirs_exist_ok=True,
        symlinks=True,
        ignore=shutil.ignore_patterns("elastic_constants.txt"),
    )
    return new_material_dir
