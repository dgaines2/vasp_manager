# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import shutil

import importlib_resources
import pytest


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
def base_calcs_dir(tmp_path_factory):
    """
    Calculations directory that needs each run type to be set up
    (no rlx*, static, bulkmod, or elastic subdirectories)
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
        ignore=shutil.ignore_patterns("rlx*", "static", "bulkmod", "elastic"),
    )
    return new_calculations_dir


@pytest.fixture(scope="session")
def stable_material_dir(tmp_path_factory):
    """
    NaCl material directory with a complete, elastically stable calculation
    """
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
    """
    CuO material directory with a complete, elastically unstable calculation
    """
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
