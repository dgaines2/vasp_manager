# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import json
from pathlib import Path

import importlib_resources
import pytest

from vasp_manager.config import (
    CalcConfig,
    ComputingConfig,
    load_calc_configs,
    load_computing_config,
)


@pytest.fixture(scope="session")
def config_dir():
    """
    Path to the shared test config directory (has calc_config.json, computing_config.json)
    """
    return Path(importlib_resources.files("vasp_manager") / "tests" / "calculations")


@pytest.fixture(autouse=True)
def clear_config_cache():
    """
    Clear lru_cache before and after each test to prevent cross-test cache pollution
    """
    load_computing_config.cache_clear()
    load_calc_configs.cache_clear()
    yield
    load_computing_config.cache_clear()
    load_calc_configs.cache_clear()


def _make_computing_config_json(tmp_path, computer="personal"):
    """
    Write a minimal but fully valid computing_config.json to tmp_path.
    Returns tmp_path for use as config_dir.
    """
    config = {
        "computer": computer,
        computer: {
            "user_id": "testuser",
            "potcar_dir": "/fake/potcars",
            "queuetype": "regular",
            "allocation": "fake_alloc",
            "constraint": "cpu",
            "vasp_module": "vasp/6.0",
            "ncore": 16,
            "ncore_per_node": 128,
        },
    }
    (tmp_path / "computing_config.json").write_text(json.dumps(config))
    return tmp_path


# ---------------------------------------------------------------------------
# load_computing_config
# ---------------------------------------------------------------------------


def test_load_computing_config_returns_computing_config(config_dir):
    """
    load_computing_config returns a ComputingConfig instance
    """
    cc = load_computing_config(config_dir)
    assert isinstance(cc, ComputingConfig)


def test_load_computing_config_computer_matches_top_level_key(config_dir):
    """
    cc.computer is populated from the top-level "computer" key in the JSON,
    not from within the nested block (verifies flattening logic)
    """
    raw = json.loads((config_dir / "computing_config.json").read_text())
    cc = load_computing_config(config_dir)
    assert cc.computer == raw["computer"]


def test_load_computing_config_potcar_dir_is_path(config_dir):
    """
    potcar_dir is coerced from str to Path by pydantic
    """
    cc = load_computing_config(config_dir)
    assert isinstance(cc.potcar_dir, Path)


def test_load_computing_config_int_fields(config_dir):
    """
    ncore and ncore_per_node are ints
    """
    cc = load_computing_config(config_dir)
    assert isinstance(cc.ncore, int)
    assert isinstance(cc.ncore_per_node, int)


def test_load_computing_config_missing_file_raises(tmp_path):
    """
    FileNotFoundError when computing_config.json does not exist in the given directory
    """
    with pytest.raises(FileNotFoundError):
        load_computing_config(tmp_path)


def test_load_computing_config_is_cached(config_dir):
    """
    Two calls with the same config_dir return the same object
    """
    cc1 = load_computing_config(config_dir)
    cc2 = load_computing_config(config_dir)
    assert cc1 is cc2


# ---------------------------------------------------------------------------
# load_calc_configs
# ---------------------------------------------------------------------------


def test_load_calc_configs_keys_match_json_modes(config_dir):
    """
    Returned dict keys match the top-level keys of the raw calc_config.json
    """
    raw = json.loads((config_dir / "calc_config.json").read_text())
    configs = load_calc_configs(config_dir)
    assert set(configs.keys()) == set(raw.keys())


def test_load_calc_configs_values_are_calc_config(config_dir):
    """
    Every value in the returned dict is a CalcConfig instance
    """
    configs = load_calc_configs(config_dir)
    for cfg in configs.values():
        assert isinstance(cfg, CalcConfig)


def test_load_calc_configs_kspacing_is_float(config_dir):
    """
    kspacing is a float across all loaded modes
    """
    configs = load_calc_configs(config_dir)
    for mode, cfg in configs.items():
        assert isinstance(cfg.kspacing, float), f"kspacing not float for mode {mode!r}"


def test_load_calc_configs_missing_file_raises(tmp_path):
    """
    FileNotFoundError when calc_config.json does not exist in the given directory
    """
    with pytest.raises(FileNotFoundError):
        load_calc_configs(tmp_path)


def test_load_calc_configs_is_cached(config_dir):
    """
    Two calls with the same config_dir return the same object
    """
    c1 = load_calc_configs(config_dir)
    c2 = load_calc_configs(config_dir)
    assert c1 is c2


# ---------------------------------------------------------------------------
# CalcConfig model
# ---------------------------------------------------------------------------


def test_calc_config_bool_flag_defaults(config_dir):
    """
    lcharge, lwave, and lvtot default to False on a freshly loaded config
    """
    configs = load_calc_configs(config_dir)
    for mode, cfg in configs.items():
        assert cfg.lcharge is False, f"lcharge not False for mode {mode!r}"
        assert cfg.lwave is False, f"lwave not False for mode {mode!r}"
        assert cfg.lvtot is False, f"lvtot not False for mode {mode!r}"


def test_calc_config_ediff_formatting(config_dir):
    """
    A string ediff value is coerced to float by pydantic
    """
    configs = load_calc_configs(config_dir)
    base_data = next(iter(configs.values())).model_dump()
    cfg = CalcConfig.model_validate({**base_data, "ediff": "1e-07"})
    assert isinstance(cfg.ediff, float)
    assert cfg.ediff == 1e-07


def test_calc_config_ediff_accepts_float(config_dir):
    """
    CalcConfig also accepts a float for ediff and preserves its type
    """
    configs = load_calc_configs(config_dir)
    base_data = next(iter(configs.values())).model_dump()
    cfg = CalcConfig.model_validate({**base_data, "ediff": 1e-5})
    assert isinstance(cfg.ediff, float)


def test_calc_config_extra_fields_accepted(config_dir):
    """
    A custom INCAR tag added to CalcConfig is stored in model_extra
    """
    configs = load_calc_configs(config_dir)
    base_data = next(iter(configs.values())).model_dump()
    cfg = CalcConfig.model_validate({**base_data, "NELMIN": 6})
    assert cfg.model_extra.get("NELMIN") == 6


def test_calc_config_elastic_optional_fields_present(config_dir):
    """
    Elastic config has write_kpoints as bool and nfree as int or None
    """
    configs = load_calc_configs(config_dir)
    if "elastic" not in configs:
        pytest.skip("No elastic mode in test calc_config.json")
    cfg = configs["elastic"]
    assert isinstance(cfg.write_kpoints, bool)
    assert cfg.nfree is None or isinstance(cfg.nfree, int)


def test_calc_config_custom_override_applied(config_dir):
    """
    A custom key in the override dict wins over the base config value.
    Uses a sentinel kspacing (0.99) that won't appear in real configs.
    """
    configs = load_calc_configs(config_dir)
    base = configs["rlx-coarse"]
    merged = CalcConfig.model_validate({**base.model_dump(), "kspacing": 0.99})
    assert merged.kspacing == 0.99
