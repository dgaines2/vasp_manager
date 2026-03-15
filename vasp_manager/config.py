import json
from functools import lru_cache
from pathlib import Path

from pydantic import BaseModel, ConfigDict


class ComputingConfig(BaseModel):
    """Typed model for a single computer's computing configuration."""

    computer: str
    user_id: str
    potcar_dir: Path
    queuetype: str
    allocation: str
    constraint: str
    vasp_module: str
    ncore: int
    ncore_per_node: int


class CalcConfig(BaseModel):
    """Typed model for a single calculation mode's configuration.

    Known fields are typed; any additional keys (custom INCAR tags) are
    accepted via extra='allow' and included in model_dump().

    Fields accepting float | str (ediff, ediffg, symprec, sigma) allow
    Fortran-compatible scientific notation strings (e.g. "1e-07") to be
    passed through as-is rather than being reformatted by Python.
    """

    model_config = ConfigDict(extra="allow")

    # INCAR fields
    prec: str
    ispin: str | int  # "auto" | 1 | 2
    kspacing: float
    symprec: float | str
    nsw: int
    ibrion: int
    isif: int
    lreal: str
    potim: float
    ediffg: float | str
    iopt: int
    nelm: int
    encut: int
    ediff: float | str
    algo: str
    ismear: int
    sigma: float | str
    kpar: int
    gga: str

    # Write flags — Python bools; converted to ".TRUE."/".FALSE." at INCAR write time
    lcharge: bool = False
    lwave: bool = False
    lvtot: bool = False

    # Meta fields (not INCAR tags)
    hubbards: str | bool | None
    walltime: str

    # Optional elastic-only fields
    write_kpoints: bool = False
    nfree: int | None = None


@lru_cache(maxsize=None)
def load_computing_config(config_dir: Path) -> ComputingConfig:
    """Load and cache computing_config.json for the given config directory.

    Flattens the nested structure by extracting the active computer's block
    and populating ComputingConfig with a validated model.
    """
    fname = "computing_config.json"
    fpath = config_dir / fname
    if not fpath.exists():
        raise FileNotFoundError(f"No {fname} found in {fpath.absolute()}")
    with open(fpath) as f:
        raw = json.load(f)
    computer = raw["computer"]
    return ComputingConfig(computer=computer, **raw[computer])


@lru_cache(maxsize=None)
def load_calc_configs(config_dir: Path) -> dict[str, CalcConfig]:
    """Load and cache calc_config.json for the given config directory.

    Returns a dict keyed by mode name (e.g. "rlx", "static").
    """
    fname = "calc_config.json"
    fpath = config_dir / fname
    if not fpath.exists():
        raise FileNotFoundError(f"No {fname} found in {fpath.absolute()}")
    with open(fpath) as f:
        raw = json.load(f)
    return {mode: CalcConfig.model_validate(cfg) for mode, cfg in raw.items()}
