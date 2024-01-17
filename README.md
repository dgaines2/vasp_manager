```text
 __      __             __  __
 \ \    / /            |  \/  |
  \ \  / /_ _ ___ _ __ | \  / | __ _ _ __   __ _  __ _  ___ _ __
   \ \/ / _` / __| '_ \| |\/| |/ _` | '_ \ / _` |/ _` |/ _ \ '__|
    \  / (_| \__ \ |_) | |  | | (_| | | | | (_| | (_| |  __/ |
     \/ \__,_|___/ .__/|_|  |_|\__,_|_| |_|\__,_|\__, |\___|_|
                 | |                              __/ |
                 |_|                             |___/ v1.1.4
```

<hr/>

<h4 align="center">

[![Tests](https://github.com/dgaines2/vasp_manager/actions/workflows/tests.yml/badge.svg)](https://github.com/dgaines2/vasp_manager/actions/workflows/tests.yml)
[![codecov](https://codecov.io/github/dgaines2/vasp_manager/graph/badge.svg?token=CQH3BRGYCR)](https://codecov.io/github/dgaines2/vasp_manager)
[![PyPI Version](https://img.shields.io/pypi/v/vasp-manager)](https://pypi.org/project/vasp-manager)
[![Requires Python 3.10+](https://img.shields.io/badge/python-3.10+-blue)](https://python.org/downloads)

Automatically run `VASP` relaxation, static, bulk moduli, or elastic constant
calculations

</h4>

## How to Install

1. Create a new environment with python version $\geq$ 3.10
2. Clone this repository
3. Run `pip install -e .` or optionally `pip install -e .[dev]` to include
packages needed for development/contribution

This package is also available on
[PyPi](https://pypi.org/project/vasp-manager/#description). To install, run
`pip install vasp-manager`.

## User Guide

This package serves to automate `VASP` calculations. `VASP` input creation is
automatic, and so is job submission, queue monitoring, calculation analysis, and
storage of the results. Simply rerun the main script and any calculations that
are ready for the next type of calculation will be created and submitted.

### VaspManager

The main class for handling all calculations is `vasp_manager.VaspManager`,
which takes in a list of calculation types and material paths. See the class
documentation for more details. By default, results are exported to
`calculations/results.json`.

The bulk moduli analysis is carried out in the backend using the open-source
`pymatgen` software to fit an EOS and elastic constant analysis using custom
scripts.

### Calculation Modes

We include calculation modes `"rlx-coarse"`, `"rlx"`, `"static"`, `"bulkmod"`,
and `"elastic"`.  The desired modes to calculate are specified when
initializing a `VaspManager` object.

* `rlx-coarse`: lower precision energy-based relaxation
* `rlx`: tighter force-based relaxation
* `static`: high accuracy static SCF calculation
* `bulkmod`: bulk modulus calculation using an Equation of State (EOS) fit to an
energy-volume curve
  * Can be run as a standalone (no relaxation required) calculation of bulk
    modulus using an EOS (although this is not recommended unless you are sure
    the cell volume is very close to the equilibrium value)
* `elastic`: Determination of elastic constants using the strain/deformation
method built into `VASP`

I generally recommend starting from `rlx-coarse`, although the functionality is
there to start a `rlx` calculation from the initially provided POSCAR.

Most users' workflows follow `rlx-coarse` &#8594; `rlx` &#8594; `static`. The
modes `static`, `bulkmod`, and `elastic` can all be run independently of each
other. For example, workflows might look like `rlx-coarse` &#8594; `rlx` &#8594;
`static` &#8594; `bulkmod`, or `rlx` &#8594; `elastic`, or simply `bulkmod`.
The `elastic` mode requires at least `rlx` preceding it in order to guarantee
converged lattice parameters and atomic positions.

## Usage Guide

In order to use this package, you MUST

1. Create a calculations folder where you'd like to run your calculations.  Each
subfolder of `calculations/` should have a unique name and contain a `POSCAR`. A
sample method of creating the calculations folder from a `json` with names and
cifs is available in `run_vasp_calculations.py`, and an example calculations
folder is provided in `calculations/`.
2. Configure `computing_config.json` and place it in the `calculations/`
directory.  You will need to specify your `user_id`, a `potcar_directory`, a
`queuetype`, your `allocation` and a `vasp_module` (VASP 6 strongly
recommended). As of now, only Perlmutter and Bridges2 at NERSC and QUEST at
Northwestern University are supported. Any other SLURM based supercomputers can
be easily added, but modifications could be made for other queue management
systems.
3. If desired, make modifications to `calc_config.json`. This must also be
placed in the `calculations/` directory. Each mode has its own configuration
settings with sensible defaults, but these can be easily customized by the user.
    * To include spin polarization, set `"ispin": "auto"` in
    `calc_config.json`; otherwise set `"ispin": 1`. With this setting, all
    elements with valence *d* or *f* electrons will start with initial magnetic
    moments of 5 and 7 $\mu_B$, respectively. `VaspManager` also accepts an
    additional argument `magmom_per_atom_cutoff` which defaults to 0. If this
    argument is passed, `rlx` calculations that finish with a magmom per atom
    less than this value with be re-run without spin polarization. This argument
    only affects `rlx` calculations, and the spin setting for following
    `static`, `bulkmod`, or `elastic` calculations is inferred from the final
    `rlx` calculation.
    * To include DFT+U for transition metal oxides, set `"hubbards": "wang"`;
    otherwise, set `"hubbards": null`.  Currently, only `"gga": "PE"` (PBE) is
    supported.

To manually stop `VaspManager` from processing a material, place a `STOP` file
in that material's directory: e.g. `calculations/NaCl/STOP`.

The module logger is also made available for information and debugging and can
be accessed through `logging.getLogger("vasp_manager")`.

### Notes

* *The current implementation has only been tested on Linux and Mac OS.*
* *At this point, KPOINT generation is handled through the KSPACING tag in the
  INCAR, but future versions will be able to specify KPPRA or a manual grid
  instead. Spin-orbit coupling calculations are also not currently supported.*
* *For those using Quest, I recommend using you own pre-compiled version of
 `VASP` rather than the Quest `VASP` module. For those in the Wolverton Group,
  consider using the `quest-hotfix` branch.*

\\\ TODO: Implement `band-structure` calculations and possibly `phonopy`
calculations
