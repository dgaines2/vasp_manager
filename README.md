# vasp_manager
Automatically run vasp relaxation and bulk moduli calculations

This package serves to automate `VASP` calculations. We include calculation
modes `"rlx-coarse"`, `"rlx-fine"`, `"bulkmod"`, `"bulkmod_rlx"`, and
`"elastic"`.  Each mode has its own configuration settings
`vasp_manager/config/calc_config.json` with sensible defaults, but these can be
easily customized by the user.

## Calculation Modes
`rlx-coarse`: (optional) lower precision energy-based relaxation

`rlx-fine`: tighter force-based relaxation

`bulkmod_standalone`: standalone (no relaxation required) bulk modulus
calculation using an Equation of State (EOS) fit to an energy-volume curve

`bulkmod`: bulk modulus calculation using the EOS but starting from the
rlx-fine output structure

`elastic`: Determination of elastic constants using the deformation method
built into VASP

\\\ TODO: Implement `static` calculations, `band-structure` calculations, and
possibly `phonopy` calculations

I generally recommend starting from `rlx-coarse`, although the functionality is
there to start a `rlx-fine` from the initially provided POSCAR.

Example workflows might look like `rlx-coarse` &#8594; `rlx` &#8594; `bulkmod`, or
`rlx` &#8594; `elastic`.

## User Info
In order to use this package, you MUST

1) Create a calculations folder. Each subfolder of `calculations/` should have a
unique name and contain a `POSCAR`. A sample method of creating the calculations
folder from a `pandas.DataFrame` is available in `run_vasp_calculations.py`.
2) Configure `vasp_manager/config/computing_config.json`. You will need to
specify your `user_id`, a `potcar directory`, a `queuetype`, and a `vasp
module`. As of now, only CORI at NERSC and QUEST from Northwestern University
are supported. Any other SLURM based supercomputers can be easily added,
but modifications could be made for other queue management systems.

Vasp input creation is automatic, and so is job submission and analysis. Simply
rerun the main script and any calculations that are ready for the next type of
calculation will be created and submitted.  The bulk moduli analysis is carried
out in the backend using the open-source `Pymatgen` software, and elastic
constant analysis is through custom personal scripts.

The main object for handling all calculations is `vasp_manager.VaspManager`,
which takes in a list of calculation modes. See the class documentation for more
details.

The module logger is also made available for information and  debugging and can
be accessed through `logging.getLogger("vasp_manager")`.

## Notes

The current implementation has only been tested on Linux and Mac OS, and it relies
on the unix utilities `cat` and `grep`.

At this point, KPOINT generation is handled through the KSPACING
tag in the INCAR, but future versions will be able to specify KPPRA or a manual
grid instead. Magnetic or spin-orbit coupling calculations are also not currently
supported.*
