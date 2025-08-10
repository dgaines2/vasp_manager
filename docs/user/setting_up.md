# Setting Up

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

For more information about `calc_config.json` and `computing_config.json`, please see

When you're done, your calculation directory should look roughly like this:
```mermaid
graph TD
  A[calculations/] --> B([calc_config.json]) & C([computing_config.json])
  A[calculations/] --> D[Material 1] --> E([POSCAR])
  A[calculations/] --> F[Material 2] --> G([POSCAR])
  A[calculations/] --> H[Material ...] --> I([POSCAR])
```
