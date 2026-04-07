# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [2.0.1] 2026-04-07

### Changed

- `LMAXMIX` is now written to all INCARs (not only DFT+U calculations); value is determined by element type (2 for s/p, 4 for d-block, 6 for f-block)

### Fixed

- Fix `STOPPED` status cascading to downstream calc types — dependency check now runs before the stopped check, so calcs with unsatisfied deps are silently skipped instead of also writing `STOPPED`


## [2.0.0] 2026-03-17

### Added

- Add typed pydantic config models (`CalcConfig`, `ComputingConfig`) with cached loaders
- Add `VaspRun` class that owns run-level concerns (error checking, magmom parsing)
- Add configurable job resource parameters (`atoms_per_node`, `rerun_increase`, `rerun_increase_factor`) to `ComputingConfig`
- Add smart rerun diagnosis: relaxations distinguish "completed all NSW without converging" (relaunch same params) from "timed out" (apply rerun strategy)

### Changed

- Add `develop` branch to CI pull request trigger
- Replace hardcoded calculation ordering with a dependency graph (`CALC_DEPENDENCIES` + topological sort)
- `VaspInputCreator` and `JobManager` now use shared cached config loaders instead of reading JSON on each instantiation
- Migrate from black/flake8/isort to ruff
- Update comprehensive test suite: `test_vasp_input_creator.py`, `test_job_manager.py`, `test_config.py`, and additional analyzer failure mode tests
- `VaspInputCreator` now accepts a pymatgen `Structure` directly instead of a `poscar_source_path`
- `BaseCalculationManager` generalized with `vasp_runs` dict to support both single-run and multi-run managers
- Job naming (`pad_string`) extracted from `VaspInputCreator` into each `CalculationManager` via `job_prefix`
- Move `job_prefix` abstract property to `BaseCalculationManager` for SLURM job naming
- Bulkmod strains are now independent jobs, each with their own `VaspRun`, `VaspInputCreator`, and `JobManager`
- Bulkmod handles per-strain failures individually instead of restarting from scratch
- Update example calculations from VASP 6.3.2 (CuO, NaCl) to VASP 6.4.3 (NaCl, Fe, NiO, BCC\_Ti)
- `make_potcar_anonymous` now produces minimal pymatgen-parseable POTCARs instead of bare TITEL-line stubs
- Update test fixtures to 4 materials covering spin-polarization, DFT+U, re-relaxation, and elastic instability
- Migrate CI workflows (tests, publish, mkdocs) from pip to uv

### Removed

- Remove bulkmod-specific SLURM templates (`vasp-bulkmod.yml`, `vasp-bulkmod-quest.yml`, `vasp-bulkmod-bridges2.yml`); bulkmod strains now reuse static templates
- Remove `bulkmod` entry from `calc_config.json`; bulkmod now uses the `static` calc config


## [1.4.2] - 2025-11-13

### Fixed

- Fix hubbard settings to only apply to oxygen-containing compounds


## [1.4.1] - 2025-10-09

### Changed

- Updated pre-commit version

### Fixed

- Fixed walltime formatting for walltimes longer than 24 hours


## [1.4.0] - 2025-08-12

### Added

- Write documentation and deploy as gh-pages
- Add type hints and type checking

### Changed

- README is simplified and points to the documentation site


## [1.3.5] - 2025-08-10

### Changed

- Enforce calculation types ordering in VaspManager


## [1.3.4] - 2025-08-03

### Changed

- Change POSCAR precision from 8 to 16 to avoid symmetry errors associated with floating point errors


## [1.3.3] - 2025-04-02

### Changed

- Use SPDX license expression for compliance with PEP 639


## [1.3.2] - 2025-04-02

### Changed

- Use dynamic versioning in pyproject.toml


## [1.3.1] - 2025-04-01

### Added

- Added publish to PyPi github workflow


## [1.3.0] - 2025-03-24

### Added

- Added optional mixing tags (AMIX, BMIX) to INCAR from calc config

### Changed

- Updated quest configuration for quest[10,11,12] nodes


## [1.2.0] - 2025-03-23

### Added

- Added oqmd psuedopotentials json, oqmd-settings branch now uses these psuedopotentials by default
- Added primitive to be passed in calculation manager kwargs
- Allow custom calculation configs for individual job managers
- Added sort\_by callable for sorting results by keys
- Allow other filenames for slurm exe and jobid files
- Added (optional) write tags LCHARG, LWAVE, and LVTOT to calc config
- Catch sbatch errors with subprocess.check\_output()
- Added ASCII logo to VaspManager init
- Added international\_monoclinic argument to get\_pmg\_structure\_from\_poscar

### Changed

- Changed vasp.q permissions to executable executable
- Added newline to end of jobid file
- Updated property setters in Analyzers
- Updated github actions versions
- [BREAKING] ElasticAnalyzer is now a better standalone analyzer, and all of the associated processing of VASP outputs is now in the from\_calc\_dir method. This makes it easier to use ElasticAnalyzer as a utility if you only have the elastic constants and a structure
- Updated managers with named loggers. Each manager should now print the material name when logging, which should make things easier to track when using multiprocessing
- Stopped pinning versions, track pymatgen master instead
- [BREAKING] Unify path naming conventions. Use \_dir for directories and \_path for filepaths for clarity

### Fixed

- Fixed rlx restart behavior if to\_rerun is False
- Fixed spacing in VASP error handling message
- Fixed github actions coverage
- Fixed static post-hoc analysis
- Fixed stopping behavior for better status tracking


## [1.1.4] - 2024-01-17

### Added

- Added ability to run static calculations without a previously existing rlx calculation
- Added capability to override job preambles/commands by placing a {computer}.yml file in the calculations folder
- Added capability to override the global calc\_config.json by placing a new calc\_config.json in a material's calculation mode folder

### Changed

- Major update to error handling, including tracking of STOPPED calculations
- Changed number of cores reserved for memory when hitting out-of-memory errors as well as for new Quest 9+ nodes
- Changed default VASP settings for elastic calculations in calc\_config.json. Convergence testing w.r.t. smearing and KPOINTS is still recommended
- Allow VaspManager to recognize previously zipped archives (for rlx-coarse or rlx calculations)
- Removed bulkmod\_standalone as a mode and incorporated the behavior into the normal bulkmod manager

### Fixed

- Ensure elastic calculations use a conventional unit cell


## [1.1.3] - 2023-10-19

### Added

- Added the ability to write a KPOINTS file by specifying write\_kpoints: true in calc\_config.json

### Changed

- Major updates to vasp error handling
- Changed default symprec in calc\_config.json
- When a job fails, double the walltime instead of doubling the number of nodes. This is typically more efficient in terms of total computational cost
- Made variable names in util functions more explicit

### Fixed

- Fixed error in increasing the walltime upon job failure
- Made the mpitasks variable uniform in the vasp.q yml files, and fix the improper reference to mpi\_tasks
- Fixed POTCAR path in computing\_config.json
- Fixed error in ElasticManager when the number of displacements >= 100
- Fixed Quest computing configuration after Quest8 nodes were retired
- Fixed automatic NBANDS for ElasticManager when the number of available compute cores is greater than 3/4 the number of electrons
- Fixed BulkmodManager symlink creation when the symlink already exists


## [1.1.2] - 2023-10-06

### Added

- Added support for bridges2 in computing\_config.json and vasp.q yml templates

### Changed

- Use only pyproject.toml instead of setup.py and pyproject.toml
- Changed LICENSE to .md instead of .txt
- Changed utility functions (phead, ptail, pgrep) to iterate over lines in the file rather than reading the entire file into memory
- pot\_dict.json is now a symbolic link to enable better support for multiple VASP versions
- Utilize yml tempaltes for vasp.q creation for easier customizability
