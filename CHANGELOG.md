# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.4] - 2024-01-17

### Added

- Added ability to run static calculations without a previously existing rlx calculation
- Added capability to override job preambles/commands by placing a {computer}.yml file in the calculations folder
- Added capability to override the global calc\_config.json by placing a new calc\_config.json in a material's calculation mode folder.

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
