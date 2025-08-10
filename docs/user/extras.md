# Extras

* To manually stop `VaspManager` from processing a material, place a `STOP` file
in that material's directory: e.g. `calculations/NaCl/STOP`.

* The module logger is also made available for information and debugging and can
be accessed through `logging.getLogger("vasp_manager")`.

* The current implementation has only been tested on Linux and Mac OS.

* At this point, KPOINT generation is handled through the `KSPACING` tag in the
  `INCAR`, but future versions will be able to specify KPPRA or a manual grid
  instead. Spin-orbit coupling calculations are also not currently supported.

\\\ TODO: Implement `band-structure` calculations and possibly `phonopy`
calculations
