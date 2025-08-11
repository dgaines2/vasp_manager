# Extras

* To manually stop `VaspManager` from processing a material, place a `STOP` file
in that material's directory: e.g. `calculations/NaCl/STOP`.
* The module logger is also made available for information and debugging and can
be accessed through `logging.getLogger("vasp_manager")`.
* The current implementation has only been tested on Linux and Mac OS.
* At this point, KPOINT generation is handled through the `KSPACING` tag in the
  `INCAR`, but future versions will be able to specify KPPRA or a manual grid
  instead. Spin-orbit coupling calculations are also not currently supported.
* On github, there's an `oqmd-settings` branch that roughly approximates the
  relaxation and static settings of the [OQMD](https://oqmd.org). But, do note
  that if you use this branch to try to compare with the OQMD, you should also
  use the corresponding potpaw.52 pseudopotentials.
