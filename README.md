# vasp_manager
Automatically run vasp relaxation and bulk moduli calculations

This package serves to automate `VASP` calculations. We include calculation modes `"rlx-coarse"`, `"rlx-fine"`, `"bulkmod"`, `"bulkmod_rlx"`, and `"elastic"`.
Each mode has its own configuration settings `vasp_manager/config/calc_config.json` with sensible defaults, but these can be easily customized by the user. 

`rlx-coarse`: (optional) lower precision energy-based relaxation  
`rlx-fine`: tighter force-based relaxation  
`bulkmod`: standalone (no relaxation required) bulk modulus calculation using an Equation of State (EOS) fit to an energy-volume curve  
`bulkmod_rlx`: bulk modulus calculation using the EOS but starting from the rlx-fine output structure  
`elastic`: Determination of elastic constants using the deformation method built into VASP  

The elastic analysis is carried out in the backend using the open-source `Pymatgen` software.   

In order to use this package, you MUST
1) Create a calculations folder. Each subfolder of `calculations/` should have a unique name and contain a `POSCAR`. A sample method of creating the calculations folder from a `pandas.DataFrame` is available in `run_vasp_calculations.py`.
2) Configure `vasp_manager/config/computing_config.json`. You will need to specify your `user_id`, a `potcar directory`, a `queuetype`, and a `vasp module`. As of now, only CORI at NERSC and QUEST from Northwestern University are supported.  

The main function is `vasp_manager.calculation_manager.manage_calculations()`

*Note:
I plan to package everything into a class soon. Additionally, At this point, KPOINT generation is handled through the KSPACING tag in the INCAR, but future versions will be able to specify KPPRA or a manual grid instead. Magnetic or spin-orbit coupling calculations are not currently supported.*
