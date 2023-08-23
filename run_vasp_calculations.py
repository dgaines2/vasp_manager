import json
import logging
from pathlib import Path

from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from vasp_manager import VaspManager


def make_calculations_folder(data_path="structures.json", calcs_path="calculations"):
    """
    Make a calculations folder
        For each material, make a named folder that contains a POSCAR
    Args:
        data_path (str | Path): path to a json with named materials and cifs
        calcs_path (str | Path): path of base calculations folder to create
    """
    calcs_path = Path(calcs_path)
    if not calcs_path.exists():
        calcs_path.mkdir()

    # This was my data file, but of course you can specify your own here
    with open(data_path) as fr:
        structure_dict = json.load(fr)

    for material, cif_string in structure_dict.items():
        print(material)
        material_path = calcs_path / material
        if not material_path.exists():
            material_path.mkdir()
        structure = Structure.from_str(cif_string, fmt="cif")
        sga = SpacegroupAnalyzer(structure)
        prim_structure = sga.get_primitive_standard_structure()
        poscar = Poscar(prim_structure)
        poscar_path = material_path / "POSCAR"
        poscar.write_file(poscar_path)


if __name__ == "__main__":
    get_logging = False
    logging_level = logging.INFO
    if get_logging:
        logging.basicConfig()
        logging.getLogger("vasp_manager").setLevel(logging_level)

    calculation_folder = Path("calculations")
    if not calculation_folder.exists():
        make_calculations_folder(calcs_path=calculation_folder)

    calc_config_path = calculation_folder / "calc_config.json"
    computing_config_path = calculation_folder / "computing_config.json"
    if not calc_config_path.exists() or not computing_config_path.exists():
        raise Exception(
            f"""Couldn't find one of the configuration files
           Be sure to set up the configuration files in {calculation_folder}
           before trying to run VaspManager"""
        )

    calculation_types = [
        "rlx-coarse",
        "rlx",
        "static",
        "bulkmod",
        "elastic",
    ]
    material_paths = [p for p in sorted(calculation_folder.glob("*")) if p.is_dir()]

    vmg = VaspManager(
        calculation_types,
        material_paths=material_paths,
        to_rerun=True,
        to_submit=True,
    )
    results = vmg.run_calculations()
    print(vmg.summary())
