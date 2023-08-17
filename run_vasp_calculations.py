import glob
import json
import logging
import os

from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from vasp_manager import VaspManager


def make_calculations_folder(data_path="structures.json", calcs_path="calculations"):
    """
    Make a calculations folder
        For each material, make a named folder that contains a POSCAR
    Args:
        data_path (str): path to a json with named materials and cifs
        calcs_path (str): path of base calculations folder to create
    """
    if not os.path.exists(calcs_path):
        os.mkdir(calcs_path)

    # This was my data file, but of course you can specify your own here
    with open(data_path) as fr:
        structure_dict = json.load(fr)

    for material, cif_string in structure_dict.items():
        print(material)
        material_path = os.path.join(calcs_path, material)
        if not os.path.exists(material_path):
            os.mkdir(material_path)
        structure = Structure.from_str(cif_string, fmt="cif")
        sga = SpacegroupAnalyzer(structure)
        prim_structure = sga.get_primitive_standard_structure()
        poscar = Poscar(prim_structure)
        poscar_path = os.path.join(material_path, "POSCAR")
        poscar.write_file(poscar_path)


if __name__ == "__main__":
    get_logging = False
    logging_level = logging.INFO
    if get_logging:
        logging.basicConfig()
        logging.getLogger("vasp_manager").setLevel(logging_level)

    calculation_folder = "calculations"
    if not os.path.exists(calculation_folder):
        make_calculations_folder(calcs_path=calculation_folder)
    calculation_types = [
        "rlx-coarse",
        "rlx",
        "static",
        "bulkmod",
        "elastic",
    ]
    calc_paths = os.path.join(calculation_folder, "*")
    material_paths = [p for p in sorted(glob.glob(calc_paths)) if os.path.isdir(p)]

    vmg = VaspManager(
        calculation_types, material_paths=material_paths, to_rerun=True, to_submit=True
    )
    results = vmg.run_calculations()
    print(vmg.summary())
