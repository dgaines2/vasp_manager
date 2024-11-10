import json
import logging
from pathlib import Path

from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from vasp_manager import VaspManager


def make_calculations_dir(data_path="structures.json", calcs_dir="calculations"):
    """
    Make a calculations directory
        For each material, make a named directory that contains a POSCAR
    Args:
        data_path (str | Path): path to a json with named materials and cifs
        calcs_dir (str | Path): path of base calculations directory to create
    """
    calcs_dir = Path(calcs_dir)
    if not calcs_dir.exists():
        calcs_dir.mkdir()

    # This was my data file, but of course you can specify your own here
    with open(data_path) as fr:
        structure_dict = json.load(fr)

    for material, cif_string in structure_dict.items():
        print(material)
        material_dir = calcs_dir / material
        if not material_dir.exists():
            material_dir.mkdir()
        structure = Structure.from_str(cif_string, fmt="cif")
        sga = SpacegroupAnalyzer(structure)
        prim_structure = sga.get_primitive_standard_structure()
        poscar = Poscar(prim_structure)
        poscar_path = material_dir / "POSCAR"
        poscar.write_file(poscar_path)


if __name__ == "__main__":
    get_logging = False
    logging_level = logging.INFO
    if get_logging:
        logging.basicConfig()
        logging.getLogger("vasp_manager").setLevel(logging_level)

    calculations_dir = Path("calculations")
    if not calculations_dir.exists():
        make_calculations_dir(calcs_dir=calculations_dir)

    calc_config_path = calculations_dir / "calc_config.json"
    computing_config_path = calculations_dir / "computing_config.json"
    if not calc_config_path.exists() or not computing_config_path.exists():
        raise Exception(
            f"""Couldn't find one of the configuration files
           Be sure to set up the configuration files in {calculations_dir}
           before trying to run VaspManager"""
        )

    calculation_types = [
        "rlx-coarse",
        "rlx",
        "static",
        "bulkmod",
        "elastic",
    ]
    material_dirs = [p for p in sorted(calculations_dir.glob("*")) if p.is_dir()]

    vmg = VaspManager(
        calculation_types,
        material_dirs=material_dirs,
        to_rerun=True,
        to_submit=True,
    )
    results = vmg.run_calculations()
    print(vmg.summary())
