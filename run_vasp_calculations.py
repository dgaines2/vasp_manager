import glob
import logging
import os
import sys

import pandas as pd
from pymatgen.io.vasp import Poscar
from vasp_manager import manage_calculations


def make_calculations_folder(data_path="structure_df.pickle.gz"):
    """
    Make a calculations folder
        For each material, make a folder named {oqmd_id} that contains a
        POSCAR
    Args:
        data_path (str): path to a pandas DF with oqmd_ids and pmg structures
    """
    calcs_path = "calculations"
    if not os.path.exists(calcs_path):
        os.mkdir(calcs_path)

    # This was my data file, but of course you can specify your own here
    df = pd.read_pickle(data_path)
    oqmd_ids = df["oqmd_id"].values
    structures = df["structure"].values
    for oqmd_id, structure in zip(oqmd_ids, structures):
        print(oqmd_id)
        oqmd_id_path = os.path.join(calcs_path, oqmd_id)
        if not os.path.exists(oqmd_id_path):
            os.mkdir(oqmd_id_path)
        poscar_path = os.path.join(oqmd_id_path, "POSCAR")
        poscar = Poscar(structure)
        poscar.write_file(poscar_path)


if __name__ == "__main__":
    get_logging = True
    logging_level = logging.DEBUG
    if get_logging:
        logging.basicConfig()
        logging.getLogger("vasp_manager").setLevel(logging_level)

    if not os.path.exists("calculations"):
        make_calculations_folder()
    calculation_types = ["rlx-coarse", "rlx-fine", "elastic"]

    manage_calculations(calculation_types)
