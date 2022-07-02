# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import json
import os
from contextlib import contextmanager

import numpy as np
import pymatgen as pmg


@contextmanager
def change_directory(new_dir):
    prev_dir = os.getcwd()
    os.chdir(os.path.expanduser(new_dir))
    try:
        yield
    finally:
        os.chdir(prev_dir)


class NumpyEncoder(json.JSONEncoder):
    """Special json encoder for numpy types"""

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def get_pmg_structure_from_poscar(
    poscar_path,
    to_process=True,
    primitive=True,
    symprec=1e-3,
    return_sg=False,
):
    """
    Args:
        poscar_path (str)
        to_process (bool): if True, get standard reduced structure
        primitive (bool): if True, get primitive structure, else get
        conventional structure
        symprec (float): symprec for SpacegroupAnalyzer
        return_sg (bool): if True, return spacegroup number
    Returns:
        structure (pmg.Structure): structure from POSCAR
    """
    structure = pmg.core.Structure.from_file(poscar_path)
    if to_process:
        sga = pmg.symmetry.analyzer.SpacegroupAnalyzer(structure, symprec=symprec)
        if primitive:
            structure = sga.get_primitive_standard_structure()
        else:
            structure = sga.get_conventional_standard_structure()
        if return_sg:
            sg = sga.get_space_group_number()
            return structure, sg
    return structure
