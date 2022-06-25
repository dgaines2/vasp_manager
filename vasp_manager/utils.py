# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import json
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


def get_primitive_structure_from_poscar(poscar_path):
    """
    Args:
        poscar_path (str)
    Returns:
        structure (pmg.Structure): primitive structure from POSCAR
    """
    structure = pmg.core.Structure.from_file(poscar_path)
    sga = pmg.symmetry.analyzer.SpacegroupAnalyzer(structure, symprec=1e-3)
    structure = sga.get_primitive_standard_structure()
    return structure
