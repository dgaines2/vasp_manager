import json
import os
import subprocess
import warnings

import numpy as np
import pymatgen as pmg


class NumpyEncoder(json.JSONEncoder):
    """ Special json encoder for numpy types """

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def get_primitive_structure_from_poscar(path):
    structure = pmg.core.Structure.from_file(path)
    sga = pmg.symmetry.analyzer.SpacegroupAnalyzer(structure, symprec=1e-3)
    structure = sga.get_primitive_standard_structure()
    return structure


def change_elastic_constants_from_vasp(vasp_elastic_tensor):
    """
    VASP ordering of elastic constants is not the same as those expected
    by my equations below
    We should expect from Voight notation:
    1: 11 or xx
    2: 22 or yy
    3: 33 or zz
    4: 23 or yz
    5: 13 or xz
    6: 12 or xy
    but VASP differs as it presents 1, 2, 3, 6 (xy), 4 (yz), 5 (xz)
    So perform swapping to match expectations
    """
    elastic_tensor = vasp_elastic_tensor.copy()
    for j in range(6):
        elastic_tensor[3, j], elastic_tensor[4, j], elastic_tensor[5, j] = (
            elastic_tensor[4, j],
            elastic_tensor[5, j],
            elastic_tensor[3, j],
        )
    for i in range(6):
        elastic_tensor[i, 3], elastic_tensor[i, 4], elastic_tensor[i, 5] = (
            elastic_tensor[i, 4],
            elastic_tensor[i, 5],
            elastic_tensor[i, 3],
        )
    return elastic_tensor


def read_stiffness_tensor(elastic_file, change_from_vasp=True):
    with open(elastic_file, "r") as fr:
        raw_elastic_data = fr.readlines()
    # Skip first 3 rows as they are just header
    # Skip the first column as it contains xx, xy, etc
    elastic_data = [line.strip().split()[1:] for line in raw_elastic_data[3:]]
    # Divide by 10 to get GPa instead of kBar
    elastic_tensor = np.array(elastic_data, dtype=float) / 10.0
    if change_from_vasp:
        elastic_tensor = change_elastic_constants_from_vasp(elastic_tensor)
    return elastic_tensor


def get_stiffness_tensor(c11, c12, c44):
    # c is stiffness tensor
    cij = np.zeros((6, 6))
    cij[0, 0], cij[1, 1], cij[2, 2] = c11 * np.ones(3)
    cij[0, 1], cij[0, 2], cij[1, 0], cij[1, 2], cij[2, 0], cij[2, 1] = c12 * np.ones(6)
    cij[3, 3], cij[4, 4], cij[5, 5] = c44 * np.ones(3)
    return cij


def get_compliance_tensor(cij):
    sij = np.linalg.inv(cij)
    return sij


def get_B_Reuss(sij):
    B_Reuss = 1 / (
        (sij[0, 0] + sij[1, 1] + sij[2, 2]) + 2 * (sij[0, 1] + sij[1, 2] + sij[2, 0])
    )
    return B_Reuss


def get_B_Voigt(cij):
    B_Voigt = (
        (cij[0, 0] + cij[1, 1] + cij[2, 2]) + 2 * (cij[0, 1] + cij[1, 2] + cij[2, 0])
    ) / 9
    return B_Voigt


def get_G_Reuss(sij):
    G_Reuss = 15 / (
        4 * (sij[0, 0] + sij[1, 1] + sij[2, 2])
        - 4 * (sij[0, 1] + sij[1, 2] + sij[2, 0])
        + 3 * (sij[3, 3] + sij[4, 4] + sij[5, 5])
    )
    return G_Reuss


def get_G_Voigt(cij):
    G_Voigt = (
        (cij[0, 0] + cij[1, 1] + cij[2, 2])
        - (cij[0, 1] + cij[1, 2] + cij[2, 0])
        + 3 * (cij[3, 3] + cij[4, 4] + cij[5, 5])
    ) / 15
    return G_Voigt


def get_VRH_average(mod1, mod2):
    return (mod1 + mod2) / 2.0


def make_elastic_constants(outcar_path):
    if not os.path.exists(outcar_path):
        print("No OUTCAR available")
    # need to get elastic dir
    outcar_parent_dir = os.path.dirname(outcar_path)
    elastic_file = os.path.join(outcar_parent_dir, "elastic_constants.txt")
    elastic_call = f'grep -A8 "TOTAL ELASTIC MOD" {outcar_path} > {elastic_file}'
    subprocess.call(elastic_call, shell=True)


def analyze_elastic(elastic_file, to_print=True):
    try:
        cij = read_stiffness_tensor(elastic_file)
        sij = get_compliance_tensor(cij)
        B_Reuss = get_B_Reuss(sij)
        B_Voigt = get_B_Voigt(cij)
        G_Reuss = get_G_Reuss(sij)
        G_Voigt = get_G_Voigt(cij)
        B_VRH = get_VRH_average(B_Reuss, B_Voigt)
        G_VRH = get_VRH_average(G_Reuss, G_Voigt)
    except Exception as e:
        warnings.warn(f"No moduli available: {e}")
        return None

    c11 = cij[0, 0]
    c12 = cij[0, 1]
    if c11 < 0 or c12 < 0 or c11 <= 1.1 * c12:
        warning = True
    else:
        warning = False

    elastic_dict = {}
    elastic_dict["B_Reuss"] = np.round(B_Reuss, 3)
    elastic_dict["B_Voigt"] = np.round(B_Voigt, 3)
    elastic_dict["B_VRH"] = np.round(B_VRH, 3)
    elastic_dict["G_Reuss"] = np.round(G_Reuss, 3)
    elastic_dict["G_Voigt"] = np.round(G_Voigt, 3)
    elastic_dict["G_VRH"] = np.round(G_VRH, 3)
    elastic_dict["warning"] = warning
    elastic_dict["elastic_tensor"] = np.round(cij, 3)

    if to_print:
        print(json.dumps(elastic_dict, cls=NumpyEncoder, indent=2))

    return elastic_dict
