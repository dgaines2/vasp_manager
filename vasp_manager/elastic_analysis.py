import json
import logging
import os
import subprocess
import warnings

import numpy as np

from .utils import NumpyEncoder

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def change_elastic_constants_from_vasp(vasp_elastic_tensor):
    """
    VASP ordering of elastic constants is not the same as those expected
    by my equations below
    We should expect from Voigt notation:
    1: 11 or xx
    2: 22 or yy
    3: 33 or zz
    4: 23 or yz
    5: 13 or xz
    6: 12 or xy
    but VASP differs as it presents 1, 2, 3, 6 (xy), 4 (yz), 5 (xz)
    So perform swapping to match expectations

    Args:
        vasp_elastic_tensor (6x6 np.array)
    Returns:
        elastic_tensor (6x6 np.array): reordered to match Voigt notation
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
    """
    Read vasp stiffness tensor from elastic_file
    """
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
    """
    cij is stiffness tensor
    Provide alternate way to construct the tensor from c11, c12, and c44

    Args:
        c11 (float)
        c12 (float)
        c44 (float)
    Returns:
        cij (6x6 np.array)
    """
    cij = np.zeros((6, 6))
    cij[0, 0], cij[1, 1], cij[2, 2] = c11 * np.ones(3)
    cij[0, 1], cij[0, 2], cij[1, 0], cij[1, 2], cij[2, 0], cij[2, 1] = c12 * np.ones(6)
    cij[3, 3], cij[4, 4], cij[5, 5] = c44 * np.ones(3)
    return cij


def get_compliance_tensor(cij):
    """
    Args:
        cij (6x6 np.arrray): stiffness tensor
    Returns:
        sij (6x6 np.arrray): compliance tensor
    """
    sij = np.linalg.inv(cij)
    return sij


def get_B_Reuss(sij):
    """
    Args:
        sij (6x6 np.arrray): compliance tensor
    Returns:
        B_Reuss (float): Reuss bulk modulus
    """
    B_Reuss = 1 / (
        (sij[0, 0] + sij[1, 1] + sij[2, 2]) + 2 * (sij[0, 1] + sij[1, 2] + sij[2, 0])
    )
    return B_Reuss


def get_B_Voigt(cij):
    """
    Args:
        cij (6x6 np.arrray): compliance_tensor
    Returns:
        B_Reuss (float): Reuss bulk modulus
    """
    B_Voigt = (
        (cij[0, 0] + cij[1, 1] + cij[2, 2]) + 2 * (cij[0, 1] + cij[1, 2] + cij[2, 0])
    ) / 9
    return B_Voigt


def get_G_Reuss(sij):
    """
    Args:
        sij (6x6 np.arrray): compliance tensor
    Returns:
        B_Reuss (float): Reuss bulk modulus
    """
    G_Reuss = 15 / (
        4 * (sij[0, 0] + sij[1, 1] + sij[2, 2])
        - 4 * (sij[0, 1] + sij[1, 2] + sij[2, 0])
        + 3 * (sij[3, 3] + sij[4, 4] + sij[5, 5])
    )
    return G_Reuss


def get_G_Voigt(cij):
    """
    Args:
        cij (6x6 np.arrray): stiffness tensor
    Returns:
        G_Reuss (float): Reuss shear modulus
    """
    G_Voigt = (
        (cij[0, 0] + cij[1, 1] + cij[2, 2])
        - (cij[0, 1] + cij[1, 2] + cij[2, 0])
        + 3 * (cij[3, 3] + cij[4, 4] + cij[5, 5])
    ) / 15
    return G_Voigt


def get_VRH_average(mod1, mod2):
    """
    Args:
        mod1 (float):  B or G Voigt
        mode2 (float): B or G Reuss
    Returns:
        VRH_average (float)
    """
    return (mod1 + mod2) / 2.0


def make_elastic_constants(outcar_path):
    """
    Utility function to scrape OUTCAR for elastic constants
    Writes to elastic_constants.txt
    """
    if not os.path.exists(outcar_path):
        raise Exception("No OUTCAR available")
    # need to get elastic dir
    outcar_parent_dir = os.path.dirname(outcar_path)
    elastic_file = os.path.join(outcar_parent_dir, "elastic_constants.txt")
    elastic_call = f'grep -A8 "TOTAL ELASTIC MOD" {outcar_path} > {elastic_file}'
    subprocess.call(elastic_call, shell=True)


def analyze_elastic_file(elastic_file, verbose=True):
    """
    Grab important quantities from the elastic calculation results

    Args:
        elastic_file (str): filepath
        to_print (bool): if True, print the result
    """
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
        logger.warnings(f"No moduli available: {e}")
        return None

    c11 = cij[0, 0]
    c12 = cij[0, 1]
    if c11 < 0 or c12 < 0 or c11 <= 1.1 * c12:
        logger.warning("-" * 10 + " WARNING: Elastically Unstable " + "-" * 10)
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

    if verbose:
        logging.info(json.dumps(elastic_dict, cls=NumpyEncoder, indent=2))

    return elastic_dict
