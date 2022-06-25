# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import glob
import json
import logging
import os
import pkgutil
import shutil
import subprocess

import numpy as np
import pymatgen as pmg
from pymatgen.analysis.eos import BirchMurnaghan

from .elastic_analysis import analyze_elastic_file, make_elastic_constants
from .utils import change_directory

logger = logging.getLogger(__name__)

computing_config_dict = json.loads(
    pkgutil.get_data("vasp_manager", "config/computing_config.json").decode("utf-8")
)
calc_config_dict = json.loads(
    pkgutil.get_data("vasp_manager", "config/calc_config.json").decode("utf-8")
)
q_mapper = json.loads(
    pkgutil.get_data("vasp_manager", "static_files/q_handles.json").decode("utf-8")
)
potcar_dict = json.loads(
    pkgutil.get_data("vasp_manager", "static_files/pot_dict.json").decode("utf-8")
)
incar_template = pkgutil.get_data("vasp_manager", "static_files/INCAR_template").decode(
    "utf-8"
)
computer = computing_config_dict["computer"]


def make_potcar(structure, write_path):
    """
    Create and write a POTCAR given a pmg structure

    Args:
        structure (pmg.Structure)
        write_path (str): destination path
    """
    potcar_dir = computing_config_dict[computer]["potcar_dir"]

    el_names = [el.name for el in structure.composition]
    logger.debug(f"{structure.composition.reduced_formula}, {el_names}")
    pot_singles = [
        os.path.join(potcar_dir, potcar_dict[el_name], "POTCAR") for el_name in el_names
    ]
    cmd = "cat " + " ".join(pot_singles) + " > " + write_path
    subprocess.call(cmd, shell=True)


def make_incar(incar_path, mode):
    """
    mode should be 'rlx-coarse', 'rlx', 'bulkmod', 'bulkmod_rlx', or 'elastic'

    Need to modify this to account for spin/magmom
    Current kpoints coming from the kspacing tag in the INCAR,
        but future versions should include ability to make kpoints from kppra
    """
    ncore = computing_config_dict[computer]["ncore"]

    # bulkmod and bulkmod_rlx have same calc configs
    if "bulkmod" in mode:
        mode = "bulkmod"
    calc_config = calc_config_dict[mode]
    # if calc_config_dict["magmom"] == "auto":
    #     calc_config_dict["magmom"] = ""

    # Add lines to the vaspq file for only elastic calculations
    incar_tmp = incar_template
    if "elastic" in mode:
        incar_tmp = incar_tmp.split("\n")
        for i, line in enumerate(incar_tmp):
            if "KSPACING" in line:
                nfree_line = "NFREE = {nfree}"
                symprec_line = "SYMPREC = {symprec}"
                incar_tmp.insert(i + 1, symprec_line)
                incar_tmp.insert(i + 1, nfree_line)
            if "NCORE" in line:
                incar_tmp[i] = "NPAR = 1"
        incar_tmp = "\n".join([line for line in incar_tmp])
    incar = incar_tmp.format(**calc_config, ncore=ncore)
    logger.debug(incar)
    with open(incar_path, "w+") as fw:
        fw.write(incar)


def make_vaspq(vaspq_path, mode, jobname=None, structure=None, increase_nodes=False):
    """
    mode should be 'rlx-coarse', 'rlx', 'bulkmod', 'bulkmod_rlx', or 'elastic'
    """
    calc_config = calc_config_dict[mode]
    walltime = calc_config["walltime"]

    parent_dir = os.path.dirname(vaspq_path)
    if structure is None:
        try:
            poscar_path = os.path.join(parent_dir, "POSCAR")
            structure = pmg.core.Structure.from_file(poscar_path)
        except Exception as e:
            raise Exception(f"Cannot load POSCAR in {parent_dir}: {e}")

    # start with 1 node per 32 atoms
    n_nodes = (len(structure) // 32) + 1
    # create pad string for job naming to differentiate in the queue
    if mode == "rlx" or mode == "rlx-coarse":
        if mode == "rlx":
            pad_string = "r"
        elif mode == "rlx-coarse":
            pad_string = "rc"
        mode = "rlx"
    elif "bulkmod" in mode:
        pad_string = "b"
        mode = "bulkmod"
    elif mode == "elastic":
        pad_string = "e"

    if computer == "quest":
        # quest has small nodes
        n_nodes *= 2
        n_procs = n_nodes * 28
    else:
        n_procs = n_nodes * 68

    if jobname is None:
        jobname = pad_string + structure.composition.reduced_formula
    else:
        jobname = pad_string + jobname

    if increase_nodes:
        n_nodes *= 2
        n_procs *= 2
        hours, minutes, seconds = walltime.split(":")
        hours = str(int(hours) * 2)
        walltime = ":".join([hours, minutes, seconds])

    computer_config = computing_config_dict[computer].copy()
    computer_config.update({"n_nodes": n_nodes, "n_procs": n_procs, "jobname": jobname})

    q_name = q_mapper[computer][mode]
    vaspq_tmp = pkgutil.get_data("vasp_manager", f"static_files/{q_name}").decode(
        "utf-8"
    )
    vaspq = vaspq_tmp.format(**computer_config, walltime=walltime)
    logger.debug(vaspq)
    with open(vaspq_path, "w+") as fw:
        fw.write(vaspq)


def make_archive(compound_path, mode):
    """
    Make an archive of a VASP calculation and copy back over relevant files

    Args:
        compound_path (str): base path of calculation
        mode (str): subdirectory of compound_path to run the calculation

    Returns:
        archive_made (bool): if job was never submitted, return False
    """
    mode_path = os.path.join(compound_path, mode)
    oqmd_id = compound_path.split("/")[1]

    # if job was never submitted, don't make an archive
    jobid_path = os.join(mode_path, "jobid")
    if not os.path.exists(jobid_path):
        return False

    with change_directory(mode_path):
        num_archives = len(glob.glob("archive*"))
        all_files = [d for d in glob.glob("*") if os.path.isfile(d)]
        archive_name = f"archive_{num_archives}"
        logger.info(f"Making {archive_name}...")
        os.mkdir(archive_name)
        for f in all_files:
            shutil.move(f, archive_name)

        contcar_path = os.path.join(archive_name, "CONTCAR")
        shutil.copy(contcar_path, "POSCAR")
        potcar_path = os.path.join(archive_name, "POTCAR")
        shutil.copy(potcar_path, "POTCAR")

        incar_path = os.path.join(mode_path, "INCAR")
        vaspq_path = os.path.join(mode_path, "vasp.q")
        make_incar(incar_path, mode=mode)
        make_vaspq(vaspq_path, mode=mode, jobname=oqmd_id)
    return True


def make_bulkmod_strains(bulkmod_path, strains):
    """
    Create a set of strain directory for fitting the E-V info

    Args:
        bulkmod_path (str)
        strains (iterable of floats)
    """
    logger.info("Making strain directories")
    for i, strain in enumerate(strains):
        middle = int(len(strains) / 2)
        strain_index = i - middle
        strain_name = f"strain_{strain_index}"
        strain_path = os.path.join(bulkmod_path, strain_name)
        logger.info(strain_path)

        if not os.path.exists(strain_path):
            os.mkdir(strain_path)
        orig_poscar_path = os.path.join(bulkmod_path, "POSCAR")
        strain_poscar_path = os.path.join(strain_path, "POSCAR")
        shutil.copy(orig_poscar_path, strain_poscar_path)

        # change second line to be {strain} rather than 1.0
        with open(strain_poscar_path, "r") as fr:
            strain_poscar = fr.read().splitlines()
        strain_poscar[1] = f"{strain}"
        logger.debug(f"{strain}")
        with open(strain_poscar_path, "w+") as fw:
            fw.write("\n".join(strain_poscar))

        with change_directory(strain_path):
            for f in ["POTCAR", "INCAR"]:
                orig_path = "../" + f
                f_path = f
                os.symlink(orig_path, f_path, target_is_directory=False)


def analyze_bulkmod(bulkmod_path, from_relax=False):
    """
    Fit an EOS to calculate the bulk modulus from a finished bulkmod calculation

    Args:
        bulkmod_path (str): bulkmod path of the calculation
        from_relax (bool): if True, copy the CONTCAR from the relaxation folder
            If False, copy the POSCAR from compound_path
    """
    strain_paths = [
        path for path in glob.glob(bulkmod_path + "/strain*") if os.path.isdir(path)
    ]
    strain_paths = sorted(strain_paths, key=lambda d: int(d.split("_")[-1]))
    volumes = []
    final_energies = []
    for i, strain_path in enumerate(strain_paths):
        poscar_path = os.path.join(strain_path, "POSCAR")
        vasprun_path = os.path.join(strain_path, "vasprun.xml")
        volume = pmg.core.Structure.from_file(poscar_path).volume
        vasprun = pmg.io.vasp.outputs.Vasprun(
            filename=vasprun_path,
            parse_dos=False,
            parse_eigen=False,
            parse_potcar_file=False,
        )
        final_energy = vasprun.final_energy
        volumes.append(volume)
        final_energies.append(final_energy)
    logger.debug("Volumes")
    logger.debug(f"{volumes}")
    logger.debug("Final Energies")
    logger.debug(f"{final_energies}")
    eos_analyzer = BirchMurnaghan(volumes, final_energies)
    eos_analyzer.fit()
    bulk_modulus = np.round(eos_analyzer.b0_GPa, 3)
    logger.info(f"{mode.upper()} Calculation: Successful")
    logger.info(f"BULK MODULUS: {bulk_modulus}")
    return bulk_modulus


def analyze_elastic(elastic_path):
    """
    Get results from elastic calculation
    """
    elastic_file = os.path.join(elastic_path, "elastic_constants.txt")
    if not os.path.exists(elastic_file):
        outcar_file = os.path.join(elastic_path, "OUTCAR")
        make_elastic_constants(outcar_file)
    results = analyze_elastic_file(elastic_file)
    return results
