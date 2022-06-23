import glob
import json
import logging
import os
import shutil
import subprocess
import sys
import warnings
from contextlib import contextmanager

import numpy as np
import pandas as pd
import pymatgen as pmg
from pymatgen.analysis.eos import BirchMurnaghan
from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

global logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()


@contextmanager
def change_directory(new_dir):
    prev_dir = os.getcwd()
    os.chdir(os.path.expanduser(new_dir))
    try:
        yield
    finally:
        os.chdir(prev_dir)


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


def get_computing_config_dict(
    computing_config_dict_path="static_files/computing_config.json",
):
    """
    config dict info about the current computer
    """
    with open(computing_config_dict_path, "r") as fr:
        computing_config_dict = json.load(fr)
    return computing_config_dict


def get_calc_config_dict(
    calc_config_dict_path="static_files/calc_config.json", mode=None
):
    """
    config dict info for calculations
    """
    with open(calc_config_dict_path, "r") as fr:
        calc_config_dict = json.load(fr)
    if mode is None:
        return calc_config_dict
    else:
        return calc_config_dict[mode]


def get_potcar_dict(potcar_dict_path="static_files/pot_dict.json"):
    """
    potcar_dict should contain the name of the psuedopotential you want to
    use for each element
    """
    with open(potcar_dict_path, "r") as fr:
        potcar_dict = json.load(fr)
    return potcar_dict


def make_potcar(structure, write_path):
    """
    Create and write a POTCAR given a pmg structure

    Args:
        structure (pmg.Structure)
        write_path (str): destination path
    """
    potcar_dict = get_potcar_dict()
    computing_config_dict = get_computing_config_dict()
    computer = computing_config_dict["computer"]
    potcar_dir = computing_config_dict[computer]["potcar_dir"]

    el_names = [el.name for el in structure.composition]
    logging.info(el_names)
    pot_singles = [
        os.path.join(potcar_dir, potcar_dict[el_name], "POTCAR") for el_name in el_names
    ]
    cmd = "cat " + " ".join(pot_singles) + " > " + write_path
    subprocess.call(cmd, shell=True)


def make_incar(incar_path, mode):
    """
    mode should be 'rlx-coarse', 'rlx', or 'bulkmod'

    Need to modify this to account for spin/magmom
    Current kpoints coming from the kspacing tag in the INCAR,
        but future versions should include ability to make kpoints from kppra
    """
    computing_config_dict = get_computing_config_dict()
    computer = computing_config_dict["computer"]
    ncore = computing_config_dict[computer]["ncore"]

    if "rlx" in mode:
        mode = "relax"
    elif "bulkmod" in mode:
        mode = "bulkmod"
    calc_config_dict = get_calc_config_dict(mode=mode)
    # if calc_config_dict["magmom"] == "auto":
    #     calc_config_dict["magmom"] = ""

    source_incar_path = f"static_files/INCAR_template"
    with open(source_incar_path, "r") as fr:
        incar_tmp = fr.read()
    incar = incar_tmp.format(**calc_config_dict, ncore=ncore)
    with open(incar_path, "w+") as fw:
        fw.write(incar)


def make_vaspq(vaspq_path, mode, jobname=None, structure=None):
    """
    mode should be 'rlx-coarse', 'rlx', or 'bulkmod'
    """
    computing_config_dict = get_computing_config_dict()
    computer = computing_config_dict["computer"]
    computer_config_dict = computing_config_dict[computer]

    calc_config_dict = get_calc_config_dict(mode=mode)
    walltime = calc_config_dict["walltime"]

    parent_dir = os.path.dirname(vaspq_path)
    if structure is None:
        try:
            poscar_path = os.path.join(parent_dir, "POSCAR")
            structure = pmg.core.Structure.from_file(poscar_path)
        except Exception as e:
            raise Exception(f"Cannot load POSCAR in {parent_dir}")

    n_nodes = (len(structure) // 32) + 1
    if "rlx" in mode:
        mode = "rlx"
    if computer == "quest":
        source_vaspq_path = f"static_files/vasp-{mode}-quest.q"
        # quest has small nodes
        n_nodes *= 2
        n_procs = n_nodes * 28
    else:
        source_vaspq_path = f"static_files/vasp-{mode}-cori.q"
        n_procs = n_nodes * 68
    if jobname is None:
        jobname = structure.composition.reduced_formula
    computer_config_dict.update(
        {"n_nodes": n_nodes, "n_procs": n_procs, "jobname": jobname}
    )

    with open(source_vaspq_path, "r") as fr:
        vaspq_tmp = fr.read()
    vaspq = vaspq_tmp.format(**computer_config_dict, walltime=walltime)
    with open(vaspq_path, "w+") as fw:
        fw.write(vaspq)


def submit_job(compound_path, mode, ignore_errors=False):
    computing_config_dict = get_computing_config_dict()
    computer = computing_config_dict["computer"]
    if "personal" in computer:
        ignore_errors = True
        error_msg = "Cannot submit job on personal computer"
        error_msg += "\n\tIgnoring job submission..."
        logger.info(error_msg)
        return True

    calc_path = os.path.join(compound_path, mode)
    jobid_path = os.path.join(calc_path, "jobid")
    if os.path.exists(jobid_path):
        logger.info("Job already exists")
        return
    submission_call = "sbatch vasp.q | awk '{ print $4 }' | tee jobid"
    with change_directory(calc_path):
        subprocess.call(submission_call, shell=True)


def check_job_complete(jobid, ignore_errors=False):
    """Return True if job done"""
    computing_config_dict = get_computing_config_dict()
    computer = computing_config_dict["computer"]
    if computer == "personal":
        ignore_errors = True
        error_msg = "Cannot check job on personal computer"
        error_msg += "\n\tIgnoring job status check..."
        logger.info(error_msg)
        return True
    else:
        user_id = computing_config_dict[computer]["user_id"]
        check_queue_call = f"squeue -u {user_id}"
        queue_call = (
            subprocess.check_output(check_queue_call, shell=True)
            .decode("utf-8")
            .splitlines()
        )
        for line in queue_call:
            line = line.strip().split()
            if jobid in line:
                return False
        return True


def make_archive(compound_path, mode):
    """
    Make an archive of a VASP calculation and copy back over relevant files

    Args:
        compound_path (str): base path of calculation
        mode (str): subdirectory of compound_path to run the calculation
    """
    mode_path = os.path.join(compound_path, mode)
    oqmd_id = compound_path.split("/")[1]
    with change_directory(mode_path):
        num_archives = len(glob.glob("archive*"))
        all_files = [d for d in glob.glob("*") if os.path.isfile(d)]
        archive_name = f"archive_{num_archives}"
        logging.info(f"Making {archive_name}...")
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


def setup_coarse_relax(
    compound_path,
    submit=False,
    from_scratch=False,
    rerun_relax=False,
):
    """
    Set up a coarse relaxation

    Args:
        compound_path (str): base path of calculation
        submit (bool): if True, submit calculation
        from_scratch (bool): if True, remove previous calculation
        rerun_relax (bool): if True, make an archive and resubmit
        to_change_incar (bool): if True, call make_incar_subsitutions
    """
    computing_config_dict = get_computing_config_dict()
    # POSCAR, POTCAR, INCAR, vasp.q
    oqmd_id = compound_path.split("/")[1]
    mode = "rlx-coarse"
    crelax_path = os.path.join(compound_path, mode)
    if from_scratch:
        if os.path.exists(crelax_path):
            shutil.rmtree(crelax_path)
    if not os.path.exists(crelax_path):
        os.mkdir(crelax_path)

    # POSCAR
    if rerun_relax:
        make_archive(compound_path, mode=mode)
    else:
        orig_poscar_path = os.path.join(compound_path, "POSCAR")
        final_poscar_path = os.path.join(crelax_path, "POSCAR")
        shutil.copy(orig_poscar_path, final_poscar_path)

        # POTCAR
        structure = pmg.core.Structure.from_file(final_poscar_path)
        potcar_path = os.path.join(crelax_path, "POTCAR")
        make_potcar(structure, potcar_path)

        # INCAR
        incar_path = os.path.join(crelax_path, "INCAR")
        make_incar(incar_path, mode=mode)

        # vasp.q
        vaspq_path = os.path.join(crelax_path, "vasp.q")
        make_vaspq(vaspq_path, mode=mode, jobname=oqmd_id)

    if submit:
        submit_job(compound_path, mode=mode)


def setup_relax(
    compound_path,
    submit=False,
    from_scratch=False,
    rerun_relax=False,
    from_coarse=False,
    change_incar=False,
):
    """
    Set up a fine relaxation

    Args:
        compound_path (str): base path of calculation
        submit (bool): if True, submit calculation
        from_scratch (bool): if True, remove previous calculation
        rerun_relax (bool): if True, make an archive and if submit, resubmit
        from_coarse (bool): if True, copy the CONTCAR from the coarse relaxation
            folder
            If False, copy the POSCAR from compound_path
        to_change_incar (bool): if True, run make_incar_substitutions()
    """
    computing_config_dict = get_computing_config_dict()
    mode = "rlx"
    if from_coarse:
        crelax_done = check_relax(compound_path, mode="rlx-coarse")
        if not crelax_done:
            logger.info("Coarse relax not finished")
            return

    oqmd_id = compound_path.split("/")[1]
    relax_path = os.path.join(compound_path, mode)
    if from_scratch:
        if os.path.exists(relax_path):
            shutil.rmtree(relax_path)
    if not os.path.exists(relax_path):
        os.mkdir(relax_path)

    if rerun_relax:
        make_archive(compound_path, mode=mode)
    else:
        # POSCAR, POTCAR, INCAR, vasp.q
        if from_coarse:
            # POSCAR
            crelax_path = os.path.join(compound_path, "rlx-coarse")
            orig_poscar_path = os.path.join(crelax_path, "CONTCAR")
            final_poscar_path = os.path.join(relax_path, "POSCAR")
            shutil.copy(orig_poscar_path, final_poscar_path)
            structure = pmg.core.Structure.from_file(final_poscar_path)

            # POTCAR
            orig_potcar_path = os.path.join(crelax_path, "POTCAR")
            final_potcar_path = os.path.join(relax_path, "POTCAR")
            shutil.copy(orig_potcar_path, final_potcar_path)
        else:
            # POSCAR
            orig_poscar_path = os.path.join(compound_path, "POSCAR")
            final_poscar_path = os.path.join(relax_path, "POSCAR")
            shutil.copy(orig_poscar_path, final_poscar_path)
            structure = pmg.core.Structure.from_file(final_poscar_path)

            # POTCAR
            potcar_path = os.path.join(relax_path, "POTCAR")
            make_potcar(structure, potcar_path)

        # INCAR
        incar_path = os.path.join(relax_path, "INCAR")
        make_incar(incar_path, mode=mode)

        # vasp.q
        vaspq_path = os.path.join(relax_path, "vasp.q")
        make_vaspq(vaspq_path, mode=mode, jobname=oqmd_id)

    if submit:
        submit_job(compound_path, mode=mode)


def check_relax(
    compound_path,
    rerun_relax=False,
    submit=False,
    tail=5,
    mode="rlx-coarse",
):
    """
    Check if calculation has finished and reached required accuracy
    (No real automatic logging or fixing of VASP errors)

    Args:
        compound_path (str): base path of calculation
        rerun_relax (bool): if True and job fails, make an archive and re-setup
            the rlx folder
        submit (bool): if True and rerun_relax, submit calculation
        tail (int): If job fails, print {tail} lines from stdout.txt
        mode (str): "rlx" or "rlx-coarse"
    Returns:
        relaxation_successful (bool): if True, relaxation completed successfully
    """
    relax_path = os.path.join(compound_path, mode)
    stdout_path = os.path.join(relax_path, "stdout.txt")
    jobdone_path = os.path.join(relax_path, "jobdone")

    if os.path.exists(stdout_path):
        jobid_path = os.path.join(relax_path, "jobid")
        # read jobid from file
        with open(jobid_path, "r") as fr:
            job_id = fr.readline().splitlines()[0]
        job_complete = check_job_complete(job_id)
        if not job_complete:
            logger.info(f"{mode} not finished")
            return False

        grep_call = f"tail -n{tail} {stdout_path}"
        grep_output = (
            subprocess.check_output(grep_call, shell=True).decode("utf-8").strip()
        )
        if "reached required accuracy" in grep_output:
            logger.info(f"{mode} reached required accuracy")
            return True
        else:
            archive_dirs = glob.glob(f"{relax_path}/archive*")
            if "coarse" in mode and len(archive_dirs) >= 3:
                logger.info("Many archives exist, suggest force based relaxation")
                if rerun_relax:
                    setup_relax(compound_path, submit=submit, rerun_relax=True)
                return True

            logger.info("FAILED")
            logger.debug(grep_output)
            if rerun_relax:
                logger.info(f"Rerunning {compound_path}")
                if "coarse" in mode:
                    setup_coarse_relax(compound_path, submit=submit, rerun_relax=True)
                else:
                    setup_relax(compound_path, submit=submit, rerun_relax=True)
            return False
    else:
        logger.info(f"{mode} not started")
        return False


def check_volume_difference(compound_path, rerun_relax=False, submit=False):
    """
    check relaxation runs for volume difference
    if abs(volume difference) is >= 5%, rerun relaxation

    Args:
        compound_path (str): base path of calculation
        rerun_relax (bool): if True, make an archive and re-setup the rlx folder
        submit (bool): if True and rerun_relax, submit calculation
    """
    relaxation_done = check_relax(
        compound_path, rerun_relax=rerun_relax, submit=submit, mode="rlx"
    )
    if not relaxation_done:
        logger.info(f"  Relaxation not finished")
        return False

    relax_path = os.path.join(compound_path, "rlx")
    poscar_path = os.path.join(relax_path, "POSCAR")
    contcar_path = os.path.join(relax_path, "CONTCAR")
    try:
        p_structure = pmg.core.Structure.from_file(poscar_path)
        c_structure = pmg.core.Structure.from_file(contcar_path)
    except Exception:
        logger.info(f"  Relaxation CONTCAR doesn't exist or is empty")
        return False

    sga1 = SpacegroupAnalyzer(p_structure, symprec=1e-3)
    sga2 = SpacegroupAnalyzer(c_structure, symprec=1e-3)
    p_spacegroup = sga1.get_space_group_number()
    c_spacegroup = sga2.get_space_group_number()
    if p_spacegroup == c_spacegroup:
        logger.info(f"  Spacegroups match {p_spacegroup}=={c_spacegroup}")
    else:
        logger.warning(
            f"   Warning: spacegroups do not match {p_spacegroup} != {c_spacegroup}"
        )

    volume_diff = (c_structure.volume - p_structure.volume) / p_structure.volume
    if np.abs(volume_diff) >= 0.05:
        logger.warning(f"  NEED TO RE-RELAX: dV = {volume_diff:.4f}")
        volume_converged = False
        if rerun_relax:
            setup_relax(compound_path, rerun_relax=True, submit=submit)
    else:
        logger.info(f"  dV = {volume_diff:.4f}")
        volume_converged = True
    return volume_converged


def setup_simple_bulkmod(
    compound_path,
    submit=True,
    from_scratch=False,
    from_relax=False,
):
    """
    Set up EOS bulkmod calculation

    Args:
        compound_path (str): base path of calculation
        submit (bool): if True, submit calculation
        from_scratch (bool): if True, remove previous calculation
        from_relax (bool): if True, copy the CONTCAR from the relaxation folder
            If False, copy the POSCAR from compound_path
    """
    computing_config_dict = get_computing_config_dict()
    oqmd_id = compound_path.split("/")[1]
    if from_relax:
        relax_done = check_relax(compound_path, mode="rlx")
        volume_converged = check_volume_difference(compound_path)
        if not relax_done or not volume_converged:
            logger.info("Not ready to set up elastic calculation")
            return
        logger.info("Relaxations successful")
        relax_path = os.path.join(compound_path, "rlx")
        mode = "bulkmod_rlx"
    else:
        mode = "bulkmod"
    bulkmod_path = os.path.join(compound_path, mode)

    if from_scratch:
        if os.path.exists(bulkmod_path):
            shutil.rmtree(bulkmod_path)
    if not os.path.exists(bulkmod_path):
        os.mkdir(bulkmod_path)

    # POSCAR
    poscar_path = os.path.join(bulkmod_path, "POSCAR")
    if from_relax:
        contcar_path = os.path.join(relax_path, "CONTCAR")
        shutil.copy(contcar_path, poscar_path)
    else:
        orig_poscar_path = os.path.join(compound_path, "POSCAR")
        shutil.copy(orig_poscar_path, poscar_path)
    structure = pmg.core.Structure.from_file(poscar_path)

    # POTCAR
    potcar_path = os.path.join(bulkmod_path, "POTCAR")
    make_potcar(structure, potcar_path)

    # INCAR
    incar_path = os.path.join(bulkmod_path, "INCAR")
    make_incar(incar_path, mode=mode)

    # vasp.q
    vaspq_path = os.path.join(bulkmod_path, "vasp.q")
    make_vaspq(vaspq_path, mode=mode, jobname=oqmd_id)

    logger.info("Making strain directories")
    strains = np.power(np.linspace(0.8, 1.2, 11), 1 / 3)
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
        sed_call_poscar = f"sed -i '2 c\{strain}' {strain_poscar_path}"
        subprocess.call(sed_call_poscar, shell=True)

        with change_directory(strain_path):
            for f in ["POTCAR", "INCAR"]:
                orig_path = "../" + f
                f_path = f
                os.symlink(orig_path, f_path, target_is_directory=False)

    if submit:
        mode = "bulkmod"
        if from_relax:
            mode += "_rlx"
        submit_job(compound_path, mode=mode)


def analyze_bulkmod(compound_path, from_relax=False):
    """
    Fit an EOS to calculate the bulk modulus from a finished bulkmod calculation

    Args:
        compound_path (str): base path of calculation
        from_relax (bool): if True, copy the CONTCAR from the relaxation folder
            If False, copy the POSCAR from compound_path
    """
    mode = "bulkmod"
    if from_relax:
        mode += "_rlx"
    bulkmod_path = os.path.join(compound_path, mode)

    jobid_path = os.path.join(bulkmod_path, "jobid")
    if os.path.exists(jobid_path):
        with open(jobid_path, "r") as fr:
            job_id = fr.readline().splitlines()[0]
        job_complete = check_job_complete(job_id)
        if not job_complete:
            logger.info(f"{mode} not finished")
            return -1
    else:
        raise Exception(f"No job id exists in {bulkmod_path}")

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
    eos_analyzer = BirchMurnaghan(volumes, final_energies)
    eos_analyzer.fit()
    bulk_modulus = np.round(eos_analyzer.b0_GPa, 3)
    logger.info(f"BULK MODULUS: {bulk_modulus}")
    return bulk_modulus


def manage_calculations(compound_path, calculation_types):
    if "rlx-coarse" in calculation_types:
        coarse_rlx_path = os.path.join(compound_path, "rlx-coarse")
        if not os.path.exists(coarse_rlx_path):
            setup_coarse_relax(compound_path, submit=True)
            return

        coarse_rlx_successful = check_relax(
            compound_path, rerun_relax=True, submit=True, mode="rlx-coarse"
        )
        if not coarse_rlx_successful:
            return
        from_coarse = True
    else:
        from_coarse = False

    if "rlx-fine" in calculation_types:
        rlx_path = os.path.join(compound_path, "rlx")
        if not os.path.exists(rlx_path):
            setup_relax(compound_path, submit=True, from_coarse=from_coarse)
            return

        rlx_successful = check_relax(
            compound_path, rerun_relax=True, submit=True, mode="rlx"
        )
        if not rlx_successful:
            return
        from_relax = True
    else:
        from_relax = False

    if "bulkmod" in calculation_types:
        if from_relax:
            bulkmod_path = os.path.join(compound_path, "bulkmod_rlx")
        else:
            bulkmod_path = os.path.join(compound_path, "bulkmod")

        if not os.path.exists(bulkmod_path):
            setup_simple_bulkmod(compound_path, submit=True, from_relax=from_relax)
            return

        analyze_bulkmod(compound_path, from_relax=from_relax)


if __name__ == "__main__":
    calculation_types = ["rlx-coarse", "rlx-fine", "bulkmod"]

    if not os.path.exists("calculations"):
        make_calculations_folder()

    compound_paths = [d for d in glob.glob("calculations/*") if os.path.isdir(d)]
    # Sort the paths by name
    compound_paths = sorted(compound_paths, key=lambda d: int(d.split("/")[1]))

    for compound_path in compound_paths:
        compound_name = compound_path.split("/")[1]
        print(compound_name)
        manage_calculations(compound_path, calculation_types)
        print("")
