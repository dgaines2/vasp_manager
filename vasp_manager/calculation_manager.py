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
from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from .utils import change_directory, get_primitive_structure_from_poscar
from .vasp_utils import (
    analyze_bulkmod,
    analyze_elastic,
    make_archive,
    make_bulkmod_strains,
    make_incar,
    make_potcar,
    make_vaspq,
)

logger = logging.getLogger(__name__)

computing_config_dict = json.loads(
    pkgutil.get_data("vasp_manager", "config/computing_config.json").decode("utf-8")
)
computer = computing_config_dict["computer"]


def submit_job(compound_path, mode, ignore_errors=False):
    """Call SLURM sbatch for the calculation and log the jobid"""
    if "personal" in computer:
        ignore_errors = True
        error_msg = f"Cannot submit {mode.upper()} job for on personal computer"
        error_msg += "\n\tIgnoring job submission..."
        logger.debug(error_msg)
        return True

    calc_path = os.path.join(compound_path, mode)
    jobid_path = os.path.join(calc_path, "jobid")
    if os.path.exists(jobid_path):
        logger.info(f"{mode.upper()} Job already exists")
        return True

    vaspq_location = os.path.join(calc_path, "vasp.q")
    if not os.path.exists(vaspq_location):
        logger.info(f"No vasp.q file in {calc_path}")
        # return False here instead of catching an exception
        # This enables job resubmission
        return False
    submission_call = "sbatch vasp.q | awk '{ print $4 }' | tee jobid"
    with change_directory(calc_path):
        subprocess.call(submission_call, shell=True)
    return True


def check_job_complete(jobid, ignore_errors=False):
    """Return True if job done"""
    if computer == "personal":
        ignore_errors = True
        error_msg = "Cannot check job on personal computer"
        error_msg += "\n\tIgnoring job status check..."
        logger.debug(error_msg)
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


def check_relax(
    compound_path,
    mode="rlx-coarse",
    rerun_relax=False,
    submit=False,
    tail=5,
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
    if mode not in ["rlx", "rlx-coarse"]:
        raise Exception(f"Mode '{mode}' is incorrect or not supported")
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
            logger.info(f"{mode.upper()} not finished")
            return False

        grep_call = f"tail -n{tail} {stdout_path}"
        grep_output = (
            subprocess.check_output(grep_call, shell=True).decode("utf-8").strip()
        )
        if "reached required accuracy" in grep_output:
            logger.info(f"{mode.upper()} Calculation: reached required accuracy")
            logger.debug(grep_output)
            return True
        else:
            archive_dirs = glob.glob(f"{relax_path}/archive*")
            if "coarse" in mode and len(archive_dirs) >= 3:
                logger.warning("Many archives exist, suggest force based relaxation")
                if rerun_relax:
                    setup_relax(compound_path, submit=submit, rerun_relax=True)
                return True

            logger.warning(f"{mode.upper()} FAILED")
            logger.debug(grep_output)
            if rerun_relax:
                logger.info(f"Rerunning {compound_path}")
                if "coarse" in mode:
                    setup_coarse_relax(compound_path, submit=submit, rerun_relax=True)
                else:
                    setup_relax(compound_path, submit=submit, rerun_relax=True)
            return False
    else:
        logger.info(f"{mode.upper()} not started")
        return False


def check_volume_difference(compound_path, rerun_relax=False, submit=False):
    """
    check relaxation runs for volume difference

    if abs(volume difference) is >= 5%, rerun relaxation
    only check for mode='rlx' as that's the structure for elastic analysis

    Args:
        compound_path (str): base path of calculation
        rerun_relax (bool): if True, make an archive and re-setup the rlx folder
        submit (bool): if True and rerun_relax, submit calculation

    Returns:
        volume_converged (bool): if True, relaxation completed successfully
    """
    relaxation_done = check_relax(
        compound_path, rerun_relax=rerun_relax, submit=submit, mode="rlx"
    )
    if not relaxation_done:
        return False

    relax_path = os.path.join(compound_path, "rlx")
    poscar_path = os.path.join(relax_path, "POSCAR")
    contcar_path = os.path.join(relax_path, "CONTCAR")
    try:
        p_structure = pmg.core.Structure.from_file(poscar_path)
        c_structure = pmg.core.Structure.from_file(contcar_path)
    except Exception as e:
        logger.info(f"  RLX CONTCAR doesn't exist or is empty: {e}")
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
        logger.info(f"  RLX volume converged")
        logger.info(f"  dV = {volume_diff:.4f}")
        volume_converged = True
    return volume_converged


def setup_coarse_relax(
    compound_path,
    submit=False,
    rerun_relax=False,
):
    """
    Set up a coarse relaxation

    Args:
        compound_path (str): base path of calculation
        submit (bool): if True, submit calculation
        rerun_relax (bool): if True, make an archive and resubmit
    """
    # POSCAR, POTCAR, INCAR, vasp.q
    oqmd_id = compound_path.split("/")[1]
    logger.info(f"{oqmd_id} Setting up rlx-coarse")
    mode = "rlx-coarse"

    crelax_path = os.path.join(compound_path, mode)
    if not os.path.exists(crelax_path):
        os.mkdir(crelax_path)

    # POSCAR
    if rerun_relax:
        archive_made = make_archive(compound_path, mode=mode)
        if not archive_made:
            # set rerun_relax to not make an achive and instead
            # continue to make the input files
            setup_coarse_relax(compound_path, submit=submit, rerun_relax=False)
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
        job_status = submit_job(compound_path, mode=mode)
        if not job_submitted:
            setup_coarse_relax(compound_path, submit=True, rerun_relax=True)


def setup_relax(
    compound_path,
    submit=False,
    rerun_relax=False,
    from_coarse=False,
):
    """
    Set up a fine relaxation

    Args:
        compound_path (str): base path of calculation
        submit (bool): if True, submit calculation
        rerun_relax (bool): if True, make an archive and if submit, resubmit
        from_coarse (bool): if True, copy the CONTCAR from the coarse relaxation
            folder
            If False, copy the POSCAR from compound_path
    """
    oqmd_id = compound_path.split("/")[1]
    logger.info(f"{oqmd_id} Setting up rlx")

    mode = "rlx"
    if from_coarse:
        crelax_done = check_relax(compound_path, mode="rlx-coarse")
        if not crelax_done:
            logger.info("RLX-COARSE not finished")
            return

    relax_path = os.path.join(compound_path, mode)
    if not os.path.exists(relax_path):
        os.mkdir(relax_path)

    if rerun_relax:
        archive_made = make_archive(compound_path, mode=mode)
        if not archive_made:
            # set rerun_relax to not make an achive and instead
            # continue to make the input files
            setup_relax(
                compound_path, submit=submit, rerun_relax=False, from_coarse=from_coarse
            )
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
        job_status = submit_job(compound_path, mode=mode)
        if not job_submitted:
            setup_relax(
                compound_path, submit=True, rerun_relax=True, from_coarse=from_coarse
            )


def setup_bulkmod(
    compound_path,
    submit=True,
    from_relax=False,
):
    """
    Set up EOS bulkmod calculation

    Args:
        compound_path (str): base path of calculation
        submit (bool): if True, submit calculation
        from_relax (bool): if True, copy the CONTCAR from the relaxation folder
            If False, copy the POSCAR from compound_path
    """
    oqmd_id = compound_path.split("/")[1]
    logger.info(f"{oqmd_id} Setting up BULKMOD Calculation")
    if from_relax:
        relax_done = check_relax(compound_path, mode="rlx")
        volume_converged = check_volume_difference(compound_path)
        mode = "bulkmod_rlx"
        if not relax_done or not volume_converged:
            logger.info("Not ready to set up {mode.upper()} calculation")
            logger.info("relax_done = {relax_done}")
            logger.info("volume_converged = {volume_converged}")
            return
        logger.info("Relaxations Successful")
        relax_path = os.path.join(compound_path, "rlx")
    else:
        mode = "bulkmod"
    bulkmod_path = os.path.join(compound_path, mode)

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

    strains = np.power(np.linspace(0.8, 1.2, 11), 1 / 3)
    make_bulkmod_strains(bulkmod_path, strains)

    if submit:
        mode = "bulkmod"
        if from_relax:
            mode += "_rlx"
        job_status = submit_job(compound_path, mode=mode)
        if not job_submitted:
            setup_bulkmod(compound_path, submit=True, rerun_relax=from_relax)


def check_bulkmod(bulkmod_path, from_relax=False):
    jobid_path = os.path.join(bulkmod_path, "jobid")
    if os.path.exists(jobid_path):
        with open(jobid_path, "r") as fr:
            job_id = fr.readline().splitlines()[0]
        job_complete = check_job_complete(job_id)
        if not job_complete:
            logger.info(f"{mode.upper()} job not finished")
        return False
    else:
        raise Exception(f"No job id exists in {bulkmod_path}")
    return True


def setup_elastic(compound_path, submit=True, increase_nodes=False):
    """
    Run elastic constants routine through VASP
    By default, requires relaxation (as the elastic constants routine needs
        the cell to be nearly at equilibrium)

    Args:
        compound_path (str): base path of calculation
        submit (bool): if True, submit calculation
        increase_nodes (bool): if True, copy the CONTCAR from the relaxation folder
    """
    oqmd_id = compound_path.split("/")[1]
    logger.info(f"{oqmd_id} Setting up ELASTIC Calculation")
    # POSCAR, POTCAR, KPOINTS, INCAR, vasp.q
    relax_done = check_relax(compound_path, mode="rlx")
    volume_converged = check_volume_difference(compound_path)
    if not relax_done or not volume_converged:
        logger.info("Not ready to set up ELASTIC calculation")
        logger.info("relax_done = {relax_done}")
        logger.info("volume_converged = {volume_converged}")
        return
    logger.info("Relaxations Successful")

    relax_path = os.path.join(compound_path, "rlx")
    elastic_path = os.path.join(compound_path, "elastic")

    if increase_nodes:
        shutil.rmtree(elastic_path)
    if not os.path.exists(elastic_path):
        os.mkdir(elastic_path)

    # POSCAR
    contcar_path = os.path.join(relax_path, "CONTCAR")
    poscar_path = os.path.join(elastic_path, "POSCAR")
    structure = get_primitive_structure_from_poscar(contcar_path)
    poscar = Poscar(structure)
    poscar.write_file(poscar_path)

    # POTCAR
    potcar_path = os.path.join(elastic_path, "POTCAR")
    make_potcar(structure, potcar_path)

    # KPOINTS
    # shutil.copy("../static_files/KPOINTS", elastic_path)

    # INCAR
    incar_path = os.path.join(elastic_path, "INCAR")
    make_incar(incar_path, mode="elastic")
    # shutil.copy("../static_files/INCAR-elastic", incar_path)

    # vasp.q
    vaspq_path = os.path.join(elastic_path, "vasp.q")
    make_vaspq(
        vaspq_path, mode="elastic", jobname=oqmd_id, increase_nodes=increase_nodes
    )

    if submit:
        job_status = submit_job(compound_path, mode=mode)
        if not job_submitted:
            setup_bulkmod(compound_path, submit=True, increase_nodes=increase_nodes)


def check_elastic(elastic_path, rerun=False, submit=False, tail=5):
    """
    Check result of elastic calculation

    Args:
        compound_path (str): base path of calculation
        rerun (bool): if True, rerun failed jobs
        submit (bool): if True, submit failed jobs
        tail (int): If job fails, print {tail} lines from stdout.txt

    Returns
        elastic_successful (bool): if True, elastic calculation completed successfully
    """
    mode = "elastic"
    jobid_path = os.path.join(elastic_path, "jobid")
    if os.path.exists(jobid_path):
        with open(jobid_path, "r") as fr:
            job_id = fr.readline().splitlines()[0]
        job_complete = check_job_complete(job_id)
        if not job_complete:
            logger.info(f"{mode.upper()} not finished")
            return False
    else:
        raise Exception(f"No job id exists in {elastic_path}")

    stdout_path = os.path.join(elastic_path, "stdout.txt")
    if os.path.exists(stdout_path):
        grep_call = f"grep 'Total' {stdout_path}"
        grep_output = (
            subprocess.check_output(grep_call, shell=True).decode("utf-8").splitlines()
        )
        last_grep_line = grep_output[-1].strip().split()
        # last grep line looks something like 'Total: 36/ 36'
        finished_deformations = int(last_grep_line[-2].replace("/", ""))
        total_deformations = int(last_grep_line[-1])
        logger.debug(last_grep_line)
        if finished_deformations == total_deformations:
            logger.info(f"{mode.upper()} Calculation: Success")
            return True
        else:
            grep_call = f"tail -n{tail} {stdout_path}"
            grep_output = (
                subprocess.check_output(grep_call, shell=True).decode("utf-8").strip()
            )
            logger.info(grep_output)
            logger.info(f"{mode.upper()} Calculation: FAILED")
            if rerun:
                setup_elastic(compound_path, submit=submit, increase_nodes=True)
            return False
    else:
        # shouldn't get here unless function was called with submit=False
        logger.info("{mode.upper()} Calculation: No stdout.txt available")
        if rerun:
            # increase nodes as its likely the calculation failed
            setup_elastic(compound_path, submit=submit, increase_nodes=True)
        return False


def manage_calculation(compound_path, calculation_types, submit=True, rerun_relax=True):
    if "rlx-coarse" in calculation_types:
        coarse_rlx_path = os.path.join(compound_path, "rlx-coarse")
        # if rlx-coarse doesn't exist, set it up
        if not os.path.exists(coarse_rlx_path):
            setup_coarse_relax(compound_path, submit=submit)
            return
        # else, check if it finished
        check_relax(compound_path, mode="rlx-coarse", rerun_relax=rerun_relax)
        # if rlx-coarse fails, it'll be caught in check_relax
        # and setup_relax won't run
        from_coarse = True
    else:
        from_coarse = False

    if "rlx-fine" in calculation_types:
        rlx_path = os.path.join(compound_path, "rlx")
        # if rlx doesn't exist, set it up
        if not os.path.exists(rlx_path):
            setup_relax(compound_path, submit=submit, from_coarse=from_coarse)
            return
        # else, check if it finished
        check_relax(compound_path, mode="rlx", rerun_relax=rerun_relax, submit=submit)
        # if rlx fails, it'll be caught in check_relax
        # and setup_bulkmod (if from_relax=True) won't run
        from_relax = True
    else:
        from_relax = False

    if "bulkmod" in calculation_types:
        if from_relax:
            bulkmod_path = os.path.join(compound_path, "bulkmod_rlx")
        else:
            bulkmod_path = os.path.join(compound_path, "bulkmod")

        # if bulkmod doesn't exist, set it up
        if not os.path.exists(bulkmod_path):
            setup_bulkmod(compound_path, submit=submit, from_relax=from_relax)
            return
        # else, check it, then analyze it
        bulkmod_successful = check_bulkmod(bulkmod_path, from_relax=from_relax)
        if not bulkmod_successful:
            return
        analyze_bulkmod(compound_path, from_relax=from_relax)
        # if it fails, don't auto resubmit so you can find the problem

    if "elastic" in calculation_types:
        # if elastic doesn't exist, set it up
        elastic_path = os.path.join(compound_path, "elastic")
        if not os.path.exists(elastic_path):
            setup_elastic(compound_path, submit=submit)
            return
        # else, check if it finished
        elastic_successful = check_elastic(elastic_path, rerun=True)
        # else, check it, then analyze it
        if not elastic_successful:
            return
        analyze_elastic(elastic_path)
        # if it fails, don't auto resubmit so you can find the problem


def manage_calculations(calculation_types):
    compound_paths = [d for d in glob.glob("calculations/*") if os.path.isdir(d)]
    # Sort the paths by name
    compound_paths = sorted(compound_paths, key=lambda d: int(d.split("/")[1]))

    for compound_path in compound_paths:
        compound_name = compound_path.split("/")[1]
        print(compound_name)
        manage_calculation(compound_path, calculation_types)
        print("\n")
