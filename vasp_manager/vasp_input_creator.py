# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import glob
import json
import logging
import os
import pkgutil
import shutil
import subprocess

import pymatgen as pmg

from vasp_manager.utils import change_directory, get_pmg_structure_from_poscar

logger = logging.getLogger(__name__)


class VaspInputCreator:
    def __init__(
        self, calc_path, mode, poscar_source_path, name=None, increase_nodes=False
    ):
        self.calc_path = calc_path
        self.poscar_source_path = poscar_source_path
        self.increase_nodes = increase_nodes
        self.name = name
        self.mode = self._get_mode(mode)

    def _get_mode(self, mode):
        # rlx-coarse, rlx, bulkmod, or elastic
        # needed to add this to ensure bulkmod or bulkmod_standalone
        # share same config
        if "bulkmod" in mode:
            mode = "bulkmod"
        return mode

    @property
    def computer(self):
        return self.computing_config_dict["computer"]

    @property
    def structure(self):
        try:
            structure = get_pmg_structure_from_poscar(self.poscar_source_path)
        except Exception as e:
            raise Exception(f"Cannot load POSCAR in {self.poscar_source_path}: {e}")
        return structure

    @property
    def calc_config_dict(self):
        calc_config_dict = json.loads(
            pkgutil.get_data("vasp_manager", "config/calc_config.json").decode("utf-8")
        )
        return calc_config_dict

    @property
    def computing_config_dict(self):
        computing_config_dict = json.loads(
            pkgutil.get_data("vasp_manager", "config/computing_config.json").decode(
                "utf-8"
            )
        )
        return computing_config_dict

    @property
    def incar_template(self):
        incar_template = pkgutil.get_data(
            "vasp_manager", "static_files/INCAR_template"
        ).decode("utf-8")
        return incar_template

    @property
    def potcar_dict(self):
        potcar_dict = json.loads(
            pkgutil.get_data("vasp_manager", "static_files/pot_dict.json").decode(
                "utf-8"
            )
        )
        return potcar_dict

    @property
    def q_mapper(self):
        q_mapper = json.loads(
            pkgutil.get_data("vasp_manager", "static_files/q_handles.json").decode(
                "utf-8"
            )
        )
        return q_mapper

    def make_poscar(self):
        """
        Create and write a POSCAR
        """
        poscar = pmg.io.vasp.Poscar(self.structure)
        poscar_path = os.path.join(self.calc_path, "POSCAR")
        poscar.write_file(poscar_path)

    @property
    def n_nodes(self):
        # start with 1 node per 32 atoms
        num_nodes = (len(self.structure) // 32) + 1
        if self.computer == "quest":
            # quest has small nodes
            num_nodes *= 2
        if self.increase_nodes:
            n_nodes *= 2
        return num_nodes

    def make_potcar(self):
        """
        Create and write a POTCAR
        """
        potcar_path = os.path.join(self.calc_path, "POTCAR")
        potcar_dir = self.computing_config_dict[self.computer]["potcar_dir"]

        el_names = [el.name for el in self.structure.composition]
        logger.debug(f"{self.structure.composition.reduced_formula}, {el_names}")
        pot_singles = [
            os.path.join(potcar_dir, self.potcar_dict[el_name], "POTCAR")
            for el_name in el_names
        ]
        for pot_single in pot_singles:
            if not os.path.exists(pot_single):
                msg = "Unable to create POTCAR"
                msg += f"\n\t POTCAR not found at path {pot_single}"
                raise Exception(msg)
        cmd = "cat " + " ".join(pot_singles) + " > " + potcar_path
        subprocess.call(cmd, shell=True)

    def make_incar(self):
        """
        Create and write an INCAR

        Need to modify this to account for spin/magmom
        Current kpoints coming from the kspacing tag in the INCAR,
            but future versions should include ability to make kpoints from kppra
        """
        incar_path = os.path.join(self.calc_path, "INCAR")
        ncore = self.computing_config_dict[self.computer]["ncore"]
        calc_config = self.calc_config_dict[self.mode]

        # Add lines to the vaspq file for only elastic calculations
        incar_tmp = self.incar_template
        if self.mode == "elastic":
            # add extra flags for elastic mode
            incar_tmp = incar_tmp.split("\n")
            for i, line in enumerate(incar_tmp):
                if "KSPACING" in line:
                    nfree_line = "NFREE = {nfree}"
                    symprec_line = "SYMPREC = {symprec}"
                    incar_tmp.insert(i + 1, symprec_line)
                    incar_tmp.insert(i + 1, nfree_line)
                if "NCORE" in line:
                    # elastic calculation won't run unless NCORE=1
                    incar_tmp[i] = f"NCORE = 1"
            incar_tmp = "\n".join([line for line in incar_tmp])
        incar = incar_tmp.format(**calc_config, ncore=ncore)
        logger.debug(incar)
        with open(incar_path, "w+") as fw:
            fw.write(incar)

    def make_vaspq(self):
        """
        Create and write vasp.q file
        """
        vaspq_path = os.path.join(self.calc_path, "vasp.q")
        calc_config = self.calc_config_dict[self.mode]
        walltime = calc_config["walltime"]

        # create pad string for job naming to differentiate in the queue
        if self.mode == "rlx" or self.mode == "rlx-coarse":
            if self.mode == "rlx":
                pad_string = "r"
            elif self.mode == "rlx-coarse":
                pad_string = "rc"
            mode = "rlx"
        elif self.mode == "static":
            pad_string = "s"
            mode = "static"
        elif self.mode == "bulkmod":
            pad_string = "b"
            mode = "bulkmod"
        elif self.mode == "elastic":
            pad_string = "e"
            mode = "elastic"

        n_procs = (
            self.n_nodes * self.computing_config_dict[self.computer]["ncore_per_node"]
        )

        if self.name is None:
            jobname = pad_string + self.structure.composition.reduced_formula
        else:
            jobname = pad_string + self.name

        if self.increase_nodes:
            n_procs *= 2
            hours, minutes, seconds = walltime.split(":")
            hours = str(int(hours) * 2)
            walltime = ":".join([hours, minutes, seconds])

        computer_config = self.computing_config_dict[self.computer].copy()
        computer_config.update(
            {"n_nodes": self.n_nodes, "n_procs": n_procs, "jobname": jobname}
        )

        q_name = self.q_mapper[self.computer][mode]
        vaspq_tmp = pkgutil.get_data("vasp_manager", f"static_files/{q_name}").decode(
            "utf-8"
        )
        vaspq = vaspq_tmp.format(**computer_config, walltime=walltime)
        logger.debug(vaspq)
        with open(vaspq_path, "w+") as fw:
            fw.write(vaspq)

    def make_archive(self):
        """
        Make an archive of a VASP calculation and copy back over relevant files

        Returns:
            archive_made (bool): if job was never submitted, return False
        """
        # if job was never submitted, don't make an archive
        jobid_path = os.path.join(self.calc_path, "jobid")
        if not os.path.exists(jobid_path):
            return False
        # check if CONTCAR is empty -- calculation failed almost immediately
        # if it is empty, don't make an archive, just recreate the files
        contcar_path = os.path.join(self.calc_path, "CONTCAR")
        if os.stat(contcar_path).st_size == 0:
            os.path.remove(jobid_path)
            return True

        with change_directory(self.calc_path):
            num_archives = len(glob.glob("archive*"))
            all_files = [d for d in glob.glob("*") if os.path.isfile(d)]
            archive_name = f"archive_{num_archives}"
            logger.info(f"Making {archive_name}...")
            os.mkdir(archive_name)
            for f in all_files:
                shutil.move(f, archive_name)

        contcar_path = os.path.join(self.calc_path, archive_name, "CONTCAR")
        self.poscar_source_path = contcar_path
        self.create()
        return True

    def create(self):
        """
        Make VASP input files
        """
        if not os.path.exists(self.calc_path):
            os.mkdir(self.calc_path)
        self.make_poscar()
        self.make_potcar()
        self.make_incar()
        self.make_vaspq()
