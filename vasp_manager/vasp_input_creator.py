# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import json
import logging
import os
import shutil
import warnings
from datetime import time, timedelta
from functools import cached_property
from pathlib import Path

import importlib_resources
import numpy as np
import yaml
from pymatgen.io.vasp import Poscar, Potcar

from vasp_manager.utils import change_directory, get_pmg_structure_from_poscar, pcat

logger = logging.getLogger(__name__)


class VaspInputCreator:
    """
    Handles VASP file creation
    """

    def __init__(
        self,
        calc_path,
        mode,
        poscar_source_path,
        config_dir=None,
        primitive=True,
        name=None,
        increase_nodes_by_factor=1,
        increase_walltime_by_factor=1,
        poscar_significant_figures=8,
        ncore_per_node_for_memory=0,
        use_spin=True,
    ):
        """
        Args:
            calc_path
            mode
            poscar_source_path
            config_dir
            primitive
            name
            increase_nodes_by_factor
            increase_walltime_by_factor
            poscar_significant_figures
            ncore_per_node_for_memory
            use_spin -- pass False to supress spin polarization
        """
        self.calc_path = Path(calc_path)
        self.poscar_source_path = Path(poscar_source_path)
        self.config_dir = Path(config_dir) if config_dir else self.calc_path.parent.parent
        self.primitive = primitive
        self.increase_nodes_by_factor = int(increase_nodes_by_factor)
        self.increase_walltime_by_factor = int(increase_walltime_by_factor)
        self.name = name
        self.mode = mode
        self.poscar_significant_figures = poscar_significant_figures
        self.ncore_per_node_for_memory = ncore_per_node_for_memory
        self.use_spin = use_spin

    @cached_property
    def calc_config(self):
        fname = "calc_config.json"
        fpath = self.config_dir / fname
        if fpath.exists():
            with open(fpath) as fr:
                calc_config_dict = json.load(fr)
                calc_config = calc_config_dict[self.mode]
        else:
            raise Exception(f"No {fname} found in path {self.config_dir.absolute()}")
        # if material_name/mode calc_config.json exists, update with that
        # This allows for custom settings for a specific material/mode to be used
        custom_calc_config_path = self.calc_path / fname
        if custom_calc_config_path.exists():
            with open(custom_calc_config_path) as fr:
                custom_calc_config = json.load(fr)
            calc_config.update(custom_calc_config)
        return calc_config

    @cached_property
    def computing_config_dict(self):
        fname = "computing_config.json"
        fpath = self.config_dir / fname
        if fpath.exists():
            with open(fpath) as fr:
                computing_config = json.load(fr)
        else:
            raise Exception(f"No {fname} found in path {self.config_dir.absolute()}")
        return computing_config

    @cached_property
    def computing_config(self):
        return self.computing_config_dict[self.computer]

    @cached_property
    def computer(self):
        return self.computing_config_dict["computer"]

    @cached_property
    def source_structure(self):
        num_archives = len(list(self.calc_path.glob("archive*")))
        if num_archives > 0:
            archive_name = f"archive_{num_archives-1}"
            self.poscar_source_path = self.calc_path / archive_name / "CONTCAR"
        try:
            structure = get_pmg_structure_from_poscar(
                self.poscar_source_path, primitive=self.primitive
            )
        except Exception as e:
            raise Exception(f"Cannot load POSCAR in {self.poscar_source_path}: {e}")
        return structure

    @cached_property
    def incar_template(self):
        incar_template = (
            importlib_resources.files("vasp_manager")
            .joinpath(str(Path("static_files") / "INCAR_template"))
            .read_text()
        )
        return incar_template

    @cached_property
    def potcar_dict(self):
        potcar_dict = json.loads(
            importlib_resources.files("vasp_manager")
            .joinpath(str(Path("static_files") / "pot_dict.json"))
            .read_text()
        )
        return potcar_dict

    @cached_property
    def q_mapper(self):
        q_mapper = json.loads(
            importlib_resources.files("vasp_manager")
            .joinpath(str(Path("static_files") / "q_handles.json"))
            .read_text()
        )
        return q_mapper

    def make_poscar(self):
        """
        Create and write a POSCAR
        """
        poscar = Poscar(self.source_structure)
        poscar_path = self.calc_path / "POSCAR"
        poscar.write_file(
            poscar_path, significant_figures=self.poscar_significant_figures
        )

    def make_potcar(self):
        """
        Create and write a POTCAR
        """
        potcar_path = self.calc_path / "POTCAR"
        potcar_dir = Path(self.computing_config["potcar_dir"])

        el_names = [el.name for el in self.source_structure.composition]
        logger.debug(f"{self.source_structure.composition.reduced_formula}, {el_names}")
        pot_singles = [
            potcar_dir / self.potcar_dict[el_name] / "POTCAR" for el_name in el_names
        ]
        for pot_single in pot_singles:
            if not pot_single.exists():
                msg = "Unable to create POTCAR"
                msg += f"\n\t POTCAR not found at path {pot_single}"
                raise Exception(msg)

        potcar = pcat(pot_singles)
        with open(potcar_path, "w+") as fw:
            fw.write(potcar)

    @cached_property
    def n_nodes(self):
        # start with 1 node per 32 atoms
        num_nodes = (len(self.source_structure) // 32) + 1
        if self.computer == "quest":
            # quest has ~4x smaller nodes than perlmutter
            num_nodes *= 4
        num_nodes *= self.increase_nodes_by_factor
        return num_nodes

    @cached_property
    def n_procs(self):
        n_procs = self.n_nodes * self.computing_config["ncore_per_node"]
        return n_procs

    @cached_property
    def n_procs_used(self):
        # typically request all processors on each node, and then
        # leave some ~4/node empty for memory
        ncore_per_node = self.computing_config["ncore_per_node"]
        if self.computer == "quest":
            self.ncore_per_node_for_memory += 4
        return self.n_nodes * (ncore_per_node - self.ncore_per_node_for_memory)

    @cached_property
    def d_f_block(self):
        d_f_block = json.loads(
            importlib_resources.files("vasp_manager")
            .joinpath(str(Path("static_files") / "d_f_block.json"))
            .read_text()
        )
        return d_f_block

    def _check_needs_spin_polarization(self, composition_dict):
        needs_spin_polarization = False
        for el_name in composition_dict:
            if (
                el_name in self.d_f_block["d_block"]
                or el_name in self.d_f_block["f_block"]
            ):
                needs_spin_polarization = True
        return needs_spin_polarization

    def _get_auto_magmom(self, composition_dict):
        magmom_string = "MAGMOM ="
        for el_name in composition_dict:
            if el_name in self.d_f_block["d_block"]:
                init_mag = "5.0"
            elif el_name in self.d_f_block["f_block"]:
                init_mag = "7.0"
            else:
                init_mag = "0.0"
            magmom_atom = str(int(composition_dict[el_name])) + "*" + init_mag
            magmom_string += f" {magmom_atom}"
        return magmom_string

    @cached_property
    def hubbards(self):
        hubbards = json.loads(
            importlib_resources.files("vasp_manager")
            .joinpath(str(Path("static_files") / "hubbards.json"))
            .read_text()
        )
        return hubbards

    def _check_needs_dftu(self, hubbards_type, composition_dict):
        if not hubbards_type:
            return False

        needs_dftu = False
        for el_name in composition_dict:
            if el_name in self.hubbards[hubbards_type]:
                needs_dftu = True
        return needs_dftu

    def _get_ldau_string(self, hubbards_type, composition_dict):
        u_string = "LDAUU ="
        j_string = "LDAUJ ="
        l_string = "LDAUL ="
        hubbards = self.hubbards[hubbards_type]

        for el_name in composition_dict:
            if el_name in hubbards:
                u_val = hubbards[el_name]["U"]
                l_val = hubbards[el_name]["L"]
            else:
                u_val = 0.0
                l_val = -1
            u_string += f" {float(u_val)}"
            j_string += " 0.0"
            l_string += f" {l_val}"

        ldau_string = "LDAU = .TRUE.\n"
        ldau_string += "LDAUPRINT = 1\n"
        ldau_string += f"{u_string}\n"
        ldau_string += f"{j_string}\n"
        ldau_string += f"{l_string}"
        return ldau_string

    def _get_lmaxmix(self, composition_dict):
        d_block_elements = self.d_f_block["d_block"]
        f_block_elements = self.d_f_block["f_block"]
        lmaxmix = 2
        for element in composition_dict:
            if element in d_block_elements:
                lmaxmix = 4
        for element in composition_dict:
            if element in f_block_elements:
                lmaxmix = 6
        return lmaxmix

    def make_incar(self):
        """
        Create and write an INCAR

        Current kpoints coming from the kspacing tag in the INCAR,
            but future versions should include ability to make kpoints from kppra
        """
        incar_path = self.calc_path / "INCAR"
        calc_config = self.calc_config.copy()

        if calc_config["ispin"] not in [1, "auto"]:
            raise RuntimeError("ISPIN must be set to 1 or auto")

        if calc_config["hubbards"] not in [None, False, "wang"]:
            raise RuntimeError("hubbards must be set to None, False, or 'wang'")
        if calc_config["hubbards"] == "wang" and calc_config["gga"] != "PE":
            raise RuntimeError("DFT+U is only support for PBE")

        if calc_config["iopt"] != 0 and calc_config["potim"] != 0:
            raise RuntimeError("To use IOPT != 0, POTIM must be set to 0")

        composition_dict = self.source_structure.composition.as_dict()
        # read POTCAR
        potcar_path = self.calc_path / "POTCAR"
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            potcar = Potcar.from_file(potcar_path)
        n_electrons = 0
        for potcar_single in potcar:
            n_electrons += (
                potcar_single.nelectrons * composition_dict[potcar_single.element]
            )
        # make n_bands divisible by NCORE (VASP INCAR tag)
        # elastic calculation won't run unless NCORE=1
        if self.mode == "elastic":
            self.computing_config["ncore"] = 1
        ncore = self.computing_config["ncore"]
        calc_config["ncore"] = ncore
        cores_per_group = int(self.n_procs_used / self.calc_config["kpar"])
        calc_config["nbands"] = np.max(
            [
                int(np.ceil(0.75 * n_electrons / ncore) * ncore),
                cores_per_group,
            ]
        )

        needs_spin_polarization = self._check_needs_spin_polarization(composition_dict)
        use_spin_polarization = (
            needs_spin_polarization and calc_config["ispin"] == "auto" and self.use_spin
        )
        if use_spin_polarization:
            ispin = 2
            magmom_line = self._get_auto_magmom(composition_dict)
        else:
            ispin = 1

        needs_dftu = self._check_needs_dftu(calc_config["hubbards"], composition_dict)
        use_dftu = needs_dftu and calc_config["hubbards"]
        if use_dftu:
            lmaxmix = self._get_lmaxmix(composition_dict)
            lmaxmix_line = f"LMAXMIX = {lmaxmix}"
            ldau_string = self._get_ldau_string(calc_config["hubbards"], composition_dict)

        # Add lines to the incar file
        incar_tmp = self.incar_template.split("\n")
        for i, line in enumerate(incar_tmp):
            # add extra flags for spin polarization
            if "ISPIN" in line:
                incar_tmp[i] = f"ISPIN = {ispin}"
                if use_spin_polarization:
                    incar_tmp.insert(i + 1, magmom_line)
            if "LVTOT" in line and use_spin_polarization:
                lorbit_line = "LORBIT = 11"
                incar_tmp.insert(i + 1, lorbit_line)
            # add extra flags for DFT+U
            if "DFT+U" in line and use_dftu:
                incar_tmp.insert(i + 1, lmaxmix_line)
                incar_tmp.insert(i + 1, ldau_string)
            # add extra flags for elastic mode
            if self.mode == "elastic":
                if "KSPACING" in line:
                    nfree_line = "NFREE = {nfree}"
                    incar_tmp.insert(i + 1, nfree_line)
        incar_tmp = "\n".join([line for line in incar_tmp])
        incar = incar_tmp.format(**calc_config)
        logger.debug(incar)
        with open(incar_path, "w+") as fw:
            fw.write(incar)

    def make_kpoints(self):
        kpoints_path = self.calc_path / "KPOINTS"
        kspacing = self.calc_config["kspacing"]
        reciprocal_lattice = self.source_structure.lattice.reciprocal_lattice
        abc = np.asarray(reciprocal_lattice.abc)
        kpoints = [int(k) for k in np.ceil(abc / kspacing)]
        kpoints_str = " ".join(str(kpoint) for kpoint in kpoints)
        kpoints_text = f"kpoints generated from kspacing={kspacing}\n"
        kpoints_text += "0\n"
        kpoints_text += "Gamma\n"
        kpoints_text += f"{kpoints_str}\n"
        with open(kpoints_path, "w+") as fw:
            fw.write(kpoints_text)

    def make_vaspq(self):
        """
        Create and write vasp.q file
        """
        vaspq_path = self.calc_path / "vasp.q"

        # create pad string for job naming to differentiate in the queue
        match self.mode:
            case "rlx-coarse" | "rlx":
                if self.mode == "rlx":
                    pad_string = "r"
                elif self.mode == "rlx-coarse":
                    pad_string = "rc"
                mode = "rlx"
            case "static":
                pad_string = "s"
                mode = "static"
            case "bulkmod":
                pad_string = "b"
                mode = "bulkmod"
            case "elastic":
                pad_string = "e"
                mode = "elastic"
            case _:
                raise ValueError(
                    "Calculation type {self.mode} not in supported calculation types"
                    "of VaspInputCreator"
                )

        if self.name is None:
            jobname = pad_string + self.source_structure.composition.reduced_formula
        else:
            jobname = pad_string + self.name

        # convert walltime into seconds for increase_walltime_by_factor
        walltime_iso = time.fromisoformat(self.calc_config["walltime"])
        # walltime_duration is in seconds
        walltime_duration = timedelta(
            hours=walltime_iso.hour,
            minutes=walltime_iso.minute,
            seconds=walltime_iso.second,
        )
        walltime_duration *= self.increase_walltime_by_factor
        # convert to HH:MM:SS
        walltime = str(walltime_duration)
        # cut walltime short by 1 minute so job metrics log properly
        timeout = walltime_duration.seconds - 60
        # quest uses mpirun which needs timeout in seconds
        # otherwise, convert it back to HH:MM:SS
        if not self.computer == "quest":
            timeout = str(timedelta(seconds=timeout))

        computing_config = self.computing_config.copy()
        ncore_per_node = self.n_procs_used // self.n_nodes
        computing_config.update(
            {
                "n_nodes": self.n_nodes,
                "n_procs": self.n_procs,
                "ncore_per_node": ncore_per_node,
                "jobname": jobname,
                "walltime": walltime,
                "timeout": timeout,
            }
        )

        vaspq_settings_path = self.q_mapper[self.computer][mode]
        vaspq_settings = yaml.load(
            importlib_resources.files("vasp_manager")
            .joinpath(str(Path("static_files") / vaspq_settings_path))
            .read_text(),
            Loader=yaml.SafeLoader,
        )
        override_vaspq_settings_path = self.config_dir / f"{self.computer}.yml"
        if override_vaspq_settings_path.exists():
            with open(override_vaspq_settings_path) as fr:
                override_vaspq_settings = yaml.load(
                    fr,
                    Loader=yaml.SafeLoader,
                )
            vaspq_settings.update(override_vaspq_settings)
        vaspq_tmp = (
            importlib_resources.files("vasp_manager")
            .joinpath(str(Path("static_files") / "vasp.q"))
            .read_text()
        )
        vaspq_tmp = vaspq_tmp.format(**vaspq_settings)
        vaspq = vaspq_tmp.format(**computing_config)
        logger.debug(vaspq)
        with open(vaspq_path, "w+") as fw:
            fw.write(vaspq)

    def make_archive_and_repopulate(self):
        """
        Make an archive of a VASP calculation and copy back over relevant files
        """
        with change_directory(self.calc_path):
            contcar_path = Path("CONTCAR")
            contcar_exists = contcar_path.exists()
            if contcar_exists:
                contcar_is_empty = contcar_path.stat().st_size == 0
            else:
                contcar_is_empty = True

            # if CONTCAR is empty, don't make an archive and clean up
            if contcar_is_empty:
                all_files = [
                    f
                    for f in Path(".").glob("*")
                    if f.is_file() and "archive" not in f.name and "json" not in f.name
                ]
                for f in all_files:
                    os.remove(f)
            # else, make the archive
            else:
                num_previous_archives = len(list(Path(".").glob("archive*")))
                archive_name = Path(f"archive_{num_previous_archives}")
                logger.info(f"Making {archive_name}...")
                archive_name.mkdir()

                all_files = [
                    f
                    for f in Path(".").glob("*")
                    if f.is_file() and "archive" not in f.name and "json" not in f.name
                ]
                for f in all_files:
                    # add if symlink for testing
                    if f.is_symlink():
                        f_links_to = f.readlink()
                        os.remove(f)
                        shutil.copy2(f_links_to, f)
                    shutil.move(f, archive_name)

        self.create()

    def create(self):
        """
        Make VASP input files

        Don't touch the order! make_incar and make_vaspq rely on the poscar and
        potcar existing already
        """
        if not self.calc_path.exists():
            self.calc_path.mkdir()
        self.make_poscar()
        self.make_potcar()
        if "write_kpoints" in self.calc_config:
            if self.calc_config["write_kpoints"]:
                self.make_kpoints()
        self.make_incar()
        self.make_vaspq()
