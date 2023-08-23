# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import gzip
import json
import os
from contextlib import contextmanager
from pathlib import Path

import numpy as np
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


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
    structure = Structure.from_file(poscar_path)
    if to_process:
        sga = SpacegroupAnalyzer(structure, symprec=symprec, angle_tolerance=-1.0)
        if primitive:
            structure = sga.get_primitive_standard_structure()
        else:
            structure = sga.get_conventional_standard_structure()
        if return_sg:
            sg = sga.get_space_group_number()
            return structure, sg
    return structure


def pcat(f_names):
    """
    Custom python-only replacement for cat

    Args:
        f_names (list): names of files to cat together
    Returns:
        catted (str)
    """
    f_contents = []
    if isinstance(f_names, (str, Path)):
        f_names = [f_names]
    for f_name in f_names:
        with open(f_name) as fr:
            f_content = fr.read()
        f_contents.append(f_content)
    catted = "".join(f_content for f_content in f_contents)
    return catted


def pgrep(
    f_name,
    str_to_grep,
    stop_after_first_match=False,
    after=None,
    as_string=False,
):
    """
    Custom python-only replacement for grep

    Args:
        f_name (str): path of file
        str_to_grep (str): target string
        stop_after_first_match (bool): if True, stop after first found instance of
            str_to_grep
        after (int): if not None, return {after} lines found after str_to_grep
        as_str (bool): if as_string, return a single string, else return splitlines
    Returns:
        matches (str | list)
    """
    opener = gzip.open if ".gz" in str(f_name) else open
    with opener(f_name, "rt") as fr:
        f_lines = [line.strip() for line in fr.readlines()]
    matches = []
    for i, line in enumerate(f_lines):
        if str_to_grep in line:
            matches.append(line)
            if after is not None:
                matches.extend(f_lines[(i + 1) : (i + after + 1)])
            if stop_after_first_match:
                break
    if as_string:
        matches = "\n".join([line for line in matches])
    return matches


def phead(f_name, n_head=1, as_string=False):
    """
    Custom python-only replacement for head

    Args:
        f_name (str): path of file
        n_tail (int): n lines to head
        as_str (bool): if as_string, return a single string, else return splitlines
    Returns:
        head (str | list)
    """
    opener = gzip.open if ".gz" in str(f_name) else open
    with opener(f_name, "rt") as fr:
        head = [line.strip() for line in fr.readlines()[:n_head]]
    if as_string:
        head = "\n".join([line for line in head])
    return head


def ptail(f_name, n_tail=1, as_string=False):
    """
    Custom python-only replacement for grep

    Args:
        f_name (str): path of file
        n_tail (int): n lines to tail
        as_str (bool): if as_string, return a single string, else return splitlines
    Returns:
        tail (str | list)
    """
    opener = gzip.open if ".gz" in str(f_name) else open
    with opener(f_name, "rt") as fr:
        tail = [line.strip() for line in fr.readlines()[-n_tail:]]
    if as_string:
        tail = "\n".join([line for line in tail])
    return tail


def make_potcar_anonymous(input_f_name, output_f_name=None):
    """
    Replace full POTCAR with only single POTCAR names

    Args:
        input_f_name (str): path of POTCAR file
        output_f_name (str): path to write anonymized POTCAR
    Returns:
        None
    """
    if output_f_name is None:
        output_f_name = input_f_name

    with open(input_f_name, "rt") as fr:
        full_potcar_text = [line.strip() for line in fr.readlines()]
    trimmed_potcar_lines = []
    for line in full_potcar_text:
        if "TITEL" in line:
            trimmed_potcar_lines.append(line.split("=")[1].strip())
    trimmed_potcar_string = "\n".join([line for line in trimmed_potcar_lines])
    with open(output_f_name, "w+") as fw:
        fw.write(trimmed_potcar_string)
