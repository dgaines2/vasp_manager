# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

import gzip
import json
import os
from collections import deque
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
    symprec=1e-03,
    return_spacegroup=False,
):
    """
    Args:
        poscar_path (str | Path)
        to_process (bool): if True, get standard reduced structure
        primitive (bool): if True, get primitive structure, else get
            conventional structure
        symprec (float): symprec for SpacegroupAnalyzer
        return_spacegroup (bool): if True, return spacegroup number
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
        if return_spacegroup:
            spacegroup = sga.get_space_group_number()
            return structure, spacegroup
    return structure


def pcat(file_names):
    """
    Custom python-only replacement for cat

    Args:
        file_names (list): names of files to cat together
    Returns:
        catted (str)
    """
    file_contents = []
    if isinstance(file_names, (str, Path)):
        file_names = [file_names]
    for file_name in file_names:
        with open(file_name) as fr:
            file_content = fr.read()
        file_contents.append(file_content)
    catted = "".join(file_content for file_content in file_contents)
    return catted


def pgrep(
    file_name,
    str_to_grep,
    stop_after_first_match=False,
    after=None,
    as_string=False,
):
    """
    Custom python-only replacement for grep

    Args:
        file_name (str | Path): path of file
        str_to_grep (str): target string
        stop_after_first_match (bool): if True, stop after first found instance of
            str_to_grep
        after (int): if not None, return {after} lines found after str_to_grep
        as_str (bool): if as_string, return a single string, else return splitlines
    Returns:
        matches (str | list)
    """
    opener = gzip.open if ".gz" in str(file_name) else open
    matches = []
    line_idx_to_include = set()
    with opener(file_name, "rt") as fr:
        for line_idx, line in enumerate(fr):
            if str_to_grep in line:
                matches.append(line.strip("\n"))
                if after is not None:
                    line_idx_to_include.update(range(line_idx + 1, line_idx + after + 1))
            if line_idx in line_idx_to_include:
                matches.append(line.strip("\n"))
                line_idx_to_include.discard(line_idx)
            if stop_after_first_match:
                if len(matches) != 0 and len(line_idx_to_include) == 0:
                    break
    if as_string:
        matches = "\n".join([line for line in matches])
    return matches


def phead(file_name, n_head=1, as_string=False):
    """
    Custom python-only replacement for head

    Args:
        file_name (str): path of file
        n_head (int): n lines to head
        as_str (bool): if as_string, return a single string, else return splitlines
    Returns:
        head (str | list)
    """
    opener = gzip.open if ".gz" in str(file_name) else open
    head = []
    with opener(file_name, "rt") as fr:
        for i, line in enumerate(fr):
            if i < n_head:
                head.append(line.strip("\n"))
    if as_string:
        head = "\n".join([line for line in head])
    return head


def ptail(file_name, n_tail=1, as_string=False):
    """
    Custom python-only replacement for grep

    Args:
        file_name (str): path of file
        n_tail (int): n lines to tail
        as_str (bool): if as_string, return a single string, else return splitlines
    Returns:
        tail (str | list)
    """
    opener = gzip.open if ".gz" in str(file_name) else open
    tail = deque(maxlen=n_tail)
    with opener(file_name, "rt") as fr:
        for line in fr:
            tail.append(line.strip("\n"))
    tail = list(tail)
    if as_string:
        tail = "\n".join([line for line in tail])
    return tail


def make_potcar_anonymous(input_file_name, output_file_name=None):
    """
    Replace full POTCAR with only single POTCAR names

    Args:
        input_file_name (str | Path): path of POTCAR file
        output_file_name (str | Path | None): path to write anonymized POTCAR
            if None, write to the location of input_f_name
    Returns:
        None
    """
    if output_file_name is None:
        output_file_name = input_file_name

    with open(input_file_name, "rt") as fr:
        full_potcar_text = [line.strip() for line in fr.readlines()]
    trimmed_potcar_lines = []
    for line in full_potcar_text:
        if "TITEL" in line:
            trimmed_potcar_lines.append(line.split("=")[1].strip())
    trimmed_potcar_string = "\n".join([line for line in trimmed_potcar_lines])
    with open(output_file_name, "w+") as fw:
        fw.write(trimmed_potcar_string)
