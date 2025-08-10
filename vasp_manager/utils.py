# Copyright (c) Dale Gaines II
# Distributed under the terms of the MIT LICENSE

from __future__ import annotations

import gzip
import json
import logging
import os
from collections import deque
from contextlib import contextmanager
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

if TYPE_CHECKING:
    from vasp_manager.types import Filepath


@contextmanager
def change_directory(new_dir: str | Path):
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


class LoggerAdapter(logging.LoggerAdapter):
    """Logging adapter to add a custom prefix to a Logger"""

    def __init__(
        self,
        logger: logging.Logger,
        prefix: str,
        separator: str = " -- ",
    ):
        super(LoggerAdapter, self).__init__(logger, {})
        self.prefix = prefix
        self.separator = separator

    def process(self, msg, kwargs):
        return f"{self.prefix}{self.separator}{msg}", kwargs


def get_pmg_structure_from_poscar(
    poscar_path: Filepath,
    to_process: bool = True,
    primitive: bool = True,
    symprec: float = 1e-03,
    angle_tolerance: float = -1.0,
    international_monoclinic: bool = False,
    return_spacegroup: bool = False,
) -> Structure | tuple[Structure, int]:
    """
    Args:
        poscar_path: path to POSCAR file
        to_process: if True, get standard reduced structure
        primitive: if True, get primitive structure, else get
            conventional structure
        symprec: symprec for SpacegroupAnalyzer
        angle_tolerance: angle tolerance for SpacegroupAnalyzer
        international_monoclinic: if True, convert to proper
            international convention such that beta is the non-right angle.
            WARNING: setting True is not compatible with pymatgen kpaths
        return_spacegroup: if True, return spacegroup number

    Returns:
        structure: structure from POSCAR
    """
    structure = Structure.from_file(poscar_path)
    if to_process:
        sga = SpacegroupAnalyzer(
            structure, symprec=symprec, angle_tolerance=angle_tolerance
        )
        if primitive:
            structure = sga.get_primitive_standard_structure(
                international_monoclinic=international_monoclinic
            )
        else:
            structure = sga.get_conventional_standard_structure(
                international_monoclinic=international_monoclinic
            )
        if return_spacegroup:
            spacegroup = sga.get_space_group_number()
            return structure, spacegroup
    return structure


def pcat(file_names: Filepath | list[Filepath]) -> str:
    """
    Custom python-only replacement for cat

    Args:
        file_names: names of files to cat together

    Returns:
        catted:
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
    file_name: Filepath,
    str_to_grep: str,
    stop_after_first_match: bool = False,
    after: int | None = None,
    as_string: bool = False,
) -> str | list[str]:
    """
    Custom python-only replacement for grep

    Args:
        file_name: path of file
        str_to_grep: target string
        stop_after_first_match: if True, stop after first found instance of
            str_to_grep
        after: if not None, return {after} lines found after str_to_grep
        as_string: if as_string, return a single string, else return splitlines

    Returns:
        matches:
    """
    opener = gzip.open if ".gz" in str(file_name) else open
    matches = []
    line_idx_to_include: set[int] = set()
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
        matches_as_str = "\n".join([line for line in matches])
    return matches_as_str if as_string else matches


def phead(
    file_name: Filepath,
    n_head: int = 1,
    as_string: bool = False,
) -> str | list[str]:
    """
    Custom python-only replacement for head

    Args:
        file_name: path of file
        n_head: n lines to head
        as_string: if as_string, return a single string, else return splitlines

    Returns:
        head:
    """
    opener = gzip.open if ".gz" in str(file_name) else open
    head = []
    with opener(file_name, "rt") as fr:
        for i, line in enumerate(fr):
            if i < n_head:
                head.append(line.strip("\n"))
    if as_string:
        head_as_str = "\n".join([line for line in head])
    return head_as_str if as_string else head


def ptail(
    file_name: Filepath,
    n_tail: int = 1,
    as_string: bool = False,
) -> str | list[str]:
    """
    Custom python-only replacement for grep

    Args:
        file_name: path of file
        n_tail: n lines to tail
        as_string: if as_string, return a single string, else return splitlines

    Returns:
        tail:
    """
    opener = gzip.open if ".gz" in str(file_name) else open
    tail: deque[str] = deque(maxlen=n_tail)
    with opener(file_name, "rt") as fr:
        for line in fr:
            tail.append(line.strip("\n"))
    if as_string:
        tail_as_str = "\n".join([line for line in tail])
    return tail_as_str if as_string else list(tail)


def make_potcar_anonymous(
    input_file_name: Filepath,
    output_file_name: Filepath | None = None,
) -> None:
    """
    Replace full POTCAR with only single POTCAR names

    Args:
        input_file_name: path of POTCAR file
        output_file_name: path to write anonymized POTCAR. If None, write to
            the location of input_f_name
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
