#!/usr/bin/env python3
"""Decompress VASP output files.

Run from within a calculations/ directory.

Performs two operations:
  1. Extract any archive_N.tar.gz archives back into archive_N/ directories
  2. Decompress individual output files (OUTCAR, vasprun.xml, etc.) in-place
"""

import gzip
import shutil
import tarfile
from pathlib import Path

OUTPUT_FILES = [
    "OUTCAR",
    "vasprun.xml",
    "EIGENVAL",
    "DOSCAR",
    "PROCAR",
    "vaspout.h5",
]


def gunzip_file(gz_path: Path) -> None:
    out_path = gz_path.parent / gz_path.stem  # strips .gz
    with gzip.open(gz_path, "rb") as f_in, open(out_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    gz_path.unlink()
    print(f"  {gz_path} -> {out_path.name}")


def extract_archive_tar(tar_path: Path) -> None:
    with tarfile.open(tar_path, "r:gz") as tar:
        tar.extractall(path=tar_path.parent)
    tar_path.unlink()
    archive_name = tar_path.name.removesuffix(".tar.gz")
    print(f"  {tar_path.name} -> {tar_path.parent / archive_name}/")


root = Path(".")

print("Extracting archive directories...")
for tar_path in sorted(root.rglob("archive_*.tar.gz")):
    extract_archive_tar(tar_path)

print("Decompressing output files...")
for filename in OUTPUT_FILES:
    for gz_path in sorted(root.rglob(f"{filename}.gz")):
        gunzip_file(gz_path)
