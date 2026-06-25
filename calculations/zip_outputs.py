#!/usr/bin/env python3
"""Compress VASP output files and archive directories.

Run from within a calculations/ directory.

Performs two operations:
  1. gzip individual output files (OUTCAR, vasprun.xml, etc.) in-place
  2. Compress any archive_N/ directories into archive_N.tar.gz and remove the originals
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


def gzip_file(path: Path) -> None:
    gz_path = path.parent / (path.name + ".gz")
    with open(path, "rb") as f_in, gzip.open(gz_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    path.unlink()
    print(f"  {path} -> {gz_path.name}")


def compress_archive_dir(archive_dir: Path) -> None:
    tar_path = archive_dir.parent / f"{archive_dir.name}.tar.gz"
    with tarfile.open(tar_path, "w:gz") as tar:
        tar.add(archive_dir, arcname=archive_dir.name)
    shutil.rmtree(archive_dir)
    print(f"  {archive_dir}/ -> {tar_path.name}")


root = Path(".")

print("Compressing output files...")
for filename in OUTPUT_FILES:
    for path in sorted(root.rglob(filename)):
        gzip_file(path)

print("Compressing archive directories...")
for archive_dir in sorted(root.rglob("archive_*")):
    if archive_dir.is_dir():
        compress_archive_dir(archive_dir)
