#!/bin/bash

files=(
    OUTCAR
    vasprun.xml
    EIGENVAL
    DOSCAR
    PROCAR
    vaspout.h5
)

for name in "${files[@]}"; do
    find . -name "$name" -exec gunzip -v {} \;
done
