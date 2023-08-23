#!/bin/bash

archive_dirs=$(find * -name 'archive*.tar.gz')
for archive_dir in $archive_dirs; do
    parent_dir=$(echo $archive_dir | cut -d '/' -f 1-2)
    archive_name=$(echo $archive_dir | rev | cut -d '/' -f 1 | rev)
    echo $archive_dir $parent_dir $archive_name
    ( cd $parent_dir &&
        tar -xzvf $archive_name &&
        rm $archive_name
    )
done

find * -name OUTCAR.gz -exec gunzip -v {} \;
find * -name vasprun.xml.gz -exec gunzip -v {} \;
find * -name EIGENVAL.gz -exec gunzip -v {} \;
find * -name PROCAR -exec gunzip -v {} \;
find * -name vaspout.h5 -exec gunzip -v {} \;
