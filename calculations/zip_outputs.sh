#!/bin/bash

find * -name OUTCAR -exec gzip -v {} \;
find * -name vasprun.xml -exec gzip -v {} \;
find * -name EIGENVAL -exec gzip -v {} \;

archive_dirs=$(find * -type d -name archive*)
for archive_dir in $archive_dirs; do
    parent_dir=$(echo $archive_dir | cut -d '/' -f 1-2)
    archive_name=$(echo $archive_dir | rev | cut -d '/' -f 1 | rev)
    echo $parent_dir $archive_name
    ( cd $parent_dir &&
        tar -cvzf $archive_name.tar.gz $archive_name && rm -r $archive_name
    )
done
