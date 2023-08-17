#!/bin/bash

calc_dirs=$(find * -type d -maxdepth 0)
for calc_dir in $calc_dirs; do
    echo $calc_dir
    tar -czvf $calc_dir.tar.gz $calc_dir && rm -r $calc_dir
done
