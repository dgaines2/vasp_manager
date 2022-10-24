#!/bin/bash

for calc_dir in */; do
    echo $calc_dir
    tar -cvf - $calc_dir | zstd --rm -z > "${calc_dir%/}.tar.zst"
done
