#!/bin/bash

compressed_calc_dirs=$(find * -name '*.tar.gz')
for compressed_calc_dir in $compressed_calc_dirs; do
    echo $compressed_calc_dir
    tar -xzvf $compressed_calc_dir && rm $compressed_calc_dir
done
