#!/bin/bash

for f in *.tar.zst; do
    zstd --rm -d $f 
    tar_name="${f%.zst}"
    tar -xvf $tar_name
    rm $tar_name
done
