#!/bin/bash

find * -name OUTCAR.zst -exec zstd --rm -d {} \;
find * -name vasprun.xml.zst -exec zstd --rm -d {} \;
