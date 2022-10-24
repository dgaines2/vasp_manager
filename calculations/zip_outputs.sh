#!/bin/bash

find * -name OUTCAR -exec zstd --rm -z {} \;
find * -name vasprun.xml -exec zstd --rm -z {} \;
