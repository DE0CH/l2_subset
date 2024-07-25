#!/bin/bash

set -e
set -eo pipefail

gcc -O3 -o l2_subset l2_subset.c dem_disc.c -lm
./l2_subset large_points.txt 2 1 1

gcc -O2 -o l2_subset l2_subset.c dem_disc.c -lm
./l2_subset large_points.txt 2 1 1