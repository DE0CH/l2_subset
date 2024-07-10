#!/bin/bash

gcc -O3 -DCOMPUTE_MODE=2 -o l2_subset l2_subset.c dem_disc.c -lm
./l2_subset large_points.txt 2 1 1