#!/bin/bash -e

make -C build CFLAGS="-DDEBUG_SLOW=1 -Wall -g -O3"

# build/gen_points 9 build/large_points.txt 50000
build/l2_subset build/large_points.txt 1000
cd ..
