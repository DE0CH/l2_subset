#!/bin/bash -e

make -C build CFLAGS="-Wall -Werror -g -Og"

valgrind --leak-check=full \
         --show-leak-kinds=all \
         --track-origins=yes \
         --verbose \
build/l2_subset build/points.txt 4
cd ..
