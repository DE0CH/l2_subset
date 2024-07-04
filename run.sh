#!/bin/bash -e

make -C build CFLAGS="-UDEBUG_SLOW -Wall -O3"

build/gen_points 3 build/large_points.txt 10000
build/l2_subset build/large_points.txt 10 | tee build/l2_subset_log.txt
cat build/l2_subset_log.txt | grep "Active points:" | sed 's/Active points: //g' > build/selected_points.txt
echo "Verifying descrepancy of selected points with python script"
python3 formula.py build/large_points.txt build/selected_points.txt

