#!/bin/bash -e

set -eo pipefail

make -C build clean
make -C build CFLAGS="-DDEBUG_SLOW -Wall -O3"
N=20
M=10
TRIALS=10
echo "Generating $N points"
build/gen_points 3 build/large_points.txt "$N"

build/l2_subset build/large_points.txt "$M" 42 "$TRIALS" | tee build/l2_subset_log.txt
cat build/l2_subset_log.txt | grep "Active points:" | sed 's/Active points: //g' > build/selected_points.txt
while read line; do
    echo "Verifying descrepancy of selected points with python script"
    python3 formula.py build/large_points.txt <(echo "$line")
done <build/selected_points.txt


