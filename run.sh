#!/bin/bash -e

set -eo pipefail

make -C build clean
make -C build CFLAGS="-UDEBUG_SLOW -DCOMPUTE_MODE=2 -Wall -O3"
N=100000
M=50
DIM=3
TRIALS=100
echo "Generating $N points"
build/gen_points "$DIM" build/large_points.txt "$N"

echo "Compiling matrix"
build/l2_subset_compile_matrix build/large_points.txt "$M" build/large_points.p
for i in $(seq 1 "$TRIALS"); do
    TEMP=`mktemp`
    echo "Starting trial $i"
    build/l2_subset_from_compiled_matrix build/large_points.p $i | tee build/l2_subset_log.txt
    cat build/l2_subset_log.txt | grep "Active points:" | sed 's/Active points: //g' > build/selected_points.txt
    echo "Verifying trial $i"
    python3 src/formula.py build/large_points.txt build/selected_points.txt
    python3 src/filter_point.py build/large_points.txt build/selected_points.txt > build/selected_pf.txt
    build/linf_subset build/selected_pf.txt "$DIM" "$M" "$M"  > $TEMP 2>&1 || true # There is some segfault at the end of the program, and we don't really care. 
    echo -n "L infinity: "
    grep -oP "Final discre:\K[0-9.]+" $TEMP 
    rm $TEMP
done
