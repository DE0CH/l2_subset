#!/bin/bash -e

set -eo pipefail

make -C build clean
make -C build CFLAGS="-UDEBUG_SLOW -DCOMPUTE_MODE=2 -Wall -O3 -g"
export N=500000
export M=50
export DIM=3
TRIALS=300
echo "Generating $N points"
build/gen_points "$DIM" build/large_points.txt "$N"

echo "Compiling matrix"
build/l2_subset_compile_matrix build/large_points.txt "$M" build/large_points.p
echo "finished compiling matrix"
# Function to execute a single trial
run_trial() {
    set -eo pipefail  # Ensure error handling is applied within the function
    
    local trial=$1
    local temp_dir=$(mktemp -d)
    local output_file=$(mktemp)
    $(
    {
        echo "Starting trial $trial"
        build/l2_subset_from_compiled_matrix build/large_points.p "$trial" | tee "$temp_dir/l2_subset_log.txt"
        grep "Active points:" "$temp_dir/l2_subset_log.txt" | sed 's/Active points: //g' > "$temp_dir/selected_points.txt"
        echo "Verifying trial $trial"
        python3 src/formula.py build/large_points.txt "$temp_dir/selected_points.txt"
        python3 src/filter_point.py build/large_points.txt "$temp_dir/selected_points.txt" > "$temp_dir/selected_pf.txt"
        build/linf_subset "$temp_dir/selected_pf.txt" "$DIM" "$M" "$M" > "$temp_dir/output.txt" 2>&1 || true # There is some segfault at the end of the program, and we don't really care.
        echo -n "Trial $trial: L infinity: "
        grep -oP "Final discre:\K[0-9.]+" "$temp_dir/output.txt" | head -1
        rm -r "$temp_dir"
    } &> "$output_file"
    )
    cat "$output_file"
    rm "$output_file"
}
export -f run_trial

# Run the trials in parallel
parallel -j $(($(nproc) - 1))  run_trial ::: $(seq 1 "$TRIALS")
