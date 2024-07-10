#!/bin/bash
#SBATCH -n 1
#SBATCH --array=1-100
#SBATCH --time=24:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null

set -eo pipefail

export N=100
export M=50
export DIM=3

generate_points() {
    echo "" > output.txt
    make -C build clean
    make -C build
    mkdir -p output
    echo "Generating $N points"
    build/gen_points "$DIM" build/large_points.txt "$N"

    echo "Compiling matrix"
    build/l2_subset_compile_matrix build/large_points.txt "$M" build/large_points.p
    echo "finished compiling matrix"
}

# Function to execute a single trial
run_trial() {
    set -e
    set -eo pipefail  # Ensure error handling is applied within the function
    local trial=$1
    local temp_dir=$(mktemp -d)
    $(
    {
        echo "Starting trial $trial"
        build/l2_subset_from_compiled_matrix build/large_points.p "$trial" | tee "$temp_dir/l2_subset_log.txt"
        grep "Active points:" "$temp_dir/l2_subset_log.txt" | sed 's/Active points: //g' > "$temp_dir/selected_points.txt"
        python3 src/filter_point.py build/large_points.txt "$temp_dir/selected_points.txt" > "$temp_dir/selected_pf.txt"
        echo -n "python code evaluation: "
        echo -n "Trial $trial: L infinity: "
        build/linf_disc "$temp_dir/selected_pf.txt"
        rm -r "$temp_dir"
    } > "output/$trial.txt"
    )
}
export -f run_trial

if [ -z ${SLURM_ARRAY_TASK_ID+x} ]; then
    generate_points
else
    run_trial $SLURM_ARRAY_TASK_ID
fi
