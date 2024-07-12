#!/bin/bash
#SBATCH -n 1
#SBATCH --array=1-100
#SBATCH --time=24:00:00
#SBATCH -o {{ OUTPUT_DIR }}/%a.out
#SBATCH -e {{ OUTPUT_DIR }}/%a.err

set -e
set -eo pipefail

export N={{ N }}
export M={{ M }}
export DIM={{ DIM }}
export BUILD_DIR={{ BUILD_DIR }}
export OUTPUT_DIR={{ OUTPUT_DIR }}

generate_points() {
    mkdir "$OUTPUT_DIR"
    echo "Generating $N points"
    "$BUILD_DIR/gen_points" "$DIM" "$OUTPUT_DIR/points.txt" "$N"

    echo "Compiling matrix"
    "$BUILD_DIR/l2_subset_compile_matrix" "$OUTPUT_DIR/points.txt" "$M" "$OUTPUT_DIR/points.p"
    echo "finished compiling matrix"
}

# Function to execute a single trial
run_trial() {
    trial=$1
    echo "Starting trial $trial"
    "$BUILD_DIR/l2_subset_from_compiled_matrix" "$OUTPUT_DIR/points.p" "$trial" 2>&1
    echo "End of trial $trial"
}
export -f run_trial

if [ -z ${SLURM_ARRAY_TASK_ID+x} ]; then
    generate_points
else
    run_trial $SLURM_ARRAY_TASK_ID
fi
