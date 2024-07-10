#!/bin/bash

set -e
set -eo pipefail

range=$(seq 1 100)
export SCRIPT="./run-single.slurm"
run_script() {
    SLURM_ARRAY_TASK_ID=$1 $SCRIPT
}

export -f run_script
parallel -j $(($(nproc) - 4)) run_script ::: $range
