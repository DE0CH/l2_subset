import argparse
from itertools import product
from subprocess import run
from os import mkdir
import json

parser = argparse.ArgumentParser(description='Run matrix configuration')
parser.add_argument('-N', nargs='+', type=int, help='N: the size of the selection set')
parser.add_argument('-M', nargs='+', type=int, help='M: the size of the subset')
parser.add_argument('-DIM', nargs='+', type=int, help='DIM: the dimension')
substitution = parser.parse_args()

build_dir = 'build'
output_dir = 'output'

mkdir('output')
mkdir('output/job-scripts')

for N, M, DIM in product(substitution.N, substitution.M, substitution.DIM):
    substitution = {'N': N, 'M': M, 'DIM': DIM, 'BUILD_DIR': build_dir, 'OUTPUT_DIR': f"{output_dir}/output-{N}-{M}-{DIM}"}
    run(['jinja2', 'run-single.slurm'], input=json.dumps(substitution).encode('utf-8'), stdout=open(f'output/job-scripts/run-single-{N}-{M}-{DIM}.slurm', 'w'), check=True)
    run(['chmod', '+x', f'output/job-scripts/run-single-{N}-{M}-{DIM}.slurm'], check=True)
    run([f'output/job-scripts/run-single-{N}-{M}-{DIM}.slurm'], check=True)
    run(['sbatch', f'output/job-scripts/run-single-{N}-{M}-{DIM}.slurm'], check=True)
