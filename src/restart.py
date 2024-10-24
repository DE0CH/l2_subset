import argparse
import subprocess
import tempfile
from os.path import join
from os import rmdir, mkdir
import random
import sys
import shutil

parser = argparse.ArgumentParser()
parser.add_argument('n', type=int, help='n')
parser.add_argument('m', type=int, help='m')
parser.add_argument('d', type=int, help='d')
parser.add_argument('alpha', type=float, help='alpha')
parser.add_argument('p', type=int, help='p')
parser.add_argument('number_local_restart', type=int, help='Number of restarts')
parser.add_argument('num_global_restart', type=int, help='Number of global restarts')

parser.add_argument('seed', type=int, help='Seed for random number generator')
args = parser.parse_args()

random.seed(args.seed)

# make a temporary directory

print(f'Running a new restart with seed {args.seed}, n = {args.n}, m = {args.m}, d = {args.d}, p = {args.p}')
temp_dir = tempfile.TemporaryDirectory()
temp_dir = temp_dir.name
mkdir(temp_dir)

try:
    print(f"generating point file with {args.n} points")
    subprocess.run(['./gen_points', str(args.d), join(temp_dir, 'points.txt'), str(args.n)], check=True)
    print(f"Compiling matrix")
    subprocess.run(['./l2_subset_compile_matrix', join(temp_dir, 'points.txt'), str(args.m), str(args.alpha), join(temp_dir, 'points.p')], check=True)
    for i in range(args.num_global_restart):
        subprocess.run([sys.executable, '-u', '../src/perturb.py', join(temp_dir, 'points.p'), str(random.randrange(0, 2**63)), str(args.p), str(args.number_local_restart)], check=True)
finally:
    print("cleaning up temporary directory")
    shutil.rmtree(temp_dir)