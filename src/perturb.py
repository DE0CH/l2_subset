import struct
import argparse
import subprocess
import random
import re
import sys

def read_numbers_from_file(filename):
    with open(filename, "rb") as file:  # Open the file in binary mode
        # Read the first 24 bytes (3 long long values, 8 bytes each)
        data = file.read(24)
        
        if len(data) < 24:
            raise ValueError("File is too small to contain 3 long long numbers.")
        
        # Unpack the data as three long long integers ('q' format specifier)
        n, m, d = struct.unpack('qqq', data)
        
        return n, m, d

parser = argparse.ArgumentParser()
parser.add_argument('point_file', type=argparse.FileType('rb'), help='Point file')
parser.add_argument('seed', type=int, help='Seed for random number generator')
parser.add_argument('p', type=int, help="number of points to pertube")
parser.add_argument('iterations', type=int, help='Number of iterations')
args = parser.parse_args()

random.seed(args.seed)

print(f"Running with local perturbation. Seed: {args.seed}, p: {args.p}, Iterations: {args.iterations}")

best_l2 = 1
best_linf = 1
n, m, d = read_numbers_from_file(args.point_file.name)
points = random.sample(list(range(n)), m)
best_points = points

for i in range(args.iterations):
    lines = []
    print("Iteration", i + 1)
    p = subprocess.Popen(["./l2_subset_from_compiled_matrix_w_starting_point", args.point_file.name, str(random.randrange(0, 2**61-1)), *map(str, points)], stdout=subprocess.PIPE)
    for line in p.stdout:
        line = line.decode('utf-8')
        print(line, end='')
        matches = re.match(r"Active points: ([\d\s]+)", line)
        if matches:
            new_points = list(map(int, matches.group(1).split()))
        matches = re.match(r"Active point sum: -?([\d\.]+)", line)
        if matches:
            new_l2 = float(matches.group(1))
        matches = re.match(r"linf discrepancy: -?([\d\.]+)", line)
        if matches:
            new_linf = float(matches.group(1))
        lines.append(line)
    p.wait()
    # check return code
    if p.returncode != 0:
        raise subprocess.CalledProcessError(p.returncode, p.args)

    print("perturbing for next iteration")
    if new_linf < best_linf:
        best_l2 = new_l2
        best_linf = new_linf
        best_points = new_points
    points = best_points.copy()
    other_points = []
    points_set = set(points)
    for i in range(n):
        if i not in points_set:
            other_points.append(i)
    for i in range(args.p):
        old = random.randrange(0, m)
        new = random.randrange(0, n-m)
        t = points[old]
        points[old] = other_points[new]
        other_points[new] = t

