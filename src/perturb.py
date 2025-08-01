import struct
import argparse
import subprocess
import random
import re
import sys
from utils import read_matrix_from_binary, read_points_from_bianry
import numpy as np
from code_snippet_FC import KSD_loss_RBF
import torch
import threading
import time

def read_numbers_from_file(filename):
    with open(filename, "rb") as file:  # Open the file in binary mode
        # Read the first 24 bytes (3 long long values, 8 bytes each)
        data = file.read(24)
        
        if len(data) < 24:
            raise ValueError("File is too small to contain 3 long long numbers.")
        
        # Unpack the data as three long long integers ('q' format specifier)
        n, m, d = struct.unpack('qqq', data)
        
        return n, m, d
    
def calc_linf(points, subset):
    sub_points = points[subset]
    v = KSD_loss_RBF(sub_points, 1, len(subset), 2)
    return v.mean(dim=(1, 2)).mean().item()
def timer():
    for i in range(1, 4):
        print(f"{i} second(s) elapsed")
        time.sleep(1)
    

parser = argparse.ArgumentParser()
parser.add_argument('compiled_point_file', type=argparse.FileType('r'), help='Compiled Point file')
parser.add_argument('raw_point_file', type=argparse.FileType('r'), help='Raw Point file')
parser.add_argument('seed', type=int, help='Seed for random number generator')
parser.add_argument('p', type=int, help="number of points to perturb")
parser.add_argument('iterations', type=int, help='Number of iterations')
parser.add_argument('initial_population_size', type=int, help='Initial population size')
args = parser.parse_args()

random.seed(args.seed)

print(f"Running with local perturbation. Seed: {args.seed}, p: {args.p}, Iterations: {args.iterations}")

best_l2 = 1
best_linf = 1
n, m, d = read_numbers_from_file(args.compiled_point_file.name)

matrix = read_matrix_from_binary(args.compiled_point_file.name)
raw_points = torch.tensor(read_points_from_bianry(args.raw_point_file.name, n), dtype=torch.float32)

best_linf = float('inf')
for i in range(args.initial_population_size):
    print("Calculating initial population", i + 1)
    points = random.sample(list(range(n)), m)
    ans = calc_linf(raw_points, points)
    if ans < best_linf:
        best_linf = ans
        best_points = points

points = best_points
print("initial best points:", *points)

for i in range(args.iterations):
    lines = []
    print("Iteration", i + 1)
    should_raise_keyboard_interrupt = False
    should_raise_none_zero_return_code = False
    try:
        p = subprocess.Popen(["./l2_subset_from_compiled_matrix_w_starting_point", args.compiled_point_file.name, str(random.randrange(0, 2**61-1)), *map(str, points)], stdout=subprocess.PIPE)
        for line in p.stdout:
            line = line.decode('utf-8')
            print(line, end='')
            matches = re.match(r"Active points: ([\d\s]+)", line)
            if matches:
                new_points = list(map(int, matches.group(1).split()))
            matches = re.match(r"Active point sum: -?([\d\.]+)", line)
            if matches:
                new_l2 = float(matches.group(1))
            lines.append(line)
    except KeyboardInterrupt:
        for line in p.stdout:
            line = line.decode('utf-8')
            print(line, end='')
            matches = re.match(r"Active points: ([\d\s]+)", line)
            if matches:
                new_points = list(map(int, matches.group(1).split()))
            matches = re.match(r"Active point sum: -?([\d\.]+)", line)
            if matches:
                new_l2 = float(matches.group(1))
            lines.append(line)
        should_raise_keyboard_interrupt = True
    finally:
        p.wait()
        if p.returncode != 0:
            should_raise_none_zero_return_code = True

    if should_raise_keyboard_interrupt:
        print("KeyboardInterrupt raised, continuing regardless with errors to attempt to calculate linf discrepancy")
    if should_raise_none_zero_return_code:
        print("Non-zero return code raised, continuing regardless with errors to attempt to calculate linf discrepancy")
    if should_raise_keyboard_interrupt or should_raise_none_zero_return_code:
        print("Starting a timer (it prints out a number every second, if it disappears without any more output, it means that the program has been forcefully exited, most likely due to a timeout on the cluster)")
        timer_thread = threading.Thread(target=timer)
        timer_thread.start()
    new_linf = calc_linf(raw_points, new_points)

    print("linf discrepancy:", new_linf)
    if new_linf < best_linf:
        best_l2 = new_l2
        best_linf = new_linf
        best_points = new_points
    if should_raise_keyboard_interrupt or should_raise_none_zero_return_code:
        break
    print("perturbing for next iteration")
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

print(f"Index of best points: {' '.join(map(str, best_points))}")
print(f"linf discrepancy: {best_linf}")
print(f"======= points below ======")

for p in raw_points[best_points]:
    print(*map(lambda x: f"{x:.8g}", p.numpy()))

if should_raise_keyboard_interrupt:
    print("KeyboardInterrupt raised, exiting with errors")
    raise KeyboardInterrupt

if should_raise_none_zero_return_code:
    print("Non-zero return code raised, exiting with errors")
    raise subprocess.CalledProcessError(p.returncode, p.args)
