import argparse
import subprocess
import tempfile
from os.path import join
from os import mkdir
import os
import random
import sys
import shutil
import signal
import time
import threading

def send_sigint_to_process_group_after_delay(delay, event):
    if not event.wait(delay):
        print("Timeout reached")
        os.killpg(os.getpgid(os.getpid()), signal.SIGINT)

parser = argparse.ArgumentParser()
parser.add_argument('n', type=int, help='n (we select m points from n points)')
parser.add_argument('m', type=int, help='m (we select m points from n points)')
parser.add_argument('d', type=int, help='dimension')
parser.add_argument('p', type=int, help='the number of points to perturb. It shouldn\'t be larger than m. Try 4 first (for no reason other that I like the number 4).')
parser.add_argument('number_local_restart', type=int, help='Number of local restarts. Between two local restart, we remember the best subset and do perturbation.')
parser.add_argument('num_global_restart', type=int, help='Number of global restarts. This forgets everything between each restart.')
parser.add_argument('initial_population_size', type=int, help='Number of initial population for each global restart. We randomly sample this many populations and select the best one to work with.')

parser.add_argument('seed', type=int, help='Seed for random number generator')
parser.add_argument('--timeout', type=float, help='Timeout for the entire process')
args = parser.parse_args()

random.seed(args.seed)

# make a temporary directory

print(f'Running a new restart with seed {args.seed}, n = {args.n}, m = {args.m}, d = {args.d}, p = {args.p}')
temp_dir = tempfile.TemporaryDirectory()
temp_dir = temp_dir.name
mkdir(temp_dir)

if args.timeout:
    stop_event = threading.Event()
    timer_thread = threading.Thread(target=send_sigint_to_process_group_after_delay, args=(args.timeout, stop_event))
    timer_thread.start()

try:
    print(f"Generating stein kernel with {args.n} points")
    subprocess.run(['../src/code_snippet_FC.py', str(args.n), str(args.d), str(args.m), join(temp_dir, 'points.p'), join(temp_dir, 'raw_points.p'), '--seed', str(args.seed)], check=True)
    for i in range(args.num_global_restart):
        try:
            p = subprocess.Popen([sys.executable, '-u', '../src/perturb.py', join(temp_dir, 'points.p'), join(temp_dir, 'raw_points.p'), str(random.randrange(0, 2**63)), str(args.p), str(args.number_local_restart), str(args.initial_population_size)])
            p.wait()

        finally:
            p.wait()
            if p.returncode != 0:
                raise subprocess.CalledProcessError(p.returncode, p.args)
finally:
    if args.timeout:
        stop_event.set()
    print("cleaning up temporary directory")
    shutil.rmtree(temp_dir)