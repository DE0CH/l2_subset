import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Select points by index')
parser.add_argument('point_file', type=str, help='Point File')
parser.add_argument('selected_points', type=str, help='Selected Points File') 
args = parser.parse_args()

with open(args.selected_points) as f:
    selected_points = np.array(list(map(int, f.readline().split())))

with open(args.point_file) as f:
    f.readline()
    points = []
    for line in f:
        points.append(line)
    points = np.array(points)
    print(''.join(points[selected_points]), end='')

