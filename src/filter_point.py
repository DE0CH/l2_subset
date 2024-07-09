import argparse

parser = argparse.ArgumentParser(description='Select points by index')
parser.add_argument('point_file', type=str, help='Point File')
parser.add_argument('selected_points', type=str, help='Selected Points File')
args = parser.parse_args()

with open(args.selected_points) as f:
    selected_points = list(map(int, f.readline().split()))

with open(args.point_file) as f:
    d, n, _ = f.readline().split()
    points = []
    for line in f:
        points.append(line)
    remaining_points = []
    for i in selected_points:
        remaining_points.append(points[i])
print(d, len(selected_points), "0.0", sep=' ')
print(''.join(remaining_points), end='')
