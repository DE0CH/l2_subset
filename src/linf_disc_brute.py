import argparse
import itertools

parser = argparse.ArgumentParser(description='Compute the L_inf using brute force')
parser.add_argument('file', type=str, help='The file containing the data')
args = parser.parse_args()

points = []
with open(args.file, 'r') as f:
    d, n, _ = f.readline().split()
    d = int(d)
    n = int(n)
    for line in f:
        points.append(list(map(float, line.split())))

index_one = range(n+1)

def get_volume(point_index):
    volume = 1.0
    for id in range(d):
        if point_index[id] == n:
            volume *= 1.0
        else:
            volume *= points[point_index[id]][id]
    return volume

def open(point_index, i): # is point i in the open box bounded by the corder point_index? 
    for id in range(d):
        if point_index[id] == n:
            corner_value = 1.0
        else:
            corner_value = points[point_index[id]][id]
        if points[i][id] >= corner_value:
            return False
    return True

def closed(point_index, i):
    for id in range(d):
        if point_index[id] == n:
            corner_value = 1.0
        else:
            corner_value = points[point_index[id]][id]
        if points[i][id] > corner_value:
            return False
    return True

def delta(point_index):
    m = 0
    for i in range(n):
        if open(point_index, i):
            m += 1
    return get_volume(point_index) - m / n

def bardelta(point_index):
    m = 0
    for i in range(n):
        if closed(point_index, i):
            m += 1
    return m / n - get_volume(point_index)

v = 0.0

for point in itertools.product(*[index_one]*d):
    point = list(point)
    for i in point:
        v = max(v, delta(point), bardelta(point))

print(v)
