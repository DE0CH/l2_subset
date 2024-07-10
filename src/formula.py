import numpy as np
import argparse

def d_star_2(X):
    n, d = X.shape
    
    # First term
    term1 = 1 / (3 ** d)
    
    # Second term
    term2_sum = 0
    for i in range(n):
        prod = 1.0
        for j in range(d):
            prod *= (1 + 2 * X[i, j] - 2 * X[i, j] ** 2) / 4
        term2_sum += prod
    term2 = -2 / n * term2_sum
    
    # Third term
    term3_sum = 0
    for i in range(n):
        for j in range(n):
            prod = 1.0
            for k in range(d):
                prod *= (1 - abs(X[i, k] - X[j, k])) / 2
            term3_sum += prod
    term3 = 1 / (n ** 2) * term3_sum
    
    # Combine terms
    if (i == j):
        d_star = term1 + term2 + term3
    else:
        d_star = term1 + term3
    
    return d_star

parser = argparse.ArgumentParser(description='Compute L2')
parser.add_argument('point_file', type=str, help='Point File')
parser.add_argument('selected_points', type=str, help='Selected Points File')
args = parser.parse_args()

with open(args.selected_points) as f:
    selected_points = np.array(list(map(int, f.readline().split())))

with open(args.point_file) as f:
    f.readline()
    points = []
    for line in f:
        point = np.array(list(map(float, line.split())))
        points.append(point)
    points = np.array(points)
    print(d_star_2(points[selected_points, :]))
