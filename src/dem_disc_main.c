#include <stdlib.h>
#include <stdio.h>
#include "dem_disc.h"
#include "l2_subset.h"

int main(int argc, char **argv)
{
    if (argc != 2) {
        printf("Usage: %s [file]\n", argv[0]);
        return 1;
    }
    long long n, d;
    double *points_store = read_points_from_file(argv[1], &d, &n);
    double **points = malloc(n * sizeof(double *));
    for (size_t i = 0; i < n; i++) {
        points[i] = points_store + i * d;
    }
    double discrepancy = oydiscr(points, d, n);
    printf("%.10lf\n", discrepancy);
    free(points);
    free(points_store);
}
