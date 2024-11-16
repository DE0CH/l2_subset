#include <stdio.h>
#include <stdlib.h>
#include "utils.h"

double *read_points_from_file(char *filename, long long *d, long long *n) { // return the points, write d and n to the pointers supplied in the arugment
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        die("Could not open file: %s", filename);
    }
    if (fscanf(file, "%lld %lld %*f", d, n) != 2) {
        die("Could not read dimensions. Hint: check if the file has the right format\n");
    }

    if (*d < 2) {
        die("The dimension must be at least 2\n");
    }

    double *points = (double *)malloc(*n * *d * sizeof(double));

    for (size_t i = 0; i < *n**d; i++) {
        if(fscanf(file, "%lf", &points[i]) != 1) {
            die("Could not read point. Hint: check if the file has the right number of points, and if the points are actually numbers\n");
        }
    }
    double temp;
    if (fscanf(file, "%lf", &temp) == 1) {
        die("Too many points. Hint: check if you have given the right number of points at the top of the file\n");
    }
    fclose(file);
    return points;
}