#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <errno.h>
#include <sys/mman.h>
#include "l2_subset.h"
#include "dem_disc.h"

double get_weight(struct weights *w, size_t i, size_t j) {
    return w_ij(w->points, i, j, w->d, w->m, w->n);
}

void replace_points(struct weights *w, size_t dest, size_t src) {
    remove_point(w, dest);
    add_point(w, src);
}

void add_point(struct weights *w, size_t src) {
    size_t j = src;
    w->total_discrepancy += w->point_weights[j];
    double old_weight = w->point_weights[j];
    for (size_t i = 0; i < w->n; ++i) {
        double change = 2*get_weight(w, i, j);
        w->point_weights[i] += change;
    }
    w->points_category[j] = ACTIVE;
    w->point_weights[j] = old_weight;
}


void remove_point(struct weights *w, size_t src) {
    size_t j = src;
    w->total_discrepancy -= w->point_weights[j];
    double old_weight = w->point_weights[j];
    for (size_t i = 0; i < w->n; i++) {
        double change = 2*get_weight(w, i, j);
        w->point_weights[i] -= change;
    }
    w->point_weights[j] = old_weight;
    w->points_category[j] = INACTIVE;
}

void recalculate_weights(struct weights *w) {
    w->total_discrepancy = 0.0;
    for (size_t i = 0; i < w->n; i++) {
        w->point_weights[i] = 0.0;
        // used to be some meaningful calculation, but no needed to reproduce the bug
    }
}

double w_ij(double* X, int i, int j, int d, int m, int n) {
    // used to be some meaningful calculation, but no needed to reproduce the bug
    return 1;
}

struct weights *weights_alloc(size_t d, size_t n) {
    struct weights *w = (struct weights *)malloc(sizeof(struct weights));
    w->n = n;
    w->points = (double *)malloc(n * d * sizeof(double));
    w->point_weights = (double *)malloc(n * sizeof(double));
    w->points_category = (bool *)malloc(n * sizeof(bool));
    return w;
}

void weights_free(struct weights *w) {
    free(w->points);
    free(w->point_weights);
    free(w->points_category);
    free(w);
}

void select_random_points(struct weights *w) {
    size_t m = w->m;
    size_t n = w->n;
    size_t resevoir[m];

    resevoir_sample(resevoir, n, m);
    array_to_mask(w->points_category, resevoir, n, m);
    recalculate_weights(w);
}

void process_points(struct weights *w) {
}

void resevoir_sample(size_t *resevoir, size_t n, size_t k) {
    for (size_t i = 0; i < k; i++) {
        resevoir[i] = i;
    }
    for (size_t i = k; i < n; i++) {
        size_t j = rand() % (i + 1);
        if (j < k) {
            resevoir[j] = i;
        }
    }
}

// length of the array is k
// length of mask is n
void array_to_mask(bool *mask, size_t *array, size_t n, size_t k) {
    // actually I could use a for loop because compiler actually optimizes it to become memset, but I want to be cool...
    memset(mask, false, n * sizeof(bool));
    for (size_t i = 0; i < k; ++i) {
        mask[array[i]] = true;
    }
}

int atoi_or_die(char *str) {
    char *endptr;
    int result = strtol(str, &endptr, 10);
    if (*str == '\0' || *endptr != '\0') {
        die("Invalid number: \"%s\"", str);
    }
    return result;
}

double *read_points_from_file(char *filename, int *d, int *n) { // return the points, write d and n to the pointers supplied in the arugment
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        die("Could not open file: %s", filename);
    }
    if (fscanf(file, "%d %d %*f", d, n) != 2) {
        die("Could not read dimensions. Hint: check if the file has the right format\n");
    }

    if (*d < 2) {
        die("The dimension must be at least 2\n");
    }

    double *points = (double *)malloc(*n * *d * sizeof(double));

    for (int i = 0; i < *n**d; i++) {
        if(fscanf(file, "%lf", &points[i]) != 1) {
            die("Could not read point. Hint: check if the file has the right number of points, and if the points are actually numbers\n");
        }
    }
    if (fscanf(file, "%lf", points) == 1) {
        die("Too many points. Hint: check if you have given the right number of points at the top of the file\n");
    }
    fclose(file);
    return points;
}

struct weights *read_point_file(struct input_data *data, int argc, char *argv[]) { // return the points

    if (argc != 5) {
        die("Usage: %s <points_file> <m> <seed> <n_trials>. Selecte m low discrepancy points.", argv[0]);
    }
    int m = atoi_or_die(argv[2]);
    int seed = atoi_or_die(argv[3]);
    int n_trials = atoi_or_die(argv[4]);

    data->n_trials = n_trials;
    data->seed = seed;

    int d, n;
    double *points = read_points_from_file(argv[1], &d, &n);
    struct weights *w = weights_alloc(d, n);
    free(w->points);
    w->points = points;
    w->m = m;
    w->n = n;
    w->d = d;
    process_points(w);
    process_points(w);

    return w;
}

size_t round_up(size_t location, size_t alignment) {
    return (location + alignment - 1) & ~(alignment - 1);
}

void main_loop(struct weights *w) {
    printf("linf %.10lf\n", linf_disc(w));
    replace_points(w, 1, 2);
}

void print_results(struct weights *w) {
    printf("Active points: ");
    for (size_t i = 0; i < w->n; ++i) {
        if (w->points_category[i] == ACTIVE) {
            printf("%zu ", i);
        }
    }
    printf("\n");

    printf("Active point sum: %.10lf\n", total_discrepancy(w));
    printf("linf discrepancy: %.10lf\n", linf_disc(w));
    fflush(stdout);
}

double total_discrepancy(struct weights *w) {
    double constant = 1.0 / (double)pow(3, w->d);
    return constant + w->total_discrepancy;
}

double linf_disc(struct weights *w) {
    int n = w->n;
    int d = w->d;
    int m = w->m;
    double **points = malloc(m * sizeof(double *));
    int j = 0;
    for (int i = 0; i < n; i++) {
        if (w->points_category[i] == ACTIVE) {
            points[j] = w->points + i * d;
            j++;
        }
    }
    double ans = oydiscr(points, d, m);
    free(points);
    return ans;
}

int main(int argc, char *argv[]) {
    struct input_data data;
    struct weights *w = read_point_file(&data, argc, argv);
    srand(data.seed);
    select_random_points(w);
    main_loop(w);
    print_results(w);

    weights_free(w);
    return 0;
}
