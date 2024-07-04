#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "l2_subset.h"
#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
#define true 1
#define false 0
#define die(fmt, ...) do { \
    fprintf(stderr, "Error: " fmt "\n", ##__VA_ARGS__); \
    exit(EXIT_FAILURE); \
} while (0)
typedef long long ll;
typedef unsigned char bool;

double get_weight(struct weights *w, size_t i, size_t j) {
    return w->entries[i * w->n + j];
}

void replace_points(struct weights *w, size_t dest, size_t src) {
    remove_point(w, dest);
    add_point(w, src);
}

void print_array(double *arr, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        printf("%lf ", arr[i]);
    }
    printf("\n");
}

void print_matrix(struct weights *w) {
    for (size_t i = 0; i < w->n; ++i) {
        for (size_t j = 0; j < w->n; ++j) {
            printf("%lf ", get_weight(w, i, j));
        }
        printf("\n");
    }
}

bool isclose(double a, double b) {
    double rel_tol = 1e-9;
    double abs_tol = 0.0;
    return (fabs(a-b) <= max(rel_tol * max(fabs(a), fabs(b)), abs_tol));
}

bool double_array_close(double *a, double *b, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        if (!isclose(a[i], b[i])) {
            return false;
        }
    }
    return true;
}

void print_bool_arr(bool *arr, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        printf("%d ", arr[i]);
    }
    printf("\n");
}

void debug_cmp(struct weights *w) {
    size_t n = w->n;
    double point_weights[n];
    memcpy(point_weights, w->point_weights, n * sizeof(double));
    double total_discrepancy = w->total_discrepancy;
    recalculate_weights(w);
    if (!double_array_close(point_weights, w->point_weights, n)) {
        fprintf(stderr, "Active points mismatch\n");
        print_matrix(w);
        print_bool_arr(w->points_category, n);
        printf("found: ");
        print_array(point_weights, n);
        printf("expected: ");
        print_array(w->point_weights, n);
        exit(EXIT_FAILURE);
    }
    if (!isclose(total_discrepancy, w->total_discrepancy)) {
        print_matrix(w);
        printf("point_weights: ");
        print_array(point_weights, n);
        print_array(w->point_weights, n);
        printf("found: %lf expected: %lf\n", total_discrepancy, w->total_discrepancy);
        fprintf(stderr, "total_discrepancy mismatch\n");
        exit(EXIT_FAILURE);
    }
    fprintf(stderr, "All good\n");
}

void add_point(struct weights *w, size_t src) {
    size_t j = src;
    w->total_discrepancy += w->point_weights[j];
    double old_weight = w->point_weights[j];
    for (size_t i = 0; i < w->n; ++i) {
        double change = get_weight(w, i, j) + get_weight(w, j, i);
        w->point_weights[i] += change;
    }
    w->points_category[j] = ACTIVE;
    w->point_weights[j] = old_weight;
#if DEBUG_SLOW
    printf("Adding point %zu\n", j);
    debug_cmp(w);
#endif
}


void remove_point(struct weights *w, size_t src) {
    size_t j = src;
    w->total_discrepancy -= w->point_weights[j];
    double old_weight = w->point_weights[j];
    for (size_t i = 0; i < w->n; i++) {
        double change = get_weight(w, i, j) + get_weight(w, j, i);
        w->point_weights[i] -= change;
    }
    w->point_weights[j] = old_weight;
    w->points_category[j] = INACTIVE;
#if DEBUG_SLOW
    printf("Removing point %zu\n", j);
    debug_cmp(w);
#endif
}

void recalculate_weights(struct weights *w) {
    w->total_discrepancy = 0.0;
    for (size_t i = 0; i < w->n; i++) {
        w->point_weights[i] = 0.0;
        for (size_t j = 0; j < w->n; j++) {
            if (w->points_category[j] == ACTIVE || i == j) {
                double change = get_weight(w, i, j) + get_weight(w, j, i);
                if (i == j) {
                    change /= 2;
                }
                w->point_weights[i] += change;
                if (w->points_category[i] == ACTIVE) {
                    w->total_discrepancy += get_weight(w, i, j);
                }
            }
        }
    }
}

double w_ij(double* X, int i, int j, int d, int m, int n) {
    if (i == j) {
        double prod1 = 1.0;
        double prod2 = 1.0;
        for (int h = 0; h < d; h++) {
            prod1 *= (1 - pow(X[i * d + h], 2));
            prod2 *= (1 - X[i * d + h]);
        }
        return - prod1 * pow(2, 1-d) / m + prod2 / (m * m);
    } else {
        double prod = 1.0;
        for (int h = 0; h < d; h++) {
            prod *= min(1 - X[i * d + h], 1 - X[j * d + h]);
        }
        return prod / (m * m);
    }
}

struct weights *weights_alloc(size_t n) {
    struct weights *w = (struct weights *)malloc(sizeof(struct weights));
    w->entries = (double *)malloc(n * n * sizeof(double));
    w->n = n;
    w->point_weights = (double *)malloc(n * sizeof(double));
    w->points_category = (bool *)malloc(n * sizeof(bool));
    return w;
}

void weights_free(struct weights *w) {
    free(w->entries);
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

void process_points(struct weights *w, double *points, int d, int m, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            w->entries[i * n + j] = w_ij(points, i, j, d, m, n);
        }
    }
    w->d = d;
    w->m = m;
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

size_t largest_active_point(struct weights *w) {
    size_t max = SIZE_MAX;
    for (size_t i = 0; i < w->n; ++i) {
        if (w->points_category[i] == ACTIVE && (max == SIZE_MAX || w->point_weights[i] >= w->point_weights[max])) {
            max = i;
        }
    }
    return max;
}

double relative_inactive_weight(struct weights *w, size_t inactive_point, size_t active_point) {
    return w->point_weights[inactive_point] - get_weight(w, inactive_point, active_point) - get_weight(w, active_point, inactive_point);
}

size_t smallest_inactive_point(struct weights *w, size_t largest_active_point) {
    size_t min = SIZE_MAX; 
    for (size_t i = 0; i < w->n; i++) {
        if (w->points_category[i] == INACTIVE && (min == SIZE_MAX || relative_inactive_weight(w, i, largest_active_point) <= relative_inactive_weight(w, min, largest_active_point))) {
            min = i;
        }
    }
    return min;
}

int atoi_or_die(char *str) {
    char *endptr;
    int result = strtol(str, &endptr, 10);
    if (*str == '\0' || *endptr != '\0') {
        die("Invalid number: \"%s\"", str);
    }
    return result;
}

struct weights *read_input(struct input_data *data, int argc, char *argv[]) {
    double *points;
    
    if (argc != 5) {
        die("Usage: %s <points_file> <m> <seed> <n_trials>. Selecte m low discrepancy points.", argv[0]);
    }

    int m = atoi_or_die(argv[2]);
    int seed = atoi_or_die(argv[3]);
    int n_trials = atoi_or_die(argv[4]);
    
    FILE *file = fopen(argv[1], "r");
    if (file == NULL) {
        die("Could not open file: %s", argv[3]);
    }
    
    int n, d;
    if (fscanf(file, "%d %d %*f", &d, &n) != 2) {
        die("Could not read dimensions. Hint: check if the file has the right format\n");
    }
    points = (double *)malloc(n * d * sizeof(double));

    for (int i = 0; i < n*d; i++) {
        if(fscanf(file, "%lf", &points[i]) != 1) {
            die("Could not read point. Hint: check if the file has the right number of points, and if the points are actually numbers\n");
        }
    }

    fclose(file);
    data->n_trials = n_trials;
    data->seed = seed;

    struct weights *w = weights_alloc(n);
    process_points(w, points, d, m, n);
    free(points);
    return w;
}

struct analytics *analytics_alloc() {
    return (struct analytics *)malloc(sizeof(struct analytics));
}

void analytics_free(struct analytics *a) {
    free(a);
}

struct analytics *main_loop(struct weights *w) {
    struct analytics *a = analytics_alloc();
    a->num_iterations = 0;
    double old_sum;
    double new_sum = w->total_discrepancy;
    if (w->n == w->m) {
        return a;
    }
    do {
        size_t i = largest_active_point(w);
        size_t j = smallest_inactive_point(w, i);
        replace_points(w, i, j);
        old_sum = new_sum;
        new_sum = w->total_discrepancy;
        if (old_sum <= new_sum) {
            replace_points(w, j, i);
        }
        a->num_iterations++;
    } while (new_sum < old_sum);
    return a;
}

void print_results(struct weights *w, struct analytics *a) {
    printf("Active points: ");
    for (size_t i = 0; i < w->n; ++i) {
        if (w->points_category[i] == ACTIVE) {
            printf("%zu ", i);
        }
    }
    printf("\n");
    double constant = 1.0 / (double)pow(3, w->d);
    printf("Number of iterations: %lld\n", a->num_iterations);
    printf("Active point sum: %lf\n", constant + w->total_discrepancy);
}

int main(int argc, char *argv[]) {
    struct input_data data;
    struct weights *w = read_input(&data, argc, argv);
    srand(data.seed);
    for (ll i = 0; i < data.n_trials; i++) {
        printf("Trial %lld\n", i);
        select_random_points(w);
        struct analytics *a = main_loop(w);
        print_results(w, a);
        analytics_free(a);
    }

    weights_free(w);
    return 0;
}
