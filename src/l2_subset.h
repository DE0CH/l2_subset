#ifndef L2_SUBSET_H
#define L2_SUBSET_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
#define true 1
#define false 0
#define ACTIVE 1
#define INACTIVE 0
#define die(fmt, ...) do { \
    fprintf(stderr, "Error: " fmt "\n", ##__VA_ARGS__); \
    exit(EXIT_FAILURE); \
} while (0)

typedef long long ll;
typedef unsigned char bool;

// Structure definitions

struct weights {
    size_t n;
    size_t d;
    double *entries;
    double *point_weights;
    double active_point_sum;
    bool *points_category;
};

struct analytics {
    ll num_iterations;
};

// Function prototypes

// Weights functions
double get_weight(struct weights *w, size_t i, size_t j);
void replace_points(struct weights *w, size_t dest, size_t src);
void add_point(struct weights *w, size_t rc);
void remove_point(struct weights *w, size_t src);
void recalculate_weights(struct weights *w);
struct weights *weights_alloc(size_t n);
void weights_free(struct weights *w);
void process_points(struct weights *w, double *points, int d, int m, int n);
size_t largest_active_point(struct weights *w);
double relative_inactive_weight(struct weights *w, size_t inactive_point, size_t active_point);
size_t smallest_inactive_point(struct weights *w, size_t largest_active_point);
double w_ij(double* X, int i, int j, int d, int m, int n);
void debug_cmp(struct weights *w);

// Utility functions
void resevoir_sample(size_t *resevoir, size_t n, size_t k);
void array_to_mask(bool *mask, size_t *array, size_t n, size_t k);
int atoi_or_die(char *str);
struct weights *read_input(int argc, char *argv[]);

// Analytics functions
struct analytics *analytics_alloc(void);
void analytics_free(struct analytics *a);
struct analytics *main_loop(struct weights *w);
void print_results(struct weights *w, struct analytics *a);

#endif // L2_SUBSET_H