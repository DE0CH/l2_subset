#ifndef L2_SUBSET_H
#define L2_SUBSET_H

#include <stdbool.h>
#include "utils.h"
#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
#define ACTIVE 1
#define INACTIVE 0
#define USE_MATRIX 1
#define USE_POINTS 2
#ifndef COMPUTE_MODE
#define COMPUTE_MODE 1
#endif
#define die(fmt, ...) do { \
    fprintf(stderr, "Error: " fmt "\n", ##__VA_ARGS__); \
    exit(EXIT_FAILURE); \
} while (0)

// Structure definitions

struct weights {
    long long n;
    long long m;
    long long d;
    char mode;
#if COMPUTE_MODE == USE_MATRIX
    double *entries;
#elif COMPUTE_MODE == USE_POINTS
#else
    #error "Invalid COMPUTE_MODE"
#endif
    double *points;
    double *point_weights;
    double total_discrepancy;
    bool *points_category;
    size_t *pp; // permuation of lengh n, used for shuffling
};

struct analytics {
    long long num_iterations;
};

struct serialize_header {
    long long n;
    long long m;
    long long d;
};

struct pair {
    size_t i;
    size_t j;
};

// Function prototypes

// Weights functions
double get_weight(struct weights *w, size_t i, size_t j);
void replace_points(struct weights *w, size_t dest, size_t src);
void add_point(struct weights *w, size_t rc);
void remove_point(struct weights *w, size_t src);
void recalculate_weights(struct weights *w);
struct weights *weights_alloc(size_t d, size_t n);
void weights_free(struct weights *w);
void process_points(struct weights *w);
size_t largest_active_point(struct weights *w);
double relative_inactive_weight(struct weights *w, size_t inactive_point, size_t active_point);
size_t smallest_inactive_point(struct weights *w, size_t largest_active_point);
double w_ij(double* X, size_t i, size_t j, long long d, long long m, long long n);
void debug_cmp(struct weights *w);
void select_random_points(struct weights *w);

// Utility functions
void resevoir_sample(size_t *resevoir, size_t n, size_t k);
void array_to_mask(bool *mask, size_t *array, size_t n, size_t k);
long long atoi_or_die(char *str);
struct weights *read_point_file(struct input_data *data, int argc, char *argv[]);
double relative_inactive_weight(struct weights *w, size_t inactive_point, size_t active_point);
bool isclose(double a, double b);
void print_matrix(struct weights *w);
void print_array(double *arr, size_t n);
void initialise_permuation(struct weights *w);
void shuffle(struct weights *w);

// Analytics functions
double total_discrepancy(struct weights *w);
struct analytics *analytics_alloc(void);
void analytics_free(struct analytics *a);
struct analytics *main_loop(struct weights *w);
void print_results(struct weights *w, struct analytics *a);

struct weights *read_point_file_and_save(struct input_data *data, int argc, char *argv[]);
size_t weight_serialized_file_size(struct serialize_header h);
int weights_serialize(struct weights *w, char *filename);
struct weights *weights_deserialize(char *filename, void **mmapedData);
struct weights *read_from_compiled_matrix(struct input_data *data, int argc, char *argv[], void **mmaped_data);
struct weights *read_from_compiled_matrix_w_starting_point(struct input_data *data, int argc, char *argv[], void **mmaped_data);
int free_mmaped_matrix(struct weights *w, void *mmaped_data);

#endif // L2_SUBSET_H
