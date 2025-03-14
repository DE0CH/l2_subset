#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <signal.h>
#include <string.h>
#include <errno.h>
#include <sys/mman.h>
#include "utils.h"
#include "l2_subset.h"
#include <time.h>
#include "mt19937-64/mt64.h"

#define _XOPEN_SOURCE 700

volatile sig_atomic_t run = true;

void signal_handler(int sig) {
    if (!run) {
        signal(sig, SIG_DFL);
        raise(sig);
    } else {
        fprintf(stderr, "Caught SIGINT, setting run to false. Going to print out the best points so far.\nSend the SIGINT again (e.g. ctrl-C) to force quit the program.\n");
        run = false;
    }
}

double get_weight(struct weights *w, size_t i, size_t j) {
#if COMPUTE_MODE == USE_MATRIX
    return w->entries[i * w->n + j];
#elif COMPUTE_MODE == USE_POINTS
    return w_ij(w->points, i, j, w->d, w->m, w->n);
#else
    #error "Invalid COMPUTE_MODE"
#endif
}

void replace_points(struct weights *w, size_t dest, size_t src) {
#if DEBUG_SLOW
    double predicted_change = -w->point_weights[dest] + relative_inactive_weight(w, src, dest);
    double old_total_discrepancy = w->total_discrepancy;
#endif
    remove_point(w, dest);
    add_point(w, src);
#if DEBUG_SLOW
    if (!isclose(predicted_change, w->total_discrepancy - old_total_discrepancy)) {
        fprintf(stderr, "Predicted change: %.10lf, actual change: %.10lf\n", predicted_change, w->total_discrepancy - old_total_discrepancy);
        fprintf(stderr, "dest: %zu, src: %zu\n", dest, src);
        exit(EXIT_FAILURE);
    }
#endif
}

void print_array(double *arr, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        printf("%.10lf ", arr[i]);
    }
    printf("\n");
}

void print_matrix(struct weights *w) {
    for (size_t i = 0; i < w->n; ++i) {
        for (size_t j = 0; j < w->n; ++j) {
            printf("%.10lf ", get_weight(w, i, j));
        }
        printf("\n");
    }
}

bool isclose(double a, double b) {
    double rel_tol = 1e-8;
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
        printf("found: %.10lf expected: %.10lf\n", total_discrepancy, w->total_discrepancy);
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
    fprintf(stderr, "Adding point %zu\n", j);
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
    fprintf(stderr, "Removing point %zu\n", j);
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

void calculate_init_weights(struct weights *w) {
    w->total_discrepancy = 0.0;
    for (size_t i = 0; i < w->n; i++) {
	    w->point_weights[i] = get_weight(w, i, i);
    }
#if DEBUG_SLOW
    bool points_category[w->n];
    memcpy(points_category, w->points_category, w->n * sizeof(bool));
    memset(w->points_category, INACTIVE, w->n * sizeof(bool));
#else
    bool *points_category = w->points_category;
#endif
    for (size_t i = 0; i < w->n; i++) {
        if (points_category[i] == ACTIVE) {
            add_point(w, i);
        }
    }
}

double w_ij(double* X, size_t i, size_t j, long long d, long long m, long long n) {
    if (i == j) {
        double prod1 = 1.0;
        double prod2 = 1.0;
        for (size_t h = 0; h < d; h++) {
            prod1 *= (3 - pow(X[i * d + h], 2));
            prod2 *= (2 - X[i * d + h]);
        }
        return - prod1 * pow(2, 1-d) / m + prod2 / (m * m);
    } else {
        double prod = 1.0;
        for (size_t h = 0; h < d; h++) {
            prod *= min(2 - X[i * d + h], 2 - X[j * d + h]);
        }
        return prod / (m * m);
    }
}

struct weights *weights_alloc(size_t d, size_t n) {
    struct weights *w = (struct weights *)malloc(sizeof(struct weights));
#if COMPUTE_MODE == USE_MATRIX
    w->entries = (double *)malloc(n * n * sizeof(double));
#elif COMPUTE_MODE == USE_POINTS
#else
    #error "Invalid COMPUTE_MODE"
#endif
    w->n = n;
    w->points = (double *)malloc(n * d * sizeof(double));
    w->point_weights = (double *)malloc(n * sizeof(double));
    w->points_category = (bool *)malloc(n * sizeof(bool));
    w->pp = (size_t *)malloc(n * sizeof(size_t));
    return w;
}

void weights_free(struct weights *w) {
#if COMPUTE_MODE == USE_MATRIX
    free(w->entries);
#elif COMPUTE_MODE == USE_POINTS
#else
    #error "Invalid COMPUTE_MODE"
#endif
    free(w->points);
    free(w->point_weights);
    free(w->points_category);
    free(w->pp);
    free(w);
}

void select_random_points(struct weights *w) {
    size_t m = w->m;
    size_t n = w->n;
    size_t resevoir[m];

    resevoir_sample(resevoir, n, m);
    array_to_mask(w->points_category, resevoir, n, m);
    calculate_init_weights(w);
}

void process_points(struct weights *w) {
#if COMPUTE_MODE == USE_MATRIX
    for (size_t i = 0; i < w->n; ++i) {
        for (size_t j = 0; j < w->n; ++j) {
            w->entries[i * w->n + j] = w_ij(w->points, i, j, w->d, w->m, w->n);
        }
    }
#elif COMPUTE_MODE == USE_POINTS
#else
    #error "Invalid COMPUTE_MODE"
#endif
}

void resevoir_sample(size_t *resevoir, size_t n, size_t k) {
    for (size_t i = 0; i < k; i++) {
        resevoir[i] = i;
    }
    for (size_t i = k; i < n; i++) {
        size_t j = genrand64_int64() % (i + 1);
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

long long atoi_or_die(char *str) {
    char *endptr;
    long long result = strtoll(str, &endptr, 10);
    if (*str == '\0' || *endptr != '\0') {
        die("Invalid number: \"%s\"", str);
    }
    return result;
}

struct weights *read_point_file(struct input_data *data, int argc, char *argv[]) { // return the points

    if (argc != 5) {
        die("Usage: %s <points_file> <m> <seed> <n_trials>. Selecte m low discrepancy points.", argv[0]);
    }
    long long m = atoi_or_die(argv[2]);
    long long seed = atoi_or_die(argv[3]);
    long long n_trials = atoi_or_die(argv[4]);

    data->n_trials = n_trials;
    data->seed = seed;

    long long d, n;
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

struct weights *read_from_compiled_matrix(struct input_data *data, int argc, char *argv[], void **mmaped_data) {
    if (argc != 3) {
        die("Usage: %s <compiled_matrix_file> <seed>. Selecte m low discrepancy points.", argv[0]);
    }
    long long seed = atoi_or_die(argv[2]);
    data->seed = seed;
    struct weights *w = weights_deserialize(argv[1], mmaped_data);
    if (w == NULL) {
        die("Could not read compiled matrix file: %s", argv[1]);
    }
    return w;
}

struct weights *read_from_compiled_matrix_w_starting_point(struct input_data *data, int argc, char *argv[], void **mmaped_data) {\
    if (argc < 3) {
        die("Usage: %s <compiled_matrix_file> <seed> <starting point #1> <starting point #2> .... Selecte m low discrepancy points.", argv[0]);
    }
    long long seed = atoi_or_die(argv[2]);
    data->seed = seed;
    struct weights *w = weights_deserialize(argv[1], mmaped_data);
    if (w == NULL) {
        die("Could not read compiled matrix file: %s", argv[1]);
    }
    if (argc != 3 + w->m) {
        die("Usage: %s <compiled_matrix_file> <seed> <starting point #1> <starting point #2> .... Selecte m low discrepancy points.", argv[0]);
    }
    memset(w->points_category, INACTIVE, w->n * sizeof(bool));
    for (size_t i = 0; i < w->m; i++) {
        if (3 + i >= argc) {
            die("%s not enough starting point is present in the argument", argv[0]);
        }
        size_t pi = atoi_or_die(argv[3 + i]);
        if (pi >= w->n) {
            die("Point index out of bounds\n");
        }
        w->points_category[pi] = ACTIVE;
    }
    calculate_init_weights(w);
    return w;
}

struct weights *read_point_file_and_save(struct input_data *data, int argc, char *argv[]) {
    if (argc != 4) {
        die("Usage: %s <points_file> <m> <output_file>. Precompute the matrix and save to output file.", argv[0]);
    }
    long long m = atoi_or_die(argv[2]);
    long long d, n;
    double *points = read_points_from_file(argv[1], &d, &n);
    data->output_filename = argv[3];
    struct weights *w = weights_alloc(d, n);
    free(w->points);
    w->points = points;
    w->m = m;
    w->n = n;
    w->d = d;
    process_points(w);
    return w;
}

size_t round_up(size_t location, size_t alignment) {
    return (location + alignment - 1) & ~(alignment - 1);
}

int weights_serialize(struct weights *w, char *filename) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        return 1;
    }
    size_t written = 0;
    size_t total = 0;
    struct serialize_header h;
    h.n = w->n;
    h.m = w->m;
    h.d = w->d;

    written += fwrite(&h, sizeof(struct serialize_header), 1, file);
    total += 1;
    // pad zero until alignment of double
    size_t diff = round_up(sizeof(struct serialize_header), sizeof(double)) - sizeof(struct serialize_header);
    char padding = '\0';
    written += fwrite(&padding, sizeof(char), diff, file);
    total += diff;
    written += fwrite(w->points, sizeof(double), w->n * w->d, file);
    total += w->n * w->d;
#if COMPUTE_MODE == USE_MATRIX
    written += fwrite(w->entries, sizeof(double), w->n * w->n, file);
    total += w->n * w->n;
#elif COMPUTE_MODE == USE_POINTS
#else
    #error "Invalid COMPUTE_MODE"
#endif
    fclose(file);
    if (written != total) {
        return 1;
    }
    return 0;
}

size_t weight_serialized_file_size(struct serialize_header h) {
    size_t padding_length = round_up(sizeof(struct serialize_header), sizeof(double)) - sizeof(struct serialize_header);
    size_t filesize = sizeof(struct serialize_header) + padding_length + h.n * h.d * sizeof(double);
#if COMPUTE_MODE == USE_MATRIX
    filesize += h.n * h.n * sizeof(double);
#elif COMPUTE_MODE == USE_POINTS
#else
    #error "Invalid COMPUTE_MODE"
#endif
    return filesize;
}

struct weights *weights_deserialize(char *filename, void **mmapedData) {

    // mmap the file to entries
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        return NULL;
    }
    struct serialize_header h;
    size_t read = fread(&h, sizeof(struct serialize_header), 1, file);
    if (read != 1) {
        fclose(file);
        return NULL;
    }
    long long n = h.n;
    long long m = h.m;
    long long d = h.d;
    size_t padding_length = round_up(sizeof(struct serialize_header), sizeof(double)) - sizeof(struct serialize_header);
    fseek(file, padding_length, SEEK_CUR);
    struct weights *w = weights_alloc(d, n);
    free(w->points);
#if COMPUTE_MODE == USE_MATRIX
    free(w->entries);
#elif COMPUTE_MODE == USE_POINTS
#else
    #error "Invalid COMPUTE_MODE"
#endif
    w->m = m;
    w->d = d;
    w->n = n;
    size_t filesize = weight_serialized_file_size(h);

    *mmapedData = mmap(NULL, filesize, PROT_READ, MAP_PRIVATE, fileno(file), 0);
    if (*mmapedData == MAP_FAILED) {
        fclose(file);
        return NULL;
    }

    size_t offset = sizeof(struct serialize_header) + padding_length;
    double *points = (double *)((char *)*mmapedData + offset);
    w->points = points;
    offset += sizeof(double) * n * d;
#if COMPUTE_MODE == USE_MATRIX
    double *entries = (double *)((char *)*mmapedData + offset);
    w->entries = entries;
#elif COMPUTE_MODE == USE_POINTS
#else
    #error "Invalid COMPUTE_MODE"
#endif
    fclose(file);
    return w;
}

int free_mmaped_matrix(struct weights *w, void *mmaped_data) {
    struct serialize_header h;
    h.n = w->n;
    h.m = w->m;
    h.d = w->d;
    size_t filesize = weight_serialized_file_size(h);
    if (munmap(mmaped_data, filesize) == -1) {
        return 1;
    }
    return 0;
}

struct analytics *analytics_alloc() {
    return (struct analytics *)malloc(sizeof(struct analytics));
}

void analytics_free(struct analytics *a) {
    free(a);
}

void shuffle(struct weights *w){
    size_t *pp = w->pp;
    size_t n = w->n;
    for (size_t i = 0; i < n - 1; i++) {
        size_t j = i + genrand64_int64() % (n - i);
        size_t t = pp[j];
        pp[j] = pp[i];
        pp[i] = t;
    }
}

struct pair most_significant_pair(struct weights *w) {
    struct pair p = {SIZE_MAX, SIZE_MAX};
    double min = 0.0;
    size_t *pp = w->pp;
    for (size_t ii = 0; ii < w->n; ii++) {
        size_t i = pp[ii];
        if (w->points_category[i] == ACTIVE) {
            for (size_t jj = 0; jj < w->n; jj++) {
                size_t j = pp[jj];
                if (w->points_category[j] == INACTIVE) {
                    double weight = -w->point_weights[i] + relative_inactive_weight(w, j, i);
                    if (weight < min) {
                        min = weight;
                        p.i = i;
                        p.j = j;
                        return p;
                    }
                }
            }
        }
    }
    return p;

}

void initialise_permuation(struct weights *w) {
    for (size_t i = 0; i < w->n; i++) {
        w->pp[i] = i;
    }
}

struct analytics *main_loop(struct weights *w) {
    struct sigaction sa;
    sa.sa_handler = signal_handler;
    sa.sa_flags = 0;
    sigemptyset(&sa.sa_mask);

    if (sigaction(SIGINT, &sa, NULL) == -1) {
        perror("Error setting SIGINT handler");
        exit(EXIT_FAILURE);
    }

    struct analytics *a = analytics_alloc();
    initialise_permuation(w);
    shuffle(w);
    a->num_iterations = 0;
    if (w->n == w->m) {
        return a;
    }
    while (run) {
        printf("l2   %.10lf\n", total_discrepancy(w));
        struct pair p = most_significant_pair(w);
        if (p.i == SIZE_MAX) {
            break;
        }
        replace_points(w, p.i, p.j);
        a->num_iterations++;
    }
    if (!run) {
        printf("Interrupted\n");
        print_results(w, a);
        raise(SIGINT);
    }
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

    printf("Number of iterations: %lld\n", a->num_iterations);
    printf("Active point sum: %.10lf\n", total_discrepancy(w));
    fflush(stdout);
}

double total_discrepancy(struct weights *w) {
    return w->total_discrepancy;
}
