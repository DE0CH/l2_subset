#ifndef _UTILS_H_
#define _UTILS_H_


#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
#define die(fmt, ...) do { \
    fprintf(stderr, "Error: " fmt "\n", ##__VA_ARGS__); \
    exit(EXIT_FAILURE); \
} while (0)

struct input_data {
    long long n_trials;
    long long seed;
    char *output_filename;
};

double *read_points_from_file(char *filename, long long *d, long long *n);

#endif