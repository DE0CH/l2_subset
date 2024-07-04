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
#define die(fmt, ...) do { \
    fprintf(stderr, "Error: " fmt "\n", ##__VA_ARGS__); \
    exit(EXIT_FAILURE); \
} while (0)

typedef long long ll;
typedef unsigned char bool;

// Structure definitions
struct node {
    size_t value;
    struct node *next;
    struct node *prev;
};

struct small_allocator {
    struct node *mem;
    struct node *mem_len;
    struct node **stack;
    struct node **stack_len;
};

struct linked_list {
    struct small_allocator *allocator;
    struct node head;
    struct node end;
};

struct weights {
    size_t n;
    double *entries;
    double *point_weights;
    double active_point_sum;
    struct linked_list *active_points;
    struct linked_list *inactive_points;
};

struct analytics {
    ll num_iterations;
};

// Function prototypes

// Linked List functions
struct linked_list *linked_list_alloc(size_t n);
void linked_list_free(struct linked_list *list);
void linked_list_add(struct linked_list *list, size_t point);
void linked_list_remove(struct linked_list *list, struct node *node);
struct node *linked_list_begin(struct linked_list *list);
struct node *linked_list_end(struct linked_list *list);

// Small Allocator functions
struct small_allocator *small_allocator_alloc(size_t n);
void small_allocator_free(struct small_allocator *alloc);
struct node *node_alloc(struct small_allocator *allocator);
void node_free(struct small_allocator *allocator, struct node *node);

// Weights functions
double get_weight(struct weights *w, size_t i, size_t j);
void replace_points(struct weights *w, struct node *dest, struct node *src);
void add_point(struct weights *w, struct node *src);
void remove_point(struct weights *w, struct node *src);
void recalculate_weights(struct weights *w);
struct weights *weights_alloc(size_t n);
void weights_free(struct weights *w);
void process_points(struct weights *w, double *points, int d, int m, int n);
struct node *largest_active_point(struct weights *w);
double relative_inactive_weight(struct weights *w, size_t inactive_point, size_t active_point);
struct node *smallest_inactive_point(struct weights *w, struct node *largest_active_point);
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