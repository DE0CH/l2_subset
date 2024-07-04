#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
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
#define DEBUG_SLOW 1

struct linked_list *linked_list_alloc(size_t n) {
    struct linked_list *list = (struct linked_list *)malloc(sizeof(struct linked_list));
    list->allocator = small_allocator_alloc(n);
    list->head.next = &list->end;
    list->end.prev = &list->head;
    return list;
}

void linked_list_free(struct linked_list *list) {
    small_allocator_free(list->allocator);
    free(list);
}

void linked_list_add(struct linked_list *list, size_t point) {
    struct node *new_node = node_alloc(list->allocator);
    new_node->value = point;
    new_node->next = list->head.next;
    new_node->prev = &list->head;
    new_node->prev->next = new_node;
    new_node->next->prev = new_node;
}

void linked_list_remove(struct linked_list *list, struct node *node) {
    node->prev->next = node->next;
    node->next->prev = node->prev;
    node_free(list->allocator, node);
}

struct node *linked_list_begin(struct linked_list *list) {
    return list->head.next;
}

struct node *linked_list_end(struct linked_list *list) {
    return &list->end;
}

struct small_allocator *small_allocator_alloc(size_t n) {
    struct small_allocator *alloc = (struct small_allocator *)malloc(sizeof(struct small_allocator));
    alloc->mem = (struct node *)malloc(n * sizeof(struct node));
    alloc->mem_len = alloc->mem;
    alloc->stack = (struct node **)malloc(n * sizeof(struct node *));
    alloc->stack_len = alloc->stack;
    return alloc;
}

void small_allocator_free(struct small_allocator *alloc) {
    free(alloc->mem);
    free(alloc->stack);
    free(alloc);
}

struct node *node_alloc(struct small_allocator *allocator) {
    if (allocator->stack_len == allocator->stack) {
        return allocator->mem_len++;
    } else {
        return *--(allocator->stack_len);
    }
}

void node_free(struct small_allocator *allocator, struct node *node) {
    *allocator->stack_len++ = node;
}

double get_weight(struct weights *w, size_t i, size_t j) {
    return w->entries[i * w->n + j];
}

void replace_points(struct weights *w, struct node *dest, struct node *src) {
    remove_point(w, dest);
    add_point(w, src);
}

void debug_cmp(struct weights *w) {
    size_t n = w->n;
    double active_points[n];
    memcpy(active_points, w->active_point_weights, n * sizeof(double));
    double inactive_points[n];
    memcpy(inactive_points, w->inactive_point_weights, n * sizeof(double));
    double active_point_sum = w->active_point_sum;
    double inactive_point_sum = w->inactive_point_sum;
    recalculate_weights(w);
    if (memcmp(active_points, w->active_point_weights, n * sizeof(double)) != 0) {
        fprintf(stderr, "Active points mismatch\n");
        exit(EXIT_FAILURE);
    }
    if (memcmp(inactive_points, w->inactive_point_weights, n * sizeof(double)) != 0) {
        fprintf(stderr, "Inactive points mismatch\n");
        exit(EXIT_FAILURE);
    }
    if (active_point_sum != w->active_point_sum) {
        fprintf(stderr, "Active point sum mismatch\n");
        exit(EXIT_FAILURE);
    }
    if (inactive_point_sum != w->inactive_point_sum) {
        fprintf(stderr, "Inactive point sum mismatch\n");
        exit(EXIT_FAILURE);
    }
}

void add_point(struct weights *w, struct node *src) {
    size_t j = src->value;
    for (struct node *node = linked_list_begin(w->active_points); node != linked_list_end(w->active_points); node = node->next) {
        size_t i = node->value;
        double change = get_weight(w, i, j) + get_weight(w, j, i);
        w->active_point_weights[i] += change;
        w->active_point_sum += change;
    }
    linked_list_add(w->active_points, src->value);
    linked_list_remove(w->inactive_points, src);
    // nothing to do for src because invariant guarantees that it has value of all the active points including itself
    for (struct node *node = linked_list_begin(w->inactive_points); node != linked_list_end(w->inactive_points); node = node->next) {
        size_t i = node->value;
        double change = get_weight(w, i, j) + get_weight(w, j, i);
        w->inactive_point_weights[i] += change;
        w->inactive_point_sum += change;
    }
#if DEBUG_SLOW
    debug_cmp(w);
#endif
}



void remove_point(struct weights *w, struct node *src) {
    size_t j = src->value;
    // need to remove itself from active list first to make sure its value doesn't get changed
    linked_list_remove(w->active_points, src);
    for (struct node *node = linked_list_begin(w->active_points); node != linked_list_end(w->active_points); node = node->next) {
        size_t i = node->value;
        double change = get_weight(w, i, j) + get_weight(w, j, i);
        w->active_point_weights[i] -= change;
        w->active_point_sum -= change;
    }
    for (struct node *node = linked_list_begin(w->inactive_points); node != linked_list_end(w->inactive_points); node = node->next) {
        size_t i = node->value;
        double change = get_weight(w, i, j) + get_weight(w, j, i);
        w->inactive_point_weights[i] -= change;
        w->inactive_point_sum -= change;
    }
    linked_list_add(w->inactive_points, src->value);
#if DEBUG_SLOW
    debug_cmp(w);
#endif
}

void recalculate_weights(struct weights *w) {
    w->active_point_sum = 0.0;
    for (struct node *node = linked_list_begin(w->active_points); node != linked_list_end(w->active_points); node = node->next) {
        size_t i = node->value;
        w->active_point_weights[i] = 0.0;
        for (struct node *node2 = linked_list_begin(w->active_points); node2 != linked_list_end(w->active_points); node2 = node2->next) {
            size_t j = node->value;
            double change = get_weight(w, i, j) + get_weight(w, j, i);
            w->active_point_weights[i] += change;
            w->active_point_sum += change;
        }
    }
    w->inactive_point_sum = 0.0;
    for (struct node *node = linked_list_begin(w->inactive_points); node != linked_list_end(w->inactive_points); node = node->next) {
        size_t i = node->value;
        w->inactive_point_weights[i] = 0.0;
        for (struct node *node2 = linked_list_begin(w->active_points); node2 != linked_list_end(w->active_points); node2 = node2->next) {
            size_t j = node->value;
            double change = get_weight(w, i, j) + get_weight(w, j, i); 
            w->inactive_point_weights[i] += change;
            w->inactive_point_sum += change;
        }
        double change = get_weight(w, i, i) + get_weight(w, i, i);
        w->inactive_point_weights[i] += change;
        w->inactive_point_sum += change;
    }
}

double w_ij(double* X, int i, int j, int d, int m, int n) {
    if (i == j) {
        double prod1 = 1.0;
        double prod2 = 1.0;
        for (int h = 0; h < d; h++) {
            prod1 *= (1 - pow(X[i * n + h], 2));
            prod2 *= (1 - X[i * n + h]);
        }
        return -pow(2, 1-d) / (2.0 * m) * prod1 + 1.0 / (2.0 * m * m) * prod2;
    } else {
        double prod = 1.0;
        for (int h = 0; h < d; h++) {
            prod *= min(1 - X[i * n + h], 1 - X[j * n + h]);
        }
        return prod / (m * m);
    }
}

struct weights *weights_alloc(size_t n) {
    struct weights *w = (struct weights *)malloc(sizeof(struct weights));
    w->entries = (double *)malloc(n * n * sizeof(double));
    w->n = n;
    w->active_point_weights = (double *)malloc(n * sizeof(double));
    w->inactive_point_weights = (double *)malloc(n * sizeof(double));
    w->active_points = linked_list_alloc(n);
    w->inactive_points = linked_list_alloc(n);
    return w;
}

void weights_free(struct weights *w) {
    free(w->entries);
    free(w->active_point_weights);
    free(w->inactive_point_weights);
    linked_list_free(w->active_points);
    linked_list_free(w->inactive_points);
    free(w);    
}

void process_points(struct weights *w, double *points, int d, int m, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            w->entries[i * n + j] = w_ij(points, i, j, d, m, n);
        }
    }
    
    size_t resevoir[m];
    bool mask[n];
    resevoir_sample(resevoir, n, m);
    array_to_mask(mask, resevoir, n, m);
    for (size_t i = 0; i < n; ++i) {
        if (mask[i]) {
            linked_list_add(w->active_points, i);
        } else {
            linked_list_add(w->inactive_points, i);
        }
    }
    recalculate_weights(w);
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
    for (size_t i = 0; i < n; ++i) {
        mask[array[i]] = true;
    }
}

struct node *largest_active_point(struct weights *w) {
    struct node *max = NULL;
    for (struct node *node = linked_list_begin(w->active_points); node != linked_list_end(w->active_points); node = node->next) {
        size_t i = node->value;
        if (max == NULL || w->active_point_weights[i] >= w->active_point_weights[node->value]) {
            max = node;
        }
    }
    return max;
}


double relative_inactive_weight(struct weights *w, size_t inactive_point, size_t active_point) {
    return w->inactive_point_weights[inactive_point] - get_weight(w, inactive_point, active_point) - get_weight(w, active_point, inactive_point);
}

struct node *smallest_inactive_point(struct weights *w, struct node *largest_active_point) {
    struct node *min = NULL;
    for (struct node *node = linked_list_begin(w->inactive_points); node != linked_list_end(w->inactive_points); node = node->next) {
        size_t i = node->value;
        if (min == NULL || relative_inactive_weight(w, i, largest_active_point->value) <= relative_inactive_weight(w, min->value, largest_active_point->value)) {
            min = node;
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

struct weights *read_input(int argc, char *argv[]) {
    double *points;
    int d, m, n;
    
    if (argc != 5) {
        die("Usage: %s <d> <n> <points_file> <m>", argv[0]);
    }

    d = atoi_or_die(argv[1]);
    n = atoi_or_die(argv[2]);
    m = atoi_or_die(argv[4]);

    points = (double *)malloc(n * d * sizeof(double));
    
    FILE *file = fopen(argv[3], "r");
    if (file == NULL) {
        die("Could not open file: %s", argv[3]);
    }

    for (int i = 0; i < n*d; i++) {
        if(fscanf(file, "%lf", &points[i]) != 1) {
            die("Could not read point. Hint: check if the file has the right number of points, and if the points are actually numbers\n");
        }
    }

    fclose(file);
    srand(42);
    
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
    double old_sum = w->active_point_sum;
    double new_sum;
    do {
        struct node *i = largest_active_point(w);
        struct node *j = smallest_inactive_point(w, i);
        replace_points(w, i, j);
        new_sum = w->active_point_sum;
        if (old_sum < new_sum) {
            replace_points(w, j, i);
            break;
        }
    } while (new_sum <= old_sum);
    return a;
}

void print_results(struct weights *w, struct analytics *a) {
    printf("Active points: ");
    for (struct node *node = linked_list_begin(w->active_points); node != linked_list_end(w->active_points); node = node->next) {
        printf("%zu ", node->value);
    }
    printf("\n");
    printf("Number of iterations: %lld\n", a->num_iterations);
}

int main(int argc, char *argv[]) {
    struct weights *w = read_input(argc, argv);
    struct analytics *a = main_loop(w);
    print_results(w, a);
    weights_free(w);
    analytics_free(a);
    return 0;
}
