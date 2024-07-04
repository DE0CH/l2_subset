#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

struct node *linked_list_find(struct linked_list *list, size_t point) {
    for (struct node *node = linked_list_begin(list); node != linked_list_end(list); node = node->next) {
        if (node->value == point) {
            return node;
        }
    }
    return linked_list_end(list);
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

void print_array(double *arr, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        printf("%lf ", arr[i]);
    }
    printf("\n");
}

void print_linked_list(struct linked_list *list) {
    for (struct node *node = linked_list_begin(list); node != linked_list_end(list); node = node->next) {
        printf("%zu ", node->value);
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

void debug_cmp(struct weights *w) {
    size_t n = w->n;
    double point_weights[n];
    memcpy(point_weights, w->point_weights, n * sizeof(double));
    double active_point_sum = w->active_point_sum;
    recalculate_weights(w);
    if (!double_array_close(point_weights, w->point_weights, n)) {
        fprintf(stderr, "Active points mismatch\n");
        print_linked_list(w->active_points);
        print_matrix(w);
        print_array(point_weights, n);
        print_array(w->point_weights, n);
        exit(EXIT_FAILURE);
    }
    if (!isclose(active_point_sum, w->active_point_sum)) {
        print_linked_list(w->active_points);
        print_matrix(w);
        print_array(point_weights, n);
        print_array(w->point_weights, n);
        printf("found: %lf expected: %lf\n", active_point_sum, w->active_point_sum);
        fprintf(stderr, "Active point sum mismatch\n");
        exit(EXIT_FAILURE);
    }
    fprintf(stderr, "All good\n");
}

void add_point(struct weights *w, struct node *src) {
    size_t j = src->value;
    for (struct node *node = linked_list_begin(w->active_points); node != linked_list_end(w->active_points); node = node->next) {
        size_t i = node->value;
        double change = get_weight(w, i, j) + get_weight(w, j, i);
        w->point_weights[i] += change;
        w->active_point_sum += change;
    }
    linked_list_add(w->active_points, src->value);
    w->active_point_sum += w->point_weights[j];
    linked_list_remove(w->inactive_points, src);
    // nothing to do for src because invariant guarantees that it has value of all the active points including itself
    for (struct node *node = linked_list_begin(w->inactive_points); node != linked_list_end(w->inactive_points); node = node->next) {
        size_t i = node->value;
        double change = get_weight(w, i, j) + get_weight(w, j, i);
        w->point_weights[i] += change;
    }
#if DEBUG_SLOW
    debug_cmp(w);
#endif
}



void remove_point(struct weights *w, struct node *src) {
    size_t j = src->value;
    // need to remove itself from active list first to make sure its value doesn't get changed
    linked_list_remove(w->active_points, src);
    w->active_point_sum -= w->point_weights[j];
    for (struct node *node = linked_list_begin(w->active_points); node != linked_list_end(w->active_points); node = node->next) {
        size_t i = node->value;
        double change = get_weight(w, i, j) + get_weight(w, j, i);
        w->point_weights[i] -= change;
        w->active_point_sum -= change;
    }
    for (struct node *node = linked_list_begin(w->inactive_points); node != linked_list_end(w->inactive_points); node = node->next) {
        size_t i = node->value;
        double change = get_weight(w, i, j) + get_weight(w, j, i);
        w->point_weights[i] -= change;
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
        w->point_weights[i] = 0.0;
        for (struct node *node2 = linked_list_begin(w->active_points); node2 != linked_list_end(w->active_points); node2 = node2->next) {
            size_t j = node2->value;
            double change = get_weight(w, i, j) + get_weight(w, j, i);
            w->point_weights[i] += change;
            w->active_point_sum += change;
        }
    }
    for (struct node *node = linked_list_begin(w->inactive_points); node != linked_list_end(w->inactive_points); node = node->next) {
        size_t i = node->value;
        w->point_weights[i] = 0.0;
        for (struct node *node2 = linked_list_begin(w->active_points); node2 != linked_list_end(w->active_points); node2 = node2->next) {
            size_t j = node2->value;
            double change = get_weight(w, i, j) + get_weight(w, j, i); 
            w->point_weights[i] += change;
        }
        double change = get_weight(w, i, i) + get_weight(w, i, i);
        w->point_weights[i] += change;
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
        return -pow(2, 1-d) / (2.0 * m) * prod1 + 1.0 / (2.0 * m * m) * prod2;
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
    w->active_points = linked_list_alloc(n);
    w->inactive_points = linked_list_alloc(n);
    return w;
}

void weights_free(struct weights *w) {
    free(w->entries);
    free(w->point_weights);
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
    for (size_t i = 0; i < k; ++i) {
        mask[array[i]] = true;
    }
}

struct node *largest_active_point(struct weights *w) {
    struct node *max = NULL;
    for (struct node *node = linked_list_begin(w->active_points); node != linked_list_end(w->active_points); node = node->next) {
        size_t i = node->value;
        if (max == NULL || w->point_weights[i] >= w->point_weights[node->value]) {
            max = node;
        }
    }
    return max;
}


double relative_inactive_weight(struct weights *w, size_t inactive_point, size_t active_point) {
    return w->point_weights[inactive_point] - get_weight(w, inactive_point, active_point) - get_weight(w, active_point, inactive_point);
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
    
    if (argc != 3) {
        die("Usage: %s <points_file> <m>", argv[0]);
    }

    int m = atoi_or_die(argv[2]);

    
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
    double old_sum;
    double new_sum = w->active_point_sum;
    do {
        struct node *i = largest_active_point(w);
        size_t ii = i->value;
        struct node *j = smallest_inactive_point(w, i);
        size_t jj = j->value;
        replace_points(w, i, j);
        old_sum = new_sum;
        new_sum = w->active_point_sum;
        if (old_sum <= new_sum) {
            struct node *iii = linked_list_find(w->inactive_points, ii);
            struct node *jjj = linked_list_find(w->active_points, jj);
            replace_points(w, jjj, iii);
        }
        a->num_iterations++;
    } while (new_sum < old_sum);
    return a;
}

void print_results(struct weights *w, struct analytics *a) {
    printf("Active points: ");
    for (struct node *node = linked_list_begin(w->active_points); node != linked_list_end(w->active_points); node = node->next) {
        printf("%zu ", node->value);
    }
    printf("\n");
    printf("Number of iterations: %lld\n", a->num_iterations);
    printf("Active point sum: %lf\n", w->active_point_sum);
}

int main(int argc, char *argv[]) {
    struct weights *w = read_input(argc, argv);
    struct analytics *a = main_loop(w);
    print_results(w, a);
    weights_free(w);
    analytics_free(a);
    return 0;
}
