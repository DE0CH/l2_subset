#include <stdlib.h>
#include <sys/mman.h>
#include "l2_subset.h"

int main(int argc, char *argv[]) {
    struct input_data data;
    void *mmaped_data;
    struct weights *w = read_from_compiled_matrix(&data, argc, argv, &mmaped_data);
    srand(data.seed);
    select_random_points(w);
    struct analytics *a = main_loop(w);
    print_results(w, a);
    analytics_free(a);
    free_mmaped_matrix(w, mmaped_data);
    w->points = malloc(0);
#if COMPUTE_MODE == USE_MATRIX
    w->entries = malloc(0);
#elif COMPUTE_MODE == USE_POINTS
#else
    #error "Invalid COMPUTE_MODE"
#endif
    weights_free(w);
    return 0;
}
