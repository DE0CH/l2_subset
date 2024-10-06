#include <stdlib.h>
#include <stdio.h>
#include <sys/mman.h>
#include "l2_subset.h"
#include "mt19937-64/mt64.h"

int main(int argc, char *argv[]) {
    setvbuf(stdout, NULL, _IOLBF, 0);
    struct input_data data;
    void *mmaped_data;
    struct weights *w = read_from_compiled_matrix(&data, argc, argv, &mmaped_data);
    init_genrand64(data.seed);
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
