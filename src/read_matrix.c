#include <stdlib.h>
#include <stdio.h>
#include "l2_subset.h"
#include "utils.h"
int main() {
    void *mmapped_data;
    struct weights *w = weights_deserialize("../src/t.bin", &mmapped_data);
    if (w == NULL) {
        die("Could not read compiled matrix file");
    }
    printf("n: %lld\n", w->n);
    printf("m: %lld\n", w->m);
    printf("d: %lld\n", w->d);
    for(size_t i = 0; i<w->n; i++) {
        for(size_t j=0; j<w->n; j++) {
            printf("(%zu, %zu): %.3lf\n", i, j, get_weight(w, i, j));
        }
    }
    free_mmaped_matrix(w, mmapped_data);
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