#include <stdio.h>
#include <stdlib.h>
#include "l2_subset.h"
#include "mt19937-64/mt64.h"

int main(int argc, char *argv[]) {
    struct input_data data;
    struct weights *w = read_point_file(&data, argc, argv);
    init_genrand64(data.seed);
    for (long long i = 0; i < data.n_trials; i++) {
        printf("Trial %lld\n", i);
        select_random_points(w);
        struct analytics *a = main_loop(w);
        print_results(w, a);
        analytics_free(a);
    }

    weights_free(w);
    return 0;
}
