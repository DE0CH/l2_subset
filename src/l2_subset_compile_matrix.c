#include <stdio.h>
#include <stdlib.h>
#include "l2_subset.h"

int main(int argc, char *argv[]) {
    struct input_data data;
    struct weights *w = read_point_file_and_save(&data, argc, argv);
    if (weights_serialize(w, data.output_filename) != 0) {
        die("Failed to serializes weights");
    }
    weights_free(w);
    return 0;
}