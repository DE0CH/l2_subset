AUTOMAKE_OPTIONS = subdir-objects

BASE_INCLUDE = src/l2_subset.c src/dem_disc.c src/TA_common.c src/mt19937-64/mt19937-64.c

bin_PROGRAMS = l2_subset gen_points l2_subset_compile_matrix l2_subset_from_compiled_matrix linf_disc l2_subset_from_compiled_matrix_w_starting_point ta_delta ta_bardelta

l2_subset_SOURCES = $(BASE_INCLUDE) src/l2_subset_independent.c
gen_points_SOURCES = src/gen_points.c
l2_subset_compile_matrix_SOURCES = $(BASE_INCLUDE) src/l2_subset_compile_matrix.c
l2_subset_from_compiled_matrix_SOURCES = $(BASE_INCLUDE) src/l2_subset_from_compiled_matrix.c
l2_subset_from_compiled_matrix_w_starting_point_SOURCES = $(BASE_INCLUDE) src/l2_subset_from_compiled_matrix_w_starting_point.c
linf_disc_SOURCES = $(BASE_INCLUDE) src/dem_disc_main.c
ta_delta_SOURCES = src/TA_common.c src/TA_shrink_delta.c
ta_bardelta_SOURCES = src/TA_common.c src/TA_shrink_bardelta.c
