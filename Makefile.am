AUTOMAKE_OPTIONS = subdir-objects

BASE_INCLUDE = src/l2_subset.c src/mt19937-64/mt19937-64.c src/utils.c

bin_PROGRAMS = l2_subset gen_points l2_subset_compile_matrix l2_subset_from_compiled_matrix linf_disc l2_subset_from_compiled_matrix_w_starting_point ta_delta ta_bardelta read_matrix

l2_subset_SOURCES = $(BASE_INCLUDE) src/l2_subset_independent.c
gen_points_SOURCES = src/gen_points.c
l2_subset_compile_matrix_SOURCES = $(BASE_INCLUDE) src/l2_subset_compile_matrix.c
l2_subset_from_compiled_matrix_SOURCES = $(BASE_INCLUDE) src/l2_subset_from_compiled_matrix.c
read_matrix_SOURCES = $(BASE_INCLUDE) src/read_matrix.c 
l2_subset_from_compiled_matrix_w_starting_point_SOURCES = $(BASE_INCLUDE) src/l2_subset_from_compiled_matrix_w_starting_point.c
linf_disc_SOURCES = src/utils.c src/dem_disc.c src/dem_disc_main.c
ta_delta_SOURCES = src/TA_common.c src/TA_shrink_delta.c
ta_bardelta_SOURCES = src/TA_common.c src/TA_shrink_bardelta.c
