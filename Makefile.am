AUTOMAKE_OPTIONS = subdir-objects

BASE_INCLUDE = src/l2_subset.c src/dem_disc.c src/mt19937-64/mt19937-64.c

bin_PROGRAMS = l2_subset gen_points l2_subset_compile_matrix l2_subset_from_compiled_matrix linf_disc l2_subset_from_compiled_matrix_w_starting_point

l2_subset_SOURCES = $(BASE_INCLUDE) src/l2_subset_independent.c
gen_points_SOURCES = src/gen_points.c
l2_subset_compile_matrix_SOURCES = $(BASE_INCLUDE) src/l2_subset_compile_matrix.c
l2_subset_from_compiled_matrix_SOURCES = $(BASE_INCLUDE) src/l2_subset_from_compiled_matrix.c
l2_subset_from_compiled_matrix_w_starting_point_SOURCES = $(BASE_INCLUDE) src/l2_subset_from_compiled_matrix_w_starting_point.c
linf_disc_SOURCES = $(BASE_INCLUDE) src/dem_disc_main.c
