#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

struct grid {
	int n_dimensions;
	int n_points;
	double *n_coords;
	double **coord;
	int **point_index;
	// stores a random permutation, used to shuffle coordinate. Used here for speedup reasons
	int *coordinate;
};

struct initial_params {
	double **pointset;
	int npoints;
	int dim;
	int mc;
	int i_tilde;
	int trials;
};

struct kmc {
	int *k;
	int mc;
};

// we want to use C library qsort to sort.
// I made a replacement stump
void quicksort(int left, int right, double *arr);

double get_coord(struct grid *grid, int d, int i);

// The following functions use an "output variable" provided by the caller
// (anything else would just mess up the C memory management)
void round_point_up(struct grid *grid, double *point, int *output);
void round_point_down(struct grid *grid, double *point, int *output);
void round_point_extradown(struct grid *grid, double *point, int *output);

void process_coord_data(struct grid *grid, double **points, int n, int d);

double volume(struct grid *grid, int *corner);

//is x in [0,z[ ?
int open(struct grid *grid, int *x, int *z);

//is x in [0,z] ?
int closed(struct grid *grid, int *x, int *z);

int count_open(struct grid *grid, int *corner);

int count_closed(struct grid *grid, int *corner);

void grow_box_randomly(struct grid *grid, int *corner);

// Rounds box down to borders given by its contained points.
void snap_box(struct grid *grid, int *corner);

//calculates delta(x)
double get_delta(struct grid *grid, int *x);

//calculates bar(delta)(x)
double get_bar_delta(struct grid *grid, int *x);

//Generate random search point xc
//Fills coordinates (indexes, not numbers) into its three arguments
void generate_xc(struct grid *grid, int *xn_plus, int *xn_minus, int *xn_extraminus);

void generate_xc_delta(struct grid *grid, int *xn_plus);

void generate_xc_bardelta(struct grid* grid, int *xn_minus, int *xn_extraminus);

//Generate a random neighbor of xc
//k[i] == range radius in component i, mc = number of components to change
//the three xn_* variables are filled with indexes
void generate_neighbor (struct grid *grid, int *xn_plus_index, int *xn_minus_index, int *xn_extraminus_index,
			int *xc_index, int *k, int mc);

void generate_neighbor_delta(struct grid *grid, int *xn_plus_index,
			     int *xc_index, int *k, int mc);

void generate_neighbor_bardelta(struct grid *grid, int *xn_minus_index, int *xn_extraminus_index,
			int *xc_index, int *k, int mc);

double best_of_rounded_delta(struct grid *grid, int *xn_plus);

double best_of_rounded_bardelta(struct grid *grid, int *xn_minus, int *xn_extraminus, int *xc_index);

void read_points(int argc, char *argv[], struct initial_params *param);

void free_initial_params(struct initial_params *param);

void get_kmc(struct grid *grid, struct kmc *kmc, int current_iteration, int total_iteration);

void ta_update_point(double fxc, double *current, double T, int *xc_index, int *xn_best_index, int d);

void free_grid(struct grid *grid);
