/* Contains common functions (= most) for TA discrepancy search */

#include "TA_common.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include <string.h>

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))
#define TRUE 1
#define FALSE 0
typedef unsigned char bool;

// for stupid speedup reasons (alternative: make it static, and take care of initialization somehow)

#ifndef MC
#define MC 2
#endif
#define I_TILDE 316 // thresholds to be calculated (sqrt(iterations)), default value 100k
#define TRIALS 1 // nbr of runs (mean and max will be calculated), default value


// we want to use C library qsort to sort.
// I made a replacement stump
int dbl_compare(const void *a, const void *b)
{
  const double *xp=a, *yp=b;
  const double x=*xp, y=*yp;
  if (x<y)
    return -1;
  if (x>y)
    return 1;
  return 0;
} // oops: "x-y" gets cast to integer, rounded wrong

void quicksort(int left, int right, double *arr) {
  qsort(&arr[left], right-left, sizeof(double), dbl_compare);
}

double get_coord(struct grid *grid, int d, int i) {
  return grid->coord[d][i];
}

int get_index_up(struct grid *grid, int d, double key) {
  int i=(int)key*(grid->n_coords[d]-1);
  while (grid->coord[d][i] < key)
    i++;
  while ((i>0) && (grid->coord[d][i-1] >= key))
    i--;
  // termination: at first coordinate which is >= key
  // (if i=0 then key==0)
  return i;
}

int get_index_down(struct grid *grid, int d, double key){
  int bound=grid->n_coords[d]-1;
  int i=key*bound;
  while (grid->coord[d][i] > key)
    i--;
  while ((i < bound) && (grid->coord[d][i+1] <= key))
    i++;
  return i;
}

// The following two functions use an "output variable" provided by the caller
// (anything else would just mess up the C memory management)
void round_point_up(struct grid *grid, double *point, int *output) {
  int i;
  for (i=0; i<grid->n_dimensions; i++)
    output[i] = get_index_up(grid, i, point[i]);
}

void round_point_down(struct grid *grid, double *point, int *output) {
  int i;
  for (i=0; i<grid->n_dimensions; i++)
    output[i] = get_index_down(grid, i, point[i]);
}

void round_point_extradown(struct grid *grid, double *point, int *output) {
  int i;
  for (i=0; i<grid->n_dimensions; i++) {
    output[i] = get_index_down(grid, i, point[i]);
    if(output[i]==0)
      output[i]=grid->n_coords[i]-1;
  }
}

void process_coord_data(struct grid *grid, double **points, int n, int d) {
  int i,j;
  double tmp_coords[n+2];
  int n_uniq, idx;
  double prev_coord;
  //initialise n_coords[], coord[][]
  grid->n_dimensions=d;
  grid->n_points=n;
  grid->coordinate=malloc(d*sizeof(int));
  for (i=0; i<d; i++)
    grid->coordinate[i]=i;
  grid->n_coords = malloc(grid->n_dimensions*sizeof(double));
  grid->coord = malloc(grid->n_dimensions*sizeof(double *));
  for (i=0; i<grid->n_dimensions; i++) {
    for (j=0; j<grid->n_points; j++)
      tmp_coords[j+1]=points[j][i];
    tmp_coords[0]=0.0;
    tmp_coords[n+1]=1.0;
    quicksort(1, n+1, tmp_coords);
    // 1. count
    n_uniq=1;
    prev_coord=0;
    for (j=1; j<=n+1; j++) { // inclusive bound: contains 1.0
      if (prev_coord==tmp_coords[j])
	      continue;
      n_uniq++;
      prev_coord=tmp_coords[j];
    }
    // 2. transfer
    grid->coord[i]=malloc(n_uniq*sizeof(double));
    idx=1;
    prev_coord=tmp_coords[0];
    grid->coord[i][0]=prev_coord;
    for (j=1; j<=n+1; j++) {
      if (prev_coord==tmp_coords[j])
	      continue;
      prev_coord=tmp_coords[j];
      grid->coord[i][idx++] = prev_coord;
    }
    grid->n_coords[i]=n_uniq;
    //    fprintf(stderr, "Coordinates %d, %d: ", i, n_uniq);
    //    for (j=0; j<n_uniq; j++)
    //      fprintf(stderr, "%g ", coord[i][j]);
    //    fprintf(stderr,"\n");
  }
  // finished setup for: n_coords[], coord[][]

  // next: transfer point set to into point_index
  grid->point_index=malloc(grid->n_points*sizeof(int *));
  for (i=0; i<grid->n_points; i++)
    grid->point_index[i]=malloc(grid->n_dimensions*sizeof(int));
  for (i=0; i<grid->n_points; i++)
    for (j=0; j<grid->n_dimensions; j++) {
      idx=get_index_up(grid, j, points[i][j]);
      if (grid->coord[j][idx] != points[i][j]) {
        fprintf(stderr, "ERROR: located incorrect coordinate (%g at %d, wanted %g).\n", grid->coord[j][idx], idx, points[i][j]);
        abort();
      }
      grid->point_index[i][j] = idx;
    }
  // setup finished.
}

double volume(struct grid *grid, int *corner)
{
  double vol=1.0;
  int i;
  for (i=0; i<grid->n_dimensions; i++)
    vol *= get_coord(grid, i, corner[i]);
  return vol;
}


//is x in [0,z[ ?
int open(struct grid *grid, int *x, int *z)
{
  int i;
  for (i=0; i<grid->n_dimensions; i++)
    if (x[i] >= z[i])
      return 0;
  return 1;
}

//is x in [0,z] ?
int closed(struct grid *grid, int *x, int *z)
{
  int i;
  for (i=0; i<grid->n_dimensions; i++)
    if (x[i] > z[i])
      return 0;
  return 1;
}

int count_open(struct grid *grid, int *corner)
{
  int n=grid->n_points;
  int i;
  int res=0;
  for (i=0; i<n; i++)
    res += open(grid, grid->point_index[i], corner);
  return res;
}

int count_closed(struct grid *grid, int *corner)
{
  int n=grid->n_points;
  int i;
  int res=0;
  for (i=0; i<n; i++)
    res += closed(grid, grid->point_index[i], corner);
  return res;
}

// group of functions for grow/"snap up"
int point_critical_dimension(struct grid *grid, int *corner, int *point)
{
  int i;
  int crit=-1;
  //  fprintf(stderr, "Point");
  //  for (i=1; i<=d; i++)
  //    fprintf(stderr, " %g", point[i]);
  for (i=0; i<grid->n_dimensions; i++)
    if (point[i] >= corner[i]) {
      if (crit>=0) {
	//	fprintf(stderr, " double out (%d,%d)\n", crit, i);
	return -1;
      }
      else
	crit=i;
    }
  //  fprintf(stderr, " crit %d\n", crit);
  return crit;
}

void grow_box_randomly(struct grid *grid, int *corner)
{
  int order[grid->n_dimensions];
  int n=grid->n_points, d=grid->n_dimensions;
  int i,j,k, swap, memo;
  int curr_d;
  int new_box[d];

  for (i=0; i<d; i++)
    order[i]=i;
  for (i=0; i<d; i++) {
    j = i + random()%(d-i);
    if (i != j) {
      swap=order[i];
      order[i]=order[j];
      order[j]=swap;
    }
  }
  for (i=0; i<d; i++)
    new_box[i] = grid->n_coords[i]-1;

  for (i=0; i<n; i++) {
    memo=-1;
    for (j=0; j<d; j++) {
      curr_d=order[j];
      k=grid->point_index[i][curr_d];
      if (k >= corner[curr_d]) {
        if (k < new_box[curr_d]) {
          if (memo<0)
            memo=curr_d;
        }
        else {
          memo=-1;
          break;
        }
      }
    }
    if (memo >= 0)
      new_box[memo] = grid->point_index[i][memo];
  }


  for (i=0; i<d; i++)
    corner[i]=new_box[i];
}

// Rounds box down to borders given by its contained points.
void snap_box(struct grid *grid, int *corner)
{
  int d=grid->n_dimensions;
  int n=grid->n_points;
  int max_idx[d];
  int i,j;
  for (i=0; i<d; i++)
    max_idx[i]=-1;
  for (i=0; i<n; i++)
    if (closed(grid, grid->point_index[i], corner))
      for (j=0; j<d; j++)
	      if (grid->point_index[i][j] > max_idx[j])
	        max_idx[j]=grid->point_index[i][j];
        for (i=0; i<d; i++)
          corner[i] = (max_idx[i] < 0) ? 0 : max_idx[i];
}



//calculates delta(x)
double get_delta(struct grid *grid, int *x)
{
  int op;
  double vol, delta;
  int n=grid->n_points;

  //calculates no of points in [0,x[
  op = count_open(grid, x);

  //calcualtes volume of box generated by x
  vol = volume(grid, x);

  //calculates delta
  delta = vol - (double)op/n;
  //  fprintf(stderr, "Stand.: vol %g pts %d -> %g error %g\n",
  //	  vol, op, (double)op/n, delta);
  return delta;
}

//calculates bar(delta)(x)
double get_bar_delta(struct grid *grid, int *x)
{
  int cl;
  int n=grid->n_points;
  double vol, bdelta;

  //calculates no of points in [0,x]
  cl = count_closed(grid, x);

  //calcualtes volume of box generated by x
  vol = volume(grid, x);

  //calculates bar(delta)
  bdelta = (double)cl/n - vol;
  //  fprintf(stderr, "Stand.: vol %g pts %d -> %g error %g\n",
  //	  vol, cl, (double)cl/n, bdelta);
  return bdelta;
}

//Generate random search point xc
//Fills coordinates (indexes, not numbers) into its three arguments
void generate_xc(struct grid *grid, int *xn_plus, int *xn_minus, int *xn_extraminus)
{
  int j, d=grid->n_dimensions;
  double xn[d];
  double temp;
  for(j=0; j<d; j++) {
    temp=(double)((double)rand()/RAND_MAX);
    xn[j]=pow(temp,(double)((double)1/(double)d));
  }
  round_point_up(grid, xn, xn_plus);
  round_point_down(grid, xn, xn_minus);
  round_point_extradown(grid, xn, xn_extraminus);
}

void generate_xc_delta(struct grid *grid, int *xn_plus)
{
  int j, d=grid->n_dimensions;
  double xn[d];
  double temp;
  for(j=0; j<d; j++) {
    temp=(double)((double)rand()/RAND_MAX);
    xn[j]=pow(temp,(double)((double)1/(double)d));
  }
  round_point_up(grid, xn, xn_plus);
}

void generate_xc_bardelta(struct grid *grid, int *xn_minus, int *xn_extraminus)
{
  int j, d=grid->n_dimensions;
  double xn[d];
  double temp;
  for(j=0; j<d; j++) {
    temp=(double)((double)rand()/RAND_MAX);
    xn[j]=pow(temp,(double)((double)1/(double)d));
  }
  round_point_down(grid, xn, xn_minus);
  round_point_extradown(grid, xn, xn_extraminus);
}


//Generate a random neighbor of xc
//k[i] == range radius in component i, mc = number of components to change
//the three xn_* variables are filled with indexes
void generate_neighbor (struct grid *grid, int *xn_plus_index, int *xn_minus_index, int *xn_extraminus_index,
			int *xc_index, int *k, int mc)
{
  int i, j, q, d=grid->n_dimensions;
  double temp, upper_bound, lower_bound;
  double xn[d];

  //First copy the values of the current search point
  for(j=0; j<d; j++) {
    xn[j] = grid->coord[j][xc_index[j]];
  }

  // find mc different coordinates to be changed
  for (j=0; j<mc; j++) {
    i = j + random()%(d-j);
    if (i != j) {
      q = grid->coordinate[j];
      grid->coordinate[j]=grid->coordinate[i];
      grid->coordinate[i] = q;
    }
  }

  //set lower and upper bound to the box from which the random neighbor will be sampled
  for(j=0; j<mc; j++){
    if (xc_index[grid->coordinate[j]]-k[grid->coordinate[j]]>=0)
      lower_bound = grid->coord[grid->coordinate[j]][xc_index[grid->coordinate[j]]-k[grid->coordinate[j]]];
    else
      lower_bound=0.0;
    if (xc_index[grid->coordinate[j]]+k[grid->coordinate[j]] < grid->n_coords[grid->coordinate[j]])
      upper_bound = grid->coord[grid->coordinate[j]][xc_index[grid->coordinate[j]]+k[grid->coordinate[j]]];
    else
      upper_bound=1.0;

    //draw a random number in [0,1]
    temp=(double)((double)rand()/RAND_MAX);
    temp=(pow(upper_bound,d)-pow(lower_bound,d))*temp + pow(lower_bound,d);
    xn[grid->coordinate[j]]=pow(temp,(double)((double)1/(double)d));
  }

  round_point_up(grid, xn, xn_plus_index);
  round_point_down(grid, xn, xn_minus_index);
  round_point_extradown(grid, xn, xn_extraminus_index);
}

void generate_neighbor_delta(struct grid *grid, int *xn_plus_index, int *xc_index, int *k, int mc) {
  int i, j, q, d=grid->n_dimensions;
  double temp, upper_bound, lower_bound;
  double xn[d];

  //First copy the values of the current search point
  for(j=0; j<d; j++) {
    xn[j] = grid->coord[j][xc_index[j]];
  }

  // find mc different coordinates to be changed
  for (j=0; j<mc; j++) {
    i = j + random()%(d-j);
    if (i != j) {
      q = grid->coordinate[j];
      grid->coordinate[j]=grid->coordinate[i];
      grid->coordinate[i] = q;
    }
  }

  //set lower and upper bound to the box from which the random neighbor will be sampled
  for(j=0; j<mc; j++){
    if (xc_index[grid->coordinate[j]]-k[grid->coordinate[j]]>=0)
      lower_bound = grid->coord[grid->coordinate[j]][xc_index[grid->coordinate[j]]-k[grid->coordinate[j]]];
    else
      lower_bound=0.0;
    if (xc_index[grid->coordinate[j]]+k[grid->coordinate[j]] < grid->n_coords[grid->coordinate[j]])
      upper_bound = grid->coord[grid->coordinate[j]][xc_index[grid->coordinate[j]]+k[grid->coordinate[j]]];
    else
      upper_bound=1.0;

    //draw a random number in [0,1]
    temp=(double)((double)rand()/RAND_MAX);
    temp=(pow(upper_bound,d)-pow(lower_bound,d))*temp + pow(lower_bound,d);
    xn[grid->coordinate[j]]=pow(temp,(double)((double)1/(double)d));
  }

  round_point_up(grid, xn, xn_plus_index);
}

void generate_neighbor_bardelta(struct grid *grid, int *xn_minus_index, int *xn_extraminus_index,
			int *xc_index, int *k, int mc) {
  int i, j, q, d=grid->n_dimensions;
  double temp, upper_bound, lower_bound;
  double xn[d];

  //First copy the values of the current search point
  for(j=0; j<d; j++) {
    xn[j] = grid->coord[j][xc_index[j]];
  }

  // find mc different coordinates to be changed
  for (j=0; j<mc; j++) {
    i = j + random()%(d-j);
    if (i != j) {
      q = grid->coordinate[j];
      grid->coordinate[j]=grid->coordinate[i];
      grid->coordinate[i] = q;
    }
  }

  //set lower and upper bound to the box from which the random neighbor will be sampled
  for(j=0; j<mc; j++){
    if (xc_index[grid->coordinate[j]]-k[grid->coordinate[j]]>=0)
      lower_bound = grid->coord[grid->coordinate[j]][xc_index[grid->coordinate[j]]-k[grid->coordinate[j]]];
    else
      lower_bound=0.0;
    if (xc_index[grid->coordinate[j]]+k[grid->coordinate[j]] < grid->n_coords[grid->coordinate[j]])
      upper_bound = grid->coord[grid->coordinate[j]][xc_index[grid->coordinate[j]]+k[grid->coordinate[j]]];
    else
      upper_bound=1.0;

    //draw a random number in [0,1]
    temp=(double)((double)rand()/RAND_MAX);
    temp=(pow(upper_bound,d)-pow(lower_bound,d))*temp + pow(lower_bound,d);
    xn[grid->coordinate[j]]=pow(temp,(double)((double)1/(double)d));
  }

  round_point_down(grid, xn, xn_minus_index);
  round_point_extradown(grid, xn, xn_extraminus_index);
}

// Computes the best of the rounded points -- basic version
// Constant neighbourhood size and mc-values.
// Does not split the search.
// Copies the appropriate "thing" into xc_index (output variable)
double best_of_rounded_delta(struct grid *grid, int *xn_plus)
{
  double fxc;
  int j, d = grid->n_dimensions;
  int xn_plus_grow[d];

  // Growing, shrinking.
  // Grower, shrinker that copy the point
  for (j = 0; j < d; j++)
    xn_plus_grow[j] = xn_plus[j];
  grow_box_randomly(grid, xn_plus_grow);

  // Now, create the official numbers.
  // official update from modified points
  fxc = get_delta(grid, xn_plus_grow);

  return fxc;
}

// Computes the best of the rounded points -- basic version
// Constant neighbourhood size and mc-values.
// Does not split the search.
// Copies the appropriate "thing" into xc_index (output variable)
double best_of_rounded_bardelta(struct grid *grid, int *xn_minus, int *xn_extraminus, int *xc_index)
{
  double fxn_extraminus;
  double fxc;
  int j, d = grid->n_dimensions;
  int use_extraminus = 0;
  int xn_minus_snap[d], xn_extraminus_snap[d];
  for (j = 0; j < d; j++)
    if (xn_minus[j] != xn_extraminus[j])
    {
      use_extraminus = 1;
      break;
    }

  // Growing, shrinking.
  // Grower, shrinker that copy the point
  for (j = 0; j < d; j++)
    xn_minus_snap[j] = xn_minus[j];
  snap_box(grid, xn_minus_snap);
  if (use_extraminus)
  {
    for (j = 0; j < d; j++)
      xn_extraminus_snap[j] = xn_extraminus[j];
    snap_box(grid, xn_extraminus_snap);
  }

  // Now, create the official numbers.
  // official update from modified points
  fxc = get_bar_delta(grid, xn_minus_snap);
  if (use_extraminus)
  {
    fxn_extraminus = get_bar_delta(grid, xn_extraminus_snap);
    fxc = max(fxc, fxn_extraminus);
  }

  // Remains only to copy the winning point to output variable xc_index.
  if (use_extraminus && (fxn_extraminus >= fxc))
  {
    for (j = 0; j < d; j++)
      xc_index[j] = xn_extraminus[j];
  }
  else
  {
    for (j = 0; j < d; j++)
      xc_index[j] = xn_minus[j];
  }

  return fxc;
}

void read_points(int argc, char *argv[], struct initial_params *param) {
  param->mc = MC;
  param->i_tilde = I_TILDE;
  param->trials = TRIALS;
  FILE *pointfile;
  bool should_close_pointfile = TRUE;
  int pos = 1;

  FILE *random;
  unsigned int seed;
  random = fopen("/dev/random", "rb");
  if (fread(&seed, sizeof(int), 1, random) != 1) {
    fprintf(stderr, "Could not read random seed\n");
    exit(EXIT_FAILURE);
  }
  fclose(random);
  srand(seed);
  while (pos < argc)
  {
    if (!strcmp(argv[pos], "-mc"))
    {
      param->mc = atoi(argv[++pos]);
      pos++;
      fprintf(stderr, "Using mc = %d\n", param->mc);
    }
    else if (!strcmp(argv[pos], "-iter"))
    {
      param->i_tilde = (int)sqrt(atoi(argv[++pos]));
      pos++;
      fprintf(stderr, "Using %d iterations (adj. for sqrt)\n",
              param->i_tilde * param->i_tilde);
    }
    else if (!strcmp(argv[pos], "-trials"))
    {
      param->trials = atoi(argv[++pos]);
      pos++;
      fprintf(stderr, "Doing %d independent trials (currently: times ten thresh. rep.)\n",
              param->trials);
    }
    else
      break;
  }
  switch (argc - pos)
  {
  case 0:
    int i = scanf("%d %d reals\n", &param->dim, &param->npoints);
    if (i != 2)
    {
      fprintf(stderr, "stdin mode and header line not present\n");
      exit(EXIT_FAILURE);
    }
    pointfile = stdin;
    should_close_pointfile = FALSE;
    break;

  case 1: // one arg, interpret as file name
    pointfile = fopen(argv[pos], "r");
    if (pointfile == NULL)
    {
      fprintf(stderr, "Could not open file %s\n", argv[pos]);
      exit(EXIT_FAILURE);
    }
    i = fscanf(pointfile, "%d %d reals\n", &param->dim, &param->npoints);
    if (i != 2)
    {
      fprintf(stderr, "stdin mode and header line not present\n");
      exit(EXIT_FAILURE);
    }
    break;

  case 2: // interpret as dim npoints args
    param->dim = atoi(argv[pos++]);
    param->npoints = atoi(argv[pos]);
    pointfile = stdin;
    should_close_pointfile = FALSE;
    break;

  case 3: // interpret as dim npoints file; file not allowed to have header
    param->dim = atoi(argv[pos++]);
    param->npoints = atoi(argv[pos++]);
    pointfile = fopen(argv[pos], "r");
    if (pointfile == NULL)
    {
      fprintf(stderr, "Could not open file %s\n", argv[pos]);
      exit(EXIT_FAILURE);
    }
    break;

  default:
    fprintf(stderr, "Usage: calc_discr [dim npoints] [file]\n\nIf file not present, read from stdin. If dim, npoints not present, \nassume header '%%dim %%npoints reals' (e.g. '2 100 reals') in file.\n");
    exit(EXIT_FAILURE);
  }

  fprintf(stderr, "Reading dim %d npoints %d\n", param->dim, param->npoints);
  param->pointset = malloc(param->npoints * sizeof(double *));
  for (int i = 0; i < param->npoints; i++)
  {
    param->pointset[i] = malloc(param->dim * sizeof(double));
    for (int j = 0; j < param->dim; j++)
    {
      if (fscanf(pointfile, "%lg ", &(param->pointset[i][j])) == EOF) {
        fprintf(stderr, "Unexpected EOF when reading the pointfile\n");
        exit(EXIT_FAILURE);
      }
      // newline counts as whitespace
    }
  }
  if (should_close_pointfile) {
    fclose(pointfile);
  }
  if (param->dim < param->mc)
    param->mc = param->dim;
}

void get_kmc(struct grid *grid, struct kmc *kmc, int current_iteration, int total_iteration) {
  // Update k-value
  for (int j = 0; j < grid->n_dimensions; j++)
  {
    int start = (int)((grid->n_coords[j] - 1) / 2);
    kmc->k[j] = start * (((double)total_iteration - current_iteration) / (total_iteration)) +
            1 * ((double)current_iteration / (total_iteration));
    //	    k[j]=start[j] - (int)((3.0/4)*(current_iteration/outerloop)*(start[j]-1));
  }

  // Update mc-value
  kmc->mc = 2 + (int)(current_iteration / total_iteration * (grid->n_dimensions - 2));
}

void ta_update_point(double fxc, double *current, double T, int *xc_index, int *xn_best_index, int d) {
  // Update of current best value if necessary
  if (fxc - *current >= T)
  {
    *current = fxc;
    for (int j = 0; j < d; j++)
    {
      xc_index[j] = xn_best_index[j];
    }
  }
}

void free_initial_params(struct initial_params *param) {
  for (int i = 0; i < param->npoints; i++)
  {
    free(param->pointset[i]);
  }
  free(param->pointset);
}

void free_grid(struct grid *gird) {
  for (int i = 0; i < gird->n_dimensions; i++)
  {
    free(gird->coord[i]);
  }
  free(gird->coord);
  for (int i = 0; i < gird->n_points; i++)
  {
    free(gird->point_index[i]);
  }
  free(gird->point_index);
  free(gird->n_coords);
  free(gird->coordinate);
}
