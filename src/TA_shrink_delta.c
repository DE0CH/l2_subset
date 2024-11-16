/*
   reads a point set and computes a lower bound on its star discrepancy
*/

/*
  i_tilde inner and outer loops (effect of alpha is considered to be negligible)
  threshold aus Paaren von Nachbarn berechnet
*/

// This time, the file is written so that first element is x[0] :-p

// compute threshold sequence once, then reuse it
// not using alpha

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include <string.h>

#include "TA_common.h"

double oldmain(struct grid *grid, double **pointset, int n, int d, int i_tilde, int trials)
{
  double thresh[i_tilde]; // Thresholdsequence

  int xc_index[d], xn_plus_index[d];
  double current, global[trials], best; // current and global best values
  struct kmc kmc;
  int _k[d];
  kmc.k = _k;

  // Sort the grid points, setup global variables
  process_coord_data(grid, pointset, n, d);

  // Algorithm starts here
  for (int t = 0; t < trials; t++)
  { // Initialization
    // Generate threshold sequence
    for (int i = 0; i < i_tilde; i++)
    {
      int current_iteration = i + 1;
      // Update k-value
      get_kmc(grid, &kmc, current_iteration, i_tilde);

      // generation of random point xc
      generate_xc_delta(grid, xc_index);

      //(Possibly) Snaps the point upwards and computes the fitness
      current = best_of_rounded_delta(grid, xc_index);

      // draw a neighbour of xc
      generate_neighbor_delta(grid, xn_plus_index, xc_index, kmc.k, kmc.mc);

      // Compute the threshold
      double fxc = best_of_rounded_delta(grid, xn_plus_index);
      thresh[i] = 0.0 - fabs(fxc - current);
    }

    // sort the thresholds in increasing order
    quicksort(0, i_tilde, thresh);

    current = 0;
    global[t] = 0;
    // draw a random initial point
    generate_xc_delta(grid, xc_index);

    //(Possibly) Snap and compute the best of the rounded points and update current value
    current = best_of_rounded_delta(grid, xc_index);

    global[t] = current;

    for (int i = 0; i < i_tilde; i++)
    {
      double T = thresh[i];

      for (int p = 0; p < i_tilde ; p++)
      {
        int current_iteration = i*i_tilde + p + 1;


        // Update k-value
        // Update mc-value
        get_kmc(grid, &kmc, current_iteration, i_tilde * i_tilde);

        // Get random neighbor
        generate_neighbor_delta(grid, xn_plus_index, xc_index, kmc.k, kmc.mc);

        //(Possibly) Snap the points and compute the best of the rounded points
        double fxc = best_of_rounded_delta(grid, xn_plus_index);
        // Global update if necessary
        if (fxc > global[t])
        {
          global[t] = fxc;
        }
        ta_update_point(fxc, &current, T, xc_index, xn_plus_index, d);
      } // innerloop
    } // outerloop
  } // trials

  // best calculated value
  best = 0;
  for (int t = 0; t < trials; t++)
  {
    if (global[t] > best)
    {
      best = global[t];
    }
  }
  return best;
}

int main(int argc, char **argv)
{
  struct grid grid;
  struct initial_params param;
  fprintf(stderr, "Calling Carola calculation\n");
  read_points(argc, argv, &param);
  printf("%g\n", oldmain(&grid, param.pointset, param.npoints, param.dim, param.i_tilde, param.trials));
  free_initial_params(&param);
  free_grid(&grid);
  return EXIT_SUCCESS;
}
