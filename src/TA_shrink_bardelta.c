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

int main(int argc, char **argv)
{
  struct initial_params param;
  fprintf(stderr, "Calling Carola calculation\n");
  read_points(argc, argv, &param);
  printf("%g\n", ta_bardelta(param.pointset, param.npoints, param.dim, param.i_tilde, param.trials));
  free_initial_params(&param);
  return EXIT_SUCCESS;
}
