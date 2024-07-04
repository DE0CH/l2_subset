#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_qrng.h>
int main(int argc, char **argv)
{
    int i, j;
    int dim;
    char *endptr;

    if (argc != 4) {
        fprintf(stderr, "Usage: %s <dim> <output_file> <nb_points>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    dim = strtol(argv[1], &endptr, 10);
    if (*argv[1] == '\0' || *endptr != '\0')
    {
        fprintf(stderr, "Invalid number: %s\n", argv[1]);
        exit(EXIT_FAILURE);
    };

    int nb;

    FILE *fp;
    fp = fopen(argv[2], "w"); // Careful with dim
    if (fp == NULL) {
        printf("Error opening file!\n");
        exit(EXIT_FAILURE);
    }

    char *endptr2;
    nb = strtol(argv[3], &endptr2, 10);
    if (*argv[3] == '\0' || *endptr2 != '\0')
    {
        fprintf(stderr, "Invalid number: %s\n", argv[3]);
        exit(EXIT_FAILURE);
    };

    gsl_qrng *q = gsl_qrng_alloc(gsl_qrng_sobol, dim);

    if (q == NULL) {
        printf("Error creating quasi-random number generator!\n");
        exit(EXIT_FAILURE);
    }

    fprintf(fp, "%d %d 0.0\n", dim, nb);
    for (i = 0; i < nb; i++)
    {
        double v[dim];
        gsl_qrng_get(q, v);
        for (j = 0; j < dim; j++)
        {
            // printf ("%.5f %.5f \n", v[0], v[1]);
            fprintf(fp, "%lf ", v[j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    gsl_qrng_free(q);
    return 0;
}