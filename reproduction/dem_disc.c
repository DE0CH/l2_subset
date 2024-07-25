// Taken from MW's dem_main.c, then adapted to avoid making a mess in the original file
// Input: file name
// In file, first line has to be dim,npoints,kpoints

// expects each point on its own line with real positions.
// expects first line to be "dim npoints reals"--from thiemard
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include "l2_subset.h"

struct global
{
    double globallower;
    double glob_bound;
};

// Split this file into main and shift parts

int cmpdbl(const void *a, const void *b)
{
    if ((*(double *)a) < (*(double *)b))
        return -1;
    else if ((*(double *)a) == (*(double *)b))
        return 0;
    return 1;
}

int cmpkeyk(const void *pt1, const void *pt2, void *arg)
{
    int comparedim = *(int *)arg;
    double a = (*(double **)pt1)[comparedim], b = (*(double **)pt2)[comparedim];
    if (a < b)
        return -1;
    else if (a > b)
        return 1;
    return 0;
}

void usage()
{
    fprintf(stderr, "Usage: dem_discr [dim npoints] [file]\n\nIf file not present, read from stdin. If dim, npoints not present, \nassume header '%%dim %%npoints reals' (e.g. '2 100 reals') in file.\n");
}

/////////////////////////////////////////////////////////

double oydiscr_cell(int npoints, int dim, int rempoints,
                    double **forced, int nforced,
                    double *lowerleft, double *upperright, struct global *g)
{
    double discr = (double)(rempoints + nforced) / npoints;
    double maxdiscr = (double)rempoints / npoints;
    if (discr > maxdiscr)
    {
        maxdiscr = discr;
    }
    if (maxdiscr < g->globallower)
    {
        return maxdiscr;
    }
    // quicker code for use in some easy cells

    maxdiscr = 0;

    return maxdiscr;
}

// "forced" points are points that are strictly between the boundaries in
// one of the cdim earlier dimensions; pointset[0--rempoints-1] are points that
// so far are at most at the lower-left corner in every dimension.
// this includes a final run with lower-left=1.
// being ON a border changes nothing:
//   ON lower-left counts as in (including when lowerleft=1)
//   ON upper-right counts as out (except if previous).
double oydiscr_int(double **pointset, int npoints, int dim, int rempoints,
                   double **forced, int nforced, int cdim,
                   double *lowerleft, double *upperright, struct global *g)
{
    double coord, forcedcoord, lowedge = 0.0, highedge;
    int newcount = 0, forcedidx, i, j, previdx = 0, newforcedidx, resort = 0, curridx;
    int newrempoints, wasfinal = 0;
    double maxdiscr = 0.0, discr;
    double **newforced = malloc((nforced + rempoints) * sizeof(double *));
    // internal vars: previdx points at first element excluded from last pass
    //                    (on or above border coordinate)
    //                forcedidx: next unused forced element (a.k.a. counter)
    //                curridx: current pointset point ("i" as loop var)
    //  coords:
    //           coord is value of current pointset-point,
    //           forcedcoord is value of next-up forced boundary,
    //           lowedge is value of last boundary we used
    // newcount counts number of pointset-points since last coord
    if (cdim == dim)
    {
        // fprintf(stderr, "Reached a cell\n");
        free(newforced);
        discr = oydiscr_cell(npoints, dim, rempoints,
                             forced, nforced,
                             lowerleft, upperright, g);
        if (discr > g->glob_bound)
        {
            g->glob_bound = discr;
        }
        return discr;
    }

    int comparedim = cdim;
    qsort_r(pointset, rempoints, sizeof(double *), cmpkeyk, &comparedim);
    qsort_r(forced, nforced, sizeof(double *), cmpkeyk, &comparedim);
    i = 0;
    forcedidx = 0;
    while ((i < rempoints) || (forcedidx < nforced))
    {
        if (i < rempoints)
            coord = pointset[i][cdim];
        else
            coord = 1.0;
        if (forcedidx < nforced)
            forcedcoord = forced[forcedidx][cdim];
        else
            forcedcoord = 1.0;
        if ((forcedcoord > coord) && (newcount <= sqrt(npoints)))
        {
            newcount++;
            i++;
            if ((i < rempoints) || (forcedidx < nforced))
                continue;
            else
            { // Add one trailing cell
                lowerleft[cdim] = lowedge;
                highedge = upperright[cdim] = 1.0;
                wasfinal = 1;
            }
        } // below: create new cell
        if (!wasfinal)
        {
            if (forcedcoord <= coord)
            {
                lowerleft[cdim] = lowedge;
                highedge = upperright[cdim] = forcedcoord;
            }
            else
            { // must be count-based border
                lowerleft[cdim] = lowedge;
                highedge = upperright[cdim] = coord;
            }
        } // end "if (!wasfinal)"
        curridx = i; // for better mnemonics
        // creating a new cell (or subslab):
        // 1. surviving forced copied
        for (j = 0; (j < nforced) && (forced[j][cdim] < highedge); j++)
            newforced[j] = forced[j];
        newforcedidx = j;
        // 2. new (strictly) internal points appended as forced
        j = previdx;
        while ((j < rempoints) && (pointset[j][cdim] <= lowedge))
            j++;
        newrempoints = j;
        for (; (j < rempoints) && (pointset[j][cdim] < highedge); j++)
            newforced[newforcedidx++] = pointset[j];
        if (j > (curridx + 1))
            resort = 1;
        // 3. make call with properly adjusted boundaries, update variables
        discr = oydiscr_int(pointset, npoints, dim, newrempoints,
                            newforced, newforcedidx, cdim + 1,
                            lowerleft, upperright, g);
        if (discr > maxdiscr)
        {
            maxdiscr = discr;
        }
        if (resort)
        {
            comparedim = cdim;
            qsort_r(pointset, rempoints, sizeof(double *), cmpkeyk, &comparedim);
            resort = 0;
        }
        while ((forcedidx < nforced) && (forced[forcedidx][cdim] == highedge))
            forcedidx++;
        while ((i < rempoints) && (pointset[i][cdim] <= highedge))
            i++;
        lowedge = highedge;
        previdx = i;
        newcount = 0;
    }
    // one final call to capture the border cases (for boxes containing a point with coord 1.0)
    // 1. new forced == old forced (copy to avoid interfering with previous stages)
    for (j = 0; j < nforced; j++)
        newforced[j] = forced[j];
    // 2. per above, we have no new internal/forced points
    // 3. make the call
    lowerleft[cdim] = lowedge;
    upperright[cdim] = 1.0;
    discr = oydiscr_int(pointset, npoints, dim, rempoints,
                        newforced, nforced, cdim + 1,
                        lowerleft, upperright, g);
    if (discr > maxdiscr)
    {
        maxdiscr = discr;
    }

    free(newforced);
    return maxdiscr;
}

double oydiscr(double **pointset, int dim, int npoints)
{
    struct global g;
    g.globallower = 0.0;
    double lowerleft[dim], upperright[dim];
    double **pre_force = malloc(2 * dim * sizeof(double *));
    double discr, *border;
    double maxcoord;
    int maxpos;
    int is_border[npoints];
    double **clone_set = malloc(npoints * sizeof(double *));
    int i, j, k;
    // fprintf(stderr,"Going to int");
    for (i = 0; i < dim; i++)
    {
        border = malloc(dim * sizeof(double));
        for (j = 0; j < dim; j++)
        {
            if (i == j)
                border[j] = 1.0;
            else
                border[j] = 0.0;
        }
        pre_force[i] = border;
    }
    for (i = 0; i < npoints; i++)
        is_border[i] = 0;
    for (i = 0; i < dim; i++)
    {
        maxcoord = -1.0;
        maxpos = -1;
        for (j = 0; j < npoints; j++)
            if (pointset[j][i] > maxcoord)
            {
                maxcoord = pointset[j][i];
                maxpos = j;
            }
        is_border[maxpos] = 1;
    }
    j = dim;
    k = 0;
    for (i = 0; i < npoints; i++)
        if (is_border[i])
            pre_force[j++] = pointset[i];
        else
            clone_set[k++] = pointset[i];
    discr = oydiscr_int(pointset, npoints, dim, npoints,
                        pre_force, 0, 0,
                        lowerleft, upperright, &g); // Careful, this is called npoints like in original code but is equal to kpoints
    for (i = 0; i < dim; i++)
        free(pre_force[i]);
    free(pre_force);
    free(clone_set);
    // fprintf(stderr,"%lf\n",discr);
    return discr;
}
