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
    double discr, maxdiscr, coordlist[dim][nforced + 1]; // actually coordlist[dim][nforced] is enough but nforced sometimes is 0 which upsets sanitizer. Should probably address the root cause, but I don't know how at the moment.
    int indexes[dim];
    int i, j, k, h, status, dimension;
    double biggest[dim][nforced + 1], smallest[dim][nforced + 1];
    double ***big_bord;
    double ***small_bord;
    big_bord = malloc(dim * sizeof(double **));
    small_bord = malloc(dim * sizeof(double **));
    for (i = 0; i < dim; i++)
    {
        big_bord[i] = malloc((nforced + 1) * sizeof(double *));
        small_bord[i] = malloc((nforced + 1) * sizeof(double *));
        for (j = 0; j < (nforced + 1); j++)
        {
            big_bord[i][j] = malloc(dim * sizeof(double));
            small_bord[i][j] = malloc(dim * sizeof(double));
            for (k = 0; k < dim; k++)
            {
                big_bord[i][j][k] = 0.0;
                small_bord[i][j][k] = 0.0;
            }
        }
    }
    int maxpoints[dim];
    // biggest[i][j]: biggest product of coords 0--i for hitting j points
    // smallest[i][j]: smallest product of coords 0--i for hitting j+1 points
    // maxpoints[i]: number of points you get in total from coords 0--i
    double vol1 = 1.0, vol2 = 1.0;
    for (i = 0; i < dim; i++)
    {
        vol1 *= lowerleft[i];
        vol2 *= upperright[i];
    }
    // fprintf(stderr, "Get in cell\n");
    maxdiscr = vol2 - (double)rempoints / npoints;
    discr = (double)(rempoints + nforced) / npoints - vol1;
    if (discr > maxdiscr)
    {
        maxdiscr = discr;
    }
    if (maxdiscr < g->globallower)
    {
        for (i = 0; i < dim; i++)
        {
            for (j = 0; j < (nforced + 1); j++)
            {
                free(small_bord[i][j]);
                free(big_bord[i][j]);
            }
            free(big_bord[i]);
            free(small_bord[i]);
        }
        free(big_bord);
        free(small_bord);
        return maxdiscr;
    }
    // quicker code for use in some easy cells
    // otherwise, get to work...
    for (i = 0; i < dim; i++)
    {
        indexes[i] = 0;
        for (j = 0; j <= nforced; j++)
        {
            smallest[i][j] = 1.0;
            biggest[i][j] = 0.0;
        }
    }
    for (i = 0; i < nforced; i++)
    {
        status = 0;
        for (j = 0; j < dim; j++)
        {
            // order is chosen to handle final box
            if (forced[i][j] <= lowerleft[j])
                continue;
            else if (forced[i][j] >= upperright[j])
                break;
            else
            { // strictly internal
                if (status)
                {
                    fprintf(stderr, "PROBLEM: Point occurs as double-internal\n");
                    fflush(stderr);
                    abort();
                }
                status = 1;
                dimension = j;
            }
        }
        if (j == dim)
        { // else: hit "break", skip
            if (status)
            {
                coordlist[dimension][indexes[dimension]] = forced[i][dimension];
                indexes[dimension] += 1;
            }
            else
            { // completely internal
                rempoints++;
            }
        }
    }

    for (i = 0; i < dim; i++)
        qsort(&(coordlist[i][0]), indexes[i], sizeof(double), cmpdbl);
#pragma GCC diagnostic push                            // save the actual diag context
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized" // disable maybe warnings
    maxpoints[0] = indexes[0];
#pragma GCC diagnostic pop // restore previous diag context
    for (i = 1; i < dim; i++)
        maxpoints[i] = maxpoints[i - 1] + indexes[i];

    // coord 0 first, since there is no recursion for that:
    smallest[0][0] = lowerleft[0];
    small_bord[0][0][0] = lowerleft[0]; // Changed here
    for (j = 0; j < indexes[0] && j < nforced; j++) // actually indexes[0] cannot be larger than nforced, but the compiler is too dumb to understand that (and I want to use -Werror)
    {
        smallest[0][j + 1] = coordlist[0][j];
        biggest[0][j] = coordlist[0][j];
        small_bord[0][j + 1][0] = coordlist[0][j]; // Changed here
        big_bord[0][j][0] = coordlist[0][j];
    }
    biggest[0][indexes[0]] = upperright[0];
    big_bord[0][indexes[0]][0] = upperright[0];

    // INIT CORRECT

    // now, use these to find lower, upper limits
    // mode: always contain "rempoints", additionally
    //         pick from smallest[dim-1], biggest[dim-1]
    maxdiscr = 0;
    // fprintf(stderr, "DynProg time\n");
    for (i = 0; i <= maxpoints[dim - 1]; i++)
    {
        discr = (double)(rempoints + i) / npoints - smallest[dim - 1][i];
    }
    // fprintf(stderr, "Get out cell2\n");
    if (maxdiscr > g->globallower)
    {
        g->globallower = maxdiscr;
    }
    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < (nforced + 1); j++)
        {
            free(small_bord[i][j]);
            free(big_bord[i][j]);
        }
        free(big_bord[i]);
        free(small_bord[i]);
    }
    free(big_bord);
    free(small_bord);

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
