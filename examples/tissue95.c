/* Issue #95: the off-diagonal entries of the problem-10 Jacobian in vecjac.c
   and errjac.c must match a central finite difference of vecfcn.c. */
#include <stdio.h>
#include <math.h>
#include <cminpack.h>
#include "vec.h"
#define N 5

int main(void)
{
    __cminpack_real__ x[N], xp[N], fp[N], fm[N], fjv[N * N], fje[N * N], fd;
    const double h = 1e-5;
    int i, j, nbad = 0;
    for (j = 0; j < N; ++j) x[j] = 0.1 * (j + 1) - 0.3;
    vecjac(N, x, fjv, N, 10);
    errjac(N, x, fje, N, 10);
    for (j = 0; j < N; ++j) {
        for (i = 0; i < N; ++i) xp[i] = x[i];
        xp[j] = x[j] + h; vecfcn(N, xp, fp, 10);
        xp[j] = x[j] - h; vecfcn(N, xp, fm, 10);
        for (i = 0; i < N; ++i) {
            fd = (fp[i] - fm[i]) / (2 * h);
            if (i != j && (fabs(fjv[i + j * N] - fd) > 1e-6 || fabs(fje[i + j * N] - fd) > 1e-6)) {
                printf("fjac(%d,%d): vecjac %g, errjac %g, finite difference %g\n",
                       i + 1, j + 1, fjv[i + j * N], fje[i + j * N], fd);
                ++nbad;
            }
        }
    }
    return nbad != 0;
}
