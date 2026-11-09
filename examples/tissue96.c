/* Issue #96: the deliberately wrong problem-6 (Watson) Jacobian of errjac.c
   must match errjac.f. Expected values from errjac.f (gfortran -O0). */
#include <stdio.h>
#include <math.h>
#include <cminpack.h>
#include "vec.h"
#define N 4

int main(void)
{
    static const double expected[2] = {-63.082997091728032, -30.014371415760280};
    __cminpack_real__ x[N], fjac[N * N];
    int i, nbad = 0;
    for (i = 0; i < N; ++i) x[i] = 0.1 * (i + 1) - 0.25;
    errjac(N, x, fjac, N, 6);
    for (i = 0; i < 2; ++i) {
        if (fabs(fjac[i] - expected[i]) > 1e-9 * fabs(expected[i])) {
            printf("errjac(%d,1) = %.10g, expected %.10g (errjac.f)\n",
                   i + 1, (double)fjac[i], expected[i]);
            ++nbad;
        }
    }
    return nbad != 0;
}
