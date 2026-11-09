/* Issue #84: with m < n, qrfac must not read past the end of a (USE_LAPACK).
   a is followed by NaN guard values; reading them shows up as a NaN in rdiag. */
#include <stdio.h>
#include <math.h>
#include <cminpack.h>
#define real __cminpack_real__

#define GUARD 16

static int run(int m, int n, int pivot)
{
    real buf[10 + GUARD], rdiag[5], acnorm[5], wa[5];
    int ipvt[5], i, j, ok = 1;

    for (i = 0; i < 10 + GUARD; ++i)
        buf[i] = (real)NAN;
    for (j = 0; j < n; ++j)
        for (i = 0; i < m; ++i)
            buf[i + j * m] = (real)(1 + i + 2 * j + (i == j ? 3 : 0));

    __cminpack_func__(qrfac)(m, n, buf, m, pivot, ipvt, n, rdiag, acnorm, wa);

    printf("m=%d n=%d pivot=%d: rdiag =", m, n, pivot);
    for (j = 0; j < n; ++j) {
        printf(" %g", (double)rdiag[j]);
        ok &= isfinite((double)rdiag[j]) != 0;
    }
    for (i = m * n; i < 10 + GUARD; ++i)
        ok &= isnan((double)buf[i]) != 0; /* guard untouched */
    printf(" %s\n", ok ? "OK" : "FAIL");
    return !ok;
}

int main(void)
{
    int pivot, failures = 0;

    for (pivot = 0; pivot <= 1; ++pivot) {
        failures += run(1, 3, pivot);
        failures += run(2, 5, pivot);
    }
    return failures != 0;
}
