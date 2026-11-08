/*
 * Regression test for cminpack issue #84:
 *   "USE_LAPACK qrfac reads out of bounds when m < n".
 *
 * With fewer rows than columns, R has only min(m,n) diagonal entries, but the
 * USE_LAPACK path of qrfac copied n "diagonal" entries a[i*lda+i] into rdiag,
 * reading past the end of a. Here a is followed by NaN guard values: any read
 * past its end shows up as a NaN in rdiag. rdiag[j] for j >= min(m,n) must be
 * 0, not an off-diagonal element of R.
 *
 * Only registered for USE_LAPACK builds (the non-LAPACK path and the f2c
 * version are fine). Self-checking: exits non-zero on failure.
 */
#include <stdio.h>
#include <math.h>
#include <cminpack.h>
#define real __cminpack_real__

#define MAXSIZE 16
#define GUARD 16

static int run(int m, int n, int pivot)
{
    real buf[MAXSIZE + GUARD];
    real *a = buf;
    real rdiag[MAXSIZE], acnorm[MAXSIZE], wa[MAXSIZE];
    int ipvt[MAXSIZE];
    int i, j, ok = 1;
    const int lda = m;
    const int ltau = m < n ? m : n;

    for (i = 0; i < MAXSIZE + GUARD; ++i) {
        buf[i] = (real)NAN;
    }
    for (j = 0; j < n; ++j) {
        for (i = 0; i < m; ++i) {
            a[i + j * lda] = (real)(1 + i + 2 * j + (i == j ? 3 : 0));
        }
    }
    for (j = 0; j < n; ++j) {
        rdiag[j] = acnorm[j] = wa[j] = 0;
    }

    __cminpack_func__(qrfac)(m, n, a, lda, pivot, ipvt, n, rdiag, acnorm, wa);

    printf("m=%d n=%d pivot=%d: rdiag =", m, n, pivot);
    for (j = 0; j < n; ++j) {
        printf(" %g", (double)rdiag[j]);
        if (!isfinite((double)rdiag[j])) {
            ok = 0;
        }
    }
    for (j = 0; j < ltau; ++j) {
        if (rdiag[j] == 0) {
            ok = 0;
        }
    }
    for (j = ltau; j < n; ++j) {
        if (rdiag[j] != 0) {
            ok = 0; /* R has no diagonal entry there: rdiag is zeroed */
        }
    }
    for (i = m * n; i < MAXSIZE + GUARD; ++i) {
        if (!isnan((double)buf[i])) {
            ok = 0; /* guard overwritten */
        }
    }
    printf(" %s\n", ok ? "OK" : "FAIL");
    return ok ? 0 : 1;
}

int main(void)
{
    int failures = 0;
    int pivot;

    for (pivot = 0; pivot <= 1; ++pivot) {
        failures += run(1, 3, pivot);
        failures += run(2, 3, pivot);
        failures += run(2, 5, pivot);
        failures += run(3, 3, pivot); /* m == n, unaffected */
    }
    if (failures) {
        printf("FAILED: %d check(s)\n", failures);
        return 1;
    }
    printf("all checks passed\n");
    return 0;
}
