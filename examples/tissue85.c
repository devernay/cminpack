/* Issue #85: with pivot=0, qrfac documents that wa may coincide with rdiag;
   the result must be the same as with separate arrays. */
#include <stdio.h>
#include <string.h>
#include <cminpack.h>
#define real __cminpack_real__

#define M 4
#define N 3

static const real a0[M * N] = {2., 1., -1., 3., 1., 4., 2., -2., 0., -3., 5., 1.};

int main(void)
{
    real a1[M * N], rdiag1[N], acnorm1[N], wa1[N];
    real a2[M * N], rdiag2[N], acnorm2[N];
    int ipvt[N], i, failures = 0;

    memcpy(a1, a0, sizeof a0);
    __cminpack_func__(qrfac)(M, N, a1, M, 0, ipvt, N, rdiag1, acnorm1, wa1);
    memcpy(a2, a0, sizeof a0);
    __cminpack_func__(qrfac)(M, N, a2, M, 0, ipvt, N, rdiag2, acnorm2, rdiag2);

    for (i = 0; i < M * N; ++i) {
        printf("a[%2d]: separate %13.6g aliased %13.6g\n", i,
               (double)a1[i], (double)a2[i]);
        failures += a1[i] != a2[i];
    }
    for (i = 0; i < N; ++i) {
        printf("rdiag[%d]: separate %13.6g aliased %13.6g\n", i,
               (double)rdiag1[i], (double)rdiag2[i]);
        failures += rdiag1[i] != rdiag2[i] || acnorm1[i] != acnorm2[i];
    }
    printf("%s\n", failures ? "FAIL" : "OK");
    return failures != 0;
}
