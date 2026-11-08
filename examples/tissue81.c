/* Issue #81: a user abort during the lmstr Jacobian must be returned as info. */
#include <stdio.h>
#include <cminpack.h>
#define real __cminpack_real__

#define M 4
#define N 2
#define ABORT_VALUE (-7)

/* Linear model f_i(x) = x_0 + t_i*x_1 - y_i. Returns ABORT_VALUE when called
   with iflag == abort_iflag. */
static int fcn(void *p, int m, int n, const real *x, real *fvec, real *fjrow,
               int iflag)
{
    static const real t[M] = {0., 1., 2., 3.};
    static const real y[M] = {1., 3., 4., 7.};
    int abort_iflag = *(const int *)p;
    int i;
    (void)n;
    if (iflag == abort_iflag) {
        return ABORT_VALUE;
    }
    if (iflag == 1) {
        for (i = 0; i < m; ++i) {
            fvec[i] = x[0] + t[i] * x[1] - y[i];
        }
    } else if (iflag >= 2) {
        fjrow[0] = 1.;
        fjrow[1] = t[iflag - 2];
    }
    return 0;
}

static int check(const char *name, int abort_iflag, int info)
{
    int ok = (info == ABORT_VALUE);
    printf("%s, abort at iflag=%d: info=%d (expected %d) %s\n", name,
           abort_iflag, info, ABORT_VALUE, ok ? "OK" : "FAIL");
    return ok ? 0 : 1;
}

int main(void)
{
    /* abort during the function evaluation (iflag=1, already correct), and
       during the first, second and last Jacobian row (iflag=2, 3, M+1). */
    static const int aborts[] = {1, 2, 3, M + 1};
    int k, failures = 0;

    for (k = 0; k < (int)(sizeof aborts / sizeof aborts[0]); ++k) {
        int abort_iflag = aborts[k];
        int ipvt[N];
        real x[N], fvec[M], fjac[N * N], wa[5 * N + M];
        real diag[N], qtf[N], wa1[N], wa2[N], wa3[N], wa4[M];
        int info, nfev, njev;
        real tol = 1e-8;

        x[0] = 0.;
        x[1] = 0.;
        info = __cminpack_func__(lmstr1)(fcn, &abort_iflag, M, N, x, fvec,
                                         fjac, N, tol, ipvt, wa, 5 * N + M);
        failures += check("lmstr1", abort_iflag, info);

        x[0] = 0.;
        x[1] = 0.;
        info = __cminpack_func__(lmstr)(fcn, &abort_iflag, M, N, x, fvec,
                                        fjac, N, tol, tol, 0., 200, diag, 1,
                                        100., 0, &nfev, &njev, ipvt, qtf,
                                        wa1, wa2, wa3, wa4);
        failures += check("lmstr ", abort_iflag, info);
    }

    if (failures) {
        printf("FAILED: %d check(s)\n", failures);
        return 1;
    }
    printf("all checks passed\n");
    return 0;
}
