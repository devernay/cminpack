/* Issue #89: fdjac1/fdjac2 must restore x when the user function aborts. */
#include <stdio.h>
#include <cminpack.h>
#include <minpack.h>
#define real __cminpack_real__
#define N 3

/* abort on every call: fdjac1/fdjac2 only call fcn with a perturbed x */
static int fcn_nn(void *p, int n, const real *x, real *f, int iflag)
{ (void)p; (void)n; (void)x; (void)f; (void)iflag; return -1; }
static void fcn_nn_(const int *n, const real *x, real *f, int *iflag)
{ (void)n; (void)x; (void)f; *iflag = -1; }
static int fcn_mn(void *p, int m, int n, const real *x, real *f, int iflag)
{ (void)p; (void)m; (void)n; (void)x; (void)f; (void)iflag; return -1; }
static void fcn_mn_(const int *m, const int *n, const real *x, real *f, int *iflag)
{ (void)m; (void)n; (void)x; (void)f; *iflag = -1; }

static int fails = 0;
static void check(const char *name, const real *x)
{
    int ok = x[0] == 1 && x[1] == 2 && x[2] == 3;
    printf("%-16s %s (x = %.10g %.10g %.10g)\n", name, ok ? "OK" : "FAIL",
           (double)x[0], (double)x[1], (double)x[2]);
    fails += !ok;
}

int main(void)
{
    real x[N], fvec[N] = {0, 0, 0}, fjac[N*N], wa1[N], wa2[N], e = 0;
    int n = N, iflag = 2, dense = N - 1, band = 0;
#define X0 (x[0] = 1, x[1] = 2, x[2] = 3)

    X0; __cminpack_func__(fdjac1)(fcn_nn, NULL, N, x, fvec, fjac, N, dense, dense,
                                  0., wa1, wa2);
    check("fdjac1 dense", x);
    X0; __cminpack_func__(fdjac1)(fcn_nn, NULL, N, x, fvec, fjac, N, band, band,
                                  0., wa1, wa2);
    check("fdjac1 banded", x);
    X0; __cminpack_func__(fdjac2)(fcn_mn, NULL, N, N, x, fvec, fjac, N, 0., wa1);
    check("fdjac2", x);
    X0; __minpack_func__(fdjac1)(fcn_nn_, &n, x, fvec, fjac, &n, &iflag, &dense,
                                 &dense, &e, wa1, wa2);
    check("fdjac1_ dense", x);
    X0; iflag = 2;
    __minpack_func__(fdjac1)(fcn_nn_, &n, x, fvec, fjac, &n, &iflag, &band, &band,
                             &e, wa1, wa2);
    check("fdjac1_ banded", x);
    X0; iflag = 2;
    __minpack_func__(fdjac2)(fcn_mn_, &n, &n, x, fvec, fjac, &n, &iflag, &e, wa1);
    check("fdjac2_", x);
    return fails != 0;
}
