/* Issue #88: ml or mu = INT_MAX (a valid "dense" value) must not overflow. */
#include <stdio.h>
#include <limits.h>
#include <cminpack.h>
#include <minpack.h>
#define real __cminpack_real__
#define N 3

/* Broyden tridiagonal function */
static void broyden(int n, const real *x, real *f)
{
    int k;
    for (k = 0; k < n; ++k)
        f[k] = (3 - 2 * x[k]) * x[k] + 1 - (k > 0 ? x[k-1] : 0)
               - (k < n - 1 ? 2 * x[k+1] : 0);
}
static int fcn(void *p, int n, const real *x, real *f, int iflag)
{ (void)p; (void)iflag; broyden(n, x, f); return 0; }
static void fcn_(const int *n, const real *x, real *f, int *iflag)
{ (void)iflag; broyden(*n, x, f); }

static int fails = 0;
static void check(const char *what, int ok)
{ printf("%-28s %s\n", what, ok ? "OK" : "FAIL"); fails += !ok; }

/* fdjac1 at x = -1 into a poisoned fjac; 1 if equal to the dense result */
static int fdjac1_ok(int f2c, int ml, int mu)
{
    real x[N] = {-1, -1, -1}, fvec[N], wa1[N], wa2[N], d[N*N], j[N*N], e = 0;
    int i, n = N, iflag = 2, m1 = N - 1, ok = 1;
    broyden(N, x, fvec);
    for (i = 0; i < N*N; ++i) d[i] = j[i] = 12345;
    if (f2c) {
        __minpack_func__(fdjac1)(fcn_, &n, x, fvec, d, &n, &iflag, &m1, &m1, &e, wa1, wa2);
        __minpack_func__(fdjac1)(fcn_, &n, x, fvec, j, &n, &iflag, &ml, &mu, &e, wa1, wa2);
    } else {
        __cminpack_func__(fdjac1)(fcn, NULL, N, x, fvec, d, N, m1, m1, 0., wa1, wa2);
        __cminpack_func__(fdjac1)(fcn, NULL, N, x, fvec, j, N, ml, mu, 0., wa1, wa2);
    }
    for (i = 0; i < N*N; ++i) ok &= (d[i] == j[i]);
    return ok;
}

/* nfev of hybrd from x = -1 */
static int hybrd_nfev(int f2c, int ml, int mu)
{
    real x[N] = {-1, -1, -1}, fvec[N], fjac[N*N], diag[N], r[N*(N+1)/2], qtf[N];
    real wa1[N], wa2[N], wa3[N], wa4[N], xtol = 1e-10, e = 0, factor = 100;
    int n = N, maxfev = 2000, mode = 1, nprint = 0, info, nfev, lr = N*(N+1)/2;
    if (f2c)
        __minpack_func__(hybrd)(fcn_, &n, x, fvec, &xtol, &maxfev, &ml, &mu, &e, diag,
                                &mode, &factor, &nprint, &info, &nfev, fjac, &n, r, &lr,
                                qtf, wa1, wa2, wa3, wa4);
    else
        __cminpack_func__(hybrd)(fcn, NULL, N, x, fvec, xtol, maxfev, ml, mu, 0., diag,
                                 1, 100., 0, &nfev, fjac, N, r, lr, qtf, wa1, wa2, wa3, wa4);
    return nfev;
}

int main(void)
{
    int f2c, dense, big;
    for (f2c = 0; f2c <= 1; ++f2c) {
        printf("%s API:\n", f2c ? "f2c" : "C");
        check("fdjac1 ml=mu=INT_MAX", fdjac1_ok(f2c, INT_MAX, INT_MAX));
        check("fdjac1 ml=INT_MAX mu=0", fdjac1_ok(f2c, INT_MAX, 0));
        dense = hybrd_nfev(f2c, N - 1, N - 1);
        big = hybrd_nfev(f2c, INT_MAX, 0);
        printf("hybrd nfev: ml=mu=n-1 %d, ml=INT_MAX mu=0 %d\n", dense, big);
        check("hybrd ml=INT_MAX mu=0", dense == big);
    }
    return fails != 0;
}
