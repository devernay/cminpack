/* Issue #87: the work-array size checks must reject huge m/n without overflow,
   and a negative lwa. */
#include <stdio.h>
#include <limits.h>
#include <cminpack.h>
#define real __cminpack_real__
#define NWA 64

static int ncalls = 0, fails = 0;

static int fcn_nn(void *p, int n, const real *x, real *fvec, int iflag)
{ (void)p; (void)n; (void)x; (void)fvec; (void)iflag; ++ncalls; return -1; }
static int fcnder_mn(void *p, int m, int n, const real *x, real *fvec,
                     real *fjac, int ldfjac, int iflag)
{ (void)p; (void)m; (void)n; (void)x; (void)fvec; (void)fjac; (void)ldfjac;
  (void)iflag; ++ncalls; return -1; }
static int fcnderstr_mn(void *p, int m, int n, const real *x, real *fvec,
                        real *fjrow, int iflag)
{ (void)p; (void)m; (void)n; (void)x; (void)fvec; (void)fjrow; (void)iflag;
  ++ncalls; return -1; }

/* lwa is too small: expect info = 0 and no call of fcn */
static void check(const char *name, int m, int n, int lwa, int info)
{
    int ok = (info == 0 && ncalls == 0);
    printf("%-6s m=%d n=%d lwa=%d: info=%d, %d call(s) of fcn %s\n",
           name, m, n, lwa, info, ncalls, ok ? "OK" : "FAIL");
    fflush(stdout);
    fails += !ok;
    ncalls = 0;
}

int main(void)
{
    real x[NWA], fvec[NWA], fjac[NWA], wa[NWA];
    int ipvt[NWA];
    int m = INT_MAX - 1, n = 1;    /* n*5 + m overflows */

    check("lmder1", m, n, NWA, __cminpack_func__(lmder1)(fcnder_mn, NULL, m, n,
          x, fvec, fjac, m, 1e-8, ipvt, wa, NWA));
    check("lmstr1", m, n, NWA, __cminpack_func__(lmstr1)(fcnderstr_mn, NULL, m,
          n, x, fvec, fjac, n, 1e-8, ipvt, wa, NWA));
    n = 2;                         /* a negative lwa must be rejected too */
    check("hybrd1", n, n, -1, __cminpack_func__(hybrd1)(fcn_nn, NULL, n, x,
          fvec, 1e-8, wa, -1));
    n = 1431655761;                /* n*3 + 13 wraps to 0 */
    check("hybrd1", n, n, NWA, __cminpack_func__(hybrd1)(fcn_nn, NULL, n, x,
          fvec, 1e-8, wa, NWA));
    return fails != 0;
}
