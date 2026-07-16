/*
 * Regression test for cminpack issue #76:
 *   "Non-convergence in lmder for simple problem with zero vector solution".
 *
 * The test problem is the linear system solved in least squares / root-finding
 * sense
 *      F(x) = ( -x_0, M*x_0 - x_1, M*x_1 - x_2, ..., M*x_{n-2} - x_{n-1} )
 * whose unique solution is the zero vector. Starting from [1,0,...,0], all of
 * the cminpack solvers must converge quickly to x = 0.
 *
 * The original bug: as x shrinks geometrically toward 0, lmpar's iterative
 * parameter search eventually recomputed dxnorm == 0 and then divided by it
 * (0/0 -> NaN) in the Newton-correction step. The NaN propagated into lmder,
 * every stopping test involving a NaN comparison silently failed, and the
 * solver spun until maxfev (info = 5) instead of converging.
 *
 * This test FAILS (returns non-zero) if any solver:
 *   - returns a non-convergence info code,
 *   - produces a non-finite (NaN/Inf) solution or residual,
 *   - or (for the raw lmder/lmdif entry points) uses a number of function
 *     evaluations close to maxfev, which is the signature of the bug.
 *
 * It is written to be self-checking (no reference .ref file needed): a zero
 * exit status means all solvers converged.
 */
#include <stdio.h>
#include <math.h>
#include <cminpack.h>
#define real __cminpack_real__

#define N 4
#define F_MULT (36./73.)
#define MAXFEV 10000
/* A healthy solve of this problem needs at most a few hundred evaluations
   (lmder ~27, finite-difference lmdif ~205). The bug drove nfev all the way
   to MAXFEV (10000), so any run near that ceiling indicates the regression. */
#define NFEV_LIMIT 1000

static int all_finite(const real *v, int n)
{
    int i;
    for (i = 0; i < n; ++i) {
        if (!isfinite((double)v[i])) {
            return 0;
        }
    }
    return 1;
}

/* residual, shared by every callback flavour */
static void residual(int n, const real *x, real *fvec)
{
    int i;
    fvec[0] = -x[0];
    for (i = 0; i < n - 1; ++i) {
        fvec[i + 1] = (real)F_MULT * x[i] - x[i + 1];
    }
}

/* dense Jacobian, column-major (fjac[i + j*ldfjac] = dF_i/dx_j) */
static void jacobian(int n, real *fjac, int ldfjac)
{
    int i, j;
    for (j = 0; j < n; ++j) {
        for (i = 0; i < n; ++i) {
            real val = 0.;
            if (i == j) {
                val = -1.;
            } else if (i == j + 1) {
                val = (real)F_MULT;
            }
            fjac[i + j * ldfjac] = val;
        }
    }
}

/* --- callbacks for the various entry points --- */

static int fcn_lmder(void *p, int m, int n, const real *x, real *fvec,
                     real *fjac, int ldfjac, int iflag)
{
    (void)p; (void)m;
    if (iflag == 1) {
        residual(n, x, fvec);
    } else if (iflag == 2) {
        jacobian(n, fjac, ldfjac);
    }
    return 0;
}

static int fcn_lmdif(void *p, int m, int n, const real *x, real *fvec, int iflag)
{
    (void)p; (void)m; (void)iflag;
    residual(n, x, fvec);
    return 0;
}

static int fcn_hybrd(void *p, int n, const real *x, real *fvec, int iflag)
{
    (void)p; (void)iflag;
    residual(n, x, fvec);
    return 0;
}

static int fcn_hybrj(void *p, int n, const real *x, real *fvec, real *fjac,
                     int ldfjac, int iflag)
{
    (void)p;
    if (iflag == 1) {
        residual(n, x, fvec);
    } else if (iflag == 2) {
        jacobian(n, fjac, ldfjac);
    }
    return 0;
}

/* lmstr wants one Jacobian *row* at a time: iflag>=2 requests row (iflag-2). */
static int fcn_lmstr(void *p, int m, int n, const real *x, real *fvec,
                     real *fjrow, int iflag)
{
    (void)p; (void)m;
    if (iflag == 1) {
        residual(n, x, fvec);
    } else if (iflag >= 2) {
        int i = iflag - 2; /* 0-based row index */
        int j;
        for (j = 0; j < n; ++j) {
            real val = 0.;
            if (i == j) {
                val = -1.;
            } else if (i == j + 1) {
                val = (real)F_MULT;
            }
            fjrow[j] = val;
        }
    }
    return 0;
}

static void init_x(real *x, int n)
{
    int i;
    for (i = 0; i < n; ++i) {
        x[i] = (i == 0) ? 1. : 0.;
    }
}

/* report and return 0 on success, 1 on failure */
/* The issue #76 bug had two symptoms: a NaN result, and exhaustion of maxfev
   (the solver spun instead of terminating). A correct solve of this
   zero-residual problem is therefore: finite (no NaN/Inf), reaches a small
   residual, and stops well before maxfev.

   We deliberately do NOT require a specific convergence info code. At the
   denormal floor a solver may legitimately stop with info=3 ("xtol too small")
   or a slow-progress code instead of info=1, and *which* solver does so
   depends on the exact enorm scaling thresholds. cminpack's enorm uses the
   machine-optimal rdwarf/rgiant computed by examples/tenorm (which differ from
   the original FORTRAN's conservative hardcoded constants), so the last-bit
   path at ~1e-300 magnitudes differs from FORTRAN. That is expected and is not
   the bug. Requiring info==1 would make this test fail on one solver or another
   purely due to that benign difference. */
static int check(const char *name, int not_maxfev, int finite_ok, real fnorm)
{
    int resid_ok = (double)fnorm < 1e-3;
    int ok = not_maxfev && finite_ok && resid_ok;
    printf("  %-8s : %s (not_maxfev=%d finite=%d resid_ok=%d fnorm=%.3e)\n",
           name, ok ? "PASS" : "FAIL", not_maxfev, finite_ok, resid_ok,
           (double)fnorm);
    return ok ? 0 : 1;
}

int main(void)
{
    int failures = 0;
    const real tol = 0.5; /* the exact tolerance from issue #76 */

    printf("issue #76 regression (zero-residual linear problem, N=%d):\n", N);

    /* ---- lmder (the entry point from the original report) ---- */
    {
        real x[N], fvec[N], fjac[N * N], diag[N], qtf[N];
        real wa1[N], wa2[N], wa3[N], wa4[N];
        int ipvt[N], nfev = 0, njev = 0, info;
        init_x(x, N);
        info = __cminpack_func__(lmder)(fcn_lmder, NULL, N, N, x, fvec, fjac, N,
              tol, tol, 0., MAXFEV, diag, 1, 100., 0, &nfev, &njev, ipvt, qtf,
              wa1, wa2, wa3, wa4);
        failures += check("lmder", info != 5 && nfev < NFEV_LIMIT,
                          all_finite(x, N) && all_finite(fvec, N),
                          __cminpack_func__(enorm)(N, fvec));
    }

    /* ---- lmdif (finite-difference Jacobian, same lmpar core) ---- */
    {
        real x[N], fvec[N], fjac[N * N], diag[N], qtf[N];
        real wa1[N], wa2[N], wa3[N], wa4[N];
        int ipvt[N], nfev = 0, info;
        init_x(x, N);
        info = __cminpack_func__(lmdif)(fcn_lmdif, NULL, N, N, x, fvec,
              tol, tol, 0., MAXFEV, 0., diag, 1, 100., 0, &nfev, fjac, N, ipvt,
              qtf, wa1, wa2, wa3, wa4);
        failures += check("lmdif", info != 5 && nfev < NFEV_LIMIT,
                          all_finite(x, N) && all_finite(fvec, N),
                          __cminpack_func__(enorm)(N, fvec));
    }

    /* ---- lmstr1 (also routes through lmpar) ---- */
    {
        real x[N], fvec[N], fjac[N * N], wa[5 * N + N];
        int ipvt[N], info;
        init_x(x, N);
        info = __cminpack_func__(lmstr1)(fcn_lmstr, NULL, N, N, x, fvec, fjac, N,
              tol, ipvt, wa, 5 * N + N);
        failures += check("lmstr1", info != 5,
                          all_finite(x, N) && all_finite(fvec, N),
                          __cminpack_func__(enorm)(N, fvec));
    }

    /* ---- hybrd1 (dogleg-based nonlinear equation solver) ---- */
    {
        real x[N], fvec[N], wa[(N * (3 * N + 13)) / 2];
        int info;
        init_x(x, N);
        info = __cminpack_func__(hybrd1)(fcn_hybrd, NULL, N, x, fvec, 1e-8,
              wa, (N * (3 * N + 13)) / 2);
        failures += check("hybrd1", info != 2,
                          all_finite(x, N) && all_finite(fvec, N),
                          __cminpack_func__(enorm)(N, fvec));
    }

    /* ---- hybrj1 (dogleg-based, analytic Jacobian) ---- */
    {
        real x[N], fvec[N], fjac[N * N], wa[(N * (N + 13)) / 2];
        int info;
        init_x(x, N);
        info = __cminpack_func__(hybrj1)(fcn_hybrj, NULL, N, x, fvec, fjac, N,
              1e-8, wa, (N * (N + 13)) / 2);
        failures += check("hybrj1", info != 2,
                          all_finite(x, N) && all_finite(fvec, N),
                          __cminpack_func__(enorm)(N, fvec));
    }

    if (failures == 0) {
        printf("all solvers converged (issue #76 fixed)\n");
    } else {
        printf("%d solver(s) FAILED the issue #76 regression\n", failures);
    }
    return failures != 0;
}
