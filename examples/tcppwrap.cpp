/*
 * Test for the header-only C++ wrapper cminpackcpp.hpp (issue #74).
 *
 * It exercises the wrapper with the kinds of callables the C API cannot take
 * directly: a capturing lambda, a stateful functor, and a std::function. Each
 * solver must converge to the known solution, which proves both that the
 * templated trampolines forward the arguments correctly and that user state is
 * carried without any global variable.
 *
 * Self-checking: exits non-zero if any solver fails to converge or lands away
 * from the expected solution.
 */
#include <cminpackcpp.hpp>

#include <cstdio>
#include <cmath>
#include <functional>
#include <vector>

typedef __cminpack_real__ real;

static int failures = 0;

static void check(const char *name, bool converged, double err)
{
    bool ok = converged && err < 1e-6;
    std::printf("  %-20s : %s (converged=%d max_err=%.3e)\n",
                name, ok ? "PASS" : "FAIL", (int)converged, err);
    if (!ok) {
        ++failures;
    }
}

static double max_abs_diff(const real *a, const double *b, int n)
{
    double e = 0.;
    for (int i = 0; i < n; ++i) {
        double d = std::fabs((double)a[i] - b[i]);
        if (d > e) {
            e = d;
        }
    }
    return e;
}

int main()
{
    /* ---------- least-squares: quadratic fit y = c0 + c1*i + c2*i^2 ----------
       The residual is linear in the parameters, so LM recovers the exact
       coefficients. The data is *captured* by the callable, which is exactly
       what the plain C function-pointer API cannot express. */
    const int m = 10, n = 3;
    const double truec[3] = { 0.5, -1.0, 2.0 };
    std::vector<real> data(m);
    for (int i = 0; i < m; ++i) {
        data[i] = (real)(truec[0] + truec[1] * i + truec[2] * i * i);
    }

    /* lmder1 with a capturing lambda (analytic Jacobian) */
    {
        auto fcn = [&](int m_, int n_, const real *x, real *fvec, real *fjac,
                       int ldfjac, int iflag) -> int {
            (void)n_;
            if (iflag == 1) {
                for (int i = 0; i < m_; ++i) {
                    fvec[i] = (x[0] + x[1] * i + x[2] * (real)i * i) - data[i];
                }
            } else if (iflag == 2) {
                for (int i = 0; i < m_; ++i) {
                    fjac[i + 0 * ldfjac] = 1.;
                    fjac[i + 1 * ldfjac] = (real)i;
                    fjac[i + 2 * ldfjac] = (real)i * i;
                }
            }
            return 0;
        };
        real x[3] = { 0., 0., 0. }, fvec[10], fjac[10 * 3], wa[5 * 3 + 10];
        int ipvt[3];
        int info = cminpack::lmder1(fcn, m, n, x, fvec, fjac, m,
                                    1e-10, ipvt, wa, 5 * 3 + 10);
        check("lmder1 (lambda)", info >= 1 && info <= 4, max_abs_diff(x, truec, n));
    }

    /* lmdif1 with a stateful functor (finite-difference Jacobian) */
    {
        struct Fitter {
            const std::vector<real> &d;
            int calls;
            explicit Fitter(const std::vector<real> &dd) : d(dd), calls(0) {}
            int operator()(int m_, int n_, const real *x, real *fvec, int iflag) {
                (void)n_; (void)iflag;
                ++calls;
                for (int i = 0; i < m_; ++i) {
                    fvec[i] = (x[0] + x[1] * i + x[2] * (real)i * i) - d[i];
                }
                return 0;
            }
        } fitter(data);
        real x[3] = { 0., 0., 0. }, fvec[10], wa[10 * 3 + 5 * 3 + 10];
        int iwa[3];
        int info = cminpack::lmdif1(fitter, m, n, x, fvec, 1e-10,
                                    iwa, wa, 10 * 3 + 5 * 3 + 10);
        check("lmdif1 (functor)", info >= 1 && info <= 4, max_abs_diff(x, truec, n));
        /* prove the functor state was actually used (no globals) */
        if (fitter.calls == 0) {
            std::printf("    (functor was never called!)\n");
            ++failures;
        }
    }

    /* ---------- nonlinear equations, root at (2, 3, 2) ---------- */
    const double rootc[3] = { 2., 3., 2. };

    /* hybrd1 with a std::function (finite-difference Jacobian) */
    {
        std::function<int(int, const real *, real *, int)> fcn =
            [](int n_, const real *x, real *fvec, int iflag) -> int {
                (void)n_; (void)iflag;
                fvec[0] = x[0] * x[0] - 4.;
                fvec[1] = x[0] * x[1] - 6.;
                fvec[2] = x[1] + x[2] - 5.;
                return 0;
            };
        real x[3] = { 1., 1., 1. }, fvec[3], wa[(3 * (3 * 3 + 13)) / 2];
        int info = cminpack::hybrd1(fcn, 3, x, fvec, 1e-10,
                                    wa, (3 * (3 * 3 + 13)) / 2);
        check("hybrd1 (std::function)", info == 1, max_abs_diff(x, rootc, 3));
    }

    /* hybrj1 with a capturing lambda (analytic Jacobian) */
    {
        int njev = 0;
        auto fcn = [&](int n_, const real *x, real *fvec, real *fjac,
                       int ldfjac, int iflag) -> int {
            (void)n_;
            if (iflag == 1) {
                fvec[0] = x[0] * x[0] - 4.;
                fvec[1] = x[0] * x[1] - 6.;
                fvec[2] = x[1] + x[2] - 5.;
            } else if (iflag == 2) {
                ++njev;
                for (int j = 0; j < 3; ++j)
                    for (int i = 0; i < 3; ++i)
                        fjac[i + j * ldfjac] = 0.;
                fjac[0 + 0 * ldfjac] = 2. * x[0];
                fjac[1 + 0 * ldfjac] = x[1];
                fjac[1 + 1 * ldfjac] = x[0];
                fjac[2 + 1 * ldfjac] = 1.;
                fjac[2 + 2 * ldfjac] = 1.;
            }
            return 0;
        };
        real x[3] = { 1., 1., 1. }, fvec[3], fjac[3 * 3], wa[(3 * (3 + 13)) / 2];
        int info = cminpack::hybrj1(fcn, 3, x, fvec, fjac, 3, 1e-10,
                                    wa, (3 * (3 + 13)) / 2);
        check("hybrj1 (lambda)", info == 1 && njev > 0, max_abs_diff(x, rootc, 3));
    }

    if (failures == 0) {
        std::printf("cminpackcpp.hpp: all callable kinds converged\n");
    } else {
        std::printf("cminpackcpp.hpp: %d test(s) FAILED\n", failures);
    }
    return failures != 0;
}
