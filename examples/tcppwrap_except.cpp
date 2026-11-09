/*
 * Test for cminpackcpp.hpp (issue #94): an exception thrown by the callable
 * reaches the caller without unwinding through the C solver (this program is
 * linked with a copy of cminpack built without unwind tables), and the
 * callable is not called again after it threw.
 */
#include <cminpackcpp.hpp>

#include <cstdio>

typedef __cminpack_real__ real;

struct UserError {};

static int failures = 0;

static void check(const char *name, bool ok)
{
    std::printf("%-42s: %s\n", name, ok ? "PASS" : "FAIL");
    if (!ok) {
        ++failures;
    }
}

int main()
{
    const int lwa = (2 * (3 * 2 + 13)) / 2;

    /* hybrd1: the exception reaches the caller */
    {
        auto fcn = [](int, const real *x, real *fvec, int) -> int {
            if (x[0] != 0.) {
                throw UserError();
            }
            fvec[0] = x[0] - 1;
            fvec[1] = x[1] - 1;
            return 0;
        };
        real x[2] = { 0., 0. }, fvec[2], wa[lwa];
        bool caught = false;
        try {
            cminpack::hybrd1(fcn, 2, x, fvec, (real)1e-6, wa, lwa);
        } catch (const UserError &) {
            caught = true;
        }
        check("hybrd1: the exception reaches the caller", caught);
    }

    /* lmdif with nprint = 1: the callable is not called again after throwing */
    {
        int calls = 0;
        auto fcn = [&](int m, int, const real *x, real *fvec, int) -> int {
            if (++calls == 2) {
                throw UserError();
            }
            for (int i = 0; i < m; ++i) {
                fvec[i] = x[0] + x[1] * i - 1;
            }
            return 0;
        };
        real x[2] = { 0., 0. }, fvec[3], diag[2], fjac[3 * 2], qtf[2];
        real wa1[2], wa2[2], wa3[2], wa4[3];
        int ipvt[2], nfev;
        bool caught = false;
        try {
            cminpack::lmdif(fcn, 3, 2, x, fvec, (real)1e-6, (real)1e-6, 0., 200,
                            0., diag, 1, 100., 1, &nfev, fjac, 3, ipvt, qtf,
                            wa1, wa2, wa3, wa4);
        } catch (const UserError &) {
            caught = true;
        }
        check("lmdif: not called again after throwing", caught && calls == 2);
    }

    /* a negative return value still gives info < 0, without throwing */
    {
        auto fcn = [](int, const real *, real *, int) -> int { return -1; };
        real x[2] = { 0., 0. }, fvec[2], wa[lwa];
        int info = 0;
        try {
            info = cminpack::hybrd1(fcn, 2, x, fvec, (real)1e-6, wa, lwa);
        } catch (...) {
            info = 1;
        }
        check("negative return: info < 0, no exception", info < 0);
    }

    return failures != 0;
}
