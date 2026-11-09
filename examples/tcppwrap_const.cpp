/*
 * Compile-only test for cminpackcpp.hpp (issue #93): a const functor whose
 * operator() is not const must be rejected, like calling a non-const member
 * function on a const object. CTest expects this file to FAIL to compile.
 */
#include <cminpackcpp.hpp>

typedef __cminpack_real__ real;

struct Residual {
    int calls;
    Residual() : calls(0) {}
    int operator()(int m, int n, const real *x, real *fvec, int iflag)
    {
        (void)n; (void)iflag;
        ++calls;
        for (int i = 0; i < m; ++i) {
            fvec[i] = x[0] - (real)i;
        }
        return 0;
    }
};

int run_const_functor()
{
    const Residual fcn;
    real x[1] = { 0. }, fvec[3], wa[3 * 1 + 5 * 1 + 3];
    int iwa[1];
    return cminpack::lmdif1(fcn, 3, 1, x, fvec, (real)1e-6, iwa, wa,
                            3 * 1 + 5 * 1 + 3);
}
