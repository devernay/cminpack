/* Issue #82: the long-double enorm must scale whatever the long double format.
   Built with -D__cminpack_long_double__ directly from the enorm sources. */
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <cminpack.h>
#include <minpack.h>

static int fails = 0;
static void check(const char *name, long double value, long double expected)
{
    int ok = isfinite((double)value) && value > 0 &&
             fabsl(value - expected) <= 16 * LDBL_EPSILON * expected;
    /* %Lg is broken on MinGW */
    printf("%-26s = %g (expected %g) %s\n", name, (double)value,
           (double)expected, ok ? "OK" : "FAIL");
    fails += !ok;
}

int main(void)
{
    static const long double big[2] = {1e200L, 1e200L};
    static const long double small[2] = {1e-200L, 1e-200L};
    const long double sqrt2 = 1.4142135623730950488016887242096981L;
    int two = 2;

    printf("LDBL_MANT_DIG = %d, LDBL_MAX_EXP = %d\n", LDBL_MANT_DIG, LDBL_MAX_EXP);
    check("ldenorm({1e200,1e200})", __cminpack_func__(enorm)(2, big), 1e200L * sqrt2);
    check("ldenorm({1e-200,1e-200})", __cminpack_func__(enorm)(2, small), 1e-200L * sqrt2);
    check("ldenorm_({1e200,1e200})", __minpack_func__(enorm)(&two, big), 1e200L * sqrt2);
    check("ldenorm_({1e-200,1e-200})", __minpack_func__(enorm)(&two, small), 1e-200L * sqrt2);
    return fails != 0;
}
