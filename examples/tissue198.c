/*
 * Regression test for cminpack issue #198: dpmpar(1) must be the machine
 * precision, 2^(1-p) for a p-bit mantissa. For the IBM double-double long
 * double of Linux ppc64le (p = 106), LDBL_EPSILON is 4.9e-324, and dpmpar(1)
 * returned it: the finite-difference steps were then wrong.
 */
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <cminpack.h>
#define real __cminpack_real__

int main(void)
{
#if defined(__cminpack_long_double__)
  const int p = LDBL_MANT_DIG;
#elif defined(__cminpack_float__)
  const int p = FLT_MANT_DIG;
#else
  const int p = DBL_MANT_DIG;
#endif
  real epsmch = __cminpack_func__(dpmpar)(1);
  real expected = (real)ldexp(1., 1 - p);

  printf("dpmpar(1) = %g, expected %g\n", (double)epsmch, (double)expected);
  return epsmch == expected ? 0 : 1;
}
