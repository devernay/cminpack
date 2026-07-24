/*      driver for enorm example. */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <minpack.h>
#define real __minpack_real__

#ifdef __cminpack_double__
#define N1 3
#define N2 367
#define SQRTFAC 1.5
#define FAC 10.
#endif
#ifdef __cminpack_long_double__
#define N1 3
#define N2 367
#define SQRTFAC 1.5L
#define FAC 10.L
#define sqrt(x) sqrtl(x)
#endif
#ifdef __cminpack_float__
#define N1 3
#define N2 367
#define SQRTFAC 1.5
#define FAC 10.
#endif
#ifdef __cminpack_half__
#define N1 2
#define N2 2
#define SQRTFAC 1.
#define FAC 2.
#endif

int main()
{
  int i;
  int n;
  real* x;
  real norm;
  real agiant;
  /* dpmpar uses the FORTRAN/f2c calling convention (__minpack_func__), which
     takes the index by pointer, so pass the address of an int -- not a literal. */
  int i2 = 2, i3 = 3;

  real rdwarf = sqrt(__minpack_func__(dpmpar)(&i2)*SQRTFAC) * FAC;
  real rgiant = sqrt(__minpack_func__(dpmpar)(&i3)) / FAC;
#ifdef __cminpack_long_double__
  printf ("dpmpar(2) = %18.7Lg, dpmpar(3) = %18.7Lg\n", __minpack_func__(dpmpar)(&i2), __minpack_func__(dpmpar)(&i3));
  printf ("rdwarf = %.19Lg, rgiant = %.19Lg\n", rdwarf, rgiant);
#else
  printf ("dpmpar(2) = %15.7g, dpmpar(3) = %15.7g\n", (double)__minpack_func__(dpmpar)(&i2), (double)__minpack_func__(dpmpar)(&i3));
  printf ("rdwarf = %.16g, rgiant = %.16g\n", (double)rdwarf, (double)rgiant);
#endif

  n = N1*N1*N2*N2*2+2;
  x = (real*)malloc(n*sizeof(real));
  for (i = 0; i < n; ++i) {
    if (i < 2) {
      x[i] = rdwarf * N1;
    } else {
      x[i] = rdwarf / N2;
    }
  }

  norm = 0;
  for (i = 0; i < n; ++i) {
    norm += x[i]*x[i];
  }
  norm = sqrt(norm) / rdwarf / N1 / 2;
#ifdef __cminpack_long_double__
  printf( "norm/rdwarf (naive) = %.18Lg\n", norm);
#else
  printf( "norm/rdwarf (naive) = %.15g\n", (double)norm);
#endif

  norm = __minpack_func__(enorm)(&n, x);
  norm = norm / rdwarf / N1 / 2;
#ifdef __cminpack_long_double__
  printf( "norm/rdwarf (enorm) = %.18Lg\n", norm);
#else
  printf( "norm/rdwarf (enorm) = %.15g\n", (double)norm);
#endif

  agiant = rgiant / (real)n;
  for (i = 0; i < n; ++i) {
    if (i < 2) {
      x[i] = agiant * N1;
    } else {
      x[i] = agiant / N2;
    }
  }

  norm = 0;
  for (i = 0; i < n; ++i) {
    norm += x[i]*x[i];
  }
  norm = sqrt(norm) / agiant / N1 / 2;
#ifdef __cminpack_long_double__
  printf( "norm/agiant (naive) = %.18Lg\n", norm);
#else
  printf( "norm/agiant (naive) = %.15g\n", (double)norm);
#endif

  norm = __minpack_func__(enorm)(&n, x);
  norm = norm / agiant / N1 / 2;
#ifdef __cminpack_long_double__
  printf( "norm/agiant (enorm) = %.18Lg\n", norm);
#else
  printf( "norm/agiant (enorm) = %.15g\n", (double)norm);
#endif

  free (x);

  return 0;
}
