/* enorm.f -- translated by f2c (version 20020621).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "minpack.h"
#include <math.h>

#include "minpackP.h"

/*
  About the values for rdwarf and rgiant.

  enorm accumulates the sum of squares in three buckets so that squaring a
  component can neither overflow nor destructively underflow: components in the
  middle band [rdwarf, rgiant/n] are summed directly as x*x, while smaller and
  larger ones are rescaled first. The only real constraints are that rdwarf**2
  must not underflow and rgiant**2 must not overflow *in the floating-point
  format actually being used*.

  The original MINPACK values, used identically in both the single- and the
  double-precision FORTRAN sources, were:
#define rdwarf 3.834e-20
#define rgiant 1.304e19
  See for example:
    http://www.netlib.org/slatec/src/denorm.f
    http://www.netlib.org/slatec/src/enorm.f
  These give rdwarf**2 ~ 1.5e-39 and rgiant**2 ~ 1.7e38, i.e. they are sized for
  a machine whose dynamic range is about 1e+-38. MINPACK dates from 1980, before
  IEEE 754 (1985). Back then the worst-case "double precision" in wide use --
  notably the DEC VAX D_floating format -- had the *same* exponent range as
  single precision (~1e+-38, with only more mantissa bits). A single conservative
  pair of constants was therefore chosen to be portable "for every known
  computer".

  On IEEE 754 those constants are safe but far too conservative:
    - in IEEE single, rdwarf = 3.834e-20 is below sqrt(FLT_MIN) = 1.084e-19, so
      rdwarf**2 ~ 1.5e-39 UNDERFLOWS, violating the constraint above;
    - in IEEE double (DBL_MIN ~ 2.2e-308, DBL_MAX ~ 1.8e308) the usable middle
      band spans ~1e+-308, but 3.834e-20/1.304e19 restrict it to ~1e+-38 and so
      needlessly rescale many components.

  cminpack therefore tunes rdwarf/rgiant to the actual range of each real type,
  following MPFIT (http://cow.physics.wisc.edu/~craigm/idl/fitting.html):
    rdwarf = sqrt(dpmpar(2)*1.5) * 10      (dpmpar(2) = smallest normal)
    rgiant = sqrt(dpmpar(3)) * 0.1         (dpmpar(3) = largest magnitude)
  Half precision does not work well with that formula, so for half we use:
    rdwarf = sqrt(dpmpar(2)) * 2
    rgiant = sqrt(dpmpar(3)) * 0.5
  (half cminpack is only a proof of concept). The values below are those
  formulas evaluated per type; see examples/tenorm*.c, which prints them.

  Consequence -- why cminpack's double enorm can differ from FORTRAN MINPACK:
  because the middle-band boundary differs, a vector with very small
  (< 3.834e-20) or very large (> 1.304e19) components is bucketed differently
  and its Euclidean norm can then differ in the last bit. Both norms are equally
  correct. cminpack's own pure-C and f2c sources are, however, bit-identical to
  each other at every precision (verified by examples/crosscheck.py).

  These constants are NOT the main reason cminpack and the original FORTRAN take
  different iteration counts on the harder test problems. Aligning them
  (restoring the Argonne values here) does not reproduce the FORTRAN runs: many
  of the difficult More/Garbow/Hillstrom problems (the intensive driver programs
  under examples/) still diverge, and enorm-insensitive drivers such as lmddrv
  diverge by the same amount either way. dpmpar returns identical machine constants on both sides, so the
  convergence tolerance is not involved. The dominant cause is floating-point
  contraction: gcc and gfortran fuse "a*b + c" into a fused multiply-add (FMA) at
  different places, so intermediate values differ by one ULP. Building both
  compilers with -ffp-contract=off removes most of the divergence; the small
  residual is chaotic amplification of last-bit differences over hundreds of
  iterations on ill-conditioned problems, where a single ULP can flip a
  trust-region accept/reject decision and send the two runs to different -- but
  equally valid -- results. Every affected problem still converges. cminpack
  keeps the IEEE-appropriate constants deliberately; see README.md, section
  "Numerical differences from FORTRAN MINPACK".
*/
#define double_dwarf (1.82691291192569e-153)
#define double_giant (1.34078079299426e+153)
#define long_double_dwarf (2.245696932951581572e-2465l)
#define long_double_giant (1.090748135619415929e+2465l)
#define float_dwarf (1.327871072777421e-18f)
#define float_giant (1.844674297419792e+18f)
#define half_dwarf (0.015625f)
#define half_giant (127.9375f)

#define dwarf(type) _dwarf(type)
#define _dwarf(type) type ## _dwarf
#define giant(type) _giant(type)
#define _giant(type) type ## _giant

#define rdwarf dwarf(realm)
#define rgiant giant(realm)

__minpack_attr__
real __minpack_func__(enorm)(const int *n, const real *x)
{
    /* System generated locals */
    int i__1;
    real ret_val, d__1;

    /* Local variables */
    int i__;
    real s1, s2, s3, xabs, x1max, x3max, agiant, floatn;

/*     ********** */

/*     function enorm */

/*     given an n-vector x, this function calculates the */
/*     euclidean norm of x. */

/*     the euclidean norm is computed by accumulating the sum of */
/*     squares in three different sums. the sums of squares for the */
/*     small and large components are scaled so that no overflows */
/*     occur. non-destructive underflows are permitted. underflows */
/*     and overflows do not occur in the computation of the unscaled */
/*     sum of squares for the intermediate components. */
/*     the definitions of small, intermediate and large components */
/*     depend on two constants, rdwarf and rgiant. the main */
/*     restrictions on these constants are that rdwarf**2 not */
/*     underflow and rgiant**2 not overflow. the constants */
/*     given here are suitable for every known computer. */

/*     the function statement is */

/*       double precision function enorm(n,x) */

/*     where */

/*       n is a positive integer input variable. */

/*       x is an input array of length n. */

/*     subprograms called */

/*       fortran-supplied ... dabs,dsqrt */

/*     argonne national laboratory. minpack project. march 1980. */
/*     burton s. garbow, kenneth e. hillstrom, jorge j. more */

/*     ********** */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    s1 = 0.;
    s2 = 0.;
    s3 = 0.;
    x1max = 0.;
    x3max = 0.;
    floatn = (real) (*n);
    agiant = rgiant / floatn;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xabs = fabs(x[i__]);
	if (xabs > rdwarf && xabs < agiant) {
	    goto L70;
	}
	if (xabs <= rdwarf) {
	    goto L30;
	}

/*              sum for large components. */

	if (xabs <= x1max) {
	    goto L10;
	}
/* Computing 2nd power */
	d__1 = x1max / xabs;
	s1 = 1 + s1 * (d__1 * d__1);
	x1max = xabs;
	goto L20;
L10:
/* Computing 2nd power */
	d__1 = xabs / x1max;
	s1 += d__1 * d__1;
L20:
	goto L60;
L30:

/*              sum for small components. */

	if (xabs <= x3max) {
	    goto L40;
	}
/* Computing 2nd power */
	d__1 = x3max / xabs;
	s3 = 1 + s3 * (d__1 * d__1);
	x3max = xabs;
	goto L50;
L40:
	if (xabs != 0.) {
/* Computing 2nd power */
	    d__1 = xabs / x3max;
	    s3 += d__1 * d__1;
	}
L50:
L60:
	goto L80;
L70:

/*           sum for intermediate components. */

/* Computing 2nd power */
	d__1 = xabs;
	s2 += d__1 * d__1;
L80:
/* L90: */
	;
    }

/*     calculation of norm. */

    if (s1 == 0.) {
	goto L100;
    }
    ret_val = x1max * sqrt(s1 + s2 / x1max / x1max);
    goto L130;
L100:
    if (s2 == 0.) {
	goto L110;
    }
    if (s2 >= x3max) {
	ret_val = sqrt(s2 * (1 + x3max / s2 * (x3max * s3)));
    } else {
	ret_val = sqrt(x3max * (s2 / x3max + x3max * s3));
    }
    goto L120;
L110:
    ret_val = x3max * sqrt(s3);
L120:
L130:
    return ret_val;

/*     last card of function enorm. */
} /* enorm_ */

