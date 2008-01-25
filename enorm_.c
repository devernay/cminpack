/* enorm.f -- translated by f2c (version 20020621).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <math.h>
#include "minpack.h"

double enorm_(const int *n, const double *x)
{
    /* Initialized data */

#define rdwarf 3.834e-20
#define rgiant 1.304e19

    /* System generated locals */
    int i__1;
    double ret_val, d__1;

    /* Local variables */
    int i__;
    double s1, s2, s3, xabs, x1max, x3max, agiant, floatn;

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
    floatn = (double) (*n);
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
	s1 = 1. + s1 * (d__1 * d__1);
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
	s3 = 1. + s3 * (d__1 * d__1);
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
	ret_val = sqrt(s2 * (1. + x3max / s2 * (x3max * s3)));
    }
    if (s2 < x3max) {
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

