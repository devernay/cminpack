/* rwupdt.f -- translated by f2c (version 20020621).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <math.h>
#include "cminpack.h"
#define abs(x) ((x) >= 0 ? (x) : -(x))

/* Subroutine */ void rwupdt(int n, double *r__, int ldr, 
	const double *w, double *b, double *alpha, double *cos__, 
	double *sin__)
{
    /* Initialized data */

#define p5 .5
#define p25 .25

    /* System generated locals */
    int r_dim1, r_offset;

    /* Local variables */
    int i__, j, jm1;
    double tan__, temp, rowj, cotan;

/*     ********** */

/*     subroutine rwupdt */

/*     given an n by n upper triangular matrix r, this subroutine */
/*     computes the qr decomposition of the matrix formed when a row */
/*     is added to r. if the row is specified by the vector w, then */
/*     rwupdt determines an orthogonal matrix q such that when the */
/*     n+1 by n matrix composed of r augmented by w is premultiplied */
/*     by (q transpose), the resulting matrix is upper trapezoidal. */
/*     the matrix (q transpose) is the product of n transformations */

/*           g(n)*g(n-1)* ... *g(1) */

/*     where g(i) is a givens rotation in the (i,n+1) plane which */
/*     eliminates elements in the (n+1)-st plane. rwupdt also */
/*     computes the product (q transpose)*c where c is the */
/*     (n+1)-vector (b,alpha). q itself is not accumulated, rather */
/*     the information to recover the g rotations is supplied. */

/*     the subroutine statement is */

/*       subroutine rwupdt(n,r,ldr,w,b,alpha,cos,sin) */

/*     where */

/*       n is a positive integer input variable set to the order of r. */

/*       r is an n by n array. on input the upper triangular part of */
/*         r must contain the matrix to be updated. on output r */
/*         contains the updated triangular matrix. */

/*       ldr is a positive integer input variable not less than n */
/*         which specifies the leading dimension of the array r. */

/*       w is an input array of length n which must contain the row */
/*         vector to be added to r. */

/*       b is an array of length n. on input b must contain the */
/*         first n elements of the vector c. on output b contains */
/*         the first n elements of the vector (q transpose)*c. */

/*       alpha is a variable. on input alpha must contain the */
/*         (n+1)-st element of the vector c. on output alpha contains */
/*         the (n+1)-st element of the vector (q transpose)*c. */

/*       cos is an output array of length n which contains the */
/*         cosines of the transforming givens rotations. */

/*       sin is an output array of length n which contains the */
/*         sines of the transforming givens rotations. */

/*     subprograms called */

/*       fortran-supplied ... dabs,dsqrt */

/*     argonne national laboratory. minpack project. march 1980. */
/*     burton s. garbow, dudley v. goetschel, kenneth e. hillstrom, */
/*     jorge j. more */

/*     ********** */
    /* Parameter adjustments */
    --sin__;
    --cos__;
    --b;
    --w;
    r_dim1 = ldr;
    r_offset = 1 + r_dim1 * 1;
    r__ -= r_offset;

    /* Function Body */

    for (j = 1; j <= n; ++j) {
	rowj = w[j];
	jm1 = j - 1;

/*        apply the previous transformations to */
/*        r(i,j), i=1,2,...,j-1, and to w(j). */

	if (jm1 >= 1) {
            for (i__ = 1; i__ <= jm1; ++i__) {
                temp = cos__[i__] * r__[i__ + j * r_dim1] + sin__[i__] * rowj;
                rowj = -sin__[i__] * r__[i__ + j * r_dim1] + cos__[i__] * rowj;
                r__[i__ + j * r_dim1] = temp;
            }
        }

/*        determine a givens rotation which eliminates w(j). */

	cos__[j] = 1.;
	sin__[j] = 0.;
	if (rowj != 0.) {
            if (fabs(r__[j + j * r_dim1]) < fabs(rowj)) {
                cotan = r__[j + j * r_dim1] / rowj;
                sin__[j] = p5 / sqrt(p25 + p25 * (cotan * cotan));
                cos__[j] = sin__[j] * cotan;
            } else {
                tan__ = rowj / r__[j + j * r_dim1];
                cos__[j] = p5 / sqrt(p25 + p25 * (tan__ * tan__));
                sin__[j] = cos__[j] * tan__;
            }

/*        apply the current transformation to r(j,j), b(j), and alpha. */

            r__[j + j * r_dim1] = cos__[j] * r__[j + j * r_dim1] + sin__[j] * rowj;
            temp = cos__[j] * b[j] + sin__[j] * *alpha;
            *alpha = -sin__[j] * b[j] + cos__[j] * *alpha;
            b[j] = temp;
        }
    }

/*     last card of subroutine rwupdt. */

} /* rwupdt_ */

