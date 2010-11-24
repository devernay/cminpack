/* covar.f -- translated by f2c (version 20050501).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include <math.h>
#include "minpack.h"

/* Subroutine */ void covar(int n, double *r__, int ldr, 
	const int *ipvt, double tol, double *wa)
{
    /* System generated locals */
    int r_dim1, r_offset;

    /* Local variables */
    int i__, j, k, l, ii, jj, km1;
    int sing;
    double temp, tolr;

/*     ********** */

/*     subroutine covar */

/*     given an m by n matrix a, the problem is to determine */
/*     the covariance matrix corresponding to a, defined as */

/*                    t */
/*           inverse(a *a) . */

/*     this subroutine completes the solution of the problem */
/*     if it is provided with the necessary information from the */
/*     qr factorization, with column pivoting, of a. that is, if */
/*     a*p = q*r, where p is a permutation matrix, q has orthogonal */
/*     columns, and r is an upper triangular matrix with diagonal */
/*     elements of nonincreasing magnitude, then covar expects */
/*     the full upper triangle of r and the permutation matrix p. */
/*     the covariance matrix is then computed as */

/*                      t     t */
/*           p*inverse(r *r)*p  . */

/*     if a is nearly rank deficient, it may be desirable to compute */
/*     the covariance matrix corresponding to the linearly independent */
/*     columns of a. to define the numerical rank of a, covar uses */
/*     the tolerance tol. if l is the largest integer such that */

/*           abs(r(l,l)) .gt. tol*abs(r(1,1)) , */

/*     then covar computes the covariance matrix corresponding to */
/*     the first l columns of r. for k greater than l, column */
/*     and row ipvt(k) of the covariance matrix are set to zero. */

/*     the subroutine statement is */

/*       subroutine covar(n,r,ldr,ipvt,tol,wa) */

/*     where */

/*       n is a positive integer input variable set to the order of r. */

/*       r is an n by n array. on input the full upper triangle must */
/*         contain the full upper triangle of the matrix r. on output */
/*         r contains the square symmetric covariance matrix. */

/*       ldr is a positive integer input variable not less than n */
/*         which specifies the leading dimension of the array r. */

/*       ipvt is an integer input array of length n which defines the */
/*         permutation matrix p such that a*p = q*r. column j of p */
/*         is column ipvt(j) of the identity matrix. */

/*       tol is a nonnegative input variable used to define the */
/*         numerical rank of a in the manner described above. */

/*       wa is a work array of length n. */

/*     subprograms called */

/*       fortran-supplied ... dabs */

/*     argonne national laboratory. minpack project. august 1980. */
/*     burton s. garbow, kenneth e. hillstrom, jorge j. more */

/*     ********** */
    /* Parameter adjustments */
    --wa;
    --ipvt;
    tolr = tol * fabs(r__[0]);
    r_dim1 = ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;

    /* Function Body */

/*     form the inverse of r in the full upper triangle of r. */

    l = 0;
    for (k = 1; k <= n; ++k) {
	if (fabs(r__[k + k * r_dim1]) <= tolr) {
	    goto TERMINATE_INVERSE;
	}
	r__[k + k * r_dim1] = 1. / r__[k + k * r_dim1];
	km1 = k - 1;
	if (km1 >= 1) {
            for (j = 1; j <= km1; ++j) {
                temp = r__[k + k * r_dim1] * r__[j + k * r_dim1];
                r__[j + k * r_dim1] = 0.;
                for (i__ = 1; i__ <= j; ++i__) {
                    r__[i__ + k * r_dim1] -= temp * r__[i__ + j * r_dim1];
                }
            }
        }
	l = k;
    }
TERMINATE_INVERSE:

/*     form the full upper triangle of the inverse of (r transpose)*r */
/*     in the full upper triangle of r. */

    if (l >= 1) {
        for (k = 1; k <= l; ++k) {
            km1 = k - 1;
            if (km1 >= 1) {
                for (j = 1; j <= km1; ++j) {
                    temp = r__[j + k * r_dim1];
                    for (i__ = 1; i__ <= j; ++i__) {
                        r__[i__ + j * r_dim1] += temp * r__[i__ + k * r_dim1];
                    }
                }
            }
            temp = r__[k + k * r_dim1];
            for (i__ = 1; i__ <= k; ++i__) {
                r__[i__ + k * r_dim1] = temp * r__[i__ + k * r_dim1];
            }
        }
    }

/*     form the full lower triangle of the covariance matrix */
/*     in the strict lower triangle of r and in wa. */

    for (j = 1; j <= n; ++j) {
	jj = ipvt[j];
	sing = j > l;
	for (i__ = 1; i__ <= j; ++i__) {
	    if (sing) {
		r__[i__ + j * r_dim1] = 0.;
	    }
	    ii = ipvt[i__];
	    if (ii > jj) {
		r__[ii + jj * r_dim1] = r__[i__ + j * r_dim1];
	    }
	    if (ii < jj) {
		r__[jj + ii * r_dim1] = r__[i__ + j * r_dim1];
	    }
	}
	wa[jj] = r__[j + j * r_dim1];
    }

/*     symmetrize the covariance matrix in r. */

    for (j = 1; j <= n; ++j) {
	for (i__ = 1; i__ <= j; ++i__) {
	    r__[i__ + j * r_dim1] = r__[j + i__ * r_dim1];
	}
	r__[j + j * r_dim1] = wa[j];
    }

/*     last card of subroutine covar. */

} /* covar_ */

