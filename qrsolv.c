/* qrsolv.f -- translated by f2c (version 20020621).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <math.h>
#include "cminpack.h"
#define abs(x) ((x) >= 0 ? (x) : -(x))

/* Subroutine */ void qrsolv(int n, double *r, int ldr, 
	const int *ipvt, const double *diag, const double *qtb, double *x, 
	double *sdiag, double *wa)
{
    /* Initialized data */

#define p5 .5
#define p25 .25

    /* System generated locals */
    int r_dim1, r_offset;

    /* Local variables */
    int i, j, k, l, jp1, kp1;
    double tan, cos, sin, sum, temp, cotan;
    int nsing;
    double qtbpj;

/*     ********** */

/*     subroutine qrsolv */

/*     given an m by n matrix a, an n by n diagonal matrix d, */
/*     and an m-vector b, the problem is to determine an x which */
/*     solves the system */

/*           a*x = b ,     d*x = 0 , */

/*     in the least squares sense. */

/*     this subroutine completes the solution of the problem */
/*     if it is provided with the necessary information from the */
/*     qr factorization, with column pivoting, of a. that is, if */
/*     a*p = q*r, where p is a permutation matrix, q has orthogonal */
/*     columns, and r is an upper triangular matrix with diagonal */
/*     elements of nonincreasing magnitude, then qrsolv expects */
/*     the full upper triangle of r, the permutation matrix p, */
/*     and the first n components of (q transpose)*b. the system */
/*     a*x = b, d*x = 0, is then equivalent to */

/*                  t       t */
/*           r*z = q *b ,  p *d*p*z = 0 , */

/*     where x = p*z. if this system does not have full rank, */
/*     then a least squares solution is obtained. on output qrsolv */
/*     also provides an upper triangular matrix s such that */

/*            t   t               t */
/*           p *(a *a + d*d)*p = s *s . */

/*     s is computed within qrsolv and may be of separate interest. */

/*     the subroutine statement is */

/*       subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa) */

/*     where */

/*       n is a positive integer input variable set to the order of r. */

/*       r is an n by n array. on input the full upper triangle */
/*         must contain the full upper triangle of the matrix r. */
/*         on output the full upper triangle is unaltered, and the */
/*         strict lower triangle contains the strict upper triangle */
/*         (transposed) of the upper triangular matrix s. */

/*       ldr is a positive integer input variable not less than n */
/*         which specifies the leading dimension of the array r. */

/*       ipvt is an integer input array of length n which defines the */
/*         permutation matrix p such that a*p = q*r. column j of p */
/*         is column ipvt(j) of the identity matrix. */

/*       diag is an input array of length n which must contain the */
/*         diagonal elements of the matrix d. */

/*       qtb is an input array of length n which must contain the first */
/*         n elements of the vector (q transpose)*b. */

/*       x is an output array of length n which contains the least */
/*         squares solution of the system a*x = b, d*x = 0. */

/*       sdiag is an output array of length n which contains the */
/*         diagonal elements of the upper triangular matrix s. */

/*       wa is a work array of length n. */

/*     subprograms called */

/*       fortran-supplied ... dabs,dsqrt */

/*     argonne national laboratory. minpack project. march 1980. */
/*     burton s. garbow, kenneth e. hillstrom, jorge j. more */

/*     ********** */
    /* Parameter adjustments */
    --wa;
    --sdiag;
    --x;
    --qtb;
    --diag;
    --ipvt;
    r_dim1 = ldr;
    r_offset = 1 + r_dim1 * 1;
    r -= r_offset;

    /* Function Body */

/*     copy r and (q transpose)*b to preserve input and initialize s. */
/*     in particular, save the diagonal elements of r in x. */

    for (j = 1; j <= n; ++j) {
	for (i = j; i <= n; ++i) {
	    r[i + j * r_dim1] = r[j + i * r_dim1];
	}
	x[j] = r[j + j * r_dim1];
	wa[j] = qtb[j];
    }

/*     eliminate the diagonal matrix d using a givens rotation. */

    for (j = 1; j <= n; ++j) {

/*        prepare the row of d to be eliminated, locating the */
/*        diagonal element using p from the qr factorization. */

	l = ipvt[j];
	if (diag[l] != 0.) {
            for (k = j; k <= n; ++k) {
                sdiag[k] = 0.;
            }
            sdiag[j] = diag[l];

/*        the transformations to eliminate the row of d */
/*        modify only a single element of (q transpose)*b */
/*        beyond the first n, which is initially zero. */

            qtbpj = 0.;
            for (k = j; k <= n; ++k) {

/*           determine a givens rotation which eliminates the */
/*           appropriate element in the current row of d. */

                if (sdiag[k] != 0.) {
                    if (fabs(r[k + k * r_dim1]) < fabs(sdiag[k])) {
                        cotan = r[k + k * r_dim1] / sdiag[k];
                        sin = p5 / sqrt(p25 + p25 * (cotan * cotan));
                        cos = sin * cotan;
                    } else {
                        tan = sdiag[k] / r[k + k * r_dim1];
                        cos = p5 / sqrt(p25 + p25 * (tan * tan));
                        sin = cos * tan;
                    }

/*           compute the modified diagonal element of r and */
/*           the modified element of ((q transpose)*b,0). */

                    r[k + k * r_dim1] = cos * r[k + k * r_dim1] + sin * sdiag[k];
                    temp = cos * wa[k] + sin * qtbpj;
                    qtbpj = -sin * wa[k] + cos * qtbpj;
                    wa[k] = temp;

/*           accumulate the tranformation in the row of s. */

                    kp1 = k + 1;
                    if (n >= kp1) {
                        for (i = kp1; i <= n; ++i) {
                            temp = cos * r[i + k * r_dim1] + sin * sdiag[i];
                            sdiag[i] = -sin * r[i + k * r_dim1] + cos * sdiag[i];
                            r[i + k * r_dim1] = temp;
                        }
                    }
                }
            }
        }

/*        store the diagonal element of s and restore */
/*        the corresponding diagonal element of r. */

	sdiag[j] = r[j + j * r_dim1];
	r[j + j * r_dim1] = x[j];
    }

/*     solve the triangular system for z. if the system is */
/*     singular, then obtain a least squares solution. */

    nsing = n;
    for (j = 1; j <= n; ++j) {
	if (sdiag[j] == 0. && nsing == n) {
	    nsing = j - 1;
	}
	if (nsing < n) {
	    wa[j] = 0.;
	}
    }
    if (nsing >= 1) {
        for (k = 1; k <= nsing; ++k) {
            j = nsing - k + 1;
            sum = 0.;
            jp1 = j + 1;
            if (nsing >= jp1) {
                for (i = jp1; i <= nsing; ++i) {
                    sum += r[i + j * r_dim1] * wa[i];
                }
            }
            wa[j] = (wa[j] - sum) / sdiag[j];
        }
    }

/*     permute the components of z back to components of x. */

    for (j = 1; j <= n; ++j) {
	l = ipvt[j];
	x[l] = wa[j];
    }
    return;

/*     last card of subroutine qrsolv. */

} /* qrsolv_ */

