#include "cminpack.h"
#include <math.h>
#ifdef USE_LAPACK
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#endif
#include "cminpackP.h"

__cminpack_attr__
void __cminpack_func__(qrfac)(int m, int n, real *a, int
	lda, int pivot, int *ipvt, int lipvt, real *rdiag,
	 real *acnorm, real *wa)
{
#ifdef USE_LAPACK
    __CLPK_integer m_ = m;
    __CLPK_integer n_ = n;
    __CLPK_integer lda_ = lda;
    __CLPK_integer *jpvt;

    int i, j, k;
    real t;
    real* tau = wa;
    const __CLPK_integer ltau = m > n ? n : m;
    __CLPK_integer lwork = -1;
    __CLPK_integer info = 0;
    real* work;

    if (pivot) {
        assert( lipvt >= n );
        if (sizeof(__CLPK_integer) != sizeof(ipvt[0])) {
            jpvt = malloc(n*sizeof(__CLPK_integer));
            if (jpvt == NULL) {
                fprintf(stderr, "cminpack qrfac: out of memory\n");
                abort();
            }
        } else {
            /* __CLPK_integer is actually an int, just do a cast */
            jpvt = (__CLPK_integer *)ipvt;
        }
        /* set all columns free. Use sizeof(*jpvt) (== sizeof(__CLPK_integer)),
           not sizeof(int): when __CLPK_integer is wider than int (an ILP64
           LAPACK) jpvt is a separately allocated __CLPK_integer array, and
           sizeof(int)*n would clear only the first int-sized bytes -- i.e. only
           the first few whole elements -- leaving the remaining elements as
           garbage that geqp3 would treat as fixed (leading) columns. */
        memset(jpvt, 0, sizeof(*jpvt)*n);
    }
    
    /* query optimal size of work */
    lwork = -1;
    if (pivot) {
        __cminpack_lapack__(geqp3_)(&m_,&n_,a,&lda_,jpvt,tau,tau,&lwork,&info);
        lwork = (__CLPK_integer)tau[0];
        assert( lwork >= 3*n+1  );
    } else {
        __cminpack_lapack__(geqrf_)(&m_,&n_,a,&lda_,tau,tau,&lwork,&info);
        lwork = (__CLPK_integer)tau[0];
        assert( lwork >= 1 && lwork >= n );
    }
    if (info != 0) {
        fprintf(stderr, "cminpack qrfac: LAPACK QR workspace query failed "
                "(info=%d)\n", (int)info);
        abort();
    }
    
    /* alloc work area.
       NOTE: qrfac() returns void, so it cannot report a LAPACK error or an
       out-of-memory condition to its callers (lmder/lmstr/hybr*). Rather than
       let a NULL work pointer or a nonzero LAPACK info silently crash or
       corrupt a release (NDEBUG) build -- where the asserts are compiled out --
       these unrecoverable conditions are reported to stderr and abort(). A
       non-fatal alternative would require adding an int error return to the
       qrfac API (a breaking change), deliberately not done here. */
    work = (real *)malloc(sizeof(real)*lwork);
    if (work == NULL) {
        fprintf(stderr, "cminpack qrfac: out of memory\n");
        abort();
    }
    
    /* set acnorm first (from the doc of qrfac, acnorm may point to the same area as rdiag) */
    if (acnorm != rdiag) {
        for (j = 0; j < n; ++j) {
            acnorm[j] = __cminpack_func__(enorm)(m, &a[j * lda]);
        }
    }
    
    /* QR decomposition */
    if (pivot) {
        __cminpack_lapack__(geqp3_)(&m_,&n_,a,&lda_,jpvt,tau,work,&lwork,&info);
    } else {
        __cminpack_lapack__(geqrf_)(&m_,&n_,a,&lda_,tau,work,&lwork,&info);
    }
    if (info != 0) {
        fprintf(stderr, "cminpack qrfac: LAPACK QR factorization failed "
                "(info=%d)\n", (int)info);
        abort();
    }
    
    /* set rdiag and rebuild the lower triangle in minpack qrfac's u-vector
       format. tau (wa) may coincide with rdiag: read tau[i] first. */
    for(i=0 ; i<ltau ; ++i) {
        k = i*lda+i;
        t = tau[i];
        rdiag[i] = a[k];
        a[k] = t;
        for(j=i+1 ; j<m ; j++) {
            k++;
            a[k] *= t;
        }
    }
    for(i=ltau ; i<n ; ++i) {
        rdiag[i] = 0.;
    }
    
    free(work);
    if (pivot) {
        /* convert back jpvt to ipvt */
        if (sizeof(__CLPK_integer) != sizeof(ipvt[0])) {
            for(i=0; i<n; ++i) {
                ipvt[i] = jpvt[i];
            }
            free(jpvt);
        }
    }
#else /* !USE_LAPACK */
    /* Initialized data */

#define p05 .05

    /* System generated locals */
    real d1;

    /* Local variables */
    int i, j, k, jp1;
    real sum;
    real temp;
    int minmn;
    real epsmch;
    real ajnorm;

/*     ********** */

/*     subroutine qrfac */

/*     this subroutine uses householder transformations with column */
/*     pivoting (optional) to compute a qr factorization of the */
/*     m by n matrix a. that is, qrfac determines an orthogonal */
/*     matrix q, a permutation matrix p, and an upper trapezoidal */
/*     matrix r with diagonal elements of nonincreasing magnitude, */
/*     such that a*p = q*r. the householder transformation for */
/*     column k, k = 1,2,...,min(m,n), is of the form */

/*                           t */
/*           i - (1/u(k))*u*u */

/*     where u has zeros in the first k-1 positions. the form of */
/*     this transformation and the method of pivoting first */
/*     appeared in the corresponding linpack subroutine. */

/*     the subroutine statement is */

/*       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa) */

/*     where */

/*       m is a positive integer input variable set to the number */
/*         of rows of a. */

/*       n is a positive integer input variable set to the number */
/*         of columns of a. */

/*       a is an m by n array. on input a contains the matrix for */
/*         which the qr factorization is to be computed. on output */
/*         the strict upper trapezoidal part of a contains the strict */
/*         upper trapezoidal part of r, and the lower trapezoidal */
/*         part of a contains a factored form of q (the non-trivial */
/*         elements of the u vectors described above). */

/*       lda is a positive integer input variable not less than m */
/*         which specifies the leading dimension of the array a. */

/*       pivot is a logical input variable. if pivot is set true, */
/*         then column pivoting is enforced. if pivot is set false, */
/*         then no column pivoting is done. */

/*       ipvt is an integer output array of length lipvt. ipvt */
/*         defines the permutation matrix p such that a*p = q*r. */
/*         column j of p is column ipvt(j) of the identity matrix. */
/*         if pivot is false, ipvt is not referenced. */

/*       lipvt is a positive integer input variable. if pivot is false, */
/*         then lipvt may be as small as 1. if pivot is true, then */
/*         lipvt must be at least n. */

/*       rdiag is an output array of length n which contains the */
/*         diagonal elements of r. */

/*       acnorm is an output array of length n which contains the */
/*         norms of the corresponding columns of the input matrix a. */
/*         if this information is not needed, then acnorm can coincide */
/*         with rdiag. */

/*       wa is a work array of length n. if pivot is false, then wa */
/*         can coincide with rdiag. */

/*       when built with USE_LAPACK, the QR is computed by LAPACK, which */
/*       needs more workspace than the length-n wa above (3*n+1 with pivoting, */
/*       n without). qrfac has no wa-length argument, so it cannot use extra */
/*       caller space; instead it obtains that workspace with an internal */
/*       malloc (once per call) and abort()s if the allocation fails. Callers */
/*       that must never abort on out-of-memory should not enable USE_LAPACK. */

/*     subprograms called */

/*       minpack-supplied ... dpmpar,enorm */

/*       fortran-supplied ... dmax1,dsqrt,min0 */

/*     argonne national laboratory. minpack project. march 1980. */
/*     burton s. garbow, kenneth e. hillstrom, jorge j. more */

/*     ********** */
    (void)lipvt;

/*     epsmch is the machine precision. */

    epsmch = __cminpack_func__(dpmpar)(1);

/*     compute the initial column norms and initialize several arrays. */

    for (j = 0; j < n; ++j) {
	acnorm[j] = __cminpack_func__(enorm)(m, &a[j * lda + 0]);
	rdiag[j] = acnorm[j];
	wa[j] = rdiag[j];
	if (pivot) {
	    ipvt[j] = j+1;
	}
    }

/*     reduce a to r with householder transformations. */

    minmn = min(m,n);
    for (j = 0; j < minmn; ++j) {
	if (pivot) {

/*        bring the column of largest norm into the pivot position. */

            int kmax = j;
            for (k = j; k < n; ++k) {
                if (rdiag[k] > rdiag[kmax]) {
                    kmax = k;
                }
            }
            if (kmax != j) {
                for (i = 0; i < m; ++i) {
                    temp = a[i + j * lda];
                    a[i + j * lda] = a[i + kmax * lda];
                    a[i + kmax * lda] = temp;
                }
                rdiag[kmax] = rdiag[j];
                wa[kmax] = wa[j];
                k = ipvt[j];
                ipvt[j] = ipvt[kmax];
                ipvt[kmax] = k;
            }
        }

/*        compute the householder transformation to reduce the */
/*        j-th column of a to a multiple of the j-th unit vector. */

	ajnorm = __cminpack_func__(enorm)(m - (j+1) + 1, &a[j + j * lda]);
	if (ajnorm != 0.) {
            if (a[j + j * lda] < 0.) {
                ajnorm = -ajnorm;
            }
            for (i = j; i < m; ++i) {
                a[i + j * lda] /= ajnorm;
            }
            a[j + j * lda] += 1;

/*        apply the transformation to the remaining columns */
/*        and update the norms. */

            jp1 = j + 1;
            if (n > jp1) {
                /* Transform CMINPACK_BLOCK columns at a time: their dot
                   products are independent, so interleaving them hides the
                   latency of the additions, and column j is read once per
                   block. Each sum is accumulated in the same order as one
                   column at a time, so the results are unchanged (unless the
                   compiler contracts expressions to FMA differently). */
                for (k = jp1; k + CMINPACK_BLOCK <= n; k += CMINPACK_BLOCK) {
                    real s[CMINPACK_BLOCK];
                    int q;
                    for (q = 0; q < CMINPACK_BLOCK; ++q) {
                        s[q] = 0.;
                    }
                    /* the loops over q must be unrolled, so that the sums
                       stay in registers; gcc does not do it by itself at -O2
                       (64: any block size up to 64 is unrolled fully) */
                    for (i = j; i < m; ++i) {
                        const real u = a[i + j * lda];
#ifdef __GNUC__
#pragma GCC unroll 64
#endif
                        for (q = 0; q < CMINPACK_BLOCK; ++q) {
                            s[q] += u * a[i + (k + q) * lda];
                        }
                    }
                    for (q = 0; q < CMINPACK_BLOCK; ++q) {
                        s[q] /= a[j + j * lda];
                    }
                    for (i = j; i < m; ++i) {
                        const real u = a[i + j * lda];
#ifdef __GNUC__
#pragma GCC unroll 64
#endif
                        for (q = 0; q < CMINPACK_BLOCK; ++q) {
                            a[i + (k + q) * lda] -= s[q] * u;
                        }
                    }
                }
                for (; k < n; ++k) {
                    sum = 0.;
                    for (i = j; i < m; ++i) {
                        sum += a[i + j * lda] * a[i + k * lda];
                    }
                    temp = sum / a[j + j * lda];
                    for (i = j; i < m; ++i) {
                        a[i + k * lda] -= temp * a[i + j * lda];
                    }
                }
                for (k = jp1; pivot && k < n; ++k) {
                    if (rdiag[k] != 0.) {
                        temp = a[j + k * lda] / rdiag[k];
                        /* Computing MAX */
                        d1 = 1 - temp * temp;
                        rdiag[k] *= sqrt((max((real)0.,d1)));
                        /* Computing 2nd power */
                        d1 = rdiag[k] / wa[k];
                        if (p05 * (d1 * d1) <= epsmch) {
                            rdiag[k] = __cminpack_func__(enorm)(m - (j+1), &a[jp1 + k * lda]);
                            wa[k] = rdiag[k];
                        }
                    }
                }
            }
        }
	rdiag[j] = -ajnorm;
    }

/*     last card of subroutine qrfac. */
#endif /* !USE_LAPACK */
} /* qrfac_ */

