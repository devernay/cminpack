/* lmstr.f -- translated by f2c (version 20020621).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <math.h>
#include "cminpack.h"
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE_ (1)
#define FALSE_ (0)

/* Subroutine */ int lmstr(minpack_funcderstr_mn fcn, void *p, int m, int n, double *x, 
	double *fvec, double *fjac, int ldfjac, double ftol,
	double xtol, double gtol, int maxfev, double *
	diag, int mode, double factor, int nprint,
	int *nfev, int *njev, int *ipvt, double *qtf, 
	double *wa1, double *wa2, double *wa3, double *wa4)
{
    /* Initialized data */

#define p1 .1
#define p5 .5
#define p25 .25
#define p75 .75
#define p0001 1e-4

    /* System generated locals */
    int fjac_dim1, fjac_offset;
    double d1, d2;

    /* Local variables */
    int i, j, l;
    double par, sum;
    int sing;
    int iter;
    double temp, temp1, temp2;
    int iflag;
    double delta;
    double ratio;
    double fnorm, gnorm, pnorm, xnorm, fnorm1, actred, dirder, 
	    epsmch, prered;
    int info;

/*     ********** */

/*     subroutine lmstr */

/*     the purpose of lmstr is to minimize the sum of the squares of */
/*     m nonlinear functions in n variables by a modification of */
/*     the levenberg-marquardt algorithm which uses minimal storage. */
/*     the user must provide a subroutine which calculates the */
/*     functions and the rows of the jacobian. */

/*     the subroutine statement is */

/*       subroutine lmstr(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol, */
/*                        maxfev,diag,mode,factor,nprint,info,nfev, */
/*                        njev,ipvt,qtf,wa1,wa2,wa3,wa4) */

/*     where */

/*       fcn is the name of the user-supplied subroutine which */
/*         calculates the functions and the rows of the jacobian. */
/*         fcn must be declared in an external statement in the */
/*         user calling program, and should be written as follows. */

/*         subroutine fcn(m,n,x,fvec,fjrow,iflag) */
/*         integer m,n,iflag */
/*         double precision x(n),fvec(m),fjrow(n) */
/*         ---------- */
/*         if iflag = 1 calculate the functions at x and */
/*         return this vector in fvec. */
/*         if iflag = i calculate the (i-1)-st row of the */
/*         jacobian at x and return this vector in fjrow. */
/*         ---------- */
/*         return */
/*         end */

/*         the value of iflag should not be changed by fcn unless */
/*         the user wants to terminate execution of lmstr. */
/*         in this case set iflag to a negative integer. */

/*       m is a positive integer input variable set to the number */
/*         of functions. */

/*       n is a positive integer input variable set to the number */
/*         of variables. n must not exceed m. */

/*       x is an array of length n. on input x must contain */
/*         an initial estimate of the solution vector. on output x */
/*         contains the final estimate of the solution vector. */

/*       fvec is an output array of length m which contains */
/*         the functions evaluated at the output x. */

/*       fjac is an output n by n array. the upper triangle of fjac */
/*         contains an upper triangular matrix r such that */

/*                t     t           t */
/*               p *(jac *jac)*p = r *r, */

/*         where p is a permutation matrix and jac is the final */
/*         calculated jacobian. column j of p is column ipvt(j) */
/*         (see below) of the identity matrix. the lower triangular */
/*         part of fjac contains information generated during */
/*         the computation of r. */

/*       ldfjac is a positive integer input variable not less than n */
/*         which specifies the leading dimension of the array fjac. */

/*       ftol is a nonnegative input variable. termination */
/*         occurs when both the actual and predicted relative */
/*         reductions in the sum of squares are at most ftol. */
/*         therefore, ftol measures the relative error desired */
/*         in the sum of squares. */

/*       xtol is a nonnegative input variable. termination */
/*         occurs when the relative error between two consecutive */
/*         iterates is at most xtol. therefore, xtol measures the */
/*         relative error desired in the approximate solution. */

/*       gtol is a nonnegative input variable. termination */
/*         occurs when the cosine of the angle between fvec and */
/*         any column of the jacobian is at most gtol in absolute */
/*         value. therefore, gtol measures the orthogonality */
/*         desired between the function vector and the columns */
/*         of the jacobian. */

/*       maxfev is a positive integer input variable. termination */
/*         occurs when the number of calls to fcn with iflag = 1 */
/*         has reached maxfev. */

/*       diag is an array of length n. if mode = 1 (see */
/*         below), diag is internally set. if mode = 2, diag */
/*         must contain positive entries that serve as */
/*         multiplicative scale factors for the variables. */

/*       mode is an integer input variable. if mode = 1, the */
/*         variables will be scaled internally. if mode = 2, */
/*         the scaling is specified by the input diag. other */
/*         values of mode are equivalent to mode = 1. */

/*       factor is a positive input variable used in determining the */
/*         initial step bound. this bound is set to the product of */
/*         factor and the euclidean norm of diag*x if nonzero, or else */
/*         to factor itself. in most cases factor should lie in the */
/*         interval (.1,100.). 100. is a generally recommended value. */

/*       nprint is an integer input variable that enables controlled */
/*         printing of iterates if it is positive. in this case, */
/*         fcn is called with iflag = 0 at the beginning of the first */
/*         iteration and every nprint iterations thereafter and */
/*         immediately prior to return, with x and fvec available */
/*         for printing. if nprint is not positive, no special calls */
/*         of fcn with iflag = 0 are made. */

/*       info is an integer output variable. if the user has */
/*         terminated execution, info is set to the (negative) */
/*         value of iflag. see description of fcn. otherwise, */
/*         info is set as follows. */

/*         info = 0  improper input parameters. */

/*         info = 1  both actual and predicted relative reductions */
/*                   in the sum of squares are at most ftol. */

/*         info = 2  relative error between two consecutive iterates */
/*                   is at most xtol. */

/*         info = 3  conditions for info = 1 and info = 2 both hold. */

/*         info = 4  the cosine of the angle between fvec and any */
/*                   column of the jacobian is at most gtol in */
/*                   absolute value. */

/*         info = 5  number of calls to fcn with iflag = 1 has */
/*                   reached maxfev. */

/*         info = 6  ftol is too small. no further reduction in */
/*                   the sum of squares is possible. */

/*         info = 7  xtol is too small. no further improvement in */
/*                   the approximate solution x is possible. */

/*         info = 8  gtol is too small. fvec is orthogonal to the */
/*                   columns of the jacobian to machine precision. */

/*       nfev is an integer output variable set to the number of */
/*         calls to fcn with iflag = 1. */

/*       njev is an integer output variable set to the number of */
/*         calls to fcn with iflag = 2. */

/*       ipvt is an integer output array of length n. ipvt */
/*         defines a permutation matrix p such that jac*p = q*r, */
/*         where jac is the final calculated jacobian, q is */
/*         orthogonal (not stored), and r is upper triangular. */
/*         column j of p is column ipvt(j) of the identity matrix. */

/*       qtf is an output array of length n which contains */
/*         the first n elements of the vector (q transpose)*fvec. */

/*       wa1, wa2, and wa3 are work arrays of length n. */

/*       wa4 is a work array of length m. */

/*     subprograms called */

/*       user-supplied ...... fcn */

/*       minpack-supplied ... dpmpar,enorm,lmpar,qrfac,rwupdt */

/*       fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod */

/*     argonne national laboratory. minpack project. march 1980. */
/*     burton s. garbow, dudley v. goetschel, kenneth e. hillstrom, */
/*     jorge j. more */

/*     ********** */
    /* Parameter adjustments */
    --wa4;
    --fvec;
    --wa3;
    --wa2;
    --wa1;
    --qtf;
    --ipvt;
    --diag;
    --x;
    fjac_dim1 = ldfjac;
    fjac_offset = 1 + fjac_dim1 * 1;
    fjac -= fjac_offset;

    /* Function Body */

/*     epsmch is the machine precision. */

    epsmch = dpmpar(1);

    info = 0;
    iflag = 0;
    *nfev = 0;
    *njev = 0;

/*     check the input parameters for errors. */

    if (n <= 0 || m < n || ldfjac < n || ftol < 0. || xtol < 0. || 
	    gtol < 0. || maxfev <= 0 || factor <= 0.) {
	goto TERMINATE;
    }
    if (mode == 2) {
        for (j = 1; j <= n; ++j) {
            if (diag[j] <= 0.) {
                goto TERMINATE;
            }
        }
    }

/*     evaluate the function at the starting point */
/*     and calculate its norm. */

    iflag = (*fcn)(p, m, n, &x[1], &fvec[1], &wa3[1], 1);
    *nfev = 1;
    if (iflag < 0) {
	goto TERMINATE;
    }
    fnorm = enorm(m, &fvec[1]);

/*     initialize levenberg-marquardt parameter and iteration counter. */

    par = 0.;
    iter = 1;

/*     beginning of the outer loop. */

    for (;;) {

/*        if requested, call fcn to enable printing of iterates. */

        if (nprint > 0) {
            iflag = 0;
            if ((iter - 1) % nprint == 0) {
                iflag = (*fcn)(p, m, n, &x[1], &fvec[1], &wa3[1], 0);
            }
            if (iflag < 0) {
                goto TERMINATE;
            }
        }

/*        compute the qr factorization of the jacobian matrix */
/*        calculated one row at a time, while simultaneously */
/*        forming (q transpose)*fvec and storing the first */
/*        n components in qtf. */

        for (j = 1; j <= n; ++j) {
            qtf[j] = 0.;
            for (i = 1; i <= n; ++i) {
                fjac[i + j * fjac_dim1] = 0.;
            }
        }
        iflag = 2;
        for (i = 1; i <= m; ++i) {
            if ((*fcn)(p, m, n, &x[1], &fvec[1], &wa3[1], iflag) < 0) {
                goto TERMINATE;
            }
            temp = fvec[i];
            rwupdt(n, &fjac[fjac_offset], ldfjac, &wa3[1], &qtf[1], &temp,
                   &wa1[1], &wa2[1]);
            ++iflag;
        }
        ++(*njev);

/*        if the jacobian is rank deficient, call qrfac to */
/*        reorder its columns and update the components of qtf. */

        sing = FALSE_;
        for (j = 1; j <= n; ++j) {
            if (fjac[j + j * fjac_dim1] == 0.) {
                sing = TRUE_;
            }
            ipvt[j] = j;
            wa2[j] = enorm(j, &fjac[j * fjac_dim1 + 1]);
        }
        if (sing) {
            qrfac(n, n, &fjac[fjac_offset], ldfjac, TRUE_, &ipvt[1], n,
                  &wa1[1], &wa2[1], &wa3[1]);
            for (j = 1; j <= n; ++j) {
                if (fjac[j + j * fjac_dim1] != 0.) {
                    sum = 0.;
                    for (i = j; i <= n; ++i) {
                        sum += fjac[i + j * fjac_dim1] * qtf[i];
                    }
                    temp = -sum / fjac[j + j * fjac_dim1];
                    for (i = j; i <= n; ++i) {
                        qtf[i] += fjac[i + j * fjac_dim1] * temp;
                    }
                }
                fjac[j + j * fjac_dim1] = wa1[j];
            }
        }

/*        on the first iteration and if mode is 1, scale according */
/*        to the norms of the columns of the initial jacobian. */

        if (iter == 1) {
            if (mode != 2) {
                for (j = 1; j <= n; ++j) {
                    diag[j] = wa2[j];
                    if (wa2[j] == 0.) {
                        diag[j] = 1.;
                    }
                }
            }

/*        on the first iteration, calculate the norm of the scaled x */
/*        and initialize the step bound delta. */

            for (j = 1; j <= n; ++j) {
                wa3[j] = diag[j] * x[j];
            }
            xnorm = enorm(n, &wa3[1]);
            delta = factor * xnorm;
            if (delta == 0.) {
                delta = factor;
            }
        }

/*        compute the norm of the scaled gradient. */

        gnorm = 0.;
        if (fnorm != 0.) {
            for (j = 1; j <= n; ++j) {
                l = ipvt[j];
                if (wa2[l] != 0.) {
                    sum = 0.;
                    for (i = 1; i <= j; ++i) {
                        sum += fjac[i + j * fjac_dim1] * (qtf[i] / fnorm);
                    }
                    /* Computing MAX */
                    d1 = fabs(sum / wa2[l]);
                    gnorm = max(gnorm,d1);
                }
            }
        }

/*        test for convergence of the gradient norm. */

        if (gnorm <= gtol) {
            info = 4;
        }
        if (info != 0) {
            goto TERMINATE;
        }

/*        rescale if necessary. */

        if (mode != 2) {
            for (j = 1; j <= n; ++j) {
                /* Computing MAX */
                d1 = diag[j], d2 = wa2[j];
                diag[j] = max(d1,d2);
            }
        }

/*        beginning of the inner loop. */

        do {

/*           determine the levenberg-marquardt parameter. */

            lmpar(n, &fjac[fjac_offset], ldfjac, &ipvt[1], &diag[1], &qtf[1], delta,
                  &par, &wa1[1], &wa2[1], &wa3[1], &wa4[1]);

/*           store the direction p and x + p. calculate the norm of p. */

            for (j = 1; j <= n; ++j) {
                wa1[j] = -wa1[j];
                wa2[j] = x[j] + wa1[j];
                wa3[j] = diag[j] * wa1[j];
            }
            pnorm = enorm(n, &wa3[1]);

/*           on the first iteration, adjust the initial step bound. */

            if (iter == 1) {
                delta = min(delta,pnorm);
            }

/*           evaluate the function at x + p and calculate its norm. */

            iflag = (*fcn)(p, m, n, &wa2[1], &wa4[1], &wa3[1], 1);
            ++(*nfev);
            if (iflag < 0) {
                goto TERMINATE;
            }
            fnorm1 = enorm(m, &wa4[1]);

/*           compute the scaled actual reduction. */

            actred = -1.;
            if (p1 * fnorm1 < fnorm) {
                /* Computing 2nd power */
                d1 = fnorm1 / fnorm;
                actred = 1. - d1 * d1;
            }

/*           compute the scaled predicted reduction and */
/*           the scaled directional derivative. */

            for (j = 1; j <= n; ++j) {
                wa3[j] = 0.;
                l = ipvt[j];
                temp = wa1[l];
                for (i = 1; i <= j; ++i) {
                    wa3[i] += fjac[i + j * fjac_dim1] * temp;
                }
            }
            temp1 = enorm(n, &wa3[1]) / fnorm;
            temp2 = (sqrt(par) * pnorm) / fnorm;
            prered = temp1 * temp1 + temp2 * temp2 / p5;
            dirder = -(temp1 * temp1 + temp2 * temp2);

/*           compute the ratio of the actual to the predicted */
/*           reduction. */

            ratio = 0.;
            if (prered != 0.) {
                ratio = actred / prered;
            }

/*           update the step bound. */

            if (ratio <= p25) {
                if (actred >= 0.) {
                    temp = p5;
                } else {
                    temp = p5 * dirder / (dirder + p5 * actred);
                }
                if (p1 * fnorm1 >= fnorm || temp < p1) {
                    temp = p1;
                }
                /* Computing MIN */
                d1 = pnorm / p1;
                delta = temp * min(delta,d1);
                par /= temp;
            } else {
                if (par == 0. || ratio >= p75) {
                    delta = pnorm / p5;
                    par = p5 * par;
                }
            }

/*           test for successful iteration. */

            if (ratio >= p0001) {

/*           successful iteration. update x, fvec, and their norms. */

                for (j = 1; j <= n; ++j) {
                    x[j] = wa2[j];
                    wa2[j] = diag[j] * x[j];
                }
                for (i = 1; i <= m; ++i) {
                    fvec[i] = wa4[i];
                }
                xnorm = enorm(n, &wa2[1]);
                fnorm = fnorm1;
                ++iter;
            }

/*           tests for convergence. */

            if (fabs(actred) <= ftol && prered <= ftol && p5 * ratio <= 1.) {
                info = 1;
            }
            if (delta <= xtol * xnorm) {
                info = 2;
            }
            if (fabs(actred) <= ftol && prered <= ftol && p5 * ratio <= 1. && info == 2) {
                info = 3;
            }
            if (info != 0) {
                goto TERMINATE;
            }

/*           tests for termination and stringent tolerances. */

            if (*nfev >= maxfev) {
                info = 5;
            }
            if (fabs(actred) <= epsmch && prered <= epsmch && p5 * ratio <= 1.) {
                info = 6;
            }
            if (delta <= epsmch * xnorm) {
                info = 7;
            }
            if (gnorm <= epsmch) {
                info = 8;
            }
            if (info != 0) {
                goto TERMINATE;
            }

/*           end of the inner loop. repeat if iteration unsuccessful. */

        } while (ratio < p0001);

/*        end of the outer loop. */

    }
TERMINATE:

/*     termination, either normal or user imposed. */

    if (iflag < 0) {
	info = iflag;
    }
    if (nprint > 0) {
	(*fcn)(p, m, n, &x[1], &fvec[1], &wa3[1], 0);
    }
    return info;

/*     last card of subroutine lmstr. */

} /* lmstr_ */

