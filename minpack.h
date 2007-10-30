#ifndef __MINPACK_H__
#define __MINPACK_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* Declarations for minpack */

/* find a zero of a system of N nonlinear functions in N variables by
   a modification of the Powell hybrid method (Jacobian calculated by
   a forward-difference approximation) */
void hybrd1_ ( void (*fcn)(const int *n, const double *x, double *fvec, int *iflag ), 
	       const int *n, double *x, double *fvec, const double *tol, int *info,
	       double *wa, const int *lwa );

/* find a zero of a system of N nonlinear functions in N variables by
   a modification of the Powell hybrid method (Jacobian calculated by
   a forward-difference approximation, more general). */
void hybrd_ ( void (*fcn)(const int *n, const double *x, double *fvec, int *iflag ), 
	      const int *n, double *x, double *fvec, const double *xtol, const int *maxfev,
	      const int *ml, const int *mu, const double *epsfcn, double *diag, const int *mode,
	      const double *factor, const int *nprint, int *info, int *nfev,
	      double *fjac, const int *ldfjac, double *r, const int *lr, double *qtf,
	      double *wa1, double *wa2, double *wa3, double *wa4);
  
/* find a zero of a system of N nonlinear functions in N variables by
   a modification of the Powell hybrid method (user-supplied Jacobian) */
void hybrj1_ ( void (*fcn)(const int *n, const double *x, double *fvec, double *fjec,
			   const int *ldfjac, int *iflag ), const int *n, double *x,
	       double *fvec, double *fjec, const int *ldfjac, const double *tol,
	       int *info, double *wa, const int *lwa );
          
/* find a zero of a system of N nonlinear functions in N variables by
   a modification of the Powell hybrid method (user-supplied Jacobian,
   more general) */
void hybrj_ ( void (*fcn)(const int *n, const double *x, double *fvec, double *fjec,
			  const int *ldfjac, int *iflag ), const int *n, double *x,
	      double *fvec, double *fjec, const int *ldfjac, const double *xtol,
	      const int *maxfev, double *diag, const int *mode, const double *factor,
	      const int *nprint, int *info, int *nfev, int *njev, double *r,
	      const int *lr, double *qtf, double *wa1, double *wa2,
	      double *wa3, double *wa4 );

/* minimize the sum of the squares of nonlinear functions in N
   variables by a modification of the Levenberg-Marquardt algorithm
   (Jacobian calculated by a forward-difference approximation) */
void lmdif1_ ( void (*fcn)(const int *m, const int *n, const double *x, double *fvec,
			   int *iflag ),
	       const int *m, const int *n, double *x, double *fvec, const double *tol,
	       int *info, int *iwa, double *wa, const int *lwa );

/* minimize the sum of the squares of nonlinear functions in N
   variables by a modification of the Levenberg-Marquardt algorithm
   (Jacobian calculated by a forward-difference approximation, more
   general) */
void lmdif_ ( void (*fcn)(const int *m, const int *n, const double *x, double *fvec,
			  int *iflag ),
	      const int *m, const int *n, double *x, double *fvec, const double *ftol,
	      const double *xtol, const double *gtol, const int *maxfev, const double *epsfcn,
	      double *diag, const int *mode, const double *factor, const int *nprint,
	      int *info, int *nfev, double *fjac, const int *ldfjac, int *ipvt,
	      double *qtf, double *wa1, double *wa2, double *wa3,
	      double *wa4 );

/* minimize the sum of the squares of nonlinear functions in N
   variables by a modification of the Levenberg-Marquardt algorithm
   (user-supplied Jacobian) */
void lmder1_ ( void (*fcn)(const int *m, const int *n, const double *x, double *fvec,
			   double *fjec, const int *ldfjac, int *iflag ),
	       const int *m, const int *n, double *x, double *fvec, double *fjac,
	       const int *ldfjac, const double *tol, int *info, int *ipvt,
	       double *wa, const int *lwa );

/* minimize the sum of the squares of nonlinear functions in N
   variables by a modification of the Levenberg-Marquardt algorithm
   (user-supplied Jacobian, more general) */
void lmder_ ( void (*fcn)(const int *m, const int *n, const double *x, double *fvec,
			  double *fjec, const int *ldfjac, int *iflag ),
	      const int *m, const int *n, double *x, double *fvec, double *fjac,
	      const int *ldfjac, const double *ftol, const double *xtol, const double *gtol,
	      const int *maxfev, double *diag, const int *mode, const double *factor,
	      const int *nprint, int *info, int *nfev, int *njev, int *ipvt,
	      double *qtf, double *wa1, double *wa2, double *wa3,
	      double *wa4 );

/* minimize the sum of the squares of nonlinear functions in N
   variables by a modification of the Levenberg-Marquardt algorithm
   (user-supplied Jacobian, minimal storage) */
void lmstr1_ ( void (*fcn)(const int *m, const int *n, const double *x, double *fvec,
			   double *fjrow, int *iflag ), const int *m, const int *n,
	       double *x, double *fvec, double *fjac, const int *ldfjac,
	       const double *tol, int *info, int *ipvt, double *wa, const int *lwa );

/* minimize the sum of the squares of nonlinear functions in N
   variables by a modification of the Levenberg-Marquardt algorithm
   (user-supplied Jacobian, minimal storage, more general) */
void lmstr_ ( void (*fcn)(const int *m, const int *n, const double *x, double *fvec,
			  double *fjrow, int *iflag ), const int *m,
	      const int *n, double *x, double *fvec, double *fjac,
	      const int *ldfjac, const double *ftol, const double *xtol, const double *gtol,
	      const int *maxfev, double *diag, const int *mode, const double *factor,
	      const int *nprint, int *info, int *nfev, int *njev, int *ipvt,
	      double *qtf, double *wa1, double *wa2, double *wa3,
	      double *wa4 );
 
void chkder_ ( const int *m, const int *n, const double *x, double *fvec, double *fjec,
	       const int *ldfjac, double *xp, double *fvecp, const int *mode,
	       double *err  );

double dpmpar_ ( const int *i );

double enorm_ ( const int *n, const double *x );

/* compute a forward-difference approximation to the m by n jacobian
   matrix associated with a specified problem of m functions in n
   variables. */
void fdjac2_(void (*fcn)(const int *m, const int *n, const double *x, double *fvec,
			 int *iflag ),
	     const int *m, const int *n, double *x, const double *fvec, double *fjac,
	     const int *ldfjac, int *iflag, const double *epsfcn, double *wa);

/* compute a forward-difference approximation to the n by n jacobian
   matrix associated with a specified problem of n functions in n
   variables. if the jacobian has a banded form, then function
   evaluations are saved by only approximating the nonzero terms. */
void fdjac1_(void (*fcn)(const int *n, const double *x, double *fvec, int *iflag ),
	     const int *n, double *x, const double *fvec, double *fjac, const int *ldfjac,
	     int *iflag, const int *ml, const int *mu, const double *epsfcn, double *wa1,
	     double *wa2);

/* internal MINPACK subroutines */
void dogleg_(const int *n, const double *r__, const int *lr, 
             const double *diag, const double *qtb, const double *delta, double *x, 
             double *wa1, double *wa2);
void qrfac_(const int *m, const int *n, double *a, const int *
            lda, const int *pivot, int *ipvt, const int *lipvt, double *rdiag,
            double *acnorm, double *wa);
void qrsolv_(const int *n, double *r__, const int *ldr, 
             const int *ipvt, const double *diag, const double *qtb, double *x, 
             double *sdiag, double *wa);
void qform_(const int *m, const int *n, double *q, const int *
            ldq, double *wa);
void r1updt_(const int *m, const int *n, double *s, const int *
             ls, const double *u, double *v, double *w, int *sing);
void r1mpyq_(const int *m, const int *n, double *a, const int *
             lda, const double *v, const double *w);
void lmpar_(const int *n, double *r__, const int *ldr, 
            const int *ipvt, const double *diag, const double *qtb, const double *delta, 
            double *par, double *x, double *sdiag, double *wa1, 
            double *wa2);
void rwupdt_(const int *n, double *r__, const int *ldr, 
             const double *w, double *b, double *alpha, double *cos__, 
             double *sin__);
void covar_(const int *n, double *r__, const int *ldr, 
           const int *ipvt, const double *tol, double *wa);
#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* __CMINPACK_H__ */
