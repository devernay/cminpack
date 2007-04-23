/* Declarations for minpack */

/* find a zero of a system of N nonlinear functions in N variables by
   a modification of the Powell hybrid method (Jacobian calculated by
   a forward-difference approximation) */
void hybrd1 ( void (*fcn)(int n, const double *x, double *fvec, int *iflag ), 
	       int n, double *x, double *fvec, double tol, int *info,
	       double *wa, int lwa );

/* find a zero of a system of N nonlinear functions in N variables by
   a modification of the Powell hybrid method (Jacobian calculated by
   a forward-difference approximation, more general). */
void hybrd ( void (*fcn)(int n, const double *x, double *fvec, int *iflag ), 
	      int n, double *x, double *fvec, double xtol, int maxfev,
	      int ml, int mu, double epsfcn, double *diag, int mode,
	      double factor, int nprint, int *info, int *nfev,
	      double *fjac, int ldfjac, double *r, int lr, double *qtf,
	      double *wa1, double *wa2, double *wa3, double *wa4);
  
/* find a zero of a system of N nonlinear functions in N variables by
   a modification of the Powell hybrid method (user-supplied Jacobian) */
void hybrj1 ( void (*fcn)(int n, const double *x, double *fvec, double *fjac,
			   int ldfjac, int *iflag ), int n, double *x,
	       double *fvec, double *fjac, int ldfjac, double tol,
	       int *info, double *wa, int lwa );
          
/* find a zero of a system of N nonlinear functions in N variables by
   a modification of the Powell hybrid method (user-supplied Jacobian,
   more general) */
void hybrj ( void (*fcn)(int n, const double *x, double *fvec, double *fjac,
			  int ldfjac, int *iflag ), int n, double *x,
	      double *fvec, double *fjac, int ldfjac, double xtol,
	      int maxfev, double *diag, int mode, double factor,
	      int nprint, int *info, int *nfev, int *njev, double *r,
	      int lr, double *qtf, double *wa1, double *wa2,
	      double *wa3, double *wa4 );

/* minimize the sum of the squares of nonlinear functions in N
   variables by a modification of the Levenberg-Marquardt algorithm
   (Jacobian calculated by a forward-difference approximation) */
void lmdif1 ( void (*fcn)(int m, int n, const double *x, double *fvec,
			   int *iflag ),
	       int m, int n, double *x, double *fvec, double tol,
	       int *info, int *iwa, double *wa, int lwa );

/* minimize the sum of the squares of nonlinear functions in N
   variables by a modification of the Levenberg-Marquardt algorithm
   (Jacobian calculated by a forward-difference approximation, more
   general) */
void lmdif ( void (*fcn)(int m, int n, const double *x, double *fvec,
			  int *iflag ),
	      int m, int n, double *x, double *fvec, double ftol,
	      double xtol, double gtol, int maxfev, double epsfcn,
	      double *diag, int mode, double factor, int nprint,
	      int *info, int *nfev, double *fjac, int ldfjac, int *ipvt,
	      double *qtf, double *wa1, double *wa2, double *wa3,
	      double *wa4 );

/* minimize the sum of the squares of nonlinear functions in N
   variables by a modification of the Levenberg-Marquardt algorithm
   (user-supplied Jacobian) */
void lmder1 ( void (*fcn)(int m, int n, const double *x, double *fvec,
			   double *fjac, int ldfjac, int *iflag ),
	       int m, int n, double *x, double *fvec, double *fjac,
	       int ldfjac, double tol, int *info, int *ipvt,
	       double *wa, int lwa );

/* minimize the sum of the squares of nonlinear functions in N
   variables by a modification of the Levenberg-Marquardt algorithm
   (user-supplied Jacobian, more general) */
void lmder ( void (*fcn)(int m, int n, const double *x, double *fvec,
			  double *fjac, int ldfjac, int *iflag ),
	      int m, int n, double *x, double *fvec, double *fjac,
	      int ldfjac, double ftol, double xtol, double gtol,
	      int maxfev, double *diag, int mode, double factor,
	      int nprint, int *info, int *nfev, int *njev, int *ipvt,
	      double *qtf, double *wa1, double *wa2, double *wa3,
	      double *wa4 );

/* minimize the sum of the squares of nonlinear functions in N
   variables by a modification of the Levenberg-Marquardt algorithm
   (user-supplied Jacobian, minimal storage) */
void lmstr1 ( void (*fcn)(int m, int n, const double *x, double *fvec,
			   double *fjrow, int *iflag ), int m, int n,
	       double *x, double *fvec, double *fjac, int ldfjac,
	       double tol, int *info, int *ipvt, double *wa, int lwa );

/* minimize the sum of the squares of nonlinear functions in N
   variables by a modification of the Levenberg-Marquardt algorithm
   (user-supplied Jacobian, minimal storage, more general) */
void lmstr ( void (*fcn)(int m, int n, const double *x, double *fvec,
			  double *fjrow, int *iflag ), int m,
	      int n, double *x, double *fvec, double *fjac,
	      int ldfjac, double ftol, double xtol, double gtol,
	      int maxfev, double *diag, int mode, double factor,
	      int nprint, int *info, int *nfev, int *njev, int *ipvt,
	      double *qtf, double *wa1, double *wa2, double *wa3,
	      double *wa4 );
 
void chkder ( int m, int n, const double *x, double *fvec, double *fjac,
	       int ldfjac, double *xp, double *fvecp, int mode,
	       double *err  );

double dpmpar ( int i );

double enorm ( int n, const double *x );

/* compute a forward-difference approximation to the m by n jacobian
   matrix associated with a specified problem of m functions in n
   variables. */
void fdjac2(void (*fcn)(int m, int n, const double *x, double *fvec,
			 int *iflag ),
	     int m, int n, double *x, const double *fvec, double *fjac,
	     int ldfjac, int *iflag, double epsfcn, double *wa);

/* compute a forward-difference approximation to the n by n jacobian
   matrix associated with a specified problem of n functions in n
   variables. if the jacobian has a banded form, then function
   evaluations are saved by only approximating the nonzero terms. */
void fdjac1(void (*fcn)(int n, const double *x, double *fvec, int *iflag ),
	     int n, double *x, const double *fvec, double *fjac, int ldfjac,
	     int *iflag, int ml, int mu, double epsfcn, double *wa1,
	     double *wa2);

/* internal MINPACK subroutines */
void dogleg(int n, const double *r__, int lr, 
             const double *diag, const double *qtb, double delta, double *x, 
             double *wa1, double *wa2);
void qrfac(int m, int n, double *a, int
            lda, int pivot, int *ipvt, int lipvt, double *rdiag,
            double *acnorm, double *wa);
void qrsolv(int n, double *r__, int ldr, 
             const int *ipvt, const double *diag, const double *qtb, double *x, 
             double *sdiag, double *wa);
void qform(int m, int n, double *q, int
            ldq, double *wa);
void r1updt(int m, int n, double *s, int
             ls, const double *u, double *v, double *w, int *sing);
void r1mpyq(int m, int n, double *a, int
             lda, const double *v, const double *w);
void lmpar(int n, double *r__, int ldr, 
            const int *ipvt, const double *diag, const double *qtb, double delta, 
            double *par, double *x, double *sdiag, double *wa1, 
            double *wa2);
void rwupdt(int n, double *r__, int ldr, 
             const double *w, double *b, double *alpha, double *cos__, 
             double *sin__);
