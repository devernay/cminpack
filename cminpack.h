#ifndef __CMINPACK_H__
#define __CMINPACK_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* Cmake will define cminpack_EXPORTS on Windows when it
configures to build a shared library. If you are going to use
another build system on windows or create the visual studio
projects by hand you need to define cminpack_EXPORTS when
building a DLL on windows.
*/
#if defined (__GNUC__)
#define CMINPACK_DECLSPEC_EXPORT  __declspec(__dllexport__)
#define CMINPACK_DECLSPEC_IMPORT  __declspec(__dllimport__)
#endif
#if defined (_MSC_VER) || defined (__BORLANDC__)
#define CMINPACK_DECLSPEC_EXPORT  __declspec(dllexport)
#define CMINPACK_DECLSPEC_IMPORT  __declspec(dllimport)
#endif
#ifdef __WATCOMC__
#define CMINPACK_DECLSPEC_EXPORT  __export
#define CMINPACK_DECLSPEC_IMPORT  __import
#endif
#ifdef __IBMC__
#define CMINPACK_DECLSPEC_EXPORT  _Export
#define CMINPACK_DECLSPEC_IMPORT  _Import
#endif

#if !defined(CMINPACK_NO_DLL) && (defined(__WIN32__) || defined(WIN32) || defined (_WIN32))
#if defined(cminpack_EXPORTS) || defined(CMINPACK_EXPORTS) || defined(CMINPACK_DLL_EXPORTS)
    #define  CMINPACK_EXPORT CMINPACK_DECLSPEC_EXPORT
  #else
    #define  CMINPACK_EXPORT CMINPACK_DECLSPEC_IMPORT
  #endif /* cminpack_EXPORTS */
#else /* defined (_WIN32) */
 #define CMINPACK_EXPORT
#endif

/* Declarations for minpack */

/* Function types: */
/* The first argument can be used to store extra function parameters, thus */
/* avoiding the use of global variables. */
/* the iflag parameter is input-only (with respect to the FORTRAN */
/*  version), the output iflag value is the return value of the function. */
/* If iflag=0, the function shoulkd just print the current values (see */
/* the nprint parameters below). */
  
/* for hybrd1 and hybrd: */
/*         calculate the functions at x and */
/*         return this vector in fvec. */
/* return a negative value to terminate hybrd1/hybrd */
typedef int (*minpack_func_nn)(void *p, int n, const double *x, double *fvec, int iflag );

/* for hybrj1 and hybrj */
/*         if iflag = 1 calculate the functions at x and */
/*         return this vector in fvec. do not alter fjac. */
/*         if iflag = 2 calculate the jacobian at x and */
/*         return this matrix in fjac. do not alter fvec. */
/* return a negative value to terminate hybrj1/hybrj */
typedef int (*minpack_funcder_nn)(void *p, int n, const double *x, double *fvec, double *fjac,
                                  int ldfjac, int iflag );

/* for lmdif1 and lmdif */
/*         calculate the functions at x and */
/*         return this vector in fvec. */
/* return a negative value to terminate lmdif1/lmdif */
typedef int (*minpack_func_mn)(void *p, int m, int n, const double *x, double *fvec,
                               int iflag );

/* for lmder1 and lmder */
/*         if iflag = 1 calculate the functions at x and */
/*         return this vector in fvec. do not alter fjac. */
/*         if iflag = 2 calculate the jacobian at x and */
/*         return this matrix in fjac. do not alter fvec. */
/* return a negative value to terminate lmder1/lmder */
typedef int (*minpack_funcder_mn)(void *p, int m, int n, const double *x, double *fvec,
                                  double *fjac, int ldfjac, int iflag );

/* for lmstr1 and lmstr */
/*         if iflag = 1 calculate the functions at x and */
/*         return this vector in fvec. */
/*         if iflag = i calculate the (i-1)-st row of the */
/*         jacobian at x and return this vector in fjrow. */
/* return a negative value to terminate lmstr1/lmstr */
typedef int (*minpack_funcderstr_mn)(void *p, int m, int n, const double *x, double *fvec,
                                     double *fjrow, int iflag );






/* MINPACK functions: */
/* the info parameter was removed from most functions: the return */
/* value of the function is used instead. */
/* The argument 'p' can be used to store extra function parameters, thus */
/* avoiding the use of global variables. You can also think of it as a */
/* 'this' pointer a la C++. */

/* find a zero of a system of N nonlinear functions in N variables by
   a modification of the Powell hybrid method (Jacobian calculated by
   a forward-difference approximation) */
int CMINPACK_EXPORT hybrd1 ( minpack_func_nn fcn, 
	       void *p, int n, double *x, double *fvec, double tol,
	       double *wa, int lwa );

/* find a zero of a system of N nonlinear functions in N variables by
   a modification of the Powell hybrid method (Jacobian calculated by
   a forward-difference approximation, more general). */
int CMINPACK_EXPORT hybrd ( minpack_func_nn fcn,
	      void *p, int n, double *x, double *fvec, double xtol, int maxfev,
	      int ml, int mu, double epsfcn, double *diag, int mode,
	      double factor, int nprint, int *nfev,
	      double *fjac, int ldfjac, double *r, int lr, double *qtf,
	      double *wa1, double *wa2, double *wa3, double *wa4);
  
/* find a zero of a system of N nonlinear functions in N variables by
   a modification of the Powell hybrid method (user-supplied Jacobian) */
int CMINPACK_EXPORT hybrj1 ( minpack_funcder_nn fcn, void *p, int n, double *x,
	       double *fvec, double *fjac, int ldfjac, double tol,
	       double *wa, int lwa );
          
/* find a zero of a system of N nonlinear functions in N variables by
   a modification of the Powell hybrid method (user-supplied Jacobian,
   more general) */
int CMINPACK_EXPORT hybrj ( minpack_funcder_nn fcn, void *p, int n, double *x,
	      double *fvec, double *fjac, int ldfjac, double xtol,
	      int maxfev, double *diag, int mode, double factor,
	      int nprint, int *nfev, int *njev, double *r,
	      int lr, double *qtf, double *wa1, double *wa2,
	      double *wa3, double *wa4 );

/* minimize the sum of the squares of nonlinear functions in N
   variables by a modification of the Levenberg-Marquardt algorithm
   (Jacobian calculated by a forward-difference approximation) */
int CMINPACK_EXPORT lmdif1 ( minpack_func_mn fcn,
	       void *p, int m, int n, double *x, double *fvec, double tol,
	       int *iwa, double *wa, int lwa );

/* minimize the sum of the squares of nonlinear functions in N
   variables by a modification of the Levenberg-Marquardt algorithm
   (Jacobian calculated by a forward-difference approximation, more
   general) */
int CMINPACK_EXPORT lmdif ( minpack_func_mn fcn,
	      void *p, int m, int n, double *x, double *fvec, double ftol,
	      double xtol, double gtol, int maxfev, double epsfcn,
	      double *diag, int mode, double factor, int nprint,
	      int *nfev, double *fjac, int ldfjac, int *ipvt,
	      double *qtf, double *wa1, double *wa2, double *wa3,
	      double *wa4 );

/* minimize the sum of the squares of nonlinear functions in N
   variables by a modification of the Levenberg-Marquardt algorithm
   (user-supplied Jacobian) */
int CMINPACK_EXPORT lmder1 ( minpack_funcder_mn fcn,
	       void *p, int m, int n, double *x, double *fvec, double *fjac,
	       int ldfjac, double tol, int *ipvt,
	       double *wa, int lwa );

/* minimize the sum of the squares of nonlinear functions in N
   variables by a modification of the Levenberg-Marquardt algorithm
   (user-supplied Jacobian, more general) */
int CMINPACK_EXPORT lmder ( minpack_funcder_mn fcn,
	      void *p, int m, int n, double *x, double *fvec, double *fjac,
	      int ldfjac, double ftol, double xtol, double gtol,
	      int maxfev, double *diag, int mode, double factor,
	      int nprint, int *nfev, int *njev, int *ipvt,
	      double *qtf, double *wa1, double *wa2, double *wa3,
	      double *wa4 );

/* minimize the sum of the squares of nonlinear functions in N
   variables by a modification of the Levenberg-Marquardt algorithm
   (user-supplied Jacobian, minimal storage) */
int CMINPACK_EXPORT lmstr1 ( minpack_funcderstr_mn fcn, void *p, int m, int n,
	       double *x, double *fvec, double *fjac, int ldfjac,
	       double tol, int *ipvt, double *wa, int lwa );

/* minimize the sum of the squares of nonlinear functions in N
   variables by a modification of the Levenberg-Marquardt algorithm
   (user-supplied Jacobian, minimal storage, more general) */
int CMINPACK_EXPORT lmstr (  minpack_funcderstr_mn fcn, void *p, int m,
	      int n, double *x, double *fvec, double *fjac,
	      int ldfjac, double ftol, double xtol, double gtol,
	      int maxfev, double *diag, int mode, double factor,
	      int nprint, int *nfev, int *njev, int *ipvt,
	      double *qtf, double *wa1, double *wa2, double *wa3,
	      double *wa4 );
 
void CMINPACK_EXPORT chkder ( int m, int n, const double *x, double *fvec, double *fjac,
	       int ldfjac, double *xp, double *fvecp, int mode,
	       double *err  );

double CMINPACK_EXPORT dpmpar ( int i );

double CMINPACK_EXPORT enorm ( int n, const double *x );

/* compute a forward-difference approximation to the m by n jacobian
   matrix associated with a specified problem of m functions in n
   variables. */
int CMINPACK_EXPORT fdjac2(minpack_func_mn fcn,
	     void *p, int m, int n, double *x, const double *fvec, double *fjac,
	     int ldfjac, double epsfcn, double *wa);

/* compute a forward-difference approximation to the n by n jacobian
   matrix associated with a specified problem of n functions in n
   variables. if the jacobian has a banded form, then function
   evaluations are saved by only approximating the nonzero terms. */
int CMINPACK_EXPORT fdjac1(minpack_func_nn fcn,
	     void *p, int n, double *x, const double *fvec, double *fjac, int ldfjac,
	     int ml, int mu, double epsfcn, double *wa1,
	     double *wa2);

/* compute inverse(JtJ) after a run of lmdif or lmder. The covariance matrix is obtained
   by scaling the result by enorm(y)**2/(m-n). If JtJ is singular and k = rank(J), the
   pseudo-inverse is computed, and the result has to be scaled by enorm(y)**2/(m-k). */
void CMINPACK_EXPORT covar(int n, double *r, int ldr, 
           const int *ipvt, double tol, double *wa);

/* covar1 estimates the variance-covariance matrix:
   C = sigma**2 (JtJ)**+
   where (JtJ)**+ is the inverse of JtJ or the pseudo-inverse of JtJ (in case J does not have full rank),
   and sigma**2 = fsumsq / (m - k)
   where fsumsq is the residual sum of squares and k is the rank of J.
   The function returns 0 if J has full rank, else the rank of J.
*/
int CMINPACK_EXPORT covar1(int m, int n, double fsumsq, double *r, int ldr, 
                           const int *ipvt, double tol, double *wa);

/* internal MINPACK subroutines */
void dogleg(int n, const double *r, int lr, 
             const double *diag, const double *qtb, double delta, double *x, 
             double *wa1, double *wa2);
void qrfac(int m, int n, double *a, int
            lda, int pivot, int *ipvt, int lipvt, double *rdiag,
            double *acnorm, double *wa);
void qrsolv(int n, double *r, int ldr, 
             const int *ipvt, const double *diag, const double *qtb, double *x, 
             double *sdiag, double *wa);
void qform(int m, int n, double *q, int
            ldq, double *wa);
void r1updt(int m, int n, double *s, int
             ls, const double *u, double *v, double *w, int *sing);
void r1mpyq(int m, int n, double *a, int
             lda, const double *v, const double *w);
void lmpar(int n, double *r, int ldr, 
            const int *ipvt, const double *diag, const double *qtb, double delta, 
            double *par, double *x, double *sdiag, double *wa1, 
            double *wa2);
void rwupdt(int n, double *r, int ldr, 
             const double *w, double *b, double *alpha, double *cos, 
             double *sin);
#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* __CMINPACK_H__ */
