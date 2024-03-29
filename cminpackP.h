/* Internal header file for cminpack, by Frederic Devernay. */
#ifndef __CMINPACKP_H__
#define __CMINPACKP_H__

#ifndef __CMINPACK_H__
#error "cminpackP.h in an internal cminpack header, and must be included after all other headers (including cminpack.h)"
#endif

#if (defined (USE_BLAS) || defined (USE_LAPACK)) && !defined (__cminpack_double__) && !defined (__cminpack_float__)
#error "cminpack can use cblas and lapack only in double or single precision mode"
#endif

#ifdef USE_BLAS
int __cminpack_blas__(dot)(
  const int N, const __cminpack_real__ *X, const int incX,
  const __cminpack_real__ *Y, const int incY);
int __cminpack_blas__(nrm2)(
  const int N, const __cminpack_real__ *X, const int incX);
int __cminpack_blas__(swap)(
  const int N, __cminpack_real__ *X, const int incX,
  __cminpack_real__ *Y, const int incY);
int __cminpack_blas__(rot)(
  const int N, __cminpack_real__ *X, const int incX,
  __cminpack_real__ *Y, const int incY, const __cminpack_real__ c, const __cminpack_real__ s);
int __cminpack_blas__(trsv)(
  const char *Uplo,
  const char *TransA, const char *Diag,
  const int N, const __cminpack_real__ *A, const int lda, __cminpack_real__ *X,
  const int incX);
#define __cminpack_enorm__(n,x) __cminpack_blas__(nrm2)(n,x,1)
#else
#define __cminpack_enorm__(n,x) __cminpack_func__(enorm)(n,x)
#endif

#ifdef USE_LAPACK
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#if defined(__LP64__) /* In LP64 match sizes with the 32 bit ABI */
typedef int 		__CLPK_integer;
typedef int 		__CLPK_logical;
typedef __CLPK_logical 	(*__CLPK_L_fp)();
typedef int 		__CLPK_ftnlen;
#else
typedef long int 	__CLPK_integer;
typedef long int 	__CLPK_logical;
typedef __CLPK_logical 	(*__CLPK_L_fp)();
typedef long int 	__CLPK_ftnlen;
#endif
int __cminpack_lapack__(lartg_)(
  __cminpack_real__ *f, __cminpack_real__ *g, __cminpack_real__ *cs,
  __cminpack_real__ *sn, __cminpack_real__ *r__);
int __cminpack_lapack__(geqp3_)(
  __CLPK_integer *m, __CLPK_integer *n, __cminpack_real__ *a, __CLPK_integer * lda,
  __CLPK_integer *jpvt, __cminpack_real__ *tau, __cminpack_real__ *work, __CLPK_integer *lwork,
  __CLPK_integer *info);
int __cminpack_lapack__(geqrf_)(
  __CLPK_integer *m, __CLPK_integer *n, __cminpack_real__ *a, __CLPK_integer * lda,
  __cminpack_real__ *tau, __cminpack_real__ *work, __CLPK_integer *lwork, __CLPK_integer *info);
#endif
#endif

#include "minpackP.h"

#endif /* !__CMINPACKP_H__ */
