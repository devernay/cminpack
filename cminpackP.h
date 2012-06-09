/* Internal header file for cminpack, by Frederic Devernay. */
#ifndef __CMINPACKP_H__
#define __CMINPACKP_H__

#ifndef __CMINPACK_H__
#error "cminpackP.h in an internal cminpack header, and must be included after all other headers (including cminpack.h)"
#endif

#define real __cminpack_real__
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define TRUE_ (1)
#define FALSE_ (0)

#if (defined (USE_CBLAS) || defined (USE_LAPACK)) && !defined (__cminpack_double__)
#error "cminpack can use cblas and lapack only in double precision mode"
#endif

#ifdef USE_CBLAS
#ifdef __APPLE__
#include <vecLib/cblas.h>
#else
#include <cblas.h>
#endif
#define __cminpack_enorm__(n,x) cblas_dnrm2(n,x,1)
#else
#define __cminpack_enorm__(n,x) __cminpack_func__(enorm)(n,x)
#endif

#endif /* !__CMINPACKP_H__ */
