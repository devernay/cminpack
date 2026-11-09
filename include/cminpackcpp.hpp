/* cminpackcpp.hpp -- a thin, header-only C++ wrapper around cminpack.

   The C API passes the user callback as a plain function pointer plus an
   opaque `void *p` for user data. That makes it awkward to use a C++ lambda,
   a functor, or a std::function, because none of those convert to a plain
   function pointer once they capture state (issue #74).

   This header adds overloads in namespace `cminpack` that accept ANY callable
   whose signature matches the corresponding cminpack callback (minus the
   leading `void *p`): a lambda, a functor, a std::function, or a plain
   function passed by name or by address. The callable is forwarded through
   the `void *p` slot and invoked by a small template trampoline, so no state
   needs to be global.

   Example:

       #include <cminpackcpp.hpp>
       ...
       std::vector<double> y = ...;                 // captured measurement data
       auto residual = [&](int m, int n, const double *x, double *fvec, int iflag) {
           for (int i = 0; i < m; ++i)
               fvec[i] = model(x, i) - y[i];
           return 0;
       };
       int info = cminpack::lmdif1(residual, m, n, x, fvec, tol, iwa, wa, lwa);

   The wrapper is header-only and precision-agnostic: it uses the same
   __cminpack_real__ / __cminpack_func__ selection as cminpack.h, so defining
   __cminpack_float__ (etc.) before including this header targets that variant.

   The callable must outlive the solver call (which it always does when passed
   as an argument). The wrapper adds no allocation and no virtual dispatch.
*/
#ifndef __CMINPACKCPP_HPP__
#define __CMINPACKCPP_HPP__

#include <cminpack.h>

#if !defined(__cplusplus)
#error "cminpackcpp.hpp is a C++ header; include cminpack.h from C code."
#endif

#include <memory>       /* std::addressof */
#include <type_traits>  /* std::remove_reference, std::is_function */

namespace cminpack {

typedef __cminpack_real__ real;

namespace detail {

/* callable<F> holds what is passed through the `void *p` slot, and `type` is
   what the trampolines cast p back to.

   For objects (lambdas, functors, std::function, function pointers) p points
   to the caller's object itself.
   cv-qualifiers are kept so a const callable is invoked as const. */
template <typename F,
          bool = std::is_function<typename std::remove_reference<F>::type>::value>
struct callable {
    typedef typename std::remove_reference<F>::type type;
    void *ptr;
    explicit callable(typename std::remove_reference<F>::type &fcn)
        : ptr(const_cast<void *>(
              static_cast<const void *>(std::addressof(fcn)))) {}
    void *get() const { return ptr; }
};

/* A plain function passed by name has a function type, and a function's
   address cannot be cast to void * (issue #92). Decay it to a function
   pointer, store that pointer here, and pass the address of the pointer. */
template <typename F>
struct callable<F, true> {
    typedef typename std::remove_reference<F>::type *type;
    type fp;
    explicit callable(type fcn) : fp(fcn) {}
    void *get() { return &fp; }
};

/* trampolines: recover the callable from p and forward the arguments */
template <typename Fn>
int fcn_nn(void *p, int n, const real *x, real *fvec, int iflag)
{
    return (*static_cast<Fn *>(p))(n, x, fvec, iflag);
}

template <typename Fn>
int fcnder_nn(void *p, int n, const real *x, real *fvec, real *fjac,
              int ldfjac, int iflag)
{
    return (*static_cast<Fn *>(p))(n, x, fvec, fjac, ldfjac, iflag);
}

template <typename Fn>
int fcn_mn(void *p, int m, int n, const real *x, real *fvec, int iflag)
{
    return (*static_cast<Fn *>(p))(m, n, x, fvec, iflag);
}

template <typename Fn>
int fcnder_mn(void *p, int m, int n, const real *x, real *fvec, real *fjac,
              int ldfjac, int iflag)
{
    return (*static_cast<Fn *>(p))(m, n, x, fvec, fjac, ldfjac, iflag);
}

template <typename Fn>
int fcnderstr_mn(void *p, int m, int n, const real *x, real *fvec, real *fjrow,
                 int iflag)
{
    return (*static_cast<Fn *>(p))(m, n, x, fvec, fjrow, iflag);
}

} /* namespace detail */

/* --- nonlinear equations: hybrd / hybrd1 (finite-difference Jacobian) --- */

template <typename F>
inline int hybrd1(F &&fcn, int n, real *x, real *fvec, real tol,
                  real *wa, int lwa)
{
    detail::callable<F> c(fcn);
    typedef typename detail::callable<F>::type Fn;
    return ::__cminpack_func__(hybrd1)(
        &detail::fcn_nn<Fn>, c.get(),
        n, x, fvec, tol, wa, lwa);
}

template <typename F>
inline int hybrd(F &&fcn, int n, real *x, real *fvec, real xtol, int maxfev,
                 int ml, int mu, real epsfcn, real *diag, int mode,
                 real factor, int nprint, int *nfev, real *fjac, int ldfjac,
                 real *r, int lr, real *qtf,
                 real *wa1, real *wa2, real *wa3, real *wa4)
{
    detail::callable<F> c(fcn);
    typedef typename detail::callable<F>::type Fn;
    return ::__cminpack_func__(hybrd)(
        &detail::fcn_nn<Fn>, c.get(),
        n, x, fvec, xtol, maxfev, ml, mu, epsfcn, diag, mode, factor, nprint,
        nfev, fjac, ldfjac, r, lr, qtf, wa1, wa2, wa3, wa4);
}

/* --- nonlinear equations: hybrj / hybrj1 (analytic Jacobian) --- */

template <typename F>
inline int hybrj1(F &&fcn, int n, real *x, real *fvec, real *fjac, int ldfjac,
                  real tol, real *wa, int lwa)
{
    detail::callable<F> c(fcn);
    typedef typename detail::callable<F>::type Fn;
    return ::__cminpack_func__(hybrj1)(
        &detail::fcnder_nn<Fn>, c.get(),
        n, x, fvec, fjac, ldfjac, tol, wa, lwa);
}

/* --- least squares, finite-difference Jacobian: lmdif / lmdif1 --- */

template <typename F>
inline int lmdif1(F &&fcn, int m, int n, real *x, real *fvec, real tol,
                  int *iwa, real *wa, int lwa)
{
    detail::callable<F> c(fcn);
    typedef typename detail::callable<F>::type Fn;
    return ::__cminpack_func__(lmdif1)(
        &detail::fcn_mn<Fn>, c.get(),
        m, n, x, fvec, tol, iwa, wa, lwa);
}

template <typename F>
inline int lmdif(F &&fcn, int m, int n, real *x, real *fvec, real ftol,
                 real xtol, real gtol, int maxfev, real epsfcn, real *diag,
                 int mode, real factor, int nprint, int *nfev, real *fjac,
                 int ldfjac, int *ipvt, real *qtf,
                 real *wa1, real *wa2, real *wa3, real *wa4)
{
    detail::callable<F> c(fcn);
    typedef typename detail::callable<F>::type Fn;
    return ::__cminpack_func__(lmdif)(
        &detail::fcn_mn<Fn>, c.get(),
        m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, diag, mode, factor,
        nprint, nfev, fjac, ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4);
}

/* --- least squares, analytic Jacobian: lmder / lmder1 --- */

template <typename F>
inline int lmder1(F &&fcn, int m, int n, real *x, real *fvec, real *fjac,
                  int ldfjac, real tol, int *ipvt, real *wa, int lwa)
{
    detail::callable<F> c(fcn);
    typedef typename detail::callable<F>::type Fn;
    return ::__cminpack_func__(lmder1)(
        &detail::fcnder_mn<Fn>, c.get(),
        m, n, x, fvec, fjac, ldfjac, tol, ipvt, wa, lwa);
}

template <typename F>
inline int lmder(F &&fcn, int m, int n, real *x, real *fvec, real *fjac,
                 int ldfjac, real ftol, real xtol, real gtol, int maxfev,
                 real *diag, int mode, real factor, int nprint, int *nfev,
                 int *njev, int *ipvt, real *qtf,
                 real *wa1, real *wa2, real *wa3, real *wa4)
{
    detail::callable<F> c(fcn);
    typedef typename detail::callable<F>::type Fn;
    return ::__cminpack_func__(lmder)(
        &detail::fcnder_mn<Fn>, c.get(),
        m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, diag, mode,
        factor, nprint, nfev, njev, ipvt, qtf, wa1, wa2, wa3, wa4);
}

/* --- least squares, row-wise analytic Jacobian: lmstr / lmstr1 --- */

template <typename F>
inline int lmstr1(F &&fcn, int m, int n, real *x, real *fvec, real *fjac,
                  int ldfjac, real tol, int *ipvt, real *wa, int lwa)
{
    detail::callable<F> c(fcn);
    typedef typename detail::callable<F>::type Fn;
    return ::__cminpack_func__(lmstr1)(
        &detail::fcnderstr_mn<Fn>, c.get(),
        m, n, x, fvec, fjac, ldfjac, tol, ipvt, wa, lwa);
}

template <typename F>
inline int lmstr(F &&fcn, int m, int n, real *x, real *fvec, real *fjac,
                 int ldfjac, real ftol, real xtol, real gtol, int maxfev,
                 real *diag, int mode, real factor, int nprint, int *nfev,
                 int *njev, int *ipvt, real *qtf,
                 real *wa1, real *wa2, real *wa3, real *wa4)
{
    detail::callable<F> c(fcn);
    typedef typename detail::callable<F>::type Fn;
    return ::__cminpack_func__(lmstr)(
        &detail::fcnderstr_mn<Fn>, c.get(),
        m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, diag, mode,
        factor, nprint, nfev, njev, ipvt, qtf, wa1, wa2, wa3, wa4);
}

} /* namespace cminpack */

#endif /* !__CMINPACKCPP_HPP__ */
