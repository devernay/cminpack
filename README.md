C/C++ Minpack [![Ubuntu CI](https://github.com/devernay/cminpack/actions/workflows/ubuntu-cminpack-install.yaml/badge.svg)](https://github.com/devernay/cminpack/actions/workflows/ubuntu-cminpack-install.yaml) [![Windows MSVC CI](https://github.com/devernay/cminpack/actions/workflows/windows-visual-studio-cminpack-install.yaml/badge.svg)](https://github.com/devernay/cminpack/actions/workflows/windows-visual-studio-cminpack-install.yaml) [![Windows MSYS2 CI](https://github.com/devernay/cminpack/actions/workflows/msys2-cminpack-install.yaml/badge.svg)](https://github.com/devernay/cminpack/actions/workflows/msys2-cminpack-install.yaml) [![AppVeyor CI](https://ci.appveyor.com/api/projects/status/github/devernay/cminpack?svg=true)](https://ci.appveyor.com/project/devernay/cminpack)
==========

This is a C version of the minpack minimization package.
It has been derived from the fortran code using f2c and
some limited manual editing. Note that you need to link
against libf2c to use this version of minpack. Extern "C"
linkage permits the package routines to be called from C++.
Check ftp://netlib.bell-labs.com/netlib/f2c for the latest
f2c version. For general minpack info and test programs, see
the accompanying readme.txt and http://www.netlib.org/minpack/.

CMake is the standard build system:

    cmake -B build && cmake --build build

and `ctest --test-dir build` runs the tests. See the project home page for
build options (precision variants, BLAS/LAPACK, shared libraries). A plain
`Makefile` is also provided for backward compatibility: type `make` to compile
and `make install` to install in /usr/local, or modify the Makefile to suit
your needs.

Manolis Lourakis -- lourakis at ics forth gr, July 2002
	Institute of Computer Science,
	Foundation for Research and Technology - Hellas
	Heraklion, Crete, Greece

Repackaging by Frederic Devernay -- frederic dot devernay at m4x dot org

The project home page is at http://devernay.github.io/cminpack

C++ bindings
------------

[`include/cminpackcpp.hpp`](include/cminpackcpp.hpp) is a header-only C++
wrapper (issue [#74](https://github.com/devernay/cminpack/issues/74)). The C API takes the user callback as a plain function
pointer plus an opaque `void *p` for user data, which a capturing lambda, a
functor, or a `std::function` cannot satisfy. The wrapper adds overloads in
namespace `cminpack` that accept **any callable** whose signature matches the C
callback minus the leading `void *p`; the callable is forwarded through the
`void *p` slot by a small template trampoline, so no global state is needed. It
adds no allocation and no virtual dispatch, and the callable only has to outlive
the solver call.

The wrapped entry points are `hybrd1`/`hybrd` and `hybrj1` (nonlinear
equations), and `lmdif1`/`lmdif`, `lmder1`/`lmder` and `lmstr1`/`lmstr` (least
squares). Every argument after the callable is identical to the C function — you
just drop the `(fcn, p)` pair and pass your callable first:

```cpp
#include <cminpackcpp.hpp>
#include <vector>

std::vector<double> y = /* measured data */;

// signature = the C callback minus the leading void *p
auto residual = [&](int m, int n, const double *x, double *fvec, int iflag) {
    for (int i = 0; i < m; ++i)
        fvec[i] = model(x, i) - y[i];   // captures y, no global needed
    return 0;
};

int info = cminpack::lmdif1(residual, m, n, x, fvec, tol, iwa, wa, lwa);
```

A **class member function** works too, as long as it is wrapped in a callable
that supplies the instance — a lambda capturing `this` (or the object) is the
simplest, and `std::bind`/`std::function` also work. A bare pointer-to-member
(`&Class::residual`) is not itself callable, so it must be wrapped:

```cpp
struct LineModel {
    std::vector<double> y;
    int residual(int m, int n, const double *x, double *fvec, int iflag) const {
        for (int i = 0; i < m; ++i) fvec[i] = model(x, i) - y[i];
        return 0;
    }
};

LineModel obj = /* ... */;
auto cb = [&obj](int m, int n, const double *x, double *fvec, int iflag) {
    return obj.residual(m, n, x, fvec, iflag);      // supplies the instance
};
int info = cminpack::lmdif1(cb, m, n, x, fvec, tol, iwa, wa, lwa);
// std::bind is equivalent:
//   auto cb = std::bind(&LineModel::residual, &obj, _1, _2, _3, _4, _5);
```

The header is precision-agnostic: like `cminpack.h` it selects the real type via
`__cminpack_real__`, so defining `__cminpack_float__` (or the long-double macro)
before including it targets that variant — link against the matching cminpack
library. A worked example exercising a capturing lambda, a stateful functor and
a `std::function` is in [`examples/tcppwrap.cpp`](examples/tcppwrap.cpp).

Testing
-------

Two test routes are available:

- **CMake / CTest** (recommended): after building, run `ctest --test-dir build`.
  This runs the standard example tests (compared against `examples/ref/*.ref`
  with the dependency-free C tool `cmpfiles`), the self-checking regression
  tests, the intensive driver programs as smoke tests (run to completion, no
  NaN), and — when Python 3 is available — the FORTRAN-reference cross-check
  (see below). Set `-DCMINPACK_CROSSCHECK=OFF` to disable that cross-check.
- **Makefile**: `make check` runs the standard tests for the double, long
  double and float builds, then `make crosscheck`. The double standard tests
  are strict (a failure fails the build); the long double and float tests are
  informational.

The cross-check (`examples/crosscheck.py`, and the CMake `crosscheck_*` tests)
is a **regression gate against the original FORTRAN MINPACK**: on the intensive
driver problems, cminpack must converge on every problem FORTRAN converges on.
It compares against committed FORTRAN reference outputs
(`examples/ref/*.fortran.ref`), so **no Fortran compiler is needed**. Where
FORTRAN itself does not converge (the problems pushed from `10*x0`/`100*x0`
starting points), a different result is accepted — pure-C and f2c, like
different compilers, take different (equally valid) iteration paths there, for
the reasons in the next section. Iteration-count and last-digit differences
never fail the build.

**Python 3 is optional.** It is not needed to build cminpack, nor for the
standard or smoke tests (which use `cmpfiles`); only the cross-check uses it
(via `examples/crosscheck.py` and `examples/driver_check.py`). When Python 3 is
missing, both build systems skip only that cross-check — printing a disclaimer
that coverage is incomplete — and run everything else.

Numerical differences from FORTRAN MINPACK
------------------------------------------

On the difficult test problems (the Moré/Garbow/Hillstrom functions exercised
by the intensive driver programs `examples/*drv*`), cminpack and the original
FORTRAN MINPACK can report different iteration, function- and
Jacobian-evaluation counts, and last-digit differences in the results. **This
is expected and harmless: the problems still converge to equally valid
solutions.** Test messages reporting such differences are normal, not a sign of
a broken build.

The difficult test problems themselves are from Moré, Garbow & Hillstrom,
*Testing Unconstrained Optimization Software*, ACM TOMS 7(1):17-41 (1981)
([doi:10.1145/355934.355936](https://doi.org/10.1145/355934.355936); free
open-access technical report ANL-AMD-TM-324:
<https://www.osti.gov/servlets/purl/6650344>), and the companion Algorithm 566,
ACM TOMS 7(1):136-140 (1981)
([doi:10.1145/355934.355943](https://doi.org/10.1145/355934.355943)). They are
defined in `examples/ssqfcn.c` (least squares) and `examples/vecfcn.c` (systems
of equations); the driver programs run each one from progressively harder
starting points (`x0`, `10*x0`, `100*x0`), so the hardest problems (e.g.
problem 10, the Meyer function) are deliberately pushed to non-convergence.

The differences do *not* come from the algorithm:

- `dpmpar` returns identical machine constants on both sides, so the convergence
  tolerance is the same.
- cminpack's own pure C (`src/`) is a cleaned-up rewrite of the f2c output
  (`src/f2c/`); the two agree on well-conditioned problems, but the cleanup
  regrouped some expressions, so the compiler can contract FMAs differently in
  each — meaning even these two can reach different (equally valid) results on
  the hardest driver problems.
- The dominant cause is floating-point contraction: `gcc` and `gfortran` fuse
  `a*b + c` into a fused multiply-add (FMA) at different places, so intermediate
  values differ by one unit in the last place (ULP). On ill-conditioned problems
  a single ULP early in the iteration can flip a trust-region accept/reject
  decision, and the two runs then follow different paths to (equally valid)
  results. Compiling both sides with `-ffp-contract=off` removes most of the
  divergence.

(The `enorm` scaling constants `rdwarf`/`rgiant`, which cminpack tunes to the
IEEE range rather than MINPACK's 1980 values, are a separate and minor effect:
they only change the norm of out-of-range vectors in the last bit. See
[`src/enorm.c`](src/enorm.c) for the full explanation.)

Because of this, `examples/crosscheck.py` and the CMake `crosscheck_*` tests are
a **regression gate against the original FORTRAN MINPACK**: on the intensive
driver problems, cminpack must converge on every problem FORTRAN converges on
(compared against the committed `examples/ref/*.fortran.ref`, so no Fortran
compiler is needed). Following the original MINPACK drivers, a problem is gated
only where FORTRAN converged (exit parameter `info < 5`); where FORTRAN itself
does not converge (`info >= 5`) a different result is accepted. On a gated
problem the build fails if either its own solver reports non-convergence
(`info >= 5`) or FORTRAN drove the residual to ~0 but this build did not (which
catches a genuine divergence the exit code alone would miss). Iteration-count
and last-digit differences are reported for information only and never fail the
build. A handful of deliberately-extreme
problems have compiler-dependent, coin-flip convergence (FORTRAN itself
converges on them only at some optimization levels); those are listed, with the
rationale, in [`examples/crosscheck_exclude.txt`](examples/crosscheck_exclude.txt)
and skipped by the gate. BLAS/LAPACK builds take a different numeric path (a
BLAS `dnrm2` in `enorm`, a LAPACK QR factorization) that can flip a different,
implementation-dependent set, so they get an additional list,
[`examples/crosscheck_exclude_blas.txt`](examples/crosscheck_exclude_blas.txt),
applied only when `USE_BLAS`/`USE_LAPACK` is on.

History
------

* version 1.3.14 (25/07/2026):
  - Fix spurious `make check` failures on the intensive driver tests ([#78](https://github.com/devernay/cminpack/issues/78)): the
    non-portable driver-vs-reference text comparison was dropped. Driver
    validation is now a regression gate against the original FORTRAN MINPACK —
    cminpack (pure C and f2c) must converge on every driver problem the committed
    reference (`examples/ref/*.fortran.ref`) converges on, so no Fortran compiler
    is needed. Iteration-count and last-digit differences (compiler FMA) never
    fail the build.
  - Add user-maintained exclusion lists for the few deliberately-extreme
    coin-flip problems (FORTRAN itself converges on them only at some
    optimization levels): `examples/crosscheck_exclude.txt`, plus
    `examples/crosscheck_exclude_blas.txt` for `USE_BLAS`/`USE_LAPACK` builds.
    Each gate failure prints the problem as `nprob/dim`, the token to add.
  - Replace `crosscheck.sh` with `crosscheck.py`; add `driver_check.py`
    (`--reference-gate`) and `crosscheck_matrix.sh` (an optimization-level
    diagnostic). Wire the cross-check into CMake/CTest (`CMINPACK_CROSSCHECK`);
    Python 3 is optional (only the cross-check needs it).
  - Fix a compile error in `examples/tenorm_.c` (pass `dpmpar_` its index by
    pointer, the FORTRAN/f2c calling convention).
  - Document the numerical differences from FORTRAN MINPACK and correct the
    `enorm.c` comments (the divergence is FMA, not the rdwarf/rgiant constants).

* version 1.3.13 (16/07/2026):
  - Fix a division by zero in `covar1` when the Jacobian rank equals the number
    of residuals (m == rank, e.g. square full-rank problems), which produced
    Inf/NaN throughout the covariance matrix
  - Make the banded finite-difference branch of `fdjac1` call the user function
    with iflag=2, like the dense branch, `fdjac2` and the FORTRAN version
  - Fix the `USE_BLAS` Newton correction in `lmpar` to use all n components (it
    was truncated to the original Jacobian rank, giving a wrong step for
    rank-deficient problems)
  - Fix the jpvt memset size in the `USE_LAPACK` `qrfac`, and make the
    work-array size checks in `hybrd`, `hybrj`, `hybrd1`, `hybrj1` and `lmdif1`
    overflow-safe (also in the f2c versions)
  - Guard a harmless transient infinity in `dogleg` for rank-deficient
    Jacobians (also in the f2c version)
  - Add a `USE_LAPACK` CMake option to build the LAPACK-based QR factorization
    (`qrfac`), which was previously only reachable through the `Makefile`
  - Exercise the BLAS and LAPACK builds in CI and run the full test suite
    (including `tlmdifc`) for them; add the intensive driver programs
    (`lmddrvc`, `lmfdrvc`, `lmsdrvc`, `hyjdrvc`, `hybdrvc`, `chkdrvc`) as smoke
    tests and run the FORTRAN example tests in CI

* version 1.3.12 (16/07/2026):
  - Fix non-convergence/NaN in `lmder`, `lmdif` and `lmstr` on problems whose
    solution is the zero vector, by guarding a 0/0 division in `lmpar` [#76](https://github.com/devernay/cminpack/issues/76)
  - Compare test output using a numeric tolerance instead of an exact text
    match, so the tests pass across compilers and math libraries [#37](https://github.com/devernay/cminpack/issues/37) [#77](https://github.com/devernay/cminpack/issues/77)
  - Fix Windows linking: correct the DLL export/import macros and document
    that static-library users define `CMINPACK_NO_DLL` [#18](https://github.com/devernay/cminpack/issues/18)
  - Detect and link the CBLAS interface when `USE_BLAS` is enabled [#12](https://github.com/devernay/cminpack/issues/12)
  - Add `cminpackcpp.hpp`, a header-only C++ wrapper so the solvers accept
    lambdas, functors and `std::function` [#74](https://github.com/devernay/cminpack/issues/74)
  - Make CMake the standard build system (the `Makefile` is kept for backward
    compatibility) and remove the unmaintained Xcode, Visual Studio and Eclipse
    project files (`cminpack.xcodeproj`, `cminpack*.vcproj`/`.vcxproj`,
    `cminpack.sln`, `.cproject`, `.project`)
  - Move CI to GitHub Actions (keeping the CMake-based AppVeyor build) and
    remove the obsolete Travis CI, Coveralls and Coverity configuration and
    badges

* version 1.3.11 (13/09/2024):
  - Bump installed version number to 1.3.11 [#75](https://github.com/devernay/cminpack/issues/75)

* version 1.3.10 (11/09/2024):
  - Disable BLAS by default [#66](https://github.com/devernay/cminpack/issues/66) [#70](https://github.com/devernay/cminpack/issues/70)
  - Fix BLAS usage (broken by [#58](https://github.com/devernay/cminpack/issues/58)) [#68](https://github.com/devernay/cminpack/issues/68)
  - Fix testing on Windows [#63](https://github.com/devernay/cminpack/issues/63)
  - Fix pkg-config files [#71](https://github.com/devernay/cminpack/issues/71)
  - Bump minimum CMake version [#69](https://github.com/devernay/cminpack/issues/69)

* version 1.3.9 (28/05/2024):
  - CMake portability fixes [#50](https://github.com/devernay/cminpack/issues/50) [#53](https://github.com/devernay/cminpack/issues/53) [#56](https://github.com/devernay/cminpack/issues/56) [#57](https://github.com/devernay/cminpack/issues/57) [#58](https://github.com/devernay/cminpack/issues/58) [#62](https://github.com/devernay/cminpack/issues/62)
  - MKL-related fixes [#51](https://github.com/devernay/cminpack/issues/51) [#52](https://github.com/devernay/cminpack/issues/52)
  - Support more CI build configurations [#59](https://github.com/devernay/cminpack/issues/59) [#61](https://github.com/devernay/cminpack/issues/61)

* version 1.3.8 (09/02/2021):
  - CMake now builds by default the single-, double-, and extended-precision versions [#45](https://github.com/devernay/cminpack/issues/45) [#48](https://github.com/devernay/cminpack/issues/48)
  - Avoid promoting to doubles in all operations for the single-precision version [#47](https://github.com/devernay/cminpack/issues/47)

* version 1.3.7 (09/12/2020):
  - Makefile cleanups [#11](https://github.com/devernay/cminpack/issues/11)
  - Cmake-related fixes [#20](https://github.com/devernay/cminpack/issues/20) [#21](https://github.com/devernay/cminpack/issues/21) [#23](https://github.com/devernay/cminpack/issues/23) [#27](https://github.com/devernay/cminpack/issues/27) [#28](https://github.com/devernay/cminpack/issues/28)
  - Add Appveyor CI [#24](https://github.com/devernay/cminpack/issues/24)
  - Add support for single-precision CBLAS and LAPACK [#40](https://github.com/devernay/cminpack/issues/40)

* version 1.3.6 (24/02/2017):
  - Fix FreeBSD build [#6](https://github.com/devernay/cminpack/issues/6)
  - CMake: install CMinpackConfig.cmake rather than FindCMinpack.cmake [#8](https://github.com/devernay/cminpack/issues/8)
  - CMake: add option USE_BLAS to compile with blas [#9](https://github.com/devernay/cminpack/issues/9)

* version 1.3.5 (28/05/2016):
  - Add support for compiling a long double version (Makefile only).
  - CMake: static libraries now have the suffix `_s`.

* version 1.3.4 (28/05/2014):
  - Add FindCMinpack.cmake cmake module. If you use the cmake install,
    finding CMinpack from your `CMakeLists.txt` is as easy as
    `find_package(CMinpack)`.

* version 1.3.3 (04/02/2014):
  - Add documentation and examples abouts how to add box constraints to the variables.
  - continuous integration https://travis-ci.org/devernay/cminpack

* version 1.3.2 (27/10/2013):
  - Minor change in the CMake build: also set `SOVERSION`.

* version 1.3.1 (02/10/2013):
  - Fix CUDA examples compilation, and remove non-free files.

* version 1.3.0 (09/06/2012):
  - Optionally use LAPACK and CBLAS in lmpar, qrfac, and qrsolv. Added
    `make lapack` to build the LAPACK-based cminpack and "make
    checklapack" to test it (results of the test may depend on the
    underlying LAPACK and BLAS implementations).
    On 64-bits architectures, the preprocessor symbol `__LP64__` must be
    defined (see [`cminpackP.h`](cminpackP.h)) if the LAPACK library uses the LP64
    interface (i.e. 32-bits integer, vhereas the ILP interface uses 64
    bits integers).

* version 1.2.2 (16/05/2012):
  - Update Makefiles and documentation (see "Using CMinpack" above) for
    easier building and testing.

* version 1.2.1 (15/05/2012):
  - The library can now be built as double, float or half
    versions. Standard tests in the "examples" directory can now be
    lauched using `make check` (to run common tests, including against
    the float version), `make checkhalf` (to test the half version) and
    `make checkfail` (to run all the tests, even those that fail).

* version 1.2.0 (14/05/2012):
  - Added original FORTRAN sources for better testing (type
    `make -C fortran`, then `make -C examples` and follow the
    instructions). Added driver tests `lmsdrv`, `chkdrv`, `hyjdrv`,
    `hybdrv`. `make -C examples alltest` will run all
    possible test combinations (make sure you have gfortran installed).

* version 1.1.5 (04/05/2012):
  - cminpack now works in CUDA, thanks to Jordi Bataller Mascarell, type
    `make -C cuda` (be careful, though: this is a
    straightforward port from C, and each problem is solved using a
    single thread). cminpack can now also be compiled with
    single-precision floating point computation (define
    `__cminpack_real__` to float when compiling and using the
   library). Fix cmake support for `CMINPACK_LIB_INSTALL_DIR`. Update the
   reference files for tests.

* version 1.1.4 (30/10/2011):
  - Translated all the Levenberg-Marquardt code (`lmder`, `lmdif`, `lmstr`,
    `lmder1`, `lmdif1`, `lmstr1`, `lmpar`, `qrfac`, `qrsolv`, `fdjac2`, `chkder`) to use
    C-style indices.

* version 1.1.3 (16/03/2011):
  - Minor fix: Change non-standard `strnstr()` to `strstr()` in
    [`genf77tests.c`](examples/genf77tests.c).

* version 1.1.2 (07/01/2011):
   - Fix Windows DLL building (David Graeff) and document covar in
     cminpack.h.

* version 1.1.1 (04/12/2010):
  - Complete rewrite of the C functions (without trailing underscore in
    the function name). Using the original FORTRAN code, the original
    algorithms structure was recovered, and many goto's were converted
    to if...then...else. The code should now be both more readable and
    easier to optimize, both for humans and for compilers. Added lmddrv
    and lmfdrv test drivers, which test a lot of difficult functions
    (these functions are explained in Testing Unconstrained Optimization
    Software by Moré et al.). Also added the pkg-config files to the
    cmake build, as well as an "uninstall" target, contributed by
    Geoffrey Biggs.

* version 1.0.4 (18/10/2010):
  - Support for shared library building using CMake, thanks to Goeffrey
    Biggs and Radu Bogdan Rusu from Willow Garage. Shared libraries can be
    enabled using cmake options, as in:
    `cmake -DUSE_FPIC=ON -DBUILD_SHARED_LIBS=ON -DBUILD_EXAMPLES=OFF path_to_sources`

* version 1.0.3 (18/03/2010):
  - Added CMake support.
  - XCode build is now Universal.
  - Added `tfdjac2_` and `tfdjac2c` examples, which test the accuracy of a
    finite-differences approximation of the Jacobian.
  - Bug fix in `tlmstr1` (signaled by Thomas Capricelli).

* version 1.0.2 (27/02/2009):
  - Added Xcode and Visual Studio project files

* version 1.0.1 (17/12/2007):
  - bug fix in `covar()` and `covar_()`, the computation of tolr caused a
    segfault (signaled by Timo Hartmann).

* version 1.0.0 (24/04/2007):
  - Added FORTRAN and C examples
  - Added documentation from Debian man pages
  - Wrote pure C version
  - Added `covar()` and `covar_()`, and use it in `tlmdef`/`tlmdif`
