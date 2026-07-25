#!/usr/bin/env python3
"""Cross-version comparison of the cminpack solvers.

Python port of the former crosscheck.sh. Keeping all the test tooling in one
language lets this orchestrator import driver_check.py directly instead of
shelling out, and drops the /bin/sh dependency. It still drives the Makefile
build (it invokes `make`), so it runs in the same Unix-like environments as
`make check`; the portable, cross-platform test route remains CMake + ctest.

It builds the intensive driver programs (the difficult More/Garbow/Hillstrom
problems) against each implementation -- pure C, f2c, and the original FORTRAN
MINPACK -- and checks, with driver_check.py, that they reach *equivalent
results*: every problem's final L2 residual norm agrees. It does NOT require
byte-identical output. pure-C (src/*.c) and f2c (src/f2c/*.c) are two
independent source trees, so a compiler may contract "a*b + c" into a fused
multiply-add at different places in each; on the ill-conditioned problems that
one-ULP difference sends the runs down different (but equally convergent)
iteration paths. The same is true across compilers and versus the original
FORTRAN MINPACK (see README.md, "Numerical differences from FORTRAN MINPACK").

  * pure C vs f2c at DOUBLE precision is the one strict check: both must
    converge to the same solution on every problem. A genuine divergence
    (one converges, the other does not) is a real bug and makes this exit 1.
  * pure C vs f2c at long double/float, and f2c vs FORTRAN MINPACK, are
    INFORMATIONAL: differences are reported but never fail the run.

Exit status (so the Makefile and CI treat non-zero as FAIL):
  0  the strict check passed (pure C and f2c converge equivalently at double);
     informational differences never change this.
  1  the strict check FAILED -- pure C and f2c reached a materially different
     result at double precision.
  2  a build/run/IO error prevented the strict check from running.
A confirmed failure (1) takes precedence over a build error (2).
"""

import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import driver_check   # noqa: E402  (sibling module)

BASECFLAGS = "-O3 -g -Wall"
FLIB_REL = "../fortran/libminpack.a"
FLIB_ABS = os.path.join(HERE, "..", "fortran", "libminpack.a")
MAKE_TIMEOUT = 1800   # seconds; a full precision build should be well under this
RUN_TIMEOUT = 300     # seconds per driver; guards against a hung solver

DRIVERS = ["lmddrv", "lmfdrv", "lmsdrv", "hyjdrv", "hybdrv", "chkdrv"]
DATA = {
    "lmddrv": "testdata/ssq.data", "lmfdrv": "testdata/ssq.data",
    "lmsdrv": "testdata/ssq.data", "hyjdrv": "testdata/neq.data",
    "hybdrv": "testdata/neq.data", "chkdrv": "testdata/chkder.data",
}
# label : top-level make target : LIBSUFFIX : precision define : fatal?
PRECISIONS = [
    ("double",      "double",     "",   "",                           True),
    ("long double", "longdouble", "ld", "-D__cminpack_long_double__",  False),
    ("float",       "float",      "s",  "-D__cminpack_float__",        False),
]


def _make(args, quiet=True):
    """Run `make ARGS` in the examples directory; return its exit code (a
    non-zero code, including a timeout, means the build step failed)."""
    kw = {"cwd": HERE}
    if quiet:
        kw["stdout"] = subprocess.DEVNULL
        kw["stderr"] = subprocess.DEVNULL
    try:
        return subprocess.run(["make"] + args, timeout=MAKE_TIMEOUT, **kw).returncode
    except subprocess.TimeoutExpired:
        return 124


def build_run(lib_rel, tag, progsuffix, libsuffix, cflags, outdir):
    """Build the driver set for one configuration and run each driver on its
    input data, writing <tag>.<driver>.out into outdir. Return the set of
    drivers whose executable was missing, crashed, or timed out."""
    _make(["clean", "LIBSUFFIX=" + libsuffix])
    progs = [libsuffix + d + progsuffix for d in DRIVERS]
    _make(progs + ["LIBSUFFIX=" + libsuffix, "CFLAGS=" + cflags,
                   "MINPACK=" + lib_rel])
    failed = set()
    for d in DRIVERS:
        prog = os.path.join(HERE, libsuffix + d + progsuffix)
        out = os.path.join(outdir, "%s.%s.out" % (tag, d))
        if not (os.path.exists(prog) and os.access(prog, os.X_OK)):
            failed.add(d)
            continue
        try:
            with open(os.path.join(HERE, DATA[d]), "rb") as fin, \
                    open(out, "wb") as fout:
                rc = subprocess.run([prog], cwd=HERE, stdin=fin, stdout=fout,
                                    stderr=subprocess.STDOUT,
                                    timeout=RUN_TIMEOUT).returncode
            if rc != 0:
                failed.add(d)
        except subprocess.TimeoutExpired:
            failed.add(d)
    return failed


def main():
    tmp = tempfile.mkdtemp(prefix="cross.", dir=HERE)
    fatal_bug = False        # pure C != f2c convergence at double -> exit 1
    fatal_build_err = False  # could not run the strict check -> exit 2
    try:
        # --- pure C vs f2c, per precision -----------------------------------
        for label, top, suf, define, fatal in PRECISIONS:
            clib_rel = "../libcminpack%s.a" % suf
            if fatal:
                print("=== pure C vs f2c [%s] (must converge equivalently) ===" % label)
            else:
                print("=== pure C vs f2c [%s] (informational: iteration paths "
                      "may differ -- FMA/compiler codegen; both still converge) ==="
                      % label)
            if _make(["-C", "..", top]) != 0:
                print("  build of the %s library failed" % label)
                if fatal:
                    fatal_build_err = True
                continue
            cflags = (BASECFLAGS + " " + define).strip()
            broken = build_run(clib_rel, "C_" + suf, "c", suf, cflags, tmp)
            broken |= build_run(clib_rel, "F2C_" + suf, "_", suf, cflags, tmp)
            for d in DRIVERS:
                a = os.path.join(tmp, "C_%s.%s.out" % (suf, d))
                b = os.path.join(tmp, "F2C_%s.%s.out" % (suf, d))
                if d in broken:
                    print("  %s: build/run FAILED" % d)
                    if fatal:
                        fatal_build_err = True
                    continue
                # Compare CONVERGENCE (final residual norms), not bytes:
                # pure C and f2c are independent sources and need not be
                # bit-identical, but they must reach equivalent solutions.
                rc = driver_check.check(a, b, atol=1e-4, rtol=1e-1, label="  " + d)
                if fatal and rc != 0:
                    fatal_bug = True
            _make(["clean", "LIBSUFFIX=" + suf])

        # --- f2c vs original FORTRAN MINPACK (double only, informational) ----
        if os.path.exists(FLIB_ABS) or shutil.which("gfortran"):
            _make(["-C", "..", "double"])
            if not os.path.exists(FLIB_ABS):
                _make(["-C", "../fortran"])
            if os.path.exists(FLIB_ABS):
                print("=== f2c vs FORTRAN MINPACK [double] (informational: "
                      "differences here are EXPECTED and harmless -- FMA/compiler"
                      " codegen; the problems still converge; see README.md) ===")
                build_run("../libcminpack.a", "F2Cd", "_", "", BASECFLAGS, tmp)
                build_run(FLIB_REL, "FORT", "_", "", BASECFLAGS, tmp)
                for d in DRIVERS:
                    a = os.path.join(tmp, "F2Cd.%s.out" % d)
                    b = os.path.join(tmp, "FORT.%s.out" % d)
                    if os.path.exists(a) and os.path.exists(b):
                        driver_check.check(a, b, atol=1e-4, rtol=1e-1,
                                           label="  " + d)
                    else:
                        print("  %s: build/run failed (informational, skipped)"
                              % d)
        else:
            print("=== f2c vs FORTRAN skipped (no Fortran compiler) ===")

        if fatal_bug:
            print("=== FAIL: pure-C and f2c reached a materially different "
                  "result at double precision -- this is a real bug (see "
                  "above). ===")
            return 1
        if fatal_build_err:
            print("=== FAIL: the strict pure-C vs f2c check could not be run "
                  "(build/run error above). ===")
            return 2
        print("=== PASS: pure-C and f2c converge equivalently at double (the "
              "strict check). All other differences reported above are expected "
              "and harmless. ===")
        return 0
    finally:
        _make(["clean", "LIBSUFFIX="])
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
