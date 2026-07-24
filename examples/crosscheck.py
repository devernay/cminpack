#!/usr/bin/env python3
"""Cross-version comparison of the cminpack solvers.

Python port of the former crosscheck.sh. Keeping all the test tooling in one
language lets this orchestrator import compare.py and driver_check.py directly
instead of shelling out to them, and drops the /bin/sh dependency. It still
drives the Makefile build (it invokes `make`), so it runs in the same Unix-like
environments as `make check`; the portable, cross-platform test route remains
CMake + ctest.

It builds the intensive driver programs (the difficult More/Garbow/Hillstrom
problems) against each implementation -- pure C, f2c, and the original FORTRAN
MINPACK -- and compares their output two ways:

  * pure C (lmddrvc, ...) vs f2c (lmddrv_, ...) at DOUBLE precision MUST be
    byte-identical. This is the one strict invariant; a difference is a real
    bug. The check is a byte compare (filecmp); compare.py is used only to
    explain a failure.

  * every other comparison -- pure C vs f2c at long double/float, and f2c vs
    FORTRAN MINPACK -- is INFORMATIONAL. There the iteration paths legitimately
    differ (see README.md, "Numerical differences from FORTRAN MINPACK"), so
    driver_check.py compares the actual convergence outcome (final residual
    norms) and reports only genuinely different results, not the last-digit
    noise.

Exit status (so the Makefile and CI treat non-zero as FAIL):
  0  the strict invariant holds (pure C == f2c at double); informational
     differences, however large, never change this.
  1  the strict invariant is VIOLATED -- a real bug.
  2  a build/run/IO error prevented the strict check from running.
A confirmed bug (1) takes precedence over a build error (2).
"""

import filecmp
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare        # noqa: E402  (sibling module)
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
    fatal_bug = False        # pure C != f2c at double -> exit 1
    fatal_build_err = False  # could not run the strict check -> exit 2
    try:
        # --- pure C vs f2c, per precision -----------------------------------
        for label, top, suf, define, fatal in PRECISIONS:
            clib_rel = "../libcminpack%s.a" % suf
            if fatal:
                print("=== pure C vs f2c [%s] (must be byte-identical) ===" % label)
            else:
                print("=== pure C vs f2c [%s] (informational: differences here "
                      "are EXPECTED and harmless -- FMA/compiler codegen; the "
                      "problems still converge) ===" % label)
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
                if fatal:
                    # the strict invariant: pure C and f2c MUST be byte-identical
                    if filecmp.cmp(a, b, shallow=False):
                        print("  %s: identical" % d)
                    else:
                        print("  %s: DIFFERENT (unexpected!)" % d)
                        for dd in compare.compare(a, b, rtol=0.0, atol=0.0,
                                                  int_tol=0)[:8]:
                            print("    " + str(dd))
                        fatal_bug = True
                else:
                    # informational: equivalent CONVERGENCE, not bit-exactness
                    driver_check.check(a, b, atol=1e-4, rtol=1e-1, label="  " + d)
            _make(["clean", "LIBSUFFIX=" + suf])

        # --- f2c vs original FORTRAN MINPACK (double only, informational) ----
        if os.path.exists(FLIB_ABS) or shutil.which("gfortran"):
            _make(["-C", "..", "double"])
            if not os.path.exists(FLIB_ABS):
                _make(["-C", "../fortran"])
            if os.path.exists(FLIB_ABS):
                print("=== f2c vs FORTRAN MINPACK [double] (informational: "
                      "differences here are EXPECTED and harmless -- FMA/compiler"
                      " codegen, not enorm; the problems still converge; see "
                      "README.md) ===")
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
            print("=== FAIL: pure-C and f2c differ at double precision -- this "
                  "is a real bug (see the diffs above). ===")
            return 1
        if fatal_build_err:
            print("=== FAIL: the strict pure-C == f2c check could not be run "
                  "(build/run error above). ===")
            return 2
        print("=== PASS: pure-C == f2c at double (the only strict invariant). "
              "All differences reported above are expected and harmless. ===")
        return 0
    finally:
        _make(["clean", "LIBSUFFIX="])
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
