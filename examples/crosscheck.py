#!/usr/bin/env python3
"""Cross-version comparison of the cminpack solvers (INFORMATIONAL).

Python port of the former crosscheck.sh. It builds the intensive driver
programs (the difficult More/Garbow/Hillstrom problems) against each
implementation -- pure C, f2c, and the original FORTRAN MINPACK -- and reports,
with driver_check.py, how their convergence compares.

This is a DIAGNOSTIC, not a pass/fail gate on agreement. pure C (src/*.c) is a
hand-cleaned-up rewrite of the f2c output (src/f2c/*.c); the cleanup regrouped
some expressions, so a compiler may contract "a*b + c" into a fused multiply-add
at different places in the two. On these deliberately-extreme problems that
one-ULP difference can send the runs down different iteration paths and even to
different results -- exactly as different compilers do, and as cminpack does
versus the original FORTRAN MINPACK. Those differences are expected and harmless
(both results are valid); see README.md, "Numerical differences from FORTRAN
MINPACK". The build-failing gates live elsewhere: the standard example tests
(well-conditioned, compared against references) and the driver smoke tests
(each implementation runs to completion without NaN).

Exit status (so the Makefile treats non-zero as a real problem):
  0  the comparison ran; any differences reported above are informational.
  2  a build or driver run failed, so the comparison could not be produced.
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
# label : top-level make target : LIBSUFFIX : precision define
PRECISIONS = [
    ("double",      "double",     "",   ""),
    ("long double", "longdouble", "ld", "-D__cminpack_long_double__"),
    ("float",       "float",      "s",  "-D__cminpack_float__"),
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
    build_err = False  # a build or driver run failed -> exit 2
    try:
        # --- pure C vs f2c, per precision (informational) -------------------
        for label, top, suf, define in PRECISIONS:
            clib_rel = "../libcminpack%s.a" % suf
            print("=== pure C vs f2c [%s] (informational: differences on the "
                  "hard problems are EXPECTED -- FMA/compiler codegen; both "
                  "results are valid) ===" % label)
            if _make(["-C", "..", top]) != 0:
                print("  build of the %s library failed" % label)
                build_err = True
                continue
            cflags = (BASECFLAGS + " " + define).strip()
            broken = build_run(clib_rel, "C_" + suf, "c", suf, cflags, tmp)
            broken |= build_run(clib_rel, "F2C_" + suf, "_", suf, cflags, tmp)
            for d in DRIVERS:
                a = os.path.join(tmp, "C_%s.%s.out" % (suf, d))
                b = os.path.join(tmp, "F2C_%s.%s.out" % (suf, d))
                if d in broken:
                    print("  %s: build/run FAILED" % d)
                    build_err = True
                    continue
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

        if build_err:
            print("=== INCOMPLETE: a build or driver run failed, so the "
                  "cross-check could not be produced. ===")
            return 2
        print("=== Cross-check complete (informational). Differences on the "
              "hard problems are expected; the pass/fail gates are the standard "
              "example tests and the driver smoke tests. ===")
        return 0
    finally:
        _make(["clean", "LIBSUFFIX="])
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
