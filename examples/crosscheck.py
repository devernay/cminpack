#!/usr/bin/env python3
"""Cross-version comparison of the cminpack solvers.

Python port of the former crosscheck.sh. It builds the intensive driver
programs (the difficult More/Garbow/Hillstrom problems) against the pure-C and
f2c cminpack implementations and compares them with driver_check.py to the
committed FORTRAN reference outputs (examples/ref/*.fortran.ref). No Fortran
compiler is needed -- the reference is checked in.

The pass/fail check (DOUBLE precision) is a regression gate against the original
FORTRAN MINPACK: cminpack must converge on every problem FORTRAN converges on.
Where FORTRAN itself does not converge -- the problems pushed from 10x/100x
starting points -- a different result is accepted. pure C is a cleaned-up
rewrite of the f2c output; like different compilers, the two contract FMAs
differently and so take different (equally valid) iteration paths on those
ill-conditioned problems (see README.md, "Numerical differences from FORTRAN
MINPACK"). The long double / float runs, and the pure-C vs f2c comparison, are
printed for information only.

Exit status (so the Makefile treats non-zero as FAIL):
  0  cminpack converges wherever the FORTRAN reference does (at double).
  1  a problem the FORTRAN reference converges on is NOT converged by cminpack.
  2  a build or driver run failed, so the gate could not be evaluated.
"""

import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import driver_check   # noqa: E402  (sibling module)

REF_DIR = os.path.join(HERE, "ref")
EXCLUDE_FILE = os.path.join(HERE, "crosscheck_exclude.txt")
BASECFLAGS = "-O3 -g -Wall"


def load_exclusions():
    """Read the user-maintained 'nprob/dim' exclusion list (crosscheck_exclude.txt)."""
    excl = set()
    try:
        with open(EXCLUDE_FILE, encoding="utf-8") as f:
            for line in f:
                line = line.split("#", 1)[0].strip()
                if line:
                    excl.add(line)
    except OSError:
        pass
    return excl
MAKE_TIMEOUT = 1800   # seconds; a full precision build should be well under this
RUN_TIMEOUT = 300     # seconds per driver; guards against a hung solver

DRIVERS = ["lmddrv", "lmfdrv", "lmsdrv", "hyjdrv", "hybdrv", "chkdrv"]
DATA = {
    "lmddrv": "testdata/ssq.data", "lmfdrv": "testdata/ssq.data",
    "lmsdrv": "testdata/ssq.data", "hyjdrv": "testdata/neq.data",
    "hybdrv": "testdata/neq.data", "chkdrv": "testdata/chkder.data",
}
# label : top-level make target : LIBSUFFIX : precision define : gate-vs-FORTRAN?
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
    exclude = load_exclusions()
    gate_fail = False   # cminpack fails a problem FORTRAN converges on -> exit 1
    build_err = False   # a build or driver run failed -> exit 2
    try:
        for label, top, suf, define, gate in PRECISIONS:
            clib_rel = "../libcminpack%s.a" % suf
            if gate:
                print("=== [%s] cminpack vs FORTRAN reference (regression gate) ===" % label)
            else:
                print("=== [%s] pure C vs f2c (informational: differences on the "
                      "hard problems are EXPECTED -- FMA/compiler codegen) ===" % label)
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
                if gate:
                    ref = os.path.join(REF_DIR, "%s.fortran.ref" % d)
                    try:
                        rc_c = driver_check.check(ref, a, atol=1e-4, rtol=1e-1,
                                                  label="  %s pure-C" % d,
                                                  gate=True, exclude=exclude)
                        rc_f = driver_check.check(ref, b, atol=1e-4, rtol=1e-1,
                                                  label="  %s f2c   " % d,
                                                  gate=True, exclude=exclude)
                    except OSError as e:
                        print("  %s: cannot read reference (%s)" % (d, e))
                        build_err = True
                        continue
                    if rc_c == 1 or rc_f == 1:
                        gate_fail = True
                    elif rc_c == 2 or rc_f == 2:
                        build_err = True
                else:
                    driver_check.check(a, b, atol=1e-4, rtol=1e-1, label="  " + d)
            _make(["clean", "LIBSUFFIX=" + suf])

        if gate_fail:
            print("=== FAIL: cminpack does not converge on a problem the FORTRAN "
                  "reference converges on -- a real regression (see above). ===")
            return 1
        if build_err:
            print("=== INCOMPLETE: a build or driver run failed, so the gate "
                  "could not be fully evaluated. ===")
            return 2
        print("=== PASS: cminpack converges wherever the FORTRAN reference does. "
              "Other differences reported above are expected and harmless. ===")
        return 0
    finally:
        _make(["clean", "LIBSUFFIX="])
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
