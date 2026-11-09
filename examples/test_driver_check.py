#!/usr/bin/env python3
"""Tests of driver_check.py on FORTRAN references with one line changed."""

import io
import os
import sys
import tempfile
from contextlib import redirect_stdout

HERE = os.path.dirname(os.path.abspath(__file__))
sys.dont_write_bytecode = True  # keep the source tree clean
sys.path.insert(0, HERE)
import driver_check  # noqa: E402


def ref(name):
    with open(os.path.join(HERE, "ref", name + ".fortran.ref")) as f:
        return f.read()


def gate(name, out):
    """Exit status of the reference gate of the output text out vs ref name."""
    with tempfile.NamedTemporaryFile("w", suffix=".out", delete=False) as f:
        f.write(out)
    try:
        with redirect_stdout(io.StringIO()):
            return driver_check.check(
                os.path.join(HERE, "ref", name + ".fortran.ref"), f.name,
                1e-4, 0.1, "", gate=True)
    finally:
        os.unlink(f.name)


def test_error_vector_read_fully():  # issue 100
    hdr, blk = driver_check._split_problems(ref("chkdrv").splitlines(True))[5]
    assert "dimension    9" in hdr
    n = len(driver_check._error_vectors(blk))
    assert n == 9, "read %d of the 9 error-vector elements" % n

def test_chkdrv_error_vector_changed():  # issue 99
    out = ref("chkdrv")
    v = out.split("error vector")[1].split()[0]  # problem 1, first element
    assert gate("chkdrv", out.replace(v, "0.5", 1)) == 1, \
        "changed error vector not detected"


if __name__ == "__main__":
    for name, test in list(globals().items()):
        if name.startswith("test_"):
            test()
            print("ok", name)
