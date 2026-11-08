#!/usr/bin/env python3
"""Unit test for driver_check.py: cross-check exclusions can target chkdrv.

Exclusions are "nprob/dim" tokens extracted from each problem header. The
solver drivers print "problem N dimension D" (or "dimensions M N"), chkdrv
prints "problem N with dimension D"; both must be recognized (issue #101).
"""

import io
import os
import sys
import unittest
from contextlib import redirect_stdout

HERE = os.path.dirname(os.path.abspath(__file__))
sys.dont_write_bytecode = True  # keep the source tree clean
sys.path.insert(0, HERE)
import driver_check  # noqa: E402

REF = os.path.join(HERE, "ref")


class ProbIdTest(unittest.TestCase):
    def test_headers(self):
        cases = [
            ("      problem    1      dimension    2", "1/2"),
            ("      problem    1      dimensions    5   10", "1/5"),
            ("      problem    6      with dimension    9   is  F", "6/9"),
        ]
        for hdr, want in cases:
            self.assertEqual(driver_check._prob_id(hdr), want, hdr)

    def test_chkdrv_exclusion_applies(self):
        ref = os.path.join(REF, "chkdrv.fortran.ref")
        with redirect_stdout(io.StringIO()) as buf:
            rc = driver_check.check(ref, ref, 1e-4, 1e-1, "chkdrv", gate=True,
                                    exclude={"6/9", "14/10"})
        out = buf.getvalue()
        print(out, end="")
        self.assertEqual(rc, 0, out)
        self.assertIn("(2 excluded)", out)


if __name__ == "__main__":
    unittest.main(verbosity=2)
