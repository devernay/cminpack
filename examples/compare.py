#!/usr/bin/env python3
"""Tolerant numeric comparison for cminpack test output.

Compares two text files token by token and reports only the differences that
actually matter. It ignores last-digit and formatting noise -- the kind that
legitimately varies between compilers, C runtimes and math libraries (see
issues #37 and #77). It is a higher-level companion to the in-tree C tool
``cmpfiles`` (which stays dependency-free): use it when you want a clear,
categorised report, especially for cross-version comparison (FORTRAN vs f2c vs
pure C), where the integer counts (iterations, nfev, njev) must be identical and
only the printed floating-point values may differ in the last digit.

Tokens are classified and compared as follows:

* integers (no '.', no exponent)  -- compared EXACTLY by default. These carry
  the structurally meaningful counts (iterations, nfev, njev, exit codes), so a
  mismatch is a real ("structural") difference. Pass --int-tol to relax.
* floating-point numbers          -- compared with a relative tolerance
  (--rtol) and an absolute tolerance (--atol). Two values are equal if they
  differ by at most --atol (which also absorbs near-zero noise); otherwise a
  relative difference above --rtol is a real ("numeric") difference.
* everything else                 -- compared as text (whitespace-insensitive).

Exit status is 0 when the files are equal up to tolerance, 1 when a difference
that matters is found, and 2 on usage/IO errors. With --quiet only the exit
status is produced.
"""

import argparse
import math
import re
import sys

# A token is a run of non-whitespace characters.
_TOKEN_RE = re.compile(r"\S+")
# An integer token: optional sign, digits only (no '.', no exponent).
_INT_RE = re.compile(r"[+-]?\d+$")
# A floating-point token (also matches Fortran 'D' exponents like 1.0D-3).
_FLOAT_RE = re.compile(r"[+-]?(\d+\.?\d*|\.\d+)([eEdD][+-]?\d+)?$")


def _as_float(tok):
    """Parse a numeric token to float, accepting Fortran 'D' exponents.

    Returns None if the token is not numeric."""
    if _FLOAT_RE.match(tok):
        try:
            return float(tok.replace("D", "e").replace("d", "e"))
        except ValueError:
            return None
    return None


class Diff:
    __slots__ = ("line", "col", "kind", "a", "b", "detail")

    def __init__(self, line, col, kind, a, b, detail):
        self.line, self.col, self.kind = line, col, kind
        self.a, self.b, self.detail = a, b, detail

    def __str__(self):
        return (f"line {self.line} token {self.col}: {self.kind}: "
                f"{self.a!r} vs {self.b!r} ({self.detail})")


def compare(path_a, path_b, rtol, atol, int_tol):
    """Compare two files, returning a list of Diff objects that matter."""
    with open(path_a, encoding="utf-8", errors="replace") as fa:
        lines_a = fa.readlines()
    with open(path_b, encoding="utf-8", errors="replace") as fb:
        lines_b = fb.readlines()

    diffs = []
    if len(lines_a) != len(lines_b):
        diffs.append(Diff(0, 0, "structural", len(lines_a), len(lines_b),
                          "different number of lines"))

    for ln, (la, lb) in enumerate(zip(lines_a, lines_b), start=1):
        ta = _TOKEN_RE.findall(la)
        tb = _TOKEN_RE.findall(lb)
        if len(ta) != len(tb):
            diffs.append(Diff(ln, 0, "structural", len(ta), len(tb),
                              "different number of tokens on line"))
            continue
        for col, (a, b) in enumerate(zip(ta, tb), start=1):
            if a == b:
                continue
            a_is_int = bool(_INT_RE.match(a))
            b_is_int = bool(_INT_RE.match(b))
            if a_is_int and b_is_int:
                if abs(int(a) - int(b)) > int_tol:
                    diffs.append(Diff(ln, col, "structural", a, b,
                                      f"integer differs by {int(a) - int(b)}"))
                continue
            fa_ = _as_float(a)
            fb_ = _as_float(b)
            if fa_ is not None and fb_ is not None:
                # non-finite (e.g. an overflowing 1e999): compare directly,
                # never via a ratio -- inf/inf is NaN and would slip through as
                # "equal".
                if not (math.isfinite(fa_) and math.isfinite(fb_)):
                    if fa_ == fb_:
                        continue
                    diffs.append(Diff(ln, col, "numeric", a, b,
                                      "non-finite value differs"))
                    continue
                # both near zero -> noise, treat as equal
                if abs(fa_) <= atol and abs(fb_) <= atol:
                    continue
                denom = max(abs(fa_), abs(fb_))
                absdiff = abs(fa_ - fb_)
                if absdiff <= atol:
                    continue
                reldiff = absdiff / denom if denom else 0.0
                if reldiff > rtol:
                    diffs.append(Diff(ln, col, "numeric", a, b,
                                      f"reldiff={reldiff:.3g} > rtol={rtol:g}"))
                continue
            # non-numeric mismatch
            diffs.append(Diff(ln, col, "text", a, b, "text token differs"))
    return diffs


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("file_a")
    p.add_argument("file_b")
    p.add_argument("--rtol", type=float, default=1e-2,
                   help="relative tolerance for floating-point tokens (default 1e-2)")
    p.add_argument("--atol", type=float, default=1e-4,
                   help="absolute tolerance for floats: values within --atol are "
                        "equal, which also absorbs near-zero noise (default 1e-4)")
    p.add_argument("--int-tol", type=int, default=0,
                   help="max allowed difference between integer tokens "
                        "(default 0 = must match exactly; these are the "
                        "iteration/nfev/njev counts)")
    p.add_argument("-q", "--quiet", action="store_true",
                   help="print nothing, only set the exit status")
    p.add_argument("--max-report", type=int, default=40,
                   help="max differences to print (default 40)")
    args = p.parse_args(argv)

    try:
        diffs = compare(args.file_a, args.file_b, args.rtol, args.atol,
                        args.int_tol)
    except OSError as e:
        print(f"compare.py: {e}", file=sys.stderr)
        return 2

    if not diffs:
        if not args.quiet:
            print(f"OK: '{args.file_a}' and '{args.file_b}' match "
                  f"(rtol={args.rtol:g} atol={args.atol:g} int-tol={args.int_tol})")
        return 0

    if not args.quiet:
        structural = sum(1 for d in diffs if d.kind == "structural")
        numeric = sum(1 for d in diffs if d.kind == "numeric")
        text = sum(1 for d in diffs if d.kind == "text")
        print(f"DIFFER: '{args.file_a}' vs '{args.file_b}': "
              f"{structural} structural, {numeric} numeric, {text} text "
              f"(rtol={args.rtol:g} atol={args.atol:g} int-tol={args.int_tol})")
        for d in diffs[:args.max_report]:
            print("  " + str(d))
        if len(diffs) > args.max_report:
            print(f"  ... and {len(diffs) - args.max_report} more")
    return 1


if __name__ == "__main__":
    sys.exit(main())
