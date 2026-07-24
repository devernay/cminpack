#!/usr/bin/env python3
"""Semantic check for cminpack intensive-driver output (examples/*drv*).

Unlike compare.py (which flags every token-level difference), this script asks
the only question that matters for the hard More/Garbow/Hillstrom problems: did
two driver runs converge to *equivalent results*?  It deliberately ignores the
iteration counts, function/Jacobian-evaluation counts, exit parameter and
last-digit noise. Those legitimately differ between compilers -- gcc and
gfortran contract "a*b + c" into a fused multiply-add at different places (see
README.md, "Numerical differences from FORTRAN MINPACK") -- so they are expected
and harmless, and printing them only alarms the user.

What it actually compares:

* Solver drivers (lmddrv, lmfdrv, lmsdrv, hyjdrv, hybdrv): the *final L2 norm
  of the residuals* of each problem -- i.e. the quality of the solution found.
  Two runs are equivalent when every problem's final norm agrees within
  tolerance. A problem where both
  norms are below --atol is treated as "converged to zero" regardless of how
  many orders of magnitude the tiny residuals differ by.

* chkder (chkdrv): the per-problem error vectors -- chkder's measure of how well
  the analytic Jacobian matches finite differences (~1.0 = consistent, ~0.0 =
  inconsistent) -- compared elementwise within tolerance.

Only genuinely significant differences are reported: a problem that converged
to a materially different residual, or runs with a different number of problems.
Everything else prints a short, reassuring one-line summary.

Exit status: 0 if the two runs are equivalent, 1 if a genuine convergence
discrepancy is found, 2 on usage/IO error or unparseable input.

Usage:
    driver_check.py [--atol A] [--rtol R] [--label NAME] FILE_A FILE_B
"""

import argparse
import re
import sys

_PROBLEM_RE = re.compile(r"^\s*problem\b", re.IGNORECASE)
_FNORM_RE = re.compile(
    r"final l2 norm of the residuals\s+([-+0-9.eEdD]+)", re.IGNORECASE)
_NUM_RE = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?")


def _to_float(tok):
    try:
        return float(tok.replace("D", "e").replace("d", "e"))
    except ValueError:
        return None


def _split_problems(lines):
    """Split a driver output into (header_line, block_text) problem blocks."""
    blocks = []
    cur_header = None
    cur = []
    for ln in lines:
        if _PROBLEM_RE.match(ln):
            if cur_header is not None:
                blocks.append((cur_header, cur))
            cur_header = ln.strip()
            cur = [ln]
        elif cur_header is not None:
            cur.append(ln)
    if cur_header is not None:
        blocks.append((cur_header, cur))
    return blocks


def _error_vectors(block):
    """Return the numbers following each 'error vector' line in a chkder block."""
    vals = []
    grab = False
    for ln in block:
        low = ln.lower()
        if "error vector" in low:
            grab = True
            continue
        if grab:
            nums = [_to_float(t) for t in _NUM_RE.findall(ln)]
            nums = [x for x in nums if x is not None]
            if nums:
                vals.extend(nums)
                grab = False  # error vector is on the first non-empty line
    return vals


def _close(a, b, atol, rtol):
    return abs(a - b) <= atol + rtol * max(abs(a), abs(b))


def check(path_a, path_b, atol, rtol, label):
    with open(path_a, encoding="utf-8", errors="replace") as f:
        blocks_a = _split_problems(f.readlines())
    with open(path_b, encoding="utf-8", errors="replace") as f:
        blocks_b = _split_problems(f.readlines())

    tag = label + ": " if label else ""

    if len(blocks_a) != len(blocks_b):
        print("%sdifferent number of problems (%d vs %d) -- structural mismatch"
              % (tag, len(blocks_a), len(blocks_b)))
        return 1
    if not blocks_a:
        print("%sno problem blocks found (nothing to check)" % tag)
        return 2

    is_chkder_a = not any(_FNORM_RE.search("".join(b)) for _, b in blocks_a)
    is_chkder_b = not any(_FNORM_RE.search("".join(b)) for _, b in blocks_b)
    if is_chkder_a != is_chkder_b:
        print("%sincompatible output types (one reports residual norms, the "
              "other does not)" % tag)
        return 1
    is_chkder = is_chkder_a

    problems = 0
    worst_rel = 0.0
    worst_abs = 0.0
    bad = []

    for i, ((hdr_a, blk_a), (_, blk_b)) in enumerate(zip(blocks_a, blocks_b), 1):
        problems += 1
        if is_chkder:
            va, vb = _error_vectors(blk_a), _error_vectors(blk_b)
            if len(va) != len(vb):
                bad.append("  problem block %d (%s): error vector length differs"
                           % (i, hdr_a))
                continue
            for xa, xb in zip(va, vb):
                d = abs(xa - xb)
                worst_abs = max(worst_abs, d)
                denom = max(abs(xa), abs(xb))
                if denom:
                    worst_rel = max(worst_rel, d / denom)
                if not _close(xa, xb, atol, rtol):
                    bad.append("  problem block %d (%s): error vector %.7g vs %.7g"
                               % (i, hdr_a, xa, xb))
                    break
        else:
            ma = _FNORM_RE.search("".join(blk_a))
            mb = _FNORM_RE.search("".join(blk_b))
            if not ma or not mb:
                bad.append("  problem block %d (%s): missing final residual norm"
                           % (i, hdr_a))
                continue
            fa, fb = _to_float(ma.group(1)), _to_float(mb.group(1))
            if fa is None or fb is None:
                bad.append("  problem block %d (%s): unparseable residual "
                           "(%r vs %r)" % (i, hdr_a, ma.group(1), mb.group(1)))
                continue
            d = abs(fa - fb)
            worst_abs = max(worst_abs, d)
            denom = max(abs(fa), abs(fb))
            if denom:
                worst_rel = max(worst_rel, d / denom)
            if not _close(fa, fb, atol, rtol):
                bad.append("  problem block %d (%s): final residual %.7g vs %.7g "
                           "(rel %.2g)" % (i, hdr_a, fa, fb,
                                           d / denom if denom else 0.0))

    kind = "error vectors" if is_chkder else "final residual norms"
    if bad:
        print("%s%d of %d problems reached a materially different result (%s):"
              % (tag, len(bad), problems, kind))
        for line in bad:
            print(line)
        return 1

    print("%sOK -- %d problems, all converged to equivalent %s "
          "(worst rel %.2g, worst abs %.2g)"
          % (tag, problems, kind, worst_rel, worst_abs))
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("file_a")
    p.add_argument("file_b")
    p.add_argument("--atol", type=float, default=1e-4,
                   help="residuals/values at or below this are treated as "
                        "equal ('converged to zero'); default 1e-4")
    p.add_argument("--rtol", type=float, default=1e-1,
                   help="relative tolerance on the compared values (default "
                        "1e-1 = 10%%). This is deliberately loose. Two runs of "
                        "the same algorithm differing only in iteration path "
                        "reach the same minimum/root, so the aim is to flag "
                        "only qualitatively different outcomes -- e.g. one run "
                        "converged, the other did not -- not last-digit noise.")
    p.add_argument("--label", default="",
                   help="prefix for the summary line (e.g. the driver name)")
    args = p.parse_args(argv)
    try:
        return check(args.file_a, args.file_b, args.atol, args.rtol, args.label)
    except OSError as e:
        print("driver_check.py: %s" % e, file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
