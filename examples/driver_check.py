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

Reference-gate mode (--reference-gate): treat FILE_A as the authoritative
FORTRAN reference and apply a one-sided regression gate instead of a symmetric
comparison. Following the original MINPACK test drivers, it uses each problem's
exit parameter (info): where FORTRAN converged (info < 5) it requires this build
to converge too (info < 5), and it ignores problems where FORTRAN itself did not
converge (info >= 5, e.g. the maximum number of evaluations reached on the
10*x0/100*x0 starts). Using info rather than the residual magnitude makes the
gate well defined even for problems with a genuinely nonzero minimum. Use
--exclude to skip the handful of deliberately-extreme coin-flip problems (see
crosscheck_exclude.txt); each failure line names the problem as "nprob/dim", the
exact token to add to that list.

Exit status: 0 if equivalent (or, in gate mode, no regression), 1 if a genuine
convergence discrepancy/regression is found, 2 on usage/IO error or unparseable
input.

Usage:
    driver_check.py [--atol A] [--rtol R] [--label NAME] FILE_A FILE_B
    driver_check.py --reference-gate [--exclude "n1/d1,n2/d2"] REF_A OUT_B
"""

import argparse
import re
import sys

_PROBLEM_RE = re.compile(r"^\s*problem\b", re.IGNORECASE)
_FNORM_RE = re.compile(
    r"final l2 norm of the residuals\s+([-+0-9.eEdD]+)", re.IGNORECASE)
_INFO_RE = re.compile(r"exit parameter\s+(\d+)", re.IGNORECASE)
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


def _prob_id(header):
    """Extract 'nprob/dim' from a problem header line, or None."""
    m = re.search(r"problem\s+(\d+)\s+dimensions?\s+(\d+)", header, re.IGNORECASE)
    return "%s/%s" % (m.group(1), m.group(2)) if m else None


def check(path_a, path_b, atol, rtol, label, gate=False, exclude=None):
    exclude = exclude or set()
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
    excluded = 0
    worst_rel = 0.0
    worst_abs = 0.0
    bad = []

    for i, ((hdr_a, blk_a), (_, blk_b)) in enumerate(zip(blocks_a, blocks_b), 1):
        problems += 1
        pid = _prob_id(hdr_a) or hdr_a.strip()
        if gate and pid in exclude:
            excluded += 1
            continue
        if is_chkder:
            va, vb = _error_vectors(blk_a), _error_vectors(blk_b)
            if len(va) != len(vb):
                bad.append("  problem %s (block %d): error vector length differs"
                           % (pid, i))
                continue
            for xa, xb in zip(va, vb):
                d = abs(xa - xb)
                worst_abs = max(worst_abs, d)
                denom = max(abs(xa), abs(xb))
                if denom:
                    worst_rel = max(worst_rel, d / denom)
                if not _close(xa, xb, atol, rtol):
                    if not gate:   # chkder is not a convergence gate
                        bad.append("  problem %s (block %d): error vector %.7g vs %.7g"
                                   % (pid, i, xa, xb))
                    break
        else:
            ma = _FNORM_RE.search("".join(blk_a))
            mb = _FNORM_RE.search("".join(blk_b))
            if not ma or not mb:
                bad.append("  problem %s (block %d): missing final residual norm"
                           % (pid, i))
                continue
            fa, fb = _to_float(ma.group(1)), _to_float(mb.group(1))
            if fa is None or fb is None:
                bad.append("  problem %s (block %d): unparseable residual "
                           "(%r vs %r)" % (pid, i, ma.group(1), mb.group(1)))
                continue
            d = abs(fa - fb)
            worst_abs = max(worst_abs, d)
            denom = max(abs(fa), abs(fb))
            if denom:
                worst_rel = max(worst_rel, d / denom)
            if gate:
                # file_a is the FORTRAN reference. Following the original
                # MINPACK test drivers, use the exit parameter (info): treat
                # info >= 5 as "FORTRAN did not converge" (e.g. the maximum
                # number of evaluations was reached) and ignore those problems.
                # Where FORTRAN converged (info < 5), require this build to
                # converge too (info < 5). This mirrors the fortran-lang/minpack
                # `info_original < 5` filter and, unlike a residual-magnitude
                # test, is well defined for problems with a nonzero minimum.
                ri = int(_INFO_RE.search("".join(blk_a)).group(1)) \
                    if _INFO_RE.search("".join(blk_a)) else None
                ci = int(_INFO_RE.search("".join(blk_b)).group(1)) \
                    if _INFO_RE.search("".join(blk_b)) else None
                if ri is None or ri >= 5:
                    continue        # FORTRAN did not converge here -> ignore
                if ci is None or ci >= 5:
                    bad.append("  problem %s (block %d): FORTRAN reference "
                               "converged (exit parameter %s) but this build did "
                               "NOT (exit parameter %s; final residual %.3g vs "
                               "%.3g)" % (pid, i, ri,
                                          "missing" if ci is None else ci, fa, fb))
            elif not _close(fa, fb, atol, rtol):
                bad.append("  problem %s (block %d): final residual %.7g vs %.7g "
                           "(rel %.2g)" % (pid, i, fa, fb,
                                           d / denom if denom else 0.0))

    kind = "error vectors" if is_chkder else "final residual norms"
    if bad:
        if gate:
            print("%s%d of %d problems converge in the FORTRAN reference but "
                  "NOT in this build:" % (tag, len(bad), problems))
        else:
            print("%s%d of %d problems reached a materially different result (%s):"
                  % (tag, len(bad), problems, kind))
        for line in bad:
            print(line)
        return 1

    if gate:
        _x = " (%d excluded)" % excluded if excluded else ""
        print("%sOK -- converges wherever the FORTRAN reference does "
              "(%d problems checked%s)" % (tag, problems - excluded, _x))
        return 0
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
    p.add_argument("--reference-gate", action="store_true",
                   help="treat FILE_A as the authoritative FORTRAN reference: "
                        "fail (exit 1) only on problems the reference converges "
                        "on (residual <= --atol) but FILE_B does not. Problems "
                        "the reference does not converge on are ignored -- a "
                        "different result there is acceptable.")
    p.add_argument("--exclude", default="",
                   help="comma-separated 'nprob/dim' problem ids to skip in "
                        "--reference-gate mode (e.g. '8/30,6/9,7/6,8/40'). Use "
                        "for the handful of deliberately-extreme problems whose "
                        "convergence is a compiler-dependent coin-flip.")
    args = p.parse_args(argv)
    excl = set(t.strip() for t in args.exclude.split(",") if t.strip())
    try:
        return check(args.file_a, args.file_b, args.atol, args.rtol, args.label,
                     gate=args.reference_gate, exclude=excl)
    except OSError as e:
        print("driver_check.py: %s" % e, file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
