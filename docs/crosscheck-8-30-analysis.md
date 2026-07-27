# Why cminpack diverges on driver problem 8/30

This note documents the root cause of the one genuine divergence that the
FORTRAN-reference cross-check gate excludes on macOS: the `hybrd`/`hybrj`
driver problem **8 at dimension 30** (`crosscheck_exclude.txt`).

## The problem

Driver problem 8 is the **variably dimensioned function** (Moré, Garbow &
Hillstrom), whose exact solution is `x* = (1, 1, ..., 1)` with residual 0. The
driver runs it at dimension 30 from the hardest starting point (`100 * x0`).

## What we observe

On Apple clang, cminpack stops at `info = 4` ("iteration not making good
progress") with a final residual of **83.5**, having stalled at the spurious
point `x = (0.5, ..., 0.5)`. The FORTRAN reference reaches `x = (1, ..., 1)`
with residual `1.2136830e-13` (`info = 1`).

## Experiment: it is FMA contraction

Building the same C source with different floating-point-contraction settings
(driver `hyjdrvc`, problem 8 dim 30):

| build | final residual | info | nfev | njev |
|---|---|---|---|---|
| pure-C, clang `-O3` | 8.3476044e+01 | 4 | 25 | 2 |
| pure-C, clang `-O0` | 8.3476044e+01 | 4 | 25 | 2 |
| pure-C, clang `-O3 -ffp-contract=off` | **1.2136830e-13** | **1** | **12** | **2** |
| f2c, clang `-O3` | 8.3476044e+01 | 4 | 25 | 2 |
| FORTRAN reference (gfortran `-O0`) | **1.2136830e-13** | **1** | **12** | **2** |

Two facts stand out:

1. **`-ffp-contract=off` makes cminpack reproduce the FORTRAN result to the last
   printed digit** (`1.2136830e-13`) *and* the identical iteration counts
   (nfev = 12, njev = 2). With contraction on, it stalls at 83.5.
2. **clang contracts even at `-O0`** (its default is `-ffp-contract=on`, unlike
   gcc which only contracts at `-O1+`), which is why the `-O0` build diverges
   too. The pure-C and f2c builds diverge *identically*, so the sensitive
   arithmetic is in the shared trust-region/dogleg code, not in the pure-C
   cleanup of the f2c output.

## Answers to the two questions

**Is it just a last-significant-digit difference?** Yes. A fused multiply-add
computes `a*b + c` with a single rounding instead of two, so a contracted
expression differs from the separately-rounded one by about one unit in the last
place (ULP). The variably dimensioned function at dim 30 from the `100*x0` start
is ill-conditioned enough that a single ULP early in the iteration flips a
dogleg accept/reject decision; the two runs then follow different paths, one to
the true solution and one to the all-halves stationary point. Removing
contraction (`-ffp-contract=off`) collapses the difference and the run matches
FORTRAN exactly.

**Did the reference converge by chance?** No -- the un-contracted arithmetic
converges *robustly*. The FORTRAN reference (gfortran, no contraction at `-O0`)
converges, and so does cminpack once FMA is disabled -- to the identical value
and iteration counts. Moreover the FORTRAN build converges on this problem at
every gfortran optimisation level (`-O0` through `-O3`; see
`examples/crosscheck_matrix.sh`). It is specifically **clang's contraction of
the C expressions** that perturbs cminpack onto the diverging side of the
knife's edge. So FORTRAN is not "lucky"; rather, the exact computation converges
and a compiler-introduced 1-ULP FMA perturbation is what tips this particular
build over.

## Consequence for the cross-check

This is exactly the compiler-dependent coin-flip the cross-check tolerates: on
GNU/Linux GCC (the CI platform) cminpack matches FORTRAN and 8/30 is not
flagged, so it is only excluded to keep `make check` green on Apple clang.
Building cminpack with `-ffp-contract=off` removes the divergence entirely.
