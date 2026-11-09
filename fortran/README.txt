This directory contains the original FORTRAN sources of the Minpack library,
downloaded from ftp://ftp.netlib.org/minpack/, with these differences:

- dpmpar.f: the IEEE machine constants were replaced by full-precision
  values (the original ones are kept as comments), because the original
  values were not accurate enough according to MACHAR.
- Makefile: not from netlib; it builds libminpack.a, used by the tests in
  ../examples (e.g. "make test MINPACK=../fortran/libminpack.a").

All the other *.f files are identical to the netlib ones.

The C and f2c versions of cminpack differ from these sources in one
deliberate way: when the user function aborts (iflag < 0) during a
finite-difference Jacobian (fdjac1, fdjac2), they restore x before
returning, whereas the FORTRAN code returns with x perturbed.

For Copyright see CopyrightMINPACK.txt.


Other resources:

http://en.wikipedia.org/wiki/MINPACK
http://www.netlib.org/minpack/
