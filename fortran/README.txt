This directory contains the original FORTRAN sources of the Minpack library,
downloaded from ftp://ftp.netlib.org/minpack/, with these differences:

- dpmpar.f: the IEEE machine constants were replaced by full-precision
  values (the original ones are kept as comments), because the original
  values were not accurate enough according to MACHAR.
- Makefile: not from netlib; it builds libminpack.a, used by the tests in
  ../examples (e.g. "make test MINPACK=../fortran/libminpack.a").

All the other *.f files are identical to the netlib ones.

For Copyright see CopyrightMINPACK.txt.


Other resources:

http://en.wikipedia.org/wiki/MINPACK
http://www.netlib.org/minpack/
