====== C language readme for minpack ======

This is a C version of the minpack minimization package.
It has been derived from the fortran code using f2c and
some limited manual editing. Note that you need to link
against libf2c to use this version of minpack. Extern "C"
linkage permits the package routines to be called from C++.
Check ftp://netlib.bell-labs.com/netlib/f2c for the latest
f2c version. For general minpack info and test programs, see
the accompanying readme.txt and http://www.netlib.org/minpack/.

Type `make` to compile and `make install` to install in /usr/local
or modify the makefile to suit your needs.

This software has been tested on a RedHat 7.3 Linux machine -
usual 'use at your own risk' warnings apply.

Manolis Lourakis -- lourakis at ics forth gr, July 2002
	Institute of Computer Science,
	Foundation for Research and Technology - Hellas
	Heraklion, Crete, Greece

Repackaging by Frederic Devernay -- frederic dot devernay at m4x dot org:
- Added fortran and C examples
- Added documentation from Debian man pages
- Wrote pure C version
- Added covar() and covar_(), and use it in tlmdef/tlmdif
