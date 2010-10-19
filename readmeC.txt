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

version 1.0.0 (24/04/2007):
- Added fortran and C examples
- Added documentation from Debian man pages
- Wrote pure C version
- Added covar() and covar_(), and use it in tlmdef/tlmdif

version 1.0.1 (17/12/2007):
- bug fix in covar() and covar_(), the computation of tolr caused a
  segfault (signaled by Timo Hartmann).

version 1.0.2 (27/02/2009):
- Added Xcode and Visual Studio project files

version 1.0.3 (18/03/2010):
- Added CMake support.
- XCode build is now Universal.
- Added tfdjac2_ and tfdjac2c examples, which test the accuracy of a
  finite-differences approximation of the Jacobian.
- Bug fix in tlmstr1 (signaled by Thomas Capricelli).

version 1.0.4 (18/10/2010):
- Support for shared library building using CMake, thanks to Goeffrey
  Biggs and Radu Bogdan Rusu from Willow Garage. Shared libraries can be
  enabled using cmake options, as in;
cmake -DUSE_FPIC=ON -DSHARED_LIBS=ON -DBUILD_EXAMPLES=OFF path_to_sources
