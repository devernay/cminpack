#!/usr/bin/make 

# pick up your FORTRAN compiler
#F77=g77
F77=gfortran

# uncomment the following for FORTRAN MINPACK
#MINPACK=-lminpack
#F77C=$(F77)
#F77CFLAGS=-g

# uncomment the following for C MINPACK
MINPACK=../libminpack.a
F77C=gcc
F77CFLAGS=-I .. -g -Wall

# uncomment the following to debug using valgrind
#VALGRIND=valgrind --tool=memcheck

# Execute this script to extract the examples from the MINPACK
# documentation (/usr/share/doc/minpack-devel-19961126/minpack-documentation.txt)
# and compile and run them.  You may compare the numerical results
# with those in the same file.

FPGM=tchkder thybrd thybrd1 thybrj thybrj1 tlmder tlmder1 tlmdif \
	tlmdif1 tlmstr tlmstr1
PGM=tchkder_ thybrd_ thybrd1_ thybrj_ thybrj1_ tlmder_ tlmder1_ tlmdif_ \
	tlmdif1_ tlmstr_ tlmstr1_
CPGM=tcchkder tchybrd tchybrd1 tchybrj tchybrj1 tclmder tclmder1 tclmdif \
	tclmdif1 tclmstr tclmstr1
FSRC=tchkder.f thybrd.f thybrd1.f thybrj.f thybrj1.f tlmder.f tlmder1.f \
	tlmdif.f tlmdif1.f tlmstr.f tlmstr1.f
SRC=tchkder_.c thybrd_.c thybrd1_.c thybrj_.c thybrj1_.c tlmder_.c tlmder1_.c \
	tlmdif_.c tlmdif1_.c tlmstr_.c tlmstr1_.c
CSRC=tchkder.c thybrd.c thybrd1.c thybrj.c thybrj1.c tlmder.c tlmder1.c \
	tlmdif.c tlmdif1.c tlmstr.c tlmstr1.c
REF=test.ref ctest.ref ftest.ref

all: ${PGM} ${CPGM} ${FPGM}
	@echo "*****************************************************"
	@echo "Please type 'make test', 'make ctest' or 'make ftest'"
	@echo "make test: calling FORTRAN MINPACK from C"
	@echo "make ftest: calling FORTRAN MINPACK from FORTRAN"
	@echo "make ctest: calling C MINPACK from C"

test: ${PGM}
	for x in ${PGM}; do echo $$x; $(VALGRIND) ./$$x; done > test.out
	diff -u test.ref test.out

ctest: ${CPGM}
	for x in ${CPGM}; do echo $$x; $(VALGRIND) ./$$x; done > ctest.out
	diff -u ctest.ref ctest.out

ftest: ${FPGM} 
	for x in ${FPGM}; do echo $$x; $(VALGRIND) ./$$x; done > ftest.out
	diff -u ftest.ref ftest.out

clean:
	-rm $(PGM) $(CPGM) $(FPGM) $(FSRC) $(REF:.ref=.out) *~ #*#

${FSRC}: /usr/share/doc/minpack-devel-19961126/minpack-documentation.txt.gz
	zcat $< | awk  '					\
	/DRIVER FOR [A-Z1]+ EXAMPLE/{				\
	    pgm=tolower($$4);					\
	    oname="t" pgm ".f";					\
	    $$0 = substr($$0,3);				\
	    print >oname;					\
	    do {						\
		getline; $$0 = substr($$0,3);			\
		if (!/^ +Page$$/) print >>oname;		\
	    }							\
	    while (!/LAST CARD OF SUBROUTINE FCN/);		\
	    getline; $$0 = substr($$0,3); print >>oname;	\
	    getline; $$0 = substr($$0,3); print >>oname;	\
	}'

tchkder: tchkder.f
	$(F77) -o $@ $< $(MINPACK)
thybrd: thybrd.f
	$(F77) -o $@ $< $(MINPACK)
thybrd1: thybrd1.f
	$(F77) -o $@ $< $(MINPACK)
thybrj: thybrj.f
	$(F77) -o $@ $< $(MINPACK)
thybrj1: thybrj1.f
	$(F77) -o $@ $< $(MINPACK)
tlmder: tlmder.f
	$(F77) -o $@ $< $(MINPACK)
tlmder1: tlmder1.f
	$(F77) -o $@ $< $(MINPACK)
tlmdif: tlmdif.f
	$(F77) -o $@ $< $(MINPACK)
tlmdif1: tlmdif1.f
	$(F77) -o $@ $< $(MINPACK)
tlmstr: tlmstr.f
	$(F77) -o $@ $< $(MINPACK)
tlmstr1: tlmstr1.f
	$(F77) -o $@ $< $(MINPACK)

tchkder_: tchkder_.c					
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
thybrd_: thybrd_.c					
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
thybrd1_: thybrd1_.c					
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
thybrj_: thybrj_.c					
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
thybrj1_: thybrj1_.c					
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
tlmder_: tlmder_.c					
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
tlmder1_: tlmder1_.c					
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
tlmdif_: tlmdif_.c
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
tlmdif1_: tlmdif1_.c
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
tlmstr_: tlmstr_.c
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
tlmstr1_: tlmstr1_.c
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm

tcchkder: tchkder.c					
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
tchybrd: thybrd.c					
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
tchybrd1: thybrd1.c					
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
tchybrj: thybrj.c					
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
tchybrj1: thybrj1.c					
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
tclmder: tlmder.c					
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
tclmder1: tlmder1.c					
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
tclmdif: tlmdif.c
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
tclmdif1: tlmdif1.c
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
tclmstr: tlmstr.c
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
tclmstr1: tlmstr1.c
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
