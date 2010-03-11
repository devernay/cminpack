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
# documentation (../doc/minpack-documentation.txt)
# and compile and run them.  You may compare the numerical results
# with those in the same file.

FPGM=tchkder thybrd thybrd1 thybrj thybrj1 tlmder tlmder1 tlmdif \
	tlmdif1 tlmstr tlmstr1
PGM=tchkder_ thybrd_ thybrd1_ thybrj_ thybrj1_ tlmder_ tlmder1_ tlmdif_ \
	tlmdif1_ tlmstr_ tlmstr1_ tfdjac2_
CPGM=tchkderc thybrdc thybrd1c thybrjc thybrj1c tlmderc tlmder1c tlmdifc \
	tlmdif1c tlmstrc tlmstr1c tfdjac2c
FSRC=tchkder.f thybrd.f thybrd1.f thybrj.f thybrj1.f tlmder.f tlmder1.f \
	tlmdif.f tlmdif1.f tlmstr.f tlmstr1.f
SRC=tchkder_.c thybrd_.c thybrd1_.c thybrj_.c thybrj1_.c tlmder_.c tlmder1_.c \
	tlmdif_.c tlmdif1_.c tlmstr_.c tlmstr1_.c tfdjac2_.c
CSRC=tchkderc.c thybrdc.c thybrd1c.c thybrjc.c thybrj1c.c tlmderc.c tlmder1c.c \
	tlmdifc.c tlmdif1c.c tlmstrc.c tlmstr1c.c tfdjac2c.c
REF=test.ref ctest.ref ftest.ref

all:
	@echo "*****************************************************"
	@echo "Please type 'make test', 'make ctest' or 'make ftest'"
	@echo "make test: calling FORTRAN MINPACK from C"
	@echo "make ftest: calling FORTRAN MINPACK from FORTRAN"
	@echo "make ctest: calling C MINPACK from C"

test: ${PGM}
	for x in ${PGM}; do echo $$x; $(VALGRIND) ./$$x > $$x.out; diff -u ref/$$x.ref $$x.out ; done

ctest: ${CPGM}
	for x in ${CPGM}; do echo $$x; $(VALGRIND) ./$$x > $$x.out; diff -u ref/$$x.ref $$x.out ; done

ftest: ${FPGM} 
	for x in ${FPGM}; do echo $$x; $(VALGRIND) ./$$x > $$x.out; diff -u ref/$$x.ref $$x.out ; done

clean:
	-rm $(PGM) $(CPGM) $(FPGM) $(FSRC) $(PGM:=.out) $(CPGM:=.out) $(FPGM:=.out) *~ #*#

${FSRC}: ../doc/minpack-documentation.txt
	cat $< | awk  '					\
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
tfdjac2_: tfdjac2_.c
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm

tchkderc: tchkderc.c					
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
thybrdc: thybrdc.c					
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
thybrd1c: thybrd1c.c					
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
thybrjc: thybrjc.c					
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
thybrj1c: thybrj1c.c					
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
tlmderc: tlmderc.c					
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
tlmder1c: tlmder1c.c					
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
tlmdifc: tlmdifc.c
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
tlmdif1c: tlmdif1c.c
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
tlmstrc: tlmstrc.c
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
tlmstr1c: tlmstr1c.c
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
tfdjac2c: tfdjac2c.c
	$(F77C) $(F77CFLAGS) -o $@ $< $(MINPACK) -lm
