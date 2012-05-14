#!/usr/bin/make 

# pick up your FORTRAN compiler
#F77=g77
F77=gfortran

# uncomment the following for FORTRAN MINPACK
#MINPACK=-lminpack
#F77C=$(F77)
#F77CFLAGS=-g

# uncomment the following for C MINPACK
MINPACK=../libcminpack.a
CC=gcc
CFLAGS=-g -Wall
CPPFLAGS=-I..

FMINPACK=../fortran/libminpack.a

# uncomment the following to debug using valgrind
#VALGRIND=valgrind --tool=memcheck

# Execute this script to extract the examples from the MINPACK
# documentation (../doc/minpack-documentation.txt)
# and compile and run them.  You may compare the numerical results
# with those in the same file.

FPGM=tchkder thybrd thybrd1 thybrj thybrj1 tlmder tlmder1 tlmdif \
	tlmdif1 tlmstr tlmstr1 ibmdpdr
PGM=tchkder_ thybrd_ thybrd1_ thybrj_ thybrj1_ tlmder_ tlmder1_ tlmdif_ \
	tlmdif1_ tlmstr_ tlmstr1_ tfdjac2_
CPGM=tchkderc thybrdc thybrd1c thybrjc thybrj1c tlmderc tlmder1c tlmdifc \
	tlmdif1c tlmstrc tlmstr1c tfdjac2c ibmdpdrc
DRVPGM=lmddrv lmfdrv lmsdrv hyjdrv hybdrv chkdrv
CDRVPGM=lmddrvc lmfdrvc lmsdrvc hyjdrvc hybdrvc chkdrvc
FSRCGEN=tchkder.f thybrd.f thybrd1.f thybrj.f thybrj1.f tlmder.f tlmder1.f \
	tlmdif.f tlmdif1.f tlmstr.f tlmstr1.f
FSRC=$(FSRCGEN) \
	lmddrv.f lmdipt.f ssqfcn.f ssqjac.f lmfdrv.f lmsdrv.f \
	hyjdrv.f hybipt.f vecfcn.f vecjac.f hybdrv.f \
	errjac.f chkdrv.f
SRC=tchkder_.c thybrd_.c thybrd1_.c thybrj_.c thybrj1_.c tlmder_.c tlmder1_.c \
	tlmdif_.c tlmdif1_.c tlmstr_.c tlmstr1_.c tfdjac2_.c
CSRC=tchkderc.c thybrdc.c thybrd1c.c thybrjc.c thybrj1c.c tlmderc.c tlmder1c.c \
	tlmdifc.c tlmdif1c.c tlmstrc.c tlmstr1c.c tfdjac2c.c \
	lmddrv.c lmdipt.c ssqfcn.c ssqjac.c lmfdrv.c lmsdrv.c \
	hyjdrv.c hybipt.c vecfcn.c vecjac.c hybdrv.c \
	errjac.c chkdrv.c

CDRVSRC=lmddrv.c lmfdrv.c lmsdrv.c hyjdrv.c hybdrv.c
REF=test.ref ctest.ref ftest.ref

all:
	@echo "*****************************************************"
	@echo "Please type 'make alltest' to run all tests at once, or:"
	@echo "make test: calling FORTRAN CMINPACK from C"
	@echo "make test MINPACK=../fortran/libminpack.a: calling FORTRAN MINPACK from C"
	@echo "make ctest: calling CMINPACK from C"
	@echo "make ftest: calling FORTRAN CMINPACK from FORTRAN"
	@echo "make ftest MINPACK=../fortran/libminpack.a: calling FORTRAN MINPACK from FORTRAN"
	@echo "Intensive (driver) tests:"
	@echo "make ctestdrv: calling CMINPACK from C"
	@echo "make ftestdrv: calling FORTRAN CMINPACK from FORTRAN"
	@echo "make ftestdrv MINPACK=../fortran/libminpack.a: calling FORTRAN MINPACK from FORTRAN"

.PHONY: test ctest ctestdrv ctestlmdrv ctesthydrv ctestchkdrv ftest ftestdrv ftestlmdrv ftesthydrv ftestchkdrv alltest
test: ${PGM}
	for x in ${PGM}; do echo $$x; $(VALGRIND) ./$$x > $$x.out; diff -u ref/`echo $$x.ref|sed -e s/_.ref/c.ref/` $$x.out ; done

ctest: ${CPGM}
	for x in ${CPGM}; do echo $$x; $(VALGRIND) ./$$x > $$x.out; diff -u ref/$$x.ref $$x.out ; done

ctestdrv: ${CDRVPGM} ctestlmdrv ctesthydrv ctestchkdrv

ctestlmdrv: lmddrvc lmfdrvc lmsdrvc
	for x in $^; do echo $$x; $(VALGRIND) ./$$x < testdata/lm.data > $$x.out; diff -u ref/$$x.ref $$x.out ; done

ctesthydrv: hyjdrvc hybdrvc
	for x in $^; do echo $$x; $(VALGRIND) ./$$x < testdata/hybrd.data > $$x.out; diff -u ref/$$x.ref $$x.out ; done

ctestchkdrv: chkdrvc
	for x in $^; do echo $$x; $(VALGRIND) ./$$x < testdata/chkder.data > $$x.out; diff -u ref/$$x.ref $$x.out ; done

ftest: ${FPGM}
	for x in ${FPGM}; do echo $$x; $(VALGRIND) ./$$x > $$x.out; diff -u ref/$$x.ref $$x.out ; done

ftestdrv: ${DRVPGM} ftestlmdrv ftesthydrv ftestchkdrv

ftestlmdrv: lmddrv lmfdrv lmsdrv
	for x in $^; do echo $$x; $(VALGRIND) ./$$x < testdata/lm.data > $$x.out; diff -u ref/$$x.ref $$x.out ; done

ftesthydrv: hyjdrv hybdrv
	for x in $^; do echo $$x; $(VALGRIND) ./$$x < testdata/hybrd.data > $$x.out; diff -u ref/$$x.ref $$x.out ; done

ftestchkdrv: chkdrv
	for x in $^; do echo $$x; $(VALGRIND) ./$$x < testdata/chkder.data > $$x.out; diff -u ref/$$x.ref $$x.out ; done

alltest: ../libcminpack.a ../fortran/libminpack.a
	$(MAKE) clean
	-$(MAKE) test
	$(MAKE) clean
	-$(MAKE) test MINPACK=../fortran/libminpack.a
	$(MAKE) clean
	-$(MAKE) ctest
	$(MAKE) clean
	-$(MAKE) ftest
	$(MAKE) clean
	-$(MAKE) ftest MINPACK=../fortran/libminpack.a
	$(MAKE) clean
	-$(MAKE) -k ctestdrv
	$(MAKE) clean
	-$(MAKE) -k ftestdrv
	$(MAKE) clean
	-$(MAKE) -k ftestdrv MINPACK=../fortran/libminpack.a
	$(MAKE) clean

../libcminpack.a:
	$(MAKE) -C ..

 ../fortran/libminpack.a:
	$(MAKE) -C ../fortran

clean:
	-rm -f $(PGM) $(DRVPGM) $(CPGM) $(CDRVPGM) $(CSRC:.c=.o) $(FPGM) $(FSRCGEN) $(FSRC:.f=.o) $(PGM:=.out) $(CPGM:=.out) $(CDRVPGM:=.out) $(FPGM:=.out) *~ #*#
	-rm -rf $(PGM:=.dSYM) $(CPGM:=.dSYM) $(CDRVPGM:=.dSYM) $(FPGM:=.dSYM)

${FSRCGEN}: ../doc/minpack-documentation.txt
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

.f: $(MINPACK)
	$(F77) -o $@ $< $(MINPACK)

.c:
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $< $(MINPACK) -lm

.c.o:
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

lmddrvc: $(MINPACK) lmddrv.o lmdipt.o ssqfcn.o ssqjac.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) -lm

lmfdrvc: $(MINPACK) lmfdrv.o lmdipt.o ssqfcn.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) -lm

lmsdrvc: $(MINPACK) lmsdrv.o lmdipt.o ssqfcn.o ssqjac.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) -lm

hyjdrvc: $(MINPACK) hyjdrv.o hybipt.o vecfcn.o vecjac.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) -lm

hybdrvc: $(MINPACK) hybdrv.o hybipt.o vecfcn.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) -lm

chkdrvc: $(MINPACK) chkdrv.o hybipt.o vecfcn.o errjac.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) -lm

ibmdpdrc: $(MINPACK) ibmdpdr.o dmchar.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) -lm

lmddrv: $(MINPACK) lmddrv.f lmdipt.f ssqfcn.f ssqjac.f
	$(F77) $(FFLAGS) -o $@ $^ $(MINPACK)

lmfdrv: $(MINPACK) lmfdrv.f lmdipt.f ssqfcn.f
	$(F77) $(FFLAGS) -o $@ $^ $(MINPACK)

lmsdrv: $(MINPACK) lmsdrv.f lmdipt.f ssqfcn.f ssqjac.f
	$(F77) $(FFLAGS) -o $@ $^ $(MINPACK)

hyjdrv: $(MINPACK) hyjdrv.f hybipt.f vecfcn.f vecjac.f
	$(F77) $(FFLAGS) -o $@ $^ $(MINPACK)

hybdrv: $(MINPACK) hybdrv.f hybipt.f vecfcn.f
	$(F77) $(FFLAGS) -o $@ $^ $(MINPACK)

chkdrv: $(MINPACK) chkdrv.f hybipt.f vecfcn.f errjac.f
	$(F77) $(FFLAGS) -o $@ $^ $(MINPACK)

ibmdpdr: $(MINPACK) ibmdpdr.f dmchar.f
	$(F77) $(FFLAGS) -o $@ $^ $(MINPACK)
