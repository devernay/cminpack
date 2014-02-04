#!/usr/bin/make 

# pick up your FORTRAN compiler
#F77=g77
F77=gfortran

# uncomment the following for FORTRAN MINPACK
#MINPACK=-lminpack
#F77C=$(F77)
#F77CFLAGS=-g

# uncomment the following for C MINPACK
MINPACK=../libcminpack$(LIBSUFFIX).a
CC=gcc
CFLAGS=-g -Wall
CPPFLAGS=-I..

FMINPACK=../fortran/libminpack.a

# uncomment the following to debug using valgrind
#VALGRIND=valgrind --tool=memcheck

### The default configuration is to compile the double precision version

### configuration for the LAPACK/BLAS (double precision) version:
## make LIBSUFFIX= CFLAGS="-O3 -g -Wall -Wextra -D__cminpack_float__"
#LIBSUFFIX=s
#CFLAGS="-O3 -g -Wall -Wextra -DUSE_CBLAS -DUSE_LAPACK"
CFLAGS_L=$(CFLAGS) -DUSE_CBLAS -DUSE_LAPACK
LDADD_L=-framework vecLib

### configuration for the float (single precision) version:
## make LIBSUFFIX=s CFLAGS="-O3 -g -Wall -Wextra -D__cminpack_float__"
#LIBSUFFIX=s
#CFLAGS="-O3 -g -Wall -Wextra -D__cminpack_float__"
CFLAGS_F=$(CFLAGS) -D__cminpack_float__

### configuration for the half (half precision) version:
## make LIBSUFFIX=h CFLAGS="-O3 -g -Wall -Wextra -I/opt/local/include -D__cminpack_half__" LDADD="-L/opt/local/lib -lHalf" CC=g++
#LIBSUFFIX=h
#CFLAGS="-O3 -g -Wall -Wextra -I/opt/local/include -D__cminpack_half__"
#LDADD="-L/opt/local/lib -lHalf"
#CC=g++
CFLAGS_H=$(CFLAGS) -I/opt/local/include -D__cminpack_half__
LDADD_H=-L/opt/local/lib -lHalf
CC_H=$(CXX)

FPGM=tchkder thybrd thybrd1 thybrj thybrj1 tlmder tlmder1 tlmdif \
	tlmdif1 tlmstr tlmstr1 ibmdpdr
PGM=tchkder_ thybrd_ thybrd1_ thybrj_ thybrj1_ tlmder_ tlmder1_ tlmdif_ \
	tlmdif1_ tlmstr_ tlmstr1_ tfdjac2_ ibmdpdr_
CPGM=tchkderc thybrdc thybrd1c thybrjc thybrj1c tlmderc tlmder1c tlmdifc \
	tlmdif1c tlmstrc tlmstr1c tfdjac2c ibmdpdrc \
	tchkderc_box tlmderc_box thybrjc_box
FDRVPGM=lmddrv_ lmfdrv_ lmsdrv_ hyjdrv_ hybdrv_ chkdrv_
DRVPGM=lmddrv lmfdrv lmsdrv hyjdrv hybdrv chkdrv
CDRVPGM=lmddrvc lmfdrvc lmsdrvc hyjdrvc hybdrvc chkdrvc
FSRCGEN=tchkder.f thybrd.f thybrd1.f thybrj.f thybrj1.f tlmder.f tlmder1.f \
	tlmdif.f tlmdif1.f tlmstr.f tlmstr1.f
FSRC=$(FSRCGEN) \
	lmddrv.f lmdipt.f ssqfcn.f ssqjac.f lmfdrv.f lmsdrv.f \
	hyjdrv.f hybipt.f vecfcn.f vecjac.f hybdrv.f \
	errjac.f chkdrv.f ibmdpdr.f machar.f
SRC=tchkder_.c thybrd_.c thybrd1_.c thybrj_.c thybrj1_.c tlmder_.c tlmder1_.c \
	tlmdif_.c tlmdif1_.c tlmstr_.c tlmstr1_.c tfdjac2_.c \
	lmddrv_.c lmfdrv_.c lmsdrv_.c \
	hyjdrv_.c hybdrv_.c \
	chkdrv_.c ibmdpdr_.c
CSRC=tchkderc.c thybrdc.c thybrd1c.c thybrjc.c thybrj1c.c tlmderc.c tlmder1c.c \
	tlmdifc.c tlmdif1c.c tlmstrc.c tlmstr1c.c tfdjac2c.c \
	lmddrv.c lmdipt.c ssqfcn.c ssqjac.c lmfdrv.c lmsdrv.c \
	hyjdrv.c hybipt.c vecfcn.c vecjac.c hybdrv.c \
	errjac.c chkdrv.c ibmdpdr.c machar.c

CDRVSRC=lmddrv.c lmfdrv.c lmsdrv.c hyjdrv.c hybdrv.c
REF=test.ref ctest.ref ftest.ref

all:
	@echo "*****************************************************"
	@echo "Please type 'make check' to run all tests at once, or:"
	@echo "*****************************************************"
	@echo "make test: calling FORTRAN CMINPACK from C"
	@echo "make test MINPACK=../fortran/libminpack.a: calling FORTRAN MINPACK from C"
	@echo "make ctest: calling CMINPACK from C"
	@echo "make ftest: calling FORTRAN CMINPACK from FORTRAN"
	@echo "make ftest MINPACK=../fortran/libminpack.a: calling FORTRAN MINPACK from FORTRAN"
	@echo "Intensive (driver) tests:"
	@echo "make testdrv: calling FORTRAN CMINPACK from C"
	@echo "make testdrv MINPACK=../fortran/libminpack.a: calling FORTRAN MINPACK from C"
	@echo "make ctestdrv: calling CMINPACK from C"
	@echo "make ftestdrv: calling FORTRAN CMINPACK from FORTRAN"
	@echo "make ftestdrv MINPACK=../fortran/libminpack.a: calling FORTRAN MINPACK from FORTRAN"

.PHONY: test ctest ctestdrv ctestlmdrv ctesthydrv ctestchkdrv ftest ftestdrv ftestlmdrv ftesthydrv ftestchkdrv check checkdoublec checkdouble checkfloatc checkfloat checkhalf checkfail
test: $(MINPACK) $(PGM)
	@echo "*** Running standard tests (calling MINPACK from C using $(MINPACK))"
	for x in $(PGM); do echo $$x; $(VALGRIND) ./$$x > $$x.out; diff -u ref/$(LIBSUFFIX)`echo $$x.ref|sed -e s/_.ref/c.ref/` $$x.out ; done

testdrv: $(MINPACK) $(FDRVPGM) testlmdrv testhydrv testchkdrv

testlmdrv: lmddrv_ lmfdrv_ lmsdrv_
	@echo "*** Running LM tests (calling MINPACK from C using $(MINPACK))"
	for x in $^; do echo $$x; $(VALGRIND) ./$$x < testdata/lm.data > $$x.out; diff -u ref/$(LIBSUFFIX)`echo $$x |sed -e s/_/c/`.ref $$x.out ; done

testhydrv: hyjdrv_ hybdrv_
	@echo "*** Running HY tests (calling MINPACK from C using $(MINPACK))"
	for x in $^; do echo $$x; $(VALGRIND) ./$$x < testdata/hybrd.data > $$x.out; diff -u ref/$(LIBSUFFIX)`echo $$x |sed -e s/_/c/`.ref $$x.out ; done

testchkdrv:chkdrv_
	@echo "*** Running CHK tests (calling MINPACK from C using $(MINPACK))"
	for x in $^; do echo $$x; $(VALGRIND) ./$$x < testdata/chkder.data > $$x.out; diff -u ref/$(LIBSUFFIX)`echo $$x |sed -e s/_/c/`.ref $$x.out ; done

ctest: $(MINPACK) $(CPGM)
	@echo "*** Running standard tests (calling CMINPACK from C using $(MINPACK))"
	for x in $(CPGM); do echo $$x; $(VALGRIND) ./$$x > $$x.out; diff -u ref/$(LIBSUFFIX)$$x.ref $$x.out ; done

ctestdrv: $(MINPACK) $(CDRVPGM) ctestlmdrv ctesthydrv ctestchkdrv

ctestlmdrv: lmddrvc lmfdrvc lmsdrvc
	@echo "*** Running LM tests (calling CMINPACK from C using $(MINPACK))"
	for x in $^; do echo $$x; $(VALGRIND) ./$$x < testdata/lm.data > $$x.out; diff -u ref/$(LIBSUFFIX)$$x.ref $$x.out ; done

ctesthydrv: hyjdrvc hybdrvc
	@echo "*** Running HY tests (calling CMINPACK from C using $(MINPACK))"
	for x in $^; do echo $$x; $(VALGRIND) ./$$x < testdata/hybrd.data > $$x.out; diff -u ref/$(LIBSUFFIX)$$x.ref $$x.out ; done

ctestchkdrv:chkdrvc
	@echo "*** Running CHK tests (calling CMINPACK from C using $(MINPACK))"
	for x in $^; do echo $$x; $(VALGRIND) ./$$x < testdata/chkder.data > $$x.out; diff -u ref/$(LIBSUFFIX)$$x.ref $$x.out ; done

ftest: $(MINPACK) $(FPGM)
	@echo "*** Running standard tests (calling MINPACK from FORTRAN using $(MINPACK))"
	for x in ${FPGM}; do echo $$x; $(VALGRIND) ./$$x > $$x.out; diff -u ref/$(LIBSUFFIX)$$x.ref $$x.out ; done

ftestdrv: $(MINPACK) $(DRVPGM) ftestlmdrv ftesthydrv ftestchkdrv

ftestlmdrv: lmddrv lmfdrv lmsdrv
	@echo "*** Running LM tests (calling MINPACK from FORTRAN using $(MINPACK))"
	for x in $^; do echo $$x; $(VALGRIND) ./$$x < testdata/lm.data > $$x.out; diff -u ref/$(LIBSUFFIX)$$x.ref $$x.out ; done

ftesthydrv: hyjdrv hybdrv
	@echo "*** Running HY tests (calling MINPACK from FORTRAN using $(MINPACK))"
	for x in $^; do echo $$x; $(VALGRIND) ./$$x < testdata/hybrd.data > $$x.out; diff -u ref/$(LIBSUFFIX)$$x.ref $$x.out ; done

ftestchkdrv: chkdrv
	@echo "*** Running CHK tests (calling MINPACK from FORTRAN using $(MINPACK))"
	for x in $^; do echo $$x; $(VALGRIND) ./$$x < testdata/chkder.data > $$x.out; diff -u ref/$(LIBSUFFIX)$$x.ref $$x.out ; done

check: checkdouble checkfloat

checkdoublec:
	$(MAKE) -C .. double 2>&1 > /dev/null
	$(MAKE) clean $(CPGM) LIBSUFFIX= 2>&1 > /dev/null
	-$(MAKE) ctest LIBSUFFIX=
	$(MAKE) clean 2>&1 > /dev/null

checkdouble: checkdoublec ../fortran/libminpack.a
	$(MAKE) clean $(PGM) LIBSUFFIX= 2>&1 > /dev/null
	-$(MAKE) test LIBSUFFIX=
	$(MAKE) clean $(PGM) LIBSUFFIX= MINPACK=../fortran/libminpack.a 2>&1 > /dev/null
	-$(MAKE) test LIBSUFFIX= MINPACK=../fortran/libminpack.a
	$(MAKE) clean $(FPGM) LIBSUFFIX= 2>&1 > /dev/null
	-$(MAKE) ftest LIBSUFFIX=
	$(MAKE) clean $(FPGM) LIBSUFFIX= MINPACK=../fortran/libminpack.a 2>&1 > /dev/null
	-$(MAKE) ftest LIBSUFFIX= MINPACK=../fortran/libminpack.a
	$(MAKE) clean $(FDRVPGM) LIBSUFFIX= 2>&1 > /dev/null
	-$(MAKE) -k testdrv LIBSUFFIX=
	$(MAKE) clean $(CDRVPGM) LIBSUFFIX= 2>&1 > /dev/null
	-$(MAKE) -k ctestdrv LIBSUFFIX=
	$(MAKE) clean 2>&1 > /dev/null

checkfloatc:
	$(MAKE) -C .. float 2>&1 > /dev/null
	$(MAKE) clean $(CPGM) LIBSUFFIX=s CFLAGS="$(CFLAGS_F)" 2>&1 > /dev/null
	-$(MAKE) ctest LIBSUFFIX=s CFLAGS="$(CFLAGS_F)" 
	$(MAKE) clean LIBSUFFIX=s 2>&1 > /dev/null

checkfloat: checkfloatc
	$(MAKE) clean $(PGM) LIBSUFFIX=s CFLAGS="$(CFLAGS_F)" 2>&1 > /dev/null
	-$(MAKE) test LIBSUFFIX=s CFLAGS="$(CFLAGS_F)"
	$(MAKE) clean $(FDRVPGM) LIBSUFFIX=s CFLAGS="$(CFLAGS_F)" 2>&1 > /dev/null
	-$(MAKE) -k testdrv LIBSUFFIX=s CFLAGS="$(CFLAGS_F)" 
	$(MAKE) clean $(CDRVPGM) LIBSUFFIX=s CFLAGS="$(CFLAGS_F)" 2>&1 > /dev/null
	-$(MAKE) -k ctestdrv LIBSUFFIX=s CFLAGS="$(CFLAGS_F)"
	$(MAKE) clean LIBSUFFIX=s 2>&1 > /dev/null

checklapack:
	$(MAKE) -C .. lapack 2>&1 > /dev/null
	$(MAKE) clean $(CPGM) LIBSUFFIX=l CFLAGS="$(CFLAGS_L)" LDADD="$(LDADD_L)" 2>&1 > /dev/null
	-$(MAKE) ctest LIBSUFFIX=l CFLAGS="$(CFLAGS_L)" LDADD="$(LDADD_L)"
	$(MAKE) clean $(CDRVPGM) LIBSUFFIX=l CFLAGS="$(CFLAGS_L)" LDADD="$(LDADD_L)" 2>&1 > /dev/null
	-$(MAKE) -k ctestdrv LIBSUFFIX=l CFLAGS="$(CFLAGS_L)" LDADD="$(LDADD_L)"

checkhalf:
	$(MAKE) -C .. half 2>&1 > /dev/null
	$(MAKE) clean $(PGM) LIBSUFFIX=h CFLAGS="$(CFLAGS_H)" LDADD="$(LDADD_H)" CC="$(CC_H)" 2>&1 > /dev/null
	-$(MAKE) test LIBSUFFIX=h CFLAGS="$(CFLAGS_H)" LDADD="$(LDADD_H)" CC="$(CC_H)"
	$(MAKE) clean $(CPGM) LIBSUFFIX=h CFLAGS="$(CFLAGS_H)" LDADD="$(LDADD_H)" CC="$(CC_H)" 2>&1 > /dev/null
	-$(MAKE) ctest LIBSUFFIX=h CFLAGS="$(CFLAGS_H)" LDADD="$(LDADD_H)" CC="$(CC_H)"
	$(MAKE) clean $(FDRVPGM) LIBSUFFIX=h CFLAGS="$(CFLAGS_H)" LDADD="$(LDADD_H)" CC="$(CC_H)" 2>&1 > /dev/null
	-$(MAKE) -k testdrv LIBSUFFIX=h CFLAGS="$(CFLAGS_H)" LDADD="$(LDADD_H)" CC="$(CC_H)"
	$(MAKE) clean $(CDRVPGM) LIBSUFFIX=h CFLAGS="$(CFLAGS_H)" LDADD="$(LDADD_H)" CC="$(CC_H)" 2>&1 > /dev/null
	-$(MAKE) -k ctestdrv LIBSUFFIX=h CFLAGS="$(CFLAGS_H)" LDADD="$(LDADD_H)" CC="$(CC_H)"
#	$(MAKE) clean LIBSUFFIX=h 2>&1 > /dev/null
#	$(MAKE) -C .. clean LIBSUFFIX=h 2>&1 > /dev/null

checkfail: ../fortran/libminpack.a check
	$(MAKE) clean $(FDRVPGM) LIBSUFFIX= MINPACK=../fortran/libminpack.a 2>&1 > /dev/null
	-$(MAKE) -k testdrv LIBSUFFIX= MINPACK=../fortran/libminpack.a
	$(MAKE) clean $(DRVPGM) LIBSUFFIX= 2>&1 > /dev/null
	-$(MAKE) -k ftestdrv LIBSUFFIX=
	$(MAKE) clean $(DRVPGM) LIBSUFFIX= MINPACK=../fortran/libminpack.a 2>&1 > /dev/null
	-$(MAKE) -k ftestdrv LIBSUFFIX= MINPACK=../fortran/libminpack.a
	$(MAKE) clean LIBSUFFIX= 2>&1 > /dev/null

../libcminpack$(LIBSUFFIX).a:
	$(MAKE) -C ..

../fortran/libminpack.a:
	$(MAKE) -C ../fortran

clean:
	-rm -f $(PGM) $(FDRVPGM) $(DRVPGM) $(CPGM) $(CDRVPGM) $(SRC:.c=.o) $(CSRC:.c=.o) $(FPGM) $(FSRCGEN) $(FSRC:.f=.o) $(PGM:=.out) $(CPGM:=.out) $(FDRVPGM:=.out) $(DRVPGM:=.out) $(CDRVPGM:=.out) $(FPGM:=.out) *~ #*#
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
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $< $(MINPACK) $(LDADD) -lm

tchkderc_box: tchkderc.c
	$(CC) $(CPPFLAGS) -DBOX_CONSTRAINTS $(CFLAGS) -o $@ $< $(MINPACK) $(LDADD) -lm

tlmderc_box: tlmderc.c
	$(CC) $(CPPFLAGS) -DBOX_CONSTRAINTS $(CFLAGS) -o $@ $< $(MINPACK) $(LDADD) -lm

thybrjc_box: thybrjc.c
	$(CC) $(CPPFLAGS) -DBOX_CONSTRAINTS $(CFLAGS) -o $@ $< $(MINPACK) $(LDADD) -lm

.c.o:
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

lmddrv_: $(MINPACK) lmddrv_.o lmdipt.o ssqfcn.o ssqjac.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) $(LDADD) -lm

lmfdrv_: $(MINPACK) lmfdrv_.o lmdipt.o ssqfcn.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) $(LDADD) -lm

lmsdrv_: $(MINPACK) lmsdrv_.o lmdipt.o ssqfcn.o ssqjac.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) $(LDADD) -lm

hyjdrv_: $(MINPACK) hyjdrv_.o hybipt.o vecfcn.o vecjac.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) $(LDADD) -lm

hybdrv_: $(MINPACK) hybdrv_.o hybipt.o vecfcn.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) $(LDADD) -lm

chkdrv_: $(MINPACK) chkdrv_.o hybipt.o vecfcn.o errjac.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) $(LDADD) -lm

ibmdpdr_: $(MINPACK) ibmdpdr_.o machar.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) $(LDADD) -lm

lmddrvc: $(MINPACK) lmddrv.o lmdipt.o ssqfcn.o ssqjac.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) $(LDADD) -lm

lmfdrvc: $(MINPACK) lmfdrv.o lmdipt.o ssqfcn.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) $(LDADD) -lm

lmsdrvc: $(MINPACK) lmsdrv.o lmdipt.o ssqfcn.o ssqjac.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) $(LDADD) -lm

hyjdrvc: $(MINPACK) hyjdrv.o hybipt.o vecfcn.o vecjac.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) $(LDADD) -lm

hybdrvc: $(MINPACK) hybdrv.o hybipt.o vecfcn.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) $(LDADD) -lm

chkdrvc: $(MINPACK) chkdrv.o hybipt.o vecfcn.o errjac.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) $(LDADD) -lm

ibmdpdrc: $(MINPACK) ibmdpdr.o machar.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) $(LDADD) -lm

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

ibmdpdr: $(MINPACK) ibmdpdr.f machar.f
	$(F77) $(FFLAGS) -o $@ $^ $(MINPACK)
