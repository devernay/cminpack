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

FPGM=\
$(LIBSUFFIX)tchkder \
$(LIBSUFFIX)thybrd \
$(LIBSUFFIX)thybrd1 \
$(LIBSUFFIX)thybrj \
$(LIBSUFFIX)thybrj1 \
$(LIBSUFFIX)tlmder \
$(LIBSUFFIX)tlmder1 \
$(LIBSUFFIX)tlmdif \
$(LIBSUFFIX)tlmdif1 \
$(LIBSUFFIX)tlmstr \
$(LIBSUFFIX)tlmstr1 \
$(LIBSUFFIX)ibmdpdr

PGM=\
$(LIBSUFFIX)tchkder_ \
$(LIBSUFFIX)thybrd_ \
$(LIBSUFFIX)thybrd1_ \
$(LIBSUFFIX)thybrj_ \
$(LIBSUFFIX)thybrj1_ \
$(LIBSUFFIX)tlmder_ \
$(LIBSUFFIX)tlmder1_ \
$(LIBSUFFIX)tlmdif_ \
$(LIBSUFFIX)tlmdif1_ \
$(LIBSUFFIX)tlmstr_ \
$(LIBSUFFIX)tlmstr1_ \
$(LIBSUFFIX)tfdjac2_ \
$(LIBSUFFIX)ibmdpdr_

CPGM=\
$(LIBSUFFIX)tchkderc \
$(LIBSUFFIX)thybrdc \
$(LIBSUFFIX)thybrd1c \
$(LIBSUFFIX)thybrjc \
$(LIBSUFFIX)thybrj1c \
$(LIBSUFFIX)tlmderc \
$(LIBSUFFIX)tlmder1c \
$(LIBSUFFIX)tlmdifc \
$(LIBSUFFIX)tlmdif1c \
$(LIBSUFFIX)tlmstrc \
$(LIBSUFFIX)tlmstr1c \
$(LIBSUFFIX)tfdjac2c \
$(LIBSUFFIX)ibmdpdrc \
$(LIBSUFFIX)tchkderc_box \
$(LIBSUFFIX)tlmderc_box \
$(LIBSUFFIX)thybrjc_box

FDRVPGM=\
$(LIBSUFFIX)lmddrv_ \
$(LIBSUFFIX)lmfdrv_ \
$(LIBSUFFIX)lmsdrv_ \
$(LIBSUFFIX)hyjdrv_ \
$(LIBSUFFIX)hybdrv_ \
$(LIBSUFFIX)chkdrv_

DRVPGM=\
$(LIBSUFFIX)lmddrv \
$(LIBSUFFIX)lmfdrv \
$(LIBSUFFIX)lmsdrv \
$(LIBSUFFIX)hyjdrv \
$(LIBSUFFIX)hybdrv \
$(LIBSUFFIX)chkdrv

CDRVPGM=\
$(LIBSUFFIX)lmddrvc \
$(LIBSUFFIX)lmfdrvc \
$(LIBSUFFIX)lmsdrvc \
$(LIBSUFFIX)hyjdrvc \
$(LIBSUFFIX)hybdrvc \
$(LIBSUFFIX)chkdrvc

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
	for x in $(PGM); do echo $$x; $(VALGRIND) ./$$x > $$x.out; diff -u ref/`echo $$x.ref|sed -e s/_.ref/c.ref/` $$x.out ; done

testdrv: $(MINPACK) $(FDRVPGM) testlmdrv testhydrv testchkdrv

testlmdrv: $(LIBSUFFIX)lmddrv_ $(LIBSUFFIX)lmfdrv_ $(LIBSUFFIX)lmsdrv_
	@echo "*** Running LM tests (calling MINPACK from C using $(MINPACK))"
	for x in $^; do echo $$x; $(VALGRIND) ./$$x < testdata/lm.data > $$x.out; diff -u ref/`echo $$x |sed -e s/_/c/`.ref $$x.out ; done

testhydrv: $(LIBSUFFIX)hyjdrv_ $(LIBSUFFIX)hybdrv_
	@echo "*** Running HY tests (calling MINPACK from C using $(MINPACK))"
	for x in $^; do echo $$x; $(VALGRIND) ./$$x < testdata/hybrd.data > $$x.out; diff -u ref/`echo $$x |sed -e s/_/c/`.ref $$x.out ; done

testchkdrv: $(LIBSUFFIX)chkdrv_
	@echo "*** Running CHK tests (calling MINPACK from C using $(MINPACK))"
	for x in $^; do echo $$x; $(VALGRIND) ./$$x < testdata/chkder.data > $$x.out; diff -u ref/`echo $$x |sed -e s/_/c/`.ref $$x.out ; done

ctest: $(MINPACK) $(CPGM)
	@echo "*** Running standard tests (calling CMINPACK from C using $(MINPACK))"
	for x in $(CPGM); do echo $$x; $(VALGRIND) ./$$x > $$x.out; diff -u ref/$$x.ref $$x.out ; done

ctestdrv: $(MINPACK) $(CDRVPGM) ctestlmdrv ctesthydrv ctestchkdrv

ctestlmdrv: $(LIBSUFFIX)lmddrvc $(LIBSUFFIX)lmfdrvc $(LIBSUFFIX)lmsdrvc
	@echo "*** Running LM tests (calling CMINPACK from C using $(MINPACK))"
	for x in $^; do echo $$x; $(VALGRIND) ./$$x < testdata/lm.data > $$x.out; diff -u ref/$$x.ref $$x.out ; done

ctesthydrv: $(LIBSUFFIX)hyjdrvc $(LIBSUFFIX)hybdrvc
	@echo "*** Running HY tests (calling CMINPACK from C using $(MINPACK))"
	for x in $^; do echo $$x; $(VALGRIND) ./$$x < testdata/hybrd.data > $$x.out; diff -u ref/$$x.ref $$x.out ; done

ctestchkdrv: $(LIBSUFFIX)chkdrvc
	@echo "*** Running CHK tests (calling CMINPACK from C using $(MINPACK))"
	for x in $^; do echo $$x; $(VALGRIND) ./$$x < testdata/chkder.data > $$x.out; diff -u ref/$$x.ref $$x.out ; done

ftest: $(MINPACK) $(FPGM)
	@echo "*** Running standard tests (calling MINPACK from FORTRAN using $(MINPACK))"
	for x in ${FPGM}; do echo $$x; $(VALGRIND) ./$$x > $$x.out; diff -u ref/$$x.ref $$x.out ; done

ftestdrv: $(MINPACK) $(DRVPGM) ftestlmdrv ftesthydrv ftestchkdrv

ftestlmdrv: $(LIBSUFFIX)lmddrv $(LIBSUFFIX)lmfdrv $(LIBSUFFIX)lmsdrv
	@echo "*** Running LM tests (calling MINPACK from FORTRAN using $(MINPACK))"
	for x in $^; do echo $$x; $(VALGRIND) ./$$x < testdata/lm.data > $$x.out; diff -u ref/$$x.ref $$x.out ; done

ftesthydrv: $(LIBSUFFIX)hyjdrv $(LIBSUFFIX)hybdrv
	@echo "*** Running HY tests (calling MINPACK from FORTRAN using $(MINPACK))"
	for x in $^; do echo $$x; $(VALGRIND) ./$$x < testdata/hybrd.data > $$x.out; diff -u ref/$$x.ref $$x.out ; done

ftestchkdrv: $(LIBSUFFIX)chkdrv
	@echo "*** Running CHK tests (calling MINPACK from FORTRAN using $(MINPACK))"
	for x in $^; do echo $$x; $(VALGRIND) ./$$x < testdata/chkder.data > $$x.out; diff -u ref/$$x.ref $$x.out ; done

check: checkdouble checkfloat

checkdoublec:
	$(MAKE) -C .. double 2>&1 > /dev/null
	$(MAKE) clean LIBSUFFIX= 2>&1 > /dev/null
	-$(MAKE) ctest LIBSUFFIX=
	$(MAKE) clean 2>&1 > /dev/null

checkdouble: checkdoublec ../fortran/libminpack.a
	$(MAKE) clean LIBSUFFIX= 2>&1 > /dev/null
	-$(MAKE) test LIBSUFFIX=
	$(MAKE) clean LIBSUFFIX= MINPACK=../fortran/libminpack.a 2>&1 > /dev/null
	-$(MAKE) test LIBSUFFIX= MINPACK=../fortran/libminpack.a
	$(MAKE) clean LIBSUFFIX= 2>&1 > /dev/null
	-$(MAKE) ftest LIBSUFFIX=
	$(MAKE) clean LIBSUFFIX= MINPACK=../fortran/libminpack.a 2>&1 > /dev/null
	-$(MAKE) ftest LIBSUFFIX= MINPACK=../fortran/libminpack.a
	$(MAKE) clean LIBSUFFIX= 2>&1 > /dev/null
	-$(MAKE) -k testdrv LIBSUFFIX=
	$(MAKE) clean LIBSUFFIX= 2>&1 > /dev/null
	-$(MAKE) -k ctestdrv LIBSUFFIX=
	$(MAKE) clean 2>&1 > /dev/null

checkfloatc:
	$(MAKE) -C .. float 2>&1 > /dev/null
	$(MAKE) clean LIBSUFFIX=s CFLAGS="$(CFLAGS_F)" 2>&1 > /dev/null
	-$(MAKE) ctest LIBSUFFIX=s CFLAGS="$(CFLAGS_F)" 
	$(MAKE) clean LIBSUFFIX=s 2>&1 > /dev/null

checkfloat: checkfloatc
	$(MAKE) clean LIBSUFFIX=s CFLAGS="$(CFLAGS_F)" 2>&1 > /dev/null
	-$(MAKE) test LIBSUFFIX=s CFLAGS="$(CFLAGS_F)"
	$(MAKE) clean LIBSUFFIX=s CFLAGS="$(CFLAGS_F)" 2>&1 > /dev/null
	-$(MAKE) -k testdrv LIBSUFFIX=s CFLAGS="$(CFLAGS_F)" 
	$(MAKE) clean LIBSUFFIX=s CFLAGS="$(CFLAGS_F)" 2>&1 > /dev/null
	-$(MAKE) -k ctestdrv LIBSUFFIX=s CFLAGS="$(CFLAGS_F)"
	$(MAKE) clean LIBSUFFIX=s 2>&1 > /dev/null

checklapack:
	$(MAKE) -C .. lapack 2>&1 > /dev/null
	$(MAKE) clean LIBSUFFIX=l CFLAGS="$(CFLAGS_L)" LDADD="$(LDADD_L)" 2>&1 > /dev/null
	-$(MAKE) ctest LIBSUFFIX=l CFLAGS="$(CFLAGS_L)" LDADD="$(LDADD_L)"
	$(MAKE) clean LIBSUFFIX=l CFLAGS="$(CFLAGS_L)" LDADD="$(LDADD_L)" 2>&1 > /dev/null
	-$(MAKE) -k ctestdrv LIBSUFFIX=l CFLAGS="$(CFLAGS_L)" LDADD="$(LDADD_L)"

checkhalf:
	$(MAKE) -C .. half 2>&1 > /dev/null
	$(MAKE) clean LIBSUFFIX=h CFLAGS="$(CFLAGS_H)" LDADD="$(LDADD_H)" CC="$(CC_H)" 2>&1 > /dev/null
	-$(MAKE) test LIBSUFFIX=h CFLAGS="$(CFLAGS_H)" LDADD="$(LDADD_H)" CC="$(CC_H)"
	$(MAKE) clean LIBSUFFIX=h CFLAGS="$(CFLAGS_H)" LDADD="$(LDADD_H)" CC="$(CC_H)" 2>&1 > /dev/null
	-$(MAKE) ctest LIBSUFFIX=h CFLAGS="$(CFLAGS_H)" LDADD="$(LDADD_H)" CC="$(CC_H)"
	$(MAKE) clean LIBSUFFIX=h CFLAGS="$(CFLAGS_H)" LDADD="$(LDADD_H)" CC="$(CC_H)" 2>&1 > /dev/null
	-$(MAKE) -k testdrv LIBSUFFIX=h CFLAGS="$(CFLAGS_H)" LDADD="$(LDADD_H)" CC="$(CC_H)"
	$(MAKE) clean LIBSUFFIX=h CFLAGS="$(CFLAGS_H)" LDADD="$(LDADD_H)" CC="$(CC_H)" 2>&1 > /dev/null
	-$(MAKE) -k ctestdrv LIBSUFFIX=h CFLAGS="$(CFLAGS_H)" LDADD="$(LDADD_H)" CC="$(CC_H)"
#	$(MAKE) clean LIBSUFFIX=h 2>&1 > /dev/null
#	$(MAKE) -C .. clean LIBSUFFIX=h 2>&1 > /dev/null

checkfail: ../fortran/libminpack.a check
	$(MAKE) clean LIBSUFFIX= MINPACK=../fortran/libminpack.a 2>&1 > /dev/null
	-$(MAKE) -k testdrv LIBSUFFIX= MINPACK=../fortran/libminpack.a
	$(MAKE) clean LIBSUFFIX= 2>&1 > /dev/null
	-$(MAKE) -k ftestdrv LIBSUFFIX=
	$(MAKE) clean LIBSUFFIX= MINPACK=../fortran/libminpack.a 2>&1 > /dev/null
	-$(MAKE) -k ftestdrv LIBSUFFIX= MINPACK=../fortran/libminpack.a
	$(MAKE) clean LIBSUFFIX= 2>&1 > /dev/null

../libcminpack$(LIBSUFFIX).a:
	$(MAKE) -C ..

../fortran/libminpack.a:
	$(MAKE) -C ../fortran

clean:
	-rm -f $(PGM) $(FDRVPGM) $(DRVPGM) $(CPGM) $(CDRVPGM) *.o *.gcno *.gcda $(PGM:=.out) $(CPGM:=.out) $(FDRVPGM:=.out) $(DRVPGM:=.out) $(CDRVPGM:=.out) $(FPGM:=.out) *~ #*#
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

$(LIBSUFFIX)%: %.f $(MINPACK)
	$(F77) -o $@ $< $(MINPACK)

$(LIBSUFFIX)%: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $< $(MINPACK) $(LDADD) -lm

$(LIBSUFFIX)%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

$(LIBSUFFIX)tchkderc_box: tchkderc.c
	$(CC) $(CPPFLAGS) -DBOX_CONSTRAINTS $(CFLAGS) -o $@ $< $(MINPACK) $(LDADD) -lm

$(LIBSUFFIX)tlmderc_box: tlmderc.c
	$(CC) $(CPPFLAGS) -DBOX_CONSTRAINTS $(CFLAGS) -o $@ $< $(MINPACK) $(LDADD) -lm

$(LIBSUFFIX)thybrjc_box: thybrjc.c
	$(CC) $(CPPFLAGS) -DBOX_CONSTRAINTS $(CFLAGS) -o $@ $< $(MINPACK) $(LDADD) -lm

$(LIBSUFFIX)lmddrv_: $(MINPACK) lmddrv_.o lmdipt.o ssqfcn.o ssqjac.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) $(LDADD) -lm

$(LIBSUFFIX)lmfdrv_: $(MINPACK) lmfdrv_.o lmdipt.o ssqfcn.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) $(LDADD) -lm

$(LIBSUFFIX)lmsdrv_: $(MINPACK) lmsdrv_.o lmdipt.o ssqfcn.o ssqjac.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) $(LDADD) -lm

$(LIBSUFFIX)hyjdrv_: $(MINPACK) hyjdrv_.o hybipt.o vecfcn.o vecjac.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) $(LDADD) -lm

$(LIBSUFFIX)hybdrv_: $(MINPACK) hybdrv_.o hybipt.o vecfcn.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) $(LDADD) -lm

$(LIBSUFFIX)chkdrv_: $(MINPACK) chkdrv_.o hybipt.o vecfcn.o errjac.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) $(LDADD) -lm

$(LIBSUFFIX)ibmdpdr_: $(MINPACK) ibmdpdr_.o machar.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) $(LDADD) -lm

$(LIBSUFFIX)lmddrvc: $(MINPACK) lmddrv.o lmdipt.o ssqfcn.o ssqjac.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) $(LDADD) -lm

$(LIBSUFFIX)lmfdrvc: $(MINPACK) lmfdrv.o lmdipt.o ssqfcn.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) $(LDADD) -lm

$(LIBSUFFIX)lmsdrvc: $(MINPACK) lmsdrv.o lmdipt.o ssqfcn.o ssqjac.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) $(LDADD) -lm

$(LIBSUFFIX)hyjdrvc: $(MINPACK) hyjdrv.o hybipt.o vecfcn.o vecjac.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) $(LDADD) -lm

$(LIBSUFFIX)hybdrvc: $(MINPACK) hybdrv.o hybipt.o vecfcn.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) $(LDADD) -lm

$(LIBSUFFIX)chkdrvc: $(MINPACK) chkdrv.o hybipt.o vecfcn.o errjac.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) $(LDADD) -lm

$(LIBSUFFIX)ibmdpdrc: $(MINPACK) ibmdpdr.o machar.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $^ $(MINPACK) $(LDADD) -lm

$(LIBSUFFIX)lmddrv: $(MINPACK) lmddrv.f lmdipt.f ssqfcn.f ssqjac.f
	$(F77) $(FFLAGS) -o $@ $^ $(MINPACK)

$(LIBSUFFIX)lmfdrv: $(MINPACK) lmfdrv.f lmdipt.f ssqfcn.f
	$(F77) $(FFLAGS) -o $@ $^ $(MINPACK)

$(LIBSUFFIX)lmsdrv: $(MINPACK) lmsdrv.f lmdipt.f ssqfcn.f ssqjac.f
	$(F77) $(FFLAGS) -o $@ $^ $(MINPACK)

$(LIBSUFFIX)hyjdrv: $(MINPACK) hyjdrv.f hybipt.f vecfcn.f vecjac.f
	$(F77) $(FFLAGS) -o $@ $^ $(MINPACK)

$(LIBSUFFIX)hybdrv: $(MINPACK) hybdrv.f hybipt.f vecfcn.f
	$(F77) $(FFLAGS) -o $@ $^ $(MINPACK)

$(LIBSUFFIX)chkdrv: $(MINPACK) chkdrv.f hybipt.f vecfcn.f errjac.f
	$(F77) $(FFLAGS) -o $@ $^ $(MINPACK)

$(LIBSUFFIX)ibmdpdr: $(MINPACK) ibmdpdr.f machar.f
	$(F77) $(FFLAGS) -o $@ $^ $(MINPACK)
