#!/usr/bin/make 

# pick up your FORTRAN compiler
#F77=g77
F77=gfortran
FFLAGS=-O3
# uncomment the following for FORTRAN MINPACK
#MINPACK=-lminpack
#F77C=$(F77)
#F77CFLAGS=-g

OBJS = \
chkder.o  enorm.o   hybrd1.o  hybrj.o   lmdif1.o  lmstr1.o  qrfac.o   r1updt.o \
dogleg.o  fdjac1.o  hybrd.o   lmder1.o  lmdif.o   lmstr.o   qrsolv.o  rwupdt.o \
dpmpar.o  fdjac2.o  hybrj1.o  lmder.o   lmpar.o   qform.o   r1mpyq.o  covar.o

# target dir for install
DESTDIR=/usr/local
#
#  Static library target
#
.PHONY: all clean veryclean install

all: libminpack.a

libminpack.a:  $(OBJS)
	ar r $@ $(OBJS); ranlib $@

%.o: %.f
	${F77} ${FFLAGS} -c -o $@ $<

install: libminpack.a
	cp libminpack.a ${DESTDIR}/lib
	chmod 644 ${DESTDIR}/lib/libminpack.a
	ranlib -t ${DESTDIR}/lib/libminpack.a # might be unnecessary

clean:
	rm -f $(OBJS) libminpack.a

veryclean:
	rm -f *.o *.gcno *.gcda *~ #*#
