CC=cc
#CC=gcc
CFLAGS= -O2 -I.
OBJS = \
chkder.o  enorm.o   hybrd1.o  hybrj.o   lmdif1.o  lmstr1.o  qrfac.o   r1updt.o\
dogleg.o  fdjac1.o  hybrd.o   lmder1.o  lmdif.o   lmstr.o   qrsolv.o  rwupdt.o\
dpmpar.o  fdjac2.o  hybrj1.o  lmder.o   lmpar.o   qform.o   r1mpyq.o

# target dir for install
DESTDIR=/usr/local
#
#  Static library target
#

libminpack.a:  $(OBJS)
	ar r $@ $(OBJS); ranlib $@

%.o: %.c
	${CC} ${CFLAGS} -c -o $@ $<

install: libminpack.a
	cp libminpack.a ${DESTDIR}/lib
	chmod 644 ${DESTDIR}/lib/libminpack.a
	ranlib -t ${DESTDIR}/lib/libminpack.a # might be unnecessary
	cp minpack.h ${DESTDIR}/include
	chmod 644 ${DESTDIR}/include/minpack.h
clean:
	rm -f *.o libminpack.a
