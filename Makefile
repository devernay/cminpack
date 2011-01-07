PACKAGE=cminpack
VERSION=1.1.3

#CC=cc
CC=gcc
CFLAGS= -O2 -g -Wall -Wextra -I.
#CFLAGS= -g -Wall -I.
OBJS = \
chkder.o  enorm.o   hybrd1.o  hybrj.o   lmdif1.o  lmstr1.o  qrfac.o   r1updt.o \
dogleg.o  fdjac1.o  hybrd.o   lmder1.o  lmdif.o   lmstr.o   qrsolv.o  rwupdt.o \
dpmpar.o  fdjac2.o  hybrj1.o  lmder.o   lmpar.o   qform.o   r1mpyq.o  covar.o covar1.o \
chkder_.o enorm_.o  hybrd1_.o hybrj_.o  lmdif1_.o lmstr1_.o qrfac_.o  r1updt_.o \
dogleg_.o fdjac1_.o hybrd_.o  lmder1_.o lmdif_.o  lmstr_.o  qrsolv_.o rwupdt_.o \
dpmpar_.o fdjac2_.o hybrj1_.o lmder_.o  lmpar_.o  qform_.o  r1mpyq_.o covar_.o

# target dir for install
DESTDIR=/usr/local
#
#  Static library target
#

all: libminpack.a

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
	cp cminpack.h ${DESTDIR}/include
	chmod 644 ${DESTDIR}/include/cminpack.h

clean:
	rm -f *.o libminpack.a *~ #*#

.PHONY: dist

# COPYFILE_DISABLE=true and COPY_EXTENDED_ATTRIBUTES_DISABLE=true are used to disable inclusion
# of file attributes (._* files) in the tar file on MacOSX
dist:
	mkdir $(PACKAGE)-$(VERSION)
	env COPYFILE_DISABLE=true COPY_EXTENDED_ATTRIBUTES_DISABLE=true tar --exclude-from dist-exclude --exclude $(PACKAGE)-$(VERSION) -cf - . | (cd $(PACKAGE)-$(VERSION); tar xf -)
	tar zcvf $(PACKAGE)-$(VERSION).tar.gz $(PACKAGE)-$(VERSION)
	rm -rf $(PACKAGE)-$(VERSION)
