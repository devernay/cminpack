PACKAGE=cminpack
VERSION=1.3.6

CC=gcc
CFLAGS= -O3 -g -Wall -Wextra

### The default configuration is to compile the double precision version

### configuration for the LAPACK/BLAS (double precision) version:
## make LIBSUFFIX= CFLAGS="-O3 -g -Wall -Wextra -D__cminpack_float__"
#LIBSUFFIX=s
#CFLAGS="-O3 -g -Wall -Wextra -DUSE_CBLAS -DUSE_LAPACK"
CFLAGS_L=$(CFLAGS) -DUSE_CBLAS -DUSE_LAPACK
LDADD_L=-framework Accelerate

### configuration for the long double version:
## make LIBSUFFIX=s CFLAGS="-O3 -g -Wall -Wextra -D__cminpack_long_double__"
#LIBSUFFIX=s
#CFLAGS="-O3 -g -Wall -Wextra -D__cminpack_long_double__"
CFLAGS_LD=$(CFLAGS) -D__cminpack_long_double__

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

RANLIB=ranlib

LIB=libcminpack$(LIBSUFFIX).a

OBJS = \
src/$(LIBSUFFIX)chkder.o  src/$(LIBSUFFIX)enorm.o   src/$(LIBSUFFIX)hybrd1.o  src/$(LIBSUFFIX)hybrj.o  \
src/$(LIBSUFFIX)lmdif1.o  src/$(LIBSUFFIX)lmstr1.o  src/$(LIBSUFFIX)qrfac.o   src/$(LIBSUFFIX)r1updt.o \
src/$(LIBSUFFIX)dogleg.o  src/$(LIBSUFFIX)fdjac1.o  src/$(LIBSUFFIX)hybrd.o   src/$(LIBSUFFIX)lmder1.o \
src/$(LIBSUFFIX)lmdif.o   src/$(LIBSUFFIX)lmstr.o   src/$(LIBSUFFIX)qrsolv.o  src/$(LIBSUFFIX)rwupdt.o \
src/$(LIBSUFFIX)dpmpar.o  src/$(LIBSUFFIX)fdjac2.o  src/$(LIBSUFFIX)hybrj1.o  src/$(LIBSUFFIX)lmder.o \
src/$(LIBSUFFIX)lmpar.o   src/$(LIBSUFFIX)qform.o   src/$(LIBSUFFIX)r1mpyq.o  src/$(LIBSUFFIX)covar.o \
src/$(LIBSUFFIX)covar1.o \
legacy_minpack/$(LIBSUFFIX)chkder_.o legacy_minpack/$(LIBSUFFIX)enorm_.o  \
legacy_minpack/$(LIBSUFFIX)hybrd1_.o legacy_minpack/$(LIBSUFFIX)hybrj_.o  \
legacy_minpack/$(LIBSUFFIX)lmdif1_.o legacy_minpack/$(LIBSUFFIX)lmstr1_.o \
legacy_minpack/$(LIBSUFFIX)qrfac_.o  legacy_minpack/$(LIBSUFFIX)r1updt_.o \
legacy_minpack/$(LIBSUFFIX)dogleg_.o legacy_minpack/$(LIBSUFFIX)fdjac1_.o \
legacy_minpack/$(LIBSUFFIX)hybrd_.o  legacy_minpack/$(LIBSUFFIX)lmder1_.o \
legacy_minpack/$(LIBSUFFIX)lmdif_.o  legacy_minpack/$(LIBSUFFIX)lmstr_.o  \
legacy_minpack/$(LIBSUFFIX)qrsolv_.o legacy_minpack/$(LIBSUFFIX)rwupdt_.o \
legacy_minpack/$(LIBSUFFIX)dpmpar_.o legacy_minpack/$(LIBSUFFIX)fdjac2_.o \
legacy_minpack/$(LIBSUFFIX)hybrj1_.o legacy_minpack/$(LIBSUFFIX)lmder_.o  \
legacy_minpack/$(LIBSUFFIX)lmpar_.o  legacy_minpack/$(LIBSUFFIX)qform_.o  \
legacy_minpack/$(LIBSUFFIX)r1mpyq_.o legacy_minpack/$(LIBSUFFIX)covar_.o

# target dir for install
DESTDIR=/usr/local
#
#  Static library target
#

all: $(LIB)

double:
	$(MAKE) LIBSUFFIX=

lapack:
	$(MAKE) LIBSUFFIX=l CFLAGS="$(CFLAGS_L)" LDADD="$(LDADD_L)"

longdouble:
	$(MAKE) LIBSUFFIX=ld CFLAGS="$(CFLAGS_LD)"

float:
	$(MAKE) LIBSUFFIX=s CFLAGS="$(CFLAGS_F)"

half:
	$(MAKE) LIBSUFFIX=h CFLAGS="$(CFLAGS_H)" LDADD="$(LDADD_H)" CC="$(CC_H)"

fortran cuda:
	$(MAKE) -C $@

check checkdouble checklapack checklongdouble checkfloat checkhalf checkfail:
	$(MAKE) -C examples $@

$(LIB):  $(OBJS)
	$(AR) r $@ $(OBJS); $(RANLIB) $@

$(LIBSUFFIX)%.o: %.c
	${CC} ${CFLAGS} -c -o $@ $<

install: $(LIB)
	cp $(LIB) ${DESTDIR}/lib
	chmod 644 ${DESTDIR}/lib/$(LIB)
	$(RANLIB) -t ${DESTDIR}/lib/$(LIB) # might be unnecessary
	cp minpack.h ${DESTDIR}/include
	chmod 644 ${DESTDIR}/include/minpack.h
	cp cminpack.h ${DESTDIR}/include
	chmod 644 ${DESTDIR}/include/cminpack.h

clean:
	rm -f $(OBJS) $(LIB)
	$(MAKE) -C examples clean
	$(MAKE) -C fortran clean

veryclean: clean
	rm -f *.o libcminpack*.a *.gcno *.gcda *~ #*#
	rm -f legacy_minpack/*.o libcminpack*.a legacy_minpack/*.gcno legacy_minpack/*.gcda *~ #*#
	rm -f src/*.o libcminpack*.a src/*.gcno src/*.gcda *~ #*#
	$(MAKE) -C examples veryclean
	$(MAKE) -C examples veryclean LIBSUFFIX=s
	$(MAKE) -C examples veryclean LIBSUFFIX=h
	$(MAKE) -C examples veryclean LIBSUFFIX=l
	$(MAKE) -C fortran veryclean

.PHONY: dist all double lapack longdouble float half fortran cuda check checkhalf checkfail clean veryclean

# COPYFILE_DISABLE=true and COPY_EXTENDED_ATTRIBUTES_DISABLE=true are used to disable inclusion
# of file attributes (._* files) in the tar file on MacOSX
dist:
	mkdir $(PACKAGE)-$(VERSION)
	env COPYFILE_DISABLE=true COPY_EXTENDED_ATTRIBUTES_DISABLE=true tar --exclude-from dist-exclude --exclude $(PACKAGE)-$(VERSION) -cf - . | (cd $(PACKAGE)-$(VERSION); tar xf -)
	tar zcvf $(PACKAGE)-$(VERSION).tar.gz $(PACKAGE)-$(VERSION)
	rm -rf $(PACKAGE)-$(VERSION)
