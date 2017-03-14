PACKAGE=cminpack
VERSION=1.3.7-dev

### Use compiler flags from environment if available, else set defaults:
ifneq ($(origin CC), environment)
CC = gcc
endif

ifneq ($(origin CFLAGS), environment)
CFLAGS = -O3 -g -Wall -Wextra
endif

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

SRCS = $(wildcard *.c)
OBJS = $(addprefix $(LIBSUFFIX),$(SRCS:.c=.o))

# target dir for install
DESTDIR ?= /usr/local
#
#  Static library target
#

all: libcminpack$(LIBSUFFIX).a

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

fortran:
	$(MAKE) -C fortran

cuda:
	$(MAKE) -C cuda

check:
	$(MAKE) -C examples check

checkdouble:
	$(MAKE) -C examples checkdouble

checklapack:
	$(MAKE) -C examples checklapack

checklongdouble:
	$(MAKE) -C examples checklongdouble

checkfloat:
	$(MAKE) -C examples checkfloat

checkhalf:
	$(MAKE) -C examples checkhalf

checkfail:
	$(MAKE) -C examples checkfail


libcminpack$(LIBSUFFIX).a:  $(OBJS)
	ar r $@ $(OBJS); ranlib $@

$(LIBSUFFIX)%.o: %.c
	${CC} ${CFLAGS} -c -o $@ $<

install: libcminpack$(LIBSUFFIX).a
	cp libcminpack$(LIBSUFFIX).a ${DESTDIR}/lib
	chmod 644 ${DESTDIR}/lib/libcminpack$(LIBSUFFIX).a
	ranlib -t ${DESTDIR}/lib/libcminpack$(LIBSUFFIX).a # might be unnecessary
	cp minpack.h ${DESTDIR}/include
	chmod 644 ${DESTDIR}/include/minpack.h
	cp cminpack.h ${DESTDIR}/include
	chmod 644 ${DESTDIR}/include/cminpack.h

clean:
	rm -f $(OBJS) libcminpack$(LIBSUFFIX).a
	make -C examples clean
	make -C fortran clean

veryclean: clean
	rm -f *.o libcminpack*.a *.gcno *.gcda *~ #*#
	make -C examples veryclean
	make -C examples veryclean LIBSUFFIX=s
	make -C examples veryclean LIBSUFFIX=h
	make -C examples veryclean LIBSUFFIX=l
	make -C fortran veryclean

.PHONY: dist all double lapack longdouble float half fortran cuda check checkhalf checkfail clean veryclean

# COPYFILE_DISABLE=true and COPY_EXTENDED_ATTRIBUTES_DISABLE=true are used to disable inclusion
# of file attributes (._* files) in the tar file on MacOSX
dist:
	mkdir $(PACKAGE)-$(VERSION)
	env COPYFILE_DISABLE=true COPY_EXTENDED_ATTRIBUTES_DISABLE=true tar --exclude-from dist-exclude --exclude $(PACKAGE)-$(VERSION) -cf - . | (cd $(PACKAGE)-$(VERSION); tar xf -)
	tar zcvf $(PACKAGE)-$(VERSION).tar.gz $(PACKAGE)-$(VERSION)
	rm -rf $(PACKAGE)-$(VERSION)
