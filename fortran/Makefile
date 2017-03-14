#!/usr/bin/make 

# pick up your FORTRAN compiler
#F77=g77

ifneq ($(origin F77), environment)
F77 = gfortran
endif

ifneq ($(origin FFLAGS), environment)
FFLAGS = -O3
endif

# uncomment the following for FORTRAN MINPACK
#MINPACK=-lminpack
#F77C=$(F77)
#F77CFLAGS=-g

LIB  = libminpack.a
SRCS = $(wildcard *.f)
OBJS = $(SRCS:.f=.o)

# target dir for install
DESTDIR ?= /usr/local
#
#  Static library target
#
.PHONY: all clean veryclean install

all: $(LIB)

$(LIB): $(OBJS)
	$(AR) r $@ $^
	ranlib $@

%.o: %.f
	$(F77) $(FFLAGS) -c -o $@ $<

install: $(LIB)
	cp $(LIB) $(DESTDIR)/lib
	chmod 644 $(DESTDIR)/lib/$(LIB)
	ranlib -t $(DESTDIR)/lib/$(LIB) # might be unnecessary

clean:
	$(RM) -f $(OBJS) $(LIB)

veryclean:
	$(RM) -f *.o *.gcno *.gcda *~ #*#
