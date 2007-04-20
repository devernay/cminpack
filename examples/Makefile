#!/usr/bin/make 

# pick up your FORTRAN compiler
#F77=g77
F77=gfortran

# uncomment the following for FORTRAN MINPACK
#MINPACK=-lminpack
#F77C=$(F77)

# uncomment the following for C MINPACK
MINPACK=../libminpack.a
F77C=gcc -I..

# Execute this script to extract the examples from the MINPACK
# documentation (/usr/share/doc/minpack-devel-19961126/minpack-documentation.txt)
# and compile and run them.  You may compare the numerical results
# with those in the same file.

PGM=tchkder thybrd thybrd1 thybrj thybrj1 tlmder tlmder1 tlmdif \
	tlmdif1 tlmstr tlmstr1
CPGM=tcchkder tchybrd tchybrd1 tchybrj tchybrj1 tclmder tclmder1 tclmdif \
	tclmdif1 tclmstr tclmstr1
SRC=tchkder.f thybrd.f thybrd1.f thybrj.f thybrj1.f tlmder.f tlmder1.f \
	tlmdif.f tlmdif1.f tlmstr.f tlmstr1.f
CSRC=tchkder.c thybrd.c thybrd1.c thybrj.c thybrj1.c tlmder.c tlmder1.c \
	tlmdif.c tlmdif1.c tlmstr.c tlmstr1.c

all: ctest ftest

ctest: ${CPGM}
	for x in ${CPGM}; do echo $$x; ./$$x; done

ftest: ${PGM} 
	for x in ${PGM}; do echo $$x; ./$$x; done

clean:
	-rm $(PGM) $(CPGM) $(SRC)

${SRC}: /usr/share/doc/minpack-devel-19961126/minpack-documentation.txt.gz
	zcat $< | awk  '					\
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

tcchkder: tchkder.c					
	$(F77C) -o $@ $< $(MINPACK) -lm
tchybrd: thybrd.c					
	$(F77C) -o $@ $< $(MINPACK) -lm
tchybrd1: thybrd1.c					
	$(F77C) -o $@ $< $(MINPACK) -lm
tchybrj: thybrj.c					
	$(F77C) -o $@ $< $(MINPACK) -lm
tchybrj1: thybrj1.c					
	$(F77C) -o $@ $< $(MINPACK) -lm
tclmder: tlmder.c					
	$(F77C) -o $@ $< $(MINPACK) -lm
tclmder1: tlmder1.c					
	$(F77C) -o $@ $< $(MINPACK) -lm
tclmdif: tlmdif.c
	$(F77C) -o $@ $< $(MINPACK) -lm
tclmdif1: tlmdif1.c
	$(F77C) -o $@ $< $(MINPACK) -lm
tclmstr: tlmstr.c
	$(F77C) -o $@ $< $(MINPACK) -lm
tclmstr1: tlmstr1.c
	$(F77C) -o $@ $< $(MINPACK) -lm
