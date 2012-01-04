CC=nvcc
CFLAGS= --compile  -I..

SOURCES = $(wildcard *.cu)
OBJECTS = $(SOURCES:.cu=.o)

all: libminpack.a

libminpack.a:  $(OBJECTS)
	@echo " **************************** "
	@echo " making compiling and library for testing "
	@echo " (cuda needs __device__ source code compiled inline with kernel) "
	@echo " **************************** "
	ar r $@ $(OBJECTS); ranlib $@

%.o: %.cu
	@echo " compiling for testing "
	${CC} ${CFLAGS}  -o $@ $<

clean:
	rm -f *.o libminpack.a *~ #*#

