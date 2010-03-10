# The name of our project is "CMINPACK". CMakeLists files in this project can
# refer to the root source directory of the project as ${CMINPACK_SOURCE_DIR} and
# to the root binary directory of the project as ${CMINPACK_BINARY_DIR}.
cmake_minimum_required (VERSION 2.6)
project (CMINPACK)

#set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/../build)

# Recurse into the "Hello" and "Demo" subdirectories. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.
add_subdirectory (examples)
#add_subdirectory (Demo) 

# Create a library called "cminpack" which includes the source file "hello.cxx".
# The extension is already found. Any number of sources could be listed here.
add_library(cminpack chkder.c  enorm.c   hybrd1.c  hybrj.c   lmdif1.c  lmstr1.c  qrfac.c   r1updt.c dogleg.c  fdjac1.c  hybrd.c   lmder1.c  lmdif.c   lmstr.c   qrsolv.c  rwupdt.c dpmpar.c  fdjac2.c  hybrj1.c  lmder.c   lmpar.c   qform.c   r1mpyq.c  covar.c chkder_.c enorm_.c  hybrd1_.c hybrj_.c  lmdif1_.c lmstr1_.c qrfac_.c  r1updt_.c dogleg_.c fdjac1_.c hybrd_.c  lmder1_.c lmdif_.c  lmstr_.c  qrsolv_.c rwupdt_.c dpmpar_.c fdjac2_.c hybrj1_.c lmder_.c  lmpar_.c  qform_.c  r1mpyq_.c covar_.c)
