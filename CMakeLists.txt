# The name of our project is "CMINPACK". CMakeLists files in this project can
# refer to the root source directory of the project as ${CMINPACK_SOURCE_DIR} and
# to the root binary directory of the project as ${CMINPACK_BINARY_DIR}.
cmake_minimum_required (VERSION 2.6)
project (CMINPACK)

string (REGEX MATCH "Linux" OS_LINUX ${CMAKE_SYSTEM_NAME})
string (REGEX MATCH "BSD" OS_BSD ${CMAKE_SYSTEM_NAME})
if (WIN32)
  SET (OS_WIN TRUE BOOL INTERNAL)
endif (WIN32)
if (OS_LINUX)
  message (STATUS "Operating system is Linux")
elseif (OS_BSD)
  message (STATUS "Operating system is BSD")
elseif (OS_WIN)
  message (STATUS "Operating system is Windows")
else (OS_LINUX)
  message (STATUS "Operating system is generic Unix")
endif (OS_LINUX)

enable_testing()

if (OS_LINUX)
  option (USE_FPIC "Use the -fPIC compiler flag." ON)
else (OS_LINUX)
  option (USE_FPIC "Use the -fPIC compiler flag." OFF)
endif (OS_LINUX)

option (SHARED_LIBS "Build shared libraries instead of static." OFF)
if (SHARED_LIBS)
  message (STATUS "Building shared libraries.")
  set (LIB_TYPE SHARED)
else (SHARED_LIBS)
  message (STATUS "Building static libraries.")
  set (LIB_TYPE STATIC)
endif (SHARED_LIBS)

#set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/../build)

# Recurse into the "Hello" and "Demo" subdirectories. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.
add_subdirectory (examples)
#add_subdirectory (Demo) 

set (cminpack_srcs
  cminpack.h
  chkder.c  enorm.c   hybrd1.c  hybrj.c   lmdif1.c  lmstr1.c  qrfac.c   r1updt.c
  dogleg.c  fdjac1.c  hybrd.c   lmder1.c  lmdif.c   lmstr.c   qrsolv.c  rwupdt.c
  dpmpar.c  fdjac2.c  hybrj1.c  lmder.c   lmpar.c   qform.c   r1mpyq.c  covar.c
  minpack.h
  chkder_.c enorm_.c  hybrd1_.c hybrj_.c  lmdif1_.c lmstr1_.c qrfac_.c  r1updt_.c
  dogleg_.c fdjac1_.c hybrd_.c  lmder1_.c lmdif_.c  lmstr_.c  qrsolv_.c rwupdt_.c
  dpmpar_.c fdjac2_.c hybrj1_.c lmder_.c  lmpar_.c  qform_.c  r1mpyq_.c covar_.c
  )

# Create a library called "cminpack" which includes the source file "hello.cxx".
# The extension is already found. Any number of sources could be listed here.
add_library (cminpack ${LIB_TYPE} ${cminpack_srcs})

if (USE_FPIC AND NOT SHARED_LIBS)
  set_target_properties (cminpack PROPERTIES COMPILE_FLAGS -fPIC)
endif (USE_FPIC AND NOT SHARED_LIBS)

