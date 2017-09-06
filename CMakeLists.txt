# The name of our project is "CMINPACK". CMakeLists files in this project can
# refer to the root source directory of the project as ${CMINPACK_SOURCE_DIR} and
# to the root binary directory of the project as ${CMINPACK_BINARY_DIR}.
cmake_minimum_required (VERSION 2.6)
project (CMINPACK)
string(TOLOWER ${PROJECT_NAME} PROJECT_NAME_LOWER)

include(${PROJECT_SOURCE_DIR}/cmake/cminpack_utils.cmake)
# Set version and OS-specific settings
set(CMINPACK_VERSION 1.3.6 CACHE STRING "CMinpack version")
set(CMINPACK_SOVERSION 1 CACHE STRING "CMinpack API version")
DISSECT_VERSION()
GET_OS_INFO()

# Add an "uninstall" target
CONFIGURE_FILE ("${PROJECT_SOURCE_DIR}/cmake/uninstall_target.cmake.in"
    "${PROJECT_BINARY_DIR}/uninstall_target.cmake" IMMEDIATE @ONLY)
ADD_CUSTOM_TARGET (uninstall "${CMAKE_COMMAND}" -P
    "${PROJECT_BINARY_DIR}/uninstall_target.cmake")

enable_testing()

if (OS_LINUX OR ${CMAKE_SYSTEM_NAME} STREQUAL "FreeBSD")
  option (USE_FPIC "Use the -fPIC compiler flag." ON)
else (OS_LINUX)
  option (USE_FPIC "Use the -fPIC compiler flag." OFF)
endif (OS_LINUX)

option (BUILD_SHARED_LIBS "Build shared libraries instead of static." OFF)
if (BUILD_SHARED_LIBS)
  message (STATUS "Building shared libraries.")
else ()
  message (STATUS "Building static libraries.")
  set(CMAKE_RELEASE_POSTFIX _s)
  set(CMAKE_RELWITHDEBINFO_POSTFIX _s)
  set(CMAKE_DEBUG_POSTFIX _s)
  set(CMAKE_MINSIZEREL_POSTFIX _s)
  if(WIN32)
    add_definitions(-DCMINPACK_NO_DLL)
  endif(WIN32)
endif ()

option(USE_BLAS "Compile cminpack using a blas library if possible" ON)

#set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/../build)

if(NOT "${CMAKE_PREFIX_PATH}" STREQUAL "")
  include_directories(${CMAKE_PREFIX_PATH}/include)
endif()

set (cminpack_srcs
  src/cminpack.h src/cminpackP.h
  src/chkder.c   src/enorm.c   src/hybrd1.c  src/hybrj.c   src/lmdif1.c  src/lmstr1.c  src/qrfac.c   src/r1updt.c
  src/dogleg.c   src/fdjac1.c  src/hybrd.c   src/lmder1.c  src/lmdif.c   src/lmstr.c   src/qrsolv.c  src/rwupdt.c
  src/dpmpar.c   src/fdjac2.c  src/hybrj1.c  src/lmder.c   src/lmpar.c   src/qform.c   src/r1mpyq.c  src/covar.c
  src/covar1.c
  legacy_minpack/minpack.h
  legacy_minpack/chkder_.c legacy_minpack/enorm_.c  legacy_minpack/hybrd1_.c legacy_minpack/hybrj_.c
  legacy_minpack/lmdif1_.c legacy_minpack/lmstr1_.c legacy_minpack/qrfac_.c  legacy_minpack/r1updt_.c
  legacy_minpack/dogleg_.c legacy_minpack/fdjac1_.c legacy_minpack/hybrd_.c  legacy_minpack/lmder1_.c
  legacy_minpack/lmdif_.c  legacy_minpack/lmstr_.c  legacy_minpack/qrsolv_.c legacy_minpack/rwupdt_.c
  legacy_minpack/dpmpar_.c legacy_minpack/fdjac2_.c legacy_minpack/hybrj1_.c legacy_minpack/lmder_.c
  legacy_minpack/lmpar_.c  legacy_minpack/qform_.c  legacy_minpack/r1mpyq_.c legacy_minpack/covar_.c
  )
set (cminpack_hdrs
    src/cminpack.h legacy_minpack/minpack.h)

add_library (cminpack ${cminpack_srcs})

if (${CMAKE_SYSTEM_NAME} STREQUAL "FreeBSD")
  TARGET_LINK_LIBRARIES(cminpack m)
endif()

# Link with a BLAS library if requested
if (USE_BLAS)
  if (NOT BUILD_SHARED_LIBS)
    set(BLA_STATIC True)
  endif()
  find_package(BLAS)
  if (BLAS_FOUND)
    target_link_libraries(cminpack PUBLIC ${BLAS_LIBRARIES})
    set_target_properties(cminpack PROPERTIES LINK_FLAGS "${BLAS_LINKER_FLAGS}")
    target_compile_definitions(cminpack PUBLIC -DUSE_CBLAS)
  endif()
endif()

install (TARGETS cminpack
   LIBRARY DESTINATION ${CMINPACK_LIB_INSTALL_DIR} COMPONENT library
   ARCHIVE DESTINATION ${CMINPACK_LIB_INSTALL_DIR} COMPONENT library
   RUNTIME DESTINATION bin COMPONENT library)
install (FILES ${cminpack_hdrs} DESTINATION ${CMINPACK_INCLUDE_INSTALL_DIR}
    COMPONENT cminpack_hdrs)

if (USE_FPIC AND NOT BUILD_SHARED_LIBS)
  set_target_properties (cminpack PROPERTIES COMPILE_FLAGS -fPIC)
endif ()

set_target_properties(cminpack PROPERTIES VERSION ${CMINPACK_VERSION} SOVERSION ${CMINPACK_SOVERSION})

add_subdirectory (cmake)
add_subdirectory (examples)