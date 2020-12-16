# The name of our project is "CMINPACK". CMakeLists files in this project can
# refer to the root source directory of the project as ${CMINPACK_SOURCE_DIR} and
# to the root binary directory of the project as ${CMINPACK_BINARY_DIR}.
cmake_minimum_required (VERSION 3.0.0)

if (NOT DEFINED CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING "Build type")
endif ()

project (CMINPACK C)
string (TOLOWER ${PROJECT_NAME} PROJECT_NAME_LOWER)

list (APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

include (${PROJECT_SOURCE_DIR}/cmake/cminpack_utils.cmake)
# Set version and OS-specific settings
set (CMINPACK_VERSION 1.3.8 STRING "CMinpack version")
set (CMINPACK_SOVERSION 1 STRING "CMinpack API version")
DISSECT_VERSION ()
GET_OS_INFO ()

# Add an "uninstall" target
configure_file ("${PROJECT_SOURCE_DIR}/cmake/uninstall_target.cmake.in"
  "${PROJECT_BINARY_DIR}/uninstall_target.cmake" IMMEDIATE @ONLY)
add_custom_target (uninstall "${CMAKE_COMMAND}" -P
  "${PROJECT_BINARY_DIR}/uninstall_target.cmake")

enable_testing ()

option (BUILD_SHARED_LIBS "Build shared libraries instead of static." OFF)
option (USE_BLAS "Compile cminpack using cblas library if possible" ON)
set (CMINPACK_PRECISION "all" CACHE STRING "Precision variants to build ('s', 'd', 'ld', 'all')")

#set (CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/../build)

if (NOT "${CMAKE_PREFIX_PATH}" STREQUAL "")
  include_directories (${CMAKE_PREFIX_PATH}/include)
endif ()

set (cminpack_srcs
  cminpack.h cminpackP.h
  chkder.c  enorm.c   hybrd1.c  hybrj.c   lmdif1.c  lmstr1.c  qrfac.c   r1updt.c
  dogleg.c  fdjac1.c  hybrd.c   lmder1.c  lmdif.c   lmstr.c   qrsolv.c  rwupdt.c
  dpmpar.c  fdjac2.c  hybrj1.c  lmder.c   lmpar.c   qform.c   r1mpyq.c  covar.c covar1.c
  minpack.h
  chkder_.c enorm_.c  hybrd1_.c hybrj_.c  lmdif1_.c lmstr1_.c qrfac_.c  r1updt_.c
  dogleg_.c fdjac1_.c hybrd_.c  lmder1_.c lmdif_.c  lmstr_.c  qrsolv_.c rwupdt_.c
  dpmpar_.c fdjac2_.c hybrj1_.c lmder_.c  lmpar_.c  qform_.c  r1mpyq_.c covar_.c
  )
set (cminpack_hdrs
  cminpack.h minpack.h)

include (CheckTypeSize)
check_type_size ("double" SIZEOF_DOUBLE)
check_type_size ("long double" SIZEOF_LONG_DOUBLE)

if (${CMINPACK_PRECISION} STREQUAL s)
  message (STATUS "Building for single precision (float).")
  set (cminpack_libs cminpacks)
elseif (${CMINPACK_PRECISION} STREQUAL d)
  message (STATUS "Building for double precision (double).")
  set (cminpack_libs cminpack)
elseif (${CMINPACK_PRECISION} STREQUAL ld)
  if (${SIZEOF_LONG_DOUBLE} GREATER ${SIZEOF_DOUBLE})
    message (STATUS "Building for extended precision (long double).")
    set (cminpack_libs cminpackld)
  else ()
    # Microsoft Visual C++ for x86 makes long double a synonym for double
    message (FATAL_ERROR "Compiler does not support long double, or long double is a synonym of double")
  endif ()
elseif (${CMINPACK_PRECISION} STREQUAL all)
  if (${SIZEOF_LONG_DOUBLE} GREATER ${SIZEOF_DOUBLE})
    message (STATUS "Building for single precision (float), double precision (double) and extended precision (long double).")
    set (cminpack_libs cminpacks;cminpack;cminpackld)
  else ()
    message (STATUS "Building for single precision (float) and double precision (double).")
    set (cminpack_libs cminpacks;cminpack)
  endif ()
endif ()

if (BUILD_SHARED_LIBS)
  message (STATUS "Building shared libraries.")
  if (NOT DEFINED CMAKE_DEBUG_POSTFIX)
    set (CMAKE_DEBUG_POSTFIX _d)
  endif ()
else ()
  message (STATUS "Building static libraries.")
  if (NOT DEFINED CMAKE_RELEASE_POSTFIX)
    set (CMAKE_RELEASE_POSTFIX _s)
  endif ()
  if (NOT DEFINED CMAKE_RELWITHDEBINFO_POSTFIX)
    set (CMAKE_RELWITHDEBINFO_POSTFIX _s)
  endif ()
  if (NOT DEFINED CMAKE_DEBUG_POSTFIX)
    set (CMAKE_DEBUG_POSTFIX _s_d)
  endif ()
  if (NOT DEFINED CMAKE_MINSIZEREL_POSTFIX)
    set (CMAKE_MINSIZEREL_POSTFIX _s)
  endif ()
endif ()

foreach (cminpack_lib ${cminpack_libs})
  add_library (${cminpack_lib} ${cminpack_srcs})
  if (${cminpack_lib} STREQUAL "cminpacks")
    target_compile_definitions (cminpacks PUBLIC __cminpack_float__)
  elseif (${cminpack_lib} STREQUAL "cminpackld")
    target_compile_definitions (cminpackld PUBLIC __cminpack_long_double__)
  endif ()

  if (NOT BUILD_SHARED_LIBS AND WIN32)
    target_compile_definitions (${cminpack_lib} PUBLIC CMINPACK_NO_DLL)
  endif ()

  target_include_directories (${cminpack_lib} INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMINPACK_INCLUDE_INSTALL_DIR}/>
    )

  if (${CMAKE_SYSTEM_NAME} STREQUAL "FreeBSD")
    target_link_libraries (${cminpack_lib} m)
  endif ()

  include (CheckLibraryExists)
  include (CheckFunctionExists)
  check_function_exists (sqrt HAVE_SQRT_NO_LIBM)
  if (NOT HAVE_SQRT_NO_LIBM)
    check_library_exists ("m" sqrt "" HAVE_LIBM)
    if (HAVE_LIBM)
      target_link_libraries (${cminpack_lib} PUBLIC m)
    endif ()
  endif ()

  # Link with CBLAS library if requested
  if (USE_BLAS)
    if (${cminpack_lib} STREQUAL cminpackld)
      message (WARNING "BLAS cannot be used for the extended precision version")
    else ()
      find_package (CBLAS)
      if (CBLAS_FOUND)
        target_link_libraries (${cminpack_lib} PUBLIC ${CBLAS_LIBRARIES})
        set_target_properties (${cminpack_lib} PROPERTIES LINK_FLAGS "${CBLAS_LINKER_FLAGS}")
        target_compile_definitions (${cminpack_lib} PUBLIC USE_CBLAS)
      endif ()
    endif ()
  endif ()

  install (TARGETS ${cminpack_lib} EXPORT CMinpackConfig
    LIBRARY DESTINATION ${CMINPACK_LIB_INSTALL_DIR} COMPONENT library
    ARCHIVE DESTINATION ${CMINPACK_LIB_INSTALL_DIR} COMPONENT library
    RUNTIME DESTINATION bin COMPONENT library)
  install (FILES ${cminpack_hdrs} DESTINATION ${CMINPACK_INCLUDE_INSTALL_DIR}
    COMPONENT cminpack_hdrs)

  set_target_properties (${cminpack_lib} PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    VERSION ${CMINPACK_VERSION}
    SOVERSION ${CMINPACK_SOVERSION}
    )

  set (pkg_conf_file ${CMAKE_CURRENT_BINARY_DIR}/${cminpack_lib}.pc)
  if (USE_BLAS AND CBLAS_FOUND)
      set (PC_CMINPACK_CFLAGS "-DUSE_CBLAS")
      if (NOT "${BLAS_LIBRARIES}" STREQUAL "")
          string (REPLACE ";" " -l" PC_CMINPACK_LIBRARIES "${BLAS_LIBRARIES}")
          set (PC_CMINPACK_LIBRARIES "-l${PC_CMINPACK_LIBRARIES}")
      endif ()
      set (PC_CMINPACK_LDFLAGS "${BLAS_LINKER_FLAGS}")
  endif ()
  configure_file (cmake/cminpack.pc.in ${pkg_conf_file} @ONLY)
  install (FILES ${pkg_conf_file}
    DESTINATION ${CMINPACK_LIB_INSTALL_DIR}/pkgconfig/ COMPONENT pkgconfig)

endforeach ()

install (EXPORT CMinpackConfig DESTINATION share/cminpack NAMESPACE cminpack::)

add_subdirectory (examples)
