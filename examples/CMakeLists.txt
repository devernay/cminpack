option (BUILD_EXAMPLES "Build the examples and tests." ON)
option (BUILD_EXAMPLES_FORTRAN "Build the FORTRAN examples and tests." OFF)

if (BUILD_EXAMPLES AND TARGET cminpack)
  # Make sure the compiler can find include files from our cminpack library.
  include_directories (${CMINPACK_SOURCE_DIR})

  # Make sure the linker can find the cminpack library once it is built.
  link_directories (${CMINPACK_BINARY_DIR})

  set (FPGM tchkder thybrd thybrd1 thybrj thybrj1 tlmder tlmder1 tlmdif tlmdif1 tlmstr tlmstr1)
  set (CPGM ${FPGM} tfdjac2)

  # cmpfiles is used by runtest.cmake to compare test output to the
  # reference with a numeric tolerance (portable across platforms).
  add_executable (cmpfiles cmpfiles.c)

  # inspired by http://www.netlib.org/clapack/clapack-3.2.1-CMAKE/TESTING/CMakeLists.txt
  # except that here we have to compare the output to a reference
  macro(add_minpack_test testname reference)
    set(TEST_OUTPUT "${CMINPACK_BINARY_DIR}/examples/${testname}.out")
    set(TEST_REFERENCE "${CMINPACK_SOURCE_DIR}/examples/ref/${reference}.ref")
    add_test(NAME ${testname} COMMAND "${CMAKE_COMMAND}"
      -DCMAKE_CROSSCOMPILING_EMULATOR=${CMAKE_CROSSCOMPILING_EMULATOR}
      -DTEST=$<TARGET_FILE:${testname}>
      -DOUTPUT=${TEST_OUTPUT} 
      -DREFERENCE=${TEST_REFERENCE} 
      -DCMPFILES=$<TARGET_FILE:cmpfiles>
      -DINTDIR=${CMAKE_CFG_INTDIR}
      -P "${CMINPACK_SOURCE_DIR}/examples/runtest.cmake")
  endmacro()

  foreach (source ${CPGM})
    add_executable (${source}_ ${source}_.c)
    target_link_libraries (${source}_ cminpack)
    add_minpack_test(${source}_ ${source}c)
    add_executable (${source}c ${source}c.c)
    target_link_libraries (${source}c cminpack)
    add_minpack_test(${source}c ${source}c)
  endforeach()

  # Self-checking regression test for issue #76 (non-convergence / NaN in
  # lmpar on a zero-residual problem). It needs no reference file: it exits
  # non-zero if any solver fails to converge, produces NaN/Inf, or runs up to
  # maxfev. See tissue76.c.
  add_executable (tissue76 tissue76.c)
  target_link_libraries (tissue76 cminpack)
  add_test (NAME tissue76
    COMMAND ${CMAKE_CROSSCOMPILING_EMULATOR} $<TARGET_FILE:tissue76>)

  # C++ wrapper test (issue #74): accept lambdas/functors/std::function.
  # Only built if a C++ compiler is available, since cminpack itself is C.
  include (CheckLanguage)
  check_language (CXX)
  if (CMAKE_CXX_COMPILER)
    enable_language (CXX)
    add_executable (tcppwrap tcppwrap.cpp)
    target_link_libraries (tcppwrap cminpack)
    add_test (NAME tcppwrap
      COMMAND ${CMAKE_CROSSCOMPILING_EMULATOR} $<TARGET_FILE:tcppwrap>)
  endif ()
endif ()

if (BUILD_EXAMPLES_FORTRAN AND TARGET cminpack)
  enable_language (Fortran OPTIONAL)
  if (CMAKE_Fortran_COMPILER_WORKS)
    foreach (FPGM_item ${FPGM})
      list (APPEND FPGM_SRCS ${CMINPACK_SOURCE_DIR}/examples/${FPGM_item}.f)
    endforeach (FPGM_item ${FPGM})

    # add the executable that will do the generation
    add_executable (genf77tests genf77tests.c)
    add_custom_command (
      OUTPUT ${FPGM_SRCS}
      COMMAND $<TARGET_FILE:genf77tests> ${CMINPACK_SOURCE_DIR}/docs/minpack-documentation.txt
      WORKING_DIRECTORY ${CMINPACK_SOURCE_DIR}/examples
      DEPENDS genf77tests
      MAIN_DEPENDENCY ${CMINPACK_SOURCE_DIR}/docs/minpack-documentation.txt
      )
    foreach (source ${FPGM})
      add_executable (${source} ${source}.f)
      target_link_libraries (${source} cminpack)
      add_minpack_test(${source} ${source})
    endforeach ()
  endif ()
endif ()

