option (BUILD_EXAMPLES "Build the examples and tests." ON)
option (BUILD_EXAMPLES_FORTRAN "Build the FORTRAN examples and tests." OFF)

if (BUILD_EXAMPLES AND TARGET cminpack)
  # Make sure the compiler can find include files from our cminpack library.
  include_directories (${CMINPACK_SOURCE_DIR}/include)

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

  # Intensive driver tests: run cminpack on the difficult problems from More,
  # Garbow and Hillstrom, reading the problem set from stdin. Their detailed
  # output is platform/compiler-dependent (iteration counts and termination
  # codes vary), so these are smoke tests -- run to completion, exit 0, no NaN
  # -- rather than reference comparisons (see runsmoke.cmake).
  macro(add_minpack_driver_smoketest testname inputfile)
    add_test(NAME ${testname} COMMAND "${CMAKE_COMMAND}"
      -DCMAKE_CROSSCOMPILING_EMULATOR=${CMAKE_CROSSCOMPILING_EMULATOR}
      -DTEST=$<TARGET_FILE:${testname}>
      -DOUTPUT=${CMINPACK_BINARY_DIR}/examples/${testname}.out
      -DINPUT=${CMINPACK_SOURCE_DIR}/examples/testdata/${inputfile}
      -DINTDIR=${CMAKE_CFG_INTDIR}
      -P "${CMINPACK_SOURCE_DIR}/examples/runsmoke.cmake")
  endmacro()

  add_executable (lmddrvc lmddrv.c ssqipt.c ssqfcn.c ssqjac.c)
  add_executable (lmfdrvc lmfdrv.c ssqipt.c ssqfcn.c)
  add_executable (lmsdrvc lmsdrv.c ssqipt.c ssqfcn.c ssqjac.c)
  add_executable (hyjdrvc hyjdrv.c neqipt.c vecfcn.c vecjac.c)
  add_executable (hybdrvc hybdrv.c neqipt.c vecfcn.c)
  add_executable (chkdrvc chkdrv.c neqipt.c vecfcn.c errjac.c)
  foreach (drv lmddrvc lmfdrvc lmsdrvc hyjdrvc hybdrvc chkdrvc)
    target_link_libraries (${drv} cminpack)
  endforeach ()
  add_minpack_driver_smoketest (lmddrvc ssq.data)
  add_minpack_driver_smoketest (lmfdrvc ssq.data)
  add_minpack_driver_smoketest (lmsdrvc ssq.data)
  add_minpack_driver_smoketest (hyjdrvc neq.data)
  add_minpack_driver_smoketest (hybdrvc neq.data)
  add_minpack_driver_smoketest (chkdrvc chkder.data)

  # Cross-check: run each cminpack driver (pure-C and f2c) and require, via
  # examples/driver_check.py, that it converges on every problem the committed
  # FORTRAN reference (examples/ref/*.fortran.ref) converges on. Where FORTRAN
  # itself does not converge -- the deliberately-extreme problems -- a different
  # result is accepted (both are valid). This is a real regression gate against
  # the original FORTRAN MINPACK, and needs NO Fortran compiler (the reference is
  # checked in). Iteration-count / last-digit differences, and differences on
  # problems FORTRAN does not converge, are ignored (see README.md, "Numerical
  # differences from FORTRAN MINPACK"). Needs Python 3; if none is found the
  # cross-check tests are skipped (the rest of the build is unaffected).
  option (CMINPACK_CROSSCHECK "Add the FORTRAN-reference cross-check tests (needs Python 3)." ON)
  find_package (Python3 COMPONENTS Interpreter QUIET)
  if (CMINPACK_CROSSCHECK AND Python3_Interpreter_FOUND)
    add_executable (lmddrv_ lmddrv_.c ssqipt.c ssqfcn.c ssqjac.c)
    add_executable (lmfdrv_ lmfdrv_.c ssqipt.c ssqfcn.c)
    add_executable (lmsdrv_ lmsdrv_.c ssqipt.c ssqfcn.c ssqjac.c)
    add_executable (hyjdrv_ hyjdrv_.c neqipt.c vecfcn.c vecjac.c)
    add_executable (hybdrv_ hybdrv_.c neqipt.c vecfcn.c)
    add_executable (chkdrv_ chkdrv_.c neqipt.c vecfcn.c errjac.c)
    foreach (drv lmddrv_ lmfdrv_ lmsdrv_ hyjdrv_ hybdrv_ chkdrv_)
      target_link_libraries (${drv} cminpack)
    endforeach ()

    # Extra exclusions applied ONLY to BLAS/LAPACK builds. Those use a different
    # numeric path (cblas_dnrm2 in enorm, LAPACK qrfac/qrsolv) that can flip a
    # different set of extreme coin-flip problems than the plain build, and the
    # exact set depends on the BLAS/LAPACK implementation. Kept separate from
    # crosscheck_exclude.txt so the plain-build gate stays strict on problems
    # that only BLAS/LAPACK flips.
    set (CROSSCHECK_BLAS_EXCLUDE "")
    if (USE_BLAS OR USE_LAPACK)
      set (CROSSCHECK_BLAS_EXCLUDE
        -DEXCLUDE_FILE_BLAS=${CMINPACK_SOURCE_DIR}/examples/crosscheck_exclude_blas.txt)
    endif ()

    # name : driver target : input data : committed FORTRAN reference
    macro (add_crosscheck name exe inputfile ref)
      add_test (NAME ${name} COMMAND "${CMAKE_COMMAND}"
        -DCMAKE_CROSSCOMPILING_EMULATOR=${CMAKE_CROSSCOMPILING_EMULATOR}
        -DDRIVER=$<TARGET_FILE:${exe}>
        -DINPUT=${CMINPACK_SOURCE_DIR}/examples/testdata/${inputfile}
        -DOUTPUT=${CMINPACK_BINARY_DIR}/examples/${name}.out
        -DREF=${CMINPACK_SOURCE_DIR}/examples/ref/${ref}
        -DPYTHON=${Python3_EXECUTABLE}
        -DDRIVER_CHECK=${CMINPACK_SOURCE_DIR}/examples/driver_check.py
        -DEXCLUDE_FILE=${CMINPACK_SOURCE_DIR}/examples/crosscheck_exclude.txt
        ${CROSSCHECK_BLAS_EXCLUDE}
        -DINTDIR=${CMAKE_CFG_INTDIR}
        -P "${CMINPACK_SOURCE_DIR}/examples/runcrosscheck.cmake")
    endmacro ()

    # Gate both the pure-C (...c) and f2c (..._) builds against the FORTRAN ref.
    add_crosscheck (crosscheck_lmddrvc lmddrvc ssq.data    lmddrv.fortran.ref)
    add_crosscheck (crosscheck_lmddrv_ lmddrv_ ssq.data    lmddrv.fortran.ref)
    add_crosscheck (crosscheck_lmfdrvc lmfdrvc ssq.data    lmfdrv.fortran.ref)
    add_crosscheck (crosscheck_lmfdrv_ lmfdrv_ ssq.data    lmfdrv.fortran.ref)
    add_crosscheck (crosscheck_lmsdrvc lmsdrvc ssq.data    lmsdrv.fortran.ref)
    add_crosscheck (crosscheck_lmsdrv_ lmsdrv_ ssq.data    lmsdrv.fortran.ref)
    add_crosscheck (crosscheck_hyjdrvc hyjdrvc neq.data    hyjdrv.fortran.ref)
    add_crosscheck (crosscheck_hyjdrv_ hyjdrv_ neq.data    hyjdrv.fortran.ref)
    add_crosscheck (crosscheck_hybdrvc hybdrvc neq.data    hybdrv.fortran.ref)
    add_crosscheck (crosscheck_hybdrv_ hybdrv_ neq.data    hybdrv.fortran.ref)
    add_crosscheck (crosscheck_chkdrvc chkdrvc chkder.data chkdrv.fortran.ref)
    add_crosscheck (crosscheck_chkdrv_ chkdrv_ chkder.data chkdrv.fortran.ref)
    message (STATUS "cminpack tests: FORTRAN-reference cross-check ENABLED (Python 3 found).")
  elseif (CMINPACK_CROSSCHECK)
    # Python 3 is NOT required to build cminpack, nor for the standard tests
    # (which use the C tool cmpfiles) or the driver smoke tests (pure CMake). It
    # is only needed for the FORTRAN-reference cross-check, so we skip just that
    # part and let everything else run.
    message (WARNING
      "Python 3 was not found: the FORTRAN-reference cross-check tests are "
      "SKIPPED. The build and the rest of the test suite (standard tests + "
      "driver smoke tests) are unaffected, but `ctest` coverage is INCOMPLETE "
      "-- it will not verify that cminpack converges wherever FORTRAN MINPACK "
      "does. Install Python 3 (or use the Makefile's `make check`) for full "
      "coverage.")
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

