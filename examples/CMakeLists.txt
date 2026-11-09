option (BUILD_EXAMPLES "Build the examples and tests." ON)
option (BUILD_EXAMPLES_FORTRAN "Build the FORTRAN examples and tests." OFF)

if (BUILD_EXAMPLES)
  # the test drivers define their own helpers and macros in each file
  set (CMAKE_UNITY_BUILD OFF)

  # Make sure the compiler can find include files from our cminpack library.
  include_directories (${CMINPACK_SOURCE_DIR}/include)

  # Make sure the linker can find the cminpack library once it is built.
  link_directories (${CMINPACK_BINARY_DIR})

  set (FPGM tchkder thybrd thybrd1 thybrj thybrj1 tlmder tlmder1 tlmdif tlmdif1 tlmstr tlmstr1)
  set (CPGM ${FPGM} tfdjac2 tenorm)

  # cmpfiles is used by runtest.cmake to compare test output to the
  # reference with a numeric tolerance (portable across platforms).
  add_executable (cmpfiles cmpfiles.c)
  # Issue #151: FORTRAN D exponents must be compared as numbers.
  add_test (NAME cmpfiles_dexp
    COMMAND ${CMAKE_CROSSCOMPILING_EMULATOR} $<TARGET_FILE:cmpfiles>
      ${CMINPACK_SOURCE_DIR}/examples/testdata/cmpfiles_d_a.txt
      ${CMINPACK_SOURCE_DIR}/examples/testdata/cmpfiles_d_b.txt)

  # Self-test of cmpfiles' default tolerances (see runcmpfilestest.cmake).
  add_test (NAME cmpfiles_selftest COMMAND "${CMAKE_COMMAND}"
    -DCMPFILES=$<TARGET_FILE:cmpfiles>
    -DDIR=${CMINPACK_BINARY_DIR}/examples
    -P "${CMINPACK_SOURCE_DIR}/examples/runcmpfilestest.cmake")

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

  # The C++ wrapper tests are only built if a C++ compiler is available,
  # since cminpack itself is C.
  include (CheckLanguage)
  check_language (CXX)
  if (CMAKE_CXX_COMPILER)
    enable_language (CXX)
    include (CheckCCompilerFlag)
    check_c_compiler_flag ("-fno-unwind-tables -fno-asynchronous-unwind-tables"
      CMINPACK_C_NO_UNWIND_TABLES)
  endif ()

  # The tests below are defined for each precision built, with the prefix of
  # the Makefile and of the reference files: tlmdifc, stlmdifc, ldtlmdifc.
  foreach (lib ${cminpack_libs})
    string (REPLACE cminpack "" p ${lib})

    foreach (source ${CPGM})
      foreach (test ${source}_ ${source}c)
        add_executable (${p}${test} ${test}.c)
        target_link_libraries (${p}${test} ${lib})
        add_minpack_test(${p}${test} ${p}${source}c)
      endforeach ()
    endforeach()

    # Self-checking regression test for issue #76 (non-convergence / NaN in
    # lmpar on a zero-residual problem). It needs no reference file: it exits
    # non-zero if any solver fails to converge, produces NaN/Inf, or runs up to
    # maxfev. See tissue76.c.
    add_executable (${p}tissue76 tissue76.c)
    target_link_libraries (${p}tissue76 ${lib})
    add_test (NAME ${p}tissue76
      COMMAND ${CMAKE_CROSSCOMPILING_EMULATOR} $<TARGET_FILE:${p}tissue76>)

    # Self-checking regression test for issue #81: a negative value returned by
    # the user function while lmstr computes a Jacobian row must be returned as
    # info (it was reported as 0). See tissue81.c.
    add_executable (${p}tissue81 tissue81.c)
    target_link_libraries (${p}tissue81 ${lib})
    add_test (NAME ${p}tissue81
      COMMAND ${CMAKE_CROSSCOMPILING_EMULATOR} $<TARGET_FILE:${p}tissue81>)

    add_executable (${p}tissue88 tissue88.c)
    target_link_libraries (${p}tissue88 ${lib})
    add_test (NAME ${p}tissue88
      COMMAND ${CMAKE_CROSSCOMPILING_EMULATOR} $<TARGET_FILE:${p}tissue88>)

    add_executable (${p}tissue87 tissue87.c)
    target_link_libraries (${p}tissue87 ${lib})
    add_test (NAME ${p}tissue87
      COMMAND ${CMAKE_CROSSCOMPILING_EMULATOR} $<TARGET_FILE:${p}tissue87>)

    # Issue #84: qrfac must not read past the end of a when m < n.
    add_executable (${p}tissue84 tissue84.c)
    target_link_libraries (${p}tissue84 ${lib})
    add_test (NAME ${p}tissue84
      COMMAND ${CMAKE_CROSSCOMPILING_EMULATOR} $<TARGET_FILE:${p}tissue84>)

    # Issue #85: qrfac with pivot=0 must allow wa to coincide with rdiag.
    add_executable (${p}tissue85 tissue85.c)
    target_link_libraries (${p}tissue85 ${lib})
    add_test (NAME ${p}tissue85
      COMMAND ${CMAKE_CROSSCOMPILING_EMULATOR} $<TARGET_FILE:${p}tissue85>)

    # The low-level routines must be exported from a Windows DLL (issue #91).
    add_executable (${p}tlowlevel tlowlevel.c)
    target_link_libraries (${p}tlowlevel ${lib})
    add_test (NAME ${p}tlowlevel
      COMMAND ${CMAKE_CROSSCOMPILING_EMULATOR} $<TARGET_FILE:${p}tlowlevel>)

    add_executable (${p}tissue89 tissue89.c)
    target_link_libraries (${p}tissue89 ${lib})
    add_test (NAME ${p}tissue89
      COMMAND ${CMAKE_CROSSCOMPILING_EMULATOR} $<TARGET_FILE:${p}tissue89>)

    if (CMAKE_CXX_COMPILER)
      # C++ wrapper test (issue #74): accept lambdas/functors/std::function.
      add_executable (${p}tcppwrap tcppwrap.cpp)
      target_link_libraries (${p}tcppwrap ${lib})
      add_test (NAME ${p}tcppwrap
        COMMAND ${CMAKE_CROSSCOMPILING_EMULATOR} $<TARGET_FILE:${p}tcppwrap>)

      # Issue #94: link this test with a copy of cminpack built without unwind
      # tables, so that an exception unwinding through C calls std::terminate.
      get_target_property (_srcs ${lib} SOURCES)
      list (TRANSFORM _srcs PREPEND ${CMINPACK_SOURCE_DIR}/)
      add_library (${p}cminpack_nounwind STATIC EXCLUDE_FROM_ALL ${_srcs})
      target_include_directories (${p}cminpack_nounwind PUBLIC
        ${CMINPACK_SOURCE_DIR}/include PRIVATE ${CMINPACK_SOURCE_DIR}/src)
      target_compile_definitions (${p}cminpack_nounwind PUBLIC CMINPACK_NO_DLL
        $<TARGET_PROPERTY:${lib},INTERFACE_COMPILE_DEFINITIONS>)
      target_link_libraries (${p}cminpack_nounwind PUBLIC
        $<TARGET_PROPERTY:${lib},LINK_LIBRARIES>)
      if (CMINPACK_C_NO_UNWIND_TABLES)
        target_compile_options (${p}cminpack_nounwind PRIVATE
          -fno-unwind-tables -fno-asynchronous-unwind-tables)
      endif ()
      add_executable (${p}tcppwrap_except tcppwrap_except.cpp)
      target_link_libraries (${p}tcppwrap_except ${p}cminpack_nounwind)
      add_test (NAME ${p}tcppwrap_except
        COMMAND ${CMAKE_CROSSCOMPILING_EMULATOR} $<TARGET_FILE:${p}tcppwrap_except>)
    endif ()

    add_executable (${p}lmddrvc lmddrv.c ssqipt.c ssqfcn.c ssqjac.c)
    add_executable (${p}lmfdrvc lmfdrv.c ssqipt.c ssqfcn.c)
    add_executable (${p}lmsdrvc lmsdrv.c ssqipt.c ssqfcn.c ssqjac.c)
    add_executable (${p}hyjdrvc hyjdrv.c neqipt.c vecfcn.c vecjac.c)
    add_executable (${p}hybdrvc hybdrv.c neqipt.c vecfcn.c)
    add_executable (${p}chkdrvc chkdrv.c neqipt.c vecfcn.c errjac.c)
    foreach (drv lmddrvc lmfdrvc lmsdrvc hyjdrvc hybdrvc chkdrvc)
      target_link_libraries (${p}${drv} ${lib})
    endforeach ()
    add_minpack_driver_smoketest (${p}lmddrvc ssq.data)
    add_minpack_driver_smoketest (${p}lmfdrvc ssq.data)
    add_minpack_driver_smoketest (${p}lmsdrvc ssq.data)
    add_minpack_driver_smoketest (${p}hyjdrvc neq.data)
    add_minpack_driver_smoketest (${p}hybdrvc neq.data)
    add_minpack_driver_smoketest (${p}chkdrvc chkder.data)

    # The box variants, and the IBM driver (both use the ibmdpdrc reference).
    add_executable (${p}tchkderc_box tchkderc_box.c)
    add_executable (${p}tlmderc_box tlmderc_box.c)
    add_executable (${p}thybrjc_box thybrjc_box.c)
    add_executable (${p}ibmdpdrc ibmdpdr.c machar.c)
    add_executable (${p}ibmdpdr_ ibmdpdr_.c machar.c)
    foreach (test tchkderc_box tlmderc_box thybrjc_box ibmdpdrc ibmdpdr_)
      target_link_libraries (${p}${test} ${lib})
      string (REPLACE ibmdpdr_ ibmdpdrc reference ${test})
      add_minpack_test (${p}${test} ${p}${reference})
    endforeach ()
  endforeach ()

  # built from the enorm sources: cminpackld is not built everywhere
  add_executable (tissue82 tissue82.c
    ${CMINPACK_SOURCE_DIR}/src/enorm.c ${CMINPACK_SOURCE_DIR}/src/f2c/enorm_.c)
  target_compile_definitions (tissue82 PRIVATE __cminpack_long_double__ CMINPACK_NO_DLL)
  target_include_directories (tissue82 PRIVATE ${CMINPACK_SOURCE_DIR}/src)
  if (HAVE_LIBM)
    target_link_libraries (tissue82 m)
  endif ()
  add_test (NAME tissue82
    COMMAND ${CMAKE_CROSSCOMPILING_EMULATOR} $<TARGET_FILE:tissue82>)

  # BLAS calls must pass __cminpack_blasint__ (here long), not int: compile only
  if (CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    file (GLOB CMINPACK_C_SOURCES ${CMINPACK_SOURCE_DIR}/src/*.c)
    if (CMAKE_OSX_SYSROOT)
      set (BLASINT_SYSROOT -isysroot ${CMAKE_OSX_SYSROOT})
    endif ()
    add_test (NAME blasint_types
      COMMAND ${CMAKE_C_COMPILER} ${BLASINT_SYSROOT} -fsyntax-only
        -DUSE_BLAS -D__cminpack_blasint__=long -Werror=incompatible-pointer-types
        -I${CMINPACK_SOURCE_DIR}/include -I${CMINPACK_SOURCE_DIR}/src
        ${CMINPACK_C_SOURCES})
  endif ()

  # checks the test helpers with double precision reference values
  if (TARGET cminpack)
    add_executable (tissue95 tissue95.c vecfcn.c vecjac.c errjac.c)
    target_link_libraries (tissue95 cminpack)
    add_test (NAME tissue95
      COMMAND ${CMAKE_CROSSCOMPILING_EMULATOR} $<TARGET_FILE:tissue95>)

    add_executable (tissue96 tissue96.c errjac.c)
    target_link_libraries (tissue96 cminpack)
    add_test (NAME tissue96
      COMMAND ${CMAKE_CROSSCOMPILING_EMULATOR} $<TARGET_FILE:tissue96>)
  endif ()

  if (CMAKE_CXX_COMPILER)
    # Issue #93: a const functor with a non-const operator() must not compile.
    add_library (tcppwrap_const OBJECT EXCLUDE_FROM_ALL tcppwrap_const.cpp)
    add_test (NAME tcppwrap_const
      COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR}
        --target tcppwrap_const --config $<CONFIG>)
    set_tests_properties (tcppwrap_const PROPERTIES WILL_FAIL TRUE)
  endif ()

  # Cross-check: run each cminpack driver (pure-C and f2c) and require, via
  # examples/driver_check.py, that it converges on every problem the committed
  # FORTRAN reference (examples/ref/*.fortran.ref) converges on. Where FORTRAN
  # itself does not converge -- the deliberately-extreme problems -- a different
  # result is accepted (both are valid). This is a real regression gate against
  # the original FORTRAN MINPACK, and needs NO Fortran compiler (the reference is
  # checked in). Iteration-count / last-digit differences, and differences on
  # problems FORTRAN does not converge, are ignored (see README.md, "Numerical
  # differences from FORTRAN MINPACK"). Needs Python 3; if none is found the
  # cross-check tests are skipped (the rest of the build is unaffected). The
  # FORTRAN reference is double precision, so only the double library is checked.
  option (CMINPACK_CROSSCHECK "Add the FORTRAN-reference cross-check tests (needs Python 3)." ON)
  find_package (Python3 COMPONENTS Interpreter QUIET)
  if (CMINPACK_CROSSCHECK AND Python3_Interpreter_FOUND AND TARGET cminpack)
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
    add_test (NAME test_driver_check COMMAND ${Python3_EXECUTABLE}
      "${CMINPACK_SOURCE_DIR}/examples/test_driver_check.py")
    message (STATUS "cminpack tests: FORTRAN-reference cross-check ENABLED (Python 3 found).")
  elseif (CMINPACK_CROSSCHECK AND NOT Python3_Interpreter_FOUND)
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
    # Pure-FORTRAN lmsdrv: lmsdrv.f with FORTRAN MINPACK
    set (_f ${CMINPACK_SOURCE_DIR}/fortran)
    add_executable (lmsdrv lmsdrv.f ssqipt.f ssqfcn.f ssqjac.f ${_f}/lmstr1.f
      ${_f}/lmstr.f ${_f}/lmpar.f ${_f}/qrfac.f ${_f}/qrsolv.f ${_f}/rwupdt.f
      ${_f}/enorm.f ${_f}/dpmpar.f)
    # start locals as NaN, so a value lost between calls cannot pass by luck
    if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options (lmsdrv PRIVATE -finit-real=snan)
    endif ()
    if (COMMAND add_crosscheck)
      add_crosscheck (crosscheck_lmsdrv lmsdrv ssq.data lmsdrv.fortran.ref)
    endif ()
  endif ()
endif ()

