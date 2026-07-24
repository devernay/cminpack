# Cross-check invoked by ctest: the pure-C driver (PURE) and the f2c driver
# (F2C) MUST produce byte-identical output at double precision -- the one strict
# cminpack invariant. Both are run on INPUT, then compared with compare.py using
# an exact tolerance (rtol=0 atol=0 int-tol=0). A difference fails the test.
#
# This mirrors the strict half of examples/crosscheck.py, but built and driven
# entirely by CMake so it runs under ctest. Differences against the original
# FORTRAN MINPACK, and single/extended-precision iteration-path differences, are
# expected and are NOT checked here (see README.md, "Numerical differences from
# FORTRAN MINPACK", and examples/crosscheck.py for the full Makefile-side run).

# Multi-config generators (e.g. Visual Studio): resolve the INTDIR placeholder,
# as runtest.cmake does.
if(NOT "${INTDIR}" STREQUAL ".")
  string(REPLACE "${INTDIR}" "$ENV{CMAKE_CONFIG_TYPE}" PURE "${PURE}")
  string(REPLACE "${INTDIR}" "$ENV{CMAKE_CONFIG_TYPE}" F2C "${F2C}")
endif()

execute_process(COMMAND ${CMAKE_CROSSCOMPILING_EMULATOR} ${PURE}
  INPUT_FILE "${INPUT}" OUTPUT_FILE "${OUT_PURE}" ERROR_FILE "${OUT_PURE}.err"
  RESULT_VARIABLE R_PURE)
if(NOT "${R_PURE}" STREQUAL "0")
  message(FATAL_ERROR "pure-C driver ${PURE} exited with status ${R_PURE}")
endif()

execute_process(COMMAND ${CMAKE_CROSSCOMPILING_EMULATOR} ${F2C}
  INPUT_FILE "${INPUT}" OUTPUT_FILE "${OUT_F2C}" ERROR_FILE "${OUT_F2C}.err"
  RESULT_VARIABLE R_F2C)
if(NOT "${R_F2C}" STREQUAL "0")
  message(FATAL_ERROR "f2c driver ${F2C} exited with status ${R_F2C}")
endif()

execute_process(COMMAND "${PYTHON}" "${COMPARE}"
  --rtol 0 --atol 0 --int-tol 0 "${OUT_PURE}" "${OUT_F2C}"
  RESULT_VARIABLE R_CMP)
if(NOT "${R_CMP}" STREQUAL "0")
  message(FATAL_ERROR
    "pure-C and f2c output differ at double precision -- this is a real bug "
    "(${PURE} vs ${F2C}; see the diff above)")
endif()

message("cross-check OK: pure-C == f2c (${OUT_PURE})")
