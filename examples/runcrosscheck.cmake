# Cross-check invoked by ctest (INFORMATIONAL). Runs the pure-C driver (PURE)
# and the f2c driver (F2C) on INPUT and reports, via driver_check.py, how their
# convergence compares.
#
# This is NOT a pass/fail gate on agreement. pure C (src/*.c) is a hand-cleaned-
# up rewrite of the f2c output (src/f2c/*.c); the cleanup regrouped some
# expressions, so a compiler may contract "a*b + c" into a fused multiply-add at
# different places in the two. On the deliberately-extreme More/Garbow/Hillstrom
# driver problems that one-ULP difference can send the runs to different results
# -- exactly as different compilers do, and as cminpack does versus the original
# FORTRAN MINPACK (see
# README.md, "Numerical differences from FORTRAN MINPACK"). Those differences
# are expected and harmless. The test therefore fails ONLY if a driver crashes
# or the comparison tool errors; convergence differences are reported for
# information. The build-failing gates are the standard example tests and the
# driver smoke tests.

# Multi-config generators (e.g. Visual Studio): resolve the INTDIR placeholder.
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

# driver_check.py exit codes: 0 = equivalent, 1 = a hard problem differs
# (informational here), 2 = tooling/parse error (a real failure).
execute_process(COMMAND "${PYTHON}" "${DRIVER_CHECK}" "${OUT_PURE}" "${OUT_F2C}"
  RESULT_VARIABLE R_CMP)
if("${R_CMP}" STREQUAL "2")
  message(FATAL_ERROR "driver_check.py could not compare ${OUT_PURE} and ${OUT_F2C}")
endif()

message("cross-check ran (informational): pure-C and f2c compared for ${OUT_PURE}")
