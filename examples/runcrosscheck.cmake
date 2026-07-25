# Cross-check invoked by ctest: the pure-C driver (PURE) and the f2c driver
# (F2C) must reach EQUIVALENT results at double precision -- every problem's
# final residual norm must agree. They need NOT be byte-identical: pure C
# (src/*.c) and f2c (src/f2c/*.c) are two independent source trees, so a
# compiler may contract "a*b + c" into a fused multiply-add at different places
# in each, sending the ill-conditioned problems down different (but equally
# convergent) iteration paths. Both drivers are run on INPUT and compared with
# driver_check.py, which ignores iteration-count and last-digit noise and flags
# only a genuine convergence difference. See README.md, "Numerical differences
# from FORTRAN MINPACK", and examples/crosscheck.py for the full Makefile run.

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

execute_process(COMMAND "${PYTHON}" "${DRIVER_CHECK}" "${OUT_PURE}" "${OUT_F2C}"
  RESULT_VARIABLE R_CMP)
if(NOT "${R_CMP}" STREQUAL "0")
  message(FATAL_ERROR
    "pure-C and f2c reached a materially different result at double precision "
    "-- this is a real bug (${PURE} vs ${F2C}; see above)")
endif()

message("cross-check OK: pure-C and f2c converge equivalently (${OUT_PURE})")
