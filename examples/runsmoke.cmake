# Smoke test for the intensive driver programs (lmddrvc, lmfdrvc, lmsdrvc,
# hyjdrvc, hybdrvc, chkdrvc). They run cminpack on the difficult test problems
# from More, Garbow and Hillstrom, reading the problem set from stdin (INPUT).
#
# Their detailed numeric output is platform- and compiler-dependent (iteration
# counts, termination codes and borderline results differ from machine to
# machine, as the documentation warns), so it is NOT compared against a
# reference. Instead this checks that the driver runs the whole problem set to
# completion, exits successfully, and produces no NaN -- catching crashes,
# hangs (via the ctest timeout) and NaN propagation regressions portably.

# Handle multi-config generators (e.g. Visual Studio): replace the INTDIR
# placeholder in the test path with the actual configuration, as runtest.cmake
# does.
if(NOT "${INTDIR}" STREQUAL ".")
  string(REPLACE "${INTDIR}" "$ENV{CMAKE_CONFIG_TYPE}" TEST "${TEST}")
endif()

if(INPUT AND NOT "${INPUT}" STREQUAL "")
  execute_process(COMMAND ${CMAKE_CROSSCOMPILING_EMULATOR} ${TEST}
    INPUT_FILE ${INPUT}
    OUTPUT_FILE ${OUTPUT}
    ERROR_FILE ${OUTPUT}.err
    RESULT_VARIABLE RET)
else()
  execute_process(COMMAND ${CMAKE_CROSSCOMPILING_EMULATOR} ${TEST}
    OUTPUT_FILE ${OUTPUT}
    ERROR_FILE ${OUTPUT}.err
    RESULT_VARIABLE RET)
endif()

if(NOT "${RET}" STREQUAL "0")
  message(FATAL_ERROR "Driver ${TEST} exited with status ${RET}")
endif()

file(READ "${OUTPUT}" _out)
string(TOLOWER "${_out}" _out)
# match a standalone nan token (optionally signed), not a substring of a word
if(_out MATCHES "(^|[^a-z])[-+]?nan([^a-z]|$)")
  message(FATAL_ERROR "Driver ${TEST} produced NaN in its output")
endif()

message("Driver ${TEST} ran to completion (no crash, no NaN)")
