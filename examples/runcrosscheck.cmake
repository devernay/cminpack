# Cross-check invoked by ctest. Runs a cminpack driver (DRIVER) on INPUT and
# requires, via driver_check.py --reference-gate, that it converges on every
# problem the committed FORTRAN reference (REF) converges on. Where FORTRAN
# itself does not converge -- the deliberately-extreme More/Garbow/Hillstrom
# problems pushed from 10x/100x starting points -- a different result is
# accepted (both are valid). This is a real regression gate: cminpack must not
# fail a problem the original FORTRAN MINPACK solves.
#
# The FORTRAN reference is checked in (examples/ref/*.fortran.ref), so NO Fortran
# compiler is needed to run this. The test fails on a driver crash, a gate
# violation, or a tooling error. Iteration-count / last-digit differences, and
# differences on problems FORTRAN does not converge, are ignored (see README.md,
# "Numerical differences from FORTRAN MINPACK").

# Multi-config generators (e.g. Visual Studio): resolve the INTDIR placeholder.
if(NOT "${INTDIR}" STREQUAL ".")
  string(REPLACE "${INTDIR}" "$ENV{CMAKE_CONFIG_TYPE}" DRIVER "${DRIVER}")
endif()

execute_process(COMMAND ${CMAKE_CROSSCOMPILING_EMULATOR} ${DRIVER}
  INPUT_FILE "${INPUT}" OUTPUT_FILE "${OUTPUT}" ERROR_FILE "${OUTPUT}.err"
  RESULT_VARIABLE R_RUN)
if(NOT "${R_RUN}" STREQUAL "0")
  message(FATAL_ERROR "driver ${DRIVER} exited with status ${R_RUN}")
endif()

# Build the comma-separated exclusion list from EXCLUDE_FILE (user-maintained
# list of coin-flip problems; see examples/crosscheck_exclude.txt).
set(_excl "")
if(EXCLUDE_FILE AND EXISTS "${EXCLUDE_FILE}")
  file(STRINGS "${EXCLUDE_FILE}" _lines)
  foreach(_l IN LISTS _lines)
    string(REGEX REPLACE "#.*" "" _l "${_l}")
    string(STRIP "${_l}" _l)
    if(_l)
      if(_excl)
        set(_excl "${_excl},${_l}")
      else()
        set(_excl "${_l}")
      endif()
    endif()
  endforeach()
endif()

# driver_check.py --reference-gate exit codes: 0 = converges wherever the
# reference does, 1 = a reference-converged problem is not converged here (a
# real regression), 2 = tooling/parse error.
set(_cmd "${PYTHON}" "${DRIVER_CHECK}" --reference-gate)
if(_excl)
  list(APPEND _cmd --exclude "${_excl}")
endif()
list(APPEND _cmd "${REF}" "${OUTPUT}")
execute_process(COMMAND ${_cmd} RESULT_VARIABLE R_GATE)
if(NOT "${R_GATE}" STREQUAL "0")
  message(FATAL_ERROR
    "${DRIVER} does not converge on a problem the FORTRAN reference "
    "(${REF}) converges on -- see the report above")
endif()

message("cross-check OK: ${DRIVER} converges wherever FORTRAN does")
