# Runs cmpfiles (CMPFILES) with its default tolerances on small output/reference
# pairs written to DIR: a changed small residual must differ (exit 1), rounding
# noise around zero must not (exit 0).
set(_cases
  "small_residual_changed|1|final l2 norm 1.192636e-08|final l2 norm 9.0e-05"
  "below_noise_floor|0|final l2 norm 3.3e-17|final l2 norm 0.0")
foreach(_c IN LISTS _cases)
  string(REPLACE "|" ";" _f "${_c}")
  list(GET _f 0 _name)
  list(GET _f 1 _want)
  list(GET _f 2 _ref)
  list(GET _f 3 _out)
  file(WRITE "${DIR}/cmpfiles_${_name}.ref" "${_ref}\n")
  file(WRITE "${DIR}/cmpfiles_${_name}.out" "${_out}\n")
  execute_process(COMMAND "${CMPFILES}"
    "${DIR}/cmpfiles_${_name}.out" "${DIR}/cmpfiles_${_name}.ref"
    RESULT_VARIABLE _ret OUTPUT_QUIET ERROR_QUIET)
  if(NOT "${_ret}" STREQUAL "${_want}")
    message(FATAL_ERROR "${_name}: '${_ref}' vs '${_out}' -> exit ${_ret}, expected ${_want}")
  endif()
endforeach()
