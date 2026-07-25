#!/bin/sh
# crosscheck_matrix.sh -- DIAGNOSTIC (not part of `make check`).
#
# Shows how the intensive-driver convergence set varies with compiler
# optimization. It generates a FORTRAN -O0 reference ON THIS MACHINE, then
# builds the driver programs against FORTRAN MINPACK and against cminpack at
# several optimization levels and, for each, reports (via
# driver_check.py --reference-gate) how many problems the FORTRAN -O0 reference
# converges on but that build does NOT -- and which ones.
#
# The point: on the deliberately-extreme More/Garbow/Hillstrom problems,
# convergence is a last-bit coin-flip that is not stable even for identical
# FORTRAN source across -O0/-O1/-O2/-O3 (see README.md, "Numerical differences
# from FORTRAN MINPACK"). The driver cross-check is therefore a gate with a
# small exclusion list (examples/crosscheck_exclude.txt): this diagnostic shows
# which problems flip on a given compiler, so you can decide what belongs there.
#
# Requires gcc, gfortran and python3. Run it from the examples/ directory:
#     sh ./crosscheck_matrix.sh
#
# Environment overrides: CC (default cc), F77 (default gfortran),
# PYTHON (default python3).

set -e
cd "$(dirname "$0")"

CC=${CC:-cc}
F77=${F77:-gfortran}
PY=${PYTHON:-python3}
DRIVERS="lmddrv lmfdrv lmsdrv hyjdrv hybdrv chkdrv"

data_for() {
    case "$1" in
        lmddrv|lmfdrv|lmsdrv) echo testdata/ssq.data ;;
        hyjdrv|hybdrv)        echo testdata/neq.data ;;
        chkdrv)               echo testdata/chkder.data ;;
    esac
}

command -v "$F77" >/dev/null 2>&1 || { echo "error: $F77 not found (need a Fortran compiler)"; exit 2; }
command -v "$PY"  >/dev/null 2>&1 || { echo "error: $PY not found (need Python 3)"; exit 2; }

TMP=$(mktemp -d)
cleanup() {
    make clean LIBSUFFIX= >/dev/null 2>&1 || true
    make -C .. clean >/dev/null 2>&1 || true
    make -C ../fortran clean >/dev/null 2>&1 || true
    rm -rf "$TMP"
}
trap cleanup EXIT

# Build the six f2c driver programs against $2 (the MINPACK library) with C
# flags $1, then run each on its input; outputs go to $3.<driver>.out
build_and_run() {
    _cflags=$1; _lib=$2; _prefix=$3
    make clean LIBSUFFIX= >/dev/null 2>&1
    _progs=""
    for d in $DRIVERS; do _progs="$_progs ${d}_"; done
    make $_progs LIBSUFFIX= CC="$CC" CFLAGS="$_cflags" MINPACK="$_lib" >/dev/null 2>&1
    for d in $DRIVERS; do
        ./${d}_ < "$(data_for "$d")" > "${_prefix}.${d}.out" 2>&1
    done
}

# Gate a run against the FORTRAN -O0 reference and print one table row.
# driver_check.py --reference-gate prints one line per regressing problem:
#   "  problem <nprob>/<dim> (block N): FORTRAN reference residual ... worse)"
gate_row() {
    _label=$1; _prefix=$2
    _tot=0; _probs=""
    for d in $DRIVERS; do
        _fails=$($PY driver_check.py --reference-gate \
                    "$TMP/ref.${d}.out" "${_prefix}.${d}.out" \
                 | grep -E '^  problem [0-9]+/[0-9]+ \(block' || true)
        [ -z "$_fails" ] && continue
        _n=$(printf '%s\n' "$_fails" | grep -c '(block')
        _tot=$((_tot + _n))
        _t=$(printf '%s\n' "$_fails" \
             | sed -E 's/^  problem ([0-9]+)\/([0-9]+) .*/p\1\/d\2/')
        _probs="$_probs $_t"
    done
    if [ "$_tot" -eq 0 ]; then _probs="-"; else
        _probs=$(printf '%s\n' $_probs | sort | uniq -c | awk '{print $2 (($1>1)?"(x"$1")":"")}' | tr '\n' ' ')
    fi
    printf "%-26s %4d   %s\n" "$_label" "$_tot" "$_probs"
}

echo "Generating FORTRAN -O0 reference on this machine ($($F77 --version | head -1)) ..."
make -C ../fortran clean >/dev/null 2>&1
make -C ../fortran F77="$F77" FFLAGS="-O0 -g" >/dev/null 2>&1
build_and_run "-O0 -g" ../fortran/libminpack.a "$TMP/ref"

echo
printf "%-26s %4s   %s\n" "configuration" "fail" "problems (FORTRAN -O0 converges, build does not)"
printf "%-26s %4s   %s\n" "-------------------------" "----" "------------------------------------------------"

# FORTRAN at several optimization levels (driver harness fixed at -O0).
for spec in \
    "FORTRAN -O0|-O0 -g" \
    "FORTRAN -O1|-O1 -g" \
    "FORTRAN -O2|-O2 -g" \
    "FORTRAN -O3|-O3 -g" \
    "FORTRAN -O3 -ffast-math|-O3 -ffast-math -g" ; do
    _lbl=${spec%%|*}; _ff=${spec#*|}
    make -C ../fortran clean >/dev/null 2>&1
    make -C ../fortran F77="$F77" FFLAGS="$_ff" >/dev/null 2>&1
    build_and_run "-O0 -g" ../fortran/libminpack.a "$TMP/t"
    gate_row "$_lbl" "$TMP/t"
done

# cminpack at several optimization levels (library and drivers both).
for spec in \
    "cminpack -O0|-O0 -g" \
    "cminpack -O1|-O1 -g" \
    "cminpack -O2|-O2 -g" \
    "cminpack -O3|-O3 -g" \
    "cminpack -O3 -ffast-math|-O3 -ffast-math -g" ; do
    _lbl=${spec%%|*}; _cf=${spec#*|}
    make -C .. clean >/dev/null 2>&1
    make -C .. LIBSUFFIX= CC="$CC" CFLAGS="$_cf" >/dev/null 2>&1
    build_and_run "$_cf" ../libcminpack.a "$TMP/t"
    gate_row "$_lbl" "$TMP/t"
done

echo
echo "(Rows > 0 show problems where FORTRAN -O0 converges but the build does not."
echo " The set varies with -On even for identical FORTRAN source: convergence on"
echo " these extreme problems is a last-bit coin-flip -- see README.md.)"
