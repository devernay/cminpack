#!/bin/sh
# Cross-version comparison of the cminpack solvers, driven by compare.py.
#
# It runs the intensive driver programs (the difficult problems from More,
# Garbow and Hillstrom) built against the different implementations and compares
# their full output -- iteration counts, function- and Jacobian-evaluation
# counts and results -- with the tolerant comparator compare.py:
#
#   * pure C (lmddrvc, hyjdrvc, ...) vs f2c (lmddrv_, hyjdrv_, ...), both linked
#     against libcminpack.a. These MUST be identical (integer counts compared
#     exactly); a difference here is a real bug and makes this script exit 1.
#
#   * f2c vs the original FORTRAN MINPACK (only if a Fortran compiler and
#     fortran/libminpack.a are available). These are expected to differ slightly
#     on hybrd/hybrj because cminpack's enorm uses IEEE-appropriate scaling
#     constants rather than MINPACK's 1980 ones (see the long comment in
#     src/enorm.c). This part is INFORMATIONAL: differences are reported but do
#     not fail the script.
#
# Environment: CC, F77 (default gcc, gfortran), PYTHON (default python3).

set -e
cd "$(dirname "$0")"

CC=${CC:-gcc}
F77=${F77:-gfortran}
PY=${PYTHON:-python3}
CLIB=../libcminpack.a
FLIB=../fortran/libminpack.a

# driver base name -> input data file
data_for() {
    case "$1" in
        lmddrv|lmfdrv|lmsdrv) echo testdata/ssq.data ;;
        hyjdrv|hybdrv)        echo testdata/neq.data ;;
        chkdrv)               echo testdata/chkder.data ;;
    esac
}
DRIVERS="lmddrv lmfdrv lmsdrv hyjdrv hybdrv chkdrv"

rm -rf cross.tmp && mkdir cross.tmp

# make sure the C library exists
[ -f "$CLIB" ] || make -C .. >/dev/null 2>&1 || true
[ -f "$CLIB" ] || { echo "crosscheck: $CLIB not found (run 'make' at the top level first)"; exit 2; }

build_run() { # $1 = library, $2 = suffix tag, builds the driver set and runs it
    lib=$1; tag=$2; progsuffix=$3
    make clean LIBSUFFIX= >/dev/null 2>&1
    progs=""
    for d in $DRIVERS; do progs="$progs ${d}${progsuffix}"; done
    make $progs LIBSUFFIX= MINPACK="$lib" >/dev/null 2>&1
    for d in $DRIVERS; do
        ./${d}${progsuffix} < "$(data_for "$d")" > "cross.tmp/${tag}.${d}.out" 2>&1
    done
}

status=0

echo "=== pure C vs f2c (must be identical) ==="
build_run "$CLIB" C c
build_run "$CLIB" F2C _
for d in $DRIVERS; do
    if $PY compare.py --int-tol 0 "cross.tmp/C.${d}.out" "cross.tmp/F2C.${d}.out" >/dev/null; then
        echo "  $d: identical"
    else
        echo "  $d: DIFFERENT (unexpected!)"
        $PY compare.py --int-tol 0 --max-report 8 "cross.tmp/C.${d}.out" "cross.tmp/F2C.${d}.out" | sed 's/^/    /'
        status=1
    fi
done

if [ -f "$FLIB" ] || $F77 --version >/dev/null 2>&1; then
    [ -f "$FLIB" ] || make -C ../fortran >/dev/null 2>&1 || true
    if [ -f "$FLIB" ]; then
        echo "=== f2c vs FORTRAN MINPACK (informational: enorm differs, see src/enorm.c) ==="
        build_run "$FLIB" FORT _
        for d in $DRIVERS; do
            if $PY compare.py --int-tol 0 "cross.tmp/F2C.${d}.out" "cross.tmp/FORT.${d}.out" >/dev/null; then
                echo "  $d: identical to FORTRAN"
            else
                summary=$($PY compare.py --int-tol 0 "cross.tmp/F2C.${d}.out" "cross.tmp/FORT.${d}.out" | head -1)
                echo "  $d: differs from FORTRAN -> ${summary#DIFFER: * vs *: }"
            fi
        done
    fi
else
    echo "=== f2c vs FORTRAN skipped (no Fortran compiler) ==="
fi

make clean LIBSUFFIX= >/dev/null 2>&1
rm -rf cross.tmp
exit $status
