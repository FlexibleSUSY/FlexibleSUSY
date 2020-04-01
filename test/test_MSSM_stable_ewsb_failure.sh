#!/bin/sh

BASEDIR=$(dirname $0)
FSCONFIG="$BASEDIR/../flexiblesusy-config"
MSSM_EXE="${BASEDIR}/../models/MSSM/run_MSSM.x"
print_block="$BASEDIR/../utils/print_slha_block.awk"

[ $("$FSCONFIG" --with-MSSM) = yes -a -x ${MSSM_EXE} ] || {
    echo "Error: MSSM needs to be build!"
    exit 1;
}

slha_in="${BASEDIR}/test_MSSM_stable_ewsb_failure.in.spc"
slha_out=$("${MSSM_EXE}" --slha-input-file="$slha_in" 2>/dev/null)

spinfo_4=$(echo "$slha_out" | \
                  awk -f "$print_block" -v block=SPINFO | \
                  awk '{ if ($1 == 4) { $1 = ""; print $0 } }')

error=0

case "$spinfo_4" in
    *no\ convergence*) error=1 ;;
esac

if [ $error -eq 0 ] ; then
    echo "Test result: OK"
else
    echo "Error: point has not converged, because"
    echo "   SPINFO[4] contains \"no convergence\""
    echo ""
    echo "Test result: FAIL"
fi

exit ${error}
