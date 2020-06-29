#!/bin/sh

# directory of this script
BASEDIR=$(dirname $0)
FSCONFIG="$BASEDIR/../flexiblesusy-config"
SM_EXE="${BASEDIR}/../models/SM/run_SM.x"
SLHA_IN="${BASEDIR}/../models/SM/LesHouches.in.SM"
SLHA_OUT="${BASEDIR}/test_SM_observable_problems.out.spc"
print_block="$BASEDIR/../utils/print_slha_block.awk"

[ $("$FSCONFIG" --with-SM) = yes -a -x ${SM_EXE} ] || {
    echo "Error: SM needs to be build!"
    exit 1;
}

### run for valid point ################################################

"${SM_EXE}" --slha-input-file="${SLHA_IN}" --slha-output-file="${SLHA_OUT}"

[ $? = 0 ] || {
    echo "Error: running valid point failed!"
    exit 1
}

### run for invalid point ##############################################

{ cat "${SLHA_IN}";
  cat <<EOF
Block FlexibleSUSY
    3   0   # calculate SM pole masses
    15  1   # calculate observables
EOF
  } | "${SM_EXE}" --slha-input-file=- --slha-output-file="${SLHA_OUT}"

[ $? = 0 ] || {
    echo "Error: running invalid point failed!"
    exit 1
}

# extract general observable problem flag
general_problem=$(cat "${SLHA_OUT}" | awk -f "$print_block" -v block=OBSINFO | awk '{ if ($1 == 0) print $2 }')

[ -n "$general_problem" ] || {
    echo "Error: expecting a general observable problem flagged!"
    exit 1
}
