#!/bin/sh

# directory of this script
BASEDIR=$(dirname $0)
CMSSM_EXE="$BASEDIR/../models/CMSSM/run_CMSSM.x"

if test ! -x "$CMSSM_EXE"; then
    echo "Error: CMSSM spectrum generator not found: $CMSSM_EXE"
    exit 1
fi

CMSSM_STD_OUTPUT=$(cat <<EOF | ${CMSSM_EXE} --slha-input-file=- --slha-output-file=/dev/null 2> /dev/null
Block MINPAR                 # Input parameters
    1   125                  # m0
    2   500                  # m12
    3   50                   # TanBeta
    4   1                    # SignMu
    5  -500                  # Azero
EOF
   )

exit_code="1"

if [ -n "$CMSSM_STD_OUTPUT" ] ; then
    echo "Error: output was written to cout:"
    echo ""
    echo "$CMSSM_STD_OUTPUT"
    echo ""
    echo "Test result: FAIL"
    exit ${exit_code}
fi

echo "Test result: OK"
