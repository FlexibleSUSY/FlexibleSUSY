#!/bin/sh

# directory of this script
BASEDIR=$(dirname $0)

VALGRIND="`command -v valgrind`"
CMSSMNoFV_EXE="$BASEDIR/../models/CMSSMNoFV/run_CMSSMNoFV.x"
CMSSMNoFV_INPUT="$BASEDIR/../model_files/CMSSMNoFV/LesHouches.in.CMSSMNoFV"
CMSSMNoFV_OUTPUT="$BASEDIR/LesHouches.out.CMSSMNoFV"
VALGRIND_OUTPUT="${BASEDIR}/test_CMSSMNoFV_profile.valgrind.out"
CALLGRIND_OUTPUT="${BASEDIR}/test_CMSSMNoFV_profile.callgrind.out"

if [ -z "$VALGRIND" ]; then
    echo "Error: valgrind not found"
    echo "Profiling will not be performed."
    exit 0
fi

if test ! -x "$CMSSMNoFV_EXE"; then
    echo "Error: CMSSMNoFV spectrum generator not found: $CMSSMNoFV_EXE"
    exit 1
fi

if test ! -f "$CMSSMNoFV_INPUT"; then
    echo "Error: CMSSMNoFV SLHA input file not found: $CMSSMNoFV_INPUT"
    exit 1
fi


valgrind \
    --tool=callgrind \
    --separate-threads=yes \
    --dump-instr=yes \
    --simulate-cache=yes \
    --collect-jumps=yes \
    --log-file=${VALGRIND_OUTPUT} \
    --callgrind-out-file=${CALLGRIND_OUTPUT} \
    ${CMSSMNoFV_EXE} \
    --slha-input-file=${CMSSMNoFV_INPUT} \
    --slha-output-file=${CMSSMNoFV_OUTPUT}

echo "Valgrind log file: ${VALGRIND_OUTPUT}"
echo "Callgring output file: ${CALLGRIND_OUTPUT}"
echo ""
cat ${VALGRIND_OUTPUT}

exit 0
