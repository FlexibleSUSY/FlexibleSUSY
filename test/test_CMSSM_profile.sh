#!/bin/sh

# directory of this script
BASEDIR=$(dirname $0)

VALGRIND="`command -v valgrind`"
CMSSM_EXE="$BASEDIR/../models/CMSSM/run_CMSSM.x"
CMSSM_INPUT="$BASEDIR/../model_files/CMSSM/LesHouches.in.CMSSM"
CMSSM_OUTPUT="$BASEDIR/LesHouches.out.CMSSM"
VALGRIND_OUTPUT="${BASEDIR}/test_CMSSM_profile.valgrind.out"
CALLGRIND_OUTPUT="${BASEDIR}/test_CMSSM_profile.callgrind.out"

if [ -z "$VALGRIND" ]; then
    echo "Error: valgrind not found"
    echo "Profiling will not be performed."
    exit 0
fi

if test ! -x "$CMSSM_EXE"; then
    echo "Error: CMSSM spectrum generator not found: $CMSSM_EXE"
    exit 1
fi

if test ! -f "$CMSSM_INPUT"; then
    echo "Error: CMSSM SLHA input file not found: $CMSSM_INPUT"
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
    ${CMSSM_EXE} \
    --slha-input-file=${CMSSM_INPUT} \
    --slha-output-file=${CMSSM_OUTPUT}

echo "Valgrind log file: ${VALGRIND_OUTPUT}"
echo "Callgring output file: ${CALLGRIND_OUTPUT}"
echo ""
cat ${VALGRIND_OUTPUT}

exit 0
