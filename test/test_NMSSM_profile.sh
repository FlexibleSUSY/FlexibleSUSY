#!/bin/sh

# directory of this script
BASEDIR=$(dirname $0)

VALGRIND="`command -v valgrind`"
NMSSM_EXE="$BASEDIR/../models/NMSSM/run_NMSSM.x"
NMSSM_INPUT="$BASEDIR/../model_files/NMSSM/LesHouches.in.NMSSM"
NMSSM_OUTPUT="$BASEDIR/LesHouches.out.NMSSM"
VALGRIND_OUTPUT="${BASEDIR}/test_NMSSM_profile.valgrind.out"
CALLGRIND_OUTPUT="${BASEDIR}/test_NMSSM_profile.callgrind.out"

if [ -z "$VALGRIND" ]; then
    echo "Error: valgrind not found"
    echo "Profiling will not be performed."
    exit 0
fi

if test ! -x "$NMSSM_EXE"; then
    echo "Error: NMSSM spectrum generator not found: $NMSSM_EXE"
    exit 1
fi

if test ! -f "$NMSSM_INPUT"; then
    echo "Error: NMSSM SLHA input file not found: $NMSSM_INPUT"
    exit 1
fi


valgrind \
    --tool=callgrind \
    --dump-instr=yes \
    --simulate-cache=yes \
    --collect-jumps=yes \
    --log-file=${VALGRIND_OUTPUT} \
    --callgrind-out-file=${CALLGRIND_OUTPUT} \
    ${NMSSM_EXE} \
    --slha-input-file=${NMSSM_INPUT} \
    --slha-output-file=${NMSSM_OUTPUT}

echo "Valgrind log file: ${VALGRIND_OUTPUT}"
echo "Callgring output file: ${CALLGRIND_OUTPUT}"
echo ""
cat ${VALGRIND_OUTPUT}

exit 0
