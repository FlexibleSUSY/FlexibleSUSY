#!/bin/sh

# directory of this script
BASEDIR=$(dirname $0)

VALGRIND="`command -v valgrind`"
MSSMEFTHiggs_EXE="$BASEDIR/../models/MSSMEFTHiggs/run_MSSMEFTHiggs.x"
MSSMEFTHiggs_INPUT="$BASEDIR/../model_files/MSSMEFTHiggs/LesHouches.in.MSSMEFTHiggs"
MSSMEFTHiggs_OUTPUT="$BASEDIR/LesHouches.out.MSSMEFTHiggs"
VALGRIND_OUTPUT="${BASEDIR}/test_MSSMEFTHiggs_profile.valgrind.out"
CALLGRIND_OUTPUT="${BASEDIR}/test_MSSMEFTHiggs_profile.callgrind.out"

if [ -z "$VALGRIND" ]; then
    echo "Error: valgrind not found"
    echo "Profiling will not be performed."
    exit 0
fi

if test ! -x "$MSSMEFTHiggs_EXE"; then
    echo "Error: MSSMEFTHiggs spectrum generator not found: $MSSMEFTHiggs_EXE"
    exit 1
fi

if test ! -f "$MSSMEFTHiggs_INPUT"; then
    echo "Error: MSSMEFTHiggs SLHA input file not found: $MSSMEFTHiggs_INPUT"
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
    ${MSSMEFTHiggs_EXE} \
    --slha-input-file=${MSSMEFTHiggs_INPUT} \
    --slha-output-file=${MSSMEFTHiggs_OUTPUT}

echo "Valgrind log file: ${VALGRIND_OUTPUT}"
echo "Callgring output file: ${CALLGRIND_OUTPUT}"
echo ""
cat ${VALGRIND_OUTPUT}

exit 0
