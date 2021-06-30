#!/bin/sh

# directory of this script
BASEDIR=$(dirname $0 | sed 's|^./||')
HOMEDIR="${BASEDIR}/.."
FSCONFIG="${HOMEDIR}/flexiblesusy-config"
CXX=$($FSCONFIG --cxx)
DEPGEN="${HOMEDIR}/config/depgen.x"
OUTPUT="${BASEDIR}/test_depgen.out"

error=0

# removes line breaks (indicated with trailing backslash)
awk_join_lines='/\\$/ { printf "%s", substr($0, 1, length($0)-1); next } { print }'

run_depgens() {
    local cmd1="$1"
    local cmd2="$2"
    local file="$3"
    local diff=

    echo "--------------------------------------"
    echo "running: $cmd1 $file"
    echo "running: $cmd2 $file"

    rm -f ${OUTPUT}.1 ${OUTPUT}.2

    $cmd1 $file | awk "${awk_join_lines}" | sed '/^\/.*/d' | sed '$!N; /^\(.*\)\n\1$/!P; D' | sed 's/ */ /g' > ${OUTPUT}.1
    $cmd2 $file | awk "${awk_join_lines}" | sed '/^\/.*/d' | sed '$!N; /^\(.*\)\n\1$/!P; D' | sed 's/ */ /g' > ${OUTPUT}.2

    diff=$(diff -u ${OUTPUT}.1 ${OUTPUT}.2)

    if test -n "$diff"; then
        echo "Error: files are not equal ${OUTPUT}.1 and ${OUTPUT}.2"
        echo "Difference:"
        echo "$diff"
        echo "Test result: FAIL"
        error=1
    else
        echo "Test result: OK"
    fi
    echo "--------------------------------------"
    echo ""
}

# test inclusion of "base.hpp" header that exists in different directories

flags="-MM"
run_depgens "$CXX $flags" "$DEPGEN $flags" "${BASEDIR}/depgen/base.cpp"
run_depgens "$CXX $flags" "$DEPGEN $flags" "${BASEDIR}/depgen/subdir/base.cpp"

flags="-I${BASEDIR}/depgen -MM"
run_depgens "$CXX $flags" "$DEPGEN $flags" "${BASEDIR}/depgen/base.cpp"
run_depgens "$CXX $flags" "$DEPGEN $flags" "${BASEDIR}/depgen/subdir/base.cpp"

flags="-I${BASEDIR}/depgen -I${BASEDIR}/depgen/subdir -MM"
run_depgens "$CXX $flags" "$DEPGEN $flags" "${BASEDIR}/depgen/base.cpp"
run_depgens "$CXX $flags" "$DEPGEN $flags" "${BASEDIR}/depgen/subdir/base.cpp"

# test inclusion of non-existing header files

flags="-MM -MG"
run_depgens "$CXX $flags" "$DEPGEN $flags" "${BASEDIR}/depgen/nonexisting.cpp"
run_depgens "$CXX $flags" "$DEPGEN $flags" "${BASEDIR}/depgen/subdir/nonexisting.cpp"

flags="-I${BASEDIR}/depgen -MM -MG"
run_depgens "$CXX $flags" "$DEPGEN $flags" "${BASEDIR}/depgen/nonexisting.cpp"
run_depgens "$CXX $flags" "$DEPGEN $flags" "${BASEDIR}/depgen/subdir/nonexisting.cpp"

flags="-I${BASEDIR}/depgen -I${BASEDIR}/depgen/subdir -MM -MG"
run_depgens "$CXX $flags" "$DEPGEN $flags" "${BASEDIR}/depgen/nonexisting.cpp"
run_depgens "$CXX $flags" "$DEPGEN $flags" "${BASEDIR}/depgen/subdir/nonexisting.cpp"

# test inclusion of header in base directory

flags="-MM -MG"
run_depgens "$CXX $flags" "$DEPGEN $flags" "${BASEDIR}/depgen/subdir/include_basedir_header.cpp"

# test creation of phony targets

flags="-MM -MP"
run_depgens "$CXX $flags" "$DEPGEN $flags" "${BASEDIR}/depgen/base.cpp"
run_depgens "$CXX $flags" "$DEPGEN $flags" "${BASEDIR}/depgen/subdir/base.cpp"

# test setting of target name

flags="-MM -MT 'X.o'"
run_depgens "$CXX $flags" "$DEPGEN $flags" "${BASEDIR}/depgen/base.cpp"
run_depgens "$CXX $flags" "$DEPGEN $flags" "${BASEDIR}/depgen/subdir/base.cpp"

# test comment after include statement

flags="-MM"
run_depgens "$CXX $flags" "$DEPGEN $flags" "${BASEDIR}/depgen/comment.cpp"

# test circular dependence
run_depgens "$CXX $flags" "$DEPGEN $flags" "${BASEDIR}/depgen/circular.cpp"

# test #ifdef (currently fails)
flags="-MM"
# run_depgens "$CXX $flags" "$DEPGEN $flags" "${BASEDIR}/depgen/ifdef.cpp"

# test -MF without file argument
flags="-MF"
run_depgens "$CXX $flags" "$DEPGEN $flags" "${BASEDIR}/depgen/base.cpp"
flags="-MF -MF"
run_depgens "$CXX $flags" "$DEPGEN $flags" "${BASEDIR}/depgen/base.cpp"

rm -f ${OUTPUT}*

echo ""
echo "======================================"
if test ${error} -eq 0 ; then
    echo "Test result: OK"
else
    echo "Test result: FAIL"
fi
echo "======================================"

exit ${error}
