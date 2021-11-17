#!/bin/sh

# directory of this script
BASEDIR=$(dirname $0)
FSCONFIG="${BASEDIR}/../flexiblesusy-config"
MODELS="SM MSSMNoFVHimalaya"
exit_code=0

exit_and_log() {
    echo "=============="
    echo "Result: FAILED"
    echo "Error: build failed in directory ${dest}"
    echo "=============="
    exit 1
}

run() {
    local cmd="$1"
    echo "$cmd"
    eval "$cmd" || exit_and_log
}

[ ! -x "$FSCONFIG" ] && { echo "Error: $FSCONFIG is not executable"; exit_and_log; }

for m in ${MODELS} ; do
    source="${BASEDIR}/../models/${m}"
    dest="${BASEDIR}/${m}"

    [ "$("$FSCONFIG" --with-${m})" = no ] && { echo "Warning: ${m} is not configured" ; continue; }

    echo "> copying ${source} -> ${dest}"

    rm -rf "${dest}"
    cp -r "${source}" "${dest}"
    rm -f "${dest}"/*.a "${dest}"/*.d "${dest}"/*.o "${dest}"/*.so "${dest}"/*.x

    echo "> compiling ${dest}"
    for f in "${dest}"/*.cpp ; do
        fo="$(echo "$f" | sed 's,\.cpp,.o,')"
        run "$(${FSCONFIG} --compile-cmd) $f -c -o $fo"
    done

    echo "> creating library"
    run "$(${FSCONFIG} --module-make-lib-cmd) ${dest}/lib${m}.a ${dest}/${m}_*.o"

    echo "> creating executable"
    run "$(${FSCONFIG} --cxx) -o ${dest}/run_${m}.x ${dest}/run_${m}.o $(${FSCONFIG} --libs)"

    echo "> running spectrum generator"
    "${dest}"/run_${m}.x --slha-input-file="${dest}/LesHouches.in.${m}"

    exit_code="$?"

    if test ${exit_code} -ne 0 ; then
        echo "=============="
        echo "Result: FAILED"
        echo "Error: build failed in directory ${dest}"
        echo "=============="
    else
        echo "=========="
        echo "Result: OK"
        echo "=========="
    fi
done

exit ${exit_code}
