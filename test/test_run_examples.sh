#!/bin/sh

# directory of this script
BASEDIR=$(dirname $0)

examples_dir="${BASEDIR}/../examples"

FSCONFIG="$BASEDIR/../flexiblesusy-config"

program_dirs=""

[ $("$FSCONFIG" --with-CMSSM) = yes ] &&
    program_dirs="$program_dirs customized-betas"
[ $("$FSCONFIG" --with-MSSMD5O) = yes -a $("$FSCONFIG" --with-MSSMRHN) = yes ] &&
    program_dirs="$program_dirs tower"

exit_code=0

for dir in ${program_dirs}
do
    echo "> cleaning: rm -rf ${BASEDIR}/${dir}"
    rm -rf ${BASEDIR}/${dir}

    echo "> copying: cp -r ${examples_dir}/${dir}/ ${BASEDIR}"
    cp -r ${examples_dir}/${dir}/ ${BASEDIR}

    echo "> building: make -C ${BASEDIR}/${dir}"
    make -C ${BASEDIR}/${dir}

    exit_code="$?"
    echo "> exit code: ${exit_code}"

    if test ${exit_code} -ne 0; then
        echo "=============="
        echo "Result: FAILED"
        echo "Error: build failed in directory ${dir}"
        echo "=============="
        break
    else
        echo "=========="
        echo "Result: OK"
        echo "=========="
    fi

    echo "> running: (cd ${BASEDIR}/${dir} && ./run_example.sh)"
    (cd ${BASEDIR}/${dir} && ./run_example.sh)

    exit_code="$?"
    echo "> exit code: ${exit_code}"

    if test ${exit_code} -ne 0; then
        echo "=============="
        echo "Result: FAILED"
        echo "Error: execution failed in directory ${dir}"
        echo "=============="
        break
    else
        echo "=========="
        echo "Result: OK"
        echo "=========="
    fi

    echo "> cleaning: rm -rf ${BASEDIR}/${dir}"
    rm -rf ${BASEDIR}/${dir}/
done

exit ${exit_code}
