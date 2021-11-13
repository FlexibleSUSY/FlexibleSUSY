#!/bin/sh
BASEDIR=$(dirname $0)
HOMEDIR="${BASEDIR}/.."

models="
THDMIIMSSMBC,THDMIIMSSMBCApprox
HGTHDMIIMSSMBC,HGTHDMIIMSSMBCApprox
"

scanout1="${BASEDIR}/test_THDMIIMSSMBCFull_approximation_1.dat"
scanout2="${BASEDIR}/test_THDMIIMSSMBCFull_approximation_2.dat"

errors=0

for m in ${models}; do
    model1=$(echo $m | awk -F , '{ print $1 }')
    model2=$(echo $m | awk -F , '{ print $2 }')

    exe1="${HOMEDIR}/models/${model1}/run_${model1}.x"
    exe2="${HOMEDIR}/models/${model2}/run_${model2}.x"

    slha1="${HOMEDIR}/models/${model1}/LesHouches.in.${model1}"
    slha2="${HOMEDIR}/models/${model2}/LesHouches.in.${model2}"

    if test ! -x "${exe1}"; then
        echo "Error: ${model1} spectrum generator not executable: ${exe1}"
        errors=$(expr $errors + 1)
        continue
    fi
    if test ! -x "${exe2}"; then
        echo "Error: ${model2} spectrum generator not executable: ${exe2}"
        errors=$(expr $errors + 1)
        continue
    fi

    printf "Scanning ${model1} ... "
    ${HOMEDIR}/utils/scan-slha.sh \
              --spectrum-generator="${exe1}" \
              --slha-input-file="${slha1}" \
              --scan-range=MINPAR[3]=1~30:20 \
              --output=MINPAR[3],MASS[25] \
              > "${scanout1}"
    printf "done\n"

    printf "Scanning ${model2} ... "
    ${HOMEDIR}/utils/scan-slha.sh \
              --spectrum-generator="${exe2}" \
              --slha-input-file="${slha2}" \
              --scan-range=MINPAR[3]=1~30:20 \
              --output=MINPAR[3],MASS[25] \
              > "${scanout2}"
    printf "done\n"

    diff=$(numdiff --relative-tolerance=1.0e-8 --brief "${scanout1}" "${scanout2}" | sed '/^+++/d;/^$/d')

    if test "x${diff}" != "x" ; then
        echo "Error: scan output files are not equal:"
        echo "${diff}"
        errors=$(expr $errors + 1)
    fi
done

if [ $errors -ne 0 ]; then
    echo ""
    echo "Error: there were differences in the scan outputs!"
    echo ""
    echo "Test result: FAIL"
else
    echo ""
    echo "Test result: OK"
fi

rm -rf "${scanout1}" "${scanout2}"

exit $errors
