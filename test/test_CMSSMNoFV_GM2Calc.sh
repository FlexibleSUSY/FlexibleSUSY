#!/bin/sh

# directory of this script
BASEDIR=$(dirname $0)
FSCONFIG="$BASEDIR/../flexiblesusy-config"
CMSSMNoFV_EXE="${BASEDIR}/../models/CMSSMNoFV/run_CMSSMNoFV.x"
SLHA_IN="${BASEDIR}/../models/CMSSMNoFV/LesHouches.in.CMSSMNoFV"
SLHA_OUT="${BASEDIR}/test_CMSSMNoFV_GM2Calc.out.spc"
print_block="$BASEDIR/../utils/print_slha_block.awk"

if [ -e GM2Calc.pc ] ; then
    # shellcheck disable=SC2046
    eval $(grep '^prefix=' GM2Calc.pc)
    # shellcheck disable=SC2154
    GM2CALC_EXE="${prefix}/bin/gm2calc.x"
elif [ -e gm2calc.pc ] ; then
    # shellcheck disable=SC2046
    eval $(grep '^prefix=' gm2calc.pc)
    # shellcheck disable=SC2154
    GM2CALC_EXE="${prefix}/bin/gm2calc.x"
fi

if [ ! -x "${GM2CALC_EXE}" ] ; then
    echo "Cannot find GM2Calc executable in ${GM2CALC_EXE}"
    echo "Skipping test."
    exit
fi

[ "$("$FSCONFIG" --with-CMSSMNoFV)" = yes ] && [ -x "${CMSSMNoFV_EXE}" ] || {
    echo "Error: CMSSMNoFV needs to be build!"
    exit 1;
}

[ "$("$FSCONFIG" --enable-gm2calc)" = yes ] || {
    echo "Error: FlexibleSUSY must be configured with GM2Calc!"
    exit 1;
}

# calculate amu at 1-loop with FS
{ cat "${SLHA_IN}";
  cat <<EOF
Block FlexibleSUSY
    15  1   # calculate observables (a_muon, ...)
    32  1   # 1-loop calculation of amu
EOF
  } | "${CMSSMNoFV_EXE}" --slha-input-file=- --slha-output-file="${SLHA_OUT}"

[ $? = 0 ] || {
    echo "Error: ${CMSSMNoFV_EXE} failed!"
    exit 1
}

alpha_em_MZ_inv=$(cat "${SLHA_OUT}" | awk -f "$print_block" -v block=SMINPUTS | awk '{ if ($1 == 1) print $2 }' | sed -e 's/[eE]/\*10\^/' | sed -e 's/\^+/\^/')

alpha_em_MZ=$(cat <<EOF | bc
scale=17
1.0/(${alpha_em_MZ_inv})
EOF
    )

alpha_em_0=$(cat "${SLHA_OUT}" | awk -f "$print_block" -v block=FlexibleSUSYInput | awk '{ if ($1 == 0) print $2 }')

amu_1l_fs=$(cat "${SLHA_OUT}" | awk -f "$print_block" -v block=FlexibleSUSYLowEnergy | awk '{ if ($1 == 0) print $2 }')
damu_1l_fs=$(cat "${SLHA_OUT}" | awk -f "$print_block" -v block=FlexibleSUSYLowEnergy | awk '{ if ($1 == 1) print $2 }')

# calculate amu at 2-loop with FS
{ cat "${SLHA_IN}";
  cat <<EOF
Block FlexibleSUSY
    15  1   # calculate observables (a_muon, ...)
    32  2   # 2-loop calculation of amu
EOF
  } | "${CMSSMNoFV_EXE}" --slha-input-file=- --slha-output-file="${SLHA_OUT}"

[ $? = 0 ] || {
    echo "Error: ${CMSSMNoFV_EXE} failed!"
    exit 1
}

amu_2l_fs=$(cat "${SLHA_OUT}" | awk -f "$print_block" -v block=FlexibleSUSYLowEnergy | awk '{ if ($1 == 0) print $2 }')
damu_2l_fs=$(cat "${SLHA_OUT}" | awk -f "$print_block" -v block=FlexibleSUSYLowEnergy | awk '{ if ($1 == 1) print $2 }')

amu_2l_gm2calc_fs=$(cat "${SLHA_OUT}" | awk -f "$print_block" -v block=FlexibleSUSYLowEnergy | awk '{ if ($1 == 2) print $2 }')

# calculate amu at 2-loop with GM2Calc
amu_2l_gm2calc=$({ cat <<EOF
Block GM2CalcConfig
     1  2  # loop order
     0  0  # minimal output
     4  0  # verbose output
Block GM2CalcInput
     1  ${alpha_em_MZ}  # alpha(MZ) [1L]
     2  ${alpha_em_0}   # alpha(0)  [2L]
EOF
  cat "${SLHA_OUT}";
      } | "${GM2CALC_EXE}" --slha-input-file=-)

[ $? = 0 ] || {
    echo "Error: ${GM2CALC_EXE} failed!"
    exit 1
}

# calculate amu at 1-loop with GM2Calc
amu_1l_gm2calc=$({ cat <<EOF
Block GM2CalcConfig
     1  1  # loop order
     0  0  # minimal output
     4  0  # verbose output
Block GM2CalcInput
     1  ${alpha_em_MZ}  # alpha(MZ) [1L]
     2  ${alpha_em_0}   # alpha(0)  [2L]
EOF
  cat "${SLHA_OUT}";
      } | "${GM2CALC_EXE}" --slha-input-file=-)

[ $? = 0 ] || {
    echo "Error: ${GM2CALC_EXE} failed!"
    exit 1
}

# convert scientific notation to bc friendly notation
amu_1l_fs=$(echo "${amu_1l_fs}" | sed -e 's/[eE]/\*10\^/' | sed -e 's/\^+/\^/')
damu_1l_fs=$(echo "${damu_1l_fs}" | sed -e 's/[eE]/\*10\^/' | sed -e 's/\^+/\^/')
amu_2l_fs=$(echo "${amu_2l_fs}" | sed -e 's/[eE]/\*10\^/' | sed -e 's/\^+/\^/')
damu_2l_fs=$(echo "${damu_2l_fs}" | sed -e 's/[eE]/\*10\^/' | sed -e 's/\^+/\^/')
amu_2l_gm2calc_fs=$(echo "${amu_2l_gm2calc_fs}" | sed -e 's/[eE]/\*10\^/' | sed -e 's/\^+/\^/')
amu_2l_gm2calc=$(echo "${amu_2l_gm2calc}" | sed -e 's/[eE]/\*10\^/' | sed -e 's/\^+/\^/')
amu_1l_gm2calc=$(echo "${amu_1l_gm2calc}" | sed -e 's/[eE]/\*10\^/' | sed -e 's/\^+/\^/')

errors=0

# compares two values $1 and $2 for relative equality with a maximum relative deviation $3
test_close() {
    local val1
    local val2
    local rel_error
    local diff

    val1="$1"
    val2="$2"
    rel_error="$3"
    diff=$(cat <<EOF | bc $BASEDIR/abs.bc
scale=100
abs((abs(${val1}) - abs(${val2})) / (${val1})) < ${rel_error}
EOF
        )

    if test $diff -ne 1 ; then
        echo "Error: relative difference between"
        echo " ${val1} and ${val2} is larger than ${rel_error}"
        errors=1
    fi
}

# compares two values $1 < $2
test_lt() {
    local val1
    local val2
    local diff

    val1="$1"
    val2="$2"
    diff=$(cat <<EOF | bc
scale=100
${val1} < ${val2}
EOF
        )

    if test $diff -ne 1 ; then
        echo "Error: ${val1} !< ${val2}"
        errors=1
    fi
}

### test 2L GM2Calc vs. embedded 2L GM2Calc

# Note: The agreement between vanilla GM2Calc and the GM2Calc version
# embedded in FlexibleSUSY is not 100% perfect.  The reason is, that
# the embedded GM2Calc version uses the exact smuon pole masses
# (model.get_physical().MSm), while the vanilla GM2Calc uses the smuon
# masses from the SLHA output.  These two pole masses are different,
# because CMSSMNoFV writes the sfermion pole masses to the MASS block
# in SLHA-1 convention, i.e. without the sfermion mixing.
#
# Example point:
#
# exact  MSm = (229.991  360.947)
# SLHA-1 MSm = (230.002  360.937)

test_close "${amu_2l_gm2calc_fs}" "${amu_2l_gm2calc}" "0.0001"

### test 1L FlexibleSUSY vs. 1L GM2Calc

test_close "${amu_1l_gm2calc}" "${amu_1l_fs}" "0.005"

### test 2L FlexibleSUSY vs. embedded 2L GM2Calc

test_close "${amu_2l_gm2calc_fs}" "${amu_2l_fs}" "0.05"

### test uncertainties

test_lt "${damu_2l_fs}" "${damu_1l_fs}"

echo "FlexibleSUSY 1L    : amu = ${amu_1l_fs} +/- ${damu_1l_fs}"
echo "FlexibleSUSY 2L    : amu = ${amu_2l_fs} +/- ${damu_2l_fs}"
echo "original GM2Calc 1L: amu = ${amu_1l_gm2calc}"
echo "original GM2Calc 2L: amu = ${amu_2l_gm2calc}"
echo "embedded GM2Calc 2L: amu = ${amu_2l_gm2calc_fs}"

if test $errors -eq 0 ; then
    echo "Test status: OK"
else
    echo "Test status: FAIL"
fi

exit ${errors}
