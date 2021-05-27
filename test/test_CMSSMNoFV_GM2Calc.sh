#!/bin/sh

# directory of this script
BASEDIR=$(dirname $0)
FSCONFIG="$BASEDIR/../flexiblesusy-config"
CMSSMNoFV_EXE="${BASEDIR}/../models/CMSSMNoFV/run_CMSSMNoFV.x"
SLHA_IN="${BASEDIR}/../models/CMSSMNoFV/LesHouches.in.CMSSMNoFV"
SLHA_OUT="${BASEDIR}/test_CMSSMNoFV_GM2Calc.out.spc"
print_block="$BASEDIR/../utils/print_slha_block.awk"

if [ -e GM2Calc.pc ] ; then
    eval $(grep '^prefix=' GM2Calc.pc)
    GM2CALC_EXE="${prefix}/bin/gm2calc.x"
fi

if [ ! -x "${GM2CALC_EXE}" ] ; then
    echo "Cannot find GM2Calc executable in ${GM2CALC_EXE}"
    echo "Skipping test."
    exit
fi

[ $("$FSCONFIG" --with-CMSSMNoFV) = yes -a -x ${CMSSMNoFV_EXE} ] || {
    echo "Error: CMSSMNoFV needs to be build!"
    exit 1;
}

[ $("$FSCONFIG" --enable-gm2calc) = yes ] || {
    echo "Error: FlexibleSUSY must be configured with GM2Calc!"
    exit 1;
}

{ cat "${SLHA_IN}";
  cat <<EOF
Block FlexibleSUSY
    15  1   # calculate observables (a_muon, ...)
EOF
  } | "${CMSSMNoFV_EXE}" --slha-input-file=- --slha-output-file="${SLHA_OUT}"

[ $? = 0 ] || {
    echo "Error: ${CMSSMNoFV_EXE} failed!"
    exit 1
}

amu_fs=$(cat "${SLHA_OUT}" | awk -f "$print_block" -v block=FlexibleSUSYLowEnergy | awk '{ if ($1 == 0) print $2 }')

amu_gm2calc_fs=$(cat "${SLHA_OUT}" | awk -f "$print_block" -v block=FlexibleSUSYLowEnergy | awk '{ if ($1 == 2) print $2 }')

amu_gm2calc=$({ cat <<EOF
Block GM2CalcConfig
     0  0  # minimal output
     4  0  # verbose output
EOF
  cat "${SLHA_OUT}";
      } | "${GM2CALC_EXE}" --slha-input-file=-)

[ $? = 0 ] || {
    echo "Error: ${GM2CALC_EXE} failed!"
    exit 1
}

# convert scientific notation to bc friendly notation
amu_fs=$(echo "${amu_fs}" | sed -e 's/[eE]/\*10\^/' | sed -e 's/\^+/\^/')
amu_gm2calc_fs=$(echo "${amu_gm2calc_fs}" | sed -e 's/[eE]/\*10\^/' | sed -e 's/\^+/\^/')
amu_gm2calc=$(echo "${amu_gm2calc}" | sed -e 's/[eE]/\*10\^/' | sed -e 's/\^+/\^/')

### test GM2Calc vs. embedded GM2Calc

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

rel_error=0.0001

diff=$(cat <<EOF | bc $BASEDIR/abs.bc
scale=100
abs((abs($amu_gm2calc_fs) - abs($amu_gm2calc)) / ($amu_gm2calc_fs)) < $rel_error
EOF
    )

errors=0

if test $diff -ne 1 ; then
    echo "Error: relative difference between"
    echo " $amu_gm2calc_fs and $amu_gm2calc is larger than $rel_error"
    errors=1
fi

### test FlexibleSUSY vs. embedded GM2Calc

rel_error=0.04

diff=$(cat <<EOF | bc $BASEDIR/abs.bc
scale=100
abs((abs($amu_gm2calc_fs) - abs($amu_fs)) / ($amu_fs)) < $rel_error
EOF
    )

if test $diff -ne 1 ; then
    echo "Error: relative difference between"
    echo " $amu_gm2calc_fs and $amu_fs is larger than $rel_error"
    errors=1
fi

echo "FlexibleSUSY 1L + 2L QED                       : amu = $amu_fs"
echo "embedded GM2Calc                               : amu = $amu_gm2calc_fs"
echo "original GM2Calc                               : amu = $amu_gm2calc"

if test $errors -eq 0 ; then
    echo "Test status: OK"
else
    echo "Test status: FAIL"
fi

exit ${errors}
