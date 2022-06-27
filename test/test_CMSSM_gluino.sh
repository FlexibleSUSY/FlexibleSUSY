#!/bin/sh

BASEDIR=$(dirname $0)
print_block="$BASEDIR/../utils/print_slha_block.awk"
print_block_entry="$BASEDIR/../utils/print_slha_block_entry.awk"

mssm_exe="$BASEDIR/../models/CMSSM/run_CMSSM.x"
softsusy_exe="$BASEDIR/../test/SOFTSUSY/run_softpoint.x"

input="$BASEDIR/test_CMSSM_gluino.spc.in"
output="$BASEDIR/test_CMSSM_gluino.spc.out"
softsusy_output="$BASEDIR/test_CMSSM_gluino.softsusy.spc.out"

rel_error="0.00001"

if test ! -x "$mssm_exe"; then
    echo "Error: CMSSM spectrum generator not found: $mssm_exe"
    exit 1
fi

if test ! -x "$softsusy_exe"; then
    echo "Error: SoftsusyNMSSM spectrum generator not found: $softsusy_exe"
    exit 1
fi

echo -n "running CMSSM point ... "
$mssm_exe --slha-input-file=$input --slha-output-file=$output > /dev/null 2>&1
[ $? -ne 0 ] && { echo "\nError: $mssm_exe quitted unexpectedly"; exit 1; }
echo "done"
echo "CMSSM SLHA input file:  $input"
echo "CMSSM SLHA output file: $output"

echo -n "running CMSSM point ... "
$softsusy_exe leshouches < $input > $softsusy_output 2> /dev/null
[ $? -ne 0 ] && { echo "\nError: $softsusy_exe quitted unexpectedly"; exit 1; }
echo "done"
echo "SoftsusyNMSSM SLHA input file:  $input"
echo "SoftsusyNMSSM SLHA output file: $softsusy_output"

mg_mssm=$(awk -f "$print_block" -v block="MASS" $output | \
          awk -f "$print_block_entry" -v entries="1000021")
mg_softsusy=$(awk -f "$print_block" -v block="MASS" $output | \
              awk -f "$print_block_entry" -v entries="1000021")

# convert scientific notation to bc friendly notation
mg_mssm="`echo ${mg_mssm} | sed -e 's/[eE]+*/\\*10\\^/'`"
mg_softsusy="`echo ${mg_softsusy} | sed -e 's/[eE]+*/\\*10\\^/'`"

diff=$(cat <<EOF | bc $BASEDIR/abs.bc
scale=10
abs((abs($mg_softsusy) - abs($mg_mssm)) / ($mg_softsusy)) < $rel_error
EOF
    )

if test $diff -ne 1 ; then
    echo "Error: relative difference between"
    echo " $mg_mssm and $mg_softsusy is larger than $rel_error"
    echo "Test status: FAIL"
    exit 1
else
    echo "FlexibleSUSY gluino pole mass: $mg_mssm"
    echo "SoftSUSY     gluino pole mass: $mg_softsusy"
    echo "Test status: OK"
fi

exit 0
