#!/bin/sh

# spectrums are the same.

BASEDIR=$(dirname $0)
print_block="$BASEDIR/../utils/print_slha_block.awk"
print_block_entry="$BASEDIR/../utils/print_slha_block_entry.awk"

cmssmcpv_input="$BASEDIR/test_CMSSMCPV_wrong_higgs_state.in.spc"
cmssmcpv_output="$BASEDIR/test_CMSSMCPV_wrong_higgs_state.out.spc"
cmssm_output="$BASEDIR/test_CMSSM_wrong_higgs_state.out.spc"
rel_error="0.0001"

cmssm_exe="$BASEDIR/../models/CMSSM/run_CMSSM.x"
cmssmcpv_exe="$BASEDIR/../models/CMSSMCPV/run_CMSSMCPV.x"

if test ! -x "$cmssmcpv_exe"; then
    echo "Error: CMSSMCPV spectrum generator not found: $cmssmcpv_exe"
    exit 1
fi

if test ! -x "$cmssm_exe"; then
    echo "Error: CMSSM spectrum generator not found: $cmssm_exe"
    exit 1
fi

echo -n "running CMSSMCPV point ... "
$cmssmcpv_exe --slha-input-file=$cmssmcpv_input --slha-output-file=$cmssmcpv_output > /dev/null 2>&1
echo "done"
echo "CMSSMCPV SLHA input file:  $cmssmcpv_input"
echo "CMSSMCPV SLHA output file: $cmssmcpv_output"

echo -n "running CMSSMCPV point ... "
$cmssm_exe --slha-input-file=$cmssmcpv_input --slha-output-file=$cmssm_output > /dev/null 2>&1
echo "done"
echo "CMSSM SLHA input file:  $cmssmcpv_input"
echo "CMSSM SLHA output file: $cmssm_output"

mh_cmssmcpv=$(awk -f "$print_block" -v block="MASS" $cmssmcpv_output | \
              awk -f "$print_block_entry" -v entries="25")
mh_cmssm=$(awk -f "$print_block" -v block="MASS" $cmssm_output | \
           awk -f "$print_block_entry" -v entries="25")

# convert scientific notation to bc friendly notation
mh_cmssmcpv="`echo ${mh_cmssmcpv} | sed -e 's/[eE]+*/\\*10\\^/'`"
mh_cmssm="`echo ${mh_cmssm} | sed -e 's/[eE]+*/\\*10\\^/'`"

diff=$(cat <<EOF | bc $BASEDIR/abs.bc
scale=10
abs($mh_cmssm - $mh_cmssmcpv) / ($mh_cmssm) < $rel_error
EOF
    )

if test $diff -ne 1 ; then
    echo "Error: relative difference between"
    echo " $mh_cmssmcpv and $mh_cmssm is larger than $rel_error"
    echo "Test status: FAIL"
    exit 1
else
    echo "Test status: OK"
fi

exit 0
