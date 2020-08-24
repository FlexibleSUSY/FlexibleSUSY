#!/bin/bash

# This test compares the real CMSSM with the complex CMSSM (cCMSSM)
# for real input parameters

BASEDIR=$(dirname $0)
UTILSDIR=${BASEDIR}/../utils

cmssm_input=$(
    cat "$BASEDIR/../model_files/CMSSM/LesHouches.in.CMSSM"
    cat <<EOF
Block FlexibleSUSY
    6   2                    # beta-functions loop order
EOF
    )
cmssm_output="$BASEDIR/test_cCMSSM_CMSSM.out.spc"
ccmssm_output="$BASEDIR/test_cCMSSM_cCMSSM.out.spc"

sed_cmd=`command -v sed`
awk_cmd=`command -v awk`
numdiff_cmd=`command -v numdiff`
cmssm_exe="$BASEDIR/../models/CMSSM/run_CMSSM.x"
ccmssm_exe="$BASEDIR/../models/cCMSSM/run_cCMSSM.x"

if [ -z "$sed_cmd" ]; then
    echo "Error: sed command not found"
    exit 1
fi
if [ -z "$awk_cmd" ]; then
    echo "Error: awk command not found"
    exit 1
fi
if [ -z "$numdiff_cmd" ]; then
    echo "Error: numdiff command not found"
    exit 1
fi
if test ! -x "$cmssm_exe"; then
    echo "Error: CMSSM spectrum generator not found: $cmssm_exe"
    exit 1
fi
if test ! -x "$ccmssm_exe"; then
    echo "Error: cCMSSM spectrum generator not found: $ccmssm_exe"
    exit 1
fi

# applies abs to every entry in every line that does not start with `B'
apply_abs() {
    local bc_friendly

    while read line ; do
        case "${line}" in
            B*)
                printf "%s\n" "${line}"
                continue ;;
            *)
                for entry in ${line} ; do
                    # convert scientific notation to to bc friendly input
                    bc_friendly=$(echo "${entry}" | sed -e 's/[eE]/\*10\^/' | sed -e 's/\^+/\^/')
                    # apply abs
                    printf "%s " $(echo  "scale=15; abs(${bc_friendly})" | bc $BASEDIR/abs.bc)
                done
                printf "\n"
                continue ;;
            esac
    done
}

remove_comments() {
    $sed_cmd -e 's/ *#\(.*\)//'
}

print_blocks_to_compare() {
    cat - | tee \
        >($awk_cmd -f $UTILSDIR/print_slha_block.awk -v block=GAUGE              -v omit_comments=1 | remove_comments) \
        >($awk_cmd -f $UTILSDIR/print_slha_block.awk -v block=YU                 -v omit_comments=1 | remove_comments) \
        >($awk_cmd -f $UTILSDIR/print_slha_block.awk -v block=YD                 -v omit_comments=1 | remove_comments) \
        >($awk_cmd -f $UTILSDIR/print_slha_block.awk -v block=YE                 -v omit_comments=1 | remove_comments) \
        >($awk_cmd -f $UTILSDIR/print_slha_block.awk -v block=TU                 -v omit_comments=1 | remove_comments) \
        >($awk_cmd -f $UTILSDIR/print_slha_block.awk -v block=TD                 -v omit_comments=1 | remove_comments) \
        >($awk_cmd -f $UTILSDIR/print_slha_block.awk -v block=TE                 -v omit_comments=1 | remove_comments) \
        >($awk_cmd -f $UTILSDIR/print_slha_block.awk -v block=MSQ2               -v omit_comments=1 | remove_comments) \
        >($awk_cmd -f $UTILSDIR/print_slha_block.awk -v block=MSU2               -v omit_comments=1 | remove_comments) \
        >($awk_cmd -f $UTILSDIR/print_slha_block.awk -v block=MSD2               -v omit_comments=1 | remove_comments) \
        >($awk_cmd -f $UTILSDIR/print_slha_block.awk -v block=MSL2               -v omit_comments=1 | remove_comments) \
        >($awk_cmd -f $UTILSDIR/print_slha_block.awk -v block=MSE2               -v omit_comments=1 | remove_comments) \
        >($awk_cmd -f $UTILSDIR/print_slha_block.awk -v block=Phases             -v omit_comments=1 | remove_comments) \
        >($awk_cmd -f $UTILSDIR/print_slha_block.awk -v block=FlexibleSUSYOutput -v omit_comments=1 | remove_comments) \
        >($awk_cmd -f $UTILSDIR/print_slha_block.awk -v block=ALPHA              -v omit_comments=1 | remove_comments) \
        >($awk_cmd -f $UTILSDIR/print_slha_block.awk -v block=HMIX               -v omit_comments=1 | remove_comments) \
        >($awk_cmd -f $UTILSDIR/print_slha_block.awk -v block=Au                 -v omit_comments=1 | remove_comments) \
        >($awk_cmd -f $UTILSDIR/print_slha_block.awk -v block=Ad                 -v omit_comments=1 | remove_comments) \
        >($awk_cmd -f $UTILSDIR/print_slha_block.awk -v block=Ae                 -v omit_comments=1 | remove_comments) \
        >($awk_cmd -f $UTILSDIR/print_slha_block.awk -v block=MSOFT              -v omit_comments=1 | remove_comments) \
        >($awk_cmd -f $UTILSDIR/print_slha_block.awk -v block=MASS               -v omit_comments=1 | remove_comments | apply_abs) \
        > /dev/null
}

# generate CMSSM point
echo -n "running CMSSM point ... "
echo "$cmssm_input" | $cmssm_exe --slha-input-file=- --slha-output-file=$cmssm_output
echo "done"
echo "CMSSM SLHA output file: $cmssm_output"

# generate cCMSSM point
echo -n "running cCMSSM point ... "
echo "$cmssm_input" | $ccmssm_exe --slha-input-file=- --slha-output-file=$ccmssm_output
echo "done"
echo "cCMSSM SLHA output file: $ccmssm_output"

echo "$(cat "$cmssm_output"  | print_blocks_to_compare)" > "$cmssm_output"
echo "$(cat "$ccmssm_output" | print_blocks_to_compare)" > "$ccmssm_output"

echo ""
echo "comparing files: $cmssm_output and $ccmssm_output ..."

$numdiff_cmd $cmssm_output $ccmssm_output
exit_code=$?

if [ "x$exit_code" != "x0" ]; then
    echo "Error: difference between $cmssm_output and $ccmssm_output"
    echo ""
    echo "Test result: FAIL"
else
    echo ""
    echo "Test result: OK"
fi

exit $exit_code
