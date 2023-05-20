#!/bin/sh

# This script generates a CNMSSM spectrum using the semi-analytic solver
# and writes it to a SLHA file.  The EWSB output parameters from this
# point are then used as inputs for the semi-constrained NMSSM solved
# using the two-scale solver.  After running the two-scale NMSSM, we
# compare that the semi-analytic and two-scale spectrums are the same.

BASEDIR=$(dirname $0)
UTILSDIR=${BASEDIR}/../utils
print_block="$UTILSDIR/print_slha_block.awk"
print_block_entry="$UTILSDIR/print_slha_block_entry.awk"
remove_block="$UTILSDIR/remove_slha_block.awk"

semi_analytic_slha="$BASEDIR/../model_files/CNMSSM/LesHouches.in.CNMSSM"
semi_analytic_input="$BASEDIR/test_CNMSSM_spectrum.in.spc"
semi_analytic_output="$BASEDIR/test_CNMSSM_spectrum.out.spc"
two_scale_input="$BASEDIR/test_CNMSSM_spectrum_NMSSM.in.spc"
two_scale_output="$BASEDIR/test_CNMSSM_spectrum_NMSSM.out.spc"
rel_error="1.7e-2"

sed_cmd=`command -v sed`
awk_cmd=`command -v awk`
numdiff_cmd=`command -v numdiff`
semi_analytic_exe="$BASEDIR/../models/CNMSSM/run_CNMSSM.x"
two_scale_exe="$BASEDIR/../models/NMSSM/run_NMSSM.x"

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

if test ! -x "$semi_analytic_exe"; then
    echo "Error: CNMSSM spectrum generator not found: $semi_analytic_exe"
    exit 1
fi

if test ! -x "$two_scale_exe"; then
    echo "Error: NMSSM spectrum generator not found: $two_scale_exe"
    exit 1
fi

remove_mixing_and_input_blocks() {
    $awk_cmd -f "$remove_block" -v block=UMIX \
        | $awk_cmd -f "$remove_block" -v block=VMIX \
        | $awk_cmd -f "$remove_block" -v block=NMAMIX \
        | $awk_cmd -f "$remove_block" -v block=DSQMIX \
        | $awk_cmd -f "$remove_block" -v block=SELMIX \
        | $awk_cmd -f "$remove_block" -v block=NMHMIX \
        | $awk_cmd -f "$remove_block" -v block=NMNMIX \
        | $awk_cmd -f "$remove_block" -v block=CHARGEMIX \
        | $awk_cmd -f "$remove_block" -v block=USQMIX \
        | $awk_cmd -f "$remove_block" -v block=SNUMIX \
        | $awk_cmd -f "$remove_block" -v block=SMINPUTS \
        | $awk_cmd -f "$remove_block" -v block=MINPAR \
        | $awk_cmd -f "$remove_block" -v block=EXTPAR \
        | $awk_cmd -f "$remove_block" -v block=SPINFO \
        | $awk_cmd -f "$remove_block" -v block=MODSEL \
        | $awk_cmd -f "$remove_block" -v block=FlexibleSUSY
}

remove_extra_blocks() {
    $awk_cmd -f "$remove_block" -v block=EWSBOutputs
}

match_input_blocks() {
    if test $# -lt 3; then
        echo "match_input_blocks: Too few arguments"
        return 1
    fi
    _input_file="$1"
    _output_file="$2"
    _m0_val="$3"

    $awk_cmd -f "$print_block" -v block=FlexibleSUSY "$_input_file" \
        | $awk_cmd -f "$remove_block" -v block=FlexibleSUSY -v entry=2 \
             > "$_output_file"
    echo "    2   0" >> "$_output_file"

    $awk_cmd -f "$remove_block" -v block=FlexibleSUSY "$_input_file" \
        | $awk_cmd -f "$remove_block" -v block=EXTPAR \
                   >> "$_output_file"

    echo "   1   $_m0_val" >> "$_output_file"

    $awk_cmd -f "$print_block" -v block=EXTPAR "$_input_file" \
             >> "$_output_file"
}

# generate point using semi-analytic solver
cp $semi_analytic_slha $semi_analytic_input
echo -n "running semi-analytic solver ..."
$semi_analytic_exe --slha-input-file=$semi_analytic_input \
                   --slha-output-file=$semi_analytic_output
echo "done"
echo "CNMSSM SLHA input file:  $semi_analytic_input"
echo "CNMSSM SLHA output file: $semi_analytic_output"

if test ! -r "$semi_analytic_output"; then
    echo "Error: generated CNMSSM spectrum not found: $semi_analytic_output"
    exit 1
fi

# remove comments from the semi-analytic output spectrum
$sed_cmd -i~ -e '/^ *#/d' $semi_analytic_output

# construct matching input point for the two-scale solver
m0Sq_sol=$(awk -f "$print_block" -v block=EWSBOutputs "$semi_analytic_output" \
               | awk -f "$print_block_entry" -v entries=3 \
               | sed -e 's/[eE]+*/*10^/')
m0_sol=$(echo "scale=20; sqrt($m0Sq_sol)" | bc)
match_input_blocks $semi_analytic_input $two_scale_input $m0_sol

# generate point using the two-scale solver
echo -n "running two-scale solver ..."
$two_scale_exe --slha-input-file=$two_scale_input --slha-output-file=$two_scale_output
echo "done"
echo "NMSSM SLHA input file:  $two_scale_input"
echo "NMSSM SLHA output file: $two_scale_output"

# remove comments from two-scale output spectrum
$sed_cmd -i~ -e '/^ *#/d' $two_scale_output

# remove mixing matrix blocks because we don't want to compare objects
# with phase ambiguities, and input blocks since these do not match
cp $semi_analytic_output $semi_analytic_output~
cp $two_scale_output $two_scale_output~

remove_mixing_and_input_blocks < $semi_analytic_output~ > $semi_analytic_output
remove_mixing_and_input_blocks < $two_scale_output~ > $two_scale_output

# also remove additional output blocks present for the semi-analytic
# solver but not for the two-scale solver
cp $semi_analytic_output $semi_analytic_output~

remove_extra_blocks < $semi_analytic_output~ > $semi_analytic_output

echo ""
echo "comparing files $semi_analytic_output and $two_scale_output ..."
echo "required relative agreement: $rel_error"

diff=`$numdiff_cmd\
 --absolute-tolerance=1.0e-12\
 --relative-tolerance=$rel_error\
 $semi_analytic_output $two_scale_output`

diff_without_comments=`echo $diff | $sed_cmd -e '/^ *#/d' | $sed_cmd -e '/^+++/d'`

exit_code=0

if [ -n "$diff_without_comments" ]; then
    echo "Error: difference between $semi_analytic_output and  $two_scale_output larger than $rel_error"
    echo "$diff"
    echo ""
    echo "Test result: FAIL"
    exit_code=1
else
    echo "$diff"
    echo ""
    echo "Test result: OK"
fi

exit $exit_code
