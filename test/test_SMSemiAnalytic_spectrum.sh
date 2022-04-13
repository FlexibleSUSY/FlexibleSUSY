#!/bin/sh

# This script generates an SM spectrum using both the semi-analytic
# and two-scale solvers, and compares that the resulting spectrums are
# the same.

BASEDIR=$(dirname $0)
UTILSDIR=${BASEDIR}/../utils
remove_block="$UTILSDIR/remove_slha_block.awk"

semi_analytic_model="SMSemiAnalytic"
two_scale_model="SM"

semi_analytic_slha="$BASEDIR/../model_files/${semi_analytic_model}/LesHouches.in.${semi_analytic_model}"
semi_analytic_input="$BASEDIR/test_${semi_analytic_model}_spectrum.in.spc"
semi_analytic_output="$BASEDIR/test_${semi_analytic_model}_spectrum.out.spc"
two_scale_slha="$BASEDIR/../model_files/${two_scale_model}/LesHouches.in.${two_scale_model}"
two_scale_input="$BASEDIR/test_${semi_analytic_model}_spectrum_${two_scale_model}.in.spc"
two_scale_output="$BASEDIR/test_${semi_analytic_model}_spectrum_${two_scale_model}.out.spc"
rel_error="1.7e-2"

sed_cmd=`command -v sed`
awk_cmd=`command -v awk`
numdiff_cmd=`command -v numdiff`
semi_analytic_exe="$BASEDIR/../models/${semi_analytic_model}/run_${semi_analytic_model}.x"
two_scale_exe="$BASEDIR/../models/${two_scale_model}/run_${two_scale_model}.x"

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
    echo "Error: $semi_analytic_model spectrum generator not found: $semi_analytic_exe"
    exit 1
fi

if test ! -x "$two_scale_exe"; then
    echo "Error: $two_scale_model spectrum generator not found: $two_scale_exe"
    exit 1
fi

remove_mixing_and_input_blocks() {
    $awk_cmd -f "$remove_block" -v block=UULMIX \
        | $awk_cmd -f "$remove_block" -v block=UDLMIX \
        | $awk_cmd -f "$remove_block" -v block=UURMIX \
        | $awk_cmd -f "$remove_block" -v block=UDRMIX \
        | $awk_cmd -f "$remove_block" -v block=UELMIX \
        | $awk_cmd -f "$remove_block" -v block=UERMIX \
        | $awk_cmd -f "$remove_block" -v block=SMINPUTS \
        | $awk_cmd -f "$remove_block" -v block=MINPAR \
        | $awk_cmd -f "$remove_block" -v block=EXTPAR \
        | $awk_cmd -f "$remove_block" -v block=SPINFO \
        | $awk_cmd -f "$remove_block" -v block=MODSEL \
        | $awk_cmd -f "$remove_block" -v block=FlexibleSUSY
}

remove_extra_blocks() {
    $awk_cmd -f "$remove_block" -v block=EWSBOutputs \
        | $awk_cmd -f "$remove_block" -v block=FlexibleSUSYOutput
}

# generate point using each solver
cp "$semi_analytic_slha" "$semi_analytic_input"
echo -n "running semi-analytic solver ..."
$semi_analytic_exe --slha-input-file="$semi_analytic_input" \
                   --slha-output-file="$semi_analytic_output"
echo "done"
echo "$semi_analytic_model SLHA input file: $semi_analytic_input"
echo "$semi_analytic_model SLHA output file: $semi_analytic_output"

if test ! -r "${semi_analytic_output}"; then
    echo "Error: generated $semi_analytic_model spectrum not found: $semi_analytic_output"
    exit 1
fi

cp "$two_scale_slha" "$two_scale_input"
echo -n "running two-scale solver ..."
$two_scale_exe --slha-input-file="$two_scale_input" \
               --slha-output-file="$two_scale_output"
echo "done"
echo "$two_scale_model SLHA input file: $two_scale_input"
echo "$two_scale_model SLHA output file: $two_scale_output"

if test ! -r "${two_scale_output}"; then
    echo "Error: generated $two_scale_model spectrum not found: $two_scale_output"
    exit 1
fi

# remove comments from the output files
$sed_cmd -i~ -e '/^ *#/d' "$semi_analytic_output"
$sed_cmd -i~ -e '/^ *#/d' "$two_scale_output"

# remove mixing matrix blocks because we don't want to compare objects
# with phase ambiguities, and input blocks
cp "$semi_analytic_output" "$semi_analytic_output~"
cp "$two_scale_output" "$two_scale_output~"

remove_mixing_and_input_blocks < "$semi_analytic_output~" > "$semi_analytic_output"
remove_mixing_and_input_blocks < "$two_scale_output~" > "$two_scale_output"

# also remove additional output blocks present for the semi-analytic
# solver but not for the two-scale solver
cp "$semi_analytic_output" "$semi_analytic_output~"
remove_extra_blocks < "$semi_analytic_output~" > "$semi_analytic_output"

echo ""
echo "comparing files $semi_analytic_output and $two_scale_output ..."
echo "required relative agreement: $rel_error"

diff=`$numdiff_cmd\
 --absolute-tolerance=1.0e-12\
 --relative-tolerance=$rel_error\
 "$semi_analytic_output" "$two_scale_output"`

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
