#!/bin/sh

# This script generates a low-energy NUHMSSM spectrum using the
# semi-analytic solver and writes it to an SLHA file.  The EWSB output
# parameters from this point are then used as inputs for the model
# solved using the two-scale solver.  After running the two-scale
# lowNUHMSSM we compare that the semi-analytic and two-scale spectrums
# are the same.

BASEDIR=$(dirname $0)
UTILSDIR="${BASEDIR}/../utils"
print_block="$UTILSDIR/print_slha_block.awk"
print_block_entry="$UTILSDIR/print_slha_block_entry.awk"
remove_block="$UTILSDIR/remove_slha_block.awk"

semi_analytic_model="lowNUHMSSMSemiAnalytic"
two_scale_model="lowNUHMSSM"

semi_analytic_slha="${BASEDIR}/../model_files/${semi_analytic_model}/LesHouches.in.${semi_analytic_model}"
semi_analytic_input="${BASEDIR}/test_${semi_analytic_model}_spectrum.in.spc"
semi_analytic_output="${BASEDIR}/test_${semi_analytic_model}_spectrum.out.spc"
two_scale_input="${BASEDIR}/test_${semi_analytic_model}_spectrum_${two_scale_model}.in.spc"
two_scale_output="${BASEDIR}/test_${semi_analytic_model}_spectrum_${two_scale_model}.out.spc"
rel_error="1.7e-2"

sed_cmd=`command -v sed`
awk_cmd=`command -v awk`
numdiff_cmd=`command -v numdiff`
semi_analytic_exe="${BASEDIR}/../models/${semi_analytic_model}/run_${semi_analytic_model}.x"
two_scale_exe="${BASEDIR}/../models/${two_scale_model}/run_${two_scale_model}.x"

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

# note: blocks for soft trilinears are also removed since, if they
# come out exactly zero, they generate spurious large relative errors
remove_mixing_and_input_blocks() {
    $awk_cmd -f "$remove_block" -v block=UMIX \
        | $awk_cmd -f "$remove_block" -v block=VMIX \
        | $awk_cmd -f "$remove_block" -v block=PSEUDOSCALARMIX \
        | $awk_cmd -f "$remove_block" -v block=DSQMIX \
        | $awk_cmd -f "$remove_block" -v block=SELMIX \
        | $awk_cmd -f "$remove_block" -v block=SCALARMIX \
        | $awk_cmd -f "$remove_block" -v block=NMIX \
        | $awk_cmd -f "$remove_block" -v block=CHARGEMIX \
        | $awk_cmd -f "$remove_block" -v block=USQMIX \
        | $awk_cmd -f "$remove_block" -v block=SNUMIX \
        | $awk_cmd -f "$remove_block" -b block=Tu \
        | $awk_cmd -f "$remove_block" -b block=Td \
        | $awk_cmd -f "$remove_block" -b block=Te \
        | $awk_cmd -f "$remove_block" -v block=SMINPUTS \
        | $awk_cmd -f "$remove_block" -v block=MINPAR \
        | $awk_cmd -f "$remove_block" -v block=EXTPAR \
        | $awk_cmd -f "$remove_block" -v block=SPINFO \
        | $awk_cmd -f "$remove_block" -v block=MODSEL \
        | $awk_cmd -f "$remove_block" -v block=FlexibleSUSY
}

sign() {
    local a
    a=$(echo "$1" | sed -e 's/[eE]+*/*10^/')
    cat <<EOF | bc -l
define sign(i) {
    if (i >= 0) return(1)
    return(-1)
}

sign($a)
EOF
}

match_input_blocks() {
    if test $# -lt 4; then
        echo "match_input_blocks: Too few arguments"
        return 1
    fi
    _input_file="$1"
    _output_file="$2"
    _mHd2In_val="$3"
    _mHu2In_val="$4"

    _mu_val=$($awk_cmd -f "$print_block" -v block=EXTPAR "$_input_file" \
                  | $awk_cmd -f "$print_block_entry" -v entries=23)
    _sign_mu=$(sign $_mu_val)

    $awk_cmd -f "$print_block" -v block=FlexibleSUSY "$_input_file" \
             > "$_output_file"

    $awk_cmd -f "$remove_block" -v block=FlexibleSUSY "$_input_file" \
        | $awk_cmd -f "$remove_block" -v block=EXTPAR \
                   >> "$_output_file"

    echo "   4   $_sign_mu" >> "$_output_file"
    echo "Block EXTPAR" >> "$_output_file"
    echo "   1   $_mHd2In_val" >> "$_output_file"
    echo "   2   $_mHu2In_val" >> "$_output_file"
}

# generate point using semi-analytic solver
cp "$semi_analytic_slha" "$semi_analytic_input"
echo -n "running semi-analytic solver ..."
$semi_analytic_exe --slha-input-file="$semi_analytic_input" \
                   --slha-output-file="$semi_analytic_output"
echo "done"
echo "$semi_analytic_model SLHA input file: $semi_analytic_input"
echo "$semi_analytic_model SLHA output file: $semi_analytic_output"

if test ! -r "$semi_analytic_output"; then
    echo "Error: generated $semi_analytic_model spectrum not found: $semi_analytic_output"
    exit 1
fi

# remove comments from the semi-analytic output spectrum
$sed_cmd -i~ -e '/^ *#/d' "$semi_analytic_output"

# construct matching input point for the two-scale solver
mHd2_sol=$(awk -f "$print_block" -v block=MSOFT "$semi_analytic_output" \
               | awk -f "$print_block_entry" -v entries=21)
mHu2_sol=$(awk -f "$print_block" -v block=MSOFT "$semi_analytic_output" \
               | awk -f "$print_block_entry" -v entries=22)
match_input_blocks "$semi_analytic_input" "$two_scale_input" "$mHd2_sol" "$mHu2_sol"

# generate point using the two-scale solver
echo -n "running two-scale solver ..."
$two_scale_exe --slha-input-file="$two_scale_input" \
               --slha-output-file="$two_scale_output"
echo "done"
echo "$two_scale_model SLHA input file: $two_scale_input"
echo "$two_scale_model SLHA output file: $two_scale_output"

# remove comments from two-scale output spectrum
$sed_cmd -i~ -e '/^ *#/d' "$two_scale_output"

# remove mixing matrix blocks because we don't want to compare objects
# with phase ambiguities, and input blocks since these do not match
cp "$semi_analytic_output" "$semi_analytic_output~"
cp "$two_scale_output" "$two_scale_output~"

remove_mixing_and_input_blocks < "$semi_analytic_output~" \
                               > "$semi_analytic_output"
remove_mixing_and_input_blocks < "$two_scale_output~" \
                               > "$two_scale_output"

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
