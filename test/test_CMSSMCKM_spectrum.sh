#!/bin/sh

BASEDIR=$(dirname $0)
UTILSDIR=${BASEDIR}/../utils
print_block="$UTILSDIR/print_slha_block.awk"
print_block_entry="$UTILSDIR/print_slha_block_entry.awk"
remove_block="$UTILSDIR/remove_slha_block.awk"

slha_in="$BASEDIR/test_CMSSMCKM_spectrum.spc.in"
FS_out="$BASEDIR/test_CMSSMCKM_spectrum.spc.out.FS"
SS_out="$BASEDIR/test_CMSSMCKM_spectrum.spc.out.SS"

sed_cmd=`command -v sed`
awk_cmd=`command -v awk`
numdiff_cmd=`command -v numdiff`
FS="$BASEDIR/../models/CMSSMCKM/run_CMSSMCKM.x"
SS="$BASEDIR/../test/SOFTSUSY/run_softpoint.x"
rel_error="3e-3"

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

if test ! -x "$FS"; then
    echo "Error: CMSSMCKM spectrum generator not found: $FS"
    exit 1
fi

if test ! -x "$SS"; then
    echo "Error: SOFTSUSY spectrum generator not found: $SS"
    exit 1
fi

remove_unused_blocks() {
    $awk_cmd -f "$print_block" -v block=VCKM
}

echo "running FS ..."
${FS} --slha-input-file="$slha_in" --slha-output-file="$FS_out"
echo "done ($?)"

echo "running SS ..."
${SS} leshouches < "$slha_in" > "$SS_out"
echo "done ($?)"

$sed_cmd -i -e 's/#.*//' $FS_out
$sed_cmd -i -e 's/#.*//' $SS_out

cp $FS_out $FS_out~
cp $SS_out $SS_out~

remove_unused_blocks < $FS_out~ > $FS_out
remove_unused_blocks < $SS_out~ > $SS_out

diff=`$numdiff_cmd\
 --absolute-tolerance=1.0e-12\
 --relative-tolerance=$rel_error\
 $FS_out $SS_out`

diff_without_comments=`echo $diff | $sed_cmd -e '/^ *#/d' | $sed_cmd -e '/^+++/d'`

exit_code=0

if [ -n "$diff_without_comments" ]; then
    echo "Error: difference between ${FS_out} and  ${SS_out} larger than ${rel_error}"
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
