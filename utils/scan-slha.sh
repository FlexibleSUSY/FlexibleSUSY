#!/bin/sh

# This script allows one to scan over one parameter given an SLHA
# input file and a spectrum generator
#
# Examples:
#
# ./utils/scan-slha.sh \
#    --spectrum-generator=models/CMSSM/run_CMSSM.x \
#    --slha-input-file=model_files/CMSSM/LesHouches.in.CMSSM \
#    --scan-range=MINPAR[1]=100~300:10 \
#    --output=MINPAR[1],MASS[25],Yu[3:3]
#
# cat model_files/CMSSM/LesHouches.in.CMSSM | \
# ./utils/scan-slha.sh \
#    --spectrum-generator=models/CMSSM/run_CMSSM.x \
#    --scan-range=MINPAR[1]=100~300:10 \
#    --output=MINPAR[1],MASS[25],Yu[3:3]
#
# Author: Alexander Voigt

output=
scan_range=
sg_type=FlexibleSUSY
slha_input=
slha_input_file=
spectrum_generator=
step_size="linear"
AWK=${AWK:-awk}

# prints SLHA block
print_slha_block_awk='
BEGIN {
   is_block = 0;
   if (block == "") {
      print "Error: block name not defined";
      print "   Please define the block name with -v block=<block-name>";
      exit 1
   }
}
{
   pattern     = "^block[ \t\n\r\f]*" tolower(block) "([^a-zA-Z0-9_].*)?$";
   not_pattern = "^block[ \t\n\r\f]*.*$";

   if (tolower($0) ~ pattern) {
      is_block = 1
   } else if (tolower($0) ~ not_pattern) {
      is_block = 0
   };

   if (is_block)
      print $0
}
'

# prints block entry
# expects block entry keys in the form x or x:y or x:y:z etc.
print_block_entry_awk='
{
  if ($0 ~ /^B/) # skip block head
     next

  len = split(keys,k,":");

  matches = 1;

  for (i in k) {
     if ($(i) != k[i])
        matches = 0
  }

  if (matches == 1)
     print $(len + 1)
}
'

# string to eval when script terminates, normally or not
actions_at_exit=""

#_____________________________________________________________________
at_exit() {
    local action
    for action in "$@"; do
	printf "%s\n" "$actions_at_exit" | grep -F -e "$action" > /dev/null ||
	actions_at_exit="$actions_at_exit
$action"
    done
}

#_____________________________________________________________________
do_actions_at_exit() {
    eval "$actions_at_exit"
}

#_____________________________________________________________________
help() {
cat <<EOF
Usage: ./`basename $0` [options]
Options:

  --output=             Comma separated list of output values
                        Syntax: <block>[<entry>] | ?<block>[<entry>]
                        Example: MINPAR[1],MASS[25],Yu[3:3]
                                 ?SPINFO[3]   # 1 if SPINFO[3] is present
  --scan-range=         Scan range
                        Syntax: <block>[<entry>]=<start>~<stop>:<steps-1>
                        Example: MINPAR[1]=100~1000:10
  --slha-input-file=    SLHA input file (optional).
                        If no SLHA input file is given, the SLHA input is
                        read from stdin .
  --spectrum-generator= Spectrum generator executable
  --step-size=          the step size (linear or log)
  --type=               Spectrum generator type (default: ${sg_type})
                        Possible values: FlexibleSUSY SOFTSUSY SPheno SuSpect
  --help,-h             Print this help message

Examples:

   $ ./utils/scan-slha.sh \\
        --spectrum-generator=models/CMSSM/run_CMSSM.x \\
        --slha-input-file=model_files/CMSSM/LesHouches.in.CMSSM \\
        --scan-range=MINPAR[3]=1~30:21 \\
        --output=MINPAR[3],MASS[25],Yu[3:3] \\
     > scan-slha.dat

   $ cat model_files/CMSSM/LesHouches.in.CMSSM | \\
     ./utils/scan-slha.sh \\
        --spectrum-generator=models/CMSSM/run_CMSSM.x \\
        --scan-range=MINPAR[3]=1~30:21 \\
        --output=MINPAR[3],MASS[25],Yu[3:3] \\
     > scan-slha.dat

   $  echo "set xlabel \"tan(beta)\"; 
            set ylabel \"mh / GeV\";
            plot 'scan-slha.dat' u 1:2 w linespoints" \\
        | gnuplot -p 
EOF
}

#_____________________________________________________________________
run_flexiblesusy() {
    local sg="$1"
    local input="$2"

    echo "$input" | "$sg" --slha-input-file=- 2>/dev/null
}

#_____________________________________________________________________
run_softsusy() {
    local sg="$1"
    local input="$2"

    echo "$input" | "$sg" leshouches 2>/dev/null
}

#_____________________________________________________________________
run_spheno() {
    local sg="$1"
    local input="$2"
    local tmp_in="$$.in"
    local tmp_out="$$.spc"

    rm -f "$tmp_out" "$tmp_in"
    echo "$input" > "$tmp_in"

    "$sg" "$tmp_in" "$tmp_out" > /dev/null 2>&1

    if test "x$?" = "x0" -a -f "$tmp_out" ; then
        cat "$tmp_out"
    fi

    rm -f "$tmp_out" "$tmp_in" Messages.out SPheno.out \
       WHIZARD.par.* effC.dat BR_t.dat BR_Hplus.dat BR_H_NP.dat \
       LEP_HpHm_CS_ratios.dat MH_GammaTot.dat MHplus_GammaTot.dat fort.10
}

#_____________________________________________________________________
run_suspect() {
    local sg="$1"
    local input="$2"
    local tmp_in="suspect2_lha.in"
    local tmp_out="suspect2_lha.out"

    rm -f "$tmp_out" "$tmp_in"
    echo "$input" > "$tmp_in"

    "$sg" "$tmp_in" "$tmp_out" > /dev/null 2>&1

    if test "x$?" = "x0" -a -f "$tmp_out" ; then
        cat "$tmp_out"
    fi

    rm -f "$tmp_out" "$tmp_in" "suspect2.out"
}

#_____________________________________________________________________
run_sg() {
    local type="$1"
    local sg="$2"
    local input="$3"
    local func=

    case "$type" in
        FlexibleSUSY) func=run_flexiblesusy ;;
        SOFTSUSY)     func=run_softsusy ;;
        SPheno)       func=run_spheno ;;
        SuSpect)      func=run_suspect ;;
        *)
            echo "Error: unknown spectrum generator type: $type"
            exit 1
    esac

    ${func} "$sg" "$input"
}

trap do_actions_at_exit 0
trap "exit 1" INT QUIT TERM

if test $# -gt 0 ; then
    while test ! "x$1" = "x" ; do
        case "$1" in
            -*=*) optarg=`echo "$1" | sed 's/[-_a-zA-Z0-9]*=//'` ;;
            *) optarg= ;;
        esac

        case $1 in
            --output=*)              output=$optarg ;;
            --scan-range=*)          scan_range=$optarg ;;
            --slha-input-file=*)     slha_input_file=$optarg ;;
            --spectrum-generator=*)  spectrum_generator=$optarg ;;
            --step-size=*)           step_size=$optarg ;;
            --type=*)                sg_type=$optarg ;;
            --help|-h)               help; exit 0 ;;
            *)  echo "Invalid option '$1'. Try $0 --help" ; exit 1 ;;
        esac
        shift
    done
fi

test -z "$slha_input_file" -o -e "$slha_input_file" || \
    { echo "Error: input file $slha_input_file not found." ; exit 1; }

slha_input=`cat ${slha_input_file}`

# substitute ~ by $HOME
spectrum_generator=$(echo "$spectrum_generator" | sed 's|~|'"${HOME}"'|g')

if test ! -x "$spectrum_generator"; then
    echo "Error: spectrum generator executable $spectrum_generator not found."
    exit 1
fi

if test -z "$output"; then
    echo "Error: no output values given.  Please provide them via --output= "
    exit 1
fi

if test -z "$scan_range"; then
    echo "Error: no scan range given.  Please provide it via --scan-range= "
    exit 1
fi

# transform scientific notation into bc syntax
start=$(echo "$scan_range" | ${AWK} -F '[=:~]' '{ print $2 }' | sed -e 's/[eE]+*/*10^/')
stop=$(echo "$scan_range"  | ${AWK} -F '[=:~]' '{ print $3 }' | sed -e 's/[eE]+*/*10^/')
steps=$(echo "$scan_range" | ${AWK} -F : '{ print $NF }'      | sed -e 's/[eE]+*/*10^/')
block=$(echo "$scan_range" | ${AWK} -F [ '{ print $1 }')
entry=$(echo "$scan_range" | ${AWK} -F '[][]' '{ print $2 }')

output_fields="$(echo $output | tr ',' ' ')"

# print comment line
if test "$steps" -gt 0; then
    comment="#"
    args=" "
    for f in $output_fields; do
        comment="$comment %16s"
        args="$args $f"
    done
    comment="$comment\n"
    printf "$comment" $args
fi

# start scan over points
for i in `seq 0 $steps`; do
    # calculate current value for the scanned variable
    case "$step_size" in
        linear)
            value=$(cat <<EOF | bc
scale=20
$start + ($stop - $start)*${i} / $steps
EOF
                 ) ;;
        log)
            value=$(cat <<EOF | bc -l
scale=20
e(l($start) + (l($stop) - l($start))*${i} / $steps)
EOF
                 ) ;;
        *) echo "Error: unknown step size: $step_size"
           exit 1 ;;
    esac

    # create input data
    slha_input=$(
    { echo "$slha_input" ; \
      cat <<EOF
Block $block # added by `basename $0`
  $entry    $value
EOF
    })

    # run the spectrum generator
    slha_output=$(run_sg "$sg_type" "$spectrum_generator" "$slha_input")

    printfstr=" "
    args=

    # get the output
    for f in $output_fields; do
        output_block=$(echo "$f" | ${AWK} -F [ '{ print $1 }')

        # do we need to check only whether the entry exists?
        echo "$f" | grep -v '?' > /dev/null
        check_present=$?

        output_block=$(echo "$output_block" | sed -e 's/?//')
        full_block=$(echo "$slha_output" | ${AWK} -v block="$output_block" "$print_slha_block_awk")
        block_entries=$(echo "$f" | ${AWK} -F '[][]' '{ print $2 }')

        # get data value
        value=$(echo "$full_block" | ${AWK} -v keys="$block_entries" "$print_block_entry_awk" | tail -n 1)

        if test "$check_present" -eq 1 ; then
            # check if entry exists
            # if entry exists, set value=1, otherwise set value=0
            test -z "$value" > /dev/null
            value=$?
        fi

        if test -z "$value"; then
            value="-"
        fi

        printfstr="$printfstr %16s"
        args="$args $value"
    done

    # print data line
    printfstr="$printfstr\n"
    printf "$printfstr" $args

done
