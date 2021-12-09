#!/bin/sh

# This script allows one to scan over one parameter given an SLHA
# input file and a spectrum generator and writing the output to a
# database.
#
# Examples:
#
# ./utils/scan-database.sh \
#    --spectrum-generator=models/CMSSM/run_CMSSM.x \
#    --slha-input-file=model_files/CMSSM/LesHouches.in.CMSSM \
#    --scan-range=MINPAR[1]=100~300:10 \
#    --database-output-file=scan.db
#
# cat model_files/CMSSM/LesHouches.in.CMSSM | \
# ./utils/scan-database.sh \
#    --spectrum-generator=models/CMSSM/run_CMSSM.x \
#    --scan-range=MINPAR[1]=100~300:10 \
#    --database-output-file=scan.db
#
# Author: Alexander Voigt

database_output_file=
scan_range=
slha_input=
slha_input_file=
spectrum_generator=
step_size="linear"

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

  --database-output-file=<filename> Name of database output file
  --scan-range=         Scan range
                        Syntax: <block>[<entry>]=<start>~<stop>:<steps-1>
                        Example: MINPAR[1]=100~1000:10
  --slha-input-file=    SLHA input file (optional).
                        If no SLHA input file is given, the SLHA input is
                        read from stdin .
  --spectrum-generator= Spectrum generator executable
  --step-size=          the step size (linear or log)
  --help,-h             Print this help message

Examples:

   $ ./utils/scan-database.sh \\
        --spectrum-generator=models/CMSSM/run_CMSSM.x \\
        --slha-input-file=model_files/CMSSM/LesHouches.in.CMSSM \\
        --scan-range=MINPAR[3]=1~30:21 \\
        --database-output-file=scan.db

   $ cat model_files/CMSSM/LesHouches.in.CMSSM | \\
     ./utils/scan-database.sh \\
        --spectrum-generator=models/CMSSM/run_CMSSM.x \\
        --scan-range=MINPAR[3]=1~30:21 \\
        --database-output-file=scan.db
EOF
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
            --database-output-file=*) database_output_file=$optarg ;;
            --scan-range=*)          scan_range=$optarg ;;
            --slha-input-file=*)     slha_input_file=$optarg ;;
            --spectrum-generator=*)  spectrum_generator=$optarg ;;
            --step-size=*)           step_size=$optarg ;;
            --help|-h)               help; exit 0 ;;
            *)  echo "Invalid option '$1'. Try $0 --help" ; exit 1 ;;
        esac
        shift
    done
fi

test -z "$slha_input_file" -o -e "$slha_input_file" || \
    { echo "Error: input file $slha_input_file not found." ; exit 1; }

slha_input=`cat ${slha_input_file}`

if test ! -x "$spectrum_generator"; then
    echo "Error: spectrum generator executable $spectrum_generator not found."
    exit 1
fi

if test -z "$database_output_file"; then
    echo "Error: no database output file given.  Please provide one via --database-output-file= "
    exit 1
fi

if test -z "$scan_range"; then
    echo "Error: no scan range given.  Please provide it via --scan-range= "
    exit 1
fi

# transform scientific notation into bc syntax
start=$(echo "$scan_range" | awk -F '[=:~]' '{ print $2 }' | sed -e 's/[eE]+*/*10^/')
stop=$(echo "$scan_range"  | awk -F '[=:~]' '{ print $3 }' | sed -e 's/[eE]+*/*10^/')
steps=$(echo "$scan_range" | awk -F : '{ print $NF }'      | sed -e 's/[eE]+*/*10^/')
block=$(echo "$scan_range" | awk -F [ '{ print $1 }')
entry=$(echo "$scan_range" | awk -F '[][]' '{ print $2 }')

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

    # run the spectrum generator
    # shellcheck disable=SC2034
    slha_output=$(
    { echo "$slha_input" ; \
      cat <<EOF
Block $block # added by `basename $0`
  $entry    $value
EOF
    } | $spectrum_generator \
        --slha-input-file=- \
        --slha-output-file= \
        --database-output-file="$database_output_file")

    echo "running ${spectrum_generator} with ${block}[${entry}] = ${value}"
done
