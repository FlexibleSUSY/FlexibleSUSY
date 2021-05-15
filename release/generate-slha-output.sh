#!/bin/sh

# creates SLHA output files for all SLHA input files in the directory
# model_files/ (and sub-directories) that begin with LesHouches.in. .
# The names of the output files will begin with LesHouches.out. .

# Author: Alexander Voigt

# directory of this script
BASEDIR=$(dirname $0)
HOMEDIR=$(readlink -f "${BASEDIR}/../")
FSCONFIG="${HOMEDIR}/flexiblesusy-config"
model_file_dir="$BASEDIR/../model_files"
directory=.

#_____________________________________________________________________
help() {
cat <<EOF
Usage: ./`basename $0` [options]
Options:

  --directory=         output directory (default: ${directory})
                       If empty, the SLHA output file will be written
                       to the same directory as the input file
  --help,-h            Print this help message
EOF
}

if test $# -gt 0 ; then
    while test ! "x$1" = "x" ; do
        case "$1" in
            -*=*) optarg=`echo "$1" | sed 's/[-_a-zA-Z0-9]*=//'` ;;
            *) optarg= ;;
        esac

        case $1 in
            --directory=*)           directory=$optarg ;;
            --help|-h)               help; exit 0 ;;
            *)  echo "Invalid option '$1'. Try $0 --help" ; exit 1 ;;
        esac
        shift
    done
fi

[ -z "${directory}" ] && directory=.

[ ! -d "${directory}" ] && mkdir -p "${directory}"

SGs=$(find $model_file_dir/ -type f -iname LesHouches.in.\* -not -iname \*~ -exec dirname {} \; | awk -F / '{ print $NF }' | uniq)

errors=0

echo "Found default SLHA input files for: $SGs"

for sg in ${SGs}
do
    if [ $("$FSCONFIG" --with-${sg}) = no ] ; then
        continue
    fi

    exe="${HOMEDIR}/models/${sg}/run_${sg}.x"

    echo "========================"
    echo "   $sg"
    echo "========================"

    input_files=$(find $model_file_dir/$sg/ -type f -iname LesHouches.in.\* -not -iname \*~)
    echo "input files: "
    echo "$input_files"

    for ifile in ${input_files}
    do
        ofile=$(echo ${ifile} | sed -e 's/LesHouches\.in\./LesHouches.out./')

        if test ! "x$directory" = "x"; then
            ofile="${directory}/`basename $ofile`"
        fi

        cmd="${exe} --slha-input-file=${ifile} --slha-output-file=${ofile} > /dev/null 2>&1"

        echo ""
        echo "> running $sg"
        echo "> input file: `basename ${ifile}`"
        echo "> output file: `basename ${ofile}`"
        echo "> command: ${cmd}"

        eval "${cmd}"
        exit_code="$?"

        if test ${exit_code} -eq 0 ; then
            echo "> OK"
        else
            echo "> FAIL"
            errors=1
        fi
    done
done

echo ""
if test ${errors} -ne 0 ; then
    echo "There were errors!"
else
    echo "All output files generated successfully!"
fi

exit $errors
