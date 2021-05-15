#!/bin/sh -e

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

  --directory=         Output directory (default: ${directory})
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

default_input_files="\
models/CMSSM/LesHouches.in.CMSSM \
models/CMSSMSemiAnalytic/LesHouches.in.CMSSMSemiAnalytic \
models/MSSM/LesHouches.in.MSSM \
models/MSSMatMGUT/LesHouches.in.MSSMatMGUT \
models/MSSMNoFV/LesHouches.in.MSSMNoFV \
models/MSSMNoFVatMGUT/LesHouches.in.MSSMNoFVatMGUT \
models/CMSSMNoFV/LesHouches.in.CMSSMNoFV \
models/NUHMSSM/LesHouches.in.NUHMSSM \
models/lowMSSM/LesHouches.in.lowMSSM \
models/MSSMRHN/LesHouches.in.MSSMRHN_generated \
models/NMSSM/LesHouches.in.NMSSM \
models/NUTNMSSM/LesHouches.in.NUTNMSSM \
models/NUTSMSSM/LesHouches.in.NUTSMSSM \
models/lowNMSSM/LesHouches.in.lowNMSSM \
models/lowNMSSMTanBetaAtMZ/LesHouches.in.lowNMSSMTanBetaAtMZ \
models/SMSSM/LesHouches.in.SMSSM \
models/UMSSM/LesHouches.in.UMSSM \
models/E6SSM/LesHouches.in.E6SSM \
models/MRSSM/LesHouches.in.MRSSM \
models/TMSSM/LesHouches.in.TMSSM \
models/SM/LesHouches.in.SM \
models/HSSUSY/LesHouches.in.HSSUSY \
models/SplitMSSM/LesHouches.in.SplitMSSM \
models/THDMII/LesHouches.in.THDMII \
models/THDMIIMSSMBC/LesHouches.in.THDMIIMSSMBC \
models/HTHDMIIMSSMBC/LesHouches.in.HTHDMIIMSSMBC \
models/HGTHDMIIMSSMBC/LesHouches.in.HGTHDMIIMSSMBC \
models/MSSMEFTHiggs/LesHouches.in.MSSMEFTHiggs \
models/NMSSMEFTHiggs/LesHouches.in.NMSSMEFTHiggs \
models/E6SSMEFTHiggs/LesHouches.in.E6SSMEFTHiggs \
models/MRSSMEFTHiggs/LesHouches.in.MRSSMEFTHiggs \
models/CNMSSM/LesHouches.in.CNMSSM \
models/CE6SSM/LesHouches.in.CE6SSM \
models/MSSMNoFVatMGUTHimalaya/LesHouches.in.MSSMNoFVatMGUTHimalaya \
models/MSSMNoFVHimalaya/LesHouches.in.MSSMNoFVHimalaya \
models/NUHMSSMNoFVHimalaya/LesHouches.in.NUHMSSMNoFVHimalaya \
"

[ -z "${directory}" ] && directory=.

[ ! -d "${directory}" ] && mkdir -p "${directory}"

models=$("$FSCONFIG" --models)

echo "Configured models: $models"
echo

for model in ${models}; do
    # collect all input files that belong to the model
    input_files=
    for dif in ${default_input_files}; do
        case "${dif}" in
            ${model}/*) input_files="${input_files} ${HOMEDIR}/${dif}" ;;
        esac
    done

    for ifile in ${input_files}; do
        ofile=$(echo "${directory}/$(basename ${ifile})" | sed -e 's/\.in\./.out./')
        sg=$(echo "${model}" | awk -F / '{ print $NF }')
        exe="${HOMEDIR}/${model}/run_${sg}.x"
        cmd="${exe} --slha-input-file=${ifile} --slha-output-file=${ofile} > /dev/null 2>&1"

        printf "%s" "${cmd}"
        eval "${cmd}" || { printf " [FAIL]\n"; exit 1; }
        printf " [OK]\n"
    done
done
