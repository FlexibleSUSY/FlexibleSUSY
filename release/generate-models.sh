#!/bin/sh -e

# creates models for public release

himalaya_lib_dir=
himalaya_inc_dir=
models=
number_of_jobs=1
directory=.
MATH=math

#_____________________________________________________________________
help() {
cat <<EOF
Usage: ./`basename $0` [options]
Options:

  --with-himalaya-libdir=   Path to search for Himalaya library
  --with-himalaya-incdir=   Path to search for Himalaya header
  --with-models=            Comma separated list of models to be generated
  --number-of-jobs=         Number of parallel makefile jobs
  --directory=              Output directory (default: ${directory})
  --with-math-cmd=          Mathematic kernel (default: $MATH)
  --help,-h                 Print this help message
EOF
}

if test $# -gt 0 ; then
    while test ! "x$1" = "x" ; do
        case "$1" in
            -*=*) optarg=`echo "$1" | sed 's/[-_a-zA-Z0-9]*=//'` ;;
            *) optarg= ;;
        esac

        case $1 in
            --with-himalaya-incdir=*) himalaya_inc_dir=$optarg ;;
            --with-himalaya-libdir=*) himalaya_lib_dir=$optarg ;;
            --with-models=*)          models=$optarg ;;
            --number-of-jobs=*)       number_of_jobs=$optarg ;;
            --directory=*)            directory=$optarg ;;
            --with-math-cmd=*)        MATH=$optarg ;;
            --help|-h)                help; exit 0 ;;
            *)  echo "Invalid option '$1'. Try $0 --help" ; exit 1 ;;
        esac
        shift
    done
fi

echo "Using $number_of_jobs parallel jobs"

# use default models if no models specified
models="${models:-\
CMSSM,\
CMSSMSemiAnalytic,\
MSSM,\
MSSMatMGUT,\
MSSMNoFV,\
MSSMNoFVatMGUT,\
CMSSMNoFV,\
NUHMSSM,\
lowMSSM,\
MSSMRHN,\
NMSSM,\
NUTNMSSM,\
NUTSMSSM,\
lowNMSSM,\
lowNMSSMTanBetaAtMZ,\
SMSSM,\
UMSSM,\
E6SSM,\
MRSSM,\
TMSSM,\
SM,\
HSSUSY,\
SplitMSSM,\
THDMII-II,\
THDMIIMSSMBC-II,\
HTHDMIIMSSMBC-II,\
HGTHDMIIMSSMBC-II,\
MSSMEFTHiggs,\
NMSSMEFTHiggs,\
E6SSMEFTHiggs,\
MRSSMEFTHiggs,\
CNMSSM,\
CE6SSM,\
MSSMNoFVatMGUTHimalaya,\
MSSMNoFVHimalaya,\
NUHMSSMNoFVHimalaya,\
}"

echo "Building models: ${models}"

models_space=$(echo $models | tr ',' ' ')

# creating models
for m in ${models_space}; do
    ./createmodel --name=${m} --force  --with-math-cmd=${MATH}
done

./configure \
    --with-himalaya-libdir="${himalaya_lib_dir}" \
    --with-himalaya-incdir="${himalaya_inc_dir}" \
    --with-models=${models} \
    --with-math-cmd=${MATH}

make showbuild

make -j${number_of_jobs}

# packing models
for m in ${models_space}; do
    echo "packing ${m} ..."
    make pack-${m}-src
done

# moving models
[ -z "${directory}" ] && directory=.

[ ! -d "${directory}" ] && mkdir -p "${directory}"

if [ "x${directory}" != "x." ]; then
    for m in ${models_space}; do
        mv ${m}.tar.gz ${directory}/
    done
fi
