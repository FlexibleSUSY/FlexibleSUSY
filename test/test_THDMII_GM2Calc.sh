#!/bin/sh

# directory of this script
BASEDIR=$(dirname $0)
FSCONFIG="$BASEDIR/../flexiblesusy-config"
THDMII_EXE="${BASEDIR}/../models/THDMII/run_THDMII.x"
SLHA_OUT="${BASEDIR}/test_THDMII_GM2Calc.out.spc"
print_block="$BASEDIR/../utils/print_slha_block.awk"

SMINPUTS="\
Block SMINPUTS               # Standard Model inputs
    1     1.27934000E+02   # alpha_em(MZ)^(-1) SM MS-bar
    3     0.1184           # alpha_s(MZ) SM MS-bar [2L]
    4     9.11876000E+01   # MZ(pole)              [1L]
    5     4.18000000E+00   # mb(mb) SM MS-bar      [2L]
    6     1.73340000E+02   # mtop(pole)            [2L]
    7     1.77700000E+00   # mtau(pole)            [2L]
    8     0.00000000E+00   # mnu3(pole)            [irrelevant]
    9     8.03773317E+01   # mW(pole)              [1L]
   11     0.000510998928   # melectron(pole)       [2L]
   12     0.00000000E+00   # mnu1(pole)            [irrelevant]
   13     0.1056583715     # mmuon(pole)           [1L]
   14     0.00000000E+00   # mnu2(pole)            [irrelevant]
   21     4.76052706E-03   # md(2 GeV)             [2L]
   22     2.40534062E-03   # mu(2 GeV)             [2L]
   23     1.04230487E-01   # ms(2 GeV)             [2L]
   24     1.27183378E+00   # mc(2 GeV)             [2L]
"

SLHA_IN="\
Block MODSEL                 # Select model
#   12    1000                # parameter output scale (GeV)
Block FlexibleSUSY
    0   1.000000000e-04      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = all, 1 = two_scale, 2 = semi_analytic)
    3   1                    # calculate SM pole masses
    4   1                    # pole mass loop order
    5   1                    # EWSB loop order
    6   2                    # beta-functions loop order
    7   2                    # threshold corrections loop order
    8   0                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   0                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   0                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   0                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   0                    # force output
   13   1                    # Top pole mass QCD corrections (0 = 1L, 1 = 2L, 2 = 3L)
   14   1.000000000e-11      # beta-function zero threshold
   15   1                    # calculate observables (a_muon, ...)
   16   0                    # force positive majorana masses
   17   0                    # pole mass renormalization scale (0 = SUSY scale)
   18   0                    # pole mass renormalization scale in the EFT (0 = min(SUSY scale, Mt))
   19   0                    # EFT matching scale (0 = SUSY scale)
   20   2                    # EFT loop order for upwards matching
   21   1                    # EFT loop order for downwards matching
   22   0                    # EFT index of SM-like Higgs in the BSM model
   23   1                    # calculate BSM pole masses
   24   123111321            # individual threshold correction loop orders
   25   0                    # ren. scheme for Higgs 3L corrections (0 = DR, 1 = MDR)
   26   1                    # Higgs 3-loop corrections O(alpha_t alpha_s^2)
   27   1                    # Higgs 3-loop corrections O(alpha_b alpha_s^2)
   28   1                    # Higgs 3-loop corrections O(alpha_t^2 alpha_s)
   29   1                    # Higgs 3-loop corrections O(alpha_t^3)
   30   1                    # Higgs 4-loop corrections O(alpha_t alpha_s^3)
   31   0                    # loop library (0 = softsusy)
Block FlexibleSUSYInput
    0   0.00729735           # alpha_em(0)
    1   125.09               # Mh pole
${SMINPUTS}
Block MINPAR                 # Input parameters
    1   2                    # Lambda1IN
    2   0.1                  # Lambda2IN
    3   0.5                  # Lambda3IN
    4   0.8                  # Lambda4IN
    5   -2                   # Lambda5IN
    6   0                    # Lambda6IN
    7   0                    # Lambda7IN
    8   1000                 # M122
    10  10                   # TanBeta
Block EXTPAR
    0   9.11876000E+01       # Qin
"

GM2CALC_IN="\
Block GM2CalcConfig
     0     0     # output format (0 = minimal, 1 = detailed,
                 #  2 = NMSSMTools, 3 = SPheno, 4 = GM2Calc)
     1     2     # loop order (0, 1 or 2)
     3     0     # force output (0 or 1)
     4     0     # verbose output (0 or 1)
     5     0     # calculate uncertainty (0 or 1)
     6     1     # running couplings in the THDM
${SMINPUTS}
Block GM2CalcInput
    33     125.09           # SM Higgs boson mass   [1L]
Block MINPAR                # model parameters
     3     10               # tan(beta)
    11     4                # lambda_1
    12     0.2              # lambda_2
    13     0.5              # lambda_3
    14     0.8              # lambda_4
    15     -2               # lambda_5
    16     0                # lambda_6
    17     0                # lambda_7
    18     1000             # m_{12}^2
    24     2                # Yukawa type (1, 2, 3, 4, 5, 6 = general)
"

if [ -e GM2Calc.pc ] ; then
    # shellcheck disable=SC2046
    eval $(grep '^prefix=' GM2Calc.pc)
    # shellcheck disable=SC2154
    GM2CALC_EXE="${prefix}/bin/gm2calc.x"
elif [ -e gm2calc.pc ] ; then
    # shellcheck disable=SC2046
    eval $(grep '^prefix=' gm2calc.pc)
    # shellcheck disable=SC2154
    GM2CALC_EXE="${prefix}/bin/gm2calc.x"
fi

if [ ! -x "${GM2CALC_EXE}" ] ; then
    echo "Cannot find GM2Calc executable in ${GM2CALC_EXE}"
    echo "Skipping test."
    exit
fi

[ "$("$FSCONFIG" --with-THDMII)" = yes ] && [ -x "${THDMII_EXE}" ] || {
    echo "Error: THDMII needs to be build!"
    exit 1;
}

[ "$("$FSCONFIG" --enable-gm2calc)" = yes ] || {
    echo "Error: FlexibleSUSY must be configured with GM2Calc!"
    exit 1;
}

{ printf "%s\n" "${SLHA_IN}";
  cat <<EOF
Block FlexibleSUSY
    15  1   # calculate observables (a_muon, ...)
EOF
  } | "${THDMII_EXE}" --slha-input-file=- --slha-output-file="${SLHA_OUT}"

[ $? = 0 ] || {
    echo "Error: ${THDMII_EXE} failed!"
    exit 1
}

# amu from FS
amu_1l_fs=$(cat "${SLHA_OUT}" | awk -f "$print_block" -v block=FlexibleSUSYLowEnergy | awk '{ if ($1 == 21) print $2 }')

# Mh from FS
mh_fs=$(cat "${SLHA_OUT}" | awk -f "$print_block" -v block=MASS | awk '{ if ($1 == 25) print $2 }')

# MW from FS
mw_fs=$(cat "${SLHA_OUT}" | awk -f "$print_block" -v block=MASS | awk '{ if ($1 == 24) print $2 }')

# amu from GM2Calc, embedded in FS
amu_2l_gm2calc_fs=$(cat "${SLHA_OUT}" | awk -f "$print_block" -v block=FlexibleSUSYLowEnergy | awk '{ if ($1 == 26) print $2 }')

# uncertainty of amu from GM2Calc, embedded in FS
damu_2l_gm2calc_fs=$(cat "${SLHA_OUT}" | awk -f "$print_block" -v block=FlexibleSUSYLowEnergy | awk '{ if ($1 == 27) print $2 }')

# amu 2-loop from vanilla GM2Calc
amu_2l_gm2calc=$({ printf "%s\n" "${GM2CALC_IN}";
                cat <<EOF
Block GM2CalcConfig
     0  0  # minimal output
     1  2  # loop order (0, 1 or 2)
     4  0  # verbose output
Block SMINPUTS
     9     ${mw_fs}   # mW(pole)              [1L]
Block GM2CalcInput
    33     ${mh_fs}   # SM Higgs boson mass   [1L]
EOF
      } | "${GM2CALC_EXE}" --thdm-input-file=-)

[ $? = 0 ] || {
    echo "Error: ${GM2CALC_EXE} failed!"
    exit 1
}

# uncertainty of amu 2-loop from vanilla GM2Calc
damu_2l_gm2calc=$({ printf "%s\n" "${GM2CALC_IN}";
                cat <<EOF
Block GM2CalcConfig
     0  0  # minimal output
     1  2  # loop order (0, 1 or 2)
     4  0  # verbose output
     5  1  # calculate uncertainty (0 or 1)
Block SMINPUTS
     9     ${mw_fs}   # mW(pole)              [1L]
Block GM2CalcInput
    33     ${mh_fs}   # SM Higgs boson mass   [1L]
EOF
      } | "${GM2CALC_EXE}" --thdm-input-file=-)

[ $? = 0 ] || {
    echo "Error: ${GM2CALC_EXE} failed!"
    exit 1
}

# amu 1-loop from vanilla GM2Calc
amu_1l_gm2calc=$({ printf "%s\n" "${GM2CALC_IN}";
                cat <<EOF
Block GM2CalcConfig
     0  0  # minimal output
     1  1  # loop order (0, 1 or 2)
     4  0  # verbose output
Block SMINPUTS
     9     ${mw_fs}   # mW(pole)              [1L]
Block GM2CalcInput
    33     ${mh_fs}   # SM Higgs boson mass   [1L]
EOF
      } | "${GM2CALC_EXE}" --thdm-input-file=-)

[ $? = 0 ] || {
    echo "Error: ${GM2CALC_EXE} failed!"
    exit 1
}

# uncertainty of amu 1-loop from vanilla GM2Calc
damu_1l_gm2calc=$({ printf "%s\n" "${GM2CALC_IN}";
                cat <<EOF
Block GM2CalcConfig
     0  0  # minimal output
     1  1  # loop order (0, 1 or 2)
     4  0  # verbose output
     5  1  # calculate uncertainty (0 or 1)
Block SMINPUTS
     9     ${mw_fs}   # mW(pole)              [1L]
Block GM2CalcInput
    33     ${mh_fs}   # SM Higgs boson mass   [1L]
EOF
      } | "${GM2CALC_EXE}" --thdm-input-file=-)

[ $? = 0 ] || {
    echo "Error: ${GM2CALC_EXE} failed!"
    exit 1
}

# convert scientific notation to bc friendly notation
amu_1l_fs=$(echo "${amu_1l_fs}" | sed -e 's/[eE]/\*10\^/' | sed -e 's/\^+/\^/')
amu_2l_gm2calc_fs=$(echo "${amu_2l_gm2calc_fs}" | sed -e 's/[eE]/\*10\^/' | sed -e 's/\^+/\^/')
damu_2l_gm2calc_fs=$(echo "${damu_2l_gm2calc_fs}" | sed -e 's/[eE]/\*10\^/' | sed -e 's/\^+/\^/')
amu_1l_gm2calc=$(echo "${amu_1l_gm2calc}" | sed -e 's/[eE]/\*10\^/' | sed -e 's/\^+/\^/')
damu_1l_gm2calc=$(echo "${damu_1l_gm2calc}" | sed -e 's/[eE]/\*10\^/' | sed -e 's/\^+/\^/')
amu_2l_gm2calc=$(echo "${amu_2l_gm2calc}" | sed -e 's/[eE]/\*10\^/' | sed -e 's/\^+/\^/')
damu_2l_gm2calc=$(echo "${damu_2l_gm2calc}" | sed -e 's/[eE]/\*10\^/' | sed -e 's/\^+/\^/')

errors=0

# compares two values $1 and $2 for relative equality with a maximum relative deviation $3
test_close() {
    local val1
    local val2
    local rel_error
    local diff

    val1="$1"
    val2="$2"
    rel_error="$3"
    diff=$(cat <<EOF | bc $BASEDIR/abs.bc
scale=100
abs((abs(${val1}) - abs(${val2})) / (${val1})) < ${rel_error}
EOF
        )

    if test $diff -ne 1 ; then
        echo "Error: relative difference between"
        echo " ${val1} and ${val2} is larger than ${rel_error}"
        errors=1
    fi
}

### test 2L GM2Calc vs. embedded 2L GM2Calc
test_close "${amu_2l_gm2calc_fs}" "${amu_2l_gm2calc}" "0.0000001"

### test uncertainty 2L GM2Calc vs. embedded 2L GM2Calc
test_close "${damu_2l_gm2calc_fs}" "${damu_2l_gm2calc}" "0.0000001"

### test 1L GM2Calc vs. 1L FS
test_close "${amu_1l_fs}" "${amu_1l_gm2calc}" "0.3"

echo "FlexibleSUSY 1L + 2L QED: amu = ${amu_1l_fs}"
echo "original GM2Calc 1L     : amu = ${amu_1l_gm2calc} +/- ${damu_1l_gm2calc}"
echo "original GM2Calc 2L     : amu = ${amu_2l_gm2calc} +/- ${damu_2l_gm2calc}"
echo "embedded GM2Calc 2L     : amu = ${amu_2l_gm2calc_fs} +/- ${damu_2l_gm2calc_fs}"

if test $errors -eq 0 ; then
    echo "Test status: OK"
else
    echo "Test status: FAIL"
fi

exit ${errors}
