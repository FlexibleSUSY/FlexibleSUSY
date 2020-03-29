#!/bin/sh

BASEDIR="$(dirname $0)"
MODELDIR="${BASEDIR}/../models"

. "$BASEDIR/test.sh"

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
   pattern     = "^block[[:blank:]]*" tolower(block) "([^[:graph:]].*)?$";
   not_pattern = "^block[[:blank:]]*.*$";

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

sminputs_tmpl="\
Block SMINPUTS               # Standard Model inputs
    1   1.279440000e+02      # alpha^(-1) SM MSbar(MZ)
    2   1.166378700e-05      # G_Fermi
    3   1.184000000e-01      # alpha_s(MZ) SM MSbar
    4   9.118760000e+01      # MZ(pole)
    5   4.180000000e+00      # mb(mb) SM MSbar
    6   1.733400000e+02      # mtop(pole)
    7   1.777000000e+00      # mtau(pole)
    8   0.000000000e+00      # mnu3(pole)
    9   80.384               # MW pole
   11   5.109989020e-04      # melectron(pole)
   12   0.000000000e+00      # mnu1(pole)
   13   1.056583570e-01      # mmuon(pole)
   14   0.000000000e+00      # mnu2(pole)
   21   4.750000000e-03      # md(2 GeV) MS-bar
   22   2.400000000e-03      # mu(2 GeV) MS-bar
   23   1.040000000e-01      # ms(2 GeV) MS-bar
   24   1.270000000e+00      # mc(mc) MS-bar
"

slha_tmpl="\
Block MODSEL                 # Select model
    6   0                    # flavour violation
Block FlexibleSUSY
    0   1.000000000e-05      # precision goal
    1   10000                # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = two_scale, 1 = lattice)
    3   1                    # calculate SM pole masses
    4   2                    # pole mass loop order
    5   2                    # EWSB loop order
    6   3                    # beta-functions loop order
    7   2                    # threshold corrections loop order
    8   1                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   1                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   1                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   1                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   0                    # force output
   13   1                    # Top quark 2-loop corrections QCD
   14   1.000000000e-11      # beta-function zero threshold
   15   0                    # calculate observables (a_muon, ...)
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
${sminputs_tmpl}
Block EXTPAR
  100   1                    # LambdaLoopOrder (HSSUSY)
  101   1                    # 2-loop at*as
  102   1                    # 2-loop ab*as
  103   1                    # 2-loop at*ab
  104   1                    # 2-loop atau*atau
  105   1                    # 2-loop at*at
"

run_sg() {
    local SG="$1"
    local eftLoopOrder="$2"
    local MS=38700000
    local MS2=$(($MS * $MS))
    local TB=1
    local Msplit=1000
    local At="$Msplit"
    local slha_input
    local slha_output
    local lambdaFull
    local lambdaEFT

    slha_input=$(
      echo "$slha_tmpl" ;
      cat <<EOF
Block MODSEL                 # Select model
   12   ${Msplit}            # DRbar parameter output scale (GeV)
Block FlexibleSUSY
   21   ${eftLoopOrder}      # EFT loop order for downwards matching
Block EXTPAR                 # Input parameters
    0   ${MS}                # MSUSY
    1   ${Msplit}            # M1(MSUSY)
    2   ${Msplit}            # M2(MSUSY)
    3   ${Msplit}            # M3(MSUSY)
    4   ${Msplit}            # Mu(MSUSY)
    5   ${MS}                # mA(MSUSY)
    6   173.34               # MEWSB
    7   ${At}                # At(MSUSY)
   25   ${TB}                # TanBeta(MSUSY)
  100   2                    # LambdaLoopOrder
Block MSQ2IN
  1  1     ${MS2}   # mq2(1,1)
  2  2     ${MS2}   # mq2(2,2)
  3  3     ${MS2}   # mq2(3,3)
Block MSE2IN
  1  1     ${MS2}   # me2(1,1)
  2  2     ${MS2}   # me2(2,2)
  3  3     ${MS2}   # me2(3,3)
Block MSL2IN
  1  1     ${MS2}   # ml2(1,1)
  2  2     ${MS2}   # ml2(2,2)
  3  3     ${MS2}   # ml2(3,3)
Block MSU2IN
  1  1     ${MS2}   # mu2(1,1)
  2  2     ${MS2}   # mu2(2,2)
  3  3     ${MS2}   # mu2(3,3)
Block MSD2IN
  1  1     ${MS2}   # md2(1,1)
  2  2     ${MS2}   # md2(2,2)
  3  3     ${MS2}   # md2(3,3)
EOF
    )

    # run the spectrum generator
    slha_output=$(echo "$slha_input" | $SG --slha-input-file=- 2>/dev/null)

    lambdaFull=$(echo "$slha_output" |
		 awk -v block="SM" "$print_slha_block_awk" |
		 awk -v keys="2" "$print_block_entry_awk" |
		 tail -n 1)
    lambdaEFT=$( echo "$slha_output" |
		 awk -v block="SMSM" "$print_slha_block_awk" |
		 awk -v keys="2" "$print_block_entry_awk" |
		 tail -n 1)
    if [ -z "$lambdaFull" -o -z "$lambdaEFT" ]; then
       lambdaFull=0
       lambdaEFT=1
    fi

    echo "$lambdaFull" "$lambdaEFT"
}

error=0

# check equality of lambdas in full theory and EFT at matching scale

set -- $(run_sg "$MODELDIR/SplitMSSMEFTHiggs/run_SplitMSSMEFTHiggs.x" 0)
CHECK_EQUAL_FRACTION "$1" "$2" "0.0001" || error=$(expr $error + 1)

# check approximate equality of lambdas in full theory and EFT at matching scale
# this would fail if difference between vev normalizations
# in full theory and SM were ignored, see commit e5473865150da98e1426f2282baf31d54541169a
set -- $(run_sg "$MODELDIR/SplitMSSMEFTHiggs/run_SplitMSSMEFTHiggs.x" 1)
CHECK_EQUAL_FRACTION "$1" "$2" "0.03" || error=$(expr $error + 1)

if [ "x$error" != "x0" ] ; then
    echo "Test FAILED: There were $error errors."
else
    echo "All tests passed."
fi

exit $error
