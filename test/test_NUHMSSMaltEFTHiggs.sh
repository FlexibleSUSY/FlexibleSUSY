#!/bin/sh

BASEDIR=$(dirname $0)
MODELDIR=${BASEDIR}/../models

# shellcheck source=test.sh
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

input="
Block MODSEL                 # Select model
#   12    1000                # DRbar parameter output scale (GeV)
Block FlexibleSUSY
    0   1.000000000e-04      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = two_scale, 1 = lattice)
    3   0                    # calculate SM pole masses
    4   2                    # pole mass loop order
    5   2                    # EWSB loop order
    6   2                    # beta-functions loop order
    7   2                    # threshold corrections loop order
    8   1                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   1                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   1                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   1                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   0                    # force output
   13   1                    # Top quark 2-loop corrections QCD
   14   1.000000000e-11      # beta-function zero threshold
   15   0                    # calculate all observables
Block SMINPUTS               # Standard Model inputs
    1   1.279340000e+02      # alpha^(-1) SM MSbar(MZ)
    2   1.166370000e-05      # G_Fermi
    3   1.176000000e-01      # alpha_s(MZ) SM MSbar
    4   9.118760000e+01      # MZ(pole)
    5   4.200000000e+00      # mb(mb) SM MSbar
    6   1.733000000e+02      # mtop(pole)
    7   1.777000000e+00      # mtau(pole)
    8   0.000000000e+00      # mnu3(pole)
    9   80.404               # MW pole
   11   5.109989020e-04      # melectron(pole)
   12   0.000000000e+00      # mnu1(pole)
   13   1.056583570e-01      # mmuon(pole)
   14   0.000000000e+00      # mnu2(pole)
   21   4.750000000e-03      # md(2 GeV) MS-bar
   22   2.400000000e-03      # mu(2 GeV) MS-bar
   23   1.040000000e-01      # ms(2 GeV) MS-bar
   24   1.270000000e+00      # mc(mc) MS-bar
Block MINPAR                 # Input parameters
    1   91                   # m0
    2   91                   # m12
    3   10                   # TanBeta
    4   1                    # SignMu
    5   0                    # Azero
    6   91                   # MuInput
    7   91                   # BInput
"

run_sg() {
    local SG="$1"
    local slha_input="$2"
    local slha_output=
    local block=
    local value=
    local output_block=MASS
    local output_entry=25

    # run the spectrum generator
    slha_output=$(echo "$slha_input" | $SG --slha-input-file=- 2>/dev/null)

    block=$(echo "$slha_output" | awk -v block="$output_block" "$print_slha_block_awk")
    value=$(echo "$block"       | awk -v keys="$output_entry" "$print_block_entry_awk" | tail -n 1)

    [ "x$value" = "x" ] && value="-"

    echo $value
}

error=0

MhNUHMSSMalt=$(run_sg "$MODELDIR/NUHMSSMalt/run_NUHMSSMalt.x" "$input")
MhNUHMSSMaltEFTHiggs=$(run_sg "$MODELDIR/NUHMSSMaltEFTHiggs/run_NUHMSSMaltEFTHiggs.x" "$input")

echo "Mh in the NUHMSSMalt     : $MhNUHMSSMalt"
echo "Mh in the NUHMSSMaltEFTHiggs: $MhNUHMSSMaltEFTHiggs"

CHECK_EQUAL_FRACTION "$MhNUHMSSMalt" "$MhNUHMSSMaltEFTHiggs" "0.008" || error=$(expr $error + 1)

if [ "x$error" != "x0" ] ; then
    echo "Test FAILED: There were $error errors."
else
    echo "All tests passed."
fi

exit $error
