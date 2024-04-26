#!/bin/sh

BASEDIR="$(dirname $0)"
MODELDIR="${BASEDIR}/../models"

MS=2000
TB=5
Xt=0
Xb=0

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
    1   1000                 # max. iterations (0 = automatic)
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
   12   1                    # force output
   13   1                    # Top quark 2-loop corrections QCD
   14   1.000000000e-11      # beta-function zero threshold
   15   0                    # calculate all observables
${sminputs_tmpl}
Block MINPAR                 # Input parameters
    4   1                    # SignMu
"

run_sg() {
    local SG="$1"
    local slha_input="$2"
    local output="$3"
    local output_block
    output_block=$(echo "${output}" | cut -d'-' -f1)
    local output_entry
    output_entry=$(echo "${output}" | cut -d'-' -f2)
    local slha_output=
    local block=
    local value=

    # run the spectrum generator
    slha_output=$(echo "$slha_input" | $SG --slha-input-file=- 2>/dev/null)

    block=$(echo "$slha_output" | awk -v block="$output_block" "$print_slha_block_awk")
    value=$(echo "$block"       | awk -v keys="$output_entry" "$print_block_entry_awk" | tail -n 1)

    [ "x$value" = "x" ] && value="-"

    echo $value
}

run_MSSMEFTHiggs() {
    local SG="$1"
    local MS2 At Ab
    MS2=$(echo "scale=5; ${MS}^2" | bc)
    At=$(echo "scale=10; (1./${TB} + ${Xt}) * ${MS}" | bc)
    Ab=$(echo "scale=10; (${TB} + ${Xb}) * ${MS}" | bc)
    local slha_input="\
${slha_tmpl}
Block EXTPAR                 # Input parameters
    0   ${MS}                # MSUSY
    1   ${MS}                # M1(MSUSY)
    2   ${MS}                # M2(MSUSY)
    3   ${MS}                # M3(MSUSY)
    4   ${MS}                # Mu(MSUSY)
    5   ${MS}                # mA(MSUSY)
   25   ${TB}                # TanBeta(MSUSY)
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
Block AUIN
  1  1     ${At} # Au(1,1)
  2  2     ${At} # Au(2,2)
  3  3     ${At} # Au(3,3)
Block ADIN
  1  1     ${Ab} # Ad(1,1)
  2  2     ${Ab} # Ad(2,2)
  3  3     ${Ab} # Ad(3,3)
Block AEIN
  1  1     ${Ab} # Ad(1,1)
  2  2     ${Ab} # Ad(2,2)
  3  3     ${Ab} # Ad(3,3)
"

    echo "$(run_sg "$SG" "$slha_input" "MASS-25")"
}

run_MSSMNoFVEFTHiggs() {
    local SG="$1"
    local MS2 At Ab
    MS2=$(echo "scale=5; ${MS}^2" | bc)
    At=$(echo "scale=10; (1./${TB} + ${Xt}) * ${MS}" | bc)
    Ab=$(echo "scale=10; (${TB} + ${Xb}) * ${MS}" | bc)
    local slha_input="\
${slha_tmpl}
Block EXTPAR
    0   ${MS}                # input scale
    1   ${MS}                # M1
    2   ${MS}                # M2
    3   ${MS}                # M3
   11   ${At}                # At
   12   ${Ab}                # Ab
   13   ${Ab}                # Atau
   14   ${At}                # Ac
   15   ${Ab}                # As
   16   ${Ab}                # Amuon
   17   ${At}                # Au
   18   ${Ab}                # Ad
   19   ${Ab}                # Ae
   23   ${MS}                # Mu
   24   ${MS2}               # mA^2
   25   ${TB}                # TanBeta
   31   ${MS}
   32   ${MS}
   33   ${MS}
   34   ${MS}
   35   ${MS}
   36   ${MS}
   41   ${MS}
   42   ${MS}
   43   ${MS}
   44   ${MS}
   45   ${MS}
   46   ${MS}
   47   ${MS}
   48   ${MS}
   49   ${MS}
"

    echo "$(run_sg "$SG" "$slha_input" "MASS-25")"
}

scan() {
    local start="$1"
    local stop="$2"
    local steps="$3"

    printf "# %14s %16s %16s\n" "MS" "MSSMEFTHiggs" "MSSMNoFVEFTHiggs"

    for i in $(seq 0 $steps); do
        MS=$(cat <<EOF | bc -l
scale=10
e(l($start) + (l($stop) - l($start))*${i} / $steps)
EOF
          )

        MhMSSMEFTHiggs=$(run_MSSMEFTHiggs "$MODELDIR/MSSMEFTHiggs/run_MSSMEFTHiggs.x")
        MhMSSMNoFVEFTHiggs=$(run_MSSMNoFVEFTHiggs "$MODELDIR/MSSMNoFVEFTHiggs/run_MSSMNoFVEFTHiggs.x")

        CHECK_EQUAL_FRACTION "$MhMSSMEFTHiggs" "$MhMSSMNoFVEFTHiggs" "0.0005" || error=$(expr $error + 1)

        printf "%16s %16s %16s\n" "$MS" "$MhMSSMEFTHiggs" "$MhMSSMNoFVEFTHiggs"
    done
}

start=91.1876
stop=100000
steps=20

error=0

scan "$start" "$stop" "$steps"

if [ "x$error" != "x0" ] ; then
    echo "Test FAILED: There were $error errors."
else
    echo "All tests passed."
fi

exit $error
