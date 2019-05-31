#!/bin/sh

BASEDIR="$(dirname $0)"
MODELDIR="${BASEDIR}/../models"

start=173.34
stop=100000
steps=60

TB=5
Xt=0

output="MASS-25"
scan_data="$BASEDIR/test_MSSMEFTHiggs.dat"

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
Block MINPAR                 # Input parameters
    4   1                    # SignMu
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
    local type="$2"
    local calcBSM=$(test "x$type" = xtower ; echo $?)
    local loops=$(if [ "x$type" = xtower ] ; then echo 1; else echo 2; fi)
    # set loops=2 to have perfect agreement at large MS
    # (but this spoils agreement at MS < Mt)
    loops=2
    local MS2=$(echo "scale=5; ${MS}^2" | bc)
    local At=$(echo "scale=10; (1./${TB} + ${Xt}) * ${MS}" | bc)
    local Au=$(echo "scale=10; (1./${TB} + 0) * ${MS}" | bc)
    local Ad=$(echo "scale=10; (${TB} + 0) * ${MS}" | bc)
    local Ae="$Ad"
    local slha_output=
    local block=
    local value=
    local slha_input=
    local output_block=$(echo "${output}" | cut -d'-' -f1)
    local output_entry=$(echo "${output}" | cut -d'-' -f2)

    slha_input=$(
    { echo "$slha_tmpl" ; \
      cat <<EOF
Block EXTPAR                 # Input parameters
    0   ${MS}                # MSUSY
    1   ${MS}                # M1(MSUSY)
    2   ${MS}                # M2(MSUSY)
    3   ${MS}                # M3(MSUSY)
    4   ${MS}                # Mu(MSUSY)
    5   ${MS}                # mA(MSUSY)
    6   173.34               # MEWSB
    7   ${At}                # At(MSUSY)
   14   ${Xt}                # Xt / Ms
   25   ${TB}                # TanBeta(MSUSY)
Block FlexibleSUSY
    4   ${loops}             # pole mass loop order
    5   ${loops}             # EWSB loop order
   23   ${calcBSM}           # calculate BSM pole masses
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
  1  1     ${Au} # Au(1,1)
  2  2     ${Au} # Au(2,2)
  3  3     ${At} # Au(3,3)
Block ADIN
  1  1     ${Ad} # Ad(1,1)
  2  2     ${Ad} # Ad(2,2)
  3  3     ${Ad} # Ad(3,3)
Block AEIN
  1  1     ${Ae} # Ad(1,1)
  2  2     ${Ae} # Ad(2,2)
  3  3     ${Ae} # Ad(3,3)
EOF
    })

    # run the spectrum generator
    slha_output=$(echo "$slha_input" | $SG --slha-input-file=- 2>/dev/null)

    block=$(echo "$slha_output" | awk -v block="$output_block" "$print_slha_block_awk")
    value=$(echo "$block"       | awk -v keys="$output_entry" "$print_block_entry_awk" | tail -n 1)

    [ "x$value" = "x" ] && value="-"

    echo $value
}

scan() {
    local start="$1"
    local stop="$2"
    local steps="$3"

    printf "# %14s %16s %16s %16s\n" "MS" "MSSMEFTHiggs" "MSSMMuBMu" "HSSUSY"

    for i in $(seq 0 $steps); do
        MS=$(cat <<EOF | bc -l
scale=10
e(l($start) + (l($stop) - l($start))*${i} / $steps)
EOF
          )

        MhMSSMEFTHiggs=$(run_sg "$MODELDIR/MSSMEFTHiggs/run_MSSMEFTHiggs.x" tower)
        MhMSSMMuBMu=$(run_sg "$MODELDIR/MSSMMuBMu/run_MSSMMuBMu.x" fixedOrder)
        MhHSSUSY=$(run_sg "$MODELDIR/HSSUSY/run_HSSUSY.x" fixedOrder)

        printf "%16s %16s %16s %16s\n" "$MS" "$MhMSSMEFTHiggs" "$MhMSSMMuBMu" "$MhHSSUSY"
    done
}

scan "$start" "$stop" "$steps" | tee "$scan_data"

cat <<EOF | gnuplot
set terminal pdf enhanced size 3in,3in
set output "$BASEDIR/test_MSSMEFTHiggs.pdf"
set size square
set logscale x
set key box bottom right width 2
set grid
set xlabel "M_S / TeV"
set ylabel "M_h / GeV"

data = "$scan_data"

plot [0.1:] \
     data u (\$1/1000):2 t "MSSMEFTHiggs" w lines dt 1 lw 2 lc rgb '#FF0000', \
     data u (\$1/1000):3 t "MSSMMuBMu" w lines dt 4 lw 2 lc rgb '#00FF00', \
     data u (\$1/1000):4 t "HSSUSY"    w lines dt 2 lw 2 lc rgb '#0000FF'
EOF

error=0

# check equality of MSSMEFTHiggs and MSSMuBMu for low MS
MS=173.34
MhMSSMEFTHiggs=$(run_sg "$MODELDIR/MSSMEFTHiggs/run_MSSMEFTHiggs.x" tower)
MhMSSMMuBMu=$(run_sg "$MODELDIR/MSSMMuBMu/run_MSSMMuBMu.x" fixedOrder)
CHECK_EQUAL_FRACTION "$MhMSSMEFTHiggs" "$MhMSSMMuBMu" "0.003" || error=$(expr $error + 1)

MS=250.0
MhMSSMEFTHiggs=$(run_sg "$MODELDIR/MSSMEFTHiggs/run_MSSMEFTHiggs.x" tower)
MhMSSMMuBMu=$(run_sg "$MODELDIR/MSSMMuBMu/run_MSSMMuBMu.x" fixedOrder)
CHECK_EQUAL_FRACTION "$MhMSSMEFTHiggs" "$MhMSSMMuBMu" "0.01" || error=$(expr $error + 1)

# check equality of MSSMEFTHiggs and HSSUSY for high MS

MS=1000
MhMSSMEFTHiggs=$(run_sg "$MODELDIR/MSSMEFTHiggs/run_MSSMEFTHiggs.x" tower)
MhHSSUSY=$(run_sg "$MODELDIR/HSSUSY/run_HSSUSY.x" fixedOrder)
CHECK_EQUAL_FRACTION "$MhMSSMEFTHiggs" "$MhHSSUSY" "0.001" || error=$(expr $error + 1)

MS=10000
MhMSSMEFTHiggs=$(run_sg "$MODELDIR/MSSMEFTHiggs/run_MSSMEFTHiggs.x" tower)
MhHSSUSY=$(run_sg "$MODELDIR/HSSUSY/run_HSSUSY.x" fixedOrder)
# shift -0.05 GeV due to 2-loop O(at*as + at^2) corrections to mt(MZ)
MhHSSUSY="($MhHSSUSY - 0.05)"
CHECK_EQUAL_FRACTION "$MhMSSMEFTHiggs" "$MhHSSUSY" "0.0001" || error=$(expr $error + 1)

if [ "x$error" != "x0" ] ; then
    echo "Test FAILED: There were $error errors."
else
    echo "All tests passed."
fi

rm -f "$scan_data"

exit $error
