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

# TP3 from arxiv:1507.05093 where SM-like Higgs is 2nd lightest
input_lowNMSSM="
Block FlexibleSUSY
    0   1.000000000e-04      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = two_scale, 1 = lattice)
    3   1                    # calculate SM pole masses
    4   2                    # pole mass loop order
    5   2                    # EWSB loop order
    6   2                    # beta-functions loop order
    7   2                    # threshold corrections loop order
    8   1                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   1                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   1                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   1                    # Higgs 2-loop corrections O(alpha_tau^2)
   13   1                    # Top pole mass QCD corrections (0 = 1L, 1 = 2L, 2 = 3L)
   14   1.000000000e-11      # beta-function zero threshold
   15   0                    # calculate all observables
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
BLOCK SMINPUTS	 	 #
1	128.962	 	 #
2	0.000011663900000000002	 	 #
3	0.1172	 	 #
4	91.1876	 	 #
5	4.2	 	 #
6	172.9	 	 #
7	1.777	 	 #
9	80.385	 	 #
11	0.00051099891	 	 #
13	1.10565836	 	 #
21	0.00495	 	 #
22	0.0025	 	 #
23	0.1	 	 #
24	1.42	 	 #
BLOCK EXTPAR	 	 #
0	1000.	 	 #
1	200.	 	 #
2	400.	 	 #
3	2000.	 	 #
11	1000.	 	 #
12	1000.	 	 #
13	1000.	 	 #
25	3.	 	 #
31	1500.	 	 #
32	1500.	 	 #
33	1500.	 	 #
34	1500.	 	 #
35	1500.	 	 #
36	1500.	 	 #
41	1500.	 	 #
42	1500.	 	 #
43	1000.	 	 #
44	1500.	 	 #
45	1500.	 	 #
46	1000.	 	 #
47	1500.	 	 #
48	1500.	 	 #
49	1500.	 	 #
61	0.67	 	 #
62	0.1	 	 #
63	650.	 	 #
64	-10.	 	 #
65	200.	 	 #
"

input_NMSSMEFTHiggsTwoScale="
Block FlexibleSUSY
    0   1.000000000e-04      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = two_scale, 1 = lattice)
    3   1                    # calculate SM pole masses
    4   2                    # pole mass loop order
    5   2                    # EWSB loop order
    6   2                    # beta-functions loop order
    7   2                    # threshold corrections loop order
    8   1                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   1                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   1                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   1                    # Higgs 2-loop corrections O(alpha_tau^2)
   13   1                    # Top pole mass QCD corrections (0 = 1L, 1 = 2L, 2 = 3L)
   14   1.000000000e-11      # beta-function zero threshold
   15   0                    # calculate all observables
   16   0                    # force positive majorana masses
   17   0                    # pole mass renormalization scale (0 = SUSY scale)
   18   0                    # pole mass renormalization scale in the EFT (0 = min(SUSY scale, Mt))
   19   0                    # EFT matching scale (0 = SUSY scale)
   20   2                    # EFT loop order for upwards matching
   21   1                    # EFT loop order for downwards matching
   22   1                    # EFT index of SM-like Higgs in the BSM model
BLOCK SMINPUTS
    1   128.962
    2   0.000011663900000000002
    3   0.1172
    4   91.1876
    5   4.2
    6   172.9
    7   1.777
    9   80.385
   11   0.00051099891
   13   1.10565836
   21   0.00495
   22   0.0025
   23   0.1
   24   1.42
Block MINPAR                 # Input parameters
    4   1                    # SignMu
Block EXTPAR
    0   1000                 # Ms
    1   200                  # M1(MSUSY)
    2   400                  # M2(MSUSY)
    3   2000                 # M3(MSUSY)
    4   200                  # Mu(MSUSY) CHECK factor SQRT[2]
   25   3                    # tan(beta) at Ms
   61   0.67                 # Lambda
   62   0.1                  # Kappa
   63   650                  # ALambda
   64  -10                   # AKappa
Block MSQ2IN
  1  1     2250000          # mq2(1,1)
  2  2     2250000          # mq2(2,2)
  3  3     1.00000000E+06   # mq2(3,3)
Block MSE2IN
  1  1     2250000          # me2(1,1)
  2  2     2250000          # me2(2,2)
  3  3     2250000          # me2(3,3)
Block MSL2IN
  1  1     2250000          # ml2(1,1)
  2  2     2250000          # ml2(2,2)
  3  3     2250000          # ml2(3,3)
Block MSU2IN
  1  1     2250000          # mu2(1,1)
  2  2     2250000          # mu2(2,2)
  3  3     1.00000000E+06   # mu2(3,3)
Block MSD2IN
  1  1     2250000          # md2(1,1)
  2  2     2250000          # md2(2,2)
  3  3     2250000          # md2(3,3)
Block AUIN
  3  3     1000             # Au(3,3)
Block ADIN
  3  3     1000             # Ad(3,3)
"

run_sg() {
    local SG="$1"
    local slha_input="$2"
    local slha_output=
    local block=
    local value=
    local output_block=MASS
    local output_entry=35

    # run the spectrum generator
    slha_output=$(echo "$slha_input" | $SG --slha-input-file=- 2>/dev/null)

    block=$(echo "$slha_output" | awk -v block="$output_block" "$print_slha_block_awk")
    value=$(echo "$block"       | awk -v keys="$output_entry" "$print_block_entry_awk" | tail -n 1)

    [ "x$value" = "x" ] && value="-"

    echo $value
}

error=0

# run TP3
MhlowNMSSM=$(run_sg "$MODELDIR/lowNMSSM/run_lowNMSSM.x" "$input_lowNMSSM")
MhNMSSMEFTHiggsTwoScale=$(run_sg "$MODELDIR/NMSSMEFTHiggsTwoScale/run_NMSSMEFTHiggsTwoScale.x" "$input_NMSSMEFTHiggsTwoScale")

echo "Mh in the lowNMSSM  : $MhlowNMSSM"
echo "Mh in the NMSSMEFTHiggsTwoScale: $MhNMSSMEFTHiggsTwoScale"

CHECK_EQUAL_FRACTION "$MhlowNMSSM" "$MhNMSSMEFTHiggsTwoScale" "0.007" || error=$(expr $error + 1)

if [ "x$error" != "x0" ] ; then
    echo "Test FAILED: There were $error errors."
else
    echo "All tests passed."
fi

exit $error
