#!/bin/bash

BASEDIR=$(dirname $0)

input="$BASEDIR/test_HSSUSY_SUSYHD.dat"
output="$BASEDIR/test_HSSUSY_SUSYHD.out.dat"
exe="$BASEDIR/../models/HSSUSY/run_HSSUSY.x"
print_block="$BASEDIR/../utils/print_slha_block.awk"

rel_error="0.00065"

if test ! -x "$exe"; then
    echo "Error: HSSUSY spectrum generator not found: $exe"
    exit 1
fi

rm -f "$output"

# to make seq use point
LANG=en_US

Xt_values=$(seq -3 0.1 3)
MS=2000
MS2=$(echo "$MS*$MS" | bc)
Mu="$MS"
TB=10
Xb=0
Xtau=0

printf "Comparison FlexibleSUSY/HSSUSY and SUSYHD\n"
printf "MS = ${MS}, Mu = ${Mu}, MSf^2 = ${MS2}, tan(beta) = ${TB}\n\n"

printf "%8s\t%16s\n" "xt/MSUSY" "MH / GeV"

for Xt in ${Xt_values}
do
    At=$(cat <<EOF | bc
${Xt}*${MS} + ${Mu}/${TB}
EOF
      )
    Ab=$(cat <<EOF | bc
${Xb}*${MS} + ${Mu}*${TB}
EOF
      )
    Atau=$(cat <<EOF | bc
${Xtau}*${MS} + ${Mu}*${TB}
EOF
      )

    MH=$({ cat <<EOF
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
   10   1                    # Higgs 2-loop corrections O(alpha_t^2 + alpha_t alpha_b + alpha_b^2)
   11   1                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   0                    # force output
   13   1                    # Top pole mass QCD corrections (0 = 1L, 1 = 2L, 2 = 3L)
Block SMINPUTS               # Standard Model inputs
    1   1.279440000e+02      # alpha^(-1) SM MSbar(MZ)
    2   1.166380000e-05      # G_Fermi
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
Block EXTPAR                 # Input parameters
    0   ${MS}                # MSUSY
    1   ${MS}                # M1(MSUSY)
    2   ${MS}                # M2(MSUSY)
    3   ${MS}                # M3(MSUSY)
    4   ${Mu}                # Mu(MSUSY)
    5   ${MS}                # mA(MSUSY)
    6   173.34               # MEWSB
    7   ${At}                # At(MSUSY)
    8   ${Ab}                # Ab(MSUSY)
    9   ${Atau}              # Atau(MSUSY)
   25   ${TB}                # TanBeta(MSUSY)
  100   2                    # LambdaLoopOrder
  101   1                    # 2-loop at*as
  102   1                    # 2-loop ab*as
  103   1                    # 2-loop at*ab
  104   1                    # 2-loop atau*atau
  105   1                    # 2-loop at*at
Block MSQ2IN
  1  1     ${MS2}            # mq2(1,1)
  2  2     ${MS2}            # mq2(2,2)
  3  3     ${MS2}            # mq2(3,3)
Block MSE2IN
  1  1     ${MS2}            # me2(1,1)
  2  2     ${MS2}            # me2(2,2)
  3  3     ${MS2}            # me2(3,3)
Block MSL2IN
  1  1     ${MS2}            # ml2(1,1)
  2  2     ${MS2}            # ml2(2,2)
  3  3     ${MS2}            # ml2(3,3)
Block MSU2IN
  1  1     ${MS2}            # mu2(1,1)
  2  2     ${MS2}            # mu2(2,2)
  3  3     ${MS2}            # mu2(3,3)
Block MSD2IN
  1  1     ${MS2}            # md2(1,1)
  2  2     ${MS2}            # md2(2,2)
  3  3     ${MS2}            # md2(3,3)
EOF
    } | $exe --slha-input-file=- 2>/dev/null | \
        awk -f "$print_block" -v block=MASS | \
        awk '{ if ($1 == 25) print $2 - 0.1 }')

    # shift -0.1 GeV due to 2-loop O(at*as + at^2) corrections to mt(MZ)

    printf "%8s\t%16s\n" ${Xt} ${MH}
    printf "%8s\t%16s\n" ${Xt} ${MH} >> "$output"
done

echo ""
echo "Testing for maximum relative deviation < $rel_error ..."

# remove comments from input file
awk '{ if ($1 != "#") print $0 }' "$input" > "$input.tmp"

diff=$(numdiff --relative-tolerance=$rel_error "$input.tmp" "$output")

diff_without_comments=`echo $diff | sed -e '/^ *#/d' | sed -e '/^+++/d'`

exit_code=0

if [ -n "$diff_without_comments" ]; then
    echo "Error: difference between SUSYHD ($input) and FlexibleSUSY"
    echo "$diff"
    echo ""
    echo "Test result: FAIL"
    exit_code=1
else
    echo "$diff"
    echo ""
    echo "Test result: OK"
fi

rm -f "$output" "$input.tmp"

exit $exit_code
