#!/bin/sh

# Checks BR(mu -> 3e) in the MSSM at a point with lepton flavour
# violation in the right-handed slepton masses (me2(1,2)) and small
# tan(beta), where the Z penguins and boxes contribute about 10%.
#
# The reference value was calculated after fixing the order in which the
# scalar penguins, vector penguins and boxes are passed to the C++ code
# (the boxes were treated as penguins and vice versa, which gave
# BR(mu -> 3e) = 3.221e-12 here instead of 2.876e-12).  In a low-scale MSSM
# setup the fixed order also agrees better with the BR(mu -> 3e)/BR(mu
# -> e gamma) of SARAH's SPheno (FlavorKit).

BASEDIR=$(dirname $0)

# shellcheck source=test.sh
. "$BASEDIR/test.sh"

exe="$BASEDIR/../models/MSSM/run_MSSM.x"
input="$BASEDIR/test_MSSM_br_l_to_3l.in.spc"
output="$BASEDIR/test_MSSM_br_l_to_3l.out.spc"

if test ! -x "$exe"; then
    echo "Error: MSSM spectrum generator not found: $exe"
    exit 1
fi

"$exe" --slha-input-file="$input" --slha-output-file="$output" > /dev/null 2>&1 || {
    echo "Error: $exe failed"
    exit 1
}

# BR(mu -> 3e) in Block FlexibleSUSYLowEnergy, entry 31
br_mu_3e=$(awk '/^Block FlexibleSUSYLowEnergy/ { b = 1; next }
                /^Block/ { b = 0 }
                b && $1 == 31 { print $2 }' "$output")

echo "BR(mu -> 3e) = $br_mu_3e"

CHECK_EQUAL_FRACTION "$br_mu_3e" "2.87613447E-12" "1e-5" || {
    echo "Test status: FAIL"
    exit 1
}

echo "Test status: OK"
