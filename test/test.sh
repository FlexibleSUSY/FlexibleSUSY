# This file contains helper functions for test scripts

# compares two floating point numbers for equality
# with a given maximum relative deviation
#
# Note: scientific notation is allowed
CHECK_EQUAL_FRACTION() {
    if test $# -lt 3 ; then
        echo "Error: CHECK_EQUAL_FRACTION: Too few arguments"
        echo "Usage: CHECK_EQUAL_FRACTION $num1 $num2 $fraction"
        exit 1
    fi

    [ -z "$1" -o "x$1" = "x-" ] && {
        echo "Error: CHECK_EQUAL_FRACTION: first argument is not a number: $1"
        exit 1
    }

    [ -z "$2" -o "x$2" = "x-" ] && {
        echo "Error: CHECK_EQUAL_FRACTION: second argument is not a number: $2"
        exit 1
    }

    [ -z "$3" -o "x$3" = "x-" ] && {
        echo "Error: CHECK_EQUAL_FRACTION: third argument is not a number: $3"
        exit 1
    }

    local num1="$(echo "$1" | sed -e 's/[eE]+*/*10^/')"
    local num2="$(echo "$2" | sed -e 's/[eE]+*/*10^/')"
    local frac="$(echo "$3" | sed -e 's/[eE]+*/*10^/')"

    local scale=15

    local error=$(cat <<EOF | bc
define abs(i) {
    if (i < 0) return (-i)
    return (i)
}

define min(i,j) {
    if (i < j) return i
    return j
}

define max(i,j) {
    if (i > j) return i
    return j
}

# precision of calculation
scale=${scale}

mmin=min($num1,$num2)
mmax=max($num1,$num2)
amax=max(abs($num1),abs($num2))

(mmax - mmin) > $frac * amax
EOF
    )

    if test "x$error" != "x0" ; then
        echo "Test failed: $num1 =r= $num2 with fraction $frac"
    fi

    return $error
}

# returns minimum of two numbers
min() {
    local a=$(echo "$1" | sed -e 's/[eE]+*/*10^/')
    local b=$(echo "$2" | sed -e 's/[eE]+*/*10^/')
    cat <<EOF | bc -l
define min(i,j) {
    if (i < j) return i
    return j
}

min($a,$b)
EOF
}

# returns maximum of two numbers
max() {
    local a=$(echo "$1" | sed -e 's/[eE]+*/*10^/')
    local b=$(echo "$2" | sed -e 's/[eE]+*/*10^/')
    cat <<EOF | bc -l
define max(i,j) {
    if (i > j) return i
    return j
}

max($a,$b)
EOF
}
