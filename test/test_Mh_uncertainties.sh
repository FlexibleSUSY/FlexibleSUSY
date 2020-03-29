#!/bin/sh

BASEDIR="$(dirname $0)"
FSCONFIG="${BASEDIR}/../flexiblesusy-config"
MATH=$(${FSCONFIG} --math-cmd)

printf "%s" "<< \"${BASEDIR}/test_Mh_uncertainties.m\"; Quit[]" | "$MATH"
error="$?"

[ "x${error}" != x0 ] && {
    echo "Test FAILED: Mathematica program exited with error code ${error}."
    exit 1
}

plot() {
    local Xt="$1"
    local term="$2"

    cat <<EOF | gnuplot
set terminal ${term}
set logscale x 10
set size square
set key box bottom right width 0
set grid
set xlabel "M_S / TeV"
set ylabel "M_h / GeV"

set style line 1 lt 1 dt 2 lw 2 lc rgb '#0000FF'
set style line 2 lt 1 dt 1 lw 2 lc rgb '#FF0000'
set style line 3 lt 1 dt 4 lw 2 lc rgb '#00CC00'
set style line 4 lt 1 dt 3 lw 2 lc rgb '#000000'

data = "${BASEDIR}/test_Mh_uncertainties_TB-5_Xt-${Xt}.dat"
set output "${BASEDIR}/test_Mh_uncertainties_TB-5_Xt-${Xt}.pdf"

plot \
     data u (\$1/1000):2 t "MSSM 2L"         w lines ls 1, \
     data u (\$1/1000):4 t "MSSM 3L"         w lines ls 2, \
     data u (\$1/1000):6 t "MSSMEFTHiggs 1L" w lines ls 3, \
     data u (\$1/1000):8 t "HSSUSY 2L"       w lines ls 4

plot \
     data u (\$1/1000):2 t "MSSM 2L"         w lines ls 1, \
     data u (\$1/1000):4 t "MSSM 3L"         w lines ls 2, \
     data u (\$1/1000):6 t "MSSMEFTHiggs 1L" w lines ls 3, \
     data u (\$1/1000):8 t "HSSUSY 2L"       w lines ls 4, \
     data u (\$1/1000):(\$2-\$3):(\$2+\$3) w filledcurves ls 1 fs transparent solid 0.3 t '', \
     data u (\$1/1000):(\$4-\$5):(\$4+\$5) w filledcurves ls 2 fs transparent solid 0.3 t '', \
     data u (\$1/1000):(\$6-\$7):(\$6+\$7) w filledcurves ls 3 fs transparent solid 0.3 t '', \
     data u (\$1/1000):(\$8-\$9):(\$8+\$9) w filledcurves ls 4 fs transparent solid 0.3 t ''

set ylabel "(M_h - M_h^{MSSMEFTHiggs}) / GeV"
set key box bottom left width 0

plot \
     data u (\$1/1000):(\$2-\$6) t "MSSM 2L"         w lines ls 1, \
     data u (\$1/1000):(\$4-\$6) t "MSSM 3L"         w lines ls 2, \
     data u (\$1/1000):(\$6-\$6) t "MSSMEFTHiggs 1L" w lines ls 3, \
     data u (\$1/1000):(\$8-\$6) t "HSSUSY 2L"       w lines ls 4
     # data u (\$1/1000):((\$2-\$3):(\$2+\$3) w filledcurves ls 1 fs transparent solid 0.3 t '', \
     # data u (\$1/1000):((\$4-\$5):(\$4+\$5) w filledcurves ls 2 fs transparent solid 0.3 t '', \
     # data u (\$1/1000):((\$6-\$7):(\$6+\$7) w filledcurves ls 3 fs transparent solid 0.3 t '', \
     # data u (\$1/1000):((\$8-\$9):(\$8+\$9) w filledcurves ls 4 fs transparent solid 0.3 t ''

set ylabel "{/Symbol D}M_h / GeV"
set key box top left width 0

plot \
     data u (\$1/1000):3 t "MSSM 2L"         w lines ls 1, \
     data u (\$1/1000):5 t "MSSM 3L"         w lines ls 2, \
     data u (\$1/1000):7 t "MSSMEFTHiggs 1L" w lines ls 3, \
     data u (\$1/1000):9 t "HSSUSY 2L"       w lines ls 4
EOF
}

for Xtt in 0 -2 ; do
    plot "${Xtt}" "pdf size 3in,3in"
done

echo "Test passed."
