#!/bin/sh

command -v wolframscript > /dev/null 2>&1 || {
    printf "%s\n" "Wolframscript not found."
    exit
}

BASEDIR=$(dirname $0)
ws="${BASEDIR}/../utils/mathws"
errors=0

### test exit code ###

cat <<EOF | ${ws} > /dev/null 2>&1
Print[1]; Quit[1]
EOF

[ "x$?" != "x1" ] && errors=$(expr $errors + 1)

### test summary ###

if [ "x$errors" = "x0" ] ; then
    printf "%s\n" "All tests passed."
else
    printf "%s\n" "There were ${errors} error(s)!"
fi

exit ${errors}
