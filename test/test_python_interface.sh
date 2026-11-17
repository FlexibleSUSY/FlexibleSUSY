#!/bin/sh

# Runs the tests of the Python interface (python/tests).  Tests that need a
# generated model are skipped if the model isn't built.  The test is skipped
# if no python3 is found.

BASEDIR=$(dirname "$0")
PYTHON="${PYTHON:-python3}"

if ! command -v "$PYTHON" > /dev/null 2>&1 ; then
    echo "Skipping: $PYTHON not found"
    exit 0
fi

cd "${BASEDIR}/../python" || exit 1

# -I: isolated mode, so that a package named "tests" installed in the user's
# site-packages can't shadow python/tests
exec "$PYTHON" -I -m unittest discover -s tests -v
