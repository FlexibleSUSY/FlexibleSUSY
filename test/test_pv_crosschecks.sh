#!/bin/sh

exit_code=0

cd test || exit 1

./test_pv_fflite.x | ./test_pv_looptools.x --log_level=test_suite ||
    exit_code=1

./test_pv_fflite.x | ./test_pv_softsusy.x  --log_level=test_suite ||
    exit_code=1

exit "$exit_code"
