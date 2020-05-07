#!/bin/sh

exit_code=0

body_test () {
   export FLEXIBLESUSY_LOOP_LIBRARY="$1"
   TEST_LOOPLIBRARY=$(./test_looplibrary_environment.x -l message 2>/dev/null | grep -oG 'lib<.*>')
   [ $TEST_LOOPLIBRARY = "lib<$2>" ] || exit_code=1
}

cd test

body_test '0' 'Softsusy'
body_test 'not an integer' 'Softsusy'
grep -q '#define ENABLE_COLLIER' ../config/config.h && body_test '1' 'Collier'
grep -q '#define ENABLE_LOOPTOOLS' ../config/config.h && body_test '2' 'Looptools'
grep -q '#define ENABLE_FFLITE' ../config/config.h && body_test '3' 'Fflite'

exit "$exit_code"
