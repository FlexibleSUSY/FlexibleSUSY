#!/bin/sh

exit_code=0

exe_test () {
   TEST_LOOPLIBRARY=$(./test_looplibrary_environment.x \
                      -l message 2>/dev/null | grep -oP 'lib<(.*)>')
}

cd test

export FLEXIBLESUSY_LOOP_LIBRARY='0'
exe_test
[ $TEST_LOOPLIBRARY = 'lib<Softsusy>' ] || exit_code=1

export FLEXIBLESUSY_LOOP_LIBRARY="not an integer"
exe_test
[ $TEST_LOOPLIBRARY = 'lib<Softsusy>' ] || exit_code=1

if grep -q '#define ENABLE_COLLIER' ../config/config.h; then
   export FLEXIBLESUSY_LOOP_LIBRARY='1'
   exe_test
   [ $TEST_LOOPLIBRARY = 'lib<Collier>' ] || exit_code=1
fi

if grep -q '#define ENABLE_LOOPTOOLS' ../config/config.h; then
   export FLEXIBLESUSY_LOOP_LIBRARY='2'
   exe_test
   [ $TEST_LOOPLIBRARY = 'lib<Looptools>' ] || exit_code=1
fi

if grep -q '#define ENABLE_FFLITE' ../config/config.h; then
   export FLEXIBLESUSY_LOOP_LIBRARY='3'
   exe_test
   [ $TEST_LOOPLIBRARY = 'lib<Fflite>' ] || exit_code=1
fi

exit "$exit_code"
