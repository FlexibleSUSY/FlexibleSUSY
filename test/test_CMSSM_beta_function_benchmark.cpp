
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSM_beta_function_benchmark

#include <boost/test/unit_test.hpp>

#include "test_CMSSM.hpp"
#include "stopwatch.hpp"
#include "softsusy.h"
#include "CMSSM_two_scale_model.hpp"

using namespace flexiblesusy;
using namespace softsusy;

BOOST_AUTO_TEST_CASE( test_CMSSM_beta_function_calculation_time )
{
   CMSSM_input_parameters input;
   CMSSM<Two_scale> m;
   MssmSoftsusy s;
   setup_CMSSM(m, s, input);

   const int N_calls = 10000;

   Stopwatch stopwatch;
   stopwatch.start();
   for (int i = 0; i < N_calls; i++) {
      volatile SoftParsMssm beta(s.beta2());
   }
   stopwatch.stop();
   const double ss_time = stopwatch.get_time_in_seconds();

   stopwatch.start();
   for (int i = 0; i < N_calls; i++) {
      volatile CMSSM_soft_parameters beta(m.calc_beta());
   }
   stopwatch.stop();
   const double fs_time = stopwatch.get_time_in_seconds();

   BOOST_TEST_MESSAGE("Calcultating the CMSSM beta-functions " << N_calls
                 << " times with\n"
                 "Softsusy    : " << ss_time << "s\n"
                 "FlexibleSUSY: " << fs_time << "s\n");

   BOOST_CHECK_GT(ss_time, fs_time);
}
