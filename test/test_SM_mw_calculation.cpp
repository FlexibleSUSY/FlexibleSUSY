
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SM_mw_calculation

#include <boost/test/unit_test.hpp>

#define protected public

#include "lowe.h"
#include "SM_input_parameters.hpp"
#include "SM_two_scale_spectrum_generator.hpp"
#include "ew_input.hpp"

using namespace flexiblesusy;

SM<Two_scale> run_SM(const SM_input_parameters& input)
{
   softsusy::QedQcd qedqcd;
   qedqcd.to(qedqcd.displayPoleMZ());

   Spectrum_generator_settings settings;
   settings.set(Spectrum_generator_settings::precision, 1.0e-5);

   SM_spectrum_generator<Two_scale> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run_except(qedqcd, input);

   return spectrum_generator.get_model();
}

BOOST_AUTO_TEST_CASE( test_consistency )
{
   SM_input_parameters input;
   input.LambdaIN = 0.13;
   input.Qin      = 91.1876;
   input.QEWSB    = 173.34;

   double mw;
   BOOST_REQUIRE_NO_THROW(mw = run_SM(input).get_physical().MVWp);
   const double mwSM = Electroweak_constants::MWSM;
   BOOST_CHECK_CLOSE_FRACTION(mw, mwSM, 1.0e-10);
}
