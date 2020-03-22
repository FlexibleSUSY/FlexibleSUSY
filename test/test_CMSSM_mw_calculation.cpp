
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSM_mw_calculation

#include <boost/test/unit_test.hpp>

#define protected public

#include "lowe.h"
#include "CMSSM_input_parameters.hpp"
#include "CMSSM_two_scale_spectrum_generator.hpp"
#include "ew_input.hpp"

using namespace flexiblesusy;

CMSSM<Two_scale> run_CMSSM(const CMSSM_input_parameters& input)
{
   softsusy::QedQcd qedqcd;
   qedqcd.to(qedqcd.displayPoleMZ());

   Spectrum_generator_settings settings;
   settings.set(Spectrum_generator_settings::precision, 1.0e-5);

   CMSSM_spectrum_generator<Two_scale> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run_except(qedqcd, input);

   return spectrum_generator.get_model();
}

BOOST_AUTO_TEST_CASE( test_decoupling )
{
   double mw1, mw2, mw5, mw10;

   CMSSM_input_parameters input;
   input.TanBeta = 10.;
   input.SignMu  = 1;
   input.Azero   = 0.;

   input.m0  = 1000.;
   input.m12 = 1000.;
   BOOST_REQUIRE_NO_THROW(mw1 = run_CMSSM(input).get_physical().MVWm);

   input.m0  = 2000.;
   input.m12 = 2000.;
   BOOST_REQUIRE_NO_THROW(mw2 = run_CMSSM(input).get_physical().MVWm);

   input.m0  = 5000.;
   input.m12 = 5000.;
   BOOST_REQUIRE_NO_THROW(mw5 = run_CMSSM(input).get_physical().MVWm);

   input.m0  = 10000.;
   input.m12 = 10000.;
   BOOST_REQUIRE_NO_THROW(mw10 = run_CMSSM(input).get_physical().MVWm);

   const double mwSM = Electroweak_constants::MWSM;

   BOOST_CHECK_GT(std::abs(mw1 / mwSM - 1.), 1.0e-5);
   BOOST_CHECK_CLOSE_FRACTION(mw1, mw2, 1.0e-4);
   BOOST_CHECK_CLOSE_FRACTION(mw2, mw5, 1.0e-5);
   BOOST_CHECK_CLOSE_FRACTION(mw5, mw10, 1.0e-6);
   BOOST_CHECK_CLOSE_FRACTION(mw10, mwSM, 5.0e-7);
}
