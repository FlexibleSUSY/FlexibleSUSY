
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSM_mw_calculation

#include <boost/test/unit_test.hpp>

#define protected public

#include "lowe.h"
#include "CMSSM_input_parameters.hpp"
#include "CMSSM_two_scale_spectrum_generator.hpp"
#include "ew_input.hpp"
#include "sm_mw.hpp"
#include <utility>
#include <tuple>

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

std::pair<double, double> calc_mw_mh_CMSSM(const CMSSM_input_parameters& input)
{
   const auto model = run_CMSSM(input);
   return std::make_pair(model.get_physical().MVWm, model.get_physical().Mhh(0));
}

double calc_mw_SM(double mh)
{
   using flexiblesusy::sm_mw::calculate_mw_pole_SM_fit_MSbar;
   softsusy::QedQcd qedqcd;

   const auto res = calculate_mw_pole_SM_fit_MSbar(
      mh,
      qedqcd.displayPoleMt(),
      qedqcd.displayAlphaSInput(),
      Electroweak_constants::delta_alpha_s_5_had);

   return res.first;
}

BOOST_AUTO_TEST_CASE( test_decoupling )
{
   double mw1, mw2, mw5, mw10;
   double mh1, mh2, mh5, mh10;

   CMSSM_input_parameters input;
   input.TanBeta = 10.;
   input.SignMu  = 1;
   input.Azero   = 0.;

   input.m0  = 1000.;
   input.m12 = 1000.;
   std::tie(mw1, mh1) = calc_mw_mh_CMSSM(input);
   BOOST_CHECK_CLOSE_FRACTION(mw1, calc_mw_SM(mh1), 1.0e-4);
   BOOST_CHECK_GT(std::abs(mw1/calc_mw_SM(mh1) - 1), 1.0e-5);

   input.m0  = 2000.;
   input.m12 = 2000.;
   std::tie(mw2, mh2) = calc_mw_mh_CMSSM(input);
   BOOST_CHECK_CLOSE_FRACTION(mw2, calc_mw_SM(mh2), 1.0e-4);

   input.m0  = 5000.;
   input.m12 = 5000.;
   std::tie(mw5, mh5) = calc_mw_mh_CMSSM(input);
   BOOST_CHECK_CLOSE_FRACTION(mw5, calc_mw_SM(mh5), 1.0e-5);

   input.m0  = 10000.;
   input.m12 = 10000.;
   std::tie(mw10, mh10) = calc_mw_mh_CMSSM(input);
   BOOST_CHECK_CLOSE_FRACTION(mw10, calc_mw_SM(mh10), 1.0e-6);
}
