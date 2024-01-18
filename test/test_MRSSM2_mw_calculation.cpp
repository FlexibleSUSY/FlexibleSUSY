
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MRSSM2_mw_calculation

#include <boost/test/unit_test.hpp>

#define protected public

#include "lowe.h"
#include "MRSSM2_input_parameters.hpp"
#include "MRSSM2_two_scale_spectrum_generator.hpp"
#include "ew_input.hpp"
#include "sm_mw.hpp"

using namespace flexiblesusy;

MRSSM2<Two_scale> run_MRSSM2(double MS)
{
   softsusy::QedQcd qedqcd;
   qedqcd.to(qedqcd.displayPoleMZ());

   Spectrum_generator_settings settings;
   settings.set(Spectrum_generator_settings::precision, 1.0e-5);

   MRSSM2_input_parameters input;
   input.TanBeta    = 10.;
   input.Ms         = MS;
   input.LamSDInput = 0.0001;
   input.LamSUInput = 0.0001;
   input.LamTDInput = 0.0001;
   input.LamTUInput = -1.2;
   input.MuDInput   = MS;
   input.MuUInput   = MS;
   input.BMuInput   = Sqr(MS);
   input.mq2Input   = Eigen::DiagonalMatrix<double,3>(Sqr(MS), Sqr(MS), Sqr(MS));
   input.ml2Input   = Eigen::DiagonalMatrix<double,3>(Sqr(MS), Sqr(MS), Sqr(MS));
   input.md2Input   = Eigen::DiagonalMatrix<double,3>(Sqr(MS), Sqr(MS), Sqr(MS));
   input.mu2Input   = Eigen::DiagonalMatrix<double,3>(Sqr(MS), Sqr(MS), Sqr(MS));
   input.me2Input   = Eigen::DiagonalMatrix<double,3>(Sqr(MS), Sqr(MS), Sqr(MS));
   input.mS2Input   = Sqr(MS);
   input.mT2Input   = Sqr(MS);
   input.moc2Input  = Sqr(MS);
   input.mRd2Input  = Sqr(MS);
   input.mRu2Input  = Sqr(MS);
   input.MDBSInput  = MS;
   input.MDWBTInput = MS;
   input.MDGocInput = MS;

   MRSSM2_spectrum_generator<Two_scale> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run_except(qedqcd, input);

   return spectrum_generator.get_model();
}

std::pair<double, double> calc_mw_mh_MRSSM2(double MS)
{
   const auto model = run_MRSSM2(MS);
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
   const auto [mw1, mh1] = calc_mw_mh_MRSSM2(1000);
   const auto [mw2, mh2] = calc_mw_mh_MRSSM2(2000);
   const auto [mw5, mh5] = calc_mw_mh_MRSSM2(5000);
   // this test fails for mh >= 9999.8 with softsusy loop library
   // it works fine with Collier or LoopTools
   const auto [mw10, mh10] = calc_mw_mh_MRSSM2(9999);

   BOOST_CHECK_GT(std::abs(mw1/calc_mw_SM(mh1) - 1), 5.0e-4);

   BOOST_CHECK_CLOSE_FRACTION(mw1,  calc_mw_SM(mh1),  7e-4);
   BOOST_CHECK_CLOSE_FRACTION(mw2,  calc_mw_SM(mh2),  2e-4);
   BOOST_CHECK_CLOSE_FRACTION(mw5,  calc_mw_SM(mh5),  3e-5);
   BOOST_CHECK_CLOSE_FRACTION(mw10, calc_mw_SM(mh10), 2e-5);
}
