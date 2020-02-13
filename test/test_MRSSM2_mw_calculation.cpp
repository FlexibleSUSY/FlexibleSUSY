
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MRSSM2_mw_calculation

#include <boost/test/unit_test.hpp>

#define protected public

#include "lowe.h"
#include "MRSSM2_input_parameters.hpp"
#include "MRSSM2_two_scale_spectrum_generator.hpp"
#include "ew_input.hpp"

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

BOOST_AUTO_TEST_CASE( test_decoupling )
{
   double mw1, mw2, mw5, mw10;

   BOOST_REQUIRE_NO_THROW(mw1 = run_MRSSM2(1000.).get_physical().MVWm);
   BOOST_REQUIRE_NO_THROW(mw2 = run_MRSSM2(2000.).get_physical().MVWm);
   BOOST_REQUIRE_NO_THROW(mw5 = run_MRSSM2(5000.).get_physical().MVWm);
   BOOST_REQUIRE_NO_THROW(mw10 = run_MRSSM2(10000.).get_physical().MVWm);

   const double mwSM = Electroweak_constants::MWSM;

   BOOST_CHECK_GT(std::abs(mw1 / mwSM - 1.), 5.0e-4);
   BOOST_CHECK_CLOSE_FRACTION(mw1, mw2, 1.0e-3);
   BOOST_CHECK_CLOSE_FRACTION(mw2, mw5, 5.0e-4);
   BOOST_CHECK_CLOSE_FRACTION(mw5, mw10, 1.0e-4);
   BOOST_CHECK_CLOSE_FRACTION(mw10, mwSM, 5.0e-5);
}
