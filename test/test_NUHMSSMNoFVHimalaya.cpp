#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_NUHMSSMNoFVHimalaya

#include <boost/test/unit_test.hpp>

#include "models/NUHMSSMNoFVHimalaya/NUHMSSMNoFVHimalaya_two_scale_spectrum_generator.hpp"

using namespace flexiblesusy;

NUHMSSMNoFVHimalaya<Two_scale> run(const NUHMSSMNoFVHimalaya_input_parameters& input)
{
   softsusy::QedQcd qedqcd;
   qedqcd.to(qedqcd.displayPoleMZ());

   Spectrum_generator_settings settings;
   settings.set(Spectrum_generator_settings::precision, 1.0e-5);

   NUHMSSMNoFVHimalaya_spectrum_generator<Two_scale> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);

   return spectrum_generator.get_model();
}

double calc_Mh(int loops, double tb, double MS, double xt)
{
   NUHMSSMNoFVHimalaya_input_parameters input;

   return run(input).get_physical().Mhh(0);
}

BOOST_AUTO_TEST_CASE( test_Mh )
{
   const double tb = 20;
   const double xt = -std::sqrt(6.0);
   const double MS_start = 500;
   const double MS_stop = 3000;
   const int N = 50;

   for (int i = 0; i <= N; i++) {
      const double MS = MS_start + i*(MS_stop - MS_start)/N;
      const std::tuple<double, double, double, double> Mh{
         calc_Mh(0, tb, MS, xt),
         calc_Mh(1, tb, MS, xt),
         calc_Mh(2, tb, MS, xt),
         calc_Mh(3, tb, MS, xt)
      };
   }

   // BOOST_CHECK_GT(std::abs(mw1/calc_mw_SM(mh1) - 1), 1.0e-5);
   // BOOST_CHECK_CLOSE_FRACTION(mw1, calc_mw_SM(mh1), 1.0e-4);
   // BOOST_CHECK_CLOSE_FRACTION(mw2, calc_mw_SM(mh2), 1.0e-4);
   // BOOST_CHECK_CLOSE_FRACTION(mw5, calc_mw_SM(mh5), 1.0e-5);
   // BOOST_CHECK_CLOSE_FRACTION(mw10, calc_mw_SM(mh10), 1.0e-6);
}
