#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_NUHMSSMNoFVHimalaya

#include <boost/test/unit_test.hpp>

#include "models/NUHMSSMNoFVHimalaya/NUHMSSMNoFVHimalaya_two_scale_spectrum_generator.hpp"

#include <cstdio>
#include <vector>
#include <tuple>

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
   const double mu = MS;
   const double Xt = xt*MS;

   input.TanBeta = tb;
   input.Qin = MS;
   input.M1 = MS;
   input.M2 = MS;
   input.M3 = MS;
   input.AtIN = Xt + mu/tb;
   input.AbIN = 0;
   input.AtauIN = 0;
   input.AcIN = 0;
   input.AsIN = 0;
   input.AmuonIN = 0;
   input.AuIN = 0;
   input.AdIN = 0;
   input.AeIN = 0;
   input.MuIN = mu;
   input.mA2IN = MS*MS;
   input.ml11IN = MS;
   input.ml22IN = MS;
   input.ml33IN = MS;
   input.me11IN = MS;
   input.me22IN = MS;
   input.me33IN = MS;
   input.mq11IN = MS;
   input.mq22IN = MS;
   input.mq33IN = MS;
   input.mu11IN = MS;
   input.mu22IN = MS;
   input.mu33IN = MS;
   input.md11IN = MS;
   input.md22IN = MS;
   input.md33IN = MS;
   input.Mlow = 0;

   return run(input).get_physical().Mhh(0);
}

struct Mh {
   double MS{}, Mh0L{}, Mh1L{}, Mh2L{}, Mh3L{};
};

BOOST_AUTO_TEST_CASE( test_Mh )
{
   std::vector<Mh> data;

   const double tb = 20;
   const double xt = -std::sqrt(6.0);
   const double MS_start = 500;
   const double MS_stop = 3000;
   const int N = 10;

   for (int i = 0; i <= N; i++) {
      const double MS = MS_start + i*(MS_stop - MS_start)/N;
      data.push_back({
            MS,
            calc_Mh(0, tb, MS, xt),
            calc_Mh(1, tb, MS, xt),
            calc_Mh(2, tb, MS, xt),
            calc_Mh(3, tb, MS, xt)
         });
   }

   std::printf("# tb = %f, xt = %f\n", tb,  xt);
   std::printf("# %5s %7s %7s %7s %7s\n", "MS",  "Mh0L", "Mh1L", "Mh2L", "Mh3L");

   for (const auto& d: data) {
      std::printf("%7.3f %7.3f %7.3f %7.3f %7.3f\n", d.MS,  d.Mh0L, d.Mh1L, d.Mh2L, d.Mh3L);
   }
}
