#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_NUHNMSSMHimalaya

#include <boost/test/unit_test.hpp>

#include "models/NUHNMSSMHimalaya/NUHNMSSMHimalaya_two_scale_spectrum_generator.hpp"
#include "threshold_corrections.hpp"
#include "scan.hpp"

#include <utility>

using namespace flexiblesusy;

constexpr double prec = 1e-5;

softsusy::QedQcd make_qedqcd()
{
   // Parameter point from arxiv:1708.05720, section 4.3
   softsusy::QedQcd qedqcd;
   qedqcd.setPoleMZ(91.1876);
   qedqcd.setPoleMW(80.384);
   qedqcd.setPoleMt(173.34);
   qedqcd.setPoleMtau(1.777);
   qedqcd.setMbMb(4.18);
   qedqcd.setAlphaEmInput(1/127.944);
   qedqcd.setAlphaSInput(0.1184);
   qedqcd.setFermiConstant(1.1663787e-5);
   qedqcd.to(qedqcd.displayPoleMZ());

   return qedqcd;
}

Spectrum_generator_settings make_settings(int loops, double scale)
{
   Spectrum_generator_settings settings;
   settings.set(Spectrum_generator_settings::precision, prec);
   settings.set(Spectrum_generator_settings::max_iterations, 10000);
   settings.set(Spectrum_generator_settings::pole_mass_loop_order, loops);
   settings.set(Spectrum_generator_settings::ewsb_loop_order, loops);
   settings.set(Spectrum_generator_settings::beta_loop_order, loops);
   settings.set(Spectrum_generator_settings::threshold_corrections_loop_order, loops);
   settings.set(Spectrum_generator_settings::force_output, true); // must force the output, because xt = sqrt(6) leads to tachyons
   settings.set(Spectrum_generator_settings::pole_mass_scale, scale);

   return settings;
}

NUHNMSSMHimalaya<Two_scale> run(int loops, double scale, const NUHNMSSMHimalaya_input_parameters& input)
{
   const softsusy::QedQcd qedqcd = make_qedqcd();
   const Spectrum_generator_settings settings = make_settings(loops, scale);

   NUHNMSSMHimalaya_spectrum_generator<Two_scale> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);

   return spectrum_generator.get_model();
}

/// calculates CP-even Higgs pole mass at given loop order for degenerate SUSY parameters
double calc_Mh(int loops, double scale, double tb, double MS, double xt, double lambda, double kappa)
{
   NUHNMSSMHimalaya_input_parameters input;
   const double MS2 = MS*MS;
   const double mu = MS;
   const double mA2 = MS2;
   const double vS = std::sqrt(2.0)*mu/lambda;
   const double sb = std::sin(std::atan(tb)); // sin(beta) = vu/v
   const double cb = std::cos(std::atan(tb)); // cos(beta) = vd/v
   const double TLambda = std::sqrt(2.0)*mA2*sb*cb/vS - mu*kappa;
   const double TKappa = TLambda;
   const double Xt = xt*MS;

   input.MSUSY = MS;
   input.M1Input = MS;
   input.M2Input = MS;
   input.M3Input = MS;
   input.MuInput = mu;
   input.TanBeta = tb;
   input.LambdaInput = lambda;
   input.KappaInput = kappa;
   input.ALambdaInput = TLambda/lambda;
   input.AKappaInput = TKappa/kappa;
   input.mq2Input << MS2, 0, 0, 0, MS2, 0, 0, 0, MS2;
   input.mu2Input << MS2, 0, 0, 0, MS2, 0, 0, 0, MS2;
   input.md2Input << MS2, 0, 0, 0, MS2, 0, 0, 0, MS2;
   input.ml2Input << MS2, 0, 0, 0, MS2, 0, 0, 0, MS2;
   input.me2Input << MS2, 0, 0, 0, MS2, 0, 0, 0, MS2;
   input.AuInput << mu/tb, 0, 0, 0, mu/tb, 0, 0, 0, Xt + mu/tb;
   input.AdInput << mu*tb, 0, 0, 0, mu*tb, 0, 0, 0, mu*tb;
   input.AeInput << mu*tb, 0, 0, 0, mu*tb, 0, 0, 0, mu*tb;

   return run(loops, scale, input).get_physical().Mhh(0);
}

std::pair<double, double> calc_Mh_DMh(int loops, double tb, double MS, double xt, double lambda, double kappa)
{
   const double Mh = calc_Mh(loops, MS, tb, MS, xt, lambda, kappa);
   const auto scales = subdivide_log(MS/2, 2*MS, 10);

   double DMh_min = std::numeric_limits<double>::max(), DMh_max = 0;

   for (const auto scale: scales) {
      const double Mh = calc_Mh(loops, scale, tb, MS, xt, lambda, kappa);
      if (Mh < DMh_min) {
         DMh_min = Mh;
      }
      if (Mh > DMh_max) {
         DMh_max = Mh;
      }
   }

   const double DMh = std::abs(DMh_max - DMh_min);

   return std::pair<double, double>(Mh, DMh);
}

/// test renormalization scale invariance of Mh pole mass
BOOST_AUTO_TEST_CASE( test_Mh_scale_invariance )
{
   const double lambda = 1e-1*prec; // MSSM-limit
   const double kappa = lambda;
   const double tb = 10;
   const double xt = -std::sqrt(6.0);

   const auto susy_scales = subdivide_log(1e3, 1e4, 10);

   for (const auto ms: susy_scales) {
      const double DMh_0l = calc_Mh_DMh(0, tb, ms, xt, lambda, kappa).second;
      const double DMh_1l = calc_Mh_DMh(1, tb, ms, xt, lambda, kappa).second;
      const double DMh_2l = calc_Mh_DMh(2, tb, ms, xt, lambda, kappa).second;
      const double DMh_3l = calc_Mh_DMh(3, tb, ms, xt, lambda, kappa).second;
      BOOST_TEST_MESSAGE("MS = " << ms << ", DMh(0l, 1l, 2l, 3l) = [" << DMh_0l << ", " << DMh_1l << ", " << DMh_2l << ", " << DMh_3l << "]");
      BOOST_CHECK_GT(DMh_1l, DMh_2l);
      BOOST_CHECK_GT(DMh_2l, DMh_3l);
   }
}
