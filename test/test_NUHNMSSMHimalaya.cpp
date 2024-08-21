#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_NUHNMSSMHimalaya

#include <boost/test/unit_test.hpp>

#include "models/NUHNMSSMHimalaya/NUHNMSSMHimalaya_two_scale_spectrum_generator.hpp"
#include "models/NUHMSSMNoFVHimalaya/NUHMSSMNoFVHimalaya_two_scale_spectrum_generator.hpp"
#include "src/threshold_corrections.hpp"

using namespace flexiblesusy;

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

Threshold_corrections make_tc(int loops)
{
   // flag setting used in arxiv:1708.05720 Fig.6
   Threshold_corrections tc(122111121);
   if (loops == 2) {
      tc.set(121111121); // only 1-loop yt threshold correction
   }

   return tc;
}

Spectrum_generator_settings make_settings(int loops)
{
   const Threshold_corrections tc = make_tc(loops);

   Spectrum_generator_settings settings;
   settings.set(Spectrum_generator_settings::precision, 1.0e-5);
   settings.set(Spectrum_generator_settings::max_iterations, 10000);
   settings.set(Spectrum_generator_settings::pole_mass_loop_order, loops);
   settings.set(Spectrum_generator_settings::ewsb_loop_order, loops);
   settings.set(Spectrum_generator_settings::beta_loop_order, 3);
   settings.set(Spectrum_generator_settings::threshold_corrections_loop_order, loops);
   settings.set(Spectrum_generator_settings::force_output, true); // must force the output, because xt = sqrt(6) leads to tachyons
   settings.set_threshold_corrections(tc);

   return settings;
}

NUHMSSMNoFVHimalaya<Two_scale> run(int loops, const NUHMSSMNoFVHimalaya_input_parameters& input)
{
   const softsusy::QedQcd qedqcd = make_qedqcd();
   const Spectrum_generator_settings settings = make_settings(loops);

   NUHMSSMNoFVHimalaya_spectrum_generator<Two_scale> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);

   return spectrum_generator.get_model();
}

NUHNMSSMHimalaya<Two_scale> run(int loops, const NUHNMSSMHimalaya_input_parameters& input)
{
   const softsusy::QedQcd qedqcd = make_qedqcd();
   const Spectrum_generator_settings settings = make_settings(loops);

   NUHNMSSMHimalaya_spectrum_generator<Two_scale> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);

   return spectrum_generator.get_model();
}

/// calculates CP-even Higgs pole mass at given loop order for degenerate SUSY parameters
double calc_Mh_MSSM(int loops, double tb, double MS, double xt)
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
   input.AbIN = mu*tb;
   input.AtauIN = mu*tb;
   input.AcIN = mu/tb;
   input.AsIN = mu*tb;
   input.AmuonIN = mu*tb;
   input.AuIN = mu/tb;
   input.AdIN = mu*tb;
   input.AeIN = mu*tb;
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

   return run(loops, input).get_physical().Mhh(0);
}

/// calculates CP-even Higgs pole mass at given loop order for degenerate SUSY parameters
double calc_Mh_NMSSM(int loops, double tb, double MS, double xt, double lambda, double kappa)
{
   NUHNMSSMHimalaya_input_parameters input;
   const double MS2 = MS*MS;
   const double mu = MS;
   const double mA2 = MS2;
   const double vS = std::sqrt(2.0)*mu/lambda;
   const double sb = std::sin(std::atan(tb)); // sin(beta) = vu/v
   const double cb = std::cos(std::atan(tb)); // cos(beta) = vd/v
   const double TLambda = std::sqrt(2.0)*mA2*sb*cb/vS - mu*kappa;
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
   input.AKappaInput = 0;
   input.mq2Input << MS2, 0, 0, 0, MS2, 0, 0, 0, MS2;
   input.mu2Input << MS2, 0, 0, 0, MS2, 0, 0, 0, MS2;
   input.md2Input << MS2, 0, 0, 0, MS2, 0, 0, 0, MS2;
   input.ml2Input << MS2, 0, 0, 0, MS2, 0, 0, 0, MS2;
   input.me2Input << MS2, 0, 0, 0, MS2, 0, 0, 0, MS2;
   input.AuInput << mu/tb, 0, 0, 0, mu/tb, 0, 0, 0, Xt + mu/tb;
   input.AdInput << mu*tb, 0, 0, 0, mu*tb, 0, 0, 0, mu*tb;
   input.AeInput << mu*tb, 0, 0, 0, mu*tb, 0, 0, 0, mu*tb;

   return run(loops, input).get_physical().Mhh(0);
}

/// test 3-loop NMSSM calculation for some MSSM points from
/// arxiv:1708.05720 Fig.6 in the MSSM limit
BOOST_AUTO_TEST_CASE( test_Mh )
{
   const double lambda = 1e-4; // MSSM-limit
   const double kappa = lambda;
   const double eps = 1e-2; // @todo(alex): increase test precision, try adding 3-loop beta functions

   BOOST_CHECK_CLOSE_FRACTION(calc_Mh_MSSM(2, 5, 1e4, -std::sqrt(6.0)), calc_Mh_NMSSM(2, 5, 1e4, -std::sqrt(6.0), lambda, kappa), eps);
   BOOST_CHECK_CLOSE_FRACTION(calc_Mh_MSSM(3, 5, 1e4, -std::sqrt(6.0)), calc_Mh_NMSSM(3, 5, 1e4, -std::sqrt(6.0), lambda, kappa), eps);
   BOOST_CHECK_CLOSE_FRACTION(calc_Mh_MSSM(3, 5, 1e4, 0), calc_Mh_NMSSM(3, 5, 1e4, 0, lambda, kappa), eps);
   BOOST_CHECK_CLOSE_FRACTION(calc_Mh_MSSM(3, 5, 1e4, std::sqrt(6.0)), calc_Mh_NMSSM(3, 5, 1e4, std::sqrt(6.0), lambda, kappa), eps);
}
