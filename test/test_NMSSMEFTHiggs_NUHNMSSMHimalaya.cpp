// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_NMSSMEFTHiggs_NUHNMSSMHimalaya

#include <boost/test/unit_test.hpp>
#include "lowe.h"
#include "scan.hpp"
#include "NMSSMEFTHiggs_shooting_spectrum_generator.hpp"
#include "NUHNMSSMHimalaya_two_scale_spectrum_generator.hpp"

#include <fstream>

using namespace flexiblesusy;


constexpr double prec = 1e-5;
constexpr double Mt = 173.34;


/// returns A_lambda from given m_A
const double calc_Alambda(double mA, double mu, double tb, double lambda, double kappa)
{
   const double mA2 = mA*mA;
   const double sb = std::sin(std::atan(tb)); // sin(beta) = vu/v
   const double cb = std::cos(std::atan(tb)); // cos(beta) = vd/v

   return mA2*sb*cb/mu - mu*kappa/lambda;
}


softsusy::QedQcd make_qedqcd()
{
   // Parameter point from arxiv:1708.05720, section 4.3
   softsusy::QedQcd qedqcd;
   qedqcd.setPoleMZ(91.1876);
   qedqcd.setPoleMW(80.384);
   qedqcd.setPoleMt(Mt);
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
   settings.set(Spectrum_generator_settings::beta_loop_order, loops + 1);
   settings.set(Spectrum_generator_settings::threshold_corrections_loop_order, loops);
   settings.set(Spectrum_generator_settings::eft_matching_loop_order_down, loops); 
   settings.set(Spectrum_generator_settings::force_output, true); // must force the output, because xt = sqrt(6) leads to tachyons
   settings.set(Spectrum_generator_settings::pole_mass_scale, scale);

   settings.set(Spectrum_generator_settings::higgs_2loop_correction_at_as, 1);
   settings.set(Spectrum_generator_settings::higgs_2loop_correction_at_at, 1);
   settings.set(Spectrum_generator_settings::higgs_3loop_correction_at_as2, 1);

   return settings;
}


/// calculate Mh with FlexibleEFTHiggs for given input
double calc_Mh(const NMSSMEFTHiggs_input_parameters& input, int loops)
{
   const auto settings = make_settings(loops, input.MSUSY);
   const auto qedqcd = make_qedqcd();

   NMSSMEFTHiggs_spectrum_generator<Shooting> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);

   // BOOST_TEST_MESSAGE(spectrum_generator.get_model());

   return spectrum_generator.get_sm().get_physical().Mhh;
}


double calc_Mh(const NUHNMSSMHimalaya_input_parameters& input, int loops)
{
   const auto settings = make_settings(loops, input.MSUSY);
   const auto qedqcd = make_qedqcd();

   NUHNMSSMHimalaya_spectrum_generator<Two_scale> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);

   // BOOST_TEST_MESSAGE(spectrum_generator.get_model());

   return spectrum_generator.get_model().get_physical().Mhh(0);
}


NMSSMEFTHiggs_input_parameters make_point_feft(double ms, double tb, double xt, double lambda, double kappa)
{
   const double ms2 = ms*ms;
   const double mu = ms;
   const double mA = ms;
   const double Alambda = calc_Alambda(mA, mu, tb, lambda, kappa);
   const double Akappa = Alambda;
   const double Xt = xt*ms;

   NMSSMEFTHiggs_input_parameters input;
   input.MSUSY = ms;
   input.M1Input = ms;
   input.M2Input = ms;
   input.M3Input = ms;
   input.MuInput = mu;
   input.TanBeta = tb;
   input.LambdaInput = lambda;
   input.KappaInput = kappa;
   input.ALambdaInput = Alambda;
   input.AKappaInput = Akappa;
   input.mq2Input << ms2, 0, 0, 0, ms2, 0, 0, 0, ms2;
   input.mu2Input << ms2, 0, 0, 0, ms2, 0, 0, 0, ms2;
   input.md2Input << ms2, 0, 0, 0, ms2, 0, 0, 0, ms2;
   input.ml2Input << ms2, 0, 0, 0, ms2, 0, 0, 0, ms2;
   input.me2Input << ms2, 0, 0, 0, ms2, 0, 0, 0, ms2;
   input.AuInput << mu/tb, 0, 0, 0, mu/tb, 0, 0, 0, Xt + mu/tb;
   input.AdInput << mu*tb, 0, 0, 0, mu*tb, 0, 0, 0, mu*tb;
   input.AeInput << mu*tb, 0, 0, 0, mu*tb, 0, 0, 0, mu*tb;

   return input;
}


NUHNMSSMHimalaya_input_parameters make_point_fo(double ms, double tb, double xt, double lambda, double kappa)
{
   const NMSSMEFTHiggs_input_parameters input_feft = make_point_feft(ms, tb, xt, lambda, kappa);

   NUHNMSSMHimalaya_input_parameters input_fo;

   input_fo.MSUSY        = input_feft.MSUSY;
   input_fo.M1Input      = input_feft.M1Input;
   input_fo.M2Input      = input_feft.M2Input;
   input_fo.M3Input      = input_feft.M3Input;
   input_fo.MuInput      = input_feft.MuInput;
   input_fo.TanBeta      = input_feft.TanBeta;
   input_fo.LambdaInput  = input_feft.LambdaInput;
   input_fo.KappaInput   = input_feft.KappaInput;
   input_fo.ALambdaInput = input_feft.ALambdaInput;
   input_fo.AKappaInput  = input_feft.AKappaInput;
   input_fo.mq2Input     = input_feft.mq2Input;
   input_fo.mu2Input     = input_feft.mu2Input;
   input_fo.md2Input     = input_feft.md2Input;
   input_fo.ml2Input     = input_feft.ml2Input;
   input_fo.me2Input     = input_feft.me2Input;
   input_fo.AuInput      = input_feft.AuInput;
   input_fo.AdInput      = input_feft.AdInput;
   input_fo.AeInput      = input_feft.AeInput;

   return input_fo;
}


struct Data {
   double Mh_fo{}; ///< fixed-order calculation
   double Mh_feft{}; ///< FlexibleEFTHiggs calculation
};


Data calc_Mh(double ms, double tb, double xt, double lambda, double kappa, int loops)
{
   return {
      .Mh_fo   = calc_Mh(make_point_fo  (ms, tb, xt, lambda, kappa), loops),
      .Mh_feft = calc_Mh(make_point_feft(ms, tb, xt, lambda, kappa), loops),
   };
}


// test low-energy limit of the FlexibleEFTHiggs calculation
BOOST_AUTO_TEST_CASE( test_EFTHiggs_low_energy_limit )
{
   const double tb = 5;
   const double xt = 0;
   const double lambda = 0.001;
   const double kappa = 0.001;

   {
      const auto data = calc_Mh(Mt, tb, xt, lambda, kappa, 2);
      BOOST_CHECK_CLOSE_FRACTION(data.Mh_feft, data.Mh_fo, 5e-3);
   }

   {
      const auto data = calc_Mh(200, tb, xt, lambda, kappa, 3);
      BOOST_CHECK_CLOSE_FRACTION(data.Mh_feft, data.Mh_fo, 5e-3);
   }
}


// create data for plotting
BOOST_AUTO_TEST_CASE( test_EFTHiggs_plot )
{
   std::ofstream fstr("test/test_NMSSMEFTHiggs_NUHNMSSMHimalaya.dat");
   fstr << "# [1] ms | [2] Mh_fo(1L) | [3] Mh_feft(1L) | [4] Mh_fo(2L) | [5] Mh_feft(2L) | [6] Mh_fo(3L) | [7] Mh_feft(3L)\n";

   const double tb = 5;
   const double xt = 0;
   const double lambda = 0.001;
   const double kappa = 0.001;

   const auto ms_values = subdivide_log(91, 1e4, 20);

   for (const auto ms: ms_values) {
      const double Mh_fo_1l   = calc_Mh(make_point_fo  (ms, tb, xt, lambda, kappa), 1);
      const double Mh_feft_1l = calc_Mh(make_point_feft(ms, tb, xt, lambda, kappa), 1);
      const double Mh_fo_2l   = calc_Mh(make_point_fo  (ms, tb, xt, lambda, kappa), 2);
      const double Mh_feft_2l = calc_Mh(make_point_feft(ms, tb, xt, lambda, kappa), 2);
      const double Mh_fo_3l   = calc_Mh(make_point_fo  (ms, tb, xt, lambda, kappa), 3);
      const double Mh_feft_3l = calc_Mh(make_point_feft(ms, tb, xt, lambda, kappa), 3);

      fstr << ms << '\t'
           << Mh_fo_1l << '\t' << Mh_feft_1l << '\t'
           << Mh_fo_2l << '\t' << Mh_feft_2l << '\t'
           << Mh_fo_3l << '\t' << Mh_feft_3l << '\t'
           << '\n';
   }
}
