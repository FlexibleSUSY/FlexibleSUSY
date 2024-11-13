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
#include "NMSSMEFTHiggs_shooting_spectrum_generator.hpp"
#include "NUHNMSSMHimalaya_two_scale_spectrum_generator.hpp"

using namespace flexiblesusy;


/// returns A_lambda from given m_A
const double calc_Alambda(double mA, double mu, double tb, double lambda, double kappa)
{
   const double mA2 = mA*mA;
   const double sb = std::sin(std::atan(tb)); // sin(beta) = vu/v
   const double cb = std::cos(std::atan(tb)); // cos(beta) = vd/v

   return mA2*sb*cb/mu - mu*kappa/lambda;
}


/// calculate Mh with FlexibleEFTHiggs for given input
double calc_Mh(const NMSSMEFTHiggs_input_parameters& input)
{
   Spectrum_generator_settings settings;
   softsusy::QedQcd qedqcd;

   settings.set(Spectrum_generator_settings::eft_matching_loop_order_down, 3); 
   settings.set(Spectrum_generator_settings::pole_mass_loop_order, 3);
   settings.set(Spectrum_generator_settings::higgs_2loop_correction_at_as, 1);
   settings.set(Spectrum_generator_settings::higgs_2loop_correction_at_at, 1);
   settings.set(Spectrum_generator_settings::higgs_3loop_correction_at_as2, 1);

   NMSSMEFTHiggs_spectrum_generator<Shooting> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);
   auto sm  = spectrum_generator.get_sm();
   return sm.get_physical().Mhh;
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
   const NMSSMEFTHiggs_input_parameters point = make_point_feft(ms, tb, xt, lambda, kappa);

   NUHNMSSMHimalaya_input_parameters input;

   input.MSUSY        = point.MSUSY;
   input.M1Input      = point.M1Input;
   input.M2Input      = point.M2Input;
   input.M3Input      = point.M3Input;
   input.MuInput      = point.MuInput;
   input.TanBeta      = point.TanBeta;
   input.LambdaInput  = point.LambdaInput;
   input.KappaInput   = point.KappaInput;
   input.ALambdaInput = point.ALambdaInput;
   input.AKappaInput  = point.AKappaInput;
   input.mq2Input     = input.mq2Input;
   input.mu2Input     = input.mu2Input;
   input.md2Input     = input.md2Input;
   input.ml2Input     = input.ml2Input;
   input.me2Input     = input.me2Input;
   input.AuInput      = input.AuInput;
   input.AdInput      = input.AdInput;
   input.AeInput      = input.AeInput;

   return input;
}


// test low-energy limit of the FlexibleEFTHiggs calculation
BOOST_AUTO_TEST_CASE( test_EFTHiggs_low_energy_limit )
{
   const double eps = 1e-5;
   const double ms = 200;

   // @todo(alex): scan over low values of ms and test against FO calculation
   const double Mh_feft = calc_Mh(make_point_feft(ms, 5, 0, 0.001, 0.001));
   BOOST_CHECK_CLOSE_FRACTION(Mh_feft, 88.485127104418012, eps);
}
