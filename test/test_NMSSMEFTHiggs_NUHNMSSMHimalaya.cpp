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

   NMSSMEFTHiggs_input_parameters input;
   input.MSUSY = ms;
   input.M1Input = ms;
   input.M2Input = ms;
   input.M3Input = ms;
   input.MuInput = ms;
   input.TanBeta = tb;
   input.LambdaInput = lambda;
   input.KappaInput = kappa;
   input.ALambdaInput = 0; // @todo(alex): make consistent with FO calculation
   input.AKappaInput = input.ALambdaInput;
   input.mq2Input << ms2, 0, 0, 0, ms2, 0, 0, 0, ms2;
   input.mu2Input << ms2, 0, 0, 0, ms2, 0, 0, 0, ms2;
   input.md2Input << ms2, 0, 0, 0, ms2, 0, 0, 0, ms2;
   input.ml2Input << ms2, 0, 0, 0, ms2, 0, 0, 0, ms2;
   input.me2Input << ms2, 0, 0, 0, ms2, 0, 0, 0, ms2;
   input.AuInput << 0, 0, 0, 0, 0, 0, 0, 0, ms*(xt + 1/tb);
   input.AdInput << 0, 0, 0, 0, 0, 0, 0, 0, 0;
   input.AeInput << 0, 0, 0, 0, 0, 0, 0, 0, 0;

   return input;
}


NUHNMSSMHimalaya_input_parameters make_point_fo(double ms, double tb, double xt, double lambda, double kappa)
{
   NUHNMSSMHimalaya_input_parameters input;
   const double ms2 = ms*ms;
   const double mu = ms;
   const double mA2 = ms2;
   const double vS = std::sqrt(2.0)*mu/lambda;
   const double sb = std::sin(std::atan(tb)); // sin(beta) = vu/v
   const double cb = std::cos(std::atan(tb)); // cos(beta) = vd/v
   const double TLambda = std::sqrt(2.0)*mA2*sb*cb/vS - mu*kappa;
   const double TKappa = TLambda;
   const double Xt = xt*ms;

   input.MSUSY = ms;
   input.M1Input = ms;
   input.M2Input = ms;
   input.M3Input = ms;
   input.MuInput = mu;
   input.TanBeta = tb;
   input.LambdaInput = lambda;
   input.KappaInput = kappa;
   input.ALambdaInput = TLambda/lambda;
   input.AKappaInput = TKappa/kappa;
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


// test low-energy limit of the FlexibleEFTHiggs calculation
BOOST_AUTO_TEST_CASE( test_EFTHiggs_low_energy_limit )
{
   const double eps = 1e-5;
   const double ms = 200;

   // @todo(alex): scan over low values of ms and test against FO calculation
   const double Mh_feft = calc_Mh(make_point_feft(ms, 5, 0, 0.001, 0.001));
   BOOST_CHECK_CLOSE_FRACTION(Mh_feft, 90.126510374311536, eps);
}
