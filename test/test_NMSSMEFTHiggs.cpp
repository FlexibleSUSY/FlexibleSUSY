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
#define BOOST_TEST_MODULE test_NMSSMEFTHiggs

#include <boost/test/unit_test.hpp>
#include <tuple>
#include "lowe.h"
#include "NMSSMEFTHiggs_shooting_spectrum_generator.hpp"
#include "NMSSMEFTHiggs_slha_io.hpp"
#include "wrappers.hpp"

using namespace flexiblesusy;

struct Output_2loop {
   double Mh_2L_at_as{};
   double Mh_2L_at_at{};
};


/// calculate Mh at the precision given in `settings'
double calc_Mh(
   const NMSSMEFTHiggs_input_parameters& input,
   const softsusy::QedQcd& qedqcd,
   const Spectrum_generator_settings& settings)
{
   NMSSMEFTHiggs_spectrum_generator<Shooting> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);
   auto sm  = spectrum_generator.get_sm();

   const double Q_pole = settings.get(Spectrum_generator_settings::eft_pole_mass_scale) != 0. ? settings.get(Spectrum_generator_settings::eft_pole_mass_scale) :  qedqcd.displayPoleMt();

   sm.run_to(Q_pole);
   sm.solve_ewsb();
   sm.calculate_Mhh_pole();
   return sm.get_physical().Mhh;
}


/// extracts SLHA input
std::tuple<Spectrum_generator_settings, softsusy::QedQcd, NMSSMEFTHiggs_input_parameters>
extract_slha_input(char const * const slha_input)
{
   std::stringstream istr(slha_input);
   NMSSMEFTHiggs_slha_io slha_io;
   slha_io.read_from_stream(istr);

   softsusy::QedQcd qedqcd;
   NMSSMEFTHiggs_input_parameters input;
   Spectrum_generator_settings settings;

   try {
      slha_io.fill(settings);
      slha_io.fill(qedqcd);
      slha_io.fill(input);
   } catch (const Error& error) {
      BOOST_TEST_MESSAGE(error.what());
      BOOST_TEST(false);
   }

   return {settings, qedqcd, input};
}


/// calculate output for test
Output_2loop calc_output_2loop(char const* const slha_input)
{
   Spectrum_generator_settings settings;
   softsusy::QedQcd qedqcd;
   NMSSMEFTHiggs_input_parameters input;

   std::tie(settings, qedqcd, input) = extract_slha_input(slha_input);

   Output_2loop results{};

   settings.set(Spectrum_generator_settings::higgs_2loop_correction_at_as, 1);
   settings.set(Spectrum_generator_settings::higgs_2loop_correction_at_at, 0);
   results.Mh_2L_at_as = calc_Mh(input, qedqcd, settings);

   settings.set(Spectrum_generator_settings::higgs_2loop_correction_at_as, 0);
   settings.set(Spectrum_generator_settings::higgs_2loop_correction_at_at, 1);
   results.Mh_2L_at_at = calc_Mh(input, qedqcd, settings);

   return results;
}


// scenario 1 degenrate case at vanishing stop mixing xt = 0
char const * const slha_input_case_1 = R"(
Block FlexibleSUSY
    0   1.000000000e-05      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = all, 1 = two_scale, 2 = semi_analytic)
    3   0                    # calculate SM pole masses
    4   2                    # pole mass loop order
    5   1                    # EWSB loop order
    6   4                    # beta-functions loop order
    7   2                    # threshold corrections loop order
    8   1                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   1                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   1                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   0                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   0                    # force output
   13   1                    # Top pole mass QCD corrections (0 = 1L, 1 = 2L, 2 = 3L)
   14   1.000000000e-11      # beta-function zero threshold
   15   0                    # calculate observables (a_muon, ...)
   16   0                    # force positive majorana masses
   17   0                    # pole mass renormalization scale (0 = SUSY scale)
   18   0                    # pole mass renormalization scale in the EFT (0 = min(SUSY scale, Mt))
   19   0                    # EFT matching scale (0 = SUSY scale)
   20   1                    # EFT loop order for upwards matching
   21   1                    # EFT loop order for downwards matching
   22   0                    # EFT index of SM-like Higgs in the BSM model
   23   0                    # calculate BSM pole masses
   24   123111321            # individual threshold correction loop orders
   25   0                    # ren. scheme for Higgs 3L corrections (0 = DR, 1 = MDR)
   26   0                    # Higgs 3-loop corrections O(alpha_t alpha_s^2)
   27   0                    # Higgs 3-loop corrections O(alpha_b alpha_s^2)
   28   0                    # Higgs 3-loop corrections O(alpha_t^2 alpha_s)
   29   0                    # Higgs 3-loop corrections O(alpha_t^3)
   30   0                    # Higgs 4-loop corrections O(alpha_t alpha_s^3)
   31   0                    # loop library (0 = softsusy)
Block SMINPUTS               # Standard Model inputs
    1   1.279160000e+02      # alpha^(-1) SM MSbar(MZ)
    2   1.166378700e-05      # G_Fermi
    3   1.184000000e-01      # alpha_s(MZ) SM MSbar
    4   9.118760000e+01      # MZ(pole)
    5   4.180000000e+00      # mb(mb) SM MSbar
    6   1.733400000e+02      # mtop(pole)
    7   1.776990000e+00      # mtau(pole)
    8   0.000000000e+00      # mnu3(pole)
    9   80.385               # MW pole
   11   5.109989020e-04      # melectron(pole)
   12   0.000000000e+00      # mnu1(pole)
   13   1.056583570e-01      # mmuon(pole)
   14   0.000000000e+00      # mnu2(pole)
   21   4.750000000e-03      # md(2 GeV) MS-bar
   22   2.400000000e-03      # mu(2 GeV) MS-bar
   23   1.040000000e-01      # ms(2 GeV) MS-bar
   24   1.270000000e+00      # mc(mc) MS-bar
Block EXTPAR
    0   2000                 # Ms
    1   2000                 # M1(MSUSY)
    2   2000                 # M2(MSUSY)
    3   2000                 # M3(MSUSY)
    4   2000                 # Mu(MSUSY)
   25   5                    # tan(beta) at Ms
   61   0.2                  # Lambda
   62   0.2                  # Kappa
   63  -1.61538462E+03       # ALambda
   64  -1.00000000E+03       # AKappa
Block MSQ2IN
  1  1     4.00000000E+06   # mq2(1,1)
  2  2     4.00000000E+06   # mq2(2,2)
  3  3     4.00000000E+06   # mq2(3,3)
Block MSE2IN
  1  1     4.00000000E+06   # me2(1,1)
  2  2     4.00000000E+06   # me2(2,2)
  3  3     4.00000000E+06   # me2(3,3)
Block MSL2IN
  1  1     4.00000000E+06   # ml2(1,1)
  2  2     4.00000000E+06   # ml2(2,2)
  3  3     4.00000000E+06   # ml2(3,3)
Block MSU2IN
  1  1     4.00000000E+06   # mu2(1,1)
  2  2     4.00000000E+06   # mu2(2,2)
  3  3     4.00000000E+06   # mu2(3,3)
Block MSD2IN
  1  1     4.00000000E+06   # md2(1,1)
  2  2     4.00000000E+06   # md2(2,2)
  3  3     4.00000000E+06   # md2(3,3)
Block AUIN
  1  1     0   # Au(1,1)
  2  2     0   # Au(2,2)
  3  3     400 # Au(3,3) xt=0
Block ADIN
  1  1     0   # Ad(1,1)
  2  2     0   # Ad(2,2)
  3  3     0   # Ad(3,3)
Block AEIN
  1  1     0   # Ad(1,1)
  2  2     0   # Ad(2,2)
  3  3     0   # Ad(3,3)
)";

// scenario 2 degenrate case with non-trivial mixing At = xt*ms + mu/tb = 0 + 2000/5 = 400 and approx MSSM-limit
char const * const slha_input_case_2 = R"(
Block MODSEL                 # Select model
#   12    1000                # DRbar parameter output scale (GeV)
Block FlexibleSUSY
    0   1.000000000e-05      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = all, 1 = two_scale, 2 = semi_analytic)
    3   0                    # calculate SM pole masses
    4   3                    # pole mass loop order
    5   3                    # EWSB loop order
    6   3                    # beta-functions loop order
    7   2                    # threshold corrections loop order
    8   1                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   1                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   1                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   0                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   0                    # force output
   13   1                    # Top pole mass QCD corrections (0 = 1L, 1 = 2L, 2 = 3L)
   14   1.000000000e-11      # beta-function zero threshold
   15   0                    # calculate observables (a_muon, ...)
   16   0                    # force positive majorana masses
   17   0                    # pole mass renormalization scale (0 = SUSY scale)
   18   0                    # pole mass renormalization scale in the EFT (0 = min(SUSY scale, Mt))
   19   0                    # EFT matching scale (0 = SUSY scale)
   20   1                    # EFT loop order for upwards matching
   21   2                    # EFT loop order for downwards matching
   22   0                    # EFT index of SM-like Higgs in the BSM model
   23   0                    # calculate BSM pole masses
   24   123111321            # individual threshold correction loop orders
   25   0                    # ren. scheme for Higgs 3L corrections (0 = DR, 1 = MDR)
   26   0                    # Higgs 3-loop corrections O(alpha_t alpha_s^2)
   27   0                    # Higgs 3-loop corrections O(alpha_b alpha_s^2)
   28   0                    # Higgs 3-loop corrections O(alpha_t^2 alpha_s)
   29   0                    # Higgs 3-loop corrections O(alpha_t^3)
   30   0                    # Higgs 4-loop corrections O(alpha_t alpha_s^3)
   31   0                    # loop library (0 = softsusy)
Block SMINPUTS               # Standard Model inputs
    1   1.279160000e+02      # alpha^(-1) SM MSbar(MZ)
    2   1.166378700e-05      # G_Fermi
    3   1.184000000e-01      # alpha_s(MZ) SM MSbar
    4   9.118760000e+01      # MZ(pole)
    5   4.180000000e+00      # mb(mb) SM MSbar
    6   1.733400000e+02      # mtop(pole)
    7   1.776990000e+00      # mtau(pole)
    8   0.000000000e+00      # mnu3(pole)
    9   80.385               # MW pole
   11   5.109989020e-04      # melectron(pole)
   12   0.000000000e+00      # mnu1(pole)
   13   1.056583570e-01      # mmuon(pole)
   14   0.000000000e+00      # mnu2(pole)
   21   4.750000000e-03      # md(2 GeV) MS-bar
   22   2.400000000e-03      # mu(2 GeV) MS-bar
   23   1.040000000e-01      # ms(2 GeV) MS-bar
   24   1.270000000e+00      # mc(mc) MS-bar
Block EXTPAR
    0   2000                 # Ms
    1   2000                 # M1(MSUSY)
    2   2000                 # M2(MSUSY)
    3   2000                 # M3(MSUSY)
    4   2000                 # Mu(MSUSY)
   25   5                    # tan(beta) at Ms
   61   0.001                  # Lambda
   62   0.001                  # Kappa
   63  -1.61538462E+03       # ALambda
   64  -1.00000000E+03       # AKappa
Block MSQ2IN
  1  1     4.00000000E+06   # mq2(1,1)
  2  2     4.00000000E+06   # mq2(2,2)
  3  3     4.00000000E+06   # mq2(3,3)
Block MSE2IN
  1  1     4.00000000E+06   # me2(1,1)
  2  2     4.00000000E+06   # me2(2,2)
  3  3     4.00000000E+06   # me2(3,3)
Block MSL2IN
  1  1     4.00000000E+06   # ml2(1,1)
  2  2     4.00000000E+06   # ml2(2,2)
  3  3     4.00000000E+06   # ml2(3,3)
Block MSU2IN
  1  1     4.00000000E+06   # mu2(1,1)
  2  2     4.00000000E+06   # mu2(2,2)
  3  3     4.00000000E+06   # mu2(3,3)
Block MSD2IN
  1  1     4.00000000E+06   # md2(1,1)
  2  2     4.00000000E+06   # md2(2,2)
  3  3     4.00000000E+06   # md2(3,3)
Block AUIN
  1  1     0   # Au(1,1)
  2  2     0   # Au(2,2)
  3  3     400   # Au(3,3) xt=0
Block ADIN
  1  1     0   # Ad(1,1)
  2  2     0   # Ad(2,2)
  3  3     0   # Ad(3,3)
Block AEIN
  1  1     0   # Ad(1,1)
  2  2     0   # Ad(2,2)
  3  3     0   # Ad(3,3)
)";

// scenario 3 degenrate case with non-trivial mixing At = xt*ms + mu/tb = (-2.2 + 1/5)*2000 = -4000 and approx MSSM-limit
char const * const slha_input_case_3 = R"(
Block MODSEL                 # Select model
#   12    1000                # DRbar parameter output scale (GeV)
Block FlexibleSUSY
    0   1.000000000e-05      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = all, 1 = two_scale, 2 = semi_analytic)
    3   0                    # calculate SM pole masses
    4   3                    # pole mass loop order
    5   3                    # EWSB loop order
    6   3                    # beta-functions loop order
    7   2                    # threshold corrections loop order
    8   1                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   1                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   1                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   0                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   0                    # force output
   13   1                    # Top pole mass QCD corrections (0 = 1L, 1 = 2L, 2 = 3L)
   14   1.000000000e-11      # beta-function zero threshold
   15   0                    # calculate observables (a_muon, ...)
   16   0                    # force positive majorana masses
   17   0                    # pole mass renormalization scale (0 = SUSY scale)
   18   0                    # pole mass renormalization scale in the EFT (0 = min(SUSY scale, Mt))
   19   0                    # EFT matching scale (0 = SUSY scale)
   20   1                    # EFT loop order for upwards matching
   21   2                    # EFT loop order for downwards matching
   22   0                    # EFT index of SM-like Higgs in the BSM model
   23   0                    # calculate BSM pole masses
   24   123111321            # individual threshold correction loop orders
   25   0                    # ren. scheme for Higgs 3L corrections (0 = DR, 1 = MDR)
   26   0                    # Higgs 3-loop corrections O(alpha_t alpha_s^2)
   27   0                    # Higgs 3-loop corrections O(alpha_b alpha_s^2)
   28   0                    # Higgs 3-loop corrections O(alpha_t^2 alpha_s)
   29   0                    # Higgs 3-loop corrections O(alpha_t^3)
   30   0                    # Higgs 4-loop corrections O(alpha_t alpha_s^3)
   31   0                    # loop library (0 = softsusy)
Block SMINPUTS               # Standard Model inputs
    1   1.279160000e+02      # alpha^(-1) SM MSbar(MZ)
    2   1.166378700e-05      # G_Fermi
    3   1.184000000e-01      # alpha_s(MZ) SM MSbar
    4   9.118760000e+01      # MZ(pole)
    5   4.180000000e+00      # mb(mb) SM MSbar
    6   1.733400000e+02      # mtop(pole)
    7   1.776990000e+00      # mtau(pole)
    8   0.000000000e+00      # mnu3(pole)
    9   80.385               # MW pole
   11   5.109989020e-04      # melectron(pole)
   12   0.000000000e+00      # mnu1(pole)
   13   1.056583570e-01      # mmuon(pole)
   14   0.000000000e+00      # mnu2(pole)
   21   4.750000000e-03      # md(2 GeV) MS-bar
   22   2.400000000e-03      # mu(2 GeV) MS-bar
   23   1.040000000e-01      # ms(2 GeV) MS-bar
   24   1.270000000e+00      # mc(mc) MS-bar
Block EXTPAR
    0   2000                 # Ms
    1   2000                 # M1(MSUSY)
    2   2000                 # M2(MSUSY)
    3   2000                 # M3(MSUSY)
    4   2000                 # Mu(MSUSY)
   25   5                    # tan(beta) at Ms
   61   0.001                # Lambda
   62   0.001                # Kappa
   63  -1.61538462E+03       # ALambda
   64  -1.00000000E+03       # AKappa
Block MSQ2IN
  1  1     4.00000000E+06   # mq2(1,1)
  2  2     4.00000000E+06   # mq2(2,2)
  3  3     4.00000000E+06   # mq2(3,3)
Block MSE2IN
  1  1     4.00000000E+06   # me2(1,1)
  2  2     4.00000000E+06   # me2(2,2)
  3  3     4.00000000E+06   # me2(3,3)
Block MSL2IN
  1  1     4.00000000E+06   # ml2(1,1)
  2  2     4.00000000E+06   # ml2(2,2)
  3  3     4.00000000E+06   # ml2(3,3)
Block MSU2IN
  1  1     4.00000000E+06   # mu2(1,1)
  2  2     4.00000000E+06   # mu2(2,2)
  3  3     4.00000000E+06   # mu2(3,3)
Block MSD2IN
  1  1     4.00000000E+06   # md2(1,1)
  2  2     4.00000000E+06   # md2(2,2)
  3  3     4.00000000E+06   # md2(3,3)
Block AUIN
  1  1     0   # Au(1,1)
  2  2     0   # Au(2,2)
  3  3 -4000   # Au(3,3) xt=-2.2
Block ADIN
  1  1     0   # Ad(1,1)
  2  2     0   # Ad(2,2)
  3  3     0   # Ad(3,3)
Block AEIN
  1  1     0   # Ad(1,1)
  2  2     0   # Ad(2,2)
  3  3     0   # Ad(3,3)
)";


BOOST_AUTO_TEST_CASE( test_top_down_EFTHiggs_2loop )
{
   const struct Data {
      char const * const slha_input = nullptr;
      Output_2loop expected_output{};
      double eps{0.0};
   } data[] = {
      {slha_input_case_1, { .Mh_2L_at_as = 108.71399186416657, .Mh_2L_at_at = 108.54840607805137 }, 2e-5}, // obtained from NMSSMEFTHiggsTwoScale in EFT parametrization
      {slha_input_case_2, { .Mh_2L_at_as = 111.23847635109462, .Mh_2L_at_at = 111.18889782272626 }, 2e-5}, // MS = 2000 GeV, tb = 5, xt = 0, MSSM-limit
      {slha_input_case_3, { .Mh_2L_at_as = 120.14850196112329, .Mh_2L_at_at = 119.79021153724737 }, 2e-5}, // MS = 2000 GeV, tb = 5, xt = 0, MSSM-limit, obtained from MSSMEFTHiggs2loop in full-model parametrization w/ only 2-loop contributions of O((at+ab)*as + (at+ab)^2), i.e. no 2-loop O(atau^2) contributions and not 3- or 4-loop contributions
   };

   for (const auto& d: data) {
      const auto output = calc_output_2loop(d.slha_input);
      BOOST_CHECK_CLOSE_FRACTION(output.Mh_2L_at_as, d.expected_output.Mh_2L_at_as, d.eps);
      BOOST_CHECK_CLOSE_FRACTION(output.Mh_2L_at_at, d.expected_output.Mh_2L_at_at, d.eps);
   }
}


Spectrum_generator_settings make_settings(int loops)
{
   Spectrum_generator_settings settings;

   settings.set(Spectrum_generator_settings::precision, 1e-5);
   settings.set(Spectrum_generator_settings::pole_mass_loop_order, loops);
   settings.set(Spectrum_generator_settings::ewsb_loop_order, loops);
   settings.set(Spectrum_generator_settings::beta_loop_order, loops + 1);
   settings.set(Spectrum_generator_settings::threshold_corrections, loops);
   settings.set(Spectrum_generator_settings::eft_pole_mass_scale, 0);
   settings.set(Spectrum_generator_settings::eft_matching_scale, 0);
   settings.set(Spectrum_generator_settings::eft_matching_loop_order_down, loops);
   settings.set(Spectrum_generator_settings::calculate_bsm_masses, 0);

   return settings;
}


NMSSMEFTHiggs_input_parameters make_point(double ms, double tb, double xt, double lambda, double kappa)
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
   input.ALambdaInput = 0;
   input.AKappaInput = 0;
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


/// calculates uncertainty of Mh
double calc_DMh(double ms, double tb, double xt, double lambda, double kappa, int loops)
{
   softsusy::QedQcd qedqcd;

   const NMSSMEFTHiggs_input_parameters input = make_point(ms, tb, xt, lambda, kappa);

   const double Mh = calc_Mh(input, qedqcd, make_settings(loops));

   // variation of matchnig scale Qmatch within [ms/2, 2*ms], Eq.(8.1) arxiv:2003.04639
   const double DMh_Qmatch = [&] () {
      const double t_min = std::log(ms/2);
      const double t_max = std::log(2*ms);
      const int N_scales = 10;
      double DMh_Q = 0;

      for (int i = 0; i <= N_scales; ++i) {
         const double Q = std::exp(t_min + i*(t_max - t_min)/N_scales);
         auto settings = make_settings(loops);
         settings.set(Spectrum_generator_settings::eft_matching_scale, Q);
         const double Mh_at_Q = calc_Mh(input, qedqcd, settings);
         const double diff = std::abs(Mh - Mh_at_Q); // Eq.(8.1) arxiv:2003.04639
         if (diff > DMh_Q) {
            DMh_Q = diff;
         }
      }

      return DMh_Q;
   }();

   return DMh_Qmatch;
}


template<class T>
void print_vector(const std::vector<T>& vec)
{
   std::stringstream ss;

   ss << "[";
   for (const auto& e: vec) {
      ss << e << ", ";
   }
   ss << "]";

   BOOST_TEST_MESSAGE(ss.str());
}


void test_uncertainty_MSUSY_scan(double tb, double xt, const std::vector<double>& scales)
{
   BOOST_TEST_MESSAGE("test_uncertainty_MSUSY_scan: tb = " << tb << ", xt = " << xt);

   std::vector<double> uncertainties_1l, uncertainties_2l, uncertainties_3l;

   for (const auto& ms: scales) {
      uncertainties_1l.push_back(calc_DMh(ms, tb, xt, 0.001, 0.001, 1));
      uncertainties_2l.push_back(calc_DMh(ms, tb, xt, 0.001, 0.001, 2));
      uncertainties_3l.push_back(calc_DMh(ms, tb, xt, 0.001, 0.001, 3));
   }

   BOOST_TEST_MESSAGE("scales  = "); print_vector(scales);
   BOOST_TEST_MESSAGE("DMh(1L) = "); print_vector(uncertainties_1l);
   BOOST_TEST_MESSAGE("DMh(2L) = "); print_vector(uncertainties_2l);
   BOOST_TEST_MESSAGE("DMh(3L) = "); print_vector(uncertainties_3l);

   // check that uncertainties become smaller if MSUSY is increased
   BOOST_CHECK(std::is_sorted(uncertainties_1l.begin(), uncertainties_1l.end(), std::greater<>{}));
   BOOST_CHECK(std::is_sorted(uncertainties_2l.begin(), uncertainties_2l.end(), std::greater<>{}));
   BOOST_CHECK(std::is_sorted(uncertainties_3l.begin(), uncertainties_3l.end(), std::greater<>{}));

   // Check that uncertainties become smaller if number of loops is
   // increased.  Note: For xt = 0 the uncertainty of the 1-loop
   // calculation is unnaturally small, so we omit it from the test.
   if (xt != 0) {
      BOOST_CHECK(std::equal(uncertainties_2l.begin(), uncertainties_2l.end(),
                             uncertainties_3l.begin(), uncertainties_3l.end(),
                             std::greater<>{}));
   }
}


BOOST_AUTO_TEST_CASE( test_uncertainty )
{
   test_uncertainty_MSUSY_scan(20, 0, {400, 700, 1e3, 3e3, 1e4});
   test_uncertainty_MSUSY_scan(20, -std::sqrt(6.0), {2e3, 5e3, 7e3, 1e4});
   test_uncertainty_MSUSY_scan(20, std::sqrt(6.0), {2e3, 5e3, 7e3, 1e4});
}
