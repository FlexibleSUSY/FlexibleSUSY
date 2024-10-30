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

using Output = std::array<double,1>;


/// calculate Mh at the precision given in `settings'
double calc_Mh(
   NMSSMEFTHiggs_input_parameters& input,
   softsusy::QedQcd& qedqcd,
   Spectrum_generator_settings& settings)
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
Output calc_output(char const* const slha_input)
{
   Spectrum_generator_settings settings;
   softsusy::QedQcd qedqcd;
   NMSSMEFTHiggs_input_parameters input;

   std::tie(settings, qedqcd, input) = extract_slha_input(slha_input);

   Output results{};

   results.at(0) = calc_Mh(input, qedqcd, settings);

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


BOOST_AUTO_TEST_CASE( test_top_down_EFTHiggs )
{
   const struct Data {
      char const * const slha_input = nullptr;
      Output expected_output{};
      double eps{0.0};
   } data[] = {
      {slha_input_case_1, {1.08559604e+02}, 2e-4}, // obtained from NMSSMEFTHiggsTwoScale in EFT parametrization
      {slha_input_case_2, {133.0357875507616}, 1e-5}, // obtained from MSSMEFTHiggs2loop in full-model parametrization w/ only 2-loop contributions of O((at+ab)*as + (at+ab)^2), i.e. no 2-loop O(atau^2) contributions and not 3- or 4-loop contributions
   };

   for (const auto& d: data) {
      const auto output = calc_output(d.slha_input);
      for (int i = 0; i < d.expected_output.size(); i++) {
         BOOST_CHECK_CLOSE_FRACTION(output[i], d.expected_output[i], d.eps);
      }
   }
}
