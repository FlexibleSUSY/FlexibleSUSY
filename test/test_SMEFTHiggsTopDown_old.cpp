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
#define BOOST_TEST_MODULE test_SMEFTHiggsTopDown

#include <boost/test/unit_test.hpp>
#include <array>
#include <cstdlib>

#include "lowe.h"
#include "wrappers.hpp"
#include "spectrum_generator_settings.hpp"

#include "SM_slha_io.hpp"
#include "SM_two_scale_spectrum_generator.hpp"
#include "SMEFTHiggsTopDown_slha_io.hpp"
#include "SMEFTHiggsTopDown_shooting_spectrum_generator.hpp"

using namespace flexiblesusy;

char const * const slha_input = R"(
Block MODSEL                 # Select model
   12   100                  # parameter output scale (GeV)
Block FlexibleSUSY
    0   1.000000000e-05      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # solver (0 = all, 1 = two_scale, 2 = semi_analytic)
    3   1                    # calculate SM pole masses
    4   4                    # pole mass loop order
    5   4                    # EWSB loop order
    6   4                    # beta-functions loop order
    7   4                    # threshold corrections loop order
    8   1                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   1                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   1                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   1                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   0                    # force output
   13   3                    # Top quark 2-loop corrections QCD
   14   1.000000000e-11      # beta-function zero threshold
   15   0                    # calculate observables (a_muon, ...)
   16   0                    # force positive majorana masses
   17   0                    # pole mass renormalization scale (0 = SUSY scale)
   18   0                    # pole mass renormalization scale in the EFT (0 = min(SUSY scale, Mt))
   19   0                    # EFT matching scale (0 = SUSY scale)
   20   0                    # EFT loop order for upwards matching
   21   0                    # EFT loop order for downwards matching
   22   0                    # EFT index of SM-like Higgs in the BSM model
   23   1                    # calculate BSM pole masses
   24   124111321            # individual threshold correction loop orders
   25   0                    # ren. scheme for Higgs 3L corrections (0 = DR', 1 = MDR', 2 = H3m)
   26   1                    # Higgs 3-loop corrections O(alpha_t alpha_s^2)
   27   1                    # Higgs 3-loop corrections O(alpha_b alpha_s^2)
   28   1                    # Higgs 3-loop corrections O(alpha_t^2 alpha_s)
   29   1                    # Higgs 3-loop corrections O(alpha_t^3)
   30   1                    # Higgs 4-loop corrections O(alpha_t alpha_s^3)
   28   1                    # Higgs 3-loop corrections O(alpha_t^2 alpha_s)
   29   1                    # Higgs 3-loop corrections O(alpha_t^3)
   30   1                    # Higgs 4-loop corrections O(alpha_t alpha_s^3)
Block SMINPUTS               # Standard Model inputs
    1   1.279340000e+02      # alpha^(-1) SM MSbar(MZ)
    2   1.166378700e-05      # G_Fermi
    3   1.176000000e-01      # alpha_s(MZ) SM MSbar
    4   9.118760000e+01      # MZ(pole)
    5   4.200000000e+00      # mb(mb) SM MSbar
    6   1.733400000e+02      # mtop(pole)
    7   1.777000000e+00      # mtau(pole)
    8   0.000000000e+00      # mnu3(pole)
    9   80.404               # MW pole
   11   5.109989020e-04      # melectron(pole)
   12   0.000000000e+00      # mnu1(pole)
   13   1.056583570e-01      # mmuon(pole)
   14   0.000000000e+00      # mnu2(pole)
   21   4.750000000e-03      # md(2 GeV) MS-bar
   22   2.400000000e-03      # mu(2 GeV) MS-bar
   23   1.040000000e-01      # ms(2 GeV) MS-bar
   24   1.270000000e+00      # mc(mc) MS-bar
Block MINPAR                 # Input parameters
    1   0.192                # Lambda(Qin)
Block EXTPAR                 # Input parameters
    0   1000                 # input scale Qin
    1   173.34               # scale QEWSB
)";

using Output_t = std::array<double, 9>;

#define DEFINE_FUNCTION(model_name, solver_type, model_type)   \
   Output_t calc_ ## model_name(char const* const slha_input)  \
   {                                                           \
      using namespace flexiblesusy;                            \
                                                               \
      std::stringstream istr(slha_input);                      \
                                                               \
      model_name ## _slha_io slha_io;                          \
      slha_io.read_from_stream(istr);                          \
                                                               \
      softsusy::QedQcd qedqcd;                                 \
      model_name ## _input_parameters input;                   \
      Spectrum_generator_settings settings;                    \
                                                               \
      try {                                                    \
         slha_io.fill(settings);                               \
         slha_io.fill(qedqcd);                                 \
         slha_io.fill(input);                                  \
      } catch (const Error& error) {                           \
         BOOST_TEST_MESSAGE(error.what());                     \
         BOOST_TEST(false);                                    \
      }                                                        \
                                                               \
      model_name ## _spectrum_generator<solver_type> sg;       \
      sg.set_settings(settings);                               \
      sg.set_parameter_output_scale(100.);                     \
                                                               \
      BOOST_REQUIRE_NO_THROW(sg.run(qedqcd, input));           \
                                                               \
      const auto models = sg.get_models();                     \
      const auto sm = std::get<model_type>(models);            \
      const auto pr = sg.get_problems();                       \
                                                               \
      BOOST_TEST(!pr.have_problem());                          \
                                                               \
      if (pr.have_problem()) {                                 \
         BOOST_TEST_MESSAGE(pr);                               \
      }                                                        \
                                                               \
      return Output_t{                                         \
         sm.get_physical().Mhh,                                \
         sm.get_Lambdax(),                                     \
         sm.get_g1(),                                          \
         sm.get_g2(),                                          \
         sm.get_g3(),                                          \
         sm.get_Yu(2,2),                                       \
         sm.get_Yd(2,2),                                       \
         sm.get_Ye(2,2),                                       \
         sm.get_v()                                            \
      };                                                       \
   }

DEFINE_FUNCTION(SM, flexiblesusy::Two_scale, flexiblesusy::SM<flexiblesusy::Two_scale>)
DEFINE_FUNCTION(SMEFTHiggsTopDown, flexiblesusy::Shooting, standard_model::StandardModel<flexiblesusy::Shooting>)

BOOST_AUTO_TEST_CASE( test_Mh )
{
   const auto bu = calc_SM(slha_input);
   const auto td = calc_SMEFTHiggsTopDown(slha_input);

   BOOST_CHECK_CLOSE_FRACTION(bu.at(0), td.at(0), 1e-10);
   BOOST_CHECK_CLOSE_FRACTION(bu.at(1), td.at(1), 1e-10);
   BOOST_CHECK_CLOSE_FRACTION(bu.at(2), td.at(2), 1e-10);
   BOOST_CHECK_CLOSE_FRACTION(bu.at(3), td.at(3), 1e-10);
   BOOST_CHECK_CLOSE_FRACTION(bu.at(4), td.at(4), 1e-10);
   BOOST_CHECK_CLOSE_FRACTION(bu.at(5), td.at(5), 1e-10);
   BOOST_CHECK_CLOSE_FRACTION(bu.at(6), td.at(6), 1e-10);
   BOOST_CHECK_CLOSE_FRACTION(bu.at(7), td.at(7), 1e-10);
   BOOST_CHECK_CLOSE_FRACTION(bu.at(8), td.at(8), 1e-10);
}
