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
#define BOOST_TEST_MODULE test_NUHMSSMNoFVHimalayaEFTHiggs

#include <boost/test/unit_test.hpp>
#include <cstdlib>
#include "lowe.h"
#include "NUHMSSMNoFVHimalayaEFTHiggs_shooting_spectrum_generator.hpp"
#include "NUHMSSMNoFVHimalayaEFTHiggs_slha_io.hpp"
#include "wrappers.hpp"
using namespace flexiblesusy;

using output = std::array<double,7>;
using output_2loop = std::array<double,2>;
using output_3loop = std::array<double,1>;

double calc_lambda( 
   NUHMSSMNoFVHimalayaEFTHiggs_input_parameters& input,
   softsusy::QedQcd& qedqcd,
   Spectrum_generator_settings& settings)
{

   NUHMSSMNoFVHimalayaEFTHiggs_spectrum_generator<Shooting> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);

   const double Qmatch = 40000; 
   auto sm  = spectrum_generator.get_sm();
   sm.run_to(Qmatch);
   return sm.get_Lambdax();
}

double calc_lambda_3loop(
   NUHMSSMNoFVHimalayaEFTHiggs_input_parameters& input,
   softsusy::QedQcd& qedqcd,
   Spectrum_generator_settings& settings)
{

   NUHMSSMNoFVHimalayaEFTHiggs_spectrum_generator<Shooting> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);

   const double Qmatch = 50000; 
   auto sm  = spectrum_generator.get_sm();
   sm.run_to(Qmatch);
   return sm.get_Lambdax();
}

double calc_Mh( 
   NUHMSSMNoFVHimalayaEFTHiggs_input_parameters& input,
   softsusy::QedQcd& qedqcd,
   Spectrum_generator_settings& settings)
{

   NUHMSSMNoFVHimalayaEFTHiggs_spectrum_generator<Shooting> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);
   auto sm  = spectrum_generator.get_sm();

   double Q_pole = settings.get(Spectrum_generator_settings::eft_pole_mass_scale) != 0. ? settings.get(Spectrum_generator_settings::eft_pole_mass_scale) :  qedqcd.displayPoleMt();

   sm.run_to(Q_pole);
   sm.solve_ewsb();
   sm.calculate_Mhh_pole();
   return sm.get_physical().Mhh;
}

output edc_output( char const* const slha_input)
{
   output results = {0.,0.,0.,0.,0.,0.,0.};

   std::stringstream istream_case_1(slha_input);
   NUHMSSMNoFVHimalayaEFTHiggs_slha_io slha_io;
   slha_io.read_from_stream(istream_case_1);

   // extract the input parameters

   softsusy::QedQcd qedqcd;
   NUHMSSMNoFVHimalayaEFTHiggs_input_parameters input;
   Spectrum_generator_settings settings;

   try {
      slha_io.fill(settings);
      slha_io.fill(qedqcd);
      slha_io.fill(input);
   } catch (const Error& error) {
      BOOST_TEST_MESSAGE(error.what());
      BOOST_TEST(false);
   }

   results[1]=calc_lambda(input, qedqcd, settings);


   settings.set(Spectrum_generator_settings::eft_matching_loop_order_down, 1);
 
   results[0]=calc_Mh(input, qedqcd, settings);// /200.;
   results[2]=calc_lambda(input, qedqcd, settings);

/*   settings.set(Spectrum_generator_settings::eft_matching_loop_order_down, 2);
   
   settings.set(Spectrum_generator_settings::pole_mass_loop_order, 2);
   settings.set(Spectrum_generator_settings::higgs_2loop_correction_at_as, 1);
   settings.set(Spectrum_generator_settings::higgs_2loop_correction_ab_as, 1);

   results[3]=calc_lambda(input, qedqcd, settings);

   settings.set(Spectrum_generator_settings::higgs_2loop_correction_at_as, 0);
   results[4]=calc_lambda(input, qedqcd, settings);
   settings.set(Spectrum_generator_settings::higgs_2loop_correction_at_at, 1);
   results[5]=calc_lambda(input, qedqcd, settings);

   settings.set(Spectrum_generator_settings::higgs_2loop_correction_at_as, 1);
   settings.set(Spectrum_generator_settings::higgs_2loop_correction_ab_as, 1);
   settings.set(Spectrum_generator_settings::higgs_2loop_correction_at_at, 1);
   settings.set(Spectrum_generator_settings::higgs_2loop_correction_atau_atau, 1);
   settings.set(Spectrum_generator_settings::top_pole_qcd_corrections, 2);
   results[5]=calc_lambda(input, qedqcd, settings);

   settings.set(Spectrum_generator_settings::threshold_corrections_loop_order,3);
   settings.set(Spectrum_generator_settings::top_pole_qcd_corrections, 2);
//Have a look how implemented
   settings.set(Spectrum_generator_settings::eft_matching_loop_order_up, 2);
   settings.set(Spectrum_generator_settings::eft_matching_loop_order_down, 3);
   settings.set(Spectrum_generator_settings::pole_mass_loop_order, 3);
   settings.set(Spectrum_generator_settings::ewsb_loop_order, 3);
   settings.set(Spectrum_generator_settings::beta_loop_order, 4);
   results[6]=calc_lambda(input, qedqcd, settings);
*/
   return results;
}


output_2loop edc_output_2loop( char const* const slha_input)
{
   output_2loop results = {0.,0.};

   std::stringstream istream_case_1(slha_input);
   NUHMSSMNoFVHimalayaEFTHiggs_slha_io slha_io;
   slha_io.read_from_stream(istream_case_1);

   // extract the input parameters

   softsusy::QedQcd qedqcd;
   NUHMSSMNoFVHimalayaEFTHiggs_input_parameters input;
   Spectrum_generator_settings settings;

   try {
      slha_io.fill(settings);
      slha_io.fill(qedqcd);
      slha_io.fill(input);
   } catch (const Error& error) {
      BOOST_TEST_MESSAGE(error.what());
      BOOST_TEST(false);
   }

   settings.set(Spectrum_generator_settings::eft_matching_loop_order_down, 2); 
   settings.set(Spectrum_generator_settings::pole_mass_loop_order, 2);
   settings.set(Spectrum_generator_settings::higgs_2loop_correction_at_as, 1);
   results[0]=calc_Mh(input, qedqcd, settings);

   settings.set(Spectrum_generator_settings::higgs_2loop_correction_at_as, 0);
   settings.set(Spectrum_generator_settings::higgs_2loop_correction_at_at, 1);
   results[1]=calc_Mh(input, qedqcd, settings);

   return results;
}


output_3loop edc_output_3loop( char const* const slha_input)
{
   output_2loop results = {0.,0.};

   std::stringstream istream_case_1(slha_input);
   NUHMSSMNoFVHimalayaEFTHiggs_slha_io slha_io;
   slha_io.read_from_stream(istream_case_1);

   // extract the input parameters

   softsusy::QedQcd qedqcd;
   NUHMSSMNoFVHimalayaEFTHiggs_input_parameters input;
   Spectrum_generator_settings settings;

   try {
      slha_io.fill(settings);
      slha_io.fill(qedqcd);
      slha_io.fill(input);
   } catch (const Error& error) {
      BOOST_TEST_MESSAGE(error.what());
      BOOST_TEST(false);
   }


   settings.set(Spectrum_generator_settings::eft_matching_loop_order_down, 3); 
   settings.set(Spectrum_generator_settings::pole_mass_loop_order, 3);
   settings.set(Spectrum_generator_settings::higgs_2loop_correction_at_as, 1);
   settings.set(Spectrum_generator_settings::higgs_2loop_correction_at_at, 1);

   return { calc_lambda_3loop(input, qedqcd, settings) };
}


// scenario 1 degenrate case with vanishing stop mixing At = (1/10 + 0)*50000 = 5000
char const * const slha_input_case_1 = R"(
Block MODSEL                 # Select model
#   12    1000                # DRbar parameter output scale (GeV)
Block FlexibleSUSY
    0   1.000000000e-05      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = all, 1 = two_scale, 2 = semi_analytic)
    3   0                    # calculate SM pole masses
    4   1                    # pole mass loop order
    5   1                    # EWSB loop order
    6   3                    # beta-functions loop order
    7   1                    # threshold corrections loop order
    8   0                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   0                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   0                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   0                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   0                    # force output
   13   0                    # Top pole mass QCD corrections (0 = 1L, 1 = 2L, 2 = 3L)
   14   1.000000000e-11      # beta-function zero threshold
   15   0                    # calculate observables (a_muon, ...)
   16   0                    # force positive majorana masses
   17   0                    # pole mass renormalization scale (0 = SUSY scale)
   18   0                    # pole mass renormalization scale in the EFT (0 = min(SUSY scale, Mt))
   19   40000                    # EFT matching scale (0 = SUSY scale)
   20   1                    # EFT loop order for upwards matching
   21   0                    # EFT loop order for downwards matching
   22   0                    # EFT index of SM-like Higgs in the BSM model
   23   0                    # calculate BSM pole masses
   24   111111111            # individual threshold correction loop orders
   25   0                    # ren. scheme for Higgs 3L corrections (0 = DR, 1 = MDR)
   26   1                    # Higgs 3-loop corrections O(alpha_t alpha_s^2)
   27   0                    # Higgs 3-loop corrections O(alpha_b alpha_s^2)
   28   0                    # Higgs 3-loop corrections O(alpha_t^2 alpha_s)
   29   0                    # Higgs 3-loop corrections O(alpha_t^3)
   30   0                    # Higgs 4-loop corrections O(alpha_t alpha_s^3)
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
    0   50000                 # Ms
    1   50000                 # M1(MSUSY)
    2   50000                 # M2(MSUSY)
    3   50000                 # M3(MSUSY)
    4   50000                 # Mu(MSUSY)
    5   50000                 # mA(MSUSY)
   25   10                    # TanBeta(MSUSY)
Block MSQ2IN
  1  1     25.000E+08   # mq2(1,1)
  2  2     25.000E+08   # mq2(2,2)
  3  3     25.000E+08   # mq2(3,3)
Block MSE2IN
  1  1     25.000E+08   # me2(1,1)
  2  2     25.000E+08   # me2(2,2)
  3  3     25.000E+08   # me2(3,3)
Block MSL2IN
  1  1     25.000E+08   # ml2(1,1)
  2  2     25.000E+08   # ml2(2,2)
  3  3     25.000E+08   # ml2(3,3)
Block MSU2IN
  1  1     25.000E+08   # mu2(1,1)
  2  2     25.000E+08   # mu2(2,2)
  3  3     25.000E+08   # mu2(3,3)
Block MSD2IN
  1  1     25.000E+08   # md2(1,1)
  2  2     25.000E+08   # md2(2,2)
  3  3     25.000E+08   # md2(3,3)
Block AUIN
  1  1     0      # Au(1,1)
  2  2     0      # Au(2,2)
  3  3     5000   # Au(3,3)
Block ADIN
  1  1     0   # Ad(1,1)
  2  2     0   # Ad(2,2)
  3  3     0   # Ad(3,3)
Block AEIN
  1  1     0   # Ad(1,1)
  2  2     0   # Ad(2,2)
  3  3     0   # Ad(3,3)
)";

// scenario 2 degenrate case with high stop mixing At = (1/10 - sqrt(6))*50000 = -117474.5
char const * const slha_input_case_2 = R"(
Block MODSEL                 # Select model
#   12    1000                # DRbar parameter output scale (GeV)
Block FlexibleSUSY
    0   1.000000000e-05      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = all, 1 = two_scale, 2 = semi_analytic)
    3   0                    # calculate SM pole masses
    4   1                    # pole mass loop order
    5   1                    # EWSB loop order
    6   3                    # beta-functions loop order
    7   1                    # threshold corrections loop order
    8   0                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   0                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   0                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   0                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   0                    # force output
   13   0                    # Top pole mass QCD corrections (0 = 1L, 1 = 2L, 2 = 3L)
   14   1.000000000e-11      # beta-function zero threshold
   15   0                    # calculate observables (a_muon, ...)
   16   0                    # force positive majorana masses
   17   0                    # pole mass renormalization scale (0 = SUSY scale)
   18   0                    # pole mass renormalization scale in the EFT (0 = min(SUSY scale, Mt))
   19   40000                    # EFT matching scale (0 = SUSY scale)
   20   1                    # EFT loop order for upwards matching
   21   0                    # EFT loop order for downwards matching
   22   0                    # EFT index of SM-like Higgs in the BSM model
   23   0                    # calculate BSM pole masses
   24   111111111            # individual threshold correction loop orders
   25   0                    # ren. scheme for Higgs 3L corrections (0 = DR, 1 = MDR)
   26   1                    # Higgs 3-loop corrections O(alpha_t alpha_s^2)
   27   0                    # Higgs 3-loop corrections O(alpha_b alpha_s^2)
   28   0                    # Higgs 3-loop corrections O(alpha_t^2 alpha_s)
   29   0                    # Higgs 3-loop corrections O(alpha_t^3)
   30   0                    # Higgs 4-loop corrections O(alpha_t alpha_s^3)
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
    0   50000                 # Ms
    1   50000                 # M1(MSUSY)
    2   50000                 # M2(MSUSY)
    3   50000                 # M3(MSUSY)
    4   50000                 # Mu(MSUSY)
    5   50000                 # mA(MSUSY)
   25   10                    # TanBeta(MSUSY)
Block MSQ2IN
  1  1     25.000E+08   # mq2(1,1)
  2  2     25.000E+08   # mq2(2,2)
  3  3     25.000E+08   # mq2(3,3)
Block MSE2IN
  1  1     25.000E+08   # me2(1,1)
  2  2     25.000E+08   # me2(2,2)
  3  3     25.000E+08   # me2(3,3)
Block MSL2IN
  1  1     25.000E+08   # ml2(1,1)
  2  2     25.000E+08   # ml2(2,2)
  3  3     25.000E+08   # ml2(3,3)
Block MSU2IN
  1  1     25.000E+08   # mu2(1,1)
  2  2     25.000E+08   # mu2(2,2)
  3  3     25.000E+08   # mu2(3,3)
Block MSD2IN
  1  1     25.000E+08   # md2(1,1)
  2  2     25.000E+08   # md2(2,2)
  3  3     25.000E+08   # md2(3,3)
Block AUIN
  1  1     0            # Au(1,1)
  2  2     0            # Au(2,2)
  3  3     -117474.5    # Au(3,3)
Block ADIN
  1  1     0   # Ad(1,1)
  2  2     0   # Ad(2,2)
  3  3     0   # Ad(3,3)
Block AEIN
  1  1     0   # Ad(1,1)
  2  2     0   # Ad(2,2)
  3  3     0   # Ad(3,3)
)";

// scenario 3 non-degenrate case with vanishing stop mixing At = (1/10 - 0)*50000 = 5000
// mQ3 = MS/2
// mU3 = MS * 1.2
// M3  = MS/4
char const * const slha_input_case_3 = R"(
Block MODSEL                 # Select model
#   12    1000                # DRbar parameter output scale (GeV)
Block FlexibleSUSY
    0   1.000000000e-05      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = all, 1 = two_scale, 2 = semi_analytic)
    3   0                    # calculate SM pole masses
    4   1                    # pole mass loop order
    5   1                    # EWSB loop order
    6   3                    # beta-functions loop order
    7   1                    # threshold corrections loop order
    8   0                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   0                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   0                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   0                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   0                    # force output
   13   0                    # Top pole mass QCD corrections (0 = 1L, 1 = 2L, 2 = 3L)
   14   1.000000000e-11      # beta-function zero threshold
   15   0                    # calculate observables (a_muon, ...)
   16   0                    # force positive majorana masses
   17   0                    # pole mass renormalization scale (0 = SUSY scale)
   18   0                    # pole mass renormalization scale in the EFT (0 = min(SUSY scale, Mt))
   19   40000                    # EFT matching scale (0 = SUSY scale)
   20   1                    # EFT loop order for upwards matching
   21   0                    # EFT loop order for downwards matching
   22   0                    # EFT index of SM-like Higgs in the BSM model
   23   0                    # calculate BSM pole masses
   24   111111111            # individual threshold correction loop orders
   25   0                    # ren. scheme for Higgs 3L corrections (0 = DR, 1 = MDR)
   26   1                    # Higgs 3-loop corrections O(alpha_t alpha_s^2)
   27   0                    # Higgs 3-loop corrections O(alpha_b alpha_s^2)
   28   0                    # Higgs 3-loop corrections O(alpha_t^2 alpha_s)
   29   0                    # Higgs 3-loop corrections O(alpha_t^3)
   30   0                    # Higgs 4-loop corrections O(alpha_t alpha_s^3)
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
    0   50000                 # Ms
    1   50000                 # M1(MSUSY)
    2   50000                 # M2(MSUSY)
    3   12500                 # M3(MSUSY)
    4   50000                 # Mu(MSUSY)
    5   50000                 # mA(MSUSY)
   25   10                    # TanBeta(MSUSY)
Block MSQ2IN
  1  1     25.000E+08   # mq2(1,1)
  2  2     25.000E+08   # mq2(2,2)
  3  3     6.250E+08    # mq2(3,3)
Block MSE2IN
  1  1     25.000E+08   # me2(1,1)
  2  2     25.000E+08   # me2(2,2)
  3  3     25.000E+08   # me2(3,3)
Block MSL2IN
  1  1     25.000E+08   # ml2(1,1)
  2  2     25.000E+08   # ml2(2,2)
  3  3     25.000E+08   # ml2(3,3)
Block MSU2IN
  1  1     25.000E+08   # mu2(1,1)
  2  2     25.000E+08   # mu2(2,2)
  3  3     36.000E+08   # mu2(3,3)
Block MSD2IN
  1  1     25.000E+08   # md2(1,1)
  2  2     25.000E+08   # md2(2,2)
  3  3     25.000E+08   # md2(3,3)
Block AUIN
  1  1     0      # Au(1,1)
  2  2     0      # Au(2,2)
  3  3     5000   # Au(3,3)
Block ADIN
  1  1     0   # Ad(1,1)
  2  2     0   # Ad(2,2)
  3  3     0   # Ad(3,3)
Block AEIN
  1  1     0   # Ad(1,1)
  2  2     0   # Ad(2,2)
  3  3     0   # Ad(3,3)
)";

// scenario 4 non-degenrate case with vanishing stop mixing At = (1/10 - sqrt(6))*50000 = -117474.5
// mQ3 = MS/2
// mU3 = MS * 1.2
// M3  = MS/4
char const * const slha_input_case_4 = R"(
Block MODSEL                 # Select model
#   12    1000                # DRbar parameter output scale (GeV)
Block FlexibleSUSY
    0   1.000000000e-05      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = all, 1 = two_scale, 2 = semi_analytic)
    3   0                    # calculate SM pole masses
    4   1                    # pole mass loop order
    5   1                    # EWSB loop order
    6   3                    # beta-functions loop order
    7   1                    # threshold corrections loop order
    8   0                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   0                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   0                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   0                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   0                    # force output
   13   0                    # Top pole mass QCD corrections (0 = 1L, 1 = 2L, 2 = 3L)
   14   1.000000000e-11      # beta-function zero threshold
   15   0                    # calculate observables (a_muon, ...)
   16   0                    # force positive majorana masses
   17   0                    # pole mass renormalization scale (0 = SUSY scale)
   18   0                    # pole mass renormalization scale in the EFT (0 = min(SUSY scale, Mt))
   19   40000                    # EFT matching scale (0 = SUSY scale)
   20   1                    # EFT loop order for upwards matching
   21   0                    # EFT loop order for downwards matching
   22   0                    # EFT index of SM-like Higgs in the BSM model
   23   0                    # calculate BSM pole masses
   24   111111111            # individual threshold correction loop orders
   25   0                    # ren. scheme for Higgs 3L corrections (0 = DR, 1 = MDR)
   26   0                    # Higgs 3-loop corrections O(alpha_t alpha_s^2)
   27   0                    # Higgs 3-loop corrections O(alpha_b alpha_s^2)
   28   0                    # Higgs 3-loop corrections O(alpha_t^2 alpha_s)
   29   0                    # Higgs 3-loop corrections O(alpha_t^3)
   30   0                    # Higgs 4-loop corrections O(alpha_t alpha_s^3)
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
    0   50000                 # Ms
    1   50000                 # M1(MSUSY)
    2   50000                 # M2(MSUSY)
    3   12500                 # M3(MSUSY)
    4   50000                 # Mu(MSUSY)
    5   50000                 # mA(MSUSY)
   25   10                    # TanBeta(MSUSY)
Block MSQ2IN
  1  1     25.000E+08   # mq2(1,1)
  2  2     25.000E+08   # mq2(2,2)
  3  3     6.250E+08    # mq2(3,3)
Block MSE2IN
  1  1     25.000E+08   # me2(1,1)
  2  2     25.000E+08   # me2(2,2)
  3  3     25.000E+08   # me2(3,3)
Block MSL2IN
  1  1     25.000E+08   # ml2(1,1)
  2  2     25.000E+08   # ml2(2,2)
  3  3     25.000E+08   # ml2(3,3)
Block MSU2IN
  1  1     25.000E+08   # mu2(1,1)
  2  2     25.000E+08   # mu2(2,2)
  3  3     36.000E+08   # mu2(3,3)
Block MSD2IN
  1  1     25.000E+08   # md2(1,1)
  2  2     25.000E+08   # md2(2,2)
  3  3     25.000E+08   # md2(3,3)
Block AUIN
  1  1     0            # Au(1,1)
  2  2     0            # Au(2,2)
  3  3     -117474.5    # Au(3,3)
Block ADIN
  1  1     0   # Ad(1,1)
  2  2     0   # Ad(2,2)
  3  3     0   # Ad(3,3)
Block AEIN
  1  1     0   # Ad(1,1)
  2  2     0   # Ad(2,2)
  3  3     0   # Ad(3,3)
)";

// 2-loop scenarios: 
// a) At = (1/5 - 3)*50000 = -140000
// b) At = (1/5 - 2)*50000 = -90000
// c) At = (1/5 + 2)*50000 =  100000
// d) At = (1/5 + 3)*50000 =  160000

// scenario (1) 2-loop degenrate case with high stop mixing At = (1/5 - 3)*50000 = -140000
char const * const slha_input_case_1_2loop = R"(
Block MODSEL                 # Select model
#   12    1000                # DRbar parameter output scale (GeV)
Block FlexibleSUSY
    0   1.000000000e-05      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = all, 1 = two_scale, 2 = semi_analytic)
    3   0                    # calculate SM pole masses
    4   1                    # pole mass loop order
    5   1                    # EWSB loop order
    6   3                    # beta-functions loop order
    7   1                    # threshold corrections loop order
    8   0                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   0                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   0                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   0                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   0                    # force output
   13   0                    # Top pole mass QCD corrections (0 = 1L, 1 = 2L, 2 = 3L)
   14   1.000000000e-11      # beta-function zero threshold
   15   0                    # calculate observables (a_muon, ...)
   16   0                    # force positive majorana masses
   17   0                    # pole mass renormalization scale (0 = SUSY scale)
   18   0                    # pole mass renormalization scale in the EFT (0 = min(SUSY scale, Mt))
   19   50000                    # EFT matching scale (0 = SUSY scale)
   20   1                    # EFT loop order for upwards matching
   21   0                    # EFT loop order for downwards matching
   22   0                    # EFT index of SM-like Higgs in the BSM model
   23   0                    # calculate BSM pole masses
   24   111111111            # individual threshold correction loop orders
   25   0                    # ren. scheme for Higgs 3L corrections (0 = DR, 1 = MDR)
   26   1                    # Higgs 3-loop corrections O(alpha_t alpha_s^2)
   27   0                    # Higgs 3-loop corrections O(alpha_b alpha_s^2)
   28   0                    # Higgs 3-loop corrections O(alpha_t^2 alpha_s)
   29   0                    # Higgs 3-loop corrections O(alpha_t^3)
   30   0                    # Higgs 4-loop corrections O(alpha_t alpha_s^3)
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
    0   50000                 # Ms
    1   50000                 # M1(MSUSY)
    2   50000                 # M2(MSUSY)
    3   50000                 # M3(MSUSY)
    4   50000                 # Mu(MSUSY)
    5   50000                 # mA(MSUSY)
   25   5                    # TanBeta(MSUSY)
Block MSQ2IN
  1  1     25.000E+08   # mq2(1,1)
  2  2     25.000E+08   # mq2(2,2)
  3  3     25.000E+08   # mq2(3,3)
Block MSE2IN
  1  1     25.000E+08   # me2(1,1)
  2  2     25.000E+08   # me2(2,2)
  3  3     25.000E+08   # me2(3,3)
Block MSL2IN
  1  1     25.000E+08   # ml2(1,1)
  2  2     25.000E+08   # ml2(2,2)
  3  3     25.000E+08   # ml2(3,3)
Block MSU2IN
  1  1     25.000E+08   # mu2(1,1)
  2  2     25.000E+08   # mu2(2,2)
  3  3     25.000E+08   # mu2(3,3)
Block MSD2IN
  1  1     25.000E+08   # md2(1,1)
  2  2     25.000E+08   # md2(2,2)
  3  3     25.000E+08   # md2(3,3)
Block AUIN
  1  1     0            # Au(1,1)
  2  2     0            # Au(2,2)
  3  3     -140000    # Au(3,3)
Block ADIN
  1  1     0   # Ad(1,1)
  2  2     0   # Ad(2,2)
  3  3     0   # Ad(3,3)
Block AEIN
  1  1     0   # Ad(1,1)
  2  2     0   # Ad(2,2)
  3  3     0   # Ad(3,3)
)";

// scenario (2) 2-loop degenrate case with high stop mixing At = (1/5 - 2)*50000 approx -90100
char const * const slha_input_case_2_2loop = R"(
Block MODSEL                 # Select model
#   12    1000                # DRbar parameter output scale (GeV)
Block FlexibleSUSY
    0   1.000000000e-05      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = all, 1 = two_scale, 2 = semi_analytic)
    3   0                    # calculate SM pole masses
    4   1                    # pole mass loop order
    5   1                    # EWSB loop order
    6   3                    # beta-functions loop order
    7   1                    # threshold corrections loop order
    8   0                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   0                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   0                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   0                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   0                    # force output
   13   0                    # Top pole mass QCD corrections (0 = 1L, 1 = 2L, 2 = 3L)
   14   1.000000000e-11      # beta-function zero threshold
   15   0                    # calculate observables (a_muon, ...)
   16   0                    # force positive majorana masses
   17   0                    # pole mass renormalization scale (0 = SUSY scale)
   18   0                    # pole mass renormalization scale in the EFT (0 = min(SUSY scale, Mt))
   19   50000                    # EFT matching scale (0 = SUSY scale)
   20   1                    # EFT loop order for upwards matching
   21   0                    # EFT loop order for downwards matching
   22   0                    # EFT index of SM-like Higgs in the BSM model
   23   0                    # calculate BSM pole masses
   24   111111111            # individual threshold correction loop orders
   25   0                    # ren. scheme for Higgs 3L corrections (0 = DR, 1 = MDR)
   26   1                    # Higgs 3-loop corrections O(alpha_t alpha_s^2)
   27   0                    # Higgs 3-loop corrections O(alpha_b alpha_s^2)
   28   0                    # Higgs 3-loop corrections O(alpha_t^2 alpha_s)
   29   0                    # Higgs 3-loop corrections O(alpha_t^3)
   30   0                    # Higgs 4-loop corrections O(alpha_t alpha_s^3)
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
    0   50000                 # Ms
    1   50000                 # M1(MSUSY)
    2   50000                 # M2(MSUSY)
    3   50000                 # M3(MSUSY)
    4   50000                 # Mu(MSUSY)
    5   50000                 # mA(MSUSY)
   25   5                    # TanBeta(MSUSY)
Block MSQ2IN
  1  1     25.000E+08   # mq2(1,1)
  2  2     25.000E+08   # mq2(2,2)
  3  3     25.000E+08   # mq2(3,3)
Block MSE2IN
  1  1     25.000E+08   # me2(1,1)
  2  2     25.000E+08   # me2(2,2)
  3  3     25.000E+08   # me2(3,3)
Block MSL2IN
  1  1     25.000E+08   # ml2(1,1)
  2  2     25.000E+08   # ml2(2,2)
  3  3     25.000E+08   # ml2(3,3)
Block MSU2IN
  1  1     25.000E+08   # mu2(1,1)
  2  2     25.000E+08   # mu2(2,2)
  3  3     25.000E+08   # mu2(3,3)
Block MSD2IN
  1  1     25.000E+08   # md2(1,1)
  2  2     25.000E+08   # md2(2,2)
  3  3     25.000E+08   # md2(3,3)
Block AUIN
  1  1     0            # Au(1,1)
  2  2     0            # Au(2,2)
  3  3     -90100    # Au(3,3)
Block ADIN
  1  1     0   # Ad(1,1)
  2  2     0   # Ad(2,2)
  3  3     0   # Ad(3,3)
Block AEIN
  1  1     0   # Ad(1,1)
  2  2     0   # Ad(2,2)
  3  3     0   # Ad(3,3)
)";

// scenario (3) 2-loop degenrate case with high stop mixing At = (1/5 + 2)*50000 approx  110100
char const * const slha_input_case_3_2loop = R"(
Block MODSEL                 # Select model
#   12    1000                # DRbar parameter output scale (GeV)
Block FlexibleSUSY
    0   1.000000000e-05      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = all, 1 = two_scale, 2 = semi_analytic)
    3   0                    # calculate SM pole masses
    4   1                    # pole mass loop order
    5   1                    # EWSB loop order
    6   3                    # beta-functions loop order
    7   1                    # threshold corrections loop order
    8   0                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   0                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   0                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   0                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   0                    # force output
   13   0                    # Top pole mass QCD corrections (0 = 1L, 1 = 2L, 2 = 3L)
   14   1.000000000e-11      # beta-function zero threshold
   15   0                    # calculate observables (a_muon, ...)
   16   0                    # force positive majorana masses
   17   0                    # pole mass renormalization scale (0 = SUSY scale)
   18   0                    # pole mass renormalization scale in the EFT (0 = min(SUSY scale, Mt))
   19   50000                    # EFT matching scale (0 = SUSY scale)
   20   1                    # EFT loop order for upwards matching
   21   0                    # EFT loop order for downwards matching
   22   0                    # EFT index of SM-like Higgs in the BSM model
   23   0                    # calculate BSM pole masses
   24   111111111            # individual threshold correction loop orders
   25   0                    # ren. scheme for Higgs 3L corrections (0 = DR, 1 = MDR)
   26   1                    # Higgs 3-loop corrections O(alpha_t alpha_s^2)
   27   0                    # Higgs 3-loop corrections O(alpha_b alpha_s^2)
   28   0                    # Higgs 3-loop corrections O(alpha_t^2 alpha_s)
   29   0                    # Higgs 3-loop corrections O(alpha_t^3)
   30   0                    # Higgs 4-loop corrections O(alpha_t alpha_s^3)
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
    0   50000                 # Ms
    1   50000                 # M1(MSUSY)
    2   50000                 # M2(MSUSY)
    3   50000                 # M3(MSUSY)
    4   50000                 # Mu(MSUSY)
    5   50000                 # mA(MSUSY)
   25   5                    # TanBeta(MSUSY)
Block MSQ2IN
  1  1     25.000E+08   # mq2(1,1)
  2  2     25.000E+08   # mq2(2,2)
  3  3     25.000E+08   # mq2(3,3)
Block MSE2IN
  1  1     25.000E+08   # me2(1,1)
  2  2     25.000E+08   # me2(2,2)
  3  3     25.000E+08   # me2(3,3)
Block MSL2IN
  1  1     25.000E+08   # ml2(1,1)
  2  2     25.000E+08   # ml2(2,2)
  3  3     25.000E+08   # ml2(3,3)
Block MSU2IN
  1  1     25.000E+08   # mu2(1,1)
  2  2     25.000E+08   # mu2(2,2)
  3  3     25.000E+08   # mu2(3,3)
Block MSD2IN
  1  1     25.000E+08   # md2(1,1)
  2  2     25.000E+08   # md2(2,2)
  3  3     25.000E+08   # md2(3,3)
Block AUIN
  1  1     0            # Au(1,1)
  2  2     0            # Au(2,2)
  3  3     110100    # Au(3,3)
Block ADIN
  1  1     0   # Ad(1,1)
  2  2     0   # Ad(2,2)
  3  3     0   # Ad(3,3)
Block AEIN
  1  1     0   # Ad(1,1)
  2  2     0   # Ad(2,2)
  3  3     0   # Ad(3,3)
)";

// scenario (4) 2-loop degenrate case with high stop mixing At = (1/5 + 3)*50000 =  160000
char const * const slha_input_case_4_2loop = R"(
Block MODSEL                 # Select model
#   12    1000                # DRbar parameter output scale (GeV)
Block FlexibleSUSY
    0   1.000000000e-05      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = all, 1 = two_scale, 2 = semi_analytic)
    3   0                    # calculate SM pole masses
    4   1                    # pole mass loop order
    5   1                    # EWSB loop order
    6   3                    # beta-functions loop order
    7   1                    # threshold corrections loop order
    8   0                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   0                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   0                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   0                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   0                    # force output
   13   0                    # Top pole mass QCD corrections (0 = 1L, 1 = 2L, 2 = 3L)
   14   1.000000000e-11      # beta-function zero threshold
   15   0                    # calculate observables (a_muon, ...)
   16   0                    # force positive majorana masses
   17   0                    # pole mass renormalization scale (0 = SUSY scale)
   18   0                    # pole mass renormalization scale in the EFT (0 = min(SUSY scale, Mt))
   19   50000                    # EFT matching scale (0 = SUSY scale)
   20   1                    # EFT loop order for upwards matching
   21   0                    # EFT loop order for downwards matching
   22   0                    # EFT index of SM-like Higgs in the BSM model
   23   0                    # calculate BSM pole masses
   24   111111111            # individual threshold correction loop orders
   25   0                    # ren. scheme for Higgs 3L corrections (0 = DR, 1 = MDR)
   26   1                    # Higgs 3-loop corrections O(alpha_t alpha_s^2)
   27   0                    # Higgs 3-loop corrections O(alpha_b alpha_s^2)
   28   0                    # Higgs 3-loop corrections O(alpha_t^2 alpha_s)
   29   0                    # Higgs 3-loop corrections O(alpha_t^3)
   30   0                    # Higgs 4-loop corrections O(alpha_t alpha_s^3)
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
    0   50000                 # Ms
    1   50000                 # M1(MSUSY)
    2   50000                 # M2(MSUSY)
    3   50000                 # M3(MSUSY)
    4   50000                 # Mu(MSUSY)
    5   50000                 # mA(MSUSY)
   25   5                    # TanBeta(MSUSY)
Block MSQ2IN
  1  1     25.000E+08   # mq2(1,1)
  2  2     25.000E+08   # mq2(2,2)
  3  3     25.000E+08   # mq2(3,3)
Block MSE2IN
  1  1     25.000E+08   # me2(1,1)
  2  2     25.000E+08   # me2(2,2)
  3  3     25.000E+08   # me2(3,3)
Block MSL2IN
  1  1     25.000E+08   # ml2(1,1)
  2  2     25.000E+08   # ml2(2,2)
  3  3     25.000E+08   # ml2(3,3)
Block MSU2IN
  1  1     25.000E+08   # mu2(1,1)
  2  2     25.000E+08   # mu2(2,2)
  3  3     25.000E+08   # mu2(3,3)
Block MSD2IN
  1  1     25.000E+08   # md2(1,1)
  2  2     25.000E+08   # md2(2,2)
  3  3     25.000E+08   # md2(3,3)
Block AUIN
  1  1     0            # Au(1,1)
  2  2     0            # Au(2,2)
  3  3     160000       # Au(3,3)
Block ADIN
  1  1     0   # Ad(1,1)
  2  2     0   # Ad(2,2)
  3  3     0   # Ad(3,3)
Block AEIN
  1  1     0   # Ad(1,1)
  2  2     0   # Ad(2,2)
  3  3     0   # Ad(3,3)
)";


BOOST_AUTO_TEST_CASE( test_top_down_EFTHiggs_1loop )
{
   const struct Data {
      char const * const slha = nullptr;
      output out{};
      double eps{0.0};
   } data[] = {
      {slha_input_case_1, {129.024202, 0.125588145, 0.125768663}, 2e-4},
      {slha_input_case_2, {132.709997, 0.125588071, 0.147703015}, 8e-5},
      {slha_input_case_3, {128.315047, 0.12558815 , 0.121657248}, 8e-5},
      {slha_input_case_4, {130.665551, 0.125588152, 0.135417609}, 8e-5},
   };

   for (const auto& d: data) {
      const auto out = edc_output(d.slha);
      for (int i = 0; i < d.out.size(); i++) {
         BOOST_CHECK_CLOSE_FRACTION(out[i], d.out[i], d.eps);
      }
   }
}


BOOST_AUTO_TEST_CASE( test_top_down_EFTHiggs_2loop )
{
   const struct Data {
      char const * const slha = nullptr;
      output_2loop out{};
      double eps{0.0};
   } data[] = {
      {slha_input_case_1_2loop, {128.16151975160173, 127.9817088038441 }, 5e-4},
      {slha_input_case_2_2loop, {130.42485668023875, 130.15819922320634}, 5e-4},
      {slha_input_case_3_2loop, {130.86500014892152, 130.66749808261707}, 5e-4},
      {slha_input_case_4_2loop, {128.87146686565686, 128.53192456935258}, 5e-4},
   };

   for (const auto& d: data) {
      const auto out = edc_output_2loop(d.slha);
      for (int i = 0; i < d.out.size(); i++) {
         BOOST_CHECK_CLOSE_FRACTION(out[i], d.out[i], d.eps);
      }
   }
}


BOOST_AUTO_TEST_CASE( test_top_down_EFTHiggs_3loop )
{
   const struct Data {
      char const * const slha = nullptr;
      output_3loop out{};
      double eps{0.0};
   } data[] = {
      {slha_input_case_1_2loop, {0.11685905941993063}, 5e-5},
      {slha_input_case_4_2loop, {0.12118666568388101}, 5e-5},
   };

   for (const auto& d: data) {
      const auto out = edc_output_3loop(d.slha);
      for (int i = 0; i < d.out.size(); i++) {
         BOOST_CHECK_CLOSE_FRACTION(out[i], d.out[i], d.eps);
      }
   }
}
