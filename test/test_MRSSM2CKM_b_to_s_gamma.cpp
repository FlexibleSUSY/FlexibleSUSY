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
#define BOOST_TEST_MODULE test_MRSSM2CKM_b_to_s_gamma

#include <boost/test/unit_test.hpp>
#include <cstdlib>
#include <complex>

#include "lowe.h"
#include "test_complex_equality.hpp"

#include "MRSSM2CKM_two_scale_spectrum_generator.hpp"

#include "wrappers.hpp"
#include "MRSSM2CKM_b_to_s_gamma.hpp"
#include "MRSSM2CKM_slha_io.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_MRSSM2CKM_b_to_s_gamma, * boost::unit_test::tolerance(1e-6) )
{
   char const * const slha_input = R"(
Block MODSEL                 # Select model
    6   1                    # flavour violation
Block FlexibleSUSY
    0   1.000000000e-04      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = all, 1 = two_scale, 2 = semi_analytic)
    3   0                    # calculate SM pole masses
    4   2                    # pole mass loop order
    5   2                    # EWSB loop order
    6   3                    # beta-functions loop order
    7   2                    # threshold corrections loop order
    8   1                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   1                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   1                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   1                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   0                    # force output
   13   1                    # Top pole mass QCD corrections (0 = 1L, 1 = 2L, 2 = 3L)
   14   1.000000000e-11      # beta-function zero threshold
   15   1                    # calculate observables (a_muon, ...)
   16   0                    # force positive majorana masses
   17   160                    # pole mass renormalization scale (0 = SUSY scale)
   18   0                    # pole mass renormalization scale in the EFT (0 = min(SUSY scale, Mt))
   19   0                    # EFT matching scale (0 = SUSY scale)
   20   2                    # EFT loop order for upwards matching
   21   1                    # EFT loop order for downwards matching
   22   0                    # EFT index of SM-like Higgs in the BSM model
   23   1                    # calculate BSM pole masses
   24   123111321            # individual threshold correction loop orders
   25   0                    # ren. scheme for Higgs 3L corrections (0 = DR, 1 = MDR)
   26   1                    # Higgs 3-loop corrections O(alpha_t alpha_s^2)
   27   1                    # Higgs 3-loop corrections O(alpha_b alpha_s^2)
   28   1                    # Higgs 3-loop corrections O(alpha_t^2 alpha_s)
   29   1                    # Higgs 3-loop corrections O(alpha_t^3)
   30   1                    # Higgs 4-loop corrections O(alpha_t alpha_s^3)
Block FlexibleSUSYInput
    0   0.00729735           # alpha_em(0)
    1   125.09               # Mh pole
Block SMINPUTS               # Standard Model inputs
    1   1.279440000e+02      # alpha^(-1) SM MSbar(MZ)
    2   1.166378700e-05      # G_Fermi
    3   1.184000000e-01      # alpha_s(MZ) SM MSbar
    4   9.118760000e+01      # MZ(pole)
    5   4.180000000e+00      # mb(mb) SM MSbar
    6   1.733400000e+02      # mtop(pole)
    7   1.777000000e+00      # mtau(pole)
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
Block VCKMIN                 # CKM matrix input (Wolfenstein parameters)
    1   2.272000000e-01      # lambda(MZ) SM DR-bar
    2   8.180000000e-01      # A(MZ) SM DR-bar
    3   2.210000000e-01      # rhobar(MZ) SM DR-bar
    4   3.400000000e-01      # etabar(MZ) SM DR-bar
Block MINPAR
   3   20         # TanBeta
Block EXTPAR
   0   1000      # SUSY scale, Ms
Block HMIXIN
   303   -0.1    # LamTDInput
   304   -0.1    # LamTUInput
   301   -0.1    # LamSDInput
   302   -0.1    # LamSUInput
   201   1000    # MuDInput
   202   1000    # MuUInput
   101   1000    # BMuInput
Block MSQ2IN
   1   1   2e6   # mq2Input(1,1)
   1   2   0     # mq2Input(1,2)
   1   3   0     # mq2Input(1,3)
   2   1   0     # mq2Input(2,1)
   2   2   1e6   # mq2Input(2,2)
   2   3   1e4     # mq2Input(2,3)
   3   1   0     # mq2Input(3,1)
   3   2   0     # mq2Input(3,2)
   3   3   1e6   # mq2Input(3,3)
Block MSL2IN
   1   1   1e6   # ml2Input(1,1)
   1   2   0     # ml2Input(1,2)
   1   3   0     # ml2Input(1,3)
   2   1   0     # ml2Input(2,1)
   2   2   1e6   # ml2Input(2,2)
   2   3   0     # ml2Input(2,3)
   3   1   0     # ml2Input(3,1)
   3   2   0     # ml2Input(3,2)
   3   3   1e6   # ml2Input(3,3)
Block MSD2IN
   1   1   2e6   # md2Input(1,1)
   1   2   0     # md2Input(1,2)
   1   3   0     # md2Input(1,3)
   2   1   0     # md2Input(2,1)
   2   2   1e6   # md2Input(2,2)
   2   3   1e4     # md2Input(2,3)
   3   1   0     # md2Input(3,1)
   3   2   0     # md2Input(3,2)
   3   3   1e6   # md2Input(3,3)
Block MSU2IN
   1   1   2e6   # mu2Input(1,1)
   1   2   0     # mu2Input(1,2)
   1   3   0     # mu2Input(1,3)
   2   1   0     # mu2Input(2,1)
   2   2   1e6   # mu2Input(2,2)
   2   3   1e4     # mu2Input(2,3)
   3   1   0     # mu2Input(3,1)
   3   2   0     # mu2Input(3,2)
   3   3   1e6   # mu2Input(3,3)
Block MSE2IN
   1   1   1e6   # me2Input(1,1)
   1   2   0     # me2Input(1,2)
   1   3   0     # me2Input(1,3)
   2   1   0     # me2Input(2,1)
   2   2   1e6   # me2Input(2,2)
   2   3   0     # me2Input(2,3)
   3   1   0     # me2Input(3,1)
   3   2   0     # me2Input(3,2)
   3   3   1e6   # me2Input(3,3)
Block NMSSMRUNIN
    10   1e6     # mS2Input
Block MSOFTIN
   110   1e6     # mT2Input
   111   1e6     # moc2Input
    50   1e6     # mRd2Input
    51   1e6     # mRu2Input
   300   1000    # MDBSInput
   301   1000    # MDWBTInput
   302   1000    # MDGocInput
)";

   std::stringstream istr(slha_input);

   MRSSM2CKM_slha_io slha_io;
   slha_io.read_from_stream(istr);

   // extract the input parameters
   softsusy::QedQcd qedqcd;
   MRSSM2CKM_input_parameters input;
   Spectrum_generator_settings settings;

   try {
      slha_io.fill(settings);
      slha_io.fill(qedqcd);
      slha_io.fill(input);
   } catch (const Error& error) {
      BOOST_TEST_MESSAGE(error.what());
      BOOST_TEST(false);
   }

   settings.set(Spectrum_generator_settings::calculate_sm_masses, 0);
   settings.set(Spectrum_generator_settings::calculate_bsm_masses, 0);
   MRSSM2CKM_spectrum_generator<Two_scale> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);

   auto model = std::get<0>(spectrum_generator.get_models_slha());

   const auto reference_value = MRSSM2CKM_b_to_s_gamma::calculate_b_to_s_gamma(model, qedqcd);
   constexpr std::complex<double> C7NP {-0.20234454161004062, -0.00369650084326421};
   constexpr std::complex<double> C7pNP {-0.00435452005573289, -0.00007954989509452};
   constexpr std::complex<double> C8NP {-0.17125493908345221, -0.00312854511269766};
   constexpr std::complex<double> C8pNP {-0.00379732860057805, -0.00006937092674948};

   TEST_COMPLEX_EQUALITY(C7NP, reference_value[0]);
   TEST_COMPLEX_EQUALITY(C7pNP, reference_value[1]);
   TEST_COMPLEX_EQUALITY(C8NP, reference_value[2]);
   TEST_COMPLEX_EQUALITY(C8pNP, reference_value[3]);
}
