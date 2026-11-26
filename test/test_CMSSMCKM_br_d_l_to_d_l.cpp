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
#define BOOST_TEST_MODULE test_CMSSMCKM_br_d_l_to_d_l

#include <boost/test/unit_test.hpp>

#include "CMSSMCKM_b_to_s_gamma.hpp"
#include "CMSSMCKM_slha_io.hpp"
#include "observables/CMSSMCKM_br_d_l_to_d_l.hpp"
#include "test_CMSSMCKM.hpp"

#include "json.hpp"
#include "loop_libraries/loop_library.hpp"
#include "lowe.h"
#include "physical_input.hpp"

#include <array>
#include <complex>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

using namespace flexiblesusy;

namespace
{

// CMSSM-like point with m0 = 500 GeV, m1/2 = 1 TeV, A0 = 0, tan(beta) = 40,
// mu > 0, where the Higgs penguins enhance BR(B_s -> mu mu)
char const* const slha_input = R"(
Block MODSEL                 # Select model
    6   1                    # flavour violation
#   12    1000                # DRbar parameter output scale (GeV)
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
   13   1                    # Top pole mass QCD corrections (0 = 1L, 1 = 2L, 2 = 3L, 3 = 4L)
   14   1.000000000e-11      # beta-function zero threshold
   15   1                    # calculate all observables
   16   0                    # force positive majorana masses
   17   0                    # pole mass renormalization scale (0 = SUSY scale)
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
   31   1                    # loop library (0 = softsusy)
   32   2                    # loop level to calculate AMM
Block SMINPUTS               # Standard Model inputs
    1   1.279340000e+02      # alpha^(-1) SM MSbar(MZ)
    2   1.166378700e-05      # G_Fermi
    3   1.176000000e-01      # alpha_s(MZ) SM MSbar
    4   9.118760000e+01      # MZ(pole)
    5   4.200000000e+00      # mb(mb) SM MSbar
    6   1.733000000e+02      # mtop(pole)
    7   1.777000000e+00      # mtau(pole)
    8   0.000000000e+00      # mnu3(pole)
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
# Block UPMNSIN                # PMNS matrix input
#     1   5.837630000e-01      # theta_12
#     2   7.695840000e-01      # theta_23
#     3   1.549480000e-01      # theta_13
#     4   0.000000000e+00      # delta
#     5   0.000000000e+00      # alpha_1
#     6   0.000000000e+00      # alpha_2
Block MINPAR                 # Input parameters
    3   40                  # TanBeta
    4   1                    # SignMu
Block EXTPAR
   21   2.500000e+05                # mHd2
   22   2.500000e+05                # mHu2
Block TEIN
    1   1   0                # T[Ye][1,1]
    2   2   0                # T[Ye][2,2]
    3   3   0                # T[Ye][3,3]
Block TDIN
    1   1   0                # T[Yd][1,1]
    2   2   0                # T[Yd][2,2]
    3   3   0                # T[Yd][3,3]
Block TUIN
    1   1   0                # T[Yu][1,1]
    2   2   0                # T[Yu][2,2]
    3   3   0                # T[Yu][3,3]
Block MSQ2IN
  1  1     2.500000e+05          # mq2(1,1)
  1  2     0.00000000E+00   # mq2(1,2)
  1  3     0.00000000E+00   # mq2(1,3)
  2  1     0.00000000E+00   # mq2(2,1)
  2  2     2.500000e+05          # mq2(2,2)
  2  3     0.00000000E+00   # mq2(2,3)
  3  1     0.00000000E+00   # mq2(3,1)
  3  2     0.00000000E+00   # mq2(3,2)
  3  3     2.500000e+05          # mq2(3,3)
Block MSE2IN
  1  1     2.500000e+05          # me2(1,1)
  1  2     0.00000000E+00   # me2(1,2)
  1  3     0.00000000E+00   # me2(1,3)
  2  1     0.00000000E+00   # me2(2,1)
  2  2     2.500000e+05          # me2(2,2)
  2  3     0.00000000E+00   # me2(2,3)
  3  1     0.00000000E+00   # me2(3,1)
  3  2     0.00000000E+00   # me2(3,2)
  3  3     2.500000e+05          # me2(3,3)
Block MSL2IN
  1  1     2.500000e+05          # ml2(1,1)
  1  2     0.00000000E+00   # ml2(1,2)
  1  3     0.00000000E+00   # ml2(1,3)
  2  1     0.00000000E+00   # ml2(2,1)
  2  2     2.500000e+05          # ml2(2,2)
  2  3     0.00000000E+00   # ml2(2,3)
  3  1     0.00000000E+00   # ml2(3,1)
  3  2     0.00000000E+00   # ml2(3,2)
  3  3     2.500000e+05          # ml2(3,3)
Block MSU2IN
  1  1     2.500000e+05          # mu2(1,1)
  1  2     0.00000000E+00   # mu2(1,2)
  1  3     0.00000000E+00   # mu2(1,3)
  2  1     0.00000000E+00   # mu2(2,1)
  2  2     2.500000e+05          # mu2(2,2)
  2  3     0.00000000E+00   # mu2(2,3)
  3  1     0.00000000E+00   # mu2(3,1)
  3  2     0.00000000E+00   # mu2(3,2)
  3  3     2.500000e+05          # mu2(3,3)
Block MSD2IN
  1  1     2.500000e+05          # md2(1,1)
  1  2     0.00000000E+00   # md2(1,2)
  1  3     0.00000000E+00   # md2(1,3)
  2  1     0.00000000E+00   # md2(2,1)
  2  2     2.500000e+05          # md2(2,2)
  2  3     0.00000000E+00   # md2(2,3)
  3  1     0.00000000E+00   # md2(3,1)
  3  2     0.00000000E+00   # md2(3,2)
  3  3     2.500000e+05          # md2(3,3)
Block MSOFTIN
     1     1000   # MassB
     2     1000   # MassWB
     3     1000   # MassG
)";

std::complex<double> read_wc(const nlohmann::json& j, const std::string& name)
{
   const auto& v = j.at("values").at(name);
   return {v.at("Re").get<double>(), v.at("Im").get<double>()};
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(test_b_to_s_mu_mu)
{
   std::stringstream istr(slha_input);

   CMSSMCKM_slha_io slha_io;
   slha_io.read_from_stream(istr);

   softsusy::QedQcd qedqcd;
   CMSSMCKM_input_parameters input;
   Spectrum_generator_settings settings;

   slha_io.fill(settings);
   slha_io.fill(qedqcd);
   slha_io.fill(input);

   settings.set(Spectrum_generator_settings::calculate_observables, 0);
   Loop_library::set(1);

   const CMSSMCKM_slha model = setup_CMSSMCKM(input, qedqcd, settings);

   (void)std::remove("WC_CMSSMCKM.json");

   // b -> s mu mu: indices of b, s, mu
   const auto res =
      CMSSMCKM_br_d_l_to_d_l::calculate_FdFe_FdFe_ListVectorsScalarsBoxes_1(
         2, 1, 1, model, qedqcd);

   // regression test of the Wilson coefficients in the flavio basis
   // (C7, C7', C8, C8', C9, C9', C10, C10', CS, CS', CP, CP'),
   // values of the SLHA output of run_CMSSMCKM.x (9 digits)
   const std::array<std::complex<double>, 12> expected = {
      std::complex<double>{1.71541231E-02, 9.70417589E-09},   // C7
      std::complex<double>{4.25103704E-04, 2.46515596E-10},   // C7'
      std::complex<double>{-2.38395981E-03, 2.19517782E-08},  // C8
      std::complex<double>{-6.89512042E-05, 5.54826803E-10},  // C8'
      std::complex<double>{1.32248449E-03, -3.17774492E-10},  // C9
      std::complex<double>{2.70644760E-04, -5.68821819E-06},  // C9'
      std::complex<double>{3.00639160E-03, -8.98078050E-10},  // C10
      std::complex<double>{6.94071385E-04, -8.04746989E-10},  // C10'
      std::complex<double>{4.00966421E-03, 5.78795006E-11},   // CS
      std::complex<double>{1.04453981E-04, 6.84333959E-13},   // CS'
      std::complex<double>{-3.98574119E-03, -6.17669772E-11}, // CP
      std::complex<double>{1.03850889E-04, 7.83180556E-13},   // CP'
   };
   for (int i = 0; i < 12; i++) {
      BOOST_CHECK_CLOSE(res[i + 1].real(), expected[i].real(), 1e-5);
      BOOST_CHECK_CLOSE(res[i + 1].imag(), expected[i].imag(), 1e-5);
   }

   // the json file contains the same semi-leptonic coefficients
   std::ifstream ifs("WC_CMSSMCKM.json");
   const auto j = nlohmann::json::parse(ifs);
   BOOST_CHECK_EQUAL(j.at("eft"), "WET");
   BOOST_CHECK_EQUAL(j.at("basis"), "flavio");
   BOOST_CHECK_CLOSE(j.at("scale").get<double>(), model.get_scale(), 1e-10);

   const std::array<const char*, 8> names = {
      "C9_bsmumu", "C9p_bsmumu", "C10_bsmumu", "C10p_bsmumu",
      "CS_bsmumu", "CSp_bsmumu", "CP_bsmumu",  "CPp_bsmumu"};
   for (int i = 0; i < 8; i++) {
      const auto c = read_wc(j, names[i]);
      BOOST_CHECK_CLOSE(c.real(), res[i + 5].real(), 1e-10);
      BOOST_CHECK_CLOSE(c.imag(), res[i + 5].imag(), 1e-10);
   }
   const auto CP = read_wc(j, "CP_bsmumu");

   // C7 and C8 agree with the b -> s gamma observable, which normalizes
   // with g2^2/(8 MW^2) instead of G_F/sqrt(2)
   const auto bsg =
      CMSSMCKM_b_to_s_gamma::calculate_b_to_s_gamma(model, qedqcd);
   BOOST_CHECK_CLOSE(res[1].real(), bsg[0].real(), 5.);
   BOOST_CHECK_CLOSE(res[3].real(), bsg[2].real(), 5.);

   // The Higgs penguins interfere constructively with the Standard Model
   // (C10_SM ~ -4.2 in flavio), which enhances BR(B_s -> mu mu) by about
   // 30% (also SARAH/SPheno's FlavorKit gives an enhancement here).
   BOOST_CHECK_LT(CP.real(), -1e-3);

   // The CP-even Higgs boson with the pole mass closest to
   // FlexibleSUSYInput[1] counts as SM Higgs: diagrams with it and only SM
   // fields are removed. With the default (125.09 GeV) that is h, so the
   // pseudoscalar coefficient CP (no h exchange) and CS stay as above.
   // Declaring the heavy H as SM-like removes its diagrams with SM loops
   // (H exchange with W/top loops), which changes CS, but not CP.
   const Physical_input pi_h;
   const auto res_h =
      CMSSMCKM_br_d_l_to_d_l::calculate_FdFe_FdFe_ListVectorsScalarsBoxes_1(
         2, 1, 1, model, qedqcd, pi_h);
   BOOST_CHECK_EQUAL(res_h[9], res[9]);
   BOOST_CHECK_EQUAL(res_h[11], res[11]);

   Physical_input pi_H;
   pi_H.set(Physical_input::mh_pole, model.get_physical().Mhh(1));
   const auto res_H =
      CMSSMCKM_br_d_l_to_d_l::calculate_FdFe_FdFe_ListVectorsScalarsBoxes_1(
         2, 1, 1, model, qedqcd, pi_H);
   BOOST_CHECK_GT(std::abs(res_H[9] - res[9]), 1e-2 * std::abs(res[9]));
   BOOST_CHECK_CLOSE(res_H[11].real(), res[11].real(), 1e-10);

   (void)std::remove("WC_CMSSMCKM.json");
}
