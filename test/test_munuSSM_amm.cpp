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
#define BOOST_TEST_MODULE test_munuSSM_amm

#include <boost/test/unit_test.hpp>
#include <cstdlib>

#include "lowe.h"

#include "test_munuSSM.hpp"

#include "wrappers.hpp"
#include "munuSSM_amm.hpp"
#include "munuSSM_slha_io.hpp"
#include "cxx_qft/munuSSM_qft.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_amu )
{
   char const * const slha_input = R"(
Block MODSEL
Block FlexibleSUSY
    0   1.000000000e-04
    1   0
    2   0
    3   0
    4   1
    5   1
    6   2
    7   2
    8   1
    9   1
   10   1
   11   1
   12   0
   13   1
   14   1.000000000e-16
   15   1
   16   0
   17   0
   18   0
   19   0
   20   2
   21   1
   22   0
   23   1
   24   123111321
   25   0
   26   1
   27   1
   28   1
   29   1
   30   1
   31   0
Block SMINPUTS
    1   1.279340000e+02
    2   1.166378700e-05
    3   1.176000000e-01
    4   9.118760000e+01
    5   4.200000000e+00
    6   1.733000000e+02
    7   1.777000000e+00
    8   0.000000000e+00
    9   80.404
   11   5.109989020e-04
   12   0.000000000e+00
   13   1.056583570e-01
   14   0.000000000e+00
   21   4.750000000e-03
   22   2.400000000e-03
   23   1.040000000e-01
   24   1.270000000e+00
Block MINPAR
    1   500
    2   300
    3   10
    5   -10
Block EXTPAR
   61   0.2
   62   0.1
   65   10
   200  10
   201  10
   202  10
Block YVIN
    1   0.01
    2   0.01
    3   0.01
)";

   std::stringstream istr(slha_input);

   munuSSM_slha_io slha_io;
   slha_io.read_from_stream(istr);

   // extract the input parameters
   softsusy::QedQcd qedqcd;
   munuSSM_input_parameters input;
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

   munuSSM_slha m = setup_munuSSM(input, qedqcd, settings);

   using munuSSM_cxx_diagrams::fields::Cha;

   // 1L + 2L QED
   settings.set(Spectrum_generator_settings::calculate_amm, 1.5);

   double ae = munuSSM_amm::calculate_amm<Cha>(m, qedqcd, settings, 0);
   BOOST_CHECK_CLOSE_FRACTION(ae, 1.9398670311763684e-14, 1e-6);

   double amu = munuSSM_amm::calculate_amm<Cha>(m, qedqcd, settings, 1);
   BOOST_CHECK_CLOSE_FRACTION(amu, 3.4829855830138246e-09, 1e-6);

   double atau = munuSSM_amm::calculate_amm<Cha>(m, qedqcd, settings, 2);
   BOOST_CHECK_CLOSE_FRACTION(atau, 1.7900498254975007e-05, 1e-6);

   // 1L + 2L QED + Barr-Zee
   settings.set(Spectrum_generator_settings::calculate_amm, 2.0);

   ae = munuSSM_amm::calculate_amm<Cha>(m, qedqcd, settings, 0);
   BOOST_CHECK_CLOSE_FRACTION(ae, 1.588748576865684e-11, 1e-6);

   amu = munuSSM_amm::calculate_amm<Cha>(m, qedqcd, settings, 1);
   BOOST_CHECK_CLOSE_FRACTION(amu, 4.2677731488826819e-09, 1e-6);

   atau = munuSSM_amm::calculate_amm<Cha>(m, qedqcd, settings, 2);
   BOOST_CHECK_CLOSE_FRACTION(atau, 1.717142461318717e-05, 1e-6);
}
