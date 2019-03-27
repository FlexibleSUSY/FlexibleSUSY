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
#define BOOST_TEST_MODULE test_munuSSM_gmm2

#include <boost/test/unit_test.hpp>
#include <cstdlib>

#include "test_munuSSM.hpp"

#include "wrappers.hpp"
#include "munuSSM_a_muon.hpp"
#include "munuSSM_slha_io.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_amu )
{
   char const * const slha_input = "\n\
Block MODSEL\n\
Block FlexibleSUSY\n\
    0   1.000000000e-04\n\
    1   0\n\
    2   0\n\
    3   0\n\
    4   1\n\
    5   1\n\
    6   2\n\
    7   2\n\
    8   1\n\
    9   1\n\
   10   1\n\
   11   1\n\
   12   0\n\
   13   1\n\
   14   1.000000000e-16\n\
   15   1\n\
   16   0\n\
   17   0\n\
   18   0\n\
   19   0\n\
   20   2\n\
   21   1\n\
   22   0\n\
   23   1\n\
   24   123111321\n\
   25   0\n\
   26   1\n\
   27   1\n\
   28   1\n\
   29   1\n\
   30   1\n\
Block SMINPUTS\n\
    1   1.279340000e+02\n\
    2   1.166378700e-05\n\
    3   1.176000000e-01\n\
    4   9.118760000e+01\n\
    5   4.200000000e+00\n\
    6   1.733000000e+02\n\
    7   1.777000000e+00\n\
    8   0.000000000e+00\n\
    9   80.404\n\
   11   5.109989020e-04\n\
   12   0.000000000e+00\n\
   13   1.056583570e-01\n\
   14   0.000000000e+00\n\
   21   4.750000000e-03\n\
   22   2.400000000e-03\n\
   23   1.040000000e-01\n\
   24   1.270000000e+00\n\
Block MINPAR\n\
    1   500\n\
    2   300\n\
    3   10\n\
    5   -10\n\
Block EXTPAR\n\
   61   0.2\n\
   62   0.1\n\
   65   10\n\
   200  10\n\
   201  10\n\
   202  10\n\
Block YVIN\n\
    1   0.01\n\
    2   0.01\n\
    3   0.01\n\
";

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

   munuSSM_slha<munuSSM<Two_scale>> m = setup_munuSSM(input, qedqcd, settings);

#if defined(__INTEL_COMPILER)
   constexpr double reference_value = 1.24070126531925e-09;
#elif defined(__clang__)
   constexpr double reference_value = 1.24070101020353e-09;
#elif defined(__GNUC__)
   constexpr double reference_value = 1.24070097343634e-09;
#else
   BOOST_TEST_MESSAGE(false, "No reference value for this compiler");
#endif

   double amu = munuSSM_a_muon::calculate_a_muon(m, qedqcd);

   BOOST_CHECK_CLOSE_FRACTION(amu, reference_value, 2e-15);
}
