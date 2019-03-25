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

#include "test_munuSSM.hpp"

#include "wrappers.hpp"
#include "munuSSM_a_muon.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_amu ) {

   munuSSM_input_parameters input;

   // Block MINPAR
   input.m0 = 500.;
   input.m12 = 300.;
   input.TanBeta = 10.;
   input.Azero = -10.;

   // Block EXTPAR
   input.LambdaInput = 0.2;
   input.KappaInput = 0.1;
   input.vRInput = 10.;
   input.vL1Input = 10.;
   input.vL2Input = 10.;
   input.vL3Input = 10.;

   // Block YVIN
   input.YvInput << 0.01, 0.01, 0.01;

   softsusy::QedQcd qedqcd;

   munuSSM_slha<munuSSM<Two_scale>> m = setup_munuSSM(input, qedqcd);

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
