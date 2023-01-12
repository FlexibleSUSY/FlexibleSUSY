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
#define BOOST_TEST_MODULE test_SM_gmm2

#include <boost/test/unit_test.hpp>

#include "test_SM.hpp"

#include "spectrum_generator_settings.hpp"
#include "lowe.h"

#include "SM_two_scale_model.hpp"
#include "SM_amm.hpp"
#include "cxx_qft/SM_qft.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_zero )
{
   SM_input_parameters input;
   input.LambdaIN = 0.25;
   input.Qin = 1000;

   SM<Two_scale> sm;
   setup_SM_const(sm, input);
   sm.calculate_DRbar_masses();

   softsusy::QedQcd qedqcd;

   Spectrum_generator_settings settings;

   using SM_cxx_diagrams::fields::Fe;

   const double ae = SM_amm::calculate_amm<Fe>(sm, qedqcd, settings, 0);
   const double amu = SM_amm::calculate_amm<Fe>(sm, qedqcd, settings, 1);
   const double atau = SM_amm::calculate_amm<Fe>(sm, qedqcd, settings, 2);

   BOOST_CHECK_SMALL(ae, 1e-15);
   BOOST_CHECK_SMALL(amu, 1e-15);
   BOOST_CHECK_SMALL(atau, 1e-15);
}
