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
#define BOOST_TEST_MODULE test_SM_yukawa_convention

#include <boost/test/unit_test.hpp>
#include "test_complex_equality.hpp"
#include "SM_two_scale_spectrum_generator.hpp"
#include "standard_model.hpp"

inline int
sgn(double v) {
    return (v > 0) - (v < 0);
}

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_yukawa_convention )
{
   static constexpr double lambda = 0.12;

   const Spectrum_generator_settings settings;
   const softsusy::QedQcd qedqcd;

   SM_input_parameters input;
   input.LambdaIN = lambda;
   SM_spectrum_generator<Two_scale> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);
   auto sm = std::get<0>(spectrum_generator.get_models_slha());
   sm.set_Lambdax(lambda);
   sm.solve_ewsb();
   sm.calculate_DRbar_masses();

   standard_model::Standard_model standard_model {};
   standard_model.initialise_from_input(qedqcd);
   standard_model.set_Lambdax(lambda);
   standard_model.solve_ewsb();
   standard_model.calculate_DRbar_masses();

   // builtin SM has an opostite sign convention for Yu
   BOOST_CHECK(
      sm.get_Yu().unaryExpr(&sgn) == -standard_model.get_Yu().unaryExpr(&sgn)
   );

   // but not for Yd and Ye
   BOOST_CHECK(
      sm.get_Yd().unaryExpr(&sgn) == standard_model.get_Yd().unaryExpr(&sgn)
   );
   BOOST_CHECK(
      sm.get_Ye().unaryExpr(&sgn) == standard_model.get_Ye().unaryExpr(&sgn)
   );

   BOOST_CHECK(sgn(sm.get_g1()) == sgn(standard_model.get_g1()));
   BOOST_CHECK(sgn(sm.get_g2()) == sgn(standard_model.get_g2()));
   BOOST_CHECK(sgn(sm.get_g3()) == sgn(standard_model.get_g3()));
   BOOST_CHECK(sgn(sm.get_v()) == sgn(standard_model.get_v()));
}
