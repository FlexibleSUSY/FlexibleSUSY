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
#define BOOST_TEST_MODULE test_SM_cxxvertices

#include <boost/test/unit_test.hpp>
#include "test_complex_equality.hpp"

#include "concatenate.hpp"

#include "SM_two_scale_spectrum_generator.hpp"
#include "cxx_qft/SM_qft.hpp"
#include "cxx_qft/standard_model_qft.hpp"
#include "SM_mass_eigenstates_running.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_sm_cxxvertices )
{
   static constexpr double lambda = 0.12;

   const Spectrum_generator_settings settings;
   const softsusy::QedQcd qedqcd;

   SM_input_parameters input;
   input.LambdaIN = lambda;
   SM_spectrum_generator<Two_scale> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);
   auto _sm = std::get<0>(spectrum_generator.get_models_slha());
   _sm.set_Lambdax(lambda);
   _sm.solve_ewsb();
   _sm.calculate_DRbar_masses();
   SM_mass_eigenstates_running sm(_sm);

   standard_model::Standard_model standard_model {};
   standard_model.initialise_from_input(qedqcd);
   standard_model.set_Lambdax(lambda);
   standard_model.set_Yu(-sm.get_Yu());
   standard_model.set_Yd(sm.get_Yd());
   standard_model.set_Ye(sm.get_Ye());
   standard_model.set_v(sm.get_v());
   standard_model.solve_ewsb();
   standard_model.calculate_DRbar_masses();

   SM_cxx_diagrams::context_base context_sm {sm};
   standard_model_cxx_diagrams::context_base context_standard_model {standard_model};

   // h-t-tbar
   const auto indices = concatenate(std::array<int,0> {}, std::array<int, 1> {2}, std::array<int, 1>{2});
   const auto httbar_sm =
      SM_cxx_diagrams::Vertex<
         SM_cxx_diagrams::fields::hh,
         SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fu>::type,
         SM_cxx_diagrams::fields::Fu
      >::evaluate(indices, context_sm);
   const auto httbar_standard_model =
      standard_model_cxx_diagrams::Vertex<
         standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fu>::type,
         standard_model_cxx_diagrams::fields::Fu,
	 standard_model_cxx_diagrams::fields::hh
      >::evaluate(indices, context_standard_model);

   TEST_COMPLEX_CLOSE_FRACTION(httbar_sm.left(),  httbar_standard_model.left(), 1e-16);
   TEST_COMPLEX_CLOSE_FRACTION(httbar_sm.right(), httbar_standard_model.right(), 1e-16);

   // h-b-bbar
   const auto hbbbar_sm =
      SM_cxx_diagrams::Vertex<
         SM_cxx_diagrams::fields::hh,
         SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fd>::type,
         SM_cxx_diagrams::fields::Fd
      >::evaluate(indices, context_sm);
   const auto hbbbar_standard_model =
      standard_model_cxx_diagrams::Vertex<
         standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fd>::type,
         standard_model_cxx_diagrams::fields::Fd,
	 standard_model_cxx_diagrams::fields::hh
      >::evaluate(indices, context_standard_model);

   TEST_COMPLEX_CLOSE_FRACTION(hbbbar_sm.left(),  hbbbar_standard_model.left(), 1e-16);
   TEST_COMPLEX_CLOSE_FRACTION(hbbbar_sm.right(), hbbbar_standard_model.right(), 1e-16);

   // h-tau-taubar
   const auto htautaubar_sm =
      SM_cxx_diagrams::Vertex<
         SM_cxx_diagrams::fields::hh,
         SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fe>::type,
         SM_cxx_diagrams::fields::Fe
      >::evaluate(indices, context_sm);
   const auto htautaubar_standard_model =
      standard_model_cxx_diagrams::Vertex<
         standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fe>::type,
         standard_model_cxx_diagrams::fields::Fe,
	 standard_model_cxx_diagrams::fields::hh
      >::evaluate(indices, context_standard_model);

   TEST_COMPLEX_CLOSE_FRACTION(htautaubar_sm.left(),  htautaubar_standard_model.left(), 1e-16);
   TEST_COMPLEX_CLOSE_FRACTION(htautaubar_sm.right(), htautaubar_standard_model.right(), 1e-16);

   // h-h-h
   const auto indices_hhh = concatenate(std::array<int,0> {}, std::array<int, 0> {}, std::array<int, 0>{});
   const auto hhh_sm =
      SM_cxx_diagrams::Vertex<
         SM_cxx_diagrams::fields::hh,
         SM_cxx_diagrams::fields::hh,
         SM_cxx_diagrams::fields::hh
      >::evaluate(indices_hhh, context_sm);
   const auto hhh_standard_model =
      standard_model_cxx_diagrams::Vertex<
         standard_model_cxx_diagrams::fields::hh,
         standard_model_cxx_diagrams::fields::hh,
         standard_model_cxx_diagrams::fields::hh
      >::evaluate(indices_hhh, context_standard_model);

   TEST_COMPLEX_CLOSE_FRACTION(hhh_sm.value(),  hhh_standard_model.value(), 1e-16);
}
