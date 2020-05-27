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
#define BOOST_TEST_MODULE test_CE6SSM_gmm2

#include <boost/test/unit_test.hpp>

#include "test_CE6SSM.hpp"

#include "CE6SSM_a_muon.hpp"
#include "CE6SSM_semi_analytic_spectrum_generator.hpp"
#include "CE6SSM_slha_io.hpp"
#include "CE6SSM_spectrum_generator.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_amu )
{
   
     char const * const slha_input = R"(
Block MODSEL
Block FlexibleSUSY
    0   1.000000000e-04
    1   0
    2   2
    3   0
    4   2
    5   2
    6   2
    7   2
    8   1
    9   1
   10   1
   11   1
   12   0
   13   1
   14   1.000000000e-11
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
   11   5.109989020e-04
   12   0.000000000e+00
   13   1.056583570e-01
   14   0.000000000e+00
   21   4.750000000e-03
   22   2.400000000e-03
   23   1.040000000e-01
   24   1.270000000e+00
Block MINPAR
    3   7.1
Block EXTPAR
   61   0.4         
   62   0.171        
   63   7276          
   64   0.26           
   65   7400            
   66   0.16             
   67   1.6e7             
   68   400                
   69   400                 
)";

   std::stringstream istr(slha_input);

   CE6SSM_slha_io slha_io;
   slha_io.read_from_stream(istr);

   // extract the input parameters
   softsusy::QedQcd qedqcd;
   CE6SSM_input_parameters input;
   Spectrum_generator_settings settings;

   try {
      slha_io.fill(qedqcd);
      slha_io.fill(input);
      slha_io.fill(settings);
   } catch (const Error& error) {
      BOOST_TEST_MESSAGE(error.what());
      BOOST_TEST(false);
   }

   CE6SSM_spectrum_generator<Semi_analytic> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.set_parameter_output_scale(slha_io.get_parameter_output_scale());

   spectrum_generator.run(qedqcd, input);

   auto models = spectrum_generator.get_models_slha();
   auto amu = CE6SSM_a_muon::calculate_a_muon(std::get<0>(models), qedqcd);

   constexpr double reference_value = 1.82135849E-11;

   BOOST_CHECK_CLOSE_FRACTION(amu, reference_value, 1e-3);

}
