#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_THDMIIEWSBAtMZSemiAnalytic_semi_analytic_solutions

#include <boost/test/unit_test.hpp>

#include "test_THDMIIEWSBAtMZSemiAnalytic.hpp"
#include "THDMIIEWSBAtMZSemiAnalytic_semi_analytic_solutions.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_THDMIIEWSBAtMZSemiAnalytic_coefficients )
{
   THDMIIEWSBAtMZSemiAnalytic_input_parameters input;
   THDMIIEWSBAtMZSemiAnalytic_mass_eigenstates model(input);
   setup_THDMIIEWSBAtMZSemiAnalytic(model, input);

   const double high_scale = model.get_input().Qin;
   model.run_to(high_scale);

   Boundary_values values;
   setup_high_scale_THDMIIEWSBAtMZSemiAnalytic(model, values);

   THDMIIEWSBAtMZSemiAnalytic_semi_analytic_solutions solns;
   solns.set_input_scale(high_scale);
   solns.set_output_scale(Electroweak_constants::MZ);

   solns.calculate_coefficients(model);

   model.run_to(Electroweak_constants::MZ);

   THDMIIEWSBAtMZSemiAnalytic_mass_eigenstates coeffs_model(model);
   solns.evaluate_solutions(coeffs_model);

   BOOST_CHECK_CLOSE_FRACTION(model.get_M112(), coeffs_model.get_M112(), 1.e-3);
   BOOST_CHECK_CLOSE_FRACTION(model.get_M122(), coeffs_model.get_M122(), 1.e-3);
   BOOST_CHECK_CLOSE_FRACTION(model.get_M222(), coeffs_model.get_M222(), 1.e-3);

}
