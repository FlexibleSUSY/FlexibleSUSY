#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SMSemiAnalytic_semi_analytic_solutions

#include <boost/test/unit_test.hpp>

#include "test_SMSemiAnalytic.hpp"
#include "SMSemiAnalytic_semi_analytic_solutions.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_SMSemiAnalytic_coefficients )
{
   SMSemiAnalytic_input_parameters input;
   SMSemiAnalytic_mass_eigenstates model(input);
   setup_SMSemiAnalytic(model, input);

   const double high_scale = 2.e16;
   model.run_to(high_scale);

   Boundary_values values;
   setup_high_scale_SMSemiAnalytic(model, values);

   SMSemiAnalytic_semi_analytic_solutions solns;
   solns.set_input_scale(high_scale);
   solns.set_output_scale(Electroweak_constants::MZ);

   solns.calculate_coefficients(model);

   model.run_to(Electroweak_constants::MZ);

   SMSemiAnalytic_mass_eigenstates coeffs_model(model);
   solns.evaluate_solutions(coeffs_model);

   BOOST_CHECK_CLOSE_FRACTION(model.get_mu2(), coeffs_model.get_mu2(), 1.0e-3);
}
