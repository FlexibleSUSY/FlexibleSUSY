#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SSMSemiAnalytic_semi_analytic_solutions

#include <boost/test/unit_test.hpp>

#include "test_SSMSemiAnalytic.hpp"
#include "SSMSemiAnalytic_semi_analytic_solutions.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_SSMSemiAnalytic_coefficients )
{
   SSMSemiAnalytic_input_parameters input;
   SSMSemiAnalytic_mass_eigenstates model(input);
   setup_SSMSemiAnalytic(model, input);

   const double high_scale = 2.e16;
   model.run_to(high_scale);

   Boundary_values values;
   setup_high_scale_SSMSemiAnalytic(model, values);

   SSMSemiAnalytic_semi_analytic_solutions solns;
   solns.set_input_scale(high_scale);
   solns.set_output_scale(Electroweak_constants::MZ);

   solns.calculate_coefficients(model);

   model.run_to(Electroweak_constants::MZ);

   SSMSemiAnalytic_mass_eigenstates coeffs_model(model);
   solns.evaluate_solutions(coeffs_model);

   BOOST_CHECK_CLOSE_FRACTION(model.get_Kappa(), coeffs_model.get_Kappa(), 1.e-3);
   BOOST_CHECK_CLOSE_FRACTION(model.get_K1(), coeffs_model.get_K1(), 1.e-3);
   BOOST_CHECK_CLOSE_FRACTION(model.get_MS(), coeffs_model.get_MS(), 1.e-3);
   BOOST_CHECK_CLOSE_FRACTION(model.get_mu2(), coeffs_model.get_mu2(), 1.e-3);
}
