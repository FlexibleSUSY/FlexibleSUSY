#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SMSemiAnalytic_ewsb

#include <boost/test/unit_test.hpp>

#include "test_SMSemiAnalytic.hpp"
#include "SMSemiAnalytic_semi_analytic_ewsb_solver.hpp"
#include "SMSemiAnalytic_semi_analytic_solutions.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_SMSemiAnalytic_ewsb_tree_level_solution )
{
   const double precision = 1.e-5;
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

   SMSemiAnalytic_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_loop_order(0);
   ewsb_solver.set_semi_analytic_solutions(&solns);
   model.set_ewsb_solver(
      std::make_shared<SMSemiAnalytic_ewsb_solver<Semi_analytic> >(ewsb_solver));

   model.set_ewsb_loop_order(0);
   model.set_ewsb_iteration_precision(precision);
   const int error = model.solve_ewsb_tree_level();
   BOOST_CHECK_EQUAL(error, 0);

   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_1(), precision);
}

BOOST_AUTO_TEST_CASE( test_SMSemiAnalytic_ewsb_one_loop_solution )
{
   const double precision = 1.e-5;
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

   SMSemiAnalytic_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_loop_order(1);
   ewsb_solver.set_semi_analytic_solutions(&solns);
   model.set_ewsb_solver(
      std::make_shared<SMSemiAnalytic_ewsb_solver<Semi_analytic> >(ewsb_solver));

   model.set_ewsb_loop_order(1);
   model.set_ewsb_iteration_precision(precision);
   const int error = model.solve_ewsb_one_loop();
   BOOST_CHECK_EQUAL(error, 0);

   const std::complex<double> tadpole_hh(model.tadpole_hh_1loop());

   BOOST_CHECK_SMALL(Im(tadpole_hh), 1.0e-12);

   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_1() - Re(tadpole_hh), 0.1);
}
