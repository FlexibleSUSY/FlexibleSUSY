#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSMCPVSemiAnalytic_ewsb

#include <boost/test/unit_test.hpp>

#include "test_CMSSMCPVSemiAnalytic.hpp"
#include "CMSSMCPVSemiAnalytic_semi_analytic_ewsb_solver.hpp"
#include "CMSSMCPVSemiAnalytic_semi_analytic_model.hpp"
#include "CMSSMCPVSemiAnalytic_semi_analytic_solutions.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_CMSSMCPVSemiAnalytic_ewsb_tree_level_solution )
{
   const double precision = 1.0e-5;
   CMSSMCPVSemiAnalytic_input_parameters input;
   CMSSMCPVSemiAnalytic_mass_eigenstates model(input);
   setup_CMSSMCPVSemiAnalytic(model, input);

   const double high_scale = 2.e16;
   model.run_to(high_scale);

   Boundary_values values;
   setup_high_scale_CMSSMCPVSemiAnalytic(model, values);

   CMSSMCPVSemiAnalytic_semi_analytic_solutions solns;
   solns.set_input_scale(high_scale);
   solns.set_output_scale(Electroweak_constants::MZ);

   solns.calculate_coefficients(model);

   model.run_to(Electroweak_constants::MZ);

   CMSSMCPVSemiAnalytic_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_loop_order(0);
   ewsb_solver.set_semi_analytic_solutions(&solns);
   model.set_ewsb_solver(
      std::make_shared<CMSSMCPVSemiAnalytic_ewsb_solver<Semi_analytic> >(ewsb_solver));

   model.set_ewsb_loop_order(0);
   model.set_ewsb_iteration_precision(precision);
   const int error = model.solve_ewsb_tree_level();
   BOOST_CHECK_EQUAL(error, 0);

   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_1(), precision);
   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_2(), precision);
   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_3(), precision);
   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_4(), precision);
}

BOOST_AUTO_TEST_CASE( test_CMSSMCPVSemiAnalytic_ewsb_one_loop_solution )
{
   const double precision = 1.0e-7;
   CMSSMCPVSemiAnalytic_input_parameters input;
   CMSSMCPVSemiAnalytic_mass_eigenstates model(input);
   setup_CMSSMCPVSemiAnalytic(model, input);

   const double high_scale = 2.e16;
   model.run_to(high_scale);

   Boundary_values values;
   setup_high_scale_CMSSMCPVSemiAnalytic(model, values);

   CMSSMCPVSemiAnalytic_semi_analytic_solutions solns;
   solns.set_input_scale(high_scale);
   solns.set_output_scale(Electroweak_constants::MZ);

   solns.calculate_coefficients(model);

   model.run_to(Electroweak_constants::MZ);
   model.calculate_DRbar_masses();

   CMSSMCPVSemiAnalytic_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_loop_order(1);
   ewsb_solver.set_semi_analytic_solutions(&solns);
   model.set_ewsb_solver(
      std::make_shared<CMSSMCPVSemiAnalytic_ewsb_solver<Semi_analytic> >(ewsb_solver));

   model.set_ewsb_loop_order(1);
   model.set_ewsb_iteration_precision(precision);

   const int neq = CMSSMCPVSemiAnalytic<Semi_analytic>::number_of_ewsb_equations;
   double tadpole[neq] = { 0. };

   model.tadpole_equations(tadpole);

   BOOST_CHECK(Abs(tadpole[0]) > 1000.);
   BOOST_CHECK(Abs(tadpole[1]) > 1000.);
   BOOST_CHECK(Abs(tadpole[2]) > 1000.);
   BOOST_CHECK(Abs(tadpole[3]) > 1000.);

   const int error = model.solve_ewsb_one_loop();
   BOOST_CHECK_EQUAL(error, 0);

   model.tadpole_equations(tadpole);

   BOOST_CHECK_SMALL(Abs(tadpole[0]), 0.5);
   BOOST_CHECK_SMALL(Abs(tadpole[1]), 3.);
   BOOST_CHECK_SMALL(Abs(tadpole[2]), 20.);
   BOOST_CHECK_SMALL(Abs(tadpole[3]), 0.1);
}
