#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_munuSSMSemiAnalytic_ewsb

#include <boost/test/unit_test.hpp>

#include "test_munuSSMSemiAnalytic.hpp"
#include "munuSSMSemiAnalytic_semi_analytic_ewsb_solver.hpp"
#include "munuSSMSemiAnalytic_semi_analytic_solutions.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_munuSSMSemiAnalytic_ewsb_tree_level_solution )
{
   const double precision = 1.0e-5;
   munuSSMSemiAnalytic_input_parameters input;
   munuSSMSemiAnalytic_mass_eigenstates model(input);
   setup_munuSSMSemiAnalytic(model, input);

   const double high_scale = 2.e16;
   model.run_to(high_scale);

   Boundary_values values;
   setup_high_scale_munuSSMSemiAnalytic(model, values);

   munuSSMSemiAnalytic_semi_analytic_solutions solns;
   solns.set_input_scale(high_scale);
   solns.set_output_scale(Electroweak_constants::MZ);

   solns.calculate_coefficients(model);

   model.run_to(Electroweak_constants::MZ);

   munuSSMSemiAnalytic_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_loop_order(0);
   ewsb_solver.set_semi_analytic_solutions(&solns);
   model.set_ewsb_solver(
      std::make_shared<munuSSMSemiAnalytic_ewsb_solver<Semi_analytic> >(ewsb_solver));

   model.set_ewsb_loop_order(0);
   model.set_ewsb_iteration_precision(precision);
   const int error = model.solve_ewsb_tree_level();
   BOOST_CHECK_EQUAL(error, 0);

   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_1(), precision);
   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_2(), precision);
   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_3(), precision);
   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_4(), precision);
   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_5(), precision);
   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_6(), precision);
}

BOOST_AUTO_TEST_CASE( test_munuSSMSemiAnalytic_ewsb_one_loop_solution )
{
   const double precision = 1.0e-5;
   munuSSMSemiAnalytic_input_parameters input;
   munuSSMSemiAnalytic_mass_eigenstates model(input);
   setup_munuSSMSemiAnalytic(model, input);

   const double high_scale = 2.e16;
   model.run_to(high_scale);

   Boundary_values values;
   setup_high_scale_munuSSMSemiAnalytic(model, values);

   munuSSMSemiAnalytic_semi_analytic_solutions solns;
   solns.set_input_scale(high_scale);
   solns.set_output_scale(Electroweak_constants::MZ);

   solns.calculate_coefficients(model);

   model.run_to(Electroweak_constants::MZ);
   model.calculate_DRbar_masses();

   munuSSMSemiAnalytic_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_loop_order(1);
   ewsb_solver.set_semi_analytic_solutions(&solns);
   model.set_ewsb_solver(
      std::make_shared<munuSSMSemiAnalytic_ewsb_solver<Semi_analytic> >(ewsb_solver));

   model.set_ewsb_loop_order(1);
   model.set_ewsb_iteration_precision(precision);
   const int error = model.solve_ewsb_one_loop();
   BOOST_CHECK_EQUAL(error, 0);

   const std::complex<double> tadpole_hh_1(model.tadpole_hh_1loop(0));
   const std::complex<double> tadpole_hh_2(model.tadpole_hh_1loop(1));
   const std::complex<double> tadpole_hh_3(model.tadpole_hh_1loop(3));
   const std::complex<double> tadpole_hh_4(model.tadpole_hh_1loop(4));
   const std::complex<double> tadpole_hh_5(model.tadpole_hh_1loop(5));
   const std::complex<double> tadpole_hh_6(model.tadpole_hh_1loop(2));

   BOOST_CHECK_SMALL(Im(tadpole_hh_1), 1.0e-12);
   BOOST_CHECK_SMALL(Im(tadpole_hh_2), 1.0e-12);
   BOOST_CHECK_SMALL(Im(tadpole_hh_3), 1.0e-12);
   BOOST_CHECK_SMALL(Im(tadpole_hh_4), 1.0e-12);
   BOOST_CHECK_SMALL(Im(tadpole_hh_5), 1.0e-12);
   BOOST_CHECK_SMALL(Im(tadpole_hh_6), 1.0e-12);

   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_1() - Re(tadpole_hh_1), 15.);
   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_2() - Re(tadpole_hh_2), 15.);
   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_3() - Re(tadpole_hh_3), 15.);
   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_4() - Re(tadpole_hh_4), 15.);
   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_5() - Re(tadpole_hh_5), 15.);
   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_6() - Re(tadpole_hh_6), 15.);
}
