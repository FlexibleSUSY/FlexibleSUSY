#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CNMSSM_ewsb

#include <boost/test/unit_test.hpp>

#include "test_CNMSSM.hpp"
#include "CNMSSM_semi_analytic_ewsb_solver.hpp"
#include "CNMSSM_semi_analytic_solutions.hpp"
#include "root_finder.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_CNMSSM_ewsb_tree_level_solution )
{
   const double precision = 1.0e-5;
   CNMSSM_input_parameters input;
   CNMSSM_mass_eigenstates model(input);
   setup_CNMSSM(model, input);

   const double high_scale = 2.0e16;
   model.run_to(high_scale);

   Boundary_values values;
   setup_high_scale_CNMSSM(model, values);

   CNMSSM_semi_analytic_solutions solns;
   solns.set_input_scale(high_scale);
   solns.set_output_scale(Electroweak_constants::MZ);

   solns.calculate_coefficients(model);

   model.run_to(Electroweak_constants::MZ);

   CNMSSM_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_loop_order(0);
   ewsb_solver.set_semi_analytic_solutions(&solns);
   model.set_ewsb_solver(
      std::make_shared<CNMSSM_ewsb_solver<Semi_analytic> >(ewsb_solver));

   model.set_ewsb_loop_order(0);
   model.set_ewsb_iteration_precision(precision);
   const int error = model.solve_ewsb_tree_level();
   BOOST_CHECK_EQUAL(error, 0);

   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_1(), precision);
   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_2(), precision);
   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_3(), precision);
}

BOOST_AUTO_TEST_CASE( test_CNMSSM_ewsb_one_loop_solution )
{
   const double precision = 1.0e-5;
   CNMSSM_input_parameters input;
   CNMSSM_mass_eigenstates model(input);
   setup_CNMSSM(model, input);

   const double high_scale = 2.e16;
   model.run_to(high_scale);

   Boundary_values values;
   setup_high_scale_CNMSSM(model, values);

   CNMSSM_semi_analytic_solutions solns;
   solns.set_input_scale(high_scale);
   solns.set_output_scale(Electroweak_constants::MZ);

   solns.calculate_coefficients(model);

   model.run_to(Electroweak_constants::MZ);
   model.calculate_DRbar_masses();

   CNMSSM_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_loop_order(1);
   ewsb_solver.set_semi_analytic_solutions(&solns);
   model.set_ewsb_solver(
      std::make_shared<CNMSSM_ewsb_solver<Semi_analytic> >(ewsb_solver));

   model.set_ewsb_loop_order(1);
   model.set_ewsb_iteration_precision(precision);
   const int error = model.solve_ewsb_one_loop();
   BOOST_CHECK_EQUAL(error, 0);

   const std::complex<double> tadpole_hh_1(model.tadpole_hh_1loop(0));
   const std::complex<double> tadpole_hh_2(model.tadpole_hh_1loop(1));
   const std::complex<double> tadpole_hh_3(model.tadpole_hh_1loop(2));

   BOOST_CHECK_SMALL(Im(tadpole_hh_1), 1.0e-12);
   BOOST_CHECK_SMALL(Im(tadpole_hh_2), 1.0e-12);
   BOOST_CHECK_SMALL(Im(tadpole_hh_3), 1.0e-12);

   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_1() - Re(tadpole_hh_1), 1.0);
   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_2() - Re(tadpole_hh_2), 1.0);
   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_3() - Re(tadpole_hh_3), 1.0);
}

BOOST_AUTO_TEST_CASE( test_CNMSSM_ewsb_two_loop_solution )
{
   const double precision = 1.0e-5;
   CNMSSM_input_parameters input;
   CNMSSM_mass_eigenstates model(input);
   setup_CNMSSM(model, input);

   const double high_scale = 2.e16;
   model.run_to(high_scale);

   Boundary_values values;
   setup_high_scale_CNMSSM(model, values);

   CNMSSM_semi_analytic_solutions solns;
   solns.set_input_scale(high_scale);
   solns.set_output_scale(Electroweak_constants::MZ);

   solns.calculate_coefficients(model);

   model.run_to(Electroweak_constants::MZ);
   model.calculate_DRbar_masses();

   CNMSSM_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_loop_order(2);
   ewsb_solver.set_semi_analytic_solutions(&solns);
   model.set_ewsb_solver(
      std::make_shared<CNMSSM_ewsb_solver<Semi_analytic> >(ewsb_solver));

   model.set_ewsb_loop_order(2);
   model.set_ewsb_iteration_precision(precision);
   const int error = model.solve_ewsb();
   BOOST_CHECK_EQUAL(error, 0);

   const Eigen::Matrix<double,3,1> tadpole_2l(model.tadpole_hh_2loop());

   const std::complex<double> tadpole_hh_1(
      model.tadpole_hh_1loop(0) + tadpole_2l(0));
   const std::complex<double> tadpole_hh_2(
      model.tadpole_hh_1loop(1) + tadpole_2l(1));
   const std::complex<double> tadpole_hh_3(
      model.tadpole_hh_1loop(2) + tadpole_2l(2));

   BOOST_CHECK_SMALL(Im(tadpole_hh_1), 1.0e-12);
   BOOST_CHECK_SMALL(Im(tadpole_hh_2), 1.0e-12);
   BOOST_CHECK_SMALL(Im(tadpole_hh_3), 1.0e-12);

   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_1() - Re(tadpole_hh_1), 1.0);
   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_2() - Re(tadpole_hh_2), 1.0);
   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_3() - Re(tadpole_hh_3), 1.0);
}
