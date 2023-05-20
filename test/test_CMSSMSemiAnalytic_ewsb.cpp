#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSMSemiAnalytic_ewsb

#include <boost/test/unit_test.hpp>

#include "test_CMSSMSemiAnalytic.hpp"
#include "CMSSMSemiAnalytic_semi_analytic_ewsb_solver.hpp"
#include "CMSSMSemiAnalytic_semi_analytic_solutions.hpp"
#include "root_finder.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_CMSSMSemiAnalytic_ewsb_tree_level_solution )
{
   const double precision = 1.0e-5;
   CMSSMSemiAnalytic_input_parameters input;
   CMSSMSemiAnalytic_mass_eigenstates model(input);
   setup_CMSSMSemiAnalytic(model, input);

   const double high_scale = 2.e16;
   model.run_to(high_scale);

   Boundary_values values;
   setup_high_scale_CMSSMSemiAnalytic(model, values);

   CMSSMSemiAnalytic_semi_analytic_solutions solns;
   solns.set_input_scale(high_scale);
   solns.set_output_scale(Electroweak_constants::MZ);

   solns.calculate_coefficients(model);

   model.run_to(Electroweak_constants::MZ);

   CMSSMSemiAnalytic_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_loop_order(0);
   ewsb_solver.set_semi_analytic_solutions(&solns);
   model.set_ewsb_solver(
      std::make_shared<CMSSMSemiAnalytic_ewsb_solver<Semi_analytic> >(ewsb_solver));

   model.set_ewsb_loop_order(0);
   model.set_ewsb_iteration_precision(precision);
   const int error = model.solve_ewsb_tree_level();
   BOOST_CHECK_EQUAL(error, 0);

   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_1(), precision);
   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_2(), precision);
}

BOOST_AUTO_TEST_CASE( test_CMSSMSemiAnalytic_ewsb_one_loop_solution )
{
   const double precision = 1.0e-5;
   CMSSMSemiAnalytic_input_parameters input;
   CMSSMSemiAnalytic_mass_eigenstates model(input);
   setup_CMSSMSemiAnalytic(model, input);

   const double high_scale = 2.e16;
   model.run_to(high_scale);

   Boundary_values values;
   setup_high_scale_CMSSMSemiAnalytic(model, values);

   CMSSMSemiAnalytic_semi_analytic_solutions solns;
   solns.set_input_scale(high_scale);
   solns.set_output_scale(Electroweak_constants::MZ);

   solns.calculate_coefficients(model);

   model.run_to(Electroweak_constants::MZ);
   model.calculate_DRbar_masses();

   CMSSMSemiAnalytic_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_loop_order(1);
   ewsb_solver.set_semi_analytic_solutions(&solns);
   model.set_ewsb_solver(
      std::make_shared<CMSSMSemiAnalytic_ewsb_solver<Semi_analytic> >(ewsb_solver));

   model.set_ewsb_loop_order(1);
   model.set_ewsb_iteration_precision(precision);
   const int error = model.solve_ewsb_one_loop();
   BOOST_CHECK_EQUAL(error, 0);

   const std::complex<double> tadpole_hh_1(model.tadpole_hh_1loop(0));
   const std::complex<double> tadpole_hh_2(model.tadpole_hh_1loop(1));

   BOOST_CHECK_SMALL(Im(tadpole_hh_1), 1.0e-12);
   BOOST_CHECK_SMALL(Im(tadpole_hh_2), 1.0e-12);

   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_1() - Re(tadpole_hh_1), 5.);
   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_2() - Re(tadpole_hh_2), 5.);
}

BOOST_AUTO_TEST_CASE( test_CMSSMSemiAnalytic_ewsb_two_loop_solution )
{
   const double precision = 1.0e-5;
   CMSSMSemiAnalytic_input_parameters input;
   CMSSMSemiAnalytic_mass_eigenstates model(input);
   setup_CMSSMSemiAnalytic(model, input);

   const double high_scale = 2.e16;
   model.run_to(high_scale);

   Boundary_values values;
   setup_high_scale_CMSSMSemiAnalytic(model, values);

   CMSSMSemiAnalytic_semi_analytic_solutions solns;
   solns.set_input_scale(high_scale);
   solns.set_output_scale(Electroweak_constants::MZ);

   solns.calculate_coefficients(model);

   model.run_to(Electroweak_constants::MZ);
   model.calculate_DRbar_masses();

   CMSSMSemiAnalytic_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_loop_order(2);
   ewsb_solver.set_semi_analytic_solutions(&solns);
   model.set_ewsb_solver(
      std::make_shared<CMSSMSemiAnalytic_ewsb_solver<Semi_analytic> >(ewsb_solver));

   model.set_ewsb_loop_order(2);
   model.set_ewsb_iteration_precision(precision);
   const int error = model.solve_ewsb();
   BOOST_CHECK_EQUAL(error, 0);

   const Eigen::Matrix<double,2,1> tadpole_2l(model.tadpole_hh_2loop());
   
   const std::complex<double> tadpole_hh_1(
      model.tadpole_hh_1loop(0) + tadpole_2l(0));
   const std::complex<double> tadpole_hh_2(
      model.tadpole_hh_1loop(1) + tadpole_2l(1));

   BOOST_CHECK_SMALL(Im(tadpole_hh_1), 1.0e-12);
   BOOST_CHECK_SMALL(Im(tadpole_hh_2), 1.0e-12);

   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_1() - Re(tadpole_hh_1), 500.);
   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_2() - Re(tadpole_hh_2), 500.);
}

class High_scale_ewsb_solver {
public:
   void set_loop_order(int l) { loop_order = l; }
   void set_number_of_iterations(int it) { number_of_iterations = it; }
   void set_precision(double p) { precision = p; }
   void set_high_scale(double s) { high_scale = s; }
   void set_ewsb_scale(double s) { ewsb_scale = s; }

   int solve_ewsb_at_high_scale(CMSSMSemiAnalytic_mass_eigenstates& model) const;

private:
   int loop_order{2};
   int number_of_iterations{100};
   double precision{1.e-5};
   double high_scale{2.e16};
   double ewsb_scale{Electroweak_constants::MZ};
};

int High_scale_ewsb_solver::solve_ewsb_at_high_scale(CMSSMSemiAnalytic_mass_eigenstates& model) const
{
   auto m = model;
   m.set_ewsb_loop_order(loop_order);

   auto tadpole_stepper = [this, m](const Eigen::Matrix<double,2,1>& ewsb_pars) mutable -> Eigen::Matrix<double,2,1> {
      auto running_model = m;

      const double BMu0 = ewsb_pars(0);
      const double m0Sq = ewsb_pars(1);

      running_model.run_to(this->high_scale);

      running_model.set_mHu2(m0Sq);
      running_model.set_mHd2(m0Sq);
      running_model.set_mq2(m0Sq * UNITMATRIX(3));
      running_model.set_mu2(m0Sq * UNITMATRIX(3));
      running_model.set_md2(m0Sq * UNITMATRIX(3));
      running_model.set_ml2(m0Sq * UNITMATRIX(3));
      running_model.set_me2(m0Sq * UNITMATRIX(3));
      running_model.set_BMu(BMu0);

      running_model.run_to(this->ewsb_scale);

      if (this->loop_order > 0)
         running_model.calculate_DRbar_masses();

      return running_model.tadpole_equations();
   };

   std::unique_ptr<EWSB_solver> solvers[] = {
      std::unique_ptr<EWSB_solver>(new Root_finder<2>(tadpole_stepper, number_of_iterations, precision, Root_finder<2>::GSLHybridS)),
      std::unique_ptr<EWSB_solver>(new Root_finder<2>(tadpole_stepper, number_of_iterations, precision, Root_finder<2>::GSLBroyden))
   };

   Eigen::Matrix<double,2,1> x_init;
   x_init(0) = model.get_BMu();
   x_init(1) = model.get_mHd2();

   int status;
   for (auto& solver: solvers) {
      status = solver->solve(x_init);
      if (status == EWSB_solver::SUCCESS) {
         const auto solution = solver->get_solution();

         const double BMu0 = solution(0);
         const double m0Sq = solution(1);

         model.set_BMu0(BMu0);
         model.set_m0Sq(m0Sq);

         model.run_to(high_scale);

         model.set_mHu2(m0Sq);
         model.set_mHd2(m0Sq);
         model.set_mq2(m0Sq * UNITMATRIX(3));
         model.set_mu2(m0Sq * UNITMATRIX(3));
         model.set_md2(m0Sq * UNITMATRIX(3));
         model.set_ml2(m0Sq * UNITMATRIX(3));
         model.set_me2(m0Sq * UNITMATRIX(3));
         model.set_BMu(BMu0);

         model.run_to(ewsb_scale);

         model.calculate_DRbar_masses();

         break;
      }
   }

   return status;
}

// this test checks that solving the EWSB equations using the
// semi-analytic solutions is equivalent to solving them by
// directly setting the parameters at the high-scale
BOOST_AUTO_TEST_CASE( test_semi_analytic_ewsb_tree_level )
{
   const double precision = 1.0e-5;
   CMSSMSemiAnalytic_input_parameters input;
   CMSSMSemiAnalytic_mass_eigenstates model(input);
   setup_CMSSMSemiAnalytic(model, input);

   const double high_scale = 2.e16;
   model.run_to(high_scale);

   Boundary_values values;
   setup_high_scale_CMSSMSemiAnalytic(model, values);

   CMSSMSemiAnalytic_mass_eigenstates high_scale_model(model);

   CMSSMSemiAnalytic_semi_analytic_solutions solns;
   solns.set_input_scale(high_scale);
   solns.set_output_scale(Electroweak_constants::MZ);

   solns.calculate_coefficients(model);

   model.run_to(Electroweak_constants::MZ);

   CMSSMSemiAnalytic_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_loop_order(0);
   ewsb_solver.set_semi_analytic_solutions(&solns);
   model.set_ewsb_solver(
      std::make_shared<CMSSMSemiAnalytic_ewsb_solver<Semi_analytic> >(ewsb_solver));

   model.set_ewsb_iteration_precision(precision);
   const int error = model.solve_ewsb_tree_level();
   BOOST_CHECK_EQUAL(error, 0);

   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_1(), precision);
   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_2(), precision);

   const double m0Sq_sol = model.get_m0Sq();
   const double BMu0_sol = model.get_BMu0();

   High_scale_ewsb_solver high_scale_solver;
   high_scale_solver.set_high_scale(high_scale);
   high_scale_solver.set_ewsb_scale(Electroweak_constants::MZ);
   high_scale_solver.set_precision(precision);
   high_scale_solver.set_loop_order(0);
   
   const int high_scale_error = high_scale_solver.solve_ewsb_at_high_scale(high_scale_model);
   BOOST_CHECK_EQUAL(high_scale_error, 0);

   BOOST_CHECK_SMALL(high_scale_model.get_ewsb_eq_hh_1(), precision);
   BOOST_CHECK_SMALL(high_scale_model.get_ewsb_eq_hh_2(), precision);

   BOOST_CHECK_CLOSE(m0Sq_sol, high_scale_model.get_m0Sq(), 1.);
   BOOST_CHECK_CLOSE(BMu0_sol, high_scale_model.get_BMu0(), 1.);
}

BOOST_AUTO_TEST_CASE( test_semi_analytic_ewsb_one_loop )
{
   const double precision = 1.0e-5;
   CMSSMSemiAnalytic_input_parameters input;
   CMSSMSemiAnalytic_mass_eigenstates model(input);
   setup_CMSSMSemiAnalytic(model, input);

   const double high_scale = 2.e16;
   model.run_to(high_scale);

   Boundary_values values;
   setup_high_scale_CMSSMSemiAnalytic(model, values);

   CMSSMSemiAnalytic_mass_eigenstates high_scale_model(model);

   CMSSMSemiAnalytic_semi_analytic_solutions solns;
   solns.set_input_scale(high_scale);
   solns.set_output_scale(Electroweak_constants::MZ);

   solns.calculate_coefficients(model);

   model.run_to(Electroweak_constants::MZ);
   model.calculate_DRbar_masses();

   CMSSMSemiAnalytic_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_loop_order(1);
   ewsb_solver.set_semi_analytic_solutions(&solns);
   model.set_ewsb_solver(
      std::make_shared<CMSSMSemiAnalytic_ewsb_solver<Semi_analytic> >(ewsb_solver));

   model.set_ewsb_iteration_precision(precision);
   const int error = model.solve_ewsb_one_loop();
   BOOST_CHECK_EQUAL(error, 0);

   const std::complex<double> tadpole_hh_1(model.tadpole_hh_1loop(0));
   const std::complex<double> tadpole_hh_2(model.tadpole_hh_1loop(1));

   BOOST_CHECK_SMALL(Im(tadpole_hh_1), 1.0e-12);
   BOOST_CHECK_SMALL(Im(tadpole_hh_2), 1.0e-12);

   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_1() - Re(tadpole_hh_1), 5.);
   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_2() - Re(tadpole_hh_2), 5.);

   const double m0Sq_sol = model.get_m0Sq();
   const double BMu0_sol = model.get_BMu0();

   High_scale_ewsb_solver high_scale_solver;
   high_scale_solver.set_high_scale(high_scale);
   high_scale_solver.set_ewsb_scale(Electroweak_constants::MZ);
   high_scale_solver.set_precision(precision);
   high_scale_solver.set_loop_order(1);
   
   const int high_scale_error = high_scale_solver.solve_ewsb_at_high_scale(high_scale_model);
   BOOST_CHECK_EQUAL(high_scale_error, 0);

   const std::complex<double> running_tadpole_hh_1(high_scale_model.tadpole_hh_1loop(0));
   const std::complex<double> running_tadpole_hh_2(high_scale_model.tadpole_hh_1loop(1));

   BOOST_CHECK_SMALL(Im(running_tadpole_hh_1), 1.0e-12);
   BOOST_CHECK_SMALL(Im(running_tadpole_hh_2), 1.0e-12);

   BOOST_CHECK_SMALL(high_scale_model.get_ewsb_eq_hh_1() - Re(running_tadpole_hh_1), 5.);
   BOOST_CHECK_SMALL(high_scale_model.get_ewsb_eq_hh_2() - Re(running_tadpole_hh_2), 5.);

   BOOST_CHECK_CLOSE(m0Sq_sol, high_scale_model.get_m0Sq(), 2.);
   BOOST_CHECK_CLOSE(BMu0_sol, high_scale_model.get_BMu0(), 2.);
}
