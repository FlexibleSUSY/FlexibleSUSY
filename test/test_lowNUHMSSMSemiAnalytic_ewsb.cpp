#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_lowNUHMSSMSemiAnalytic_ewsb

#include <boost/test/unit_test.hpp>

#include "test_lowNUHMSSMSemiAnalytic.hpp"
#include "lowNUHMSSMSemiAnalytic_semi_analytic_ewsb_solver.hpp"
#include "lowNUHMSSMSemiAnalytic_semi_analytic_solutions.hpp"
#include "root_finder.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_lowNUHMSSMSemiAnalytic_ewsb_tree_level_solution )
{
   const double precision = 1.0e-5;
   lowNUHMSSMSemiAnalytic_input_parameters input;
   lowNUHMSSMSemiAnalytic_mass_eigenstates model(input);
   setup_lowNUHMSSMSemiAnalytic(model, input);

   const double susy_scale = 1000.;
   model.run_to(susy_scale);

   Boundary_values values;
   setup_susy_scale_lowNUHMSSMSemiAnalytic(model, values);

   lowNUHMSSMSemiAnalytic_semi_analytic_solutions solns;
   solns.set_input_scale(susy_scale);
   solns.set_output_scale(Electroweak_constants::MZ);

   solns.calculate_coefficients(model);

   model.run_to(Electroweak_constants::MZ);

   lowNUHMSSMSemiAnalytic_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_loop_order(0);
   ewsb_solver.set_semi_analytic_solutions(&solns);
   model.set_ewsb_solver(
      std::make_shared<lowNUHMSSMSemiAnalytic_ewsb_solver<Semi_analytic> >(ewsb_solver));

   model.set_ewsb_loop_order(0);
   model.set_ewsb_iteration_precision(precision);
   const int error = model.solve_ewsb_tree_level();
   BOOST_CHECK_EQUAL(error, 0);

   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_1(), precision);
   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_2(), precision);
}

BOOST_AUTO_TEST_CASE( test_lowNUHMSSMSemiAnalytic_ewsb_one_loop_solution )
{
   const double precision = 1.0e-5;
   lowNUHMSSMSemiAnalytic_input_parameters input;
   lowNUHMSSMSemiAnalytic_mass_eigenstates model(input);
   setup_lowNUHMSSMSemiAnalytic(model, input);

   const double susy_scale = 1000.;
   model.run_to(susy_scale);

   Boundary_values values;
   setup_susy_scale_lowNUHMSSMSemiAnalytic(model, values);

   lowNUHMSSMSemiAnalytic_semi_analytic_solutions solns;
   solns.set_input_scale(susy_scale);
   solns.set_output_scale(Electroweak_constants::MZ);

   solns.calculate_coefficients(model);

   model.run_to(Electroweak_constants::MZ);
   model.calculate_DRbar_masses();

   lowNUHMSSMSemiAnalytic_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_loop_order(1);
   ewsb_solver.set_semi_analytic_solutions(&solns);
   model.set_ewsb_solver(
      std::make_shared<lowNUHMSSMSemiAnalytic_ewsb_solver<Semi_analytic> >(ewsb_solver));

   model.set_ewsb_loop_order(1);
   model.set_ewsb_iteration_precision(precision);
   const int error = model.solve_ewsb_one_loop();
   BOOST_CHECK_EQUAL(error, 0);

   const std::complex<double> tadpole_hh_1(model.tadpole_hh_1loop(0));
   const std::complex<double> tadpole_hh_2(model.tadpole_hh_1loop(1));

   BOOST_CHECK_SMALL(Im(tadpole_hh_1), 1.0e-12);
   BOOST_CHECK_SMALL(Im(tadpole_hh_2), 1.0e-12);

   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_1() - Re(tadpole_hh_1), 1.);
   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_2() - Re(tadpole_hh_2), 1.);
}

BOOST_AUTO_TEST_CASE( test_lowNUHMSSMSemiAnalytic_ewsb_two_loop_solution )
{
   const double precision = 1.0e-5;
   lowNUHMSSMSemiAnalytic_input_parameters input;
   lowNUHMSSMSemiAnalytic_mass_eigenstates model(input);
   setup_lowNUHMSSMSemiAnalytic(model, input);

   const double susy_scale = 1000.;
   model.run_to(susy_scale);

   Boundary_values values;
   setup_susy_scale_lowNUHMSSMSemiAnalytic(model, values);

   lowNUHMSSMSemiAnalytic_semi_analytic_solutions solns;
   solns.set_input_scale(susy_scale);
   solns.set_output_scale(Electroweak_constants::MZ);

   solns.calculate_coefficients(model);

   model.run_to(Electroweak_constants::MZ);
   model.calculate_DRbar_masses();

   lowNUHMSSMSemiAnalytic_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_loop_order(2);
   ewsb_solver.set_semi_analytic_solutions(&solns);
   model.set_ewsb_solver(
      std::make_shared<lowNUHMSSMSemiAnalytic_ewsb_solver<Semi_analytic> >(ewsb_solver));

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

   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_1() - Re(tadpole_hh_1), 1.);
   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_2() - Re(tadpole_hh_2), 1.);
}

class SUSY_scale_ewsb_solver {
public:
   void set_loop_order(int l) { loop_order = l; }
   void set_number_of_iterations(int it) { number_of_iterations = it; }
   void set_precision(double p) { precision = p; }
   void set_susy_scale(double s) { susy_scale = s; }
   void set_ewsb_scale(double s) { ewsb_scale = s; }

   int solve_ewsb_at_susy_scale(lowNUHMSSMSemiAnalytic_mass_eigenstates& model) const;

private:
   int loop_order{2};
   int number_of_iterations{100};
   double precision{1.e-5};
   double susy_scale{1000.};
   double ewsb_scale{Electroweak_constants::MZ};
};

int SUSY_scale_ewsb_solver::solve_ewsb_at_susy_scale(lowNUHMSSMSemiAnalytic_mass_eigenstates& model) const
{
   auto m = model;
   m.set_ewsb_loop_order(loop_order);

   auto tadpole_stepper = [this, m](const Eigen::Matrix<double,2,1>& ewsb_pars) mutable -> Eigen::Matrix<double,2,1> {
      auto running_model = m;

      const double mHd20 = ewsb_pars(0);
      const double mHu20 = ewsb_pars(1);

      running_model.run_to(this->susy_scale);

      running_model.set_mHu2(mHu20);
      running_model.set_mHd2(mHd20);

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
   x_init(0) = model.get_mHd2();
   x_init(1) = model.get_mHu2();

   int status;
   for (auto& solver: solvers) {
      status = solver->solve(x_init);
      if (status == EWSB_solver::SUCCESS) {
         const auto solution = solver->get_solution();

         const double mHd20 = solution(0);
         const double mHu20 = solution(1);

         model.set_mHd20(mHd20);
         model.set_mHu20(mHu20);

         model.run_to(susy_scale);

         model.set_mHu2(mHu20);
         model.set_mHd2(mHd20);

         model.run_to(ewsb_scale);

         model.calculate_DRbar_masses();

         break;
      }
   }

   return status;
}

// this test checks that solving the EWSB equations using the
// semi-analytic solutions is equivalent to solving them by
// directly setting the parameters at the SUSY-scale
BOOST_AUTO_TEST_CASE( test_semi_analytic_ewsb_tree_level )
{
   const double precision = 1.0e-5;
   lowNUHMSSMSemiAnalytic_input_parameters input;
   lowNUHMSSMSemiAnalytic_mass_eigenstates model(input);
   setup_lowNUHMSSMSemiAnalytic(model, input);

   const double susy_scale = 1000.;
   model.run_to(susy_scale);

   Boundary_values values;
   setup_susy_scale_lowNUHMSSMSemiAnalytic(model, values);

   lowNUHMSSMSemiAnalytic_mass_eigenstates susy_scale_model(model);

   lowNUHMSSMSemiAnalytic_semi_analytic_solutions solns;
   solns.set_input_scale(susy_scale);
   solns.set_output_scale(Electroweak_constants::MZ);

   solns.calculate_coefficients(model);

   model.run_to(Electroweak_constants::MZ);

   lowNUHMSSMSemiAnalytic_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_loop_order(0);
   ewsb_solver.set_semi_analytic_solutions(&solns);
   model.set_ewsb_solver(
      std::make_shared<lowNUHMSSMSemiAnalytic_ewsb_solver<Semi_analytic> >(ewsb_solver));

   model.set_ewsb_iteration_precision(precision);
   const int error = model.solve_ewsb_tree_level();
   BOOST_CHECK_EQUAL(error, 0);

   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_1(), precision);
   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_2(), precision);

   const double mHd20_sol = model.get_mHd20();
   const double mHu20_sol = model.get_mHu20();

   SUSY_scale_ewsb_solver susy_scale_solver;
   susy_scale_solver.set_susy_scale(susy_scale);
   susy_scale_solver.set_ewsb_scale(Electroweak_constants::MZ);
   susy_scale_solver.set_precision(precision);
   susy_scale_solver.set_loop_order(0);

   const int susy_scale_error = susy_scale_solver.solve_ewsb_at_susy_scale(susy_scale_model);
   BOOST_CHECK_EQUAL(susy_scale_error, 0);

   BOOST_CHECK_SMALL(susy_scale_model.get_ewsb_eq_hh_1(), precision);
   BOOST_CHECK_SMALL(susy_scale_model.get_ewsb_eq_hh_2(), precision);

   BOOST_CHECK_CLOSE(mHd20_sol, susy_scale_model.get_mHd20(), 1.);
   BOOST_CHECK_CLOSE(mHu20_sol, susy_scale_model.get_mHu20(), 1.);
}

BOOST_AUTO_TEST_CASE( test_semi_analytic_ewsb_one_loop )
{
   const double precision = 1.0e-5;
   lowNUHMSSMSemiAnalytic_input_parameters input;
   lowNUHMSSMSemiAnalytic_mass_eigenstates model(input);
   setup_lowNUHMSSMSemiAnalytic(model, input);

   const double susy_scale = 1000.;
   model.run_to(susy_scale);

   Boundary_values values;
   setup_susy_scale_lowNUHMSSMSemiAnalytic(model, values);

   lowNUHMSSMSemiAnalytic_mass_eigenstates susy_scale_model(model);

   lowNUHMSSMSemiAnalytic_semi_analytic_solutions solns;
   solns.set_input_scale(susy_scale);
   solns.set_output_scale(Electroweak_constants::MZ);

   solns.calculate_coefficients(model);

   model.run_to(Electroweak_constants::MZ);
   model.calculate_DRbar_masses();

   lowNUHMSSMSemiAnalytic_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_loop_order(1);
   ewsb_solver.set_semi_analytic_solutions(&solns);
   model.set_ewsb_solver(
      std::make_shared<lowNUHMSSMSemiAnalytic_ewsb_solver<Semi_analytic> >(ewsb_solver));

   model.set_ewsb_iteration_precision(precision);
   const int error = model.solve_ewsb_one_loop();
   BOOST_CHECK_EQUAL(error, 0);

   const std::complex<double> tadpole_hh_1(model.tadpole_hh_1loop(0));
   const std::complex<double> tadpole_hh_2(model.tadpole_hh_1loop(1));

   BOOST_CHECK_SMALL(Im(tadpole_hh_1), 1.0e-12);
   BOOST_CHECK_SMALL(Im(tadpole_hh_2), 1.0e-12);

   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_1() - Re(tadpole_hh_1), 1.);
   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_2() - Re(tadpole_hh_2), 1.);

   const double mHd20_sol = model.get_mHd20();
   const double mHu20_sol = model.get_mHu20();

   SUSY_scale_ewsb_solver susy_scale_solver;
   susy_scale_solver.set_susy_scale(susy_scale);
   susy_scale_solver.set_ewsb_scale(Electroweak_constants::MZ);
   susy_scale_solver.set_precision(precision);
   susy_scale_solver.set_loop_order(1);

   const int susy_scale_error = susy_scale_solver.solve_ewsb_at_susy_scale(susy_scale_model);
   BOOST_CHECK_EQUAL(susy_scale_error, 0);

   const std::complex<double> running_tadpole_hh_1(susy_scale_model.tadpole_hh_1loop(0));
   const std::complex<double> running_tadpole_hh_2(susy_scale_model.tadpole_hh_1loop(1));

   BOOST_CHECK_SMALL(Im(running_tadpole_hh_1), 1.0e-12);
   BOOST_CHECK_SMALL(Im(running_tadpole_hh_2), 1.0e-12);

   BOOST_CHECK_SMALL(susy_scale_model.get_ewsb_eq_hh_1() - Re(running_tadpole_hh_1), 5.);
   BOOST_CHECK_SMALL(susy_scale_model.get_ewsb_eq_hh_2() - Re(running_tadpole_hh_2), 5.);

   BOOST_CHECK_CLOSE(mHd20_sol, susy_scale_model.get_mHd20(), 2.);
   BOOST_CHECK_CLOSE(mHu20_sol, susy_scale_model.get_mHu20(), 2.);
}
