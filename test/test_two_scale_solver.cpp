#include "convergence_tester.hpp"
#include "single_scale_constraint.hpp"
#include "single_scale_matching.hpp"
#include "two_scale_solver.hpp"
#include "error.hpp"

#include "mock_convergence_testers.hpp"
#include "mock_models.hpp"
#include "mock_single_scale_constraints.hpp"
#include "mock_single_scale_matchings.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_two_scale_solver

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_trival_matching )
{
   Eigen::VectorXd parameters(10);
   const Eigen::VectorXd zeros(10);
   for (int i = 0; i < 10; ++i)
      parameters(i) = i;

   Counting_convergence_tester ccc(10);

   // init two models:
   // * low-energy model with parameters
   // * high-energy model with zeros
   Static_model model1(parameters), model2(zeros);
   Counting_constraint cc(1000);

   // this trivial matching condition simply forwards the parameters
   // of one model to the other
   Trivial_matching_condition mc;
   mc.set_models(&model1, &model2);

   RGFlow<Two_scale> solver;
   solver.add(&cc, &model1);
   solver.add(&mc, &model1, &model2);
   solver.add(&cc, &model2);
   solver.set_convergence_tester(&ccc);

   BOOST_CHECK_THROW(solver.solve(), NoConvergenceError);

   // the high scale parameters should be the same as the low scale
   // parameters
   BOOST_CHECK_EQUAL(model1.get_parameters(), parameters);
   BOOST_CHECK_EQUAL(model2.get_parameters(), parameters);
}

BOOST_AUTO_TEST_CASE( test_count_method_calls )
{
   for (int number_of_iterations = 0;
        number_of_iterations < 10; ++number_of_iterations ) {
      Counting_model model1, model2;
      Counting_constraint model1_c1(1000), model1_c2(2000);
      Counting_constraint model2_c1(4000), model2_c2(5000);
      Counting_matching_condition mc(3000);
      Counting_convergence_tester ccc(number_of_iterations);

      RGFlow<Two_scale> solver;
      solver.set_convergence_tester(&ccc);
      solver.add(&model1_c1, &model1);
      solver.add(&model1_c2, &model1);
      solver.add(&mc, &model1, &model2);
      solver.add(&model2_c1, &model2);
      solver.add(&model2_c2, &model2);
      solver.add(&model2_c1, &model2);
      solver.add(&mc, &model2, &model1);
      solver.add(&model1_c2, &model1);

      if (number_of_iterations == 0) {
         try {
            solver.solve();
         } catch (Error& e) {
            BOOST_ERROR(e.what());
         }
      } else {
         // expect NoConvergenceError because the accuracy_goal_reached()
         // function of the convergence tester returns always false
         BOOST_CHECK_THROW(solver.solve(), NoConvergenceError);
      }

      // check that all iterations were done
      BOOST_CHECK_EQUAL(solver.number_of_iterations_done(),
                        number_of_iterations);
      // check how often the matching is appied
      BOOST_CHECK_EQUAL(mc.get_number_of_matches(),
                        2*number_of_iterations);
      // lowest constraint applied first time and once in each iteration
      BOOST_CHECK_EQUAL(model1_c1.get_number_of_apply_calls(),
                        number_of_iterations);
      // highest constraint applied once in each iteration
      BOOST_CHECK_EQUAL(model2_c2.get_number_of_apply_calls(),
                        number_of_iterations);
      // intermediate constraints applied twice in each iteration
      BOOST_CHECK_EQUAL(model1_c2.get_number_of_apply_calls(),
                        2 * number_of_iterations);
      BOOST_CHECK_EQUAL(model2_c1.get_number_of_apply_calls(),
                        2 * number_of_iterations);
   }
}

BOOST_AUTO_TEST_CASE( test_run_to_with_zero_models )
{
   RGFlow<Two_scale> solver;

   int status = 0;
   try { solver.run_to(1000.); } catch (Error&) { status = 1; }

   BOOST_CHECK_EQUAL(status, 0);
   BOOST_CHECK_EQUAL(solver.get_model(), (void*)NULL);
}

BOOST_AUTO_TEST_CASE( test_run_to_with_one_model )
{
   Static_model model(Eigen::VectorXd(10));
   Counting_constraint cc(1000);
   RGFlow<Two_scale> solver;
   solver.add(&cc, &model);

   int status = 0;
   try { solver.run_to(1000.); } catch (Error&) { status = 1; }

   BOOST_CHECK_EQUAL(status, 0);
   BOOST_CHECK_EQUAL(solver.get_model(), &model);
}

BOOST_AUTO_TEST_CASE( test_run_to_with_two_models )
{
   Static_model model1(Eigen::VectorXd(10)), model2(Eigen::VectorXd(10));
   Counting_constraint c1(50), c2(200);
   Trivial_matching_condition mc(100);
   const double mc_scale = mc.get_scale();

   RGFlow<Two_scale> solver;
   solver.add(&c1, &model1);
   solver.add(&mc, &model1, &model2);
   solver.add(&c2, &model2);
   solver.add(&mc, &model2, &model1);

   solver.run_to(c1.get_scale());
   BOOST_CHECK_EQUAL(solver.get_model(), &model1);

   solver.run_to((c1.get_scale() + mc_scale) / 2.);
   BOOST_CHECK_EQUAL(solver.get_model(), &model1);

   solver.run_to(mc_scale);
   BOOST_CHECK_EQUAL(solver.get_model(), &model1);

   solver.run_to((mc_scale + c2.get_scale()) / 2.);
   BOOST_CHECK_EQUAL(solver.get_model(), &model2);

   solver.run_to(c2.get_scale());
   BOOST_CHECK_EQUAL(solver.get_model(), &model2);
}
