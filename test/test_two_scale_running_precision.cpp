
#include "convergence_tester.hpp"
#include "model.hpp"
#include "single_scale_constraint.hpp"
#include "two_scale_solver.hpp"
#include "two_scale_running_precision.hpp"
#include "error.hpp"

#include "mock_convergence_testers.hpp"
#include "mock_models.hpp"
#include "mock_single_scale_constraints.hpp"

#include <cmath>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_two_scale_running_precision

#include <boost/test/unit_test.hpp>

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_constant_running_precision )
{
   for (int i = 0; i < 5; ++i) {
      const double default_precison = static_cast<double>(i);
      Two_scale_constant_precision const_prec(default_precison);

      // test that all iterations yield the same precison
      for (int iteration = 0; iteration < 10; ++iteration) {
         double precision = const_prec.get_precision(iteration);
         BOOST_CHECK_EQUAL(precision, default_precison);
      }
   }
}

BOOST_AUTO_TEST_CASE( test_increasing_running_precision )
{
   for (int i = 1; i < 10; ++i) {
      const double increasing_factor = static_cast<double>(i);
      const double min_precision = 1.0e-12;
      Two_scale_increasing_precision incr_prec(increasing_factor, min_precision);

      // test that all iterations yield expected precision
      for (int iteration = 0; iteration < 20; ++iteration) {
         double precision = incr_prec.get_precision(iteration);
         double expected_precision = std::pow(1./increasing_factor, iteration + 1);
         if (expected_precision < min_precision)
            expected_precision = min_precision;
         BOOST_CHECK_CLOSE(precision, expected_precision, 1.0e-5);
      }
   }
}

class Test_increasing_precision : public Two_scale_running_precision {
public:
   Test_increasing_precision() : call(0) {}
   virtual ~Test_increasing_precision() {}
   virtual double get_precision(int it) {
      // Check that every time this function is called the argument is
      // increased by 1, beginning with 0.
      BOOST_CHECK_EQUAL(call, it);
      ++call;
      return 1.0e-3;
   }
private:
   int call; ///< stores number of calls of get_precision()
};

BOOST_AUTO_TEST_CASE( test_increasing_iteration_number )
{
   // This checks that during each iteration (within
   // RGFlow<Two_scale>) the iteration number is increased and passed
   // to Two_scale_running_precision.
   //
   // The test in Test_increasing_precision ensures that
   // get_precision() is called exactly 1 time during each iteration.

   Static_model model;

   Static_constraint c1(100), c2(200);
   Static_convergence_tester ccc(10);
   Test_increasing_precision incr_prec;

   RGFlow<Two_scale> solver;
   solver.add(&c1, &model);
   solver.add(&c2, &model);
   solver.set_convergence_tester(&ccc);
   solver.set_running_precision(&incr_prec);

   BOOST_CHECK_THROW(solver.solve(), NoConvergenceError);
}
