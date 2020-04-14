
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_fixed_point_iterator

#include <boost/test/unit_test.hpp>

#define ENABLE_VERBOSE 1
#undef ENABLE_DEBUG

#include "fixed_point_iterator.hpp"
#include "wrappers.hpp"

using namespace flexiblesusy;

typedef Eigen::Matrix<double,2,1> EV2_t;

static int number_of_calls = 0;
void reset() { number_of_calls = 0; }
int get_number_of_calls() { return number_of_calls; }

/**
 * Finding root of f(x,y) = ((x-5)^2, (y-1)^2) ,
 *
 * => Update steps
 *
 * (x,y) = (-25/(x-10), -1/(y-2))
 *
 * @param x touple (x,y)
 *
 * @return fixed point iteration update steps
 */
auto parabola1 = [](EV2_t x) -> EV2_t {
   const double y = x(0);
   const double z = x(1);
   EV2_t f;
   f(0) = -25./(y - 10.);
   f(1) = -1./(z - 2.);
   number_of_calls++;
   return f;
};

/**
 * Finding root of f(x,y) = ((x-5)^2, (y-1)^2) ,
 *
 * => Update steps in the problematic form
 *
 * (x,y) = ((y^2 + 25)/10, (z^2+1)/2)
 *
 * @param x touple (x,y)
 *
 * @return fixed point iteration update steps
 */
auto parabola2 = [](EV2_t x) -> EV2_t {
   const double y = x(0);
   const double z = x(1);
   EV2_t f;
   f(0) = (y*y + 25.) / 10.;
   f(1) = (z*z + 1.) / 2.;
   number_of_calls++;
   return f;
};

BOOST_AUTO_TEST_CASE( test_parabola1 )
{
   const double precision = 1.0e-4;
   Eigen::Matrix<double,2,1> start;
   start << 9, 9;
   Fixed_point_iterator<2> fpi(parabola1, 1000, fixed_point_iterator::Convergence_tester_relative<2>(precision));

   reset();

   int status = fpi.find_fixed_point(start);
   const auto fixed_point = fpi.get_solution();

   const double residual_1 = MaxRelDiff(5.0, fixed_point(0));
   const double residual_2 = MaxRelDiff(1.0, fixed_point(1));

   // Note: The convergence criterion
   // MaxRelDiff(x_{n+1}, x_{n}) < precision
   // is not very good: The method converges slowly.  This means
   // subsequent steps are very close to each other, but x_n might not
   // be close to the true fixed point.

   BOOST_REQUIRE(status == EWSB_solver::SUCCESS);
   BOOST_CHECK_LT(residual_1, 200*precision);
   BOOST_CHECK_LT(residual_2, 200*precision);
   BOOST_TEST_MESSAGE("fixed point iterator used " << get_number_of_calls() << " calls");
}

BOOST_AUTO_TEST_CASE( test_parabola2 )
{
   const double precision = 1.0e-4;
   Eigen::Matrix<double,2,1> start;
   start << 9, 9;
   Fixed_point_iterator<2> fpi(parabola2, 1000, fixed_point_iterator::Convergence_tester_relative<2>(precision));

   reset();

   int status = fpi.find_fixed_point(start);

   // The form of the update steps in Parabola2 is problematic,
   // because they are quadratic in the variables and are therefore
   // not small.

   BOOST_REQUIRE(status != EWSB_solver::SUCCESS);
}

/**
 * Update function which has the form of a constant plus a
 * perturbation term.
 */
auto perturbation = [](EV2_t x) -> EV2_t {
   const double y = x(0);
   const double z = x(1);
   const double f1 = 1 + (y - z*z)/(16.*Sqr(Pi));
   const double f2 = 2 + (y*y - z)/(16.*Sqr(Pi));
   EV2_t f;
   f(0) = f1;
   f(1) = f2;
   number_of_calls++;
   return f;
};

BOOST_AUTO_TEST_CASE( test_perturbation )
{
   const double precision = 1.0e-4;
   Eigen::Matrix<double,2,1> start;
   start << 10, 10;
   Fixed_point_iterator<2> fpi(perturbation, 1000, fixed_point_iterator::Convergence_tester_relative<2>(precision));

   reset();

   int status = fpi.find_fixed_point(start);
   const auto fixed_point = fpi.get_solution();

   BOOST_REQUIRE(status == EWSB_solver::SUCCESS);
   BOOST_CHECK_CLOSE_FRACTION(fixed_point(0), 1.0, 0.02);
   BOOST_CHECK_CLOSE_FRACTION(fixed_point(1), 2.0, 0.04);

   BOOST_TEST_MESSAGE("fixed point iterator used " << get_number_of_calls() << " calls");
   BOOST_CHECK(get_number_of_calls() < 6);
}
