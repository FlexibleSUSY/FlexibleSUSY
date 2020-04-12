
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_ewsb_solver

#include <boost/test/unit_test.hpp>

#include "ewsb_solver.hpp"
#include "root_finder.hpp"
#include "minimizer.hpp"
#include "fixed_point_iterator.hpp"
#include "wrappers.hpp"

using namespace flexiblesusy;

static int number_of_calls = 0;
static const double small = 1/(16 * Pi*Pi);

void reset() { number_of_calls = 0; }
int get_number_of_calls() { return number_of_calls; }

typedef Eigen::Matrix<double,2,1> EV2_t;

BOOST_AUTO_TEST_CASE( test_parabola_2dim )
{
   /**
    * Function f(x,y) = ((x-5)^2, (y-1)^2) ,
    *
    * @param x touple (x,y)
    * @return f(x,y)
    */
   auto func = [](const EV2_t& x) -> EV2_t {
      const double y = x(0);
      const double z = x(1);
      EV2_t f;
      f(0) = (y - 5.)*(y - 5.);
      f(1) = (z - 1.)*(z - 1.);
      number_of_calls++;
      return f;
   };

   /**
    * Function f(x,y) = ((x-5)^2, (y-1)^2) ,
    *
    * @param x touple (x,y)
    * @return f_1(x,y) + f_2(x,y)
    */
   auto chi2 = [](const EV2_t& x) -> double {
      const double y = x(0);
      const double z = x(1);
      const double chi2 = (y - 5.)*(y - 5.) + (z - 1.)*(z - 1.);
      number_of_calls++;
      return chi2;
   };

   /**
    * Update steps of f(x,y) = ((x-5)^2, (y-1)^2) :
    *
    * (x,y) = (-25/(x-10), -1/(y-2))
    *
    * @param x touple (x,y)
    *
    * @return fixed point iteration update steps
    */
   auto update = [](const EV2_t& x) -> EV2_t {
      const double y = x(0);
      const double z = x(1);
      EV2_t f;
      f(0) = -25./(y - 10.);
      f(1) = -1./(z - 2.);
      number_of_calls++;
      return f;
   };

   const int dimension = 2;
   const double precision = 1.0e-4;
   const int max_iterations = 1000;

   Eigen::Matrix<double,2,1> start;
   start << 9, 9;

   EWSB_solver* solvers[] = {
      new Root_finder<dimension>(func, max_iterations, precision, Root_finder<dimension>::GSLHybrid),
      new Root_finder<dimension>(func, max_iterations, precision, Root_finder<dimension>::GSLHybridS),
      new Root_finder<dimension>(func, max_iterations, precision, Root_finder<dimension>::GSLBroyden),
      new Root_finder<dimension>(func, max_iterations, precision, Root_finder<dimension>::GSLNewton),
      new Minimizer<dimension>(chi2, max_iterations, precision, Minimizer<dimension>::GSLSimplex),
      new Minimizer<dimension>(chi2, max_iterations, precision, Minimizer<dimension>::GSLSimplex2),
      new Minimizer<dimension>(chi2, max_iterations, precision, Minimizer<dimension>::GSLSimplex2Rand),
      new Fixed_point_iterator<dimension,fixed_point_iterator::Convergence_tester_relative<dimension>>(update, max_iterations, fixed_point_iterator::Convergence_tester_relative<dimension>(precision)),
      new Fixed_point_iterator<dimension,fixed_point_iterator::Convergence_tester_absolute<dimension>>(update, max_iterations, fixed_point_iterator::Convergence_tester_absolute<dimension>(precision))
   };

   for (int i = 0; i < sizeof(solvers)/sizeof(solvers[0]); i++) {
      reset();

      const int status = solvers[i]->solve(start);

      BOOST_CHECK(status == EWSB_solver::SUCCESS);

      BOOST_CHECK_CLOSE_FRACTION(solvers[i]->get_solution()(0), 5., 0.01);
      BOOST_CHECK_CLOSE_FRACTION(solvers[i]->get_solution()(1), 1., 0.01);

      BOOST_TEST_MESSAGE("solver " << i << ": "
                    << (status == EWSB_solver::SUCCESS ?
                        "success" : "fail")
                    << " ("
                    << get_number_of_calls() << " function calls)");
   }

   for (int i = 0; i < sizeof(solvers)/sizeof(solvers[0]); i++) {
      delete solvers[i];
   }
}

/**
 * More realistic example closer to EWSB eqs (including tadpoles)
 */

BOOST_AUTO_TEST_CASE( test_perturbation )
{
   auto func = [](const EV2_t& x) -> EV2_t {
      const double y = x(0);
      const double z = x(1);
      EV2_t f;
      f(0) =  y + z + small  *(y*y + z*z) + 1;
      f(1) = -y + z + small*2*(y*y + z*z) + 2;
      number_of_calls++;
      return f;
   };

   /**
    * Function f(x,y) = ((x-5)^2, (y-1)^2) ,
    *
    * @param x touple (x,y)
    * @return f_1(x,y) + f_2(x,y)
    */
   auto chi2 = [func](const EV2_t& x) -> double {
      EV2_t f = func(x);
      const double chi2 = Sqr(f(0)) + Sqr(f(1));
      number_of_calls++;
      return chi2;
   };

   /**
    * Update steps of f(x,y) = ((x-5)^2, (y-1)^2) :
    *
    * (x,y) = (-25/(x-10), -1/(y-2))
    *
    * @param x touple (x,y)
    *
    * @return fixed point iteration update steps
    */
   auto update = [](const EV2_t& x) -> EV2_t {
      const double y = x(0);
      const double z = x(1);
      EV2_t f;
      f(0) =  0.5 + small*0.5*(y*y + z*z);
      f(1) = -1.5 - small*1.5*(y*y + z*z);
      number_of_calls++;
      return f;
   };

   const int dimension = 2;
   const double precision = 1.0e-4;
   const int max_iterations = 1000;

   Eigen::Matrix<double,2,1> start;
   start << 10, 10;

   EWSB_solver* solvers[] = {
      new Root_finder<dimension>(func, max_iterations, precision, Root_finder<dimension>::GSLHybrid),
      new Root_finder<dimension>(func, max_iterations, precision, Root_finder<dimension>::GSLHybridS),
      new Root_finder<dimension>(func, max_iterations, precision, Root_finder<dimension>::GSLBroyden),
      new Root_finder<dimension>(func, max_iterations, precision, Root_finder<dimension>::GSLNewton),
      new Minimizer<dimension>(chi2, max_iterations, precision, Minimizer<dimension>::GSLSimplex),
      new Minimizer<dimension>(chi2, max_iterations, precision, Minimizer<dimension>::GSLSimplex2),
      new Minimizer<dimension>(chi2, max_iterations, precision, Minimizer<dimension>::GSLSimplex2Rand),
      new Fixed_point_iterator<dimension,fixed_point_iterator::Convergence_tester_relative<dimension>>(update, max_iterations, fixed_point_iterator::Convergence_tester_relative<dimension>(precision)),
      new Fixed_point_iterator<dimension,fixed_point_iterator::Convergence_tester_absolute<dimension>>(update, max_iterations, fixed_point_iterator::Convergence_tester_absolute<dimension>(precision)),
      new Fixed_point_iterator<dimension,fixed_point_iterator::Convergence_tester_tadpole<dimension>>(update, max_iterations, fixed_point_iterator::Convergence_tester_tadpole<dimension>(precision, func))
   };

   for (int i = 0; i < sizeof(solvers)/sizeof(solvers[0]); i++) {
      reset();

      const int status = solvers[i]->solve(start);

      BOOST_CHECK(status == EWSB_solver::SUCCESS);

      const auto solution = solvers[i]->get_solution();

      // true solution calculated with Mathematica 7
      BOOST_CHECK_CLOSE_FRACTION(solution(0), 0.508177, 0.01);
      BOOST_CHECK_CLOSE_FRACTION(solution(1), -1.52453, 0.01);

      BOOST_TEST_MESSAGE("solver " << i << ": "
                    << (status == EWSB_solver::SUCCESS ?
                        "success" : "fail")
                    << ", solution = (" << solution(0) << ", " << solution(1)
                    << ") ("
                    << get_number_of_calls() << " function calls)");
   }

   for (int i = 0; i < sizeof(solvers)/sizeof(solvers[0]); i++) {
      delete solvers[i];
   }
}
