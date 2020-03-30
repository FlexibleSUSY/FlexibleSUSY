
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_shooting_solver

#include <boost/test/unit_test.hpp>

#include "shooting_solver.hpp"
#include "wrappers.hpp"

#include "mock_models.hpp"
#include "mock_single_scale_constraints.hpp"
#include "mock_single_scale_matchings.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_trival_matching )
{
   static const int N = 10;
   using Vec_t = Shooting_solver<N>::Vec_t;
   using Pre_t = Shooting_solver<N>::Pre_t;
   using Set_t = Shooting_solver<N>::Set_t;

   // parameters to be predicted
   Vec_t prediction_goal(Vec_t::Zero());
   for (int i = 0; i < N; ++i)
      prediction_goal(i) = i + 1;

   Static_model eft(Vec_t::Zero()), bsm(Vec_t::Zero());
   Counting_constraint ceft(100), cbsm(1000);

   // this trivial matching condition simply forwards the parameters
   // of one model to the other
   Trivial_matching_condition mc;

   Pre_t pred = [&] () {
      return Vec_t(prediction_goal - eft.get_parameters());
   };

   Set_t sett = [&] (const Vec_t& v) {
      bsm.set_parameters(v);
   };

   Shooting_solver<N> solver;
   solver.add(&cbsm, &bsm);
   solver.add(&mc, &bsm, &eft);
   solver.add(&ceft, &eft);
   solver.add(pred);
   solver.add(sett);

   const Vec_t start(Vec_t::Zero());

   BOOST_CHECK_NO_THROW(solver.solve(start));

   // the high scale parameters should be the same as the low scale
   // parameters
   BOOST_CHECK_LT((eft.get_parameters() - prediction_goal).cwiseAbs().sum(), 1e-10);
   BOOST_CHECK_LT((bsm.get_parameters() - prediction_goal).cwiseAbs().sum(), 1e-10);
}

