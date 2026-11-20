#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_derivative

#include <boost/test/unit_test.hpp>
#include <iomanip>

#include "derivative.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_derivative_x_sqr )
{
   auto f = [](double x) { return x*x; };
   auto df = [](double x) { return 2*x; };

   std::vector<double> x_values =
      {0., 1., 2., 3., 1e1, 1e2, 1e8, 1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12};

   for (auto x: x_values) {
      if (x == 0.) {
         BOOST_CHECK_SMALL(derivative_forward<0>(f,x), 1e-8);
         BOOST_CHECK_SMALL(derivative_forward<1>(f,x), 1e-9);
         BOOST_CHECK_SMALL(derivative_forward<2>(f,x), 1e-9);
         BOOST_CHECK_SMALL(derivative_forward<3>(f,x), 1e-9);
         BOOST_CHECK_SMALL(derivative_forward<4>(f,x), 1e-9);
         BOOST_CHECK_SMALL(derivative_forward<5>(f,x), 1e-9);
         BOOST_CHECK_SMALL(derivative_forward<6>(f,x), 1e-10);
         BOOST_CHECK_SMALL(derivative_forward<7>(f,x), 1e-10);
         BOOST_CHECK_SMALL(derivative_backward<0>(f,x), 1e-8);
         BOOST_CHECK_SMALL(derivative_backward<1>(f,x), 1e-9);
         BOOST_CHECK_SMALL(derivative_backward<2>(f,x), 1e-9);
         BOOST_CHECK_SMALL(derivative_backward<3>(f,x), 1e-9);
         BOOST_CHECK_SMALL(derivative_backward<4>(f,x), 1e-9);
         BOOST_CHECK_SMALL(derivative_backward<5>(f,x), 1e-9);
         BOOST_CHECK_SMALL(derivative_backward<6>(f,x), 1e-10);
         BOOST_CHECK_SMALL(derivative_backward<7>(f,x), 1e-10);
         BOOST_CHECK_SMALL(derivative_central<0>(f,x), 2e-9);
         BOOST_CHECK_SMALL(derivative_central<1>(f,x), 2e-9);
         BOOST_CHECK_SMALL(derivative_central<2>(f,x), 2e-9);
         BOOST_CHECK_SMALL(derivative_central<3>(f,x), 2e-9);
      } else {
         BOOST_CHECK_CLOSE_FRACTION(derivative_forward<0>(f,x), df(x), 1e-7);
         BOOST_CHECK_CLOSE_FRACTION(derivative_forward<1>(f,x), df(x), 1e-7);
         BOOST_CHECK_CLOSE_FRACTION(derivative_forward<2>(f,x), df(x), 1e-7);
         BOOST_CHECK_CLOSE_FRACTION(derivative_forward<3>(f,x), df(x), 1e-7);
         BOOST_CHECK_CLOSE_FRACTION(derivative_forward<4>(f,x), df(x), 1e-7);
         BOOST_CHECK_CLOSE_FRACTION(derivative_forward<5>(f,x), df(x), 1e-7);
         BOOST_CHECK_CLOSE_FRACTION(derivative_forward<6>(f,x), df(x), 2e-7);
         BOOST_CHECK_CLOSE_FRACTION(derivative_forward<7>(f,x), df(x), 3e-7);
         BOOST_CHECK_CLOSE_FRACTION(derivative_backward<0>(f,x), df(x), 1e-7);
         BOOST_CHECK_CLOSE_FRACTION(derivative_backward<1>(f,x), df(x), 1e-7);
         BOOST_CHECK_CLOSE_FRACTION(derivative_backward<2>(f,x), df(x), 1e-7);
         BOOST_CHECK_CLOSE_FRACTION(derivative_backward<3>(f,x), df(x), 1e-7);
         BOOST_CHECK_CLOSE_FRACTION(derivative_backward<4>(f,x), df(x), 1e-7);
         BOOST_CHECK_CLOSE_FRACTION(derivative_backward<5>(f,x), df(x), 1e-7);
         BOOST_CHECK_CLOSE_FRACTION(derivative_backward<6>(f,x), df(x), 2e-7);
         BOOST_CHECK_CLOSE_FRACTION(derivative_backward<7>(f,x), df(x), 2e-7);
         BOOST_CHECK_CLOSE_FRACTION(derivative_central<0>(f,x), df(x), 1e-8);
         BOOST_CHECK_CLOSE_FRACTION(derivative_central<1>(f,x), df(x), 1e-8);
         BOOST_CHECK_CLOSE_FRACTION(derivative_central<2>(f,x), df(x), 1e-8);
         BOOST_CHECK_CLOSE_FRACTION(derivative_central<3>(f,x), df(x), 1e-8);
      }
   }
}

BOOST_AUTO_TEST_CASE( test_eps_dependence )
{
   auto f = [](double x) { return x*x; };
   auto df = [](double x) { return 2*x; };

   const double x = 1.;

   std::vector<double> eps_values =
      {1e-16, 1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2};

   for (auto eps: eps_values) {
      const double df_true = df(x);
      const double df_forw0 = derivative_forward<0>(f,x,eps);
      const double df_forw1 = derivative_forward<1>(f,x,eps);
      const double df_forw2 = derivative_forward<2>(f,x,eps);
      const double df_forw3 = derivative_forward<3>(f,x,eps);
      const double df_forw4 = derivative_forward<4>(f,x,eps);
      const double df_forw5 = derivative_forward<5>(f,x,eps);
      const double df_forw6 = derivative_forward<6>(f,x,eps);
      const double df_forw7 = derivative_forward<7>(f,x,eps);
      const double df_back0 = derivative_backward<0>(f,x,eps);
      const double df_back1 = derivative_backward<1>(f,x,eps);
      const double df_back2 = derivative_backward<2>(f,x,eps);
      const double df_back3 = derivative_backward<3>(f,x,eps);
      const double df_back4 = derivative_backward<4>(f,x,eps);
      const double df_back5 = derivative_backward<5>(f,x,eps);
      const double df_back6 = derivative_backward<6>(f,x,eps);
      const double df_back7 = derivative_backward<7>(f,x,eps);
      const double df_cent0 = derivative_central<0>(f,x,eps);
      const double df_cent1 = derivative_central<1>(f,x,eps);
      const double df_cent2 = derivative_central<2>(f,x,eps);
      const double df_cent3 = derivative_central<3>(f,x,eps);

      printf("\n");
      printf("===========================================\n");
      printf("eps = % 18.16f, df = % 18.16f\n", eps, df_true);
      printf("-------------------------------------------\n");
      printf("forward[0]  - df = % 18.16f\n", df_forw0 - df_true);
      printf("forward[1]  - df = % 18.16f\n", df_forw1 - df_true);
      printf("forward[2]  - df = % 18.16f\n", df_forw2 - df_true);
      printf("forward[3]  - df = % 18.16f\n", df_forw3 - df_true);
      printf("forward[4]  - df = % 18.16f\n", df_forw4 - df_true);
      printf("forward[5]  - df = % 18.16f\n", df_forw5 - df_true);
      printf("forward[6]  - df = % 18.16f\n", df_forw6 - df_true);
      printf("forward[7]  - df = % 18.16f\n", df_forw7 - df_true);
      printf("backward[0] - df = % 18.16f\n", df_back0 - df_true);
      printf("backward[1] - df = % 18.16f\n", df_back1 - df_true);
      printf("backward[2] - df = % 18.16f\n", df_back2 - df_true);
      printf("backward[3] - df = % 18.16f\n", df_back3 - df_true);
      printf("backward[4] - df = % 18.16f\n", df_back4 - df_true);
      printf("backward[5] - df = % 18.16f\n", df_back5 - df_true);
      printf("backward[6] - df = % 18.16f\n", df_back6 - df_true);
      printf("backward[7] - df = % 18.16f\n", df_back7 - df_true);
      printf("central[0]  - df = % 18.16f\n", df_cent0 - df_true);
      printf("central[1]  - df = % 18.16f\n", df_cent1 - df_true);
      printf("central[2]  - df = % 18.16f\n", df_cent2 - df_true);
      printf("central[3]  - df = % 18.16f\n", df_cent3 - df_true);

      BOOST_CHECK_CLOSE_FRACTION(df_forw0, df_true, 6e-2);
      BOOST_CHECK_CLOSE_FRACTION(df_forw1, df_true, 5e-7);
      BOOST_CHECK_CLOSE_FRACTION(df_forw2, df_true, 5e-7);
      BOOST_CHECK_CLOSE_FRACTION(df_forw3, df_true, 5e-7);
      BOOST_CHECK_CLOSE_FRACTION(df_forw4, df_true, 5e-7);
      BOOST_CHECK_CLOSE_FRACTION(df_forw5, df_true, 5e-7);
      BOOST_CHECK_CLOSE_FRACTION(df_forw6, df_true, 5e-7);
      BOOST_CHECK_CLOSE_FRACTION(df_forw7, df_true, 5e-7);
      BOOST_CHECK_CLOSE_FRACTION(df_back0, df_true, 6e-2);
      BOOST_CHECK_CLOSE_FRACTION(df_back1, df_true, 5e-7);
      BOOST_CHECK_CLOSE_FRACTION(df_back2, df_true, 5e-7);
      BOOST_CHECK_CLOSE_FRACTION(df_back3, df_true, 5e-7);
      BOOST_CHECK_CLOSE_FRACTION(df_back4, df_true, 5e-7);
      BOOST_CHECK_CLOSE_FRACTION(df_back5, df_true, 5e-7);
      BOOST_CHECK_CLOSE_FRACTION(df_back6, df_true, 5e-7);
      BOOST_CHECK_CLOSE_FRACTION(df_back7, df_true, 1e-6);
      BOOST_CHECK_CLOSE_FRACTION(df_cent0, df_true, 2e-8);
      BOOST_CHECK_CLOSE_FRACTION(df_cent1, df_true, 2e-8);
      BOOST_CHECK_CLOSE_FRACTION(df_cent2, df_true, 1e-8);
      BOOST_CHECK_CLOSE_FRACTION(df_cent3, df_true, 1e-8);
   }
}
