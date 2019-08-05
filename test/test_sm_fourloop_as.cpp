#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_sm_fourloop_as

#include <boost/test/unit_test.hpp>
#include "sm_fourloop_as.hpp"

BOOST_AUTO_TEST_CASE( test_fixed_values )
{
   using namespace flexiblesusy::sm_fourloop_as;

   const double eps = 1e-5;
   const double pi = 3.1415926535897932384626433832795;
   const double mt = 173.34;
   const double alpha_s_SM5 = pi;
   const double nl = 5;
   const double nl2 = nl*nl;

   Parameters pars;
   pars.as = alpha_s_SM5;
   pars.mt = mt;
   pars.Q  = mt;

   const auto d1 = delta_alpha_s_1loop_as(pars);
   const auto d2 = delta_alpha_s_2loop_as_as(pars);
   const auto d3 = delta_alpha_s_3loop_as_as_as(pars);
   const auto d4 = delta_alpha_s_4loop_as_as_as_as(pars);

   // [hep-ph/0512060] Eq.(44)
   BOOST_CHECK_CLOSE_FRACTION(d1,  0.0, eps);
   BOOST_CHECK_CLOSE_FRACTION(d2, -0.152778, eps);
   BOOST_CHECK_CLOSE_FRACTION(d3, -0.972057 + 0.0846515*nl, eps);
   BOOST_CHECK_CLOSE_FRACTION(d4, -5.10032 + 1.00993*nl + 0.0219784*nl2, eps);
}
