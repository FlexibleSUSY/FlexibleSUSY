#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_sm_fourloop_as

#include <boost/test/unit_test.hpp>
#include "sm_fourloop_as.hpp"
#include "wrappers.hpp"

BOOST_AUTO_TEST_CASE( test_fixed_values )
{
   using namespace flexiblesusy;
   using namespace flexiblesusy::sm_fourloop_as;

   const double eps = 5e-5;
   const double pi = 3.1415926535897932384626433832795;
   const double mt = 173.34;
   const double alpha_s_SM5 = pi;
   const double nl = 5;
   const double nl2 = nl*nl;

   Parameters pars;
   pars.as = alpha_s_SM5;
   pars.mt = mt;
   pars.Q  = mt;

   const auto dp1 = delta_alpha_s_1loop_as(pars);
   const auto dp2 = delta_alpha_s_2loop_as_as(pars);
   const auto dp3 = delta_alpha_s_3loop_as_as_as(pars);
   const auto dp4 = delta_alpha_s_4loop_as_as_as_as(pars);

   const auto d1 = - dp1;
   const auto d2 = - dp2 + 2.*Sqr(dp1);
   const auto d3 = - dp3 - 5.*Power3(dp1) + 5.*dp1*dp2;
   const auto d4 = - dp4 + 6.*dp1*dp3 + 3.*Power2(dp2) - 21.*Power2(dp1)*dp2 + 14.*Power4(dp1);

   // [hep-ph/0512060] Eq.(43)
   // Note: expressed in terms of as(nf), not as'(nl)
   BOOST_CHECK_CLOSE_FRACTION(d1, 0.0, eps);
   BOOST_CHECK_CLOSE_FRACTION(d2, 0.152778, eps);
   BOOST_CHECK_CLOSE_FRACTION(d3, 0.972057 - 0.0846515*nl, eps);
   BOOST_CHECK_CLOSE_FRACTION(d4, 5.17035 - 1.00993*nl - 0.0219784*nl2, eps);
}

BOOST_AUTO_TEST_CASE( test_fixed_values_inverse )
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

   const auto as6 = calc_alpha_s(pars, 4);

   BOOST_CHECK_CLOSE_FRACTION(as6, alpha_s_SM5 * (1.0 + d1 + d2 + d3 + d4), eps);
}

BOOST_AUTO_TEST_CASE( test_alternatives )
{
   using namespace flexiblesusy::sm_fourloop_as;

   const double eps = 1e-5;
   const double pi = 3.1415926535897932384626433832795;
   const double mt = 173.34;
   const double alpha_s_SM5 = 0.1187;

   Parameters pars;
   pars.as = alpha_s_SM5;
   pars.mt = mt;
   pars.Q  = mt;

   const auto as6     = calc_alpha_s(pars, 4);
   const auto as6_alt = calc_alpha_s_alternative(pars, 4);

   BOOST_TEST_MESSAGE("alpha_s       = " << as6);
   BOOST_TEST_MESSAGE("alpha_s (alt) = " << as6_alt);

   BOOST_CHECK_LT(std::abs(as6 - as6_alt)/as6, 1.0e-7);
}
