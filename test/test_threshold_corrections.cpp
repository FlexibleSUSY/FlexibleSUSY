#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_threshold_corrections

#include <boost/test/unit_test.hpp>
#include "threshold_corrections.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_get )
{
   Threshold_corrections tc;

   tc.alpha_em    = 1;
   tc.sin_theta_w = 2;
   tc.alpha_s     = 3;
   tc.mz          = 4;
   tc.mw          = 5;
   tc.mh          = 6;
   tc.mt          = 7;
   tc.mb          = 8;
   tc.mtau        = 9;

   BOOST_CHECK_EQUAL(tc.get(), 987654321);
}

BOOST_AUTO_TEST_CASE( test_get_negative )
{
   Threshold_corrections tc;

   tc.alpha_em    = -0;
   tc.sin_theta_w = -1;
   tc.alpha_s     = -2;
   tc.mz          = -3;
   tc.mw          = -4;
   tc.mh          = -5;
   tc.mt          = -6;
   tc.mb          = -7;
   tc.mtau        = -8;

   BOOST_CHECK_EQUAL(tc.get(), 0);
}

BOOST_AUTO_TEST_CASE( test_set )
{
   Threshold_corrections tc;
   tc.set(987654321);

   BOOST_CHECK_EQUAL(tc.alpha_em    , 1);
   BOOST_CHECK_EQUAL(tc.sin_theta_w , 2);
   BOOST_CHECK_EQUAL(tc.alpha_s     , 3);
   BOOST_CHECK_EQUAL(tc.mz          , 4);
   BOOST_CHECK_EQUAL(tc.mw          , 5);
   BOOST_CHECK_EQUAL(tc.mh          , 6);
   BOOST_CHECK_EQUAL(tc.mt          , 7);
   BOOST_CHECK_EQUAL(tc.mb          , 8);
   BOOST_CHECK_EQUAL(tc.mtau        , 9);
}

BOOST_AUTO_TEST_CASE( test_set_negative )
{
   Threshold_corrections tc;
   tc.set(-987654321);

   BOOST_CHECK_EQUAL(tc.alpha_em    , 1);
   BOOST_CHECK_EQUAL(tc.sin_theta_w , 2);
   BOOST_CHECK_EQUAL(tc.alpha_s     , 3);
   BOOST_CHECK_EQUAL(tc.mz          , 4);
   BOOST_CHECK_EQUAL(tc.mw          , 5);
   BOOST_CHECK_EQUAL(tc.mh          , 6);
   BOOST_CHECK_EQUAL(tc.mt          , 7);
   BOOST_CHECK_EQUAL(tc.mb          , 8);
   BOOST_CHECK_EQUAL(tc.mtau        , 9);
}
