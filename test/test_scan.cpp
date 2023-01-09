#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_scan

#include <boost/test/unit_test.hpp>
#include <cmath>
#include "scan.hpp"

BOOST_AUTO_TEST_CASE( test_float_range )
{
   {
      // empty range
      const auto range = flexiblesusy::float_range(0.0, 1.0, 0);
      BOOST_CHECK(range.empty());
   }
   {
      // range with equal interval boundaries
      const auto range = flexiblesusy::float_range(1.0, 1.0, 4);
      BOOST_CHECK_EQUAL(range.at(0), 1.0);
      BOOST_CHECK_EQUAL(range.at(1), 1.0);
      BOOST_CHECK_EQUAL(range.at(2), 1.0);
      BOOST_CHECK_EQUAL(range.at(3), 1.0);
      BOOST_CHECK_EQUAL(range.size(), 4);
   }
   {
      // range with inverted interval boundaries
      const auto range = flexiblesusy::float_range(1.0, 0.0, 4);
      BOOST_CHECK_EQUAL(range.at(0), 1.0);
      BOOST_CHECK_EQUAL(range.at(1), 0.75);
      BOOST_CHECK_EQUAL(range.at(2), 0.5);
      BOOST_CHECK_EQUAL(range.at(3), 0.25);
      BOOST_CHECK_EQUAL(range.size(), 4);
   }
   {
      // range with sorted interval boundaries
      const auto range = flexiblesusy::float_range(0.0, 1.0, 4);
      BOOST_CHECK_EQUAL(range.at(0), 0.0);
      BOOST_CHECK_EQUAL(range.at(1), 0.25);
      BOOST_CHECK_EQUAL(range.at(2), 0.5);
      BOOST_CHECK_EQUAL(range.at(3), 0.75);
      BOOST_CHECK_EQUAL(range.size(), 4);
   }
   {
      // range with a negative interval boundary
      const auto range = flexiblesusy::float_range(-1.0, 0.0, 4);
      BOOST_CHECK_EQUAL(range.at(0), -1.0);
      BOOST_CHECK_EQUAL(range.at(1), -0.75);
      BOOST_CHECK_EQUAL(range.at(2), -0.5);
      BOOST_CHECK_EQUAL(range.at(3), -0.25);
      BOOST_CHECK_EQUAL(range.size(), 4);
   }
   {
      // range with a negative interval boundary (unsorted)
      const auto range = flexiblesusy::float_range(0.0, -1.0, 4);
      BOOST_CHECK_EQUAL(range.at(0), 0.0);
      BOOST_CHECK_EQUAL(range.at(1), -0.25);
      BOOST_CHECK_EQUAL(range.at(2), -0.5);
      BOOST_CHECK_EQUAL(range.at(3), -0.75);
      BOOST_CHECK_EQUAL(range.size(), 4);
   }
}

BOOST_AUTO_TEST_CASE( test_float_range_log )
{
   static constexpr double e = 2.7182818284590452;
   {
      // empty range
      const auto range = flexiblesusy::float_range_log(1.0, 1.0, 0);
      BOOST_CHECK(range.empty());
   }
   {
      // range with equal interval boundaries
      const auto range = flexiblesusy::float_range_log(1.0, 1.0, 4);
      BOOST_CHECK_EQUAL(range.at(0), 1.0);
      BOOST_CHECK_EQUAL(range.at(1), 1.0);
      BOOST_CHECK_EQUAL(range.at(2), 1.0);
      BOOST_CHECK_EQUAL(range.at(3), 1.0);
      BOOST_CHECK_EQUAL(range.size(), 4);
   }
   {
      // range with sorted interval boundaries
      const auto range = flexiblesusy::float_range_log(1.0, e, 4);
      const double step_size = 0.25; // (log(e) - log(1))/4
      const double eps = 1e-15;
      BOOST_CHECK_CLOSE_FRACTION(std::log(range.at(0)), 0*step_size, eps);
      BOOST_CHECK_CLOSE_FRACTION(std::log(range.at(1)), 1*step_size, eps);
      BOOST_CHECK_CLOSE_FRACTION(std::log(range.at(2)), 2*step_size, eps);
      BOOST_CHECK_CLOSE_FRACTION(std::log(range.at(3)), 3*step_size, eps);
      BOOST_CHECK_EQUAL(range.size(), 4);
   }
   {
      // range with inverted interval boundaries
      const auto range = flexiblesusy::float_range_log(e, 1.0, 4);
      const double step_size = 0.25; // (log(e) - log(1))/4
      const double eps = 1e-15;
      BOOST_CHECK_CLOSE_FRACTION(std::log(range.at(0)), 4*step_size, eps);
      BOOST_CHECK_CLOSE_FRACTION(std::log(range.at(1)), 3*step_size, eps);
      BOOST_CHECK_CLOSE_FRACTION(std::log(range.at(2)), 2*step_size, eps);
      BOOST_CHECK_CLOSE_FRACTION(std::log(range.at(3)), 1*step_size, eps);
      BOOST_CHECK_EQUAL(range.size(), 4);
   }
}
