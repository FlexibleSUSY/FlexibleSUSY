// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_string_format

#include <boost/test/unit_test.hpp>
#include <limits>
#include "stopwatch.hpp"
#include "string_format.hpp"


BOOST_AUTO_TEST_CASE(test_int)
{
   BOOST_CHECK_EQUAL(flexiblesusy::to_string( 0),  "0");
   BOOST_CHECK_EQUAL(flexiblesusy::to_string( 1),  "1");
   BOOST_CHECK_EQUAL(flexiblesusy::to_string(-1), "-1");
}


BOOST_AUTO_TEST_CASE(test_double)
{
   const double qnan = std::numeric_limits<double>::quiet_NaN();
   const double snan = std::numeric_limits<double>::signaling_NaN();
   const double inf  = std::numeric_limits<double>::infinity();

   BOOST_CHECK_EQUAL(flexiblesusy::to_string( 0.0),   "0");
   BOOST_CHECK_EQUAL(flexiblesusy::to_string( 1.0),   "1");
   BOOST_CHECK_EQUAL(flexiblesusy::to_string(-1.0),  "-1");
   BOOST_CHECK_EQUAL(flexiblesusy::to_string(qnan), "nan");
   BOOST_CHECK_EQUAL(flexiblesusy::to_string(snan), "nan");
   BOOST_CHECK_EQUAL(flexiblesusy::to_string( inf), "inf");
}
