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
#define BOOST_TEST_MODULE test_non_finite_width

#include <boost/test/unit_test.hpp>

#include "decays/decay.hpp"

#include <cmath>
#include <limits>
#include <string>

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE(test_nan_width)
{
   Decays_list decays(36);
   decays.set_decay(1., {5, -5}, "A -> b bbar");
   decays.set_decay(std::numeric_limits<double>::quiet_NaN(), {25, 23},
                    "A -> h Z");
   decays.set_decay(0.5, {6, -6}, "A -> t tbar");

   BOOST_CHECK(std::isnan(decays.get_total_width()));

   const auto problems = decays.remove_non_finite_widths();
   BOOST_REQUIRE_EQUAL(problems.size(), 1);
   BOOST_CHECK(problems[0].find("A -> h Z") != std::string::npos);
   BOOST_CHECK(problems[0].find("nan") != std::string::npos);

   // the channel is set to zero and the total width is recalculated
   BOOST_CHECK_EQUAL(decays.get_decay({25, 23}).get_width(), 0.);
   BOOST_CHECK_EQUAL(decays.get_decay({5, -5}).get_width(), 1.);
   BOOST_CHECK_CLOSE_FRACTION(decays.get_total_width(), 1.5, 1e-15);

   // nothing left to remove
   BOOST_CHECK(decays.remove_non_finite_widths().empty());
}

BOOST_AUTO_TEST_CASE(test_infinite_widths)
{
   Decays_list decays(36);
   decays.set_decay(std::numeric_limits<double>::infinity(), {25, 23},
                    "A -> h Z");
   decays.set_decay(-std::numeric_limits<double>::infinity(), {22, 22},
                    "A -> gamma gamma");
   decays.set_decay(-0.1, {21, 21}, "A -> g g");

   const auto problems = decays.remove_non_finite_widths();
   BOOST_REQUIRE_EQUAL(problems.size(), 2);
   BOOST_CHECK(problems[0].find("inf") != std::string::npos);
   BOOST_CHECK(problems[1].find("inf") != std::string::npos);

   // the finite negative channel is left for
   // remove_unphysical_negative_widths()
   BOOST_CHECK_CLOSE_FRACTION(decays.get_total_width(), 0.1, 1e-15);
   BOOST_CHECK_EQUAL(decays.remove_unphysical_negative_widths().size(), 1);
   BOOST_CHECK_SMALL(decays.get_total_width(), 1e-30);
}

BOOST_AUTO_TEST_CASE(test_finite_widths_unchanged)
{
   Decays_list decays(25);
   decays.set_decay(1., {5, -5}, "h -> b bbar");
   decays.set_decay(-1e-15, {22, 22}, "h -> gamma gamma");
   const double total = decays.get_total_width();

   BOOST_CHECK(decays.remove_non_finite_widths().empty());
   BOOST_CHECK_EQUAL(decays.get_total_width(), total);
   BOOST_CHECK_EQUAL(decays.get_decay({22, 22}).get_width(), -1e-15);
}
