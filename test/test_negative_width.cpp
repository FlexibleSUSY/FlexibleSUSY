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
#define BOOST_TEST_MODULE test_negative_width

#include <boost/test/unit_test.hpp>

#include "decays/decay.hpp"

#include <string>

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE(test_unphysical_negative_width)
{
   Decays_list decays(25);
   decays.set_decay(1., {5, -5}, "h -> b bbar");
   decays.set_decay(-1e-15, {22, 22}, "h -> gamma gamma"); // rounding noise
   decays.set_decay(-0.1, {23, 22}, "h -> Z gamma");       // unphysical

   // the total width is the sum of the absolute partial widths
   BOOST_CHECK_CLOSE_FRACTION(decays.get_total_width(), 1.1 + 1e-15, 1e-15);

   const double total = decays.get_total_width();
   BOOST_CHECK(
      !has_unphysical_negative_width(decays.get_decay({5, -5}), total));
   BOOST_CHECK(
      !has_unphysical_negative_width(decays.get_decay({22, 22}), total));
   BOOST_CHECK(
      has_unphysical_negative_width(decays.get_decay({23, 22}), total));

   const auto warnings = decays.remove_unphysical_negative_widths();
   BOOST_REQUIRE_EQUAL(warnings.size(), 1);
   BOOST_CHECK(warnings[0].find("h -> Z gamma") != std::string::npos);
   BOOST_CHECK(warnings[0].find("negative") != std::string::npos);

   // the unphysical channel is set to zero and removed from the total
   // width, the rounding noise is kept
   BOOST_CHECK_EQUAL(decays.get_decay({23, 22}).get_width(), 0.);
   BOOST_CHECK_EQUAL(decays.get_decay({22, 22}).get_width(), -1e-15);
   BOOST_CHECK_CLOSE_FRACTION(decays.get_total_width(), 1. + 1e-15, 1e-14);

   // nothing left to remove
   BOOST_CHECK(decays.remove_unphysical_negative_widths().empty());
}

BOOST_AUTO_TEST_CASE(test_only_unphysical_channel)
{
   Decays_list decays(1000025);
   decays.set_decay(-1e-8, {1000023, 22}, "Chi(3) -> Chi(2) VP");

   BOOST_CHECK_EQUAL(decays.remove_unphysical_negative_widths().size(), 1);
   BOOST_CHECK_SMALL(decays.get_total_width(), 1e-30);
}

BOOST_AUTO_TEST_CASE(test_negative_width_zero_total)
{
   const Decay decay(1000023, {1000022, 22}, -1e-20, "Chi(2) -> Chi(1) VP");

   // without a positive total width every negative width is unphysical
   BOOST_CHECK(has_unphysical_negative_width(decay, 0.));
}
