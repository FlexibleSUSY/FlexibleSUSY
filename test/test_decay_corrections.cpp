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
#define BOOST_TEST_MODULE test_decay_corrections

#include <boost/test/unit_test.hpp>

#include "decays/decay_corrections.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE(test_default)
{
   const Decay_corrections dc;

   BOOST_CHECK_EQUAL(dc.sm, 4);
   BOOST_CHECK_EQUAL(dc.bsm, 1);
   BOOST_CHECK_EQUAL(dc.get(), 41);
}

BOOST_AUTO_TEST_CASE(test_set_get)
{
   for (int sm = 0; sm <= 4; sm++) {
      for (int bsm = 0; bsm <= 1; bsm++) {
         const Decay_corrections::Flags_t flags = 10 * sm + bsm;
         const Decay_corrections dc(flags);

         BOOST_CHECK_EQUAL(dc.sm, sm);
         BOOST_CHECK_EQUAL(dc.bsm, bsm);
         BOOST_CHECK_EQUAL(dc.get(), flags);
      }
   }
}
