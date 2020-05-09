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

#define BOOST_TEST_MODULE "Set looplibrary via environment"

#include <boost/test/unit_test.hpp>
#include "loop_libraries/loop_library.hpp"

namespace flexiblesusy
{

BOOST_AUTO_TEST_CASE(set_library) {
   Loop_library::set(-1);
   bool predicted_behavior = true;

   switch (Loop_library::get_type()) {
   case Loop_library::Library::Softsusy:
      BOOST_TEST_MESSAGE("lib<Softsusy>");
      break;
   case Loop_library::Library::Collier:
      BOOST_TEST_MESSAGE("lib<Collier>");
      break;
   case Loop_library::Library::Looptools:
      BOOST_TEST_MESSAGE("lib<Looptools>");
      break;
   case Loop_library::Library::Fflite:
      BOOST_TEST_MESSAGE("lib<Fflite>");
      break;
   default:
      BOOST_TEST_MESSAGE("lib<Oops>");
      predicted_behavior = false;
   }

   BOOST_CHECK(predicted_behavior);
}

} // flexiblesusy
