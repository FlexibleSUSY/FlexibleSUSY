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
#define BOOST_TEST_MODULE test_two_loop_functions

#include <boost/test/unit_test.hpp>
#include "two_loop_functions.hpp"

using namespace flexiblesusy;
using namespace flexiblesusy::two_loop_functions;

BOOST_AUTO_TEST_CASE( test_TfSS )
{
   BOOST_CHECK_EQUAL(TfSS(0,0,0,1), 0.);
}

BOOST_AUTO_TEST_CASE( test_call )
{
   volatile double out;

   out = TfSS(0, 0, 0, 0);
   out = TfSSS(0, 0, 0, 0);
   out = TfSSSS(0, 0, 0, 0, 0);
   out = TfSSFF(0, 0, 0, 0, 0);
   out = TfSSFbFb(0, 0, 0, 0, 0);
   out = TfFFFbS(0, 0, 0, 0, 0);
   out = TfFFbFS(0, 0, 0, 0, 0);
   out = TfFbFbFbS(0, 0, 0, 0, 0);
   out = TfSV(0, 0);
   out = TfFV(0, 0);

   out = WfSSSS(0, 0, 0, 0, 0, 0);
   out = XfSSS(0, 0, 0, 0, 0);
   out = YfSSSS(0, 0, 0, 0, 0, 0);
   out = ZfSSSS(0, 0, 0, 0, 0, 0);
   out = SfSSS(0, 0, 0, 0, 0);
   out = UfSSSS(0, 0, 0, 0, 0, 0);
   out = VfSSSSS(0, 0, 0, 0, 0, 0, 0);
   out = MfSSSSS(0, 0, 0, 0, 0, 0, 0);
   out = WfSSSV(0, 0, 0);
   out = MfSSSSV(0, 0, 0, 0);
   out = WfSSFF(0, 0, 0, 0, 0, 0);
   out = WfSSFbFb(0, 0, 0, 0, 0, 0);
   out = MfFbFbFbFbS(0, 0, 0, 0, 0, 0, 0);
   out = MfFFbFbFS(0, 0, 0, 0, 0, 0, 0);
   out = MfFFbFFbS(0, 0, 0, 0, 0, 0, 0);
   out = MfFFFbFbS(0, 0, 0, 0, 0, 0, 0);
   out = MfFFFFS(0, 0, 0, 0, 0, 0, 0);
   out = MfSFbSFbFb(0, 0, 0, 0, 0, 0, 0);
   out = MfSFSFbF(0, 0, 0, 0, 0, 0, 0);
   out = MfSFSFFb(0, 0, 0, 0, 0, 0, 0);
   out = VfSSSFbFb(0, 0, 0, 0, 0, 0, 0);
   out = VfSSSFF(0, 0, 0, 0, 0, 0, 0);
   out = VfFbFbFbFbS(0, 0, 0, 0, 0, 0, 0);
   out = VfFbFFbFS(0, 0, 0, 0, 0, 0, 0);
   out = VfFbFFFbS(0, 0, 0, 0, 0, 0, 0);
   out = VfFFbFbFS(0, 0, 0, 0, 0, 0, 0);
   out = VfFFFbFbS(0, 0, 0, 0, 0, 0, 0);
   out = VfFFFFS(0, 0, 0, 0, 0, 0, 0);
   out = GfFFV(0, 0, 0, 0);
   out = GfFbFbV(0, 0, 0, 0);
}
