#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_standard_model_mt_calculation

#include <boost/test/unit_test.hpp>

#include "lowe.h"
#include "standard_model.hpp"

using namespace flexiblesusy;

// tests (or estimates) the precision of the calculation
// mt -> yt -> mt
BOOST_AUTO_TEST_CASE( test_consistency )
{
   softsusy::QedQcd qedqcd;
   const double mt_pole_input = qedqcd.displayPoleMt();

   standard_model::Standard_model sm;
   sm.initialise_from_input(qedqcd);

   // #1 calculate mt_pole_1
   sm.calculate_MFu_pole();
   const double mt_pole_1 = sm.get_physical().MFu(2);

   // #2 calculate yt_msbar from mt_pole_1 in #1
   const double mt_msbar = sm.calculate_MFu_DRbar(mt_pole_1, 2);
   const double yt_msbar = std::sqrt(2.0)*mt_msbar/sm.get_v();
   sm.set_Yu(2, 2, yt_msbar);

   // #3 re-calculate mt_pole_2 from the new yt_msbar
   sm.calculate_MFu_pole();
   const double mt_pole_2 = sm.get_physical().MFu(2);

   // test equality between pole masses
   BOOST_CHECK_CLOSE_FRACTION(mt_pole_1, mt_pole_2, 0.01);
}
