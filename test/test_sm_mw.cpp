#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_sm_mw

#include <boost/test/unit_test.hpp>
#include "sm_mw.hpp"

BOOST_AUTO_TEST_CASE( test_sm_mw_fit_MSbar )
{
   using flexiblesusy::sm_mw::calculate_mw_pole_SM_fit_MSbar;

   const double mz = 91.1876; // GeV (Z boson pole mass)
   const double mh = 125.15; // GeV (SM Higgs boson pole mass)
   const double mt = 173.34; // GeV
   const double as = 0.1184;
   const double da5had = 0.02750;

   const auto res = calculate_mw_pole_SM_fit_MSbar(mh, mt, as, da5had);
   const double mw = res.first;
   const double dmw = res.second;

   BOOST_TEST_MESSAGE("mw = (" << mw << " +- " << dmw << ") GeV");

   BOOST_CHECK_CLOSE_FRACTION(mw, 80.357, dmw);
}

BOOST_AUTO_TEST_CASE( test_sm_mw_fit_OS )
{
   using flexiblesusy::sm_mw::calculate_mw_pole_SM_fit_OS;

   const double mz = 91.1875; // GeV (Z boson pole mass)
   const double mh = 114.4; // GeV (SM Higgs boson pole mass)
   const double mt = 173.3; // GeV
   const double as = 0.119;
   const double Da = 0.05907;

   const auto res = calculate_mw_pole_SM_fit_OS(mz, mh, mt, as, Da);
   const double mw = res.first;
   const double dmw = res.second;

   BOOST_TEST_MESSAGE("mw = (" << mw << " +- " << dmw << ") GeV");

   BOOST_CHECK_CLOSE_FRACTION(mw, 80.357, dmw);
}
