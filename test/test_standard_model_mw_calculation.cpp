
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_standard_model_mw_calculation

#include <boost/test/unit_test.hpp>

#define protected public

#include "lowe.h"
#include "ew_input.hpp"
#include "standard_model.hpp"
#include "sm_mw.hpp"

using namespace flexiblesusy;

double calc_mw_SM(softsusy::QedQcd const& qedqcd, double mh)
{
   using flexiblesusy::sm_mw::calculate_mw_pole_SM_fit_MSbar;

   const auto res = calculate_mw_pole_SM_fit_MSbar(
      mh,
      qedqcd.displayPoleMt(),
      qedqcd.displayAlphaSInput(),
      Electroweak_constants::delta_alpha_s_5_had);

   return res.first;
}

BOOST_AUTO_TEST_CASE( test_consistency )
{
   softsusy::QedQcd qedqcd;
   standard_model::Standard_model sm{};
   sm.initialise_from_input(qedqcd);
   sm.calculate_pole_masses();
   const double mw = sm.get_physical().MVWp;
   const double mh = sm.get_physical().Mhh;
   const double mwSM = calc_mw_SM(qedqcd, mh);
   BOOST_CHECK_CLOSE_FRACTION(mw, mwSM, 1.0e-16);
}
