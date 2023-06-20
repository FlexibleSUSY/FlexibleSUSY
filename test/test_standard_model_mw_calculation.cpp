
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SM_mw_calculation

#include <boost/test/unit_test.hpp>

#define protected public

#include "lowe.h"
#include "ew_input.hpp"
#include "standard_model.hpp"
//#include "SM_input_parameters.hpp"
//#include "SM_two_scale_spectrum_generator.hpp"
//#include "ew_input.hpp"
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
   sm.solve_ewsb_tree_level();
   sm.calculate_DRbar_masses();
   sm.solve_ewsb();
   sm.calculate_pole_masses();
   const double mw = sm.get_physical().MVWp;
   const double mh = sm.get_physical().Mhh;
   const double mwSM = calc_mw_SM(qedqcd, mh);
   BOOST_CHECK_CLOSE_FRACTION(mw, mwSM, 1.0e-10);
}
