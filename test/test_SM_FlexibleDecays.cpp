
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SM_FlexibleDecays

#include <iostream>

#include <boost/test/unit_test.hpp>

#include "test_SM.hpp"
#include "SM_two_scale_model.hpp"
#include "SM_decays.hpp"

// #include "wrappers.hpp"
#include "lowe.h"
// #include "standard_model.hpp"
// TODO: remove before release
#include <iomanip>

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_SM_FlexibleDecays )
{
   
   SM_input_parameters input;
   input.LambdaIN = 0.285;
   SM<Two_scale> m;
   setup_SM_const(m, input);

   m.calculate_DRbar_masses();

   m.set_pole_mass_loop_order(1);
   m.do_calculate_sm_pole_masses(true);
   m.solve_ewsb_one_loop();
   m.calculate_pole_masses();

   if (m.get_problems().have_problem()) {
      std::ostringstream ostr;
      m.get_problems().print_problems(ostr);
      BOOST_FAIL(ostr.str());
   }

   softsusy::QedQcd qedqcd;
   SM_decays decays(m, qedqcd, input, true);

   // ------------ tree-level decays ----------------------------------------

   // h -> b bbar
   // no QED corrections
   BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_barFdFd(&m, 2, 2),
                              2.6118180765322455E-003, 2e-15);
   // QED corrections
   BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_barFdFd(&m, 2, 2),
                              2.6059181498481999E-003, 5e-15);
   // h -> tau+ tau-
   BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_barFeFe(&m, 2, 2),
                              2.6800741077127161E-004, 1e-15);
   // h -> W+ W-
   BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VWpconjVWp(&m),
                              8.4705126919250480E-004, 1e-14);
   // h -> Z Z
   BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VZVZ(&m),
                              9.4231401208556083E-005, 2e-14);

   // ------------ loop-induces decays --------------------------------------

   // h -> gluon gluon
   // BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VGVG(m), , 1e-15);
   // h -> gamma gamma
   // BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VPVP(m), , 1e-15);
   // h -> Z gamma
   // BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VZVP(m), , 1e-15);
}
