
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SM_FlexibleDecays

#include <iostream>

#include <boost/test/unit_test.hpp>

#include "test_SM.hpp"
#include "SM_two_scale_model.hpp"
#include "decays/SM_decays.hpp"

// #include "wrappers.hpp"
#include "lowe.h"
// #include "standard_model.hpp"
// TODO: remove before release
#include <iomanip>
#include "loop_libraries/loop_library.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_SM_FlexibleDecays )
{

   Loop_library::set(-1);

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

   // -----------------------------------------------------
   // decays with higher-order SM corrections

   SM_decays decays_HO = SM_decays(m, qedqcd, SM_higher_order_corrections::enable);

   // ------------ tree-level decays ------------

   // h -> b bbar
   BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_barFdFd(&m, 2, 2),
                              0.0026142076883103791, 2e-15);
   // h -> c cbar
   BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_barFuFu(&m, 1, 1),
                              7.9883035101449685e-06, 1e-16);
   // h -> tau+ tau-
   BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_barFeFe(&m, 2, 2),
                              0.00027084177013049477, 1e-15);
   // h -> W+ W-
   BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_conjVWpVWp(&m),
                              0.00088266545237025511, 1e-14);
   // h -> Z Z
   BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_VZVZ(&m),
                              8.4959557304996447e-05, 3e-14);

   // ------------ loop-induces decays_HO ------------

   // h -> gluon gluon
   BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_VGVG(&m), 0.0003937446137366738, 5e-13);
   // h -> gamma gamma
   BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_VPVP(&m), 1.0566577997631895e-05, 2e-13);
   // h -> gamma Z
   BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_VPVZ(&m), 6.9031590236381169e-06, 2e-13);

   // -----------------------------------------------------
   // decays without higher-order SM corrections

   SM_decays decays_no_HO = SM_decays(m, qedqcd, SM_higher_order_corrections::disable);

   // ------------ tree-level decays ------------

   // h -> b bbar
   BOOST_CHECK_CLOSE_FRACTION(decays_no_HO.partial_width_hh_to_barFdFd(&m, 2, 2),
                              0.00207153178400001, 2e-15);
   // h -> c cbar
   BOOST_CHECK_CLOSE_FRACTION(decays_no_HO.partial_width_hh_to_barFuFu(&m, 1, 1),
                              7.9883035101449685e-06, 1e-16);
   // h -> tau+ tau-
   BOOST_CHECK_CLOSE_FRACTION(decays_no_HO.partial_width_hh_to_barFeFe(&m, 2, 2),
                              0.00026800741537194096, 1e-15);

   // ------------ loop-induces decays_no_HO ------------

   // h -> gluon gluon
   BOOST_CHECK_CLOSE_FRACTION(decays_no_HO.partial_width_hh_to_VGVG(&m), 0.00019973935357968308, 5e-13);
   // h -> gamma gamma
   BOOST_CHECK_CLOSE_FRACTION(decays_no_HO.partial_width_hh_to_VPVP(&m), 1.0385346434891975e-05, 2e-13);
   // h -> gamma Z
   BOOST_CHECK_CLOSE_FRACTION(decays_no_HO.partial_width_hh_to_VPVZ(&m), 6.9031590236381169e-06, 2e-13);
}
