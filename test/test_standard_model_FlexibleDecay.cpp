
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SM_FlexibleDecay

#include <iostream>

#include <boost/test/unit_test.hpp>

#include "test_SM.hpp"
#include "decays/standard_model_decays.hpp"

// #include "wrappers.hpp"
#include "lowe.h"
// #include "standard_model.hpp"
#include "loop_libraries/loop_library.hpp"

using namespace flexiblesusy;
using namespace standard_model;

BOOST_AUTO_TEST_CASE( test_SM_FlexibleDecay )
{

   Loop_library::set(-1);

   softsusy::QedQcd qedqcd;
   Physical_input physical_input;
   FlexibleDecay_settings flexibledecay_settings;

   standard_model::Standard_model sm {};
   sm.initialise_from_input(qedqcd);
   sm.set_physical_input(physical_input);

   sm.calculate_pole_masses();

   // -----------------------------------------------------
   // decays with higher-order SM corrections

   Standard_model_decays decays_HO = Standard_model_decays(sm, qedqcd, physical_input, flexibledecay_settings);

   // ------------ tree-level decays ------------

   // h -> b bbar
   BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_barFdFd(sm, 2, 2),
                              0.0023811031255194888, 2e-15);
   // h -> c cbar
   BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_barFuFu(sm, 1, 1),
                              0.00011737301687946969, 2e-16);
   // h -> tau+ tau-
   BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_barFeFe(sm, 2, 2),
                              0.00026184531343741851, 1e-15);
   // h -> W+ W-
   // BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_conjVWpVWp(sm),
   //                            0.00088266545237025511, 1e-14);, 1e-14);
   BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_conjVWpVWp(sm),
                              0.00096256841980060836, 1e-3);
   // h -> Z Z
   // BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_VZVZ(sm),
   //                            8.4959557304996447e-05, 3e-14);
   BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_VZVZ(sm),
                              0.00010568799794141996, 1e-3);

   // ------------ loop-induces decays_HO ------------

   // h -> gluon gluon
   BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_VGVG(sm), 0.00035462447439152465, 5e-13);
   // h -> gamma gamma
   BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_VPVP(sm), 9.2117697375801348e-06, 2e-13);
   // h -> gamma Z
   BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_VPVZ(sm), 6.3322476114788634e-06, 2e-13);

   // -----------------------------------------------------
   // decays without higher-order SM corrections

   flexibledecay_settings.set(FlexibleDecay_settings::include_higher_order_corrections, 0.0);
   Standard_model_decays decays_no_HO = Standard_model_decays(sm, qedqcd, physical_input, flexibledecay_settings);

   // ------------ tree-level decays ------------

   // h -> b bbar
   BOOST_CHECK_CLOSE_FRACTION(decays_no_HO.partial_width_hh_to_barFdFd(sm, 2, 2),
                              0.00207153178400001, 2e-15);
   // h -> c cbar
   BOOST_CHECK_CLOSE_FRACTION(decays_no_HO.partial_width_hh_to_barFuFu(sm, 1, 1),
                              7.4579299427674489e-06, 1e-16);
   // h -> tau+ tau-
   BOOST_CHECK_CLOSE_FRACTION(decays_no_HO.partial_width_hh_to_barFeFe(sm, 2, 2),
                              0.00025910510645313849, 1e-15);

   // ------------ loop-induces decays_no_HO ------------

   // h -> gluon gluon
   BOOST_CHECK_CLOSE_FRACTION(decays_no_HO.partial_width_hh_to_VGVG(sm), 0.00019973935357968308, 5e-13);
   // h -> gamma gamma
   BOOST_CHECK_CLOSE_FRACTION(decays_no_HO.partial_width_hh_to_VPVP(sm), 9.0492996379713826e-06, 2e-13);
   // h -> gamma Z
   BOOST_CHECK_CLOSE_FRACTION(decays_no_HO.partial_width_hh_to_VPVZ(sm), 6.305199201144327e-06, 2e-13);
}
