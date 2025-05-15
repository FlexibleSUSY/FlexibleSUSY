
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_LQS1_FlexibleDecay

#include <iostream>

#include <boost/test/unit_test.hpp>

#define private public

#include <iomanip>

#include "test_LQS1.hpp"
#include "LQS1_up_basis_two_scale_model.hpp"
#include "decays/LQS1_up_basis_decays.hpp"
#include "decays/standard_model_decays.hpp"
#include "LQS1_up_basis_mass_eigenstates_running.hpp"
#include "decays/experimental_constraints.hpp"
#include "LQS1_up_basis_two_scale_spectrum_generator.hpp"

#include "lowe.h"
#include "loop_libraries/loop_library.hpp"

#include "LQS1_up_basis_mass_eigenstates.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_SM_FlexibleDecay )
{

   Loop_library::set(-1);

   LQS1_up_basis_input_parameters input;

   input.LambdaIN = 0.285;
   input.MS12Input = 2250000;
   input.gHS1Input = 0.1;

   LQS1_up_basis<Two_scale> _m;
   setup_LQS1_up_basis_const(_m, input);

   _m.calculate_DRbar_masses();

   _m.set_pole_mass_loop_order(1);
   _m.do_calculate_sm_pole_masses(true);
   _m.solve_ewsb_one_loop();
   _m.calculate_pole_masses();

   if (_m.get_problems().have_problem()) {
      std::ostringstream ostr;
      _m.get_problems().print_problems(ostr);
      BOOST_FAIL(ostr.str());
   }

   softsusy::QedQcd qedqcd;
   Physical_input physical_input;
   FlexibleDecay_settings flexibledecay_settings;
   LQS1_up_basis_mass_eigenstates_running m(_m);

   m.set_scale(100.0);

   const double p = 100.0;

   std::cout << std::setprecision(20) << '\n';
   std::cout << "Scale: " << m.get_scale() << '\n';

   std::cout << "MS1c: " << m.get_MS1c() << '\n';
   std::cout << "gHS1: " << m.get_gHS1() << '\n';
   std::cout << "lamQL: \n" << m.get_lamQL() << '\n';
   std::cout << "lamql: \n" << m.get_lamql() << '\n';


   // SM  printouts  of  physical  masses  and  tree - level  masses  and  parameters
   std::cout  << " Physical Masses \n";
   std::cout  << " mh: " << m.get_physical().Mhh << '\n';
   std::cout  << " mVZ: " << m.get_physical().MVZ << '\n';
   std::cout  << " mVW: " << m.get_physical().MVWp << '\n';
   std::cout  << " electron mass: " << m.get_physical().MFe(0) << '\n';
   std::cout  << " muon mass:  "  << m.get_physical().MFe(1) << '\n';
   std::cout  << " tau mass:  "  << m.get_physical().MFe(2) << '\n';
   std::cout  << "\n Tree-Level Masses \n";
   std::cout  << " Mh: " << m.get_Mhh() << '\n';
   std::cout  << " MZ: " << m.get_MVZ() << '\n';
   std::cout  << " MVW: " << m.get_MVWp() << '\n';
   std::cout  << " Ml: \n" << m.get_MFe() << '\n';
   // std::cout  << " Mv: \n" << dec_model.get() -> get_MFv() << ' \n';
   std::cout  << " MFu: \n" << m.get_MFu() << '\n';
   std::cout  << " MFd: \n" << m.get_MFd() << '\n';
   // std::cout  << " Ye: \n"  << dec_model.get() -> get_Ye() << ' \n';
   std::cout  << "\n Parameters \n";
   std::cout  << " g1: " << m.get_g1() << '\n';
   std::cout  << " g2: " << m.get_g2() << '\n';
   std::cout  << " v: " << m.get_v() << '\n';
   std::cout  << " TW: " << m.ThetaW() << '\n';

   // -----------------------------------------------------
   // decays with higher-order SM corrections

   LQS1_up_basis_decays decays_HO = LQS1_up_basis_decays(m, qedqcd, physical_input, flexibledecay_settings);
   decays_HO.sm.initialise_from_input(qedqcd);
   decays_HO.sm.set_loops(m.get_loops());
   decays_HO.sm.solve_ewsb_tree_level();
   decays_HO.sm.calculate_DRbar_masses();
   decays_HO.sm.calculate_pole_masses();
   decays_HO.sm_decays = Standard_model_decays(decays_HO.sm, qedqcd, physical_input, flexibledecay_settings);

   // ------------ tree-level decays ------------

   // h -> b bbar
   //BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_barFdFd(&m, 2, 2), 0.0023811031255194888, 2e-15);
   // h -> c cbar
   //BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_barFuFu(&m, 1, 1), 0.00011737301687946969, 2e-16);
   // h -> mu+ mu-
   BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_barFeFe(&m, 1, 1), 5.6230325675374227851e-7, 1e-15);
   // h -> W+ W-
   //BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_VWpconjVWp(&m), 0.00096256841980060836, 1e-3);
   // h -> Z Z
   //BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_VZVZ(&m), 0.00010568799794141996, 1e-3);

   // ------------ loop-induces decays_HO ------------

   // h -> gluon gluon
   //BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_VGVG(&m), 0.00035462447439152465, 5e-13);
   // h -> gamma gamma
   //BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_VPVP(&m), 9.2117697375801348e-06, 2e-13);
   // h -> gamma Z
   //BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_VPVZ(&m), 6.106402229544854e-06, 2e-13);

   // -----------------------------------------------------
   // decays without higher-order SM corrections

   flexibledecay_settings.set(FlexibleDecay_settings::include_higher_order_corrections, 0.0);
   LQS1_up_basis_decays decays_no_HO = LQS1_up_basis_decays(m, qedqcd, physical_input, flexibledecay_settings);

   // ------------ tree-level decays ------------

   // h -> b bbar
   //BOOST_CHECK_CLOSE_FRACTION(decays_no_HO.partial_width_hh_to_barFdFd(&m, 2, 2), 0.00207153178400001, 2e-15);
   // h -> c cbar
   //BOOST_CHECK_CLOSE_FRACTION(decays_no_HO.partial_width_hh_to_barFuFu(&m, 1, 1), 7.4579299427674489e-06, 1e-16);
   // h -> mu+ mu-
   BOOST_CHECK_CLOSE_FRACTION(decays_no_HO.partial_width_hh_to_barFeFe(&m, 1, 1), 9.077062783486042e-7, 1e-15);

   // ------------ loop-induces decays_no_HO ------------

   // h -> gluon gluon
   //BOOST_CHECK_CLOSE_FRACTION(decays_no_HO.partial_width_hh_to_VGVG(&m), 0.00019973935357968308, 5e-13);
   // h -> gamma gamma
   //BOOST_CHECK_CLOSE_FRACTION(decays_no_HO.partial_width_hh_to_VPVP(&m), 9.0492996379713826e-06, 2e-13);
   // h -> gamma Z
   //BOOST_CHECK_CLOSE_FRACTION(decays_no_HO.partial_width_hh_to_VPVZ(&m), 6.0803965974442744e-06, 2e-13);
}
