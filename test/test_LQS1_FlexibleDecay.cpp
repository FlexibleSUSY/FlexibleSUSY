
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

   LQS1_up_basis<Two_scale> m_;
   setup_LQS1_up_basis_const(m_, input);

   m_.calculate_DRbar_masses();

   m_.set_pole_mass_loop_order(1);
   m_.do_calculate_sm_pole_masses(true);
   m_.solve_ewsb_one_loop();
   m_.calculate_pole_masses();

   if (m_.get_problems().have_problem()) {
      std::ostringstream ostr;
      m_.get_problems().print_problems(ostr);
      BOOST_FAIL(ostr.str());
   }

   softsusy::QedQcd qedqcd;
   Physical_input physical_input;
   FlexibleDecay_settings flexibledecay_settings;
   LQS1_up_basis_mass_eigenstates_running m(m_);

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

   // tree-level
   flexibledecay_settings.set(FlexibleDecay_settings::include_higher_order_corrections, 0.0);
   LQS1_up_basis_decays decays_noHO = LQS1_up_basis_decays(m, qedqcd, physical_input, flexibledecay_settings);
   BOOST_CHECK_CLOSE_FRACTION(decays_noHO.partial_width_hh_to_barFeFe(&m, 1, 1), 5.6230325675374227851e-7, 1e-15);

   // 1l QED
   flexibledecay_settings.set(FlexibleDecay_settings::include_higher_order_corrections, 10.0);
   LQS1_up_basis_decays decays_HO_SM = LQS1_up_basis_decays(m, qedqcd, physical_input, flexibledecay_settings);
   BOOST_CHECK_CLOSE_FRACTION(decays_HO_SM.partial_width_hh_to_barFeFe(&m, 1, 1), 5.6230325675374227851e-7, 1e-15);

   // 1l BSM
   flexibledecay_settings.set(FlexibleDecay_settings::include_higher_order_corrections, 1.0);
   LQS1_up_basis_decays decays_HO_BSM = LQS1_up_basis_decays(m, qedqcd, physical_input, flexibledecay_settings);
   decays_HO_BSM.sm.initialise_from_input(qedqcd);
   decays_HO_BSM.sm.set_loops(m.get_loops());
   decays_HO_BSM.sm.solve_ewsb_tree_level();
   decays_HO_BSM.sm.calculate_DRbar_masses();
   decays_HO_BSM.sm.calculate_pole_masses();
   decays_HO_BSM.sm_decays = Standard_model_decays(decays_HO_BSM.sm, qedqcd, physical_input, flexibledecay_settings);
   BOOST_CHECK_CLOSE_FRACTION(decays_HO_BSM.partial_width_hh_to_barFeFe(&m, 1, 1), 5.6230325675374227851e-7, 1e-15);

   // 1l SM+BSM
   flexibledecay_settings.set(FlexibleDecay_settings::include_higher_order_corrections, 11.0);
   LQS1_up_basis_decays decays_HO_SM_BSM = LQS1_up_basis_decays(m, qedqcd, physical_input, flexibledecay_settings);
   decays_HO_SM_BSM.sm.initialise_from_input(qedqcd);
   decays_HO_SM_BSM.sm.set_loops(m.get_loops());
   decays_HO_SM_BSM.sm.solve_ewsb_tree_level();
   decays_HO_SM_BSM.sm.calculate_DRbar_masses();
   decays_HO_SM_BSM.sm.calculate_pole_masses();
   decays_HO_SM_BSM.sm_decays = Standard_model_decays(decays_HO_SM_BSM.sm, qedqcd, physical_input, flexibledecay_settings);
   BOOST_CHECK_CLOSE_FRACTION(decays_HO_SM_BSM.partial_width_hh_to_barFeFe(&m, 1, 1), 5.6230325675374227851e-7, 1e-15);
}
