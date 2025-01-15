
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SM_FlexibleDecay

#include <iostream>

#include <boost/test/unit_test.hpp>

#include "test_SM.hpp"
#include "SM_two_scale_model.hpp"
#include "decays/SM_decays.hpp"
#include "decays/standard_model_decays.hpp"
#include "decays/experimental_constraints.hpp"
#include "SM_two_scale_spectrum_generator.hpp"

#include "lowe.h"
#include "loop_libraries/loop_library.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_SM_FlexibleDecay )
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
   Physical_input physical_input;
   FlexibleDecay_settings flexibledecay_settings;

   // -----------------------------------------------------
   // decays with higher-order SM corrections

   SM_decays decays_HO = SM_decays(m, qedqcd, physical_input, flexibledecay_settings);

   // ------------ tree-level decays ------------

   // h -> b bbar
   BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_barFdFd(&m, 2, 2),
                              0.0023811031255194888, 2e-15);
   // h -> c cbar
   BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_barFuFu(&m, 1, 1),
                              0.00011737301687946969, 2e-16);
   // h -> tau+ tau-
   BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_barFeFe(&m, 2, 2),
                              0.00026184531343741851, 1e-15);
   // h -> W+ W-
   BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_VWpconjVWp(&m),
                              0.00096256841980060836, 1e-3);
   // h -> Z Z
   BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_VZVZ(&m),
                              0.00010568799794141996, 1e-3);

   // ------------ loop-induces decays_HO ------------

   // h -> gluon gluon
   BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_VGVG(&m), 0.00035462447439152465, 5e-13);
   // h -> gamma gamma
   BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_VPVP(&m), 9.2117697375801348e-06, 2e-13);
   // h -> gamma Z
   BOOST_CHECK_CLOSE_FRACTION(decays_HO.partial_width_hh_to_VPVZ(&m), 6.106402229544854e-06, 2e-13);

   // -----------------------------------------------------
   // decays without higher-order SM corrections

   flexibledecay_settings.set(FlexibleDecay_settings::include_higher_order_corrections, 0.0);
   SM_decays decays_no_HO = SM_decays(m, qedqcd, physical_input, flexibledecay_settings);

   // ------------ tree-level decays ------------

   // h -> b bbar
   BOOST_CHECK_CLOSE_FRACTION(decays_no_HO.partial_width_hh_to_barFdFd(&m, 2, 2),
                              0.00207153178400001, 2e-15);
   // h -> c cbar
   BOOST_CHECK_CLOSE_FRACTION(decays_no_HO.partial_width_hh_to_barFuFu(&m, 1, 1),
                              7.4579299427674489e-06, 1e-16);
   // h -> tau+ tau-
   BOOST_CHECK_CLOSE_FRACTION(decays_no_HO.partial_width_hh_to_barFeFe(&m, 2, 2),
                              0.00025910510645313849, 1e-15);

   // ------------ loop-induces decays_no_HO ------------

   // h -> gluon gluon
   BOOST_CHECK_CLOSE_FRACTION(decays_no_HO.partial_width_hh_to_VGVG(&m), 0.00019973935357968308, 5e-13);
   // h -> gamma gamma
   BOOST_CHECK_CLOSE_FRACTION(decays_no_HO.partial_width_hh_to_VPVP(&m), 9.0492996379713826e-06, 2e-13);
   // h -> gamma Z
   BOOST_CHECK_CLOSE_FRACTION(decays_no_HO.partial_width_hh_to_VPVZ(&m), 6.0803965974442744e-06, 2e-13);
}

BOOST_AUTO_TEST_CASE( test_SM_vs_standard_model_FlexibleDecay )
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

   standard_model::Standard_model sm{};
   setup_SM_const(sm, input);
   sm.calculate_DRbar_masses();
   sm.set_pole_mass_loop_order(1);
   sm.solve_ewsb_one_loop();
   sm.calculate_pole_masses();

   if (sm.get_problems().have_problem()) {
      std::ostringstream ostr;
      sm.get_problems().print_problems(ostr);
      BOOST_FAIL(ostr.str());
   }

   softsusy::QedQcd qedqcd;
   Physical_input physical_input;
   FlexibleDecay_settings flexibledecay_settings;

   SM_decays decays_SM_HO = SM_decays(m, qedqcd, physical_input, flexibledecay_settings);
   Standard_model_decays decays_sm_HO = Standard_model_decays(sm, qedqcd, physical_input, flexibledecay_settings);

   // ------------ tree-level decays ------------

   // h -> b bbar
   BOOST_CHECK_CLOSE_FRACTION(decays_SM_HO.partial_width_hh_to_barFdFd(&m, 2, 2),
                              decays_sm_HO.partial_width_hh_to_barFdFd(sm, 2, 2), 1e-16);
   // h -> c cbar
   BOOST_CHECK_CLOSE_FRACTION(decays_SM_HO.partial_width_hh_to_barFuFu(&m, 1, 1),
                              decays_sm_HO.partial_width_hh_to_barFuFu(sm, 1, 1), 1e-16);
   // h -> tau+ tau-
   BOOST_CHECK_CLOSE_FRACTION(decays_SM_HO.partial_width_hh_to_barFeFe(&m, 2, 2),
                              decays_sm_HO.partial_width_hh_to_barFeFe(sm, 2, 2), 1e-16);
   // h -> W+ W-
   BOOST_CHECK_CLOSE_FRACTION(decays_SM_HO.partial_width_hh_to_VWpconjVWp(&m),
                              decays_sm_HO.partial_width_hh_to_conjVWpVWp(sm), 1e-16);
   // h -> Z Z
   BOOST_CHECK_CLOSE_FRACTION(decays_SM_HO.partial_width_hh_to_VZVZ(&m),
                              decays_sm_HO.partial_width_hh_to_VZVZ(sm), 1e-16);

   // ------------ loop-induces decays_HO ------------

   // h -> gluon gluon
   BOOST_CHECK_CLOSE_FRACTION(decays_SM_HO.partial_width_hh_to_VGVG(&m), decays_sm_HO.partial_width_hh_to_VGVG(sm), 1e-16);
   // h -> gamma gamma
   BOOST_CHECK_CLOSE_FRACTION(decays_SM_HO.partial_width_hh_to_VPVP(&m), decays_sm_HO.partial_width_hh_to_VPVP(sm), 1e-16);
   // h -> gamma Z
   BOOST_CHECK_CLOSE_FRACTION(decays_SM_HO.partial_width_hh_to_VPVZ(&m), decays_sm_HO.partial_width_hh_to_VPVZ(sm), 1e-16);

   // -----------------------------------------------------
   // decays without higher-order SM corrections

   flexibledecay_settings.set(FlexibleDecay_settings::include_higher_order_corrections, 0.0);
   SM_decays decays_SM_no_HO = SM_decays(m, qedqcd, physical_input, flexibledecay_settings);
   Standard_model_decays decays_sm_no_HO = Standard_model_decays(sm, qedqcd, physical_input, flexibledecay_settings);

   // ------------ tree-level decays ------------

   // h -> b bbar
   BOOST_CHECK_CLOSE_FRACTION(decays_SM_no_HO.partial_width_hh_to_barFdFd(&m, 2, 2),
                              decays_sm_no_HO.partial_width_hh_to_barFdFd(sm, 2, 2), 1e-16);
   // h -> c cbar
   BOOST_CHECK_CLOSE_FRACTION(decays_SM_no_HO.partial_width_hh_to_barFuFu(&m, 1, 1),
                              decays_sm_no_HO.partial_width_hh_to_barFuFu(sm, 1, 1), 1e-16);
   // h -> tau+ tau-
   BOOST_CHECK_CLOSE_FRACTION(decays_SM_no_HO.partial_width_hh_to_barFeFe(&m, 2, 2),
                              decays_sm_no_HO.partial_width_hh_to_barFeFe(sm, 2, 2), 1e-16);

   // ------------ loop-induces decays_no_HO ------------

   // h -> gluon gluon
   BOOST_CHECK_CLOSE_FRACTION(decays_SM_no_HO.partial_width_hh_to_VGVG(&m), decays_sm_no_HO.partial_width_hh_to_VGVG(sm), 1e-16);
   // h -> gamma gamma
   BOOST_CHECK_CLOSE_FRACTION(decays_SM_no_HO.partial_width_hh_to_VPVP(&m), decays_sm_no_HO.partial_width_hh_to_VPVP(sm), 1e-16);
   // h -> gamma Z
   BOOST_CHECK_CLOSE_FRACTION(decays_SM_no_HO.partial_width_hh_to_VPVZ(&m), decays_sm_no_HO.partial_width_hh_to_VPVZ(sm), 1e-16);
}

BOOST_AUTO_TEST_CASE( test_SM_normalized_effective_couplings )
{

   Loop_library::set(-1);

   static constexpr double lambda = 0.194;

   SM_input_parameters input;
   input.LambdaIN = lambda;
   input.Qin = 1000;
   input.QEWSB = 173;

   Spectrum_generator_settings settings;
   softsusy::QedQcd qedqcd;

   SM_spectrum_generator<Two_scale> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);
   auto m = std::get<0>(spectrum_generator.get_models_slha());

   m.do_calculate_sm_pole_masses(true);
   m.solve_ewsb_tree_level();
   m.calculate_DRbar_masses();
   m.solve_ewsb();
   m.calculate_pole_masses();

   if (m.get_problems().have_problem()) {
      std::ostringstream ostr;
      m.get_problems().print_problems(ostr);
      BOOST_FAIL(ostr.str());
   }

   Physical_input physical_input;
   FlexibleDecay_settings flexibledecay_settings;
   flexibledecay_settings.set(FlexibleDecay_settings::calculate_normalized_effc, 1);

   // -----------------------------------------------------
   // decays with higher-order SM corrections

   SM_decays decays = SM_decays(m, qedqcd, physical_input, flexibledecay_settings);
   decays.calculate_decays();
   const auto effc = get_normalized_effective_couplings(decays.get_neutral_higgs_effc(), physical_input, qedqcd, settings, flexibledecay_settings);

   // tolerance in %
   BOOST_CHECK_CLOSE(effc[0].gg.second, 1, 0.005);
   BOOST_CHECK_CLOSE(effc[0].gamgam.second, 1, 0.01);
   BOOST_CHECK_CLOSE(effc[0].Zgam.second, 1, 0.02);

   BOOST_CHECK_CLOSE(effc[0].ZZ.second, 1, 0.04);
   BOOST_CHECK_CLOSE(effc[0].WW.second, 1, 0.007);

   BOOST_CHECK_CLOSE(std::real(effc[0].ee.second),     1, 0.002);
   BOOST_CHECK_CLOSE(std::real(effc[0].mumu.second),   1, 0.002);
   BOOST_CHECK_CLOSE(std::real(effc[0].tautau.second), 1, 0.002);

   BOOST_CHECK_CLOSE(std::real(effc[0].bb.second), 1, 0.009);
   BOOST_CHECK_CLOSE(std::real(effc[0].cc.second), 1, 0.0010);
   BOOST_CHECK_CLOSE(std::real(effc[0].ss.second), 1, 0.0020);
   BOOST_CHECK_CLOSE(std::real(effc[0].dd.second), 1, 0.002);
   BOOST_CHECK_CLOSE(std::real(effc[0].uu.second), 1, 0.002);
}
