
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SM_FlexibleDecay

#include <iostream>

#include <boost/test/unit_test.hpp>

#include "test_SM.hpp"
#include "SM_two_scale_model.hpp"
#include "decays/SM_decays.hpp"
#include "decays/standard_model_decays.hpp"
#include "SM_mass_eigenstates_running.hpp"
#include "decays/experimental_constraints.hpp"
#include "SM_two_scale_spectrum_generator.hpp"

#include "lowe.h"
#include "loop_libraries/loop_library.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_SM_self_energies )
{

   Loop_library::set(-1);

   SM_input_parameters input;
   input.LambdaIN = 0.285;
   SM<Two_scale> m_;
   setup_SM_const(m_, input);

   m_.set_scale(100.0);

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
   SM_mass_eigenstates_running m(m_);

   const double p = 100.0;

   // test self energies and their derivatives agains FeynArts calculation

   // self energies

   // Higgs self energy
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_hh_1loop(p).real(), 1599.01623458672952438, 3e-16);

   // lepton self energy
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_Fe_1loop_PL(p, 1, 1).real(),  0.00034722428545620260187,  2e-15);
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_Fe_1loop_PR(p, 1, 1).real(),  0.000947893495480889715192,  5e-16);
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_Fe_1loop_1(p, 1, 1).real(),  -0.0002141475084938143376, 1e-15);

   // gauge field self energies
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_VZ_1loop(p).real(),  177.638580157661465364, 3e-15);
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_VP_1loop(p).real(),  62.905967975444248452, 6e-14);
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_VWp_1loop(p).real(), 38.2177737084397790568, 7e-15);

   // self energy derivatives

   // Higgs self energy derivative
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_hh_1loop_deriv_p2(p).real(), -0.0251636601702882253817, 3e-15);

   // lepton self energy derivative
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_Fe_1loop_PL_deriv_p2(p, 1, 1).real(), -2.06380698930894585634e-8, 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_Fe_1loop_PR_deriv_p2(p, 1, 1).real(),  2.03919861173351302228e-8, 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_Fe_1loop_1_deriv_p2(p, 1, 1).real(),   2.94359502557603049982e-8, 1e-15);

   // gauge field self energy derivatives
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_VZ_1loop_deriv_p2(p).real(),  -0.00615148335775731521652, 4e-15);
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_VP_1loop_deriv_p2(p).real(),  0.0000232991422811751857, 5e-13);
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_VWp_1loop_deriv_p2(p).real(), -0.00814803907885960389257, 2e-15);
}
