#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_NMSSM_self_energies

#include <boost/test/unit_test.hpp>

#include "LQS1_up_basis_mass_eigenstates.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_LQS1_self_energies )
{
   LQS1_up_basis_mass_eigenstates m {};

   m.do_force_output(true);

   m.set_g1(0.36);
   m.set_g2(0.65);
   m.set_g3(1.22);

   m.set_Lambdax(0.192);
   m.set_v(246);

   m.set_Ye(0, 0, 2.9e-6);
   m.set_Ye(1, 1, 3.6e-4);
   m.set_Ye(2, 2, 1.02e-2);

   m.set_Yu(0, 0, 7.9e-6);
   m.set_Yu(1, 1, 3.6e-3);
   m.set_Yu(2, 2, 9.9e-1);

   m.set_Yd(0, 0, 1.5e-5);
   m.set_Yd(1, 1, 3.4e-4);
   m.set_Yd(2, 2, 1.6e-2);

   for(int i=1; i<=2; i++){
       for(int j=0; j<i; j++){
            m.set_Ye(i, j, 0);
            m.set_Ye(j, i, 0);
            m.set_Yu(i, j, 0);
            m.set_Yu(j, i, 0);
            m.set_Yd(i, j, 0);
            m.set_Yd(j, i, 0);
            m.set_lamQL(i, j, 0);
            m.set_lamql(j, i, 0);
       }
   }

   m.set_MS12(2250000);
   m.set_gHS1(0.1);
   m.set_lamQL(1, 1, 0.01);
   m.set_lamql(1, 1, 0.01);

   m.set_lamQL(0, 0, 0);
   m.set_lamql(2, 2, 0);

   m.set_scale(100.0);

   m.solve_ewsb_tree_level();
   m.calculate_DRbar_masses();

   const double p = 100.0;

   // self energies

   // Higgs self energy
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_hh_1loop(p).real(), -16740.274395672364335, 3e-16);

   // lepton self energy
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_Fe_1loop_PL(p, 1, 1).real(),  0.00022814295048975526061,  2e-15);
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_Fe_1loop_PR(p, 1, 1).real(),  0.00090682235298948090600,  3e-16);
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_Fe_1loop_1(p, 1, 1).real(),  -0.000074186442764721124651, 4e-16);

   // gauge field self energies
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_VZ_1loop(p).real(),  186.26861711224057672, 3e-15);
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_VP_1loop(p).real(),  38.769814264189520259, 6e-14);
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_VWp_1loop(p).real(), 64.570698773836241458, 2e-15);

   // self energy derivatives

   // Higgs self energy derivative
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_hh_1loop_deriv_p2(p).real(), -0.030816481481987057639, 3e-15);

   // lepton self energy derivative
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_Fe_1loop_PL_deriv_p2(p, 1, 1).real(), -1.72285615708972490894e-8, 2e-15);
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_Fe_1loop_PR_deriv_p2(p, 1, 1).real(), 3.2990715812435773255e-8, 2e-15);
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_Fe_1loop_1_deriv_p2(p, 1, 1).real(), 9.81945896772624064422e-9, 2e-15);

   // gauge field self energy derivatives
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_VZ_1loop_deriv_p2(p).real(),  -0.0072128881287968869829, 4e-15);
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_VP_1loop_deriv_p2(p).real(),  -0.0003227074659516424894, 5e-13);
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_VWp_1loop_deriv_p2(p).real(), -0.0084718594643131729660, 2e-15);
}
