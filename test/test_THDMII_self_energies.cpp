#include "iomanip"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_NMSSM_self_energies

#include <boost/test/unit_test.hpp>

#include "THDMII_mass_eigenstates.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_THDMII_self_energies )
{
   THDMII_mass_eigenstates m {};

   m.do_force_output(true);

   m.set_g1(0.36);
   m.set_g2(0.65);
   m.set_g3(1.22);

   const double tanb = 10;
   const double vev = 246;

   m.set_v1(vev/std::sqrt(1 + tanb*tanb));
   m.set_v2(vev*tanb/std::sqrt(1 + tanb*tanb));

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
       }
   }

   m.set_M122(1e+4);

   m.set_Lambda1(0);
   m.set_Lambda2(0);
   m.set_Lambda3(0);
   m.set_Lambda4(0);
   m.set_Lambda5(0);
   m.set_Lambda6(0);
   m.set_Lambda7(0);

   m.set_scale(100.0);

   m.solve_ewsb_tree_level();
   m.calculate_DRbar_masses();

   const double p = 100.0;

   std::cout << std::setprecision (20) << '\n';
   std::cout << " Scale: " << m.get_scale() << '\n';

   std::cout << "\n Higgs Mass: " << m.get_Mhh() << '\n';
   std::cout << " Ah Mass: " << m.get_MAh() << '\n';
   std::cout << " Hpm Mass: " << m.get_MHm() << '\n';

   std::cout << "\n Z Mass: " << m.get_MVZ() << '\n';
   std::cout << " W Mass: " << m.get_MVWm() << '\n';

   std::cout << "\n vev 1: " << m.get_v1() << '\n';
   std::cout << "\n vev 2: " << m.get_v2() << '\n';
   std::cout << "\n Electron Mass: " << m.get_MFe(0) << '\n';
   std::cout << " Muon Mass: " << m.get_MFe(1) << '\n';
   std::cout << " Tau Mass: " << m.get_MFe(2) << '\n';

   std::cout << "\n up Mass: " << m.get_MFu(0) << '\n';
   std::cout << " charm Mass: " << m.get_MFu(1) << '\n';
   std::cout << " top Mass: " << m.get_MFu(2) << '\n';

   std::cout << "\n down Mass: " << m.get_MFd(0) << '\n';
   std::cout << " strange Mass: " << m.get_MFd(1) << '\n';
   std::cout << " botton Mass: " << m.get_MFd(2) << '\n';


   std::cout << "\n thetaW: " << m.ThetaW() << '\n';

   std::cout << " g1: " << m.get_g1() << '\n';
   std::cout << " g2: " << m.get_g2() << '\n';

   // self energies

   // Higgs Tadpole
   BOOST_CHECK_CLOSE_FRACTION(m.tadpole_hh_1loop(0).real(), -4611359.147766369394958019, 3e-15);

   // Higgs self energy
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_hh_1loop(p, 0, 0).real(), -16740.274395672364335, 3e-16);

   // lepton self energy
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_Fe_1loop_PL(p, 1, 1).real(),  0.00022814295048975526061,  2e-15);
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_Fe_1loop_PR(p, 1, 1).real(),  0.00090682235298948090600,  3e-16);
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_Fe_1loop_1(p, 1, 1).real(),  -0.000074186442764721124651, 4e-16);

   // gauge field self energies
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_VZ_1loop(p).real(),  186.26861711224057672, 3e-15);
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_VP_1loop(p).real(),  38.769814264189520259, 6e-14);
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_VWm_1loop(p).real(), 64.570698773836241458, 2e-15);

   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_VPVZ_1loop(p).real(),  -8.4300615408824981500, 3e-15);

   // self energy derivatives

   // Higgs self energy derivative
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_hh_1loop_deriv_p2(p, 0, 0).real(), -0.030816481481987057639, 3e-15);

   // lepton self energy derivative
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_Fe_1loop_PL_deriv_p2(p, 1, 1).real(), -1.72285615708972490894e-8, 2e-15);
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_Fe_1loop_PR_deriv_p2(p, 1, 1).real(), 3.2990715812435773255e-8, 2e-15);
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_Fe_1loop_1_deriv_p2(p, 1, 1).real(), 9.81945896772624064422e-9, 2e-15);

   // gauge field self energy derivatives
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_VZ_1loop_deriv_p2(p).real(),  -0.0072128881287968869829, 4e-15);
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_VP_1loop_deriv_p2(p).real(),  -0.0003227074659516424894, 5e-13);
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_VWm_1loop_deriv_p2(p).real(), -0.0084718594643131729660, 2e-15);
}
