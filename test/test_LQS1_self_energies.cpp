#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_NMSSM_self_energies

#include <boost/test/unit_test.hpp>

#include "LQS1_up_basis_mass_eigenstates.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_LQS1_self_energies )
{
   LQS1_up_basis_mass_eigenstates m {};

   m.set_g1(0.3);
   m.set_g2(0.3);
   m.set_g3(0.3);

   m.solve_ewsb_tree_level();
   m.calculate_DRbar_masses();

   const double p = 100.; // Jonas, set the momentum to a correct value

   // test self energies and their derivatives agains FeynArts calculation

   // self energies
   BOOST_CHECK_CLOSE_FRACTION(m.self_energy_Hp_1loop(p).real(), 1, 1e-16);

   // self energy derivatives
}
