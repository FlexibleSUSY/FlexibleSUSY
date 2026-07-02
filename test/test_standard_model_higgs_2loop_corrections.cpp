
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SM_higgs_loop_corrections

#include <boost/test/unit_test.hpp>

#include "standard_model.hpp"
#include "wrappers.hpp"
#include "sm_twoloophiggs.hpp"

using namespace flexiblesusy;
using namespace flexiblesusy::sm_twoloophiggs;

BOOST_AUTO_TEST_CASE( test_standard_model_2l_as_at )
{
   standard_model::Standard_model m;

   const double Q = 173.34;
   const double vev = 246.;
   const double g1 = 0.4;
   const double g2 = 0.5;
   const double g3 = 1.1;

   m.set_scale(Q);
   m.set_v(vev);
   m.set_mu2(vev*vev);
   m.set_g1(g1);
   m.set_g2(g2);
   m.set_g3(g3);
   m.set_Yu(0, 0, 0.001   * Sqrt(2.) / vev);
   m.set_Yu(1, 1, 0.010   * Sqrt(2.) / vev);
   m.set_Yu(2, 2, 165.0   * Sqrt(2.) / vev);
   m.set_Yd(0, 0, 0.001   * Sqrt(2.) / vev);
   m.set_Yd(1, 1, 0.010   * Sqrt(2.) / vev);
   m.set_Yd(2, 2, 2.9     * Sqrt(2.) / vev);
   m.set_Ye(0, 0, 0.001   * Sqrt(2.) / vev);
   m.set_Ye(1, 1, 0.010   * Sqrt(2.) / vev);
   m.set_Ye(2, 2, 1.77699 * Sqrt(2.) / vev);
   m.set_Lambdax(0.2);

   const int loops = 2;
   m.set_loops(loops);
   m.set_pole_mass_loop_order(loops);
   m.set_ewsb_loop_order(loops);
   m.solve_ewsb_tree_level();
   m.calculate_DRbar_masses();
   m.solve_ewsb();
   m.calculate_pole_masses();

   const double yt = m.get_Yu(2,2);
   const double mt = m.get_MFu(2);
   const double p = 125.;

   const double se_smh = -delta_mh_1loop_at_sm(
      p, Q, mt, yt);

   BOOST_CHECK_CLOSE_FRACTION(se_smh, 0., 1.0e-10);
}
