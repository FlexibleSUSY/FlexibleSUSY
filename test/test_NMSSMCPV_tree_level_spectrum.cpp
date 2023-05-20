
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_NMSSMCPV_tree_level_spectrum

#include <boost/test/unit_test.hpp>

#include "test.hpp"
#include "test_NMSSMCPV.hpp"
#include "NMSSM_two_scale_model.hpp"
#include "NMSSMCPV_two_scale_model.hpp"

#define COMPARE_MASS(p,dev) TEST_CLOSE(a.get_##p(), b.get_##p(), dev);

using namespace flexiblesusy;

void test_spectrum_equality(const NMSSMCPV<Two_scale>& a, const NMSSM<Two_scale>& b)
{
   COMPARE_MASS(MHpm, 1e-10);
   COMPARE_MASS(MChi, 1e-10);
   COMPARE_MASS(MCha, 1e-10);
   COMPARE_MASS(MGlu, 1e-10);
   COMPARE_MASS(MVP , 1e-10);
   COMPARE_MASS(MVZ , 1e-10);
   COMPARE_MASS(MVWm, 1e-10);
   COMPARE_MASS(MVG , 1e-10);
   COMPARE_MASS(MSu , 1e-10);
   COMPARE_MASS(MSd , 1e-10);
   COMPARE_MASS(MSe , 1e-10);
   COMPARE_MASS(MSv , 1e-10);
   COMPARE_MASS(MFu , 1e-10);
   COMPARE_MASS(MFd , 1e-10);
   COMPARE_MASS(MFe , 1e-10);
   COMPARE_MASS(MFv , 1e-10);

   // Higgses
   const Eigen::Array<double,3,1> Mhh(b.get_Mhh());
   const Eigen::Array<double,3,1> MAh(b.get_MAh());
   const Eigen::Array<double,6,1> Mhh_full(a.get_Mhh());
   Eigen::Array<double,6,1> Mhh_mixed;
   Mhh_mixed << Mhh(0), Mhh(1), Mhh(2), MAh(0), MAh(1), MAh(2);

   std::sort(Mhh_mixed.data(), Mhh_mixed.data() + Mhh_mixed.size());

   TEST_CLOSE(Mhh_full, Mhh_mixed, 1e-10);
}

void copy_parameters(const NMSSMCPV<Two_scale>& a, NMSSM<Two_scale>& b)
{
   b.set_scale (a.get_scale ());
   b.set_loops (a.get_loops ());
   b.set_g1    (a.get_g1    ());
   b.set_g2    (a.get_g2    ());
   b.set_g3    (a.get_g3    ());
   b.set_Yu    (a.get_Yu    ().real());
   b.set_Yd    (a.get_Yd    ().real());
   b.set_Ye    (a.get_Ye    ().real());
   b.set_MassB (a.get_MassB ().real());
   b.set_MassG (a.get_MassG ().real());
   b.set_MassWB(a.get_MassWB().real());
   b.set_mq2   (a.get_mq2   ().real());
   b.set_ml2   (a.get_ml2   ().real());
   b.set_md2   (a.get_md2   ().real());
   b.set_mu2   (a.get_mu2   ().real());
   b.set_me2   (a.get_me2   ().real());
   b.set_mHd2  (a.get_mHd2  ());
   b.set_mHu2  (a.get_mHu2  ());
   b.set_ms2   (a.get_ms2   ());
   b.set_TYu   (a.get_TYu   ().real());
   b.set_TYd   (a.get_TYd   ().real());
   b.set_TYe   (a.get_TYe   ().real());
   b.set_vu    (a.get_vu    ());
   b.set_vd    (a.get_vd    ());
   b.set_vS    (a.get_vS    ());
   b.set_TLambdax(a.get_TLambdax().real());
   b.set_TKappa  (a.get_TKappa  ().real());
   b.set_Lambdax (a.get_Lambdax ().real());
   b.set_Kappa   (a.get_Kappa   ().real());
}

BOOST_AUTO_TEST_CASE( test_NMSSMCPV_tree_level_spectrum )
{
   NMSSMCPV_input_parameters input;
   input.m0 = 250.;
   input.m12 = 200.;
   input.TanBeta = 10.;
   input.Azero = -500.;
   input.LambdaInput = 0.1;
   input.KappaInput = 0.1;
   input.LambdaTimesvSInput = 100;
   // eta and etaS must not be both equal to zero, because this will
   // lead to a division by zero in the tree-level EWSB eqs.
   input.etaInput = 1e-8;
   input.etaSInput = 1e-10;

   NMSSMCPV<Two_scale> m1(input);
   NMSSM<Two_scale> m2;

   setup_NMSSMCPV(m1, input);

   m1.solve_ewsb_tree_level();
   copy_parameters(m1, m2);
   m1.calculate_DRbar_masses();
   m2.calculate_DRbar_masses();

   test_spectrum_equality(m1, m2);

   BOOST_REQUIRE(gErrors == 0);
   if (gErrors) {
      BOOST_FAIL("spectra are not equal");
      gErrors = 0;
   }
}

void check_goldstone_masses(const NMSSMCPV_input_parameters& input)
{
   NMSSMCPV<Two_scale> m;
   setup_NMSSMCPV(m, input);

   m.solve_ewsb_tree_level();
   m.calculate_DRbar_masses();
   m.reorder_DRbar_masses();

   BOOST_CHECK_CLOSE(m.get_Mhh(0) , m.get_MVZ() , 1e-10);
   BOOST_CHECK_CLOSE(m.get_MHpm(0), m.get_MVWm(), 1e-10);
}

BOOST_AUTO_TEST_CASE( test_NMSSMCPV_goldstone_boson_masses )
{
   for (int m = 0; m < 10; m++) {
      for (int e = 0; e < 10; e++) {
         NMSSMCPV_input_parameters input;
         input.m0 = 250.;
         input.m12 = 200.;
         input.TanBeta = 10.;
         input.Azero = -500.;
         input.LambdaInput = 0.1;
         input.KappaInput = 0.1;
         input.LambdaTimesvSInput = 100;
         // eta and etaS must not be both equal, because this will
         // lead to a division by zero in the tree-level EWSB eqs.
         input.etaInput = 2*Pi*e/10 + 1e-8;
         input.etaSInput = 2*Pi*m/10 + 1e-10;

         check_goldstone_masses(input);
      }
   }
}
