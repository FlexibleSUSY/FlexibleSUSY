
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SM_tree_level_spectrum

#include <boost/test/unit_test.hpp>

#include "test_SM.hpp"
#include "wrappers.hpp"
#include "SM_two_scale_model.hpp"
#include "standard_model.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_SM_tree_level_Higgs_masses )
{
   SM_input_parameters input;
   input.LambdaIN = 0.25;
   SM<Two_scale> m;
   setup_SM_const(m, input);

   // set mu2 which is not in agreement with EWSB
   const double mu2 = 100.;
   m.set_mu2(mu2);

   m.calculate_DRbar_masses();

   // check that mu2 was not changed
   BOOST_CHECK_CLOSE(m.get_mu2(), mu2, 1.0e-10);

   // neutral CP even Higgs
   const double hh(m.get_Mhh());
   const double lambda = m.get_Lambdax();
   const double v = m.get_v();
   const double hh_tree = Sqrt(lambda * Sqr(v));
   BOOST_CHECK_CLOSE(hh, hh_tree, 1.0e-12);

   // check Goldstone boson masses
   BOOST_CHECK_EQUAL(m.get_MHp(), m.get_MVWp());
   BOOST_CHECK_EQUAL(m.get_MAh(), m.get_MVZ());
}

#define CHECK_CLOSE_1(mass,eps)                                         \
   do {                                                                 \
      BOOST_CHECK_CLOSE(m.get_##mass(), sm.get_##mass(), (eps));        \
      BOOST_TEST_MESSAGE(#mass " = " << m.get_##mass() << " = " << sm.get_##mass()); \
   } while (false);

#define CHECK_CLOSE_N(mass,eps,N)                                       \
   do {                                                                 \
      for (int i = 0; i < N; i++) {                                     \
         BOOST_CHECK_CLOSE(m.get_##mass(i), sm.get_##mass(i), (eps));   \
         BOOST_TEST_MESSAGE(#mass "(" << i << ") = " << m.get_##mass(i) << " = " << sm.get_##mass(i)); \
      }                                                                 \
   } while (false);

BOOST_AUTO_TEST_CASE( test_SM_tree_level_masses )
{
   SM_input_parameters input;
   input.LambdaIN = 0.25;
   SM<Two_scale> m;
   standard_model::Standard_model sm;
   setup_SM_const(m, input);
   setup_SM_const(sm, input);

   m.calculate_DRbar_masses();
   sm.calculate_DRbar_masses();

   const double eps = 1e-15;

   CHECK_CLOSE_1(MVP, eps);
   CHECK_CLOSE_1(MVG, eps);
   CHECK_CLOSE_1(MVZ, eps);
   CHECK_CLOSE_1(MVWp, eps);
   CHECK_CLOSE_1(MAh, eps);
   CHECK_CLOSE_1(MHp, eps);
   CHECK_CLOSE_1(Mhh, eps);
   CHECK_CLOSE_N(MFv, eps, 3);
   CHECK_CLOSE_N(MFu, eps, 3);
   CHECK_CLOSE_N(MFd, eps, 3);
   CHECK_CLOSE_N(MFe, eps, 3);
}
