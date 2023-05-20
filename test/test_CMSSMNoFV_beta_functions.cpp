
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSMNoFV_beta_functions

#include <boost/test/unit_test.hpp>

#include "test.hpp"
#include "test_CMSSMNoFV.hpp"
#include "CMSSM_two_scale_model.hpp"
#include "CMSSMNoFV_two_scale_model.hpp"

#define COMPARE_PARAMETERS(p) TEST_EQUALITY(a.get_##p(), b.get_##p());

using namespace flexiblesusy;

template <class T1, class T2>
void test_parameter_equality(const T1& a, const T2& b)
{
   COMPARE_PARAMETERS(loops);
   COMPARE_PARAMETERS(scale);
   COMPARE_PARAMETERS(thresholds);

   COMPARE_PARAMETERS(g1);
   COMPARE_PARAMETERS(g2);
   COMPARE_PARAMETERS(g3);

   COMPARE_PARAMETERS(Yu);
   COMPARE_PARAMETERS(Yd);
   COMPARE_PARAMETERS(Ye);

   COMPARE_PARAMETERS(MassB);
   COMPARE_PARAMETERS(MassWB);
   COMPARE_PARAMETERS(MassG);

   COMPARE_PARAMETERS(mHd2);
   COMPARE_PARAMETERS(mHu2);
   COMPARE_PARAMETERS(mq2);
   COMPARE_PARAMETERS(mu2);
   COMPARE_PARAMETERS(md2);
   COMPARE_PARAMETERS(ml2);
   COMPARE_PARAMETERS(me2);

   COMPARE_PARAMETERS(TYu);
   COMPARE_PARAMETERS(TYd);
   COMPARE_PARAMETERS(TYe);

   COMPARE_PARAMETERS(Mu);
   COMPARE_PARAMETERS(BMu);

   COMPARE_PARAMETERS(vu);
   COMPARE_PARAMETERS(vd);
}

template <class T1, class T2>
void test_beta_function_equality(const T1& a, const T2& b)
{
   T1 beta_a(a.calc_beta());
   T2 beta_b(b.calc_beta());

   test_parameter_equality(beta_a, beta_b);
}

BOOST_AUTO_TEST_CASE( test_CMSSMNoFV_beta_functions )
{
   CMSSMNoFV_input_parameters input;
   input.TanBeta = 10.;
   input.m0 = 125.;
   input.m12 = 200.;
   input.SignMu = 1;
   input.Azero = 0.;

   CMSSMNoFV<Two_scale> m1;
   CMSSM<Two_scale> m2;
   setup_CMSSM_models(m1, m2, input);

   test_parameter_equality(static_cast<CMSSMNoFV_soft_parameters>(m1),
                           static_cast<CMSSM_soft_parameters>(m2));
   BOOST_REQUIRE(gErrors == 0);
   if (gErrors) {
      BOOST_FAIL("parameters are not equal");
      gErrors = 0;
   }

   test_beta_function_equality(static_cast<CMSSMNoFV_soft_parameters>(m1),
                               static_cast<CMSSM_soft_parameters>(m2));
   BOOST_CHECK(gErrors == 0);
   if (gErrors) {
      BOOST_FAIL("beta functions are not equal");
      gErrors = 0;
   }
}
