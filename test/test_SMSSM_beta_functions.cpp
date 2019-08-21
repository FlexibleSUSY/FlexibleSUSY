
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SMSSM_beta_functions

#include <boost/test/unit_test.hpp>

#include "test_legacy.hpp"
#include "test_SMSSM.hpp"
#include "wrappers.hpp"
#include "ew_input.hpp"
#include "nmssmsoftsusy.h"
#include "SMSSM_two_scale_model.hpp"

using namespace flexiblesusy;
using namespace softsusy;

void test_parameter_equality(const SoftParsNmssm& a, const SMSSM_soft_parameters& b)
{
   TEST_EQUALITY(a.displayLoops(), b.get_loops());
   TEST_EQUALITY(a.displayMu(), b.get_scale());
   TEST_EQUALITY(a.displayThresholds(), b.get_thresholds());

   TEST_EQUALITY(a.displayGaugeCoupling(1), b.get_g1());
   TEST_EQUALITY(a.displayGaugeCoupling(2), b.get_g2());
   TEST_EQUALITY(a.displayGaugeCoupling(3), b.get_g3());

   TEST_EQUALITY(a.displayYukawaMatrix(YU), b.get_Yu());
   TEST_EQUALITY(a.displayYukawaMatrix(YD), b.get_Yd());
   TEST_EQUALITY(a.displayYukawaMatrix(YE), b.get_Ye());

   TEST_EQUALITY(a.displayGaugino(1), b.get_MassB());
   TEST_EQUALITY(a.displayGaugino(2), b.get_MassWB());
   TEST_EQUALITY(a.displayGaugino(3), b.get_MassG());

   TEST_EQUALITY(a.displayMh1Squared(), b.get_mHd2());
   TEST_EQUALITY(a.displayMh2Squared(), b.get_mHu2());
   TEST_EQUALITY(a.displaySoftMassSquared(mQl), b.get_mq2());
   TEST_EQUALITY(a.displaySoftMassSquared(mUr), b.get_mu2());
   TEST_EQUALITY(a.displaySoftMassSquared(mDr), b.get_md2());
   TEST_EQUALITY(a.displaySoftMassSquared(mLl), b.get_ml2());
   TEST_EQUALITY(a.displaySoftMassSquared(mEr), b.get_me2());

   TEST_EQUALITY(a.displayTrilinear(UA), b.get_TYu());
   TEST_EQUALITY(a.displayTrilinear(DA), b.get_TYd());
   TEST_EQUALITY(a.displayTrilinear(EA), b.get_TYe());

   TEST_EQUALITY(a.displaySusyMu(), b.get_Mu());
   TEST_EQUALITY(a.displayM3Squared(), b.get_BMu());

   TEST_EQUALITY(a.displayLambda(), b.get_Lambdax());
   TEST_EQUALITY(a.displayKappa(), b.get_Kappa());
   TEST_EQUALITY(a.displaySvev(), b.get_vS());

   TEST_EQUALITY(a.displayTrialambda(), b.get_TLambdax());
   TEST_EQUALITY(a.displayTriakappa(), b.get_TKappa());
   TEST_EQUALITY(a.displayMsSquared(), b.get_ms2());

   // Z3 violating terms
   TEST_EQUALITY(a.displayMspSquared(), b.get_BMS());
   TEST_EQUALITY(a.displayXiF()       , b.get_L1());
   TEST_EQUALITY(a.displayXiS()       , b.get_LL1());
   TEST_EQUALITY(a.displayMupr()      , b.get_MS());

   const double vu = b.get_vu(), vd = b.get_vd();
   double tanBeta;
   if (is_zero(vu))
      tanBeta = 0.;
   else if (is_zero(vd))
      tanBeta = std::numeric_limits<double>::infinity();
   else
      tanBeta = vu/vd;
   const double vev = sqrt(sqr(b.get_vu()) + sqr(b.get_vd()));
   TEST_EQUALITY(a.displayTanb(), tanBeta);
   TEST_EQUALITY(a.displayHvev(), vev);
}

void test_anomalous_dimensions_equality(const SoftParsNmssm& a, const SMSSM_soft_parameters& b)
{
  DoubleMatrix gEE(3,3),gLL(3,3),gQQ(3,3),gDD(3,3),gUU(3,3);
  double gH1H1 = 0.0, gH2H2 = 0.0, gSS = 0.0;
  DoubleVector dg(1,3);
  nmsBrevity brevity;
  a.anomalousDimension(gEE, gLL, gQQ, gUU, gDD, dg, gH1H1, gH2H2, gSS, brevity);

  TEST_EQUALITY(a.displayLoops(), b.get_loops());
  TEST_EQUALITY(a.displayMu(), b.get_scale());
  TEST_EQUALITY(a.displayThresholds(), b.get_thresholds());

  TEST_EQUALITY(gEE, b.get_SeRSeR());
  TEST_EQUALITY(gLL, b.get_SlSl());
  TEST_EQUALITY(gQQ, b.get_SqSq());
  TEST_EQUALITY(gUU, b.get_SuRSuR());
  TEST_EQUALITY(gDD, b.get_SdRSdR());
  TEST_EQUALITY(gH1H1, b.get_SHdSHd());
  TEST_EQUALITY(gH2H2, b.get_SHuSHu());
  TEST_EQUALITY(gSS, b.get_SsRSsR());
}

void test_beta_function_equality(const SoftParsNmssm& a, const SMSSM_soft_parameters& b)
{
   SoftParsNmssm beta_a(a.beta2());
   SMSSM_soft_parameters beta_b(b.calc_beta());

   TEST_EQUALITY(beta_a.displayLoops(), beta_b.get_loops());
   TEST_EQUALITY(beta_a.displayMu(), beta_b.get_scale());
   TEST_EQUALITY(beta_a.displayThresholds(), beta_b.get_thresholds());

   TEST_EQUALITY(beta_a.displayGaugeCoupling(1), beta_b.get_g1());
   TEST_EQUALITY(beta_a.displayGaugeCoupling(2), beta_b.get_g2());
   TEST_EQUALITY(beta_a.displayGaugeCoupling(3), beta_b.get_g3());

   TEST_EQUALITY(beta_a.displayYukawaMatrix(YU), beta_b.get_Yu());
   TEST_EQUALITY(beta_a.displayYukawaMatrix(YD), beta_b.get_Yd());
   TEST_EQUALITY(beta_a.displayYukawaMatrix(YE), beta_b.get_Ye());

   TEST_EQUALITY(beta_a.displayLambda(), beta_b.get_Lambdax());
   TEST_EQUALITY(beta_a.displayKappa(), beta_b.get_Kappa());

   TEST_EQUALITY(beta_a.displayGaugino(1), beta_b.get_MassB());
   TEST_EQUALITY(beta_a.displayGaugino(2), beta_b.get_MassWB());
   TEST_EQUALITY(beta_a.displayGaugino(3), beta_b.get_MassG());

   TEST_EQUALITY(beta_a.displayMh1Squared(), beta_b.get_mHd2());
   TEST_CLOSE(beta_a.displayMh2Squared(), beta_b.get_mHu2(), 2.0e-12);
   TEST_EQUALITY(beta_a.displayMsSquared(), beta_b.get_ms2());
   TEST_CLOSE(beta_a.displaySoftMassSquared(mQl), beta_b.get_mq2(), 2.0e-12);
   TEST_CLOSE(beta_a.displaySoftMassSquared(mUr), beta_b.get_mu2(), 4.0e-12);
   TEST_CLOSE(beta_a.displaySoftMassSquared(mDr), beta_b.get_md2(), 2.0e-12);
   TEST_EQUALITY(beta_a.displaySoftMassSquared(mLl), beta_b.get_ml2());
   TEST_EQUALITY(beta_a.displaySoftMassSquared(mEr), beta_b.get_me2());

   TEST_EQUALITY(beta_a.displayTrilinear(UA), beta_b.get_TYu());
   TEST_EQUALITY(beta_a.displayTrilinear(DA), beta_b.get_TYd());
   TEST_EQUALITY(beta_a.displayTrilinear(EA), beta_b.get_TYe());
   TEST_EQUALITY(beta_a.displayTrialambda() , beta_b.get_TLambdax());
   TEST_EQUALITY(beta_a.displayTriakappa()  , beta_b.get_TKappa());

   TEST_EQUALITY(beta_a.displayMspSquared(), beta_b.get_BMS());
   TEST_EQUALITY(beta_a.displayXiF()       , beta_b.get_L1());
   TEST_EQUALITY(beta_a.displayXiS()       , beta_b.get_LL1());
   TEST_EQUALITY(beta_a.displayMupr()      , beta_b.get_MS());

   TEST_EQUALITY(beta_a.displaySusyMu()   , beta_b.get_Mu());
   TEST_EQUALITY(beta_a.displayM3Squared(), beta_b.get_BMu());

   const double vu = b.get_vu();
   const double vd = b.get_vd();
   const double tanBeta = vu / vd;
   const double beta_tanBeta = tanBeta * (beta_b.get_vu()/vu - beta_b.get_vd() / vd);
   const double vev = sqrt(sqr(vu) + sqr(vd));
   const double beta_vev = (vu * beta_b.get_vu() + vd * beta_b.get_vd()) / vev;

   TEST_EQUALITY(a.displayTanb(), tanBeta);
   TEST_EQUALITY(beta_a.displayTanb(), beta_tanBeta);
   TEST_EQUALITY(a.displayHvev(), vev);
   TEST_EQUALITY(beta_a.displayHvev(), beta_vev);
   TEST_EQUALITY(beta_a.displaySvev(), beta_b.get_vS());
}

BOOST_AUTO_TEST_CASE( test_SMSSM_beta_functions )
{
   SMSSM_input_parameters input;
   input.m0 = 200.;
   input.m12 = 300.;
   input.TanBeta = 10;
   input.SignMu = 1;
   input.Azero = 200;
   input.LambdaInput = 0.1;
   input.KappaInput = 0.2;
   input.LambdaSInput = 100;
   input.L1Input = 400.;
   input.MSInput = 500.;
   input.BMSInput = 600.;
   SMSSM<Two_scale> m;
   NmssmSoftsusy s;
   setup_SMSSM(m, s, input);

   m.set_loops(2);
   s.setLoops(2);

   test_parameter_equality(s, m);
   BOOST_REQUIRE(gErrors == 0);
   if (gErrors) {
      BOOST_FAIL("parameters are not equal");
      gErrors = 0;
   }

   test_anomalous_dimensions_equality(s, m);
   BOOST_CHECK(gErrors == 0);
   if (gErrors) {
      BOOST_FAIL("anomalous dimensions are not equal");
      gErrors = 0;
   }

   test_beta_function_equality(s, m);
   BOOST_CHECK(gErrors == 0);
   if (gErrors) {
      BOOST_FAIL("beta functions are not equal");
      gErrors = 0;
   }
}
