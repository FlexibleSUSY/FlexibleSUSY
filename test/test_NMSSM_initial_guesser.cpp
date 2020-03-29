
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_NMSSM_initial_guesser

#include <boost/test/unit_test.hpp>
#include "test_legacy.hpp"
#include "test_NMSSM.hpp"

#include "NMSSM_two_scale_model.hpp"
#include "NMSSM_two_scale_initial_guesser.hpp"
#include "SoftsusyNMSSM_parameter_point.hpp"
#include "SoftsusyNMSSM_two_scale.hpp"
#include "SoftsusyNMSSM_two_scale_initial_guesser.hpp"
#include "SoftsusyNMSSM_two_scale_sugra_constraint.hpp"
#include "SoftsusyNMSSM_two_scale_susy_scale_constraint.hpp"
#include "SoftsusyNMSSM_two_scale_low_scale_constraint.hpp"

softsusy::QedQcd convert(const softsusy::QedQcd_legacy& ql)
{
   softsusy::QedQcd qn;

   qn.setAlphas(flexiblesusy::ToEigenArray(ql.displayAlphas()));
   qn.setMasses(flexiblesusy::ToEigenArray(ql.displayMass()));
   qn.set_input(ql.display_input());
   qn.setPoleMb(ql.displayPoleMb());
   qn.setCKM(ql.displayCKM());
   qn.setPMNS(ql.displayPMNS());
   qn.set_number_of_parameters(ql.howMany());
   qn.set_scale(ql.displayMu());
   qn.set_loops(ql.displayLoops());
   qn.set_thresholds(ql.displayThresholds());

   return qn;
}

BOOST_AUTO_TEST_CASE( test_initial_guess )
{
   softsusy::TOLERANCE = 1.0e-3;
   softsusy::GUTlambda = true;
   NMSSM<Two_scale> m;
   m.set_loops(2);
   m.set_thresholds(2);
   SoftsusyNMSSM<Two_scale> snmssm;

   // create NMSSM initial guesser
   NMSSM_input_parameters input;
   input.m0 = 250.; // avoids tree-level tachyons
   input.m12 = 200.;
   input.TanBeta = 10.;
   input.Azero = -500.;
   input.LambdaInput = 0.1;
   input.SignvS = 1;
   QedQcd_legacy qedqcd;

   m.set_input_parameters(input);
   NMSSM_low_scale_constraint<Two_scale>  low_constraint(&m, convert(qedqcd));
   NMSSM_susy_scale_constraint<Two_scale> susy_constraint(&m, convert(qedqcd));
   NMSSM_high_scale_constraint<Two_scale> high_constraint(&m);

   NMSSM_initial_guesser<Two_scale> guesser(&m, convert(qedqcd), low_constraint,
                                            susy_constraint, high_constraint);

   // create SoftsusyNMSSM initial guesser
   SoftsusyNMSSM_parameter_point pp;
   pp.m0 = input.m0;
   pp.m12 = input.m12;
   pp.a0 = input.Azero;
   pp.mxGuess = high_constraint.get_scale();
   pp.tanBeta = input.TanBeta;
   pp.lambda = input.LambdaInput;
   pp.kappa = 0.1;  // initial guess at the low-scale
   pp.svev = 1000.; // initial guess at the low-scale
   SoftsusyNMSSM_sugra_constraint mssm_sugra_constraint(pp);
   SoftsusyNMSSM_low_scale_constraint mssm_mz_constraint(pp);
   SoftsusyNMSSM_susy_scale_constraint mssm_msusy_constraint(pp);
   SoftsusyNMSSM_initial_guesser initial_guesser(&snmssm, pp, mssm_mz_constraint,
                                          mssm_msusy_constraint,
                                          mssm_sugra_constraint);
   initial_guesser.set_QedQcd(qedqcd);

   // guess both models
   guesser.guess();
   initial_guesser.guess();

   BOOST_CHECK_EQUAL(snmssm.displayLoops()     , m.get_loops());
   BOOST_CHECK_EQUAL(snmssm.displayMu()        , m.get_scale());
   // BOOST_CHECK_EQUAL(snmssm.displayThresholds(), m.get_thresholds());

   BOOST_CHECK_CLOSE_FRACTION(snmssm.displayGaugeCoupling(1), m.get_g1(), 0.003);
   BOOST_CHECK_CLOSE_FRACTION(snmssm.displayGaugeCoupling(2), m.get_g2(), 0.01);
   BOOST_CHECK_CLOSE_FRACTION(snmssm.displayGaugeCoupling(3), m.get_g3(), 1.0e-5);

   // test off-diagonal elements
   BOOST_TEST_MESSAGE("testing off-diagonal yukawa elements");
   for (int i = 1; i <= 3; i++) {
      for (int k = 1; k <= 3; k++) {
         if (i == k)
            continue;
         BOOST_TEST_MESSAGE("testing yukawa elements " << i << ", " << k);
         BOOST_CHECK_CLOSE_FRACTION(m.get_Yu()(i-1,k-1), snmssm.displayYukawaMatrix(YU)(i,k), 1.0e-10);
         BOOST_CHECK_CLOSE_FRACTION(m.get_Yd()(i-1,k-1), snmssm.displayYukawaMatrix(YD)(i,k), 1.0e-10);
         BOOST_CHECK_CLOSE_FRACTION(m.get_Ye()(i-1,k-1), snmssm.displayYukawaMatrix(YE)(i,k), 1.0e-10);
      }
   }

   BOOST_TEST_MESSAGE("testing diagonal yukawa elements");
   BOOST_CHECK_CLOSE_FRACTION(m.get_Yu()(0,0), snmssm.displayYukawaMatrix(YU)(1,1), 1.0e-4);
   BOOST_CHECK_CLOSE_FRACTION(m.get_Yd()(0,0), snmssm.displayYukawaMatrix(YD)(1,1), 1.0e-4);
   BOOST_CHECK_CLOSE_FRACTION(m.get_Ye()(0,0), snmssm.displayYukawaMatrix(YE)(1,1), 1.0e-4);

   BOOST_CHECK_CLOSE_FRACTION(m.get_Yu()(1,1), snmssm.displayYukawaMatrix(YU)(2,2), 1.0e-4);
   BOOST_CHECK_CLOSE_FRACTION(m.get_Yd()(1,1), snmssm.displayYukawaMatrix(YD)(2,2), 1.0e-4);
   BOOST_CHECK_CLOSE_FRACTION(m.get_Ye()(1,1), snmssm.displayYukawaMatrix(YE)(2,2), 1.0e-4);

   BOOST_CHECK_CLOSE_FRACTION(m.get_Yu()(2,2), snmssm.displayYukawaMatrix(YU)(3,3), 1.0e-4);
   BOOST_CHECK_CLOSE_FRACTION(m.get_Yd()(2,2), snmssm.displayYukawaMatrix(YD)(3,3), 1.0e-4);
   BOOST_CHECK_CLOSE_FRACTION(m.get_Ye()(2,2), snmssm.displayYukawaMatrix(YE)(3,3), 1.0e-4);

   const double vu = m.get_vu(), vd = m.get_vd(), vS = m.get_vS();
   double tanBeta;
   if (is_zero(vu))
      tanBeta = 0.;
   else if (is_zero(vd))
      tanBeta = std::numeric_limits<double>::infinity();
   else
      tanBeta = vu/vd;
   const double vev = Sqrt(Sqr(m.get_vu()) + Sqr(m.get_vd()));
   BOOST_CHECK_CLOSE_FRACTION(snmssm.displayTanb(), tanBeta, 3.0e-5);
   BOOST_CHECK_CLOSE_FRACTION(snmssm.displayHvev(), vev    , 2.0e-5);
   BOOST_CHECK_CLOSE_FRACTION(snmssm.displaySvev(), vS     , 0.01);

   BOOST_CHECK_CLOSE_FRACTION(snmssm.displayGaugino(1), m.get_MassB() , 0.01);
   BOOST_CHECK_CLOSE_FRACTION(snmssm.displayGaugino(2), m.get_MassWB(), 0.01);
   BOOST_CHECK_CLOSE_FRACTION(snmssm.displayGaugino(3), m.get_MassG() , 0.001);

   BOOST_CHECK_CLOSE_FRACTION(snmssm.displayMh1Squared(), m.get_mHd2(), 0.05);
   BOOST_CHECK_CLOSE_FRACTION(snmssm.displayMh2Squared(), m.get_mHu2(), 5.0e-4);

   // BOOST_CHECK_CLOSE_FRACTION(snmssm.displaySoftMassSquared(mQl), m.get_mq2(), 1.0e-5);
   // BOOST_CHECK_CLOSE_FRACTION(snmssm.displaySoftMassSquared(mUr), m.get_mu2(), 1.0e-5);
   // BOOST_CHECK_CLOSE_FRACTION(snmssm.displaySoftMassSquared(mDr), m.get_md2(), 1.0e-5);
   // BOOST_CHECK_CLOSE_FRACTION(snmssm.displaySoftMassSquared(mLl), m.get_ml2(), 1.0e-5);
   // BOOST_CHECK_CLOSE_FRACTION(snmssm.displaySoftMassSquared(mEr), m.get_me2(), 1.0e-5);

   const DoubleMatrix Au(snmssm.displayTrilinear(UA)),
      Ad(snmssm.displayTrilinear(DA)),
      Ae(snmssm.displayTrilinear(EA));
   const Eigen::Matrix<double,3,3> TYu(m.get_TYu()),
      TYd(m.get_TYd()), TYe(m.get_TYe());

   BOOST_CHECK_CLOSE_FRACTION(Au(1,1), TYu(0,0), 1.5e-3);
   BOOST_CHECK_CLOSE_FRACTION(Au(1,2), TYu(0,1), 1.0e-5);
   BOOST_CHECK_CLOSE_FRACTION(Au(1,3), TYu(0,2), 1.0e-5);

   BOOST_CHECK_CLOSE_FRACTION(Au(2,1), TYu(1,0), 1.0e-5);
   BOOST_CHECK_CLOSE_FRACTION(Au(2,2), TYu(1,1), 1.5e-3);
   BOOST_CHECK_CLOSE_FRACTION(Au(2,3), TYu(1,2), 1.0e-5);

   BOOST_CHECK_CLOSE_FRACTION(Au(3,1), TYu(2,0), 1.0e-5);
   BOOST_CHECK_CLOSE_FRACTION(Au(3,2), TYu(2,1), 1.0e-5);
   BOOST_CHECK_CLOSE_FRACTION(Au(3,3), TYu(2,2), 3.0e-3);

   BOOST_CHECK_CLOSE_FRACTION(Ad(1,1), TYd(0,0), 2.5e-3);
   BOOST_CHECK_CLOSE_FRACTION(Ad(1,2), TYd(0,1), 1.0e-5);
   BOOST_CHECK_CLOSE_FRACTION(Ad(1,3), TYd(0,2), 1.0e-5);

   BOOST_CHECK_CLOSE_FRACTION(Ad(2,1), TYd(1,0), 1.0e-5);
   BOOST_CHECK_CLOSE_FRACTION(Ad(2,2), TYd(1,1), 2.5e-3);
   BOOST_CHECK_CLOSE_FRACTION(Ad(2,3), TYd(1,2), 1.0e-5);

   BOOST_CHECK_CLOSE_FRACTION(Ad(3,1), TYd(2,0), 1.0e-5);
   BOOST_CHECK_CLOSE_FRACTION(Ad(3,2), TYd(2,1), 1.0e-5);
   BOOST_CHECK_CLOSE_FRACTION(Ad(3,3), TYd(2,2), 2.0e-3);

   BOOST_CHECK_CLOSE_FRACTION(Ae(1,1), TYe(0,0), 2.5e-3);
   BOOST_CHECK_CLOSE_FRACTION(Ae(1,2), TYe(0,1), 1.0e-5);
   BOOST_CHECK_CLOSE_FRACTION(Ae(1,3), TYe(0,2), 1.0e-5);

   BOOST_CHECK_CLOSE_FRACTION(Ae(2,1), TYe(1,0), 1.0e-5);
   BOOST_CHECK_CLOSE_FRACTION(Ae(2,2), TYe(1,1), 2.5e-3);
   BOOST_CHECK_CLOSE_FRACTION(Ae(2,3), TYe(1,2), 1.0e-5);

   BOOST_CHECK_CLOSE_FRACTION(Ae(3,1), TYe(2,0), 1.0e-5);
   BOOST_CHECK_CLOSE_FRACTION(Ae(3,2), TYe(2,1), 1.0e-5);
   BOOST_CHECK_CLOSE_FRACTION(Ae(3,3), TYe(2,2), 2.5e-3);
}
