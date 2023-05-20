
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSM_high_scale_constraint

#include <boost/test/unit_test.hpp>
#include "test.hpp"
#include <functional>
#include <Eigen/Dense>

#define private public
#define protected public

#include "softsusy.h"
#include "SoftsusyMSSM_parameter_point.hpp"
#include "SoftsusyMSSM_two_scale_sugra_constraint.hpp"
#include "test_CMSSM.hpp"
#include "CMSSM_two_scale_model.hpp"
#include "CMSSM_two_scale_high_scale_constraint.hpp"
#include "wrappers.hpp"

CMSSM<Two_scale> setup_CMSSM()
{
   CMSSM_input_parameters input;

   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5 * ALPHAMZ / (3 * (1 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = sqrt(4 * Pi * alpha1);
   const double g2 = sqrt(4 * Pi * alpha2);
   const double g3 = sqrt(4 * Pi * ALPHASMZ);
   const double tanBeta = 10;
   const double sinBeta = sin(atan(tanBeta));
   const double cosBeta = cos(atan(tanBeta));
   const double M12 = 100.0;
   const double m0 = 250.0;
   const double a0 = 50.0;
   const double root2 = sqrt(2.0);
   const double vev = 246.0;
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double susyMu = 120.0;
   const double BMu = Sqr(2.0 * susyMu);

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero()),
      Yd(Eigen::Matrix<double,3,3>::Zero()),
      Ye(Eigen::Matrix<double,3,3>::Zero()),
      mm0(Eigen::Matrix<double,3,3>::Zero());
   Yu(2,2) = 165.0   * root2 / (vev * sinBeta);
   Yd(2,2) = 2.9     * root2 / (vev * cosBeta);
   Ye(2,2) = 1.77699 * root2 / (vev * cosBeta);
   mm0 = Sqr(m0) * Eigen::Matrix<double,3,3>::Identity();

   CMSSM<Two_scale> m(input);
   m.set_scale(91);
   m.set_loops(1);
   m.set_g1(g1);
   m.set_g2(g2);
   m.set_g3(g3);
   m.set_Yu(Yu);
   m.set_Yd(Yd);
   m.set_Ye(Ye);
   m.set_MassB(M12);
   m.set_MassG(M12);
   m.set_MassWB(M12);
   m.set_mq2(mm0);
   m.set_ml2(mm0);
   m.set_md2(mm0);
   m.set_mu2(mm0);
   m.set_me2(mm0);
   m.set_mHd2(Sqr(m0));
   m.set_mHu2(Sqr(m0));
   m.set_TYu(a0 * Yu);
   m.set_TYd(a0 * Yd);
   m.set_TYe(a0 * Ye);
   m.set_Mu(susyMu);
   m.set_BMu(BMu);
   m.set_vu(vu);
   m.set_vd(vd);

   return m;
}

BOOST_AUTO_TEST_CASE( test_unification_condition )
{
   CMSSM<Two_scale> m(setup_CMSSM());

   CMSSM_input_parameters input;
   input.m0 = 100;
   input.m12 = 200;
   input.TanBeta = m.get_vd() / m.get_vu();
   input.Azero = 300;

   m.set_input_parameters(input);
   CMSSM_high_scale_constraint<Two_scale> constraint(&m);

   double mgut = constraint.get_scale(); // initial guess
   double mgut_new = mgut;
   double diff;
   int iteration = 0;

   do {
      m.run_to(mgut);
      constraint.apply();
      mgut_new = constraint.get_scale();
      diff = std::fabs((mgut - mgut_new)/mgut);
      mgut = mgut_new;
      iteration++;
   } while (diff > 1.0e-7 && iteration < 20);

   BOOST_CHECK_GT(1.0e-7, diff);
   BOOST_CHECK_GT(20, iteration);

   m.run_to(mgut_new);

   BOOST_CHECK_CLOSE_FRACTION(m.get_g1(), m.get_g2(), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(m.get_mHd2(), Sqr(input.m0), 1.0e-9);
   BOOST_CHECK_CLOSE_FRACTION(m.get_mHu2(), Sqr(input.m0), 1.0e-10);
   TEST_CLOSE(m.get_mq2(), Sqr(input.m0) * UNITMATRIX(3), 1.0e-6);
   TEST_CLOSE(m.get_ml2(), Sqr(input.m0) * UNITMATRIX(3), 1.0e-6);
   TEST_CLOSE(m.get_md2(), Sqr(input.m0) * UNITMATRIX(3), 1.0e-6);
   TEST_CLOSE(m.get_mu2(), Sqr(input.m0) * UNITMATRIX(3), 1.0e-6);
   TEST_CLOSE(m.get_me2(), Sqr(input.m0) * UNITMATRIX(3), 1.0e-6);
   BOOST_CHECK_CLOSE_FRACTION(m.get_MassB() , input.m12, 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(m.get_MassG() , input.m12, 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(m.get_MassWB(), input.m12, 1.0e-10);
   TEST_CLOSE(m.get_TYu(), input.Azero * m.get_Yu(), 1.0e-6);
   TEST_CLOSE(m.get_TYd(), input.Azero * m.get_Yd(), 1.0e-6);
   TEST_CLOSE(m.get_TYe(), input.Azero * m.get_Ye(), 1.0e-6);
}

BOOST_AUTO_TEST_CASE( test_mx_calculation )
{
   CMSSM<Two_scale> m; SoftsusyMSSM<Two_scale> s;
   CMSSM_input_parameters input;
   setup_CMSSM(m, s, input);
   SoftsusyMSSM_parameter_point pp;
   pp.tanBeta = input.TanBeta;
   pp.a0 = input.Azero;
   pp.m12 = input.m12;
   pp.m0 = input.m0;
   pp.signMu = input.SignMu;

   CMSSM_high_scale_constraint<Two_scale> CMSSM_sugra_constraint(&m);
   SoftsusyMSSM_sugra_constraint mssm_sugra_constraint(pp);

   mssm_sugra_constraint.set_model((Model*)&s);

   CMSSM_sugra_constraint.update_scale();
   mssm_sugra_constraint.update_scale();

   BOOST_CHECK_CLOSE_FRACTION(CMSSM_sugra_constraint.get_scale(),
                              mssm_sugra_constraint.get_scale(), 1.0e-13);
}
