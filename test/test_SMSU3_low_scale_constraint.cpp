
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SMSU3_low_scale_constraint

#include <boost/test/unit_test.hpp>
#include <functional>
#include <Eigen/Dense>

#define private public

#include "SMSU3_two_scale_model.hpp"
#include "SMSU3_two_scale_low_scale_constraint.hpp"
#include "wrappers.hpp"
#include "ew_input.hpp"

using namespace flexiblesusy;

using DV3 = Eigen::Matrix<double,3,1>;

SMSU3<Two_scale> setup_SMSU3(const SMSU3_input_parameters& input)
{
   SMSU3<Two_scale> m;
   m.set_thresholds(1);

   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5 * ALPHAMZ / (3 * (1 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = sqrt(4 * Pi * alpha1);
   const double g2 = sqrt(4 * Pi * alpha2);
   const double g3 = sqrt(4 * Pi * ALPHASMZ);
   const double root2 = sqrt(2.0);
   const double vev = 246.0;
   const double scale = Electroweak_constants::MZ;

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero()),
      Yd(Eigen::Matrix<double,3,3>::Zero()),
      Ye(Eigen::Matrix<double,3,3>::Zero());
   Yu(2,2) = 165.0   * root2 / vev;
   Yd(2,2) = 2.9     * root2 / vev;
   Ye(2,2) = 1.77699 * root2 / vev;

   m.set_input_parameters(input);
   m.set_scale(scale);
   m.set_loops(1);
   m.set_thresholds(3);
   m.set_g1(g1);
   m.set_g2(g2);
   m.set_g3(g3);
   m.set_Yu(Yu);
   m.set_Yd(Yd);
   m.set_Ye(Ye);
   m.set_mu2(1);
   m.set_v(vev);
   m.set_Lambdax(input.LambdaIN);
   m.set_mF1(input.mFIN);
   m.set_mF2(input.mFIN);
   m.set_mF3(input.mFIN);

   m.solve_ewsb_tree_level();
   m.calculate_DRbar_masses();

   return m;
}

DV3 calculate_gauge_couplings(
   SMSU3<Two_scale> model,
   SMSU3_low_scale_constraint<Two_scale> constraint,
   double scale)
{
   model.set_thresholds(1);
   model.set_scale(scale);
   constraint.set_model(&model);
   constraint.apply();

   DV3 g;
   g(0) = model.get_g1();
   g(1) = model.get_g2();
   g(2) = model.get_g3();

   return g;
}

BOOST_AUTO_TEST_CASE( test_threshold_corrections )
{
   SMSU3_input_parameters input;
   input.LambdaIN = 0.2;
   input.mFIN = 100.;
   input.Qin = 1000.;
   input.QEWSB = 173.34;
   softsusy::QedQcd qedqcd;
   SMSU3<Two_scale> m = setup_SMSU3(input);

   SMSU3_low_scale_constraint<Two_scale> constraint(&m, qedqcd);

   const double Q1 = constraint.get_scale();
   const double Q2 = 1.1 * Q1;
   DV3 g_old;
   g_old(0) = Electroweak_constants::g1;
   g_old(1) = Electroweak_constants::g2;
   g_old(2) = Electroweak_constants::g3;
   DV3 prefactor;
   for (int i = 0; i < 3; i++)
      prefactor(i) = 1. / (oneOver16PiSqr * Power(g_old(i),3));

   const DV3 g_Q1 = calculate_gauge_couplings(m, constraint, Q1);
   const DV3 g_Q2 = calculate_gauge_couplings(m, constraint, Q2);

   const DV3 beta_numeric = (g_Q1 - g_Q2).cwiseProduct(prefactor) * (1. / log(Q1/Q2));
   DV3 beta_SM;
   beta_SM(0) = 41./6.;
   beta_SM(1) = -19./6.;
   beta_SM(2) = -7. - 2./3; // remove top quark
   DV3 beta_SMSU3;
   beta_SMSU3(0) = 11.3;
   beta_SMSU3(1) = -3.1666666666666665;
   beta_SMSU3(2) = -5.;

   // BOOST_CHECK_CLOSE_FRACTION(beta_numeric(0), beta_SMSU3(0) - beta_SM(0), 0.075);
   // BOOST_CHECK_CLOSE_FRACTION(beta_numeric(1), beta_SMSU3(1) - beta_SM(1), 0.075);
   BOOST_CHECK_CLOSE_FRACTION(beta_numeric(2), beta_SMSU3(2) - beta_SM(2), 0.03);
}
