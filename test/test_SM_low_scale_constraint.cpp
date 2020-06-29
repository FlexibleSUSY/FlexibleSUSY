
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SM_low_scale_constraint

#include <boost/test/unit_test.hpp>

#define private public

#include "config.h"
#include "SM_two_scale_model.hpp"
#include "SM_two_scale_low_scale_constraint.hpp"
#include "wrappers.hpp"
#include "ew_input.hpp"
#include "test_SM.hpp"

#include "standard_model.hpp"
#include "standard_model_two_scale_model.hpp"
#include "standard_model_two_scale_low_scale_constraint.hpp"

#define SARAH_VERSION_AT_LEAST(x,y,z) (                                 \
      (SARAH_MAJOR == 0 && SARAH_MINOR == 0 && SARAH_PATCH == 0) ||     \
      (SARAH_MAJOR > x || (SARAH_MAJOR >= x &&                          \
                           (SARAH_MINOR > y || (SARAH_MINOR >= y &&     \
                                                SARAH_PATCH >= z))))    \
      )

using namespace flexiblesusy;
using namespace softsusy;

BOOST_AUTO_TEST_CASE( test_model_properties )
{
   BOOST_CHECK(SM_info::is_supersymmetric_model == false);
}

BOOST_AUTO_TEST_CASE( test_delta_alpha )
{
   SM<Two_scale> m;
   SM_input_parameters input;
   input.LambdaIN = 0.1;
   QedQcd qedqcd;

   const double vev = 246.;

   m.set_scale(91.);
   m.set_v(246.);
   m.set_Yu(2, 2, 165.0   * Sqrt(2.) / vev);
   m.set_Yd(2, 2, 2.9     * Sqrt(2.) / vev);
   m.set_Ye(2, 2, 1.77699 * Sqrt(2.) / vev);

   m.calculate_DRbar_masses();

   SM_low_scale_constraint<Two_scale> constraint(&m, qedqcd);

   const double alpha_em = qedqcd.displayAlpha(ALPHA);
   const double alpha_s  = qedqcd.displayAlpha(ALPHAS);
   const double scale = m.get_scale();

   const double delta_alpha_em_fs = constraint.calculate_delta_alpha_em(alpha_em);
   const double delta_alpha_s_fs  = constraint.calculate_delta_alpha_s(alpha_s);

   const double Mtop = m.get_MFu(2);

   // the extra singlet field does not couple electromagnetically and
   // does thus not contribute here (or below)
   const double delta_alpha_em =
      alpha_em / (2 * Pi) * (- 16./9. * log(Mtop/scale));

   // no MS-bar DR-bar conversion term appears here, because the SM
   // is renormalized in MS-bar
   const double delta_alpha_s  =
      alpha_s / (2 * Pi) * (- 2./3. * log(Mtop/scale));

   BOOST_CHECK_CLOSE_FRACTION(delta_alpha_em_fs, delta_alpha_em, 1.0e-12);
   BOOST_CHECK_CLOSE_FRACTION(delta_alpha_s_fs , delta_alpha_s , 1.0e-12);
}

BOOST_AUTO_TEST_CASE( test_delta_spectrum )
{
   SM<Two_scale> m;
   SM_input_parameters input;
   input.LambdaIN = 0.1;
   QedQcd qedqcd;

   const double vev = 246.;
   const double g1 = 0.2;
   const double g2 = 0.3;
   const double g3 = 0.4;

   m.set_scale(91.);
   m.set_v(246.);
   m.set_g1(g1);
   m.set_g2(g2);
   m.set_g3(g3);
   m.set_Yu(2, 2, 165.0   * Sqrt(2.) / vev);
   m.set_Yd(2, 2, 2.9     * Sqrt(2.) / vev);
   m.set_Ye(2, 2, 1.77699 * Sqrt(2.) / vev);

   m.calculate_DRbar_masses();

   const double MVZ = 0.5 * vev * Sqrt(0.6*g1*g1 + g2*g2);

   BOOST_CHECK_CLOSE_FRACTION(m.get_MVZ(), MVZ, 1.0e-10);
}

BOOST_AUTO_TEST_CASE( test_low_scale_constraint )
{
   QedQcd qedqcd;
   SM_input_parameters input;
   input.LambdaIN = 0.25;
   SM<Two_scale> m;
   standard_model::StandardModel<Two_scale> sm;
   setup_SM_const(m, input);
   setup_SM_const(sm, input);

   m.set_thresholds(0);
   sm.set_thresholds(0);

   m.calculate_DRbar_masses();
   sm.calculate_DRbar_masses();

   SM_low_scale_constraint<Two_scale> c_m(&m, qedqcd);
   c_m.apply();
   standard_model::Standard_model_low_scale_constraint<Two_scale> c_sm(&sm, qedqcd);
   c_sm.apply();

   const double eps = 1e-14;

   BOOST_CHECK_CLOSE_FRACTION(m.get_g1(), sm.get_g1(), 1e-4); // deviation due to differen Weinberg angle class
   BOOST_CHECK_CLOSE_FRACTION(m.get_g2(), sm.get_g2(), 1e-4); // deviation due to differen Weinberg angle class
   BOOST_CHECK_CLOSE_FRACTION(m.get_g3(), sm.get_g3(), eps);
   BOOST_CHECK_CLOSE_FRACTION(m.get_Lambdax(), sm.get_Lambdax(), eps);
   BOOST_CHECK_CLOSE_FRACTION(m.get_v(), sm.get_v(), eps);

   // Note: in SARAH v4.13.0 the definition of mu2 was changed by an
   // overall sign
#if SARAH_VERSION_AT_LEAST(4,13,0)
   BOOST_CHECK_CLOSE_FRACTION(-m.get_mu2(), sm.get_mu2(), eps);
#else
   BOOST_CHECK_CLOSE_FRACTION(m.get_mu2(), sm.get_mu2(), eps);
#endif

   for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
         BOOST_CHECK_CLOSE_FRACTION(m.get_Yu(i,k), -sm.get_Yu(i,k), eps);
         BOOST_CHECK_CLOSE_FRACTION(m.get_Yd(i,k), sm.get_Yd(i,k), eps);
         BOOST_CHECK_CLOSE_FRACTION(m.get_Ye(i,k), sm.get_Ye(i,k), eps);
      }
   }
}

BOOST_AUTO_TEST_CASE( test_initialise_from_input )
{
   QedQcd qedqcd;
   SM_input_parameters input;
   input.LambdaIN = 0.25;
   standard_model::StandardModel<Two_scale> m1, m2;
   setup_SM_const(m1, input);
   setup_SM_const(m2, input);

   m1.set_thresholds(3);
   m2.set_thresholds(3);

   m1.set_scale(qedqcd.displayPoleMZ());
   m2.set_scale(qedqcd.displayPoleMZ());

   m1.calculate_DRbar_masses();
   m2.calculate_DRbar_masses();

   m2.initialise_from_input(qedqcd);

   // initialize after m2.initialise_from_input()
   qedqcd.to(qedqcd.displayPoleMZ());

   for (int i = 0; i < 10; i++) {
      standard_model::Standard_model_low_scale_constraint<Two_scale> c_m1(&m1, qedqcd);
      c_m1.apply();
      m1.solve_ewsb();
      m1.calculate_Lambdax_DRbar();
   }

   const double eps = 1e-4;

   BOOST_CHECK_CLOSE_FRACTION(m1.get_scale(), m2.get_scale(), 1e-14);

   BOOST_CHECK_CLOSE_FRACTION(m1.get_g1(), m2.get_g1(), 1e-5);
   BOOST_CHECK_CLOSE_FRACTION(m1.get_g2(), m2.get_g2(), eps);
   BOOST_CHECK_CLOSE_FRACTION(m1.get_g3(), m2.get_g3(), 1e-6);
   BOOST_CHECK_CLOSE_FRACTION(m1.get_Lambdax(), m2.get_Lambdax(), 1e-3);
   BOOST_CHECK_CLOSE_FRACTION(m1.get_v(), m2.get_v(), 1e-3);
   BOOST_CHECK_CLOSE_FRACTION(m1.get_mu2(), m2.get_mu2(), 1e-4);

   BOOST_CHECK_CLOSE_FRACTION(m1.get_Yu(0,0), m2.get_Yu(0,0), eps);
   BOOST_CHECK_CLOSE_FRACTION(m1.get_Yd(0,0), m2.get_Yd(0,0), eps);
   BOOST_CHECK_CLOSE_FRACTION(m1.get_Ye(0,0), m2.get_Ye(0,0), eps);

   BOOST_CHECK_CLOSE_FRACTION(m1.get_Yu(1,1), m2.get_Yu(1,1), eps);
   BOOST_CHECK_CLOSE_FRACTION(m1.get_Yd(1,1), m2.get_Yd(1,1), eps);
   BOOST_CHECK_CLOSE_FRACTION(m1.get_Ye(1,1), m2.get_Ye(1,1), eps);

   BOOST_CHECK_CLOSE_FRACTION(m1.get_Yu(2,2), m2.get_Yu(2,2), eps);
   BOOST_CHECK_CLOSE_FRACTION(m1.get_Yd(2,2), m2.get_Yd(2,2), eps);
   BOOST_CHECK_CLOSE_FRACTION(m1.get_Ye(2,2), m2.get_Ye(2,2), eps);
}
