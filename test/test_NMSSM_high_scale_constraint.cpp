
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_NMSSM_high_scale_constraint

#include <boost/test/unit_test.hpp>

#define private public
#define protected public

#include "test_legacy.hpp"
#include "test_NMSSM.hpp"
#include "wrappers.hpp"
#include "ew_input.hpp"
#include "nmssmsoftsusy.h"
#include "NMSSM_two_scale_model.hpp"
#include "NMSSM_two_scale_high_scale_constraint.hpp"
#include "SoftsusyNMSSM_two_scale.hpp"
#include "SoftsusyNMSSM_parameter_point.hpp"
#include "SoftsusyNMSSM_two_scale_sugra_constraint.hpp"

using namespace flexiblesusy;
using namespace softsusy;

BOOST_AUTO_TEST_CASE( test_unification_condition )
{
   NMSSM_input_parameters input;
   input.m0 = 250.; // avoids tree-level tachyons
   NMSSM<Two_scale> m;
   NmssmSoftsusy s;
   setup_NMSSM(m, s, input);

   m.set_input_parameters(input);
   NMSSM_high_scale_constraint<Two_scale> constraint(&m);

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
   BOOST_CHECK_CLOSE_FRACTION(m.get_Lambdax(), input.LambdaInput, 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(m.get_mHd2(), Sqr(input.m0), 1.0e-09);
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
   BOOST_CHECK_CLOSE_FRACTION(m.get_TLambdax(), input.Azero * m.get_Lambdax(), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(m.get_TKappa()  , input.Azero * m.get_Kappa()  , 1.0e-10);
}

BOOST_AUTO_TEST_CASE( test_mx_calculation )
{
   NMSSM<Two_scale> m;
   NmssmSoftsusy softSusy;
   NMSSM_input_parameters input;
   setup_NMSSM(m, softSusy, input);
   SoftsusyNMSSM<Two_scale> s(softSusy);

   SoftsusyNMSSM_parameter_point pp;
   pp.tanBeta = input.TanBeta;
   pp.a0 = input.Azero;
   pp.m12 = input.m12;
   pp.m0 = input.m0;
   pp.lambda = input.LambdaInput;

   NMSSM_high_scale_constraint<Two_scale> NMSSM_sugra_constraint(&m);
   SoftsusyNMSSM_sugra_constraint snmssm_sugra_constraint(pp);

   snmssm_sugra_constraint.set_model((Model*)&s);

   NMSSM_sugra_constraint.update_scale();
   snmssm_sugra_constraint.update_scale();

   BOOST_CHECK_CLOSE_FRACTION(NMSSM_sugra_constraint.get_scale(),
                              snmssm_sugra_constraint.get_scale(), 1.0e-13);
}
