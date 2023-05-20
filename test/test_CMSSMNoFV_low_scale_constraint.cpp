
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSMNoFV_low_scale_constraint

#include <boost/test/unit_test.hpp>

#define private public

#include "test_CMSSMNoFV.hpp"
#include "CMSSM_two_scale_model.hpp"
#include "CMSSM_two_scale_low_scale_constraint.hpp"
#include "CMSSMNoFV_two_scale_model.hpp"
#include "CMSSMNoFV_two_scale_low_scale_constraint.hpp"
#include "wrappers.hpp"
#include "ew_input.hpp"
#include "lowe.h"

using namespace softsusy;

BOOST_AUTO_TEST_CASE( test_delta_alpha )
{
   CMSSM<Two_scale> mssm;
   CMSSMNoFV<Two_scale> mssmnofv;
   CMSSM_input_parameters input_mssm;
   input_mssm.TanBeta = 10.;
   input_mssm.m0 = 125.;
   input_mssm.m12 = 200.;
   input_mssm.SignMu = 1;
   input_mssm.Azero = 0.;

   CMSSMNoFV_input_parameters input_mssmnofv;
   input_mssmnofv.TanBeta = 10.;
   input_mssmnofv.m0 = 125.;
   input_mssmnofv.m12 = 200.;
   input_mssmnofv.SignMu = 1;
   input_mssmnofv.Azero = 0.;

   QedQcd qedqcd;
   setup_CMSSM_models(mssm, mssmnofv, input_mssmnofv);

   mssm.calculate_DRbar_masses();
   mssmnofv.calculate_DRbar_masses();

   CMSSM_low_scale_constraint<Two_scale> constraint_mssm(&mssm, qedqcd);

   CMSSMNoFV_low_scale_constraint<Two_scale> constraint_mssmnofv(&mssmnofv,qedqcd);

   const double alpha_em = qedqcd.displayAlpha(ALPHA);
   const double alpha_s  = qedqcd.displayAlpha(ALPHAS);
   const double scale = mssm.get_scale();

   const double delta_alpha_em_mssm = constraint_mssm.calculate_delta_alpha_em(alpha_em);
   const double delta_alpha_s_mssm  = constraint_mssm.calculate_delta_alpha_s(alpha_s);

   const double delta_alpha_em_mssmnofv = constraint_mssmnofv.calculate_delta_alpha_em(alpha_em);
   const double delta_alpha_s_mssmnofv  = constraint_mssmnofv.calculate_delta_alpha_s(alpha_s);

   BOOST_CHECK_CLOSE_FRACTION(delta_alpha_em_mssm, delta_alpha_em_mssmnofv, 1.0e-12);
   BOOST_CHECK_CLOSE_FRACTION(delta_alpha_s_mssm , delta_alpha_s_mssmnofv , 1.0e-12);
}

BOOST_AUTO_TEST_CASE( test_delta_Yf )
{
   CMSSM<Two_scale> mssm;
   mssm.set_thresholds(1);
   CMSSMNoFV<Two_scale> mssmnofv;
   mssmnofv.set_thresholds(1);

   CMSSM_input_parameters input_mssm;
   input_mssm.TanBeta = 10.;
   input_mssm.m0 = 125.;
   input_mssm.m12 = 200.;
   input_mssm.SignMu = 1;
   input_mssm.Azero = 0.;

   CMSSMNoFV_input_parameters input_mssmnofv;
   input_mssmnofv.TanBeta = 10.;
   input_mssmnofv.m0 = 125.;
   input_mssmnofv.m12 = 200.;
   input_mssmnofv.SignMu = 1;
   input_mssmnofv.Azero = 0.;

   QedQcd qedqcd;
   setup_CMSSM_models(mssm, mssmnofv, input_mssmnofv);

   mssm.calculate_DRbar_masses();
   mssmnofv.calculate_DRbar_masses();

   CMSSM_low_scale_constraint<Two_scale> constraint_mssm(&mssm, qedqcd);

   CMSSMNoFV_low_scale_constraint<Two_scale> constraint_mssmnofv(&mssmnofv,qedqcd);

   constraint_mssmnofv.calculate_DRbar_yukawa_couplings();
   constraint_mssm.calculate_DRbar_yukawa_couplings();

   for (int i = 0; i < 3; i++) {
      BOOST_CHECK_CLOSE_FRACTION(mssm.get_Yu(i,i), mssmnofv.get_Yu(i,i), 1.0e-12);
      BOOST_CHECK_CLOSE_FRACTION(mssm.get_Yd(i,i), mssmnofv.get_Yd(i,i), 1.0e-12);
      BOOST_CHECK_CLOSE_FRACTION(mssm.get_Ye(i,i), mssmnofv.get_Ye(i,i), 1.0e-12);
   }
}
