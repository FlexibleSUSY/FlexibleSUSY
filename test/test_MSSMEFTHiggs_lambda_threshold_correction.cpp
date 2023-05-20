#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSMEFTHiggs_lambda_threshold_correction

#include <boost/test/unit_test.hpp>

#include "wrappers.hpp"
#include "standard_model.hpp"
#include "splitmssm_thresholds.hpp"
#include "models/MSSMEFTHiggs/MSSMEFTHiggs_standard_model_matching.hpp"

using namespace flexiblesusy;
using namespace flexiblesusy::standard_model;

Standard_model create_sm(double Q)
{
   Eigen::Matrix<double,3,3> Yu, Yd, Ye;
   Yu.setZero();
   Yd.setZero();
   Ye.setZero();
   Yu(2,2) = 0.9;

   Standard_model sm;
   sm.set_scale(Q);
   sm.set_g1(0.3);
   sm.set_g2(0.4);
   sm.set_g3(0.9);
   sm.set_Yu(Yu);
   sm.set_Yd(Yd);
   sm.set_Ye(Ye);
   sm.set_Lambdax(0.1);
   sm.set_v(0.1);

   sm.calculate_DRbar_masses();
   sm.solve_ewsb();

   return sm;
}

MSSMEFTHiggs_mass_eigenstates create_mssm(double MS, double Xt, double tb, double Q)
{
   Eigen::Matrix<double,3,3> Yu, Yd, Ye, Id;
   Id.setIdentity();
   Yu.setZero();
   Yd.setZero();
   Ye.setZero();
   Yu(2,2) = 0.9;

   const double v = 0.1;
   const double beta = ArcTan(tb);
   const double Mu = MS;
   const double Xb = 0;
   const double Xl = 0;
   const double At = Xt*MS + Mu/tb;
   const double Ab = Xb*MS + Mu*tb;
   const double Al = Xl*MS + Mu*tb;

   MSSMEFTHiggs_mass_eigenstates mssm;
   mssm.set_scale(Q);
   mssm.set_g1(0.3);
   mssm.set_g2(0.4);
   mssm.set_g3(0.9);
   mssm.set_Yu(Yu);
   mssm.set_Yd(Yd);
   mssm.set_Ye(Ye);
   mssm.set_vu(v * Sin(beta));
   mssm.set_vd(v * Cos(beta));
   mssm.set_Mu(Mu);
   mssm.set_BMu(Sqr(MS)*(Sin(beta)*Cos(beta)));
   mssm.set_TYu(At * Yu);
   mssm.set_TYd(Ab * Yd);
   mssm.set_TYe(Al * Ye);
   mssm.set_mq2(Sqr(MS) * Id);
   mssm.set_ml2(Sqr(MS) * Id);
   mssm.set_mu2(Sqr(MS) * Id);
   mssm.set_md2(Sqr(MS) * Id);
   mssm.set_me2(Sqr(MS) * Id);
   mssm.set_MassB(MS);
   mssm.set_MassWB(MS);
   mssm.set_MassG(MS);

   mssm.calculate_DRbar_masses();
   mssm.solve_ewsb();

   return mssm;
}

double calculate_delta_lambda_eft(const Standard_model& sm, const MSSMEFTHiggs_mass_eigenstates& mssm, unsigned loops)
{
   using namespace splitmssm_thresholds;

   Parameters pars;
   pars.g1 = sm.get_g1();
   pars.g2 = sm.get_g2();
   pars.g3 = sm.get_g3();
   pars.gt = sm.get_Yu(2,2);
   pars.At = mssm.get_TYu(2,2) / mssm.get_Yu(2,2);
   pars.mu = mssm.get_Mu();
   pars.mA = mssm.get_MAh(1);
   pars.m1 = mssm.get_MassB();
   pars.m2 = mssm.get_MassWB();
   pars.tan_beta = mssm.get_vu() / mssm.get_vd();
   pars.scale = mssm.get_scale();
   pars.mq2 = mssm.get_mq2();
   pars.mu2 = mssm.get_mu2();
   pars.md2 = mssm.get_md2();
   pars.ml2 = mssm.get_ml2();
   pars.me2 = mssm.get_me2();

   double lambda = lambda_tree_level(pars);

   if (loops > 0) {
      lambda +=
         + delta_lambda_1loop_reg(pars)
         + delta_lambda_1loop_phi(pars)
         + delta_lambda_1loop_chi_1(pars)
         + delta_lambda_1loop_chi_2(pars);
   }

   return lambda;
}

double calculate_delta_lambda_fs(const Standard_model& sm, const MSSMEFTHiggs_mass_eigenstates& mssm, unsigned loops)
{
   using namespace MSSMEFTHiggs_standard_model_matching;

   auto tmp_sm = sm;
   auto tmp_mssm = mssm;

   if (loops == 0) {
      match_low_to_high_scale_model_tree_level(tmp_mssm, tmp_sm);
      match_high_to_low_scale_model_tree_level(tmp_sm, tmp_mssm, 0);
   } else {
      match_low_to_high_scale_model(tmp_mssm, tmp_sm, 1, 0);
      match_high_to_low_scale_model(tmp_sm, tmp_mssm, 1, 0);
   }

   return tmp_sm.get_Lambdax();
}

BOOST_AUTO_TEST_CASE( test_threshold_0L )
{
   Standard_model sm = create_sm(3000.);
   MSSMEFTHiggs_mass_eigenstates mssm = create_mssm(2000., -2., 5, 3000.);

   const double dl_eft = calculate_delta_lambda_eft(sm, mssm, 0);
   const double dl_fs  = calculate_delta_lambda_fs(sm, mssm, 0);

   BOOST_TEST_MESSAGE("Delta lambda^(0) EFT             : " << dl_eft);
   BOOST_TEST_MESSAGE("Delta lambda^(0) FlexibleEFTHiggs: " << dl_fs);

   BOOST_CHECK_CLOSE_FRACTION(dl_eft, dl_fs, 1e-6);
}

BOOST_AUTO_TEST_CASE( test_threshold_1L )
{
   Standard_model sm = create_sm(2000.);
   MSSMEFTHiggs_mass_eigenstates mssm = create_mssm(2000., -2., 5, 2000.);

   const double dl_eft = calculate_delta_lambda_eft(sm, mssm, 1);
   const double dl_fs  = calculate_delta_lambda_fs(sm, mssm, 1);

   BOOST_TEST_MESSAGE("Delta lambda^(1) EFT             : " << dl_eft);
   BOOST_TEST_MESSAGE("Delta lambda^(1) FlexibleEFTHiggs: " << dl_fs);

   BOOST_CHECK_CLOSE_FRACTION(dl_eft, dl_fs, 1e-2);
}
