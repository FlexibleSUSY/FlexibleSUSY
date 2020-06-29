#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_sfermions

#include <boost/test/unit_test.hpp>

#include "CMSSM_two_scale_model.hpp"
#include "softsusy.h"
#include "sfermions.hpp"
#include "wrappers.hpp"
#include "test_CMSSM.hpp"

using namespace flexiblesusy;

void set_mssm_parameters(CMSSM<Two_scale>& m, MssmSoftsusy& softSusy)
{
   const int loopLevel = 1;
   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5 * ALPHAMZ / (3 * (1 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = sqrt(4 * PI * alpha1);
   const double g2 = sqrt(4 * PI * alpha2);
   const double g3 = sqrt(4 * PI * ALPHASMZ);
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
   const double BMu = sqr(2.0 * susyMu);
   DoubleMatrix Yu_SS(3,3), Yd_SS(3,3), Ye_SS(3,3);
   Yu_SS(3,3) = 165.0   * root2 / (vev * sinBeta);
   Yd_SS(3,3) = 2.9     * root2 / (vev * cosBeta);
   Ye_SS(3,3) = 1.77699 * root2 / (vev * cosBeta);
   DoubleMatrix ID(3, 3), mm0_SS(3, 3);
   for (int i=1; i<=3; i++) ID(i, i) = 1.0;
   mm0_SS = ID * sqr(m0);

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero()),
      Yd(Eigen::Matrix<double,3,3>::Zero()),
      Ye(Eigen::Matrix<double,3,3>::Zero()),
      mm0(Eigen::Matrix<double,3,3>::Zero());
   Yu(2,2) = 165.0   * root2 / (vev * sinBeta);
   Yd(2,2) = 2.9     * root2 / (vev * cosBeta);
   Ye(2,2) = 1.77699 * root2 / (vev * cosBeta);
   mm0 = sqr(m0) * Eigen::Matrix<double,3,3>::Identity();

   m.set_scale(91);
   m.set_loops(loopLevel);
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
   m.set_mHd2(sqr(m0));
   m.set_mHu2(sqr(m0));
   m.set_TYu(a0 * Yu);
   m.set_TYd(a0 * Yd);
   m.set_TYe(a0 * Ye);
   m.set_Mu(susyMu);
   m.set_BMu(BMu);
   m.set_vu(vu);
   m.set_vd(vd);

   softSusy.setMu(91);
   softSusy.setLoops(loopLevel);
   softSusy.setGaugeCoupling(1, g1);
   softSusy.setGaugeCoupling(2, g2);
   softSusy.setGaugeCoupling(3, g3);
   softSusy.setYukawaMatrix(YU, Yu_SS);
   softSusy.setYukawaMatrix(YD, Yd_SS);
   softSusy.setYukawaMatrix(YE, Ye_SS);
   softSusy.setGauginoMass(1, M12);
   softSusy.setGauginoMass(2, M12);
   softSusy.setGauginoMass(3, M12);
   softSusy.setSoftMassMatrix(mQl, mm0_SS);
   softSusy.setSoftMassMatrix(mUr, mm0_SS);
   softSusy.setSoftMassMatrix(mDr, mm0_SS);
   softSusy.setSoftMassMatrix(mLl, mm0_SS);
   softSusy.setSoftMassMatrix(mEr, mm0_SS);
   softSusy.setMh1Squared(sqr(m0));
   softSusy.setMh2Squared(sqr(m0));
   softSusy.setTrilinearMatrix(UA, a0 * Yu_SS);
   softSusy.setTrilinearMatrix(DA, a0 * Yd_SS);
   softSusy.setTrilinearMatrix(EA, a0 * Ye_SS);
   softSusy.setSusyMu(susyMu);
   softSusy.setM3Squared(BMu);
   softSusy.setHvev(vev);
   softSusy.setTanb(tanBeta);
}

void setup_mssm_models(CMSSM<Two_scale>& m, MssmSoftsusy& softSusy)
{
   set_mssm_parameters(m, softSusy);

   ensure_tree_level_ewsb(m);
   m.calculate_DRbar_masses();
   softSusy.calcDrBarPars();
}

BOOST_AUTO_TEST_CASE( test_stop )
{
   CMSSM<Two_scale> m;
   MssmSoftsusy s;
   setup_mssm_models(m, s);

   DoubleVector mf(2);
   mf(1) = s.displayDrBarPars().mu(1,3);
   mf(2) = s.displayDrBarPars().mu(2,3);
   const double thetat = s.displayDrBarPars().thetat;

   Eigen::Array<double,2,1> mstop;
   double theta_stop;

   m.calculate_MSu_3rd_generation(mstop(0), mstop(1), theta_stop);

   BOOST_CHECK_CLOSE_FRACTION(mstop(0), mf(1), 1.0e-9);
   BOOST_CHECK_CLOSE_FRACTION(mstop(1), mf(2), 1.0e-9);
   BOOST_CHECK_CLOSE_FRACTION(theta_stop, thetat, 1.0e-9);
}

BOOST_AUTO_TEST_CASE( test_sbottom )
{
   CMSSM<Two_scale> m;
   MssmSoftsusy s;
   setup_mssm_models(m, s);

   DoubleVector mf(2);
   mf(1) = s.displayDrBarPars().md(1,3);
   mf(2) = s.displayDrBarPars().md(2,3);
   double thetab = s.displayDrBarPars().thetab;

   Eigen::Array<double,2,1> msbottom;
   double theta_sbottom;

   m.calculate_MSd_3rd_generation(msbottom(0), msbottom(1), theta_sbottom);

   BOOST_CHECK_CLOSE_FRACTION(msbottom(0), mf(1), 1.0e-9);
   BOOST_CHECK_CLOSE_FRACTION(msbottom(1), mf(2), 1.0e-9);
   BOOST_CHECK_CLOSE_FRACTION(theta_sbottom, thetab, 1.0e-9);
}

BOOST_AUTO_TEST_CASE( test_stau )
{
   CMSSM<Two_scale> m;
   MssmSoftsusy s;
   setup_mssm_models(m, s);

   DoubleVector mf(2);
   mf(1) = s.displayDrBarPars().me(1,3);
   mf(2) = s.displayDrBarPars().me(2,3);
   double thetatau = s.displayDrBarPars().thetatau;

   Eigen::Array<double,2,1> mstau;
   double theta_stau;

   m.calculate_MSe_3rd_generation(mstau(0), mstau(1), theta_stau);

   BOOST_CHECK_CLOSE_FRACTION(mstau(0), mf(1), 1.0e-9);
   BOOST_CHECK_CLOSE_FRACTION(mstau(1), mf(2), 1.0e-9);
   BOOST_CHECK_CLOSE_FRACTION(theta_stau, thetatau, 1.0e-9);
}

BOOST_AUTO_TEST_CASE( test_snu )
{
   CMSSM<Two_scale> m;
   MssmSoftsusy s;
   setup_mssm_models(m, s);

   double mf = s.displayDrBarPars().msnu(3);

   Eigen::Array<double,2,1> msnu;
   double theta_snu;

   m.calculate_MSv_3rd_generation(msnu(0), msnu(1), theta_snu);

   BOOST_CHECK_CLOSE_FRACTION(msnu(0), 0., 1.0e-9);
   BOOST_CHECK_CLOSE_FRACTION(msnu(1), mf, 1.0e-9);
   BOOST_CHECK_CLOSE_FRACTION(theta_snu, 0.5 * Pi, 1.0e-9);
}

// calculates theta according to Eq. (19) of arxiv:hep-ph/0105096
double calculate_theta_stop_field_independent(const CMSSM<Two_scale>& m)
{
   double mst1, mst2, theta_stop;
   m.calculate_MSu_3rd_generation(mst1, mst2, theta_stop);

   const double mt = m.get_MFu(2);
   const double At = m.get_TYu(2,2) / m.get_Yu(2,2);
   const double Mu = m.get_Mu();
   const double tanb = m.get_vu() / m.get_vd();
   const double signMu = -1;
   const double Xt = At + signMu * Mu / tanb;
   const double s2t = 2 * mt * Xt / (Sqr(mst1) - Sqr(mst2));

   return 0.5 * ArcSin(s2t);
}

BOOST_AUTO_TEST_CASE( test_stop_different_sign )
{
   CMSSM<Two_scale> m;
   MssmSoftsusy s;
   set_mssm_parameters(m, s);
   m.set_Mu(1200.0);
   s.setSusyMu(1200.0);
   ensure_tree_level_ewsb(m);
   m.calculate_DRbar_masses();
   s.calcDrBarPars();

   DoubleVector mf(2);
   mf(1) = s.displayDrBarPars().mu(1,3);
   mf(2) = s.displayDrBarPars().mu(2,3);
   const double thetat = s.displayDrBarPars().thetat;

   Eigen::Array<double,2,1> mstop;
   double theta_stop;

   m.calculate_MSu_3rd_generation(mstop(0), mstop(1), theta_stop);

   // calculate theta according to Eq. (19) of arxiv:hep-ph/0105096
   const double theta_stop_2 = calculate_theta_stop_field_independent(m);

   BOOST_CHECK_CLOSE_FRACTION(mstop(0), mf(1), 1.0e-9);
   BOOST_CHECK_CLOSE_FRACTION(mstop(1), mf(2), 1.0e-9);
   BOOST_CHECK_CLOSE_FRACTION(theta_stop, thetat, 1.0e-9);
   BOOST_CHECK_CLOSE_FRACTION(theta_stop, theta_stop_2, 1.0e-9);
}
