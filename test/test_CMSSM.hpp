
#ifndef TEST_CMSSM_H
#define TEST_CMSSM_H

#include "test_legacy.hpp"
#include "softsusy.h"
#include "CMSSM_two_scale_model.hpp"
#include "wrappers.hpp"
#include "ew_input.hpp"
#include "test_CMSSM_like.hpp"

using namespace flexiblesusy;
using namespace softsusy;

void ensure_tree_level_ewsb(CMSSM_mass_eigenstates& m)
{
   // ensure that the EWSB eqs. are satisfied (Drees p.222)
   const double vu = m.get_vu();
   const double vd = m.get_vd();
   const double gY = m.get_g1() * sqrt(0.6);
   const double g2 = m.get_g2();
   const double Mu = m.get_Mu();
   const double BMu = m.get_BMu();
   const double mHd2 = BMu*vu/vd - (sqr(gY) + sqr(g2))*(sqr(vd) - sqr(vu))/8. - sqr(Mu);
   const double mHu2 = BMu*vd/vu + (sqr(gY) + sqr(g2))*(sqr(vd) - sqr(vu))/8. - sqr(Mu);
   m.set_mHd2(mHd2);
   m.set_mHu2(mHu2);
}

void ensure_tree_level_ewsb(MssmSoftsusy& softSusy)
{
   const double Mu = softSusy.displaySusyMu();
   // const int signMu = Mu >= 0.0 ? 1 : -1;
   const double vev = softSusy.displayHvev();
   const double tanBeta = softSusy.displayTanb();
   const double beta = atan(tanBeta);
   const double sinBeta = sin(beta);
   const double cosBeta = cos(beta);
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double g1 = softSusy.displayGaugeCoupling(1);
   const double gY = g1 * sqrt(0.6);
   const double g2 = softSusy.displayGaugeCoupling(2);
   const double BMu = softSusy.displayM3Squared();
   const double mHd2 = BMu*vu/vd - (sqr(gY) + sqr(g2))*(sqr(vd) - sqr(vu))/8. - sqr(Mu);
   const double mHu2 = BMu*vd/vu + (sqr(gY) + sqr(g2))*(sqr(vd) - sqr(vu))/8. - sqr(Mu);
   const double MZrun = 0.5 * vev * sqrt(sqr(gY) + sqr(g2));

   softSusy.setMh1Squared(mHd2);
   softSusy.setMh2Squared(mHu2);

   TEST_CLOSE(MZrun, softSusy.displayMzRun(), 1.0e-10);
   TEST_CLOSE(-2 * BMu, (mHd2 - mHu2) * tan(2*beta) + sqr(MZrun) * sin(2*beta), 1.0e-10);
}

void ensure_n_loop_ewsb(CMSSM<Two_scale>& m, int loop_level)
{
   ensure_tree_level_ewsb(m);

   const double precision = m.get_ewsb_iteration_precision();
   m.set_ewsb_loop_order(loop_level);
   m.solve_ewsb();

   if (loop_level == 1) {
      TEST_CLOSE(m.get_ewsb_eq_hh_1() - m.tadpole_hh_1loop(0).real(), 0.0, 0.05);
      TEST_CLOSE(m.get_ewsb_eq_hh_2() - m.tadpole_hh_1loop(1).real(), 0.0, 0.15);
   }
}

void ensure_one_loop_ewsb(CMSSM<Two_scale>& m)
{
   ensure_n_loop_ewsb(m, 1);
}

void ensure_n_loop_ewsb(MssmSoftsusy& s, int loop_level)
{
   ensure_tree_level_ewsb(s);

   // s.rewsbTreeLevel(m.get_Mu() >= 0.0 ? 1 : -1);
   const int signMu = s.displaySusyMu() >= 0.0 ? 1 : -1;
   const double mtrun = s.displayDrBarPars().mt;
   const DoubleVector pars(3);
   softsusy::numRewsbLoops = loop_level;
   s.rewsb(signMu, mtrun, pars);
}

<<<<<<< HEAD
void setup_CMSSM_const(CMSSM<Two_scale>& m, 
                       const CMSSM_input_parameters& input)
=======

CMSSM_mass_eigenstates setup_CMSSM(const CMSSM_input_parameters& input)
>>>>>>> upstream/feature-mass_eigenstates_interface
{
   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5 * ALPHAMZ / (3 * (1 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = sqrt(4 * Pi * alpha1);
   const double g2 = sqrt(4 * Pi * alpha2);
   const double g3 = sqrt(4 * Pi * ALPHASMZ);
   const double tanBeta = input.TanBeta;
   const double sinBeta = sin(atan(tanBeta));
   const double cosBeta = cos(atan(tanBeta));
   const double M12 = input.m12;
   const double m0 = input.m0;
   const double a0 = input.Azero;
   const double root2 = sqrt(2.0);
   const double vev = 246.0;
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double susyMu = input.SignMu * 120.0;
   const double BMu = Sqr(2.0 * susyMu);
   const double scale = Electroweak_constants::MZ;

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero()),
      Yd(Eigen::Matrix<double,3,3>::Zero()),
      Ye(Eigen::Matrix<double,3,3>::Zero()),
      mm0(Eigen::Matrix<double,3,3>::Zero());
<<<<<<< HEAD
   Yu(0,0) = 0.002   * root2 / (vev * sinBeta);
   Yu(1,1) = 1.28   * root2 / (vev * sinBeta);
   Yu(2,2) = 165.0   * root2 / (vev * sinBeta);
   Yd(0,0) = 0.1     * root2 / (vev * cosBeta);
   Yd(1,1) = 4.2     * root2 / (vev * cosBeta);
   Yd(2,2) = 2.9     * root2 / (vev * cosBeta);
   Ye(0,0) = 5.11e-4 * root2 / (vev * cosBeta);
   Ye(1,1) = 0.105 * root2 / (vev * cosBeta);
   Ye(2,2) = 1.77699 * root2 / (vev * cosBeta);
   mm0 = Sqr(m0) * Eigen::Matrix<double,3,3>::Identity();

=======
   Yu(2,2) = 165.0   * root2 / (vev * sinBeta);
   Yd(2,2) = 2.9     * root2 / (vev * cosBeta);
   Ye(2,2) = 1.77699 * root2 / (vev * cosBeta);
   mm0 = Sqr(m0) * Eigen::Matrix<double,3,3>::Identity();

   CMSSM_mass_eigenstates m(input);

>>>>>>> upstream/feature-mass_eigenstates_interface
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

   ensure_tree_level_ewsb(m);
   m.calculate_DRbar_masses();
<<<<<<< HEAD
}

=======

   return m;
}


>>>>>>> upstream/feature-mass_eigenstates_interface
void setup_CMSSM_const(CMSSM<Two_scale>& m, MssmSoftsusy& s,
                       const CMSSM_input_parameters& input)
{
   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5 * ALPHAMZ / (3 * (1 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = sqrt(4 * Pi * alpha1);
   const double g2 = sqrt(4 * Pi * alpha2);
   const double g3 = sqrt(4 * Pi * ALPHASMZ);
   const double tanBeta = input.TanBeta;
   const double sinBeta = sin(atan(tanBeta));
   const double cosBeta = cos(atan(tanBeta));
   const double M12 = input.m12;
   const double m0 = input.m0;
   const double a0 = input.Azero;
   const double root2 = sqrt(2.0);
   const double vev = 246.0;
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double susyMu = input.SignMu * 120.0;
   const double BMu = Sqr(2.0 * susyMu);
   const double scale = Electroweak_constants::MZ;

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
   mm0 = Sqr(m0) * Eigen::Matrix<double,3,3>::Identity();

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

   s.setMu(scale);
   s.setLoops(1);
   s.setThresholds(3);
   s.setGaugeCoupling(1, g1);
   s.setGaugeCoupling(2, g2);
   s.setGaugeCoupling(3, g3);
   s.setYukawaMatrix(YU, Yu_SS);
   s.setYukawaMatrix(YD, Yd_SS);
   s.setYukawaMatrix(YE, Ye_SS);
   s.setGauginoMass(1, M12);
   s.setGauginoMass(2, M12);
   s.setGauginoMass(3, M12);
   s.setSoftMassMatrix(mQl, mm0_SS);
   s.setSoftMassMatrix(mUr, mm0_SS);
   s.setSoftMassMatrix(mDr, mm0_SS);
   s.setSoftMassMatrix(mLl, mm0_SS);
   s.setSoftMassMatrix(mEr, mm0_SS);
   s.setMh1Squared(sqr(m0));
   s.setMh2Squared(sqr(m0));
   s.setTrilinearMatrix(UA, a0 * Yu_SS);
   s.setTrilinearMatrix(DA, a0 * Yd_SS);
   s.setTrilinearMatrix(EA, a0 * Ye_SS);
   s.setSusyMu(susyMu);
   s.setM3Squared(BMu);
   s.setHvev(vev);
   s.setTanb(tanBeta);
   s.setMw(s.displayMwRun());

   ensure_tree_level_ewsb(m);
   m.calculate_DRbar_masses();

   ensure_tree_level_ewsb(s);
   s.calcDrBarPars();
}

void setup_CMSSM(CMSSM<Two_scale>& m, MssmSoftsusy& s,
                 CMSSM_input_parameters& input)
{
   input.m0 = 125.;
   input.m12 = 500.;
   input.TanBeta = 10.;
   input.SignMu = 1;
   input.Azero = 0.;

   setup_CMSSM_const(m, s, input);
}

void setup_CMSSM_const_non_3rd_gen(CMSSM<Two_scale>& m, MssmSoftsusy& s,
                                   const CMSSM_input_parameters& input)
{
   setup_CMSSM_const(m, s, input);

   const double ymu = 0.1;

   s.setYukawaElement(YE, 2, 2, ymu);

   Eigen::Matrix<double,3,3> Ye(m.get_Ye());
   Ye(1,1) = ymu;
   m.set_Ye(Ye);
}

void test_parameter_equality(const SoftParsMssm& a, const CMSSM_soft_parameters& b)
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

#endif
