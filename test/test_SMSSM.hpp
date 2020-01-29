
#ifndef TEST_SMSSM_H
#define TEST_SMSSM_H

#include "wrappers.hpp"
#include "ew_input.hpp"
#include "nmssmsoftsusy.h"
#include "SMSSM_two_scale_model.hpp"

using namespace flexiblesusy;
using namespace softsusy;

void setup_SMSSM_const(SMSSM<Two_scale>& m, NmssmSoftsusy& s, const SMSSM_input_parameters& input)
{
   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5 * ALPHAMZ / (3 * (1 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = sqrt(4 * Pi * alpha1);
   const double g2 = sqrt(4 * Pi * alpha2);
   const double g3 = sqrt(4 * Pi * ALPHASMZ);
   const double lambda = input.LambdaInput;
   const double kappa = input.KappaInput;
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
   const double vS = input.LambdaSInput / lambda;
   const double susyMu = MZ;
   const double BMu = MZ * MZ;
   const double scale = Electroweak_constants::MZ;
   // Z3 violating terms
   const double MS = input.MSInput;
   const double BMS = input.BMSInput;
   const double L1 = input.L1Input;
   const double LL1 = m0 * L1;

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

   m.set_scale(scale);
   m.set_loops(1);
   m.set_g1(g1);
   m.set_g2(g2);
   m.set_g3(g3);
   m.set_Yu(Yu);
   m.set_Yd(Yd);
   m.set_Ye(Ye);
   m.set_Lambdax(lambda);
   m.set_Kappa(kappa);
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
   m.set_ms2(Sqr(m0));
   m.set_TYu(a0 * Yu);
   m.set_TYd(a0 * Yd);
   m.set_TYe(a0 * Ye);
   m.set_TLambdax(a0 * lambda);
   m.set_TKappa(a0 * kappa);
   m.set_Mu(susyMu); // EWSB eqs. output
   m.set_BMu(BMu);   // EWSB eqs. output
   m.set_vu(vu);
   m.set_vd(vd);
   m.set_vS(vS);
   // Z3 violating terms
   m.set_MS(MS);
   m.set_BMS(BMS);
   m.set_L1(L1);
   m.set_LL1(LL1);

   s.setMu(scale);
   s.setLoops(1);
   s.setGaugeCoupling(1, g1);
   s.setGaugeCoupling(2, g2);
   s.setGaugeCoupling(3, g3);
   s.setYukawaMatrix(YU, Yu_SS);
   s.setYukawaMatrix(YD, Yd_SS);
   s.setYukawaMatrix(YE, Ye_SS);
   s.setLambda(lambda);
   s.setKappa(kappa);
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
   s.setMsSquared(sqr(m0));
   s.setTrilinearMatrix(UA, a0 * Yu_SS);
   s.setTrilinearMatrix(DA, a0 * Yd_SS);
   s.setTrilinearMatrix(EA, a0 * Ye_SS);
   s.setTrialambda(a0 * lambda);
   s.setTriakappa(a0 * kappa);
   s.setSusyMu(susyMu); // EWSB eqs. output
   s.setM3Squared(BMu); // EWSB eqs. output
   s.setHvev(vev);
   s.setTanb(tanBeta);
   s.setSvev(vS);
   // Z3 violating terms
   s.setMupr(MS);
   s.setMspSquared(BMS);
   s.setXiF(L1);
   s.setXiS(LL1);

   s.setMw(s.displayMwRun());
}

void setup_SMSSM(SMSSM<Two_scale>& m, NmssmSoftsusy& s, SMSSM_input_parameters& input)
{
   input.m0 = 200.;
   input.m12 = 200.;
   input.TanBeta = 10.;
   input.SignMu = 1;
   input.Azero = 0.;
   input.LambdaInput = 0.1;
   input.KappaInput = 0.1;
   input.LambdaSInput = 100.;
   input.L1Input = 0.;
   input.MSInput = 0.;
   input.BMSInput = 0.;

   setup_SMSSM_const(m, s, input);
}

#endif
