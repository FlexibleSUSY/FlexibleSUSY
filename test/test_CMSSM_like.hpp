
#ifndef TEST_CMSSM_LIKE_H
#define TEST_CMSSM_LIKE_H

#include "test_legacy.hpp"
#include "softsusy.h"
#include "wrappers.hpp"
#include "ew_input.hpp"

using namespace flexiblesusy;
using namespace softsusy;

template <class MSSMModel>
void copy_parameters(const MSSMModel& mssm, MssmSoftsusy& softsusy)
{
   // copy base class parameters
   softsusy.setLoops(mssm.get_loops());
   softsusy.setMu(mssm.get_scale());
   softsusy.setThresholds(mssm.get_thresholds());

   // copy susy parameters
   softsusy.setGaugeCoupling(1, mssm.get_g1());
   softsusy.setGaugeCoupling(2, mssm.get_g2());
   softsusy.setGaugeCoupling(3, mssm.get_g3());

   const double vu = mssm.get_vu();
   const double vd = mssm.get_vd();
   const double tanBeta = vu / vd;
   const double vev = Sqrt(Sqr(vu) + Sqr(vd));

   softsusy.setSusyMu(mssm.get_Mu());
   softsusy.setTanb(tanBeta);
   softsusy.setHvev(vev);

   for (int i = 1; i <= 3; i++) {
      for (int k = 1; k <= 3; k++) {
         softsusy.setYukawaElement(YU, i, k, mssm.get_Yu()(i-1,k-1));
         softsusy.setYukawaElement(YD, i, k, mssm.get_Yd()(i-1,k-1));
         softsusy.setYukawaElement(YE, i, k, mssm.get_Ye()(i-1,k-1));
      }
   }

   // copy soft parameters
   softsusy.setGauginoMass(1, mssm.get_MassB());
   softsusy.setGauginoMass(2, mssm.get_MassWB());
   softsusy.setGauginoMass(3, mssm.get_MassG());

   softsusy.setM3Squared(mssm.get_BMu());
   softsusy.setMh1Squared(mssm.get_mHd2());
   softsusy.setMh2Squared(mssm.get_mHu2());

   for (int i = 1; i <= 3; i++) {
      for (int k = 1; k <= 3; k++) {
         softsusy.setSoftMassElement(mQl, i, k,  mssm.get_mq2()(i-1,k-1));
         softsusy.setSoftMassElement(mUr, i, k,  mssm.get_mu2()(i-1,k-1));
         softsusy.setSoftMassElement(mDr, i, k,  mssm.get_md2()(i-1,k-1));
         softsusy.setSoftMassElement(mLl, i, k,  mssm.get_ml2()(i-1,k-1));
         softsusy.setSoftMassElement(mEr, i, k,  mssm.get_me2()(i-1,k-1));

         softsusy.setTrilinearElement(UA, i, k, mssm.get_TYu()(i-1,k-1));
         softsusy.setTrilinearElement(DA, i, k, mssm.get_TYd()(i-1,k-1));
         softsusy.setTrilinearElement(EA, i, k, mssm.get_TYe()(i-1,k-1));
      }
   }

   softsusy.setMw(softsusy.displayMwRun());
}

#endif
