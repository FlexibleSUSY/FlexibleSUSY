// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#include "SoftsusyMSSM_two_scale_initial_guesser.hpp"
#include "SoftsusyMSSM_two_scale.hpp"
#include "softsusy.h"
#include "logger.hpp"

#include <cassert>

using namespace softsusy;

namespace flexiblesusy {

SoftsusyMSSM_initial_guesser::SoftsusyMSSM_initial_guesser(SoftsusyMSSM<Two_scale>* mssm_,
                                           const SoftsusyMSSM_parameter_point& pp_,
                                           const SoftsusyMSSM_low_scale_constraint&,
                                           const SoftsusyMSSM_susy_scale_constraint&,
                                           const SoftsusyMSSM_sugra_constraint&)
   : Initial_guesser()
   , mssm(mssm_)
   , oneset()
   , pp(pp_)
   , ewsbBCscale(false)
{
   assert(mssm && "SoftsusyMSSM_initial_guesser: Error: pointer to Mssm"
          " cannot be zero");

   const double alphasMZ = 0.1187, mtop = 173.4, mbmb = 4.2;
   oneset.setAlpha(legacy::ALPHAS, alphasMZ);
   oneset.setPoleMt(mtop);
   oneset.setMass(legacy::mBottom, mbmb);
   oneset.toMz();
}

SoftsusyMSSM_initial_guesser::~SoftsusyMSSM_initial_guesser()
{
}

void SoftsusyMSSM_initial_guesser::guess()
{
   const static MssmSoftsusy empty;

   const double mx = pp.mxGuess;
   const double muFirst = mssm->displaySusyMu(); /// Remember initial values
   const bool setTbAtMXflag = mssm->displaySetTbAtMX();
   const bool altFlag = mssm->displayAltEwsb();
   const double m32 = mssm->displayGravitino();
   const double muCondFirst = mssm->displayMuCond();
   const double maCondFirst = mssm->displayMaCond();

   mssm->setSoftsusy(empty); /// Always starts from an empty object
   mssm->setSetTbAtMX(setTbAtMXflag);
   if (altFlag)
      mssm->useAlternativeEwsb();
   mssm->setData(oneset);
   mssm->setMw(MW);
   mssm->setM32(m32);
   mssm->setMuCond(muCondFirst);
   mssm->setMaCond(maCondFirst);

   double mz = mssm->displayMz();

   if (oneset.displayMu() != mz) {
      WARNING("lowOrg in softsusy.cpp called with oneset at scale\n"
              << oneset.displayMu() << "instead of " << mz);
   }

   MssmSusy t(mssm->guessAtSusyMt(pp.tanBeta, oneset));

   t.setLoops(2); /// 2 loops should protect against ht Landau pole
   int err = t.runto(mx);
   if (err)
      ERROR("<SoftsusyMSSM_initial_guesser::guess()>: non-perturbative running to"
            " mx = " << mx);

   mssm->setSusy(t);

   /// Initial guess: B=0, mu=1st parameter, need better guesses
   const int sgnMu = pp.signMu;
   mssm->standardSugra(pp.m0, pp.m12, pp.a0);

   if ((sgnMu == 1 || sgnMu == -1) && !ewsbBCscale) {
      mssm->setSusyMu(sgnMu * 1.0);
      mssm->setM3Squared(0.);
   } else {
      if (mssm->displayAltEwsb()) {
         mssm->setSusyMu(mssm->displayMuCond());
         mssm->setM3Squared(mssm->displayMaCond());
      } else {
         mssm->setSusyMu(muFirst);
         mssm->setM3Squared(muFirst);
      }
   }

   err = mssm->run(mx, mz);
   if (err)
      ERROR("<SoftsusyMSSM_initial_guesser::guess()>: non-perturbative running from"
            " mx = " << mx << " to mz = " << mz);

   if (sgnMu == 1 || sgnMu == -1)
      mssm->rewsbTreeLevel(sgnMu);

   mssm->physical(0);
   mssm->setThresholds(3);
   mssm->setLoops(2);

   // test for problems
   sProblem problem = mssm->displayProblem();
   if (problem.test()) {
      if (problem.testSeriousProblem())
         WARNING("<SoftsusyMSSM_initial_guesser::guess()>: serious problem: " << problem);
      else
         WARNING("<SoftsusyMSSM_initial_guesser::guess()>: problem: " << problem);
   }
}

}
