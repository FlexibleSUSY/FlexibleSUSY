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

#include "SoftsusyNMSSM_two_scale_initial_guesser.hpp"
#include "SoftsusyNMSSM_two_scale.hpp"
#include "nmssmsoftsusy.h"
#include "logger.hpp"

#include <cassert>

using namespace softsusy;

namespace flexiblesusy {

SoftsusyNMSSM_initial_guesser::SoftsusyNMSSM_initial_guesser(SoftsusyNMSSM<Two_scale>* nmssm_,
                                               const SoftsusyNMSSM_parameter_point& pp_,
                                               const SoftsusyNMSSM_low_scale_constraint&,
                                               const SoftsusyNMSSM_susy_scale_constraint&,
                                               const SoftsusyNMSSM_sugra_constraint&)
   : Initial_guesser()
   , nmssm(nmssm_)
   , oneset()
   , pp(pp_)
   , ewsbBCscale(false)
{
   assert(nmssm && "SoftsusyNMSSM_initial_guesser: Error: pointer to SoftsusyNMSSM"
          " cannot be zero");

   const double alphasMZ = 0.1187, mtop = 173.4, mbmb = 4.2;
   oneset.setAlpha(legacy::ALPHAS, alphasMZ);
   oneset.setPoleMt(mtop);
   oneset.setMass(legacy::mBottom, mbmb);
   oneset.toMz();
}

SoftsusyNMSSM_initial_guesser::~SoftsusyNMSSM_initial_guesser()
{
}

void SoftsusyNMSSM_initial_guesser::guess()
{
   const static NmssmSoftsusy empty;

   const double mx = pp.mxGuess;
   const double muFirst = nmssm->displaySusyMu(); /// Remember initial values
   const bool setTbAtMXflag = nmssm->displaySetTbAtMX();

   nmssm->setSoftsusy(empty); /// Always starts from an empty object
   nmssm->setSetTbAtMX(setTbAtMXflag);
   nmssm->setData(oneset);
   nmssm->setMw(MW);

   double mz = nmssm->displayMz();

   if (oneset.displayMu() != mz) {
      WARNING("NmssmSoftsusy::lowOrg called with oneset at scale\n"
              << oneset.displayMu() << "instead of " << mz);
   }

   NmssmSusy t(nmssm->guessAtSusyMt(pp.tanBeta, pp.get_nmpars(), oneset));

   t.setLoops(2); /// 2 loops should protect against ht Landau pole
   int err = t.runto(mx);
   if (err)
      ERROR("<SoftsusyNMSSM_initial_guesser::guess()>: non-perturbative running to"
            " mx = " << mx);

   nmssm->setSusy(t);

   const int sgnMu = pp.signMu;
   nmssm->standardSugra(pp.m0, pp.m12, pp.a0);
   if (softsusy::GUTlambda)
      nmssm->setLambda(pp.lambda);

   if ((sgnMu == 1 || sgnMu == -1) && !ewsbBCscale) {
      if (Z3) {
         nmssm->setSusyMu(0.0);
         nmssm->setM3Squared(0.0);
      } else {
         nmssm->setSusyMu(sgnMu * MZ);
         nmssm->setM3Squared(1.0e6);
      }
   } else {
      nmssm->setSusyMu(muFirst);
      nmssm->setM3Squared(muFirst);
   }

   err = nmssm->run(mx, mz);
   if (err)
      ERROR("<SoftsusyNMSSM_initial_guesser::guess()>: non-perturbative running from"
            " mx = " << mx << " to mz = " << mz);

   if (sgnMu == 1 || sgnMu == -1)
      nmssm->rewsbTreeLevel(sgnMu);

   nmssm->physical(0);
   nmssm->setThresholds(3);
   nmssm->setLoops(2);

   // test for problems
   sProblem problem = nmssm->displayProblem();
   if (problem.test()) {
      if (problem.testSeriousProblem())
         WARNING("<SoftsusyNMSSM_initial_guesser::guess()>: serious problem: " << problem);
      else
         WARNING("<SoftsusyNMSSM_initial_guesser::guess()>: problem: " << problem);
   }
}

} // namespace flexiblesusy
