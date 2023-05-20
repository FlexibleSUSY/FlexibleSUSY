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

#ifndef MSSM_SOLVER_H
#define MSSM_SOLVER_H

#include "rg_flow.hpp"
#include "softsusy.h"

namespace flexiblesusy {

class SoftSusy_t;

template<>
class RGFlow<SoftSusy_t> {
public:
   typedef void (*THighScaleBoundaryCondition)(softsusy::MssmSoftsusy&, const softsusy::DoubleVector&);

   RGFlow();
   ~RGFlow();

   void solve();

   const softsusy::sPhysical& displayPhys() const;
   const softsusy::sProblem& displayProblem() const;

   void setHighScaleBoundaryCondition(THighScaleBoundaryCondition);
   void setMxGuess(double);
   void setSoftHighScalePars(const softsusy::DoubleVector&);
   void setSignMu(int);
   void setTanBeta(double);
   void setLowScaleBoundaryConditions(const softsusy::QedQcd_legacy&);
   void setGaugeUnification(bool);

private:
   softsusy::MssmSoftsusy       fMssmSoftSusy;
   THighScaleBoundaryCondition  fHighScaleBoundaryCondition;
   double                       fMxGuess;
   softsusy::DoubleVector       fSoftHighScalePars;
   int                          fSignMu;
   double                       fTanBeta;
   softsusy::QedQcd_legacy      fLowScaleBoundaryContitions;
   bool                         fGaugeUnification;
};

RGFlow<SoftSusy_t>::RGFlow()
   : fMssmSoftSusy()
   , fHighScaleBoundaryCondition(softsusy::sugraBcs)
   , fMxGuess(1.9e16)
   , fSoftHighScalePars(3)
   , fSignMu(1)
   , fTanBeta(10)
   , fLowScaleBoundaryContitions()
   , fGaugeUnification(true)
{
}

RGFlow<SoftSusy_t>::~RGFlow()
{
}

void RGFlow<SoftSusy_t>::solve()
{
   fMssmSoftSusy.lowOrg(fHighScaleBoundaryCondition, fMxGuess, fSoftHighScalePars,
                        fSignMu, fTanBeta, fLowScaleBoundaryContitions,
                        fGaugeUnification);
}

const softsusy::sPhysical& RGFlow<SoftSusy_t>::displayPhys() const
{
   return fMssmSoftSusy.displayPhys();
}

const softsusy::sProblem& RGFlow<SoftSusy_t>::displayProblem() const
{
   return fMssmSoftSusy.displayProblem();
}

void RGFlow<SoftSusy_t>::setHighScaleBoundaryCondition(THighScaleBoundaryCondition bc)
{
   fHighScaleBoundaryCondition = bc;
}

void RGFlow<SoftSusy_t>::setMxGuess(double mxGuess)
{
   fMxGuess = mxGuess;
}

void RGFlow<SoftSusy_t>::setSoftHighScalePars(const softsusy::DoubleVector& softHighScalePars)
{
   fSoftHighScalePars = softHighScalePars;
}

void RGFlow<SoftSusy_t>::setSignMu(int signMu)
{
   fSignMu = signMu;
}

void RGFlow<SoftSusy_t>::setTanBeta(double tanBeta)
{
   fTanBeta = tanBeta;
}

void RGFlow<SoftSusy_t>::setLowScaleBoundaryConditions(const softsusy::QedQcd_legacy& lowScaleBC)
{
   fLowScaleBoundaryContitions = lowScaleBC;
}

void RGFlow<SoftSusy_t>::setGaugeUnification(bool gaugeUnification)
{
   fGaugeUnification = gaugeUnification;
}

}

#endif
