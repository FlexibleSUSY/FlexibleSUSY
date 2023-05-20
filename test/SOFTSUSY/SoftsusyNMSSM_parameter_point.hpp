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

#ifndef SoftsusyNMSSM_PARAMETER_POINT_H
#define SoftsusyNMSSM_PARAMETER_POINT_H

#include "linalg.h"

namespace flexiblesusy {

struct SoftsusyNMSSM_parameter_point {
   SoftsusyNMSSM_parameter_point()
      : m0(125.0)
      , m12(500.0)
      , a0(0.0)
      , mxGuess(1.9e16)
      , msGuess(1000.)
      , signMu(1)
      , tanBeta(10.0)
      , lambda(0.1)
      , kappa(0.)
      , svev(1000.)
      , xiF(0.)
      , muPrime(0.)
   {}
   softsusy::DoubleVector get_soft_pars() const {
      softsusy::DoubleVector highScaleSoftPars(3);
      highScaleSoftPars(1) = m0;
      highScaleSoftPars(2) = m12;
      highScaleSoftPars(3) = a0;
      return highScaleSoftPars;
   }
   softsusy::DoubleVector get_nmpars() const {
      softsusy::DoubleVector nmpars(5);
      nmpars(1) = lambda;
      nmpars(2) = kappa;
      nmpars(3) = svev;
      nmpars(4) = xiF;
      nmpars(5) = muPrime;
      return nmpars;
   }
   friend std::ostream& operator<<(std::ostream& os, const SoftsusyNMSSM_parameter_point& pp) {
      os << "CNMSSM parameter point:"
         << " m0=" << pp.m0
         << ", m12=" << pp.m12
         << ", a0=" << pp.a0
         << ", mxGuess=" << pp.mxGuess
         << ", signMu=" << pp.signMu
         << ", tanBeta=" << pp.tanBeta
         << ", lambda=" << pp.lambda
         << ", kappa=" << pp.kappa
         << ", svev=" << pp.svev
         << ", xiF=" << pp.xiF
         << ", muPrime=" << pp.muPrime
         << '\n';
      return os;
   }

   double m0, m12, a0, mxGuess, msGuess, signMu, tanBeta;
   double lambda, kappa, svev, xiF, muPrime;
};

}

#endif
