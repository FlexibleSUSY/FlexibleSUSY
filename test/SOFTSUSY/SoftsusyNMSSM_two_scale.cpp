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

#include "SoftsusyNMSSM_two_scale.hpp"

namespace flexiblesusy {

SoftsusyNMSSM<Two_scale>::SoftsusyNMSSM()
   : Model()
   , precision(1.0e-5)
{
}

SoftsusyNMSSM<Two_scale>::SoftsusyNMSSM(const softsusy::SoftParsNmssm& softPars)
{
   setSoftPars(softPars);
}

SoftsusyNMSSM<Two_scale>::~SoftsusyNMSSM()
{
}

void SoftsusyNMSSM<Two_scale>::calculate_spectrum()
{
   run_to(maximum(displayMsusy(), softsusy::MZ));
   physical(3);
   runto(displayMz());
}

void SoftsusyNMSSM<Two_scale>::run_to(double scale, double eps)
{
   if (eps < 0.0)
      eps = precision;
   runto(scale, eps);
}

}
