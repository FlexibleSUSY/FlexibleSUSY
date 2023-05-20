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

#include "SoftsusyMSSM_two_scale_convergence_tester.hpp"

namespace flexiblesusy {

SoftsusyMSSM_convergence_tester::SoftsusyMSSM_convergence_tester(SoftsusyMSSM<Two_scale>* mssm_, double accuracy_goal_)
   : Convergence_tester_DRbar<SoftsusyMSSM<Two_scale> >(mssm_, accuracy_goal_)
{
}

SoftsusyMSSM_convergence_tester::~SoftsusyMSSM_convergence_tester()
{
}

double SoftsusyMSSM_convergence_tester::max_rel_diff() const
{
   return sumTol(get_current_iteration_model(), get_last_iteration_model());
}

double SoftsusyMSSM_convergence_tester::sumTol(const SoftsusyMSSM<Two_scale>& in, const SoftsusyMSSM<Two_scale>& out) const
{
  softsusy::drBarPars inforLoops(in.displayDrBarPars()),
    outforLoops(out.displayDrBarPars());
  softsusy::DoubleVector sT(32);

  ::sumTol(inforLoops, outforLoops, sT);

  return sT.max();
}

}
