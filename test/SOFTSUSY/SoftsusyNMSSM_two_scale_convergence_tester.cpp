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

#include "SoftsusyNMSSM_two_scale_convergence_tester.hpp"

namespace flexiblesusy {

SoftsusyNMSSM_convergence_tester::SoftsusyNMSSM_convergence_tester(SoftsusyNMSSM<Two_scale>* snmssm_, double accuracy_goal_)
   : Convergence_tester_DRbar<SoftsusyNMSSM<Two_scale> >(snmssm_, accuracy_goal_)
{
}

SoftsusyNMSSM_convergence_tester::~SoftsusyNMSSM_convergence_tester()
{
}

double SoftsusyNMSSM_convergence_tester::max_rel_diff() const
{
   return sumTol(get_current_iteration_model(), get_last_iteration_model());
}

double SoftsusyNMSSM_convergence_tester::sumTol(const SoftsusyNMSSM<Two_scale>& in, const SoftsusyNMSSM<Two_scale>& out) const
{
   return softsusy::sumTol(in, out, 100);
}

}
