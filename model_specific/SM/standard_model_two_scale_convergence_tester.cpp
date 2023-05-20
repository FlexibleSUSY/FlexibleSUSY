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

// File generated at Fri 28 Oct 2016 12:11:14

#include "standard_model_two_scale_convergence_tester.hpp"
#include <cmath>
#include <algorithm>
#include "wrappers.hpp"

namespace flexiblesusy {
namespace standard_model {

Standard_model_convergence_tester<Two_scale>::Standard_model_convergence_tester(
   StandardModel<Two_scale>* model, double accuracy_goal, const Scale_getter& sg)
   : Convergence_tester_DRbar<StandardModel<Two_scale> >(model, accuracy_goal, sg)
{
}

double Standard_model_convergence_tester<Two_scale>::max_rel_diff() const
{
   const StandardModel<Two_scale>& ol = get_last_iteration_model();
   const StandardModel<Two_scale>& ne = get_current_iteration_model();

   return ne.max_rel_diff(ol);
}

} // namespace standard_model
} // namespace flexiblesusy
