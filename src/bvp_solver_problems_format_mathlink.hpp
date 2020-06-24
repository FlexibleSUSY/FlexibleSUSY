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

#ifndef BVP_SOLVER_PROBLEMS_FORMAT_MATHLINK_H
#define BVP_SOLVER_PROBLEMS_FORMAT_MATHLINK_H

#include "mathlink_utils.hpp"
#include "bvp_solver_problems.hpp"

namespace flexiblesusy {

/// format BVP solver problems to MathLink output
inline void mathlink_format_problems(MLINK link, const BVP_solver_problems& sp)
{
   MLPutFunction(link, "List", sp.number_of_problems());

   if (sp.no_convergence()) {
      MLPutRuleTo(link, "True", "NoConvergence");
   }
}

/// format BVP solver warnings to MathLink output
inline void mathlink_format_warnings(MLINK link, const BVP_solver_problems& sp)
{
   MLPutFunction(link, "List", 0);
}

} // namespace flexiblesusy

#endif
