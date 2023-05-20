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

#include "composite_convergence_tester.hpp"
#include <algorithm>

namespace flexiblesusy {

/**
 * Calls the Convergence_tester::accuracy_goal_reached()
 * functions of all convergence testers.  If at least one of them
 * returns false, then false is returned.  Othewise (all convergence
 * testers yield true) true is returned.
 *
 * @note An alternative implementation would be: When the first
 * convergence tester returns false, then false is returned.
 * Otherwise (all convergence testers yield true) true is returned.
 * This implementation is currently not realized, because it is not
 * ensured that all
 * Convergence_tester::accuracy_goal_reached() functions
 * are called.
 *
 * @return true if and only if all accuracy_goal_reached() function
 * calls of all convergence testers return true.  Otherewise false is
 * returned.
 */
bool Composite_convergence_tester::accuracy_goal_reached()
{
   bool precision_reached = true;

   for (auto ct: testers) {
      const bool tester_result = ct->accuracy_goal_reached();
      precision_reached = precision_reached && tester_result;
   }

   return precision_reached;
}

int Composite_convergence_tester::max_iterations() const
{
   if (testers.empty())
      return 0;

   return (*std::max_element(testers.begin(), testers.end(),
                             [](const Convergence_tester* a,
                                const Convergence_tester* b) {
                                return a->max_iterations() < b->max_iterations();
                             }))->max_iterations();
}

void Composite_convergence_tester::restart()
{
   for (auto ct: testers)
      ct->restart();
}

void Composite_convergence_tester::add_convergence_tester(Convergence_tester* t)
{
   if (t)
      testers.push_back(t);
}

} // namespace flexiblesusy
