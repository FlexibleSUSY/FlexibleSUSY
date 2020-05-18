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

#ifndef OBSERVABLE_PROBLEMS_FORMAT_MATHLINK_H
#define OBSERVABLE_PROBLEMS_FORMAT_MATHLINK_H

#include "mathlink_utils.hpp"
#include "observables.hpp"
#include "observable_problems.hpp"

#include <string>

namespace flexiblesusy {

namespace observable_problems {

void mathlink_format_problems(MLINK link, const Problem_general& problems)
{
   MLPutRule(link, "general");
   MLPutFunction(link, "List", problems.number_of_problems());

   if (problems.have_non_perturbative_running())
      MLPutRuleTo(link, "True", "NonPerturbative");
   if (problems.have_thrown())
      MLPutRuleTo(link, "True", "Exceptions");
}



void mathlink_format_problems(MLINK link, const Problem_a_muon& problems)
{
   MLPutRule(link, observables::observable_names[observables::a_muon]);
   MLPutFunction(link, "List", problems.number_of_problems());

   if (problems.have_non_perturbative_running())
      MLPutRuleTo(link, "True", "NonPerturbative");
}

} // namespace observable_problems


void mathlink_format_problems(MLINK link, const Observable_problems& op)
{
   if (op.have_problem()) {
      MLPutFunction(link, "List", 2);
      mathlink_format_problems(link, op.general);
      mathlink_format_problems(link, op.a_muon);
   } else {
      MLPutFunction(link, "List", 0);
   }
}

} // namespace flexiblesusy

#endif
