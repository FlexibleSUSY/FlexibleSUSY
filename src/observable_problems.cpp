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

#include "observable_problems.hpp"

namespace flexiblesusy {

namespace observable_problems {

void Problem_a_muon::clear()
{
   *this = Problem_a_muon();
}

bool Problem_a_muon::have_problem() const
{
   return non_perturbative_running;
}

void Problem_a_muon::flag_non_perturbative_running(double scale)
{
   non_perturbative_running = true;
   non_perturbative_running_to_scale = scale;
}

} // namespace observable_problems

void Observable_problems::clear()
{
   a_muon.clear();
}

bool Observable_problems::have_problem() const
{
   return a_muon.have_problem();
}

} // namespace flexiblesusy
