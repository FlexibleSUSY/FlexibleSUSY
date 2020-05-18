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

// general /////////////////////////////////////////////////////////////

void Problem_general::clear()
{
   *this = Problem_general();
}

bool Problem_general::have_problem() const
{
   return size() > 0;
}

unsigned Problem_general::size() const
{
   unsigned count = 0;

   if (have_non_perturbative_running()) {
      count++;
   }

   if (have_thrown()) {
      count++;
   }

   return count;
}

void Problem_general::flag_non_perturbative_running(double scale)
{
   non_perturbative_running = true;
   non_perturbative_running_to_scale = scale;
}

bool Problem_general::have_non_perturbative_running() const
{
   return non_perturbative_running;
}

double Problem_general::get_non_perturbative_running_scale() const
{
   return non_perturbative_running_to_scale;
}

void Problem_general::flag_thrown(const char* msg)
{
   thrown = true;
   thrown_msg = msg;
}

bool Problem_general::have_thrown() const
{
   return thrown;
}

const char* Problem_general::get_thrown_message() const
{
   return thrown_msg;
}

// a_muon //////////////////////////////////////////////////////////////

void Problem_a_muon::clear()
{
   *this = Problem_a_muon();
}

bool Problem_a_muon::have_problem() const
{
   return size() > 0;
}

unsigned Problem_a_muon::size() const
{
   unsigned count = 0;

   if (have_non_perturbative_running()) {
      count++;
   }

   return count;
}

void Problem_a_muon::flag_non_perturbative_running(double scale)
{
   non_perturbative_running = true;
   non_perturbative_running_to_scale = scale;
}

bool Problem_a_muon::have_non_perturbative_running() const
{
   return non_perturbative_running;
}

double Problem_a_muon::get_non_perturbative_running_scale() const
{
   return non_perturbative_running_to_scale;
}

} // namespace observable_problems

// wrapper class ///////////////////////////////////////////////////////

void Observable_problems::clear()
{
   general.clear();
   a_muon.clear();
}

bool Observable_problems::have_problem() const
{
   return size() > 0;
}

unsigned Observable_problems::size() const
{
   unsigned count = 0;
   count += general.size();
   count += a_muon.size();
   return count;
}

} // namespace flexiblesusy
