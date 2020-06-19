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

#include "bvp_solver_problems.hpp"
#include "string_utils.hpp"

#include <algorithm>
#include <iostream>

namespace flexiblesusy {

BVP_solver_problems::BVP_solver_problems(const std::string& name_)
   : name(name_)
{
}

void BVP_solver_problems::clear()
{
   failed_convergence = false;
}

bool BVP_solver_problems::have_problem() const
{
   return no_convergence();
}

bool BVP_solver_problems::have_warning() const
{
   return false; // no warnings yet
}

std::vector<std::string> BVP_solver_problems::get_problem_strings() const
{
   std::vector<std::string> result;

   if (no_convergence())
      result.emplace_back(name + " no convergence");

   result.shrink_to_fit();

   return result;
}

std::vector<std::string> BVP_solver_problems::get_warning_strings() const
{
   return {};
}

std::string BVP_solver_problems::get_problem_string() const
{
   return concat(get_problem_strings(), '\n');
}

std::string BVP_solver_problems::get_warning_string() const
{
   return concat(get_warning_strings(), '\n');
}

void BVP_solver_problems::print_problems() const
{
   print_problems(std::cerr);
}

void BVP_solver_problems::print_problems(std::ostream& ostr) const
{
   if (!have_problem())
      return;

   ostr << get_problem_string();
}

void BVP_solver_problems::print_warnings() const
{
   print_warnings(std::cerr);
}

void BVP_solver_problems::print_warnings(std::ostream& ostr) const
{
   if (!have_warning())
      return;

   ostr << get_warning_string();
}

const std::string& BVP_solver_problems::get_solver_name() const
{
   return name;
}

void BVP_solver_problems::flag_no_convergence()
{
   failed_convergence = true;
}

void BVP_solver_problems::unflag_no_convergence()
{
   failed_convergence = false;
}

bool BVP_solver_problems::no_convergence() const
{
   return failed_convergence;
}

std::ostream& operator<<(std::ostream& ostr, const BVP_solver_problems& problems)
{
   problems.print_problems(ostr);
   problems.print_warnings(ostr);
   return ostr;
}

} // namespace flexiblesusy
