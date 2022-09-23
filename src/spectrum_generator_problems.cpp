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

#include "spectrum_generator_problems.hpp"
#include "string_utils.hpp"

#include <algorithm>
#include <iostream>

namespace flexiblesusy {

Spectrum_generator_problems::Spectrum_generator_problems(const std::vector<Problems>& problems_)
   : problems(problems_)
{
}

Spectrum_generator_problems::Spectrum_generator_problems(std::vector<Problems>&& problems_)
   : problems(std::move(problems_))
{
}

Spectrum_generator_problems::Spectrum_generator_problems(const std::vector<BVP_solver_problems>& solver_problems_)
   : solver_problems(solver_problems_)
{
}

Spectrum_generator_problems::Spectrum_generator_problems(std::vector<BVP_solver_problems>&& solver_problems_)
   : solver_problems(std::move(solver_problems_))
{
}

void Spectrum_generator_problems::set_model_problems(const std::vector<Problems>& p)
{
   problems = p;
}

void Spectrum_generator_problems::set_model_problems(std::vector<Problems>&& p)
{
   problems = std::move(p);
}

void Spectrum_generator_problems::set_bvp_solver_problems(const std::vector<BVP_solver_problems>& p)
{
   solver_problems = p;
}

void Spectrum_generator_problems::set_bvp_solver_problems(std::vector<BVP_solver_problems>&& p)
{
   solver_problems = std::move(p);
}

void Spectrum_generator_problems::clear()
{
   for (auto& p: problems)
      p.clear();
   for (auto& p: solver_problems)
      p.clear();
}

bool Spectrum_generator_problems::have_problem() const
{
   return
      std::any_of(problems.cbegin(), problems.cend(),
                  [] (const Problems& p) { return p.have_problem(); }) ||
      std::any_of(solver_problems.cbegin(), solver_problems.cend(),
                  [] (const BVP_solver_problems& p) { return p.have_problem(); });
}

bool Spectrum_generator_problems::have_warning() const
{
   return
      std::any_of(problems.cbegin(), problems.cend(),
                  [] (const Problems& p) { return p.have_warning(); }) ||
      std::any_of(solver_problems.cbegin(), solver_problems.cend(),
                  [] (const BVP_solver_problems& p) { return p.have_warning(); });
}

std::vector<std::string> Spectrum_generator_problems::get_problem_strings() const
{
   std::vector<std::string> result;

   for (const auto& p: problems) {
      const auto strings = p.get_problem_strings();
      result.insert(result.end(), strings.cbegin(), strings.cend());
   }

   for (const auto& p: solver_problems) {
      const auto strings = p.get_problem_strings();
      result.insert(result.end(), strings.cbegin(), strings.cend());
   }

   return result;
}

std::vector<std::string> Spectrum_generator_problems::get_warning_strings() const
{
   std::vector<std::string> result;

   for (const auto& p: problems) {
      const auto strings = p.get_warning_strings();
      result.insert(result.end(), strings.cbegin(), strings.cend());
   }

   for (const auto& p: solver_problems) {
      const auto strings = p.get_warning_strings();
      result.insert(result.end(), strings.cbegin(), strings.cend());
   }

   return result;
}

std::string Spectrum_generator_problems::get_problem_string(const std::string& sep) const
{
   return concat(get_problem_strings(), sep);
}

std::string Spectrum_generator_problems::get_warning_string(const std::string& sep) const
{
   return concat(get_warning_strings(), sep);
}

void Spectrum_generator_problems::print_problems() const
{
   print_problems(std::cerr);
}

void Spectrum_generator_problems::print_problems(std::ostream& ostr) const
{
   if (!have_problem())
      return;

   ostr << get_problem_string();
}

void Spectrum_generator_problems::print_warnings() const
{
   print_warnings(std::cerr);
}

void Spectrum_generator_problems::print_warnings(std::ostream& ostr) const
{
   if (!have_warning())
      return;

   ostr << get_warning_string();
}

const std::vector<Problems>& Spectrum_generator_problems::get_model_problems() const
{
   return problems;
}

std::vector<Problems>& Spectrum_generator_problems::get_model_problems()
{
   return problems;
}

const std::vector<BVP_solver_problems>& Spectrum_generator_problems::get_bvp_solver_problems() const
{
   return solver_problems;
}

std::vector<BVP_solver_problems>& Spectrum_generator_problems::get_bvp_solver_problems()
{
   return solver_problems;
}

void Spectrum_generator_problems::flag_no_convergence()
{
   for (auto& p: solver_problems)
      p.flag_no_convergence();
}

void Spectrum_generator_problems::unflag_no_convergence()
{
   for (auto& p: solver_problems)
      p.unflag_no_convergence();
}

bool Spectrum_generator_problems::no_convergence() const
{
   return std::any_of(solver_problems.cbegin(), solver_problems.cend(),
                      [] (const BVP_solver_problems& p) { return p.no_convergence(); });
}

std::ostream& operator<<(std::ostream& ostr, const Spectrum_generator_problems& problems)
{
   problems.print_problems(ostr);
   problems.print_warnings(ostr);
   return ostr;
}

} // namespace flexiblesusy
