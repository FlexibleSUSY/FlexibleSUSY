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

#ifndef BVP_SOLVER_PROBLEMS_H
#define BVP_SOLVER_PROBLEMS_H

#include "problems.hpp"

namespace flexiblesusy {

class BVP_solver_problems {
public:
   BVP_solver_problems(const std::string&);

   void clear();                      ///< clear all problems
   bool have_problem() const;         ///< problems which yield invalid spectrum
   bool have_warning() const;         ///< warnings
   std::vector<std::string> get_problem_strings() const;
   std::vector<std::string> get_warning_strings() const;
   std::string get_problem_string() const;
   std::string get_warning_string() const;
   void print_problems() const;
   void print_problems(std::ostream&) const;
   void print_warnings() const;
   void print_warnings(std::ostream&) const;
   const std::string& get_solver_name() const;

   void flag_no_convergence();
   void unflag_no_convergence();
   bool no_convergence() const;

private:
   std::string name;               ///< BVP solver name
   bool failed_convergence{false}; ///< no convergence
};

std::ostream& operator<<(std::ostream&, const BVP_solver_problems&);

} // namespace flexiblesusy

#endif
