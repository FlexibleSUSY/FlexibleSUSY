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

#ifndef SPECTRUM_GENERATOR_PROBLEMS_H
#define SPECTRUM_GENERATOR_PROBLEMS_H

#include "problems.hpp"
#include "bvp_solver_problems.hpp"

namespace flexiblesusy {

class Spectrum_generator_problems {
public:
   Spectrum_generator_problems() = default;
   explicit Spectrum_generator_problems(const std::vector<Problems>&);
   explicit Spectrum_generator_problems(std::vector<Problems>&&);
   explicit Spectrum_generator_problems(const std::vector<BVP_solver_problems>&);
   explicit Spectrum_generator_problems(std::vector<BVP_solver_problems>&&);

   void clear();                      ///< clear all problems
   bool have_problem() const;         ///< problems which yield invalid spectrum
   bool have_warning() const;         ///< warnings
   std::vector<std::string> get_problem_strings() const;
   std::vector<std::string> get_warning_strings() const;
   std::string get_problem_string(const std::string& sep = "\n") const;
   std::string get_warning_string(const std::string& sep = "\n") const;
   void print_problems() const;
   void print_problems(std::ostream&) const;
   void print_warnings() const;
   void print_warnings(std::ostream&) const;

   void set_model_problems(const std::vector<Problems>&); ///< (re)set all model problems
   void set_model_problems(std::vector<Problems>&&);      ///< (re)set all model problems
   const std::vector<Problems>& get_model_problems() const;
   std::vector<Problems>& get_model_problems();
   int get_number_of_models() const { return problems.size(); }

   void set_bvp_solver_problems(const std::vector<BVP_solver_problems>&); ///< (re)set all BVP solver problems
   void set_bvp_solver_problems(std::vector<BVP_solver_problems>&&);      ///< (re)set all BVP solver problems
   const std::vector<BVP_solver_problems>& get_bvp_solver_problems() const;
   std::vector<BVP_solver_problems>& get_bvp_solver_problems();
   int get_number_of_bvp_solvers() const { return solver_problems.size(); }

   void flag_no_convergence();
   void unflag_no_convergence();
   bool no_convergence() const;

private:
   std::vector<Problems> problems; ///< model problems
   std::vector<BVP_solver_problems> solver_problems; ///< BVP solver problems
};

std::ostream& operator<<(std::ostream&, const Spectrum_generator_problems&);

} // namespace flexiblesusy

#endif
