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

#ifndef PROBLEMS_H
#define PROBLEMS_H

#include <iosfwd>
#include <map>
#include <string>
#include <vector>

namespace flexiblesusy {

class Names;

/**
 * @class Problems
 * @brief stores problem flags for the spectrum generator
 */
class Problems {
public:
   Problems(const std::string& model_name_, const Names* particle_names_, const Names* parameter_names_);

   void flag_bad_mass(int particle, bool flag = true);
   void flag_running_tachyon(int particle, bool flag = true);
   void flag_pole_tachyon(int particle, bool flag = true);
   void flag_thrown(const std::string& msg = "unknown");
   void flag_no_ewsb();
   void flag_no_ewsb_tree_level();
   void flag_no_perturbative();
   void flag_no_pole_mass_convergence(int particle);
   void flag_non_perturbative_parameter(int parameter, double value, double scale, double threshold = 0.);
   void flag_no_sinThetaW_convergence();

   void unflag_bad_mass(int particle);
   void unflag_all_bad_masses();
   void unflag_running_tachyon(int particle);
   void unflag_pole_tachyon(int particle);
   void unflag_all_tachyons();
   void unflag_thrown();
   void unflag_no_ewsb();
   void unflag_no_ewsb_tree_level();
   void unflag_no_perturbative();
   void unflag_no_pole_mass_convergence(int particle);
   void unflag_non_perturbative_parameter(int parameter);
   void unflag_all_non_perturbative_parameters();
   void unflag_no_sinThetaW_convergence();

   bool is_bad_mass(int particle) const;
   bool is_running_tachyon(int particle) const;
   bool is_pole_tachyon(int particle) const;
   bool have_bad_mass() const;
   bool have_running_tachyon() const;
   bool have_pole_tachyon() const;
   bool have_tachyon() const;
   bool have_thrown() const;
   bool have_non_perturbative_parameter() const;
   bool have_failed_pole_mass_convergence() const;
   bool no_ewsb() const;
   bool no_ewsb_tree_level() const;
   bool no_perturbative() const;
   bool no_sinThetaW_convergence() const;

   void add(const Problems&);         ///< add problems from other class
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
   const std::string& get_model_name() const;

   std::vector<int> get_bad_masses() const;
   std::vector<int> get_running_tachyons() const;
   std::vector<int> get_pole_tachyons() const;
   std::vector<int> get_failed_pole_mass_convergence() const;

private:
   struct NonPerturbativeValue {
      NonPerturbativeValue() = default;
      NonPerturbativeValue(double value_, double scale_, double threshold_)
         : value(value_), scale(scale_), threshold(threshold_) {}
      double value{0.}, scale{0.}, threshold{0.};
   };

   std::string model_name;             ///< model name
   const Names* particle_names;        ///< access to particle names
   const Names* parameter_names;       ///< access to parameter names
   std::vector<int> bad_masses;        ///< imprecise mass eigenvalues
   std::vector<int> running_tachyons;  ///< tachyonic particles (running mass)
   std::vector<int> pole_tachyons;     ///< tachyonic particles (pole mass)
   std::vector<int> failed_pole_mass_convergence; ///< no convergence during pole mass calculation
   std::map<int, NonPerturbativeValue> non_pert_pars; ///< non-perturbative parmeters
   std::string exception_msg;          ///< exception message
   bool failed_ewsb{false};            ///< no EWSB
   bool failed_ewsb_tree_level{false}; ///< no tree-level EWSB
   bool non_perturbative{false};       ///< non-perturbative running
   bool failed_sinThetaW_convergence{false}; ///< sinThetaW-parameter not converged

   std::string get_parameter_name(int) const; ///< returns parameter name
};

std::ostream& operator<<(std::ostream&, const Problems&);

} // namespace flexiblesusy

#endif
