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

#include "problems.hpp"
#include "config.h"
#include "logger.hpp"
#include "names.hpp"
#include "string_format.hpp"
#include "string_utils.hpp"

#include <algorithm>
#include <iostream>
#include <vector>

namespace flexiblesusy {
namespace {

template <class T>
void vector_or(std::vector<T>& v1, const std::vector<T>& v2)
{
   if (v1.size() != v2.size()) {
      ERROR("Cannot combine vectors of incompatible size.");
      return;
   }

   for (typename std::vector<T>::size_type i = 0; i < v2.size(); i++)
      if (v2[i]) v1[i] = v2[i];
}

template <class T1, class T2>
void map_or(std::map<T1,T2>& m1, const std::map<T1,T2>& m2)
{
   m1.insert(m2.cbegin(), m2.cend());
}

} // anonymous namespace

Problems::Problems(const std::string& model_name_,
                   const Names* particle_names_, const Names* parameter_names_)
   : model_name(model_name_)
   , particle_names(particle_names_)
   , parameter_names(parameter_names_)
   , bad_masses(particle_names_->size())
   , running_tachyons(particle_names_->size())
   , pole_tachyons(particle_names_->size())
   , failed_pole_mass_convergence(particle_names_->size())
{
}

void Problems::clear()
{
   std::fill(bad_masses.begin(), bad_masses.end(), 0);
   std::fill(running_tachyons.begin(), running_tachyons.end(), 0);
   std::fill(pole_tachyons.begin(), pole_tachyons.end(), 0);
   std::fill(failed_pole_mass_convergence.begin(), failed_pole_mass_convergence.end(), 0);
   non_pert_pars.clear();
   exception_msg = "";
   failed_ewsb = false;
   failed_ewsb_tree_level = false;
   non_perturbative = false;
   failed_sinThetaW_convergence = false;
}

void Problems::add(const Problems& other)
{
   if (model_name != other.get_model_name())
      WARNING("Adding problems from " << other.get_model_name() << " to " << model_name);

   vector_or(bad_masses, other.bad_masses);
   vector_or(running_tachyons, other.running_tachyons);
   vector_or(pole_tachyons, other.pole_tachyons);
   vector_or(failed_pole_mass_convergence, other.failed_pole_mass_convergence);
   map_or(non_pert_pars, other.non_pert_pars);

   if (exception_msg.empty() && !other.exception_msg.empty())
      exception_msg = other.exception_msg;

   failed_ewsb = failed_ewsb || other.failed_ewsb;
   failed_ewsb_tree_level = failed_ewsb_tree_level || other.failed_ewsb_tree_level;
   non_perturbative = non_perturbative || other.non_perturbative;
   failed_sinThetaW_convergence = failed_sinThetaW_convergence || other.failed_sinThetaW_convergence;
}

bool Problems::have_problem() const
{
   return have_tachyon() || failed_ewsb || failed_ewsb_tree_level
      || non_perturbative || failed_sinThetaW_convergence
      || have_thrown()
      || have_failed_pole_mass_convergence()
      || have_non_perturbative_parameter();
}

bool Problems::have_warning() const
{
   return have_bad_mass();
}

std::string Problems::get_parameter_name(int idx) const
{
   if (idx == -1)
      return "Q";

   return parameter_names->get(idx);
}

std::vector<std::string> Problems::get_problem_strings() const
{
   std::vector<std::string> strings;
   const auto n_particles = particle_names->size();

   for (int i = 0; i < n_particles; ++i) {
      if (running_tachyons[i])
         strings.emplace_back("running tachyon " + particle_names->get(i));
   }
   for (int i = 0; i < n_particles; ++i) {
      if (pole_tachyons[i])
         strings.emplace_back("pole tachyon " + particle_names->get(i));
   }
   if (failed_ewsb)
      strings.emplace_back("no ewsb");
   if (failed_ewsb_tree_level)
      strings.emplace_back("no ewsb at tree-level");
   if (non_perturbative)
      strings.emplace_back("non-perturbative");
   if (failed_sinThetaW_convergence)
      strings.emplace_back("no sinThetaW convergence");
   if (have_thrown())
      strings.emplace_back("exception thrown(" + exception_msg + ")");
   for (int i = 0; i < n_particles; ++i) {
      if (failed_pole_mass_convergence[i])
         strings.emplace_back("no M" + particle_names->get(i) + " pole convergence");
   }

   for (const auto& par: non_pert_pars) {
      const std::string par_name = get_parameter_name(par.first);
      std::string str("non-perturbative " + par_name);
      if (par.second.threshold > 0) {
         str += " [|" + par_name + "|(" +
                flexiblesusy::to_string(par.second.scale) + ") = " +
                flexiblesusy::to_string(par.second.value) +
                " > " + flexiblesusy::to_string(par.second.threshold) + "]";
      } else {
         str += " [" + par_name + "(" +
                flexiblesusy::to_string(par.second.scale) +
                ") = " + flexiblesusy::to_string(par.second.value) + "]";
      }
      strings.emplace_back(str);
   }

   strings.shrink_to_fit();

   return strings;
}

std::vector<std::string> Problems::get_warning_strings() const
{
   std::vector<std::string> strings;
   const auto n_particles = particle_names->size();

   for (int i = 0; i < n_particles; ++i) {
      if (bad_masses[i])
         strings.emplace_back("Warning: imprecise M" + particle_names->get(i));
   }

   return strings;
}

std::string Problems::get_problem_string(const std::string& sep) const
{
   return concat(get_problem_strings(), sep);
}

std::string Problems::get_warning_string(const std::string& sep) const
{
   return concat(get_warning_strings(), sep);
}

void Problems::print_problems() const
{
   print_problems(std::cerr);
}

void Problems::print_problems(std::ostream& ostr) const
{
   if (!have_problem())
      return;

   ostr << get_problem_string();
}

void Problems::print_warnings() const
{
   print_warnings(std::cerr);
}

void Problems::print_warnings(std::ostream& ostr) const
{
   if (!have_warning())
      return;

   ostr << get_warning_string();
}

const std::string& Problems::get_model_name() const
{
   return model_name;
}

void Problems::flag_bad_mass(int particle, bool flag)
{
   bad_masses.at(particle) = flag;
}

void Problems::flag_running_tachyon(int particle, bool flag)
{
   running_tachyons.at(particle) = flag;
#if defined(ENABLE_VERBOSE) || defined(ENABLE_DEBUG)
   if (flag)
      WARNING("running " << particle_names->get(particle) << " tachyon");
#endif
}

void Problems::flag_pole_tachyon(int particle, bool flag)
{
   pole_tachyons.at(particle) = flag;
#if defined(ENABLE_VERBOSE) || defined(ENABLE_DEBUG)
   if (flag)
      WARNING("pole " << particle_names->get(particle) << " tachyon");
#endif
}

void Problems::flag_thrown(const std::string& msg)
{
   exception_msg = msg;
}

void Problems::flag_no_ewsb()
{
   failed_ewsb = true;
}

void Problems::flag_no_ewsb_tree_level()
{
   failed_ewsb_tree_level = true;
}

void Problems::flag_no_perturbative()
{
   non_perturbative = true;
}

void Problems::flag_no_pole_mass_convergence(int particle)
{
   failed_pole_mass_convergence.at(particle) = true;
}

void Problems::flag_non_perturbative_parameter(
   int parameter, double value, double scale, double threshold)
{
   const auto n_parameters = parameter_names->size();

   if (parameter < -1 || parameter >= n_parameters)
      ERROR("Parameter index " << parameter << " out of range [" << -1 << ", "
            << (n_parameters - 1) << "]");

   VERBOSE_MSG("Problem: " << get_parameter_name(parameter) << "(Q = "
               << scale << ") = " << value
               << " is non-perturbative (threshold = " << threshold << ").");

   non_pert_pars[parameter] = NonPerturbativeValue(value, scale, threshold);
}

void Problems::flag_no_sinThetaW_convergence()
{
   failed_sinThetaW_convergence = true;
}

void Problems::unflag_bad_mass(int particle)
{
   bad_masses.at(particle) = false;
}

void Problems::unflag_all_bad_masses()
{
   std::fill(bad_masses.begin(), bad_masses.end(), 0);
}

void Problems::unflag_running_tachyon(int particle)
{
   running_tachyons.at(particle) = false;
}

void Problems::unflag_pole_tachyon(int particle)
{
   pole_tachyons.at(particle) = false;
}

void Problems::unflag_all_tachyons()
{
   std::fill(running_tachyons.begin(), running_tachyons.end(), 0);
   std::fill(pole_tachyons.begin(), pole_tachyons.end(), 0);
}

void Problems::unflag_thrown()
{
   exception_msg = "";
}

void Problems::unflag_no_ewsb()
{
   failed_ewsb = false;
}

void Problems::unflag_no_ewsb_tree_level()
{
   failed_ewsb_tree_level = false;
}

void Problems::unflag_no_perturbative()
{
   non_perturbative = false;
}

void Problems::unflag_no_pole_mass_convergence(int particle)
{
   failed_pole_mass_convergence.at(particle) = false;
}

void Problems::unflag_non_perturbative_parameter(int parameter)
{
   non_pert_pars.erase(parameter);
}

void Problems::unflag_all_non_perturbative_parameters()
{
   non_pert_pars.clear();
}

void Problems::unflag_no_sinThetaW_convergence()
{
   failed_sinThetaW_convergence = false;
}

bool Problems::is_bad_mass(int particle) const
{
   return bad_masses.at(particle);
}

bool Problems::is_running_tachyon(int particle) const
{
   return running_tachyons.at(particle);
}

bool Problems::is_pole_tachyon(int particle) const
{
   return pole_tachyons.at(particle);
}

bool Problems::have_bad_mass() const
{
   return std::any_of(bad_masses.cbegin(), bad_masses.cend(),
                      [](bool x) { return x; });
}

bool Problems::have_running_tachyon() const
{
   return std::any_of(running_tachyons.cbegin(), running_tachyons.cend(),
                      [](bool x) { return x; });
}

bool Problems::have_pole_tachyon() const
{
   return std::any_of(pole_tachyons.cbegin(), pole_tachyons.cend(),
                      [](bool x) { return x; });
}

bool Problems::have_tachyon() const
{
   return have_running_tachyon() || have_pole_tachyon();
}

bool Problems::have_thrown() const
{
   return !exception_msg.empty();
}

bool Problems::have_non_perturbative_parameter() const
{
   return !non_pert_pars.empty();
}

bool Problems::have_failed_pole_mass_convergence() const
{
   return std::any_of(failed_pole_mass_convergence.cbegin(),
                      failed_pole_mass_convergence.cend(),
                      [](bool x) { return x; });
}

bool Problems::no_ewsb() const
{
   return failed_ewsb;
}

bool Problems::no_ewsb_tree_level() const
{
   return failed_ewsb_tree_level;
}

bool Problems::no_perturbative() const
{
   return non_perturbative;
}

bool Problems::no_sinThetaW_convergence() const
{
   return failed_sinThetaW_convergence;
}

std::vector<int> Problems::get_bad_masses() const
{
   return bad_masses;
}

std::vector<int> Problems::get_running_tachyons() const
{
   return running_tachyons;
}

std::vector<int> Problems::get_pole_tachyons() const
{
   return pole_tachyons;
}

std::vector<int> Problems::get_failed_pole_mass_convergence() const
{
   return failed_pole_mass_convergence;
}

std::ostream& operator<<(std::ostream& ostr, const Problems& problems)
{
   problems.print_problems(ostr);
   problems.print_warnings(ostr);
   return ostr;
}

} // namespace flexiblesusy
