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

#ifndef PROBLEMS_FORMAT_MATHLINK_H
#define PROBLEMS_FORMAT_MATHLINK_H

#include "mathlink_utils.hpp"
#include "problems.hpp"
#include <vector>

namespace flexiblesusy {

namespace {

template <class ParticleNameGetter>
void put_masses(MLINK link, const std::vector<int>& flags,
                const std::vector<std::string>& heads,
                const ParticleNameGetter& particle_name)
{
   for (std::size_t i = 0; i < flags.size(); i++) {
      if (flags[i]) {
         MLPutHeads(link, heads);
         MLPutSymbol(link, particle_name(i).c_str());
      }
   }
}

template <class ParticleNameGetter>
void put_masses(MLINK link, const std::string& rule,
                const std::vector<int>& flags,
                const std::vector<std::string>& heads,
                const ParticleNameGetter& particle_name)
{
   const auto n_masses = std::count(flags.cbegin(), flags.cend(), true);

   MLPutRule(link, rule);
   MLPutFunction(link, "List", n_masses);
   put_masses(link, flags, heads, particle_name);
}

template <class ParticleNameGetter>
void put_masses(MLINK link, const std::string& rule,
                const std::vector<int>& flags1,
                const std::vector<std::string>& heads1,
                const std::vector<int>& flags2,
                const std::vector<std::string>& heads2,
                const ParticleNameGetter& particle_name)
{
   const auto n_masses = std::count(flags1.cbegin(), flags1.cend(), true)
                       + std::count(flags2.cbegin(), flags2.cend(), true);

   MLPutRule(link, rule);
   MLPutFunction(link, "List", n_masses);
   put_masses(link, flags1, heads1, particle_name);
   put_masses(link, flags2, heads2, particle_name);
}

} // anonymous namespace

/// format problems to MathLink output
inline void mathlink_format_problems(MLINK link, const Problems& pr)
{
   const auto pn = [&pr] (int i) { return pr.get_particle_name(i); };

   MLPutFunction(link, "List", pr.number_of_problems());

   if (pr.have_tachyon()) {
      put_masses(link, "Tachyons",
                 pr.get_running_tachyons(), {"M"},
                 pr.get_pole_tachyons(), {"Pole", "M"},
                 pn);
   }
   if (pr.no_ewsb()) {
      MLPutRuleTo(link, "True", "NoEWSB");
   }
   if (pr.no_perturbative() || pr.have_non_perturbative_parameter()) {
      MLPutRuleTo(link, "True", "NonPerturbative");
   }
   if (pr.no_sinThetaW_convergence()) {
      MLPutRuleTo(link, "True", "NoSinThetaWConvergence");
   }
   if (pr.have_thrown()) {
      MLPutRuleTo(link, "True", "Exceptions");
   }
   if (pr.have_failed_pole_mass_convergence()) {
      put_masses(link, "NoPoleMassConvergence",
                 pr.get_failed_pole_mass_convergence(), {"Pole", "M"},
                 pn);
   }
}

/// format warnings to MathLink output
inline void mathlink_format_warnings(MLINK link, const Problems& pr)
{
   const auto pn = [&pr] (int i) { return pr.get_particle_name(i); };

   MLPutFunction(link, "List", pr.number_of_warnings());

   if (pr.have_bad_mass()) {
      put_masses(link, "ImpreciseMasses", pr.get_bad_masses(), {"M"}, pn);
   }
}

} // namespace flexiblesusy

#endif
