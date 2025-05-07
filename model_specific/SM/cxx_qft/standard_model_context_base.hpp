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


/**
 * @file cxx_qft/standatd_model_context_base.hpp
 *
 * This file was generated with FlexibleSUSY 2.7.1 and SARAH 4.14.5 .
 */

#ifndef STANDARDMODEL_CXXQFT_CONTEXT_BASE_H
#define STANDARDMODEL_CXXQFT_CONTEXT_BASE_H

#include "standard_model_fields.hpp"
#include "standard_model.hpp"
#include "always_false.hpp"

namespace flexiblesusy {
namespace standard_model_cxx_diagrams {

using cxx_diagrams::field_indices;
using standard_model::Standard_model;

   struct context_base {
      Standard_model const& model; ///< The model object.

      template <class Field>
      double mass(const typename cxx_diagrams::field_indices<Field>::type& indices) const
      {
         using CleanField =
            typename cxx_diagrams::fields::remove_lorentz_conjugation<Field>::type;
         return mass_impl<CleanField>(indices);
      }
      template<class Field>
      double physical_mass(const typename cxx_diagrams::field_indices<Field>::type& indices) const
      {
         using CleanField =
            typename cxx_diagrams::fields::remove_lorentz_conjugation<Field>::type;
         return physical_mass_impl<CleanField>(indices);
      }

      context_base(Standard_model const& m) : model(m) {}
      context_base(Standard_model const* const m) : model(*m) {}

      context_base(const context_base&) = default;
      context_base(context_base&&) = default;

      virtual ~context_base(void) = default;

   private:
      template <class Field>
      double
      mass_impl(const typename cxx_diagrams::field_indices<Field>::type& indices) const;
      template<class Field>
      double physical_mass_impl(const typename cxx_diagrams::field_indices<Field>::type& indices) const;
   };

template<> inline
double context_base::mass_impl<fields::VG>(const std::array<int, 0>& indices) const
{ return model.get_MVG(); }

template<> inline
double context_base::mass_impl<fields::gG>(const std::array<int, 0>& indices) const
{ return model.get_MVG(); }

template<> inline
double context_base::mass_impl<fields::Hp>(const std::array<int, 0>& indices) const
{ return model.get_MHp(); }

template<> inline
double context_base::mass_impl<fields::Fv>(const std::array<int, 1>& indices) const
{ return model.get_MFv(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Ah>(const std::array<int, 0>& indices) const
{ return model.get_MAh(); }

template<> inline
double context_base::mass_impl<fields::hh>(const std::array<int, 0>& indices) const
{ return model.get_Mhh(); }

template<> inline
double context_base::mass_impl<fields::VP>(const std::array<int, 0>& indices) const
{ return model.get_MVP(); }

template<> inline
double context_base::mass_impl<fields::VZ>(const std::array<int, 0>& indices) const
{ return model.get_MVZ(); }

template<> inline
double context_base::mass_impl<fields::gP>(const std::array<int, 0>& indices) const
{ return model.get_MVP(); }

template<> inline
double context_base::mass_impl<fields::gZ>(const std::array<int, 0>& indices) const
{ return model.get_MVZ(); }

template<> inline
double context_base::mass_impl<fields::gWp>(const std::array<int, 0>& indices) const
{ return model.get_MVWp(); }

template<> inline
double context_base::mass_impl<fields::gWpC>(const std::array<int, 0>& indices) const
{ return model.get_MVWp(); }

template<> inline
double context_base::mass_impl<fields::Fd>(const std::array<int, 1>& indices) const
{ return model.get_MFd(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Fu>(const std::array<int, 1>& indices) const
{ return model.get_MFu(indices[0]); }

template<> inline
double context_base::mass_impl<fields::Fe>(const std::array<int, 1>& indices) const
{ return model.get_MFe(indices[0]); }

template<> inline
double context_base::mass_impl<fields::VWp>(const std::array<int, 0>& indices) const
{ return model.get_MVWp(); }

template<> inline
double context_base::physical_mass_impl<fields::VG>(const std::array<int, 0>& indices) const
{ return model.get_physical().MVG; }

template<> inline
double context_base::physical_mass_impl<fields::gG>(const std::array<int, 0>& indices) const
{ return model.get_physical().MVG; }

template<> inline
double context_base::physical_mass_impl<fields::Hp>(const std::array<int, 0>& indices) const
{ return model.get_physical().MHp; }

template<> inline
double context_base::physical_mass_impl<fields::Fv>(const std::array<int, 1>& indices) const
{ return model.get_physical().MFv[indices[0]]; }

template<> inline
double context_base::physical_mass_impl<fields::Ah>(const std::array<int, 0>& indices) const
{ return model.get_physical().MAh; }

template<> inline
double context_base::physical_mass_impl<fields::hh>(const std::array<int, 0>& indices) const
{ return model.get_physical().Mhh; }

template<> inline
double context_base::physical_mass_impl<fields::VP>(const std::array<int, 0>& indices) const
{ return model.get_physical().MVP; }

template<> inline
double context_base::physical_mass_impl<fields::VZ>(const std::array<int, 0>& indices) const
{ return model.get_physical().MVZ; }

template<> inline
double context_base::physical_mass_impl<fields::gP>(const std::array<int, 0>& indices) const
{ return model.get_physical().MVP; }

template<> inline
double context_base::physical_mass_impl<fields::gZ>(const std::array<int, 0>& indices) const
{ return model.get_physical().MVZ; }

template<> inline
double context_base::physical_mass_impl<fields::gWp>(const std::array<int, 0>& indices) const
{ return model.get_physical().MVWp; }

template<> inline
double context_base::physical_mass_impl<fields::gWpC>(const std::array<int, 0>& indices) const
{ return model.get_physical().MVWp; }

template<> inline
double context_base::physical_mass_impl<fields::Fd>(const std::array<int, 1>& indices) const
{ return model.get_physical().MFd[indices[0]]; }

template<> inline
double context_base::physical_mass_impl<fields::Fu>(const std::array<int, 1>& indices) const
{ return model.get_physical().MFu[indices[0]]; }

template<> inline
double context_base::physical_mass_impl<fields::Fe>(const std::array<int, 1>& indices) const
{ return model.get_physical().MFe[indices[0]]; }

template<> inline
double context_base::physical_mass_impl<fields::VWp>(const std::array<int, 0>& indices) const
{ return model.get_physical().MVWp; }

// self energies
template<typename Field>
auto self_energy_1loop(const context_base& context, double p) {
   static_assert(always_false<Field>);
   return 0.;
}

template<typename Field>
auto self_energy_1loop_1(const context_base& context, double p) {
   static_assert(always_false<Field>);
   return 0.;
}

template<typename Field>
auto self_energy_1loop_PL(const context_base& context, double p) {
   static_assert(always_false<Field>);
   return 0.;
}

template<typename Field>
auto self_energy_1loop_PR(const context_base& context, double p) {
   static_assert(always_false<Field>);
   return 0.;
}

// self-energy derivatives w.r.t. p2
template<typename Field>
auto self_energy_1loop_deriv_p2(const context_base& context, double p) {
   static_assert(always_false<Field>);
   return 0.;
}

template<typename Field>
auto self_energy_1loop_1_deriv_p2(const context_base& context, double p) {
   static_assert(always_false<Field>);
   return 0.;
}

template<typename Field>
auto self_energy_1loop_PL_deriv_p2(const context_base& context, double p) {
   static_assert(always_false<Field>);
   return 0.;
}

template<typename Field>
auto self_energy_1loop_PR_deriv_p2(const context_base& context, double p) {
   static_assert(always_false<Field>);
   return 0.;
}

template<> inline
auto self_energy_1loop<fields::VG>(const context_base& context, double p) {
   return context.model.self_energy_VG_1loop(p);
}

template<> inline
auto self_energy_1loop<fields::Hp>(const context_base& context, double p) {
   return context.model.self_energy_Hp_1loop(p);
}

template<> inline
auto self_energy_1loop_1<fields::Fv>(const context_base& context, double p) {
   return context.model.self_energy_Fv_1loop_1(p);
}

template<> inline
auto self_energy_1loop_PL<fields::Fv>(const context_base& context, double p) {
   return context.model.self_energy_Fv_1loop_PL(p);
}

template<> inline
auto self_energy_1loop_PR<fields::Fv>(const context_base& context, double p) {
   return context.model.self_energy_Fv_1loop_PR(p);
}

template<> inline
auto self_energy_1loop<fields::Ah>(const context_base& context, double p) {
   return context.model.self_energy_Ah_1loop(p);
}

template<> inline
auto self_energy_1loop<fields::hh>(const context_base& context, double p) {
   return context.model.self_energy_hh_1loop(p);
}

template<> inline
auto self_energy_1loop<fields::VP>(const context_base& context, double p) {
   return context.model.self_energy_VP_1loop(p);
}

template<> inline
auto self_energy_1loop<fields::VZ>(const context_base& context, double p) {
   return context.model.self_energy_VZ_1loop(p);
}

template<> inline
auto self_energy_1loop_1<fields::Fd>(const context_base& context, double p) {
   return context.model.self_energy_Fd_1loop_1(p);
}

template<> inline
auto self_energy_1loop_PL<fields::Fd>(const context_base& context, double p) {
   return context.model.self_energy_Fd_1loop_PL(p);
}

template<> inline
auto self_energy_1loop_PR<fields::Fd>(const context_base& context, double p) {
   return context.model.self_energy_Fd_1loop_PR(p);
}

template<> inline
auto self_energy_1loop_1<fields::Fu>(const context_base& context, double p) {
   return context.model.self_energy_Fu_1loop_1(p);
}

template<> inline
auto self_energy_1loop_PL<fields::Fu>(const context_base& context, double p) {
   return context.model.self_energy_Fu_1loop_PL(p);
}

template<> inline
auto self_energy_1loop_PR<fields::Fu>(const context_base& context, double p) {
   return context.model.self_energy_Fu_1loop_PR(p);
}

template<> inline
auto self_energy_1loop_1<fields::Fe>(const context_base& context, double p) {
   return context.model.self_energy_Fe_1loop_1(p);
}

template<> inline
auto self_energy_1loop_PL<fields::Fe>(const context_base& context, double p) {
   return context.model.self_energy_Fe_1loop_PL(p);
}

template<> inline
auto self_energy_1loop_PR<fields::Fe>(const context_base& context, double p) {
   return context.model.self_energy_Fe_1loop_PR(p);
}

template<> inline
auto self_energy_1loop<fields::VWp>(const context_base& context, double p) {
   return context.model.self_energy_VWp_1loop(p);
}

template<> inline
auto self_energy_1loop_deriv_p2<fields::VG>(const context_base& context, double p) {
   return context.model.self_energy_VG_1loop_deriv_p2(p);
}

template<> inline
auto self_energy_1loop_deriv_p2<fields::Hp>(const context_base& context, double p) {
   return context.model.self_energy_Hp_1loop_deriv_p2(p);
}

template<> inline
auto self_energy_1loop_1_deriv_p2<fields::Fv>(const context_base& context, double p) {
   return context.model.self_energy_Fv_1loop_1_deriv_p2(p);
}

template<> inline
auto self_energy_1loop_PL_deriv_p2<fields::Fv>(const context_base& context, double p) {
   return context.model.self_energy_Fv_1loop_PL_deriv_p2(p);
}

template<> inline
auto self_energy_1loop_PR_deriv_p2<fields::Fv>(const context_base& context, double p) {
   return context.model.self_energy_Fv_1loop_PR_deriv_p2(p);
}

template<> inline
auto self_energy_1loop_deriv_p2<fields::Ah>(const context_base& context, double p) {
   return context.model.self_energy_Ah_1loop_deriv_p2(p);
}

template<> inline
auto self_energy_1loop_deriv_p2<fields::hh>(const context_base& context, double p) {
   return context.model.self_energy_hh_1loop_deriv_p2(p);
}

template<> inline
auto self_energy_1loop_deriv_p2<fields::VP>(const context_base& context, double p) {
   return context.model.self_energy_VP_1loop_deriv_p2(p);
}

template<> inline
auto self_energy_1loop_deriv_p2<fields::VZ>(const context_base& context, double p) {
   return context.model.self_energy_VZ_1loop_deriv_p2(p);
}

template<> inline
auto self_energy_1loop_1_deriv_p2<fields::Fd>(const context_base& context, double p) {
   return context.model.self_energy_Fd_1loop_1_deriv_p2(p);
}

template<> inline
auto self_energy_1loop_PL_deriv_p2<fields::Fd>(const context_base& context, double p) {
   return context.model.self_energy_Fd_1loop_PL_deriv_p2(p);
}

template<> inline
auto self_energy_1loop_PR_deriv_p2<fields::Fd>(const context_base& context, double p) {
   return context.model.self_energy_Fd_1loop_PR_deriv_p2(p);
}

template<> inline
auto self_energy_1loop_1_deriv_p2<fields::Fu>(const context_base& context, double p) {
   return context.model.self_energy_Fu_1loop_1_deriv_p2(p);
}

template<> inline
auto self_energy_1loop_PL_deriv_p2<fields::Fu>(const context_base& context, double p) {
   return context.model.self_energy_Fu_1loop_PL_deriv_p2(p);
}

template<> inline
auto self_energy_1loop_PR_deriv_p2<fields::Fu>(const context_base& context, double p) {
   return context.model.self_energy_Fu_1loop_PR_deriv_p2(p);
}

template<> inline
auto self_energy_1loop_1_deriv_p2<fields::Fe>(const context_base& context, double p) {
   return context.model.self_energy_Fe_1loop_1_deriv_p2(p);
}

template<> inline
auto self_energy_1loop_PL_deriv_p2<fields::Fe>(const context_base& context, double p) {
   return context.model.self_energy_Fe_1loop_PL_deriv_p2(p);
}

template<> inline
auto self_energy_1loop_PR_deriv_p2<fields::Fe>(const context_base& context, double p) {
   return context.model.self_energy_Fe_1loop_PR_deriv_p2(p);
}

template<> inline
auto self_energy_1loop_deriv_p2<fields::VWp>(const context_base& context, double p) {
   return context.model.self_energy_VWp_1loop_deriv_p2(p);
}

} // namespace standard_model_cxx_diagrams
} // namespace flexiblesusy

#include "standard_model_vertices.hpp"

#endif
