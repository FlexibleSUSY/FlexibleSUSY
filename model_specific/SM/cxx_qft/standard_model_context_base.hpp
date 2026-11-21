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
 * @file cxx_qft/SM_context_base.hpp
 *
 * This file was generated with FlexibleSUSY 2.7.1 and SARAH 4.14.5 .
 */

#ifndef STANDARDMODEL_CXXQFT_CONTEXT_BASE_H
#define STANDARDMODEL_CXXQFT_CONTEXT_BASE_H

#include "standard_model_fields.hpp"
#include "standard_model.hpp"

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

} // namespace standard_model_cxx_diagrams
} // namespace flexiblesusy

#include "standard_model_vertices.hpp"

#endif
