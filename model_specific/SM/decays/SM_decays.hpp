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
 * @file SM_decays.hpp
 *
 * @brief contains class for calculating particle decays
 *
 * This file was generated with FlexibleSUSY 2.6.0 and SARAH 4.14.3 .
 */

#ifndef SM_DECAYS_H
#define SM_DECAYS_H

#include "SM_decay_table.hpp"
#include "SM_mass_eigenstates.hpp"
#include "SM_mass_eigenstates_decoupling_scheme.hpp"
#include "cxx_qft/SM_qft.hpp"
#include "SM_decay_amplitudes.hpp"
#include "decays/flexibledecay_problems.hpp"
#include "lowe.h"
#include "wrappers.hpp"
#include "error.hpp"
#include "physical_input.hpp"
#include "decays/flexibledecay_settings.hpp"

namespace flexiblesusy {

template <typename Field1, typename Field2>
constexpr std::enable_if_t<!std::is_same<Field1, Field2>::value, double>
final_state_symmetry_factor(typename SM_cxx_diagrams::field_indices<Field1>::type const& idx1,
                            typename SM_cxx_diagrams::field_indices<Field2>::type const& idx2)
{
   return 1.;
}

template <typename Field1, typename Field2>
std::enable_if_t<std::is_same<Field1, Field2>::value, double>
final_state_symmetry_factor(typename SM_cxx_diagrams::field_indices<Field1>::type const& idx1,
                            typename SM_cxx_diagrams::field_indices<Field2>::type const& idx2)
{
   if (boost::range::equal(idx1, idx2)) {
      return 0.5;
   }
   else {
      return 1.;
   }
}

class SM_decays {
public:
   SM_decays() = default;
   SM_decays(SM_mass_eigenstates model_, softsusy::QedQcd const& qedqcd_,
         Physical_input const& physical_input_,
         FlexibleDecay_settings const& flexibledecay_settings_)
      : model(model_)
      , qedqcd(qedqcd_)
      , physical_input(physical_input_)
      , flexibledecay_settings(flexibledecay_settings_)
      {}
   SM_decays(const SM_decays&) = default;
   SM_decays(SM_decays&&) = default;
   ~SM_decays() = default;
   SM_decays& operator=(const SM_decays&) = default;
   SM_decays& operator=(SM_decays&&) = default;

   const SM_decay_table& get_decay_table() const;
   const FlexibleDecay_problems& get_problems() const;

   void clear();
   void clear_problems();
   void calculate_decays();

   const Decays_list& get_hh_decays() const { return decay_table.get_hh_decays();
      }
   void calculate_hh_decays();

double partial_width_hh_to_VGVG(SM_mass_eigenstates_interface*) const;
double partial_width_hh_to_VPVP(SM_mass_eigenstates_interface*) const;
double partial_width_hh_to_VPVZ(SM_mass_eigenstates_interface*) const;
double partial_width_hh_to_VZVZ(SM_mass_eigenstates_interface*) const;
double partial_width_hh_to_conjVWpVWp(SM_mass_eigenstates_interface*) const;
double partial_width_hh_to_barFvFv(SM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_barFdFd(SM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_barFuFu(SM_mass_eigenstates_interface*, int, int) const;
double partial_width_hh_to_barFeFe(SM_mass_eigenstates_interface*, int, int) const;

private:
   SM_mass_eigenstates model{};
   softsusy::QedQcd qedqcd{};
   Physical_input physical_input;
   FlexibleDecay_settings flexibledecay_settings {};
   bool run_to_decay_particle_scale {true};
   SM_decay_table decay_table{};
   FlexibleDecay_problems problems{};

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   typename Decay_amplitude_type<FieldIn, FieldOut1, FieldOut2>::type
   calculate_amplitude(
      const SM_cxx_diagrams::context_base&,
      const typename SM_cxx_diagrams::field_indices<FieldIn>::type&,
      const typename SM_cxx_diagrams::field_indices<FieldOut1>::type&,
      const typename SM_cxx_diagrams::field_indices<FieldOut2>::type&) const;

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   double amplitude_squared(SM_cxx_diagrams::context_base const& context,
                  typename SM_cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
                  typename SM_cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
                  typename SM_cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2) const;

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   double get_partial_width(
      const SM_cxx_diagrams::context_base&,
      typename SM_cxx_diagrams::field_indices<FieldIn>::type const&,
      typename SM_cxx_diagrams::field_indices<FieldOut1>::type const&,
      typename SM_cxx_diagrams::field_indices<FieldOut2>::type const&) const;
};

template<>
Decay_amplitude_SVV SM_decays::calculate_amplitude<SM_cxx_diagrams::fields::hh, SM_cxx_diagrams::fields::VG, SM_cxx_diagrams::fields::VG>(const SM_cxx_diagrams::context_base&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::hh >::type&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::VG >::type&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SVV SM_decays::calculate_amplitude<SM_cxx_diagrams::fields::hh, SM_cxx_diagrams::fields::VP, SM_cxx_diagrams::fields::VP>(const SM_cxx_diagrams::context_base&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::hh >::type&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::VP >::type&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SVV SM_decays::calculate_amplitude<SM_cxx_diagrams::fields::hh, SM_cxx_diagrams::fields::VP, SM_cxx_diagrams::fields::VZ>(const SM_cxx_diagrams::context_base&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::hh >::type&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::VP >::type&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV SM_decays::calculate_amplitude<SM_cxx_diagrams::fields::hh, SM_cxx_diagrams::fields::VZ, SM_cxx_diagrams::fields::VZ>(const SM_cxx_diagrams::context_base&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::hh >::type&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::VZ >::type&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV SM_decays::calculate_amplitude<SM_cxx_diagrams::fields::hh, typename SM_cxx_diagrams::fields::conj<SM_cxx_diagrams::fields::VWp>::type, SM_cxx_diagrams::fields::VWp>(const SM_cxx_diagrams::context_base&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::hh >::type&, const typename SM_cxx_diagrams::field_indices<typename SM_cxx_diagrams::fields::conj<SM_cxx_diagrams::fields::VWp>::type >::type&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::VWp >::type&) const;

template<>
Decay_amplitude_SFF SM_decays::calculate_amplitude<SM_cxx_diagrams::fields::hh, typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fv>::type, SM_cxx_diagrams::fields::Fv>(const SM_cxx_diagrams::context_base&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::hh >::type&, const typename SM_cxx_diagrams::field_indices<typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fv>::type >::type&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::Fv >::type&) const;

template<>
Decay_amplitude_SFF SM_decays::calculate_amplitude<SM_cxx_diagrams::fields::hh, typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fd>::type, SM_cxx_diagrams::fields::Fd>(const SM_cxx_diagrams::context_base&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::hh >::type&, const typename SM_cxx_diagrams::field_indices<typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fd>::type >::type&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SFF SM_decays::calculate_amplitude<SM_cxx_diagrams::fields::hh, typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fu>::type, SM_cxx_diagrams::fields::Fu>(const SM_cxx_diagrams::context_base&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::hh >::type&, const typename SM_cxx_diagrams::field_indices<typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fu>::type >::type&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::Fu >::type&) const;

template<>
Decay_amplitude_SFF SM_decays::calculate_amplitude<SM_cxx_diagrams::fields::hh, typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fe>::type, SM_cxx_diagrams::fields::Fe>(const SM_cxx_diagrams::context_base&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::hh >::type&, const typename SM_cxx_diagrams::field_indices<typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fe>::type >::type&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::Fe >::type&) const;


template<typename FieldIn, typename FieldOut1, typename FieldOut2>
double
SM_decays::amplitude_squared(SM_cxx_diagrams::context_base const& context,
                  typename SM_cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
                  typename SM_cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
                  typename SM_cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2) const
{

   const auto mat_elem = calculate_amplitude<FieldIn, FieldOut1, FieldOut2>(
      context, indexIn, indexOut1, indexOut2);
   return mat_elem.square();
}

// returns a squared color generator for a 3 point amplitude with FieldIn, FieldOut1 and FieldOut2
// averaged over inital state colors
// the generator is guessed from color representations of FieldIn, FieldOut1 and FieldOut2
// This is not a bulletproof solution and might fail in general but is enough for
// decays of color singlets

// 1 -> 1, 1
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
SM_cxx_diagrams::fields::is_singlet_v<FieldIn> &&
SM_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
SM_cxx_diagrams::fields::is_singlet_v<FieldOut2>, double>
squared_color_generator() {return 1.;}

// 1 -> 3, 3bar
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
   SM_cxx_diagrams::fields::is_singlet_v<FieldIn>
   &&
   (
   (SM_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
   SM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>)
   ||
   (SM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
   SM_cxx_diagrams::fields::is_triplet_v<FieldOut2>)
   ), double
>
squared_color_generator() {return 3.;}

// 1 -> 8, 8
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
SM_cxx_diagrams::fields::is_singlet_v<FieldIn> &&
SM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
SM_cxx_diagrams::fields::is_octet_v<FieldOut2>, double>
squared_color_generator() {return 8.;}

// 3 -> 3, 1; 3bar -> 3bar, 1
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
(
SM_cxx_diagrams::fields::is_triplet_v<FieldIn> &&
((SM_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
SM_cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(SM_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
SM_cxx_diagrams::fields::is_triplet_v<FieldOut2>))
)
||
(
SM_cxx_diagrams::fields::is_anti_triplet_v<FieldIn> &&
((SM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
SM_cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(SM_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
SM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>))
)
, double>
squared_color_generator() {return 1.;}

// 3 -> 3, 8; 3bar -> 3bar, 8
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
(
SM_cxx_diagrams::fields::is_triplet_v<FieldIn> &&
((SM_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
SM_cxx_diagrams::fields::is_octet_v<FieldOut2>) ||
(SM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
SM_cxx_diagrams::fields::is_triplet_v<FieldOut2>))
)
||
(
SM_cxx_diagrams::fields::is_anti_triplet_v<FieldIn> &&
((SM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
SM_cxx_diagrams::fields::is_octet_v<FieldOut2>) ||
(SM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
SM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>))
)
, double>
squared_color_generator() {return 4.;}

// 8 -> 8, 1
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
SM_cxx_diagrams::fields::is_octet_v<FieldIn> &&
((SM_cxx_diagrams::fields::is_octet_v<FieldOut1> &&
SM_cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(SM_cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
SM_cxx_diagrams::fields::is_octet_v<FieldOut2>))
, double>
squared_color_generator() {return 1.;}

// 8 -> 3, 3bar
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
   SM_cxx_diagrams::fields::is_octet_v<FieldIn>
   &&
   (
   (SM_cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
   SM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>)
   ||
   (SM_cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
   SM_cxx_diagrams::fields::is_triplet_v<FieldOut2>)
   ), double
>
squared_color_generator() {return 1./2.;}

// generic decay of FieldIn -> FieldOut1 FieldOut2
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
double SM_decays::get_partial_width(
   const SM_cxx_diagrams::context_base& context,
   typename SM_cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
   typename SM_cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
   typename SM_cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2
   ) const
{

   // on-shell masses
   const double mIn = context.physical_mass<FieldIn>(indexIn);
   const double mOut1 = context.physical_mass<FieldOut1>(indexOut1);
   const double mOut2 = context.physical_mass<FieldOut2>(indexOut2);

   // is decay kinematically allowed?
   if(mIn < mOut1 + mOut2) {
      WARNING("Called kinematically forbidden decay");
      return 0.;
   }

   // phase space without symmetry factor
   const double ps = 1./(8.*Pi) * std::sqrt(KallenLambda(1., Sqr(mOut1/mIn), Sqr(mOut2/mIn)));

   // phase space symmetry factor
   const double ps_symmetry =
      final_state_symmetry_factor<FieldOut1, FieldOut2>(indexOut1, indexOut2);

   // color factor
   constexpr double color_factor = squared_color_generator<FieldIn, FieldOut1, FieldOut2>();

   // matrix element squared
   const auto mat_elem_sq = amplitude_squared<FieldIn, FieldOut1, FieldOut2>(
      context, indexIn, indexOut1, indexOut2);

   // flux * phase space factor * symmetry factor * color factor * |matrix element|^2
   const auto result = 0.5/mIn * ps * ps_symmetry * color_factor * mat_elem_sq;

   return result;
}

template <>
double SM_decays::get_partial_width<SM_cxx_diagrams::fields::hh,SM_cxx_diagrams::fields::VZ,SM_cxx_diagrams::fields::VZ >(const SM_cxx_diagrams::context_base&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::hh >::type&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::VZ >::type&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::VZ >::type&) const;
template <>
double SM_decays::get_partial_width<SM_cxx_diagrams::fields::hh,typename SM_cxx_diagrams::fields::conj<SM_cxx_diagrams::fields::VWp>::type,SM_cxx_diagrams::fields::VWp >(const SM_cxx_diagrams::context_base&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::hh >::type&, const typename SM_cxx_diagrams::field_indices<typename SM_cxx_diagrams::fields::conj<SM_cxx_diagrams::fields::VWp>::type >::type&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::VWp >::type&) const;
template <>
double SM_decays::get_partial_width<SM_cxx_diagrams::fields::hh,SM_cxx_diagrams::fields::VG,SM_cxx_diagrams::fields::VG >(const SM_cxx_diagrams::context_base&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::hh >::type&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::VG >::type&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::VG >::type&) const;
template <>
double SM_decays::get_partial_width<SM_cxx_diagrams::fields::hh,SM_cxx_diagrams::fields::VP,SM_cxx_diagrams::fields::VP >(const SM_cxx_diagrams::context_base&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::hh >::type&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::VP >::type&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::VP >::type&) const;
template <>
double SM_decays::get_partial_width<SM_cxx_diagrams::fields::hh,SM_cxx_diagrams::fields::VP,SM_cxx_diagrams::fields::VZ >(const SM_cxx_diagrams::context_base&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::hh >::type&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::VP >::type&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::VZ >::type&) const;
template <>
double SM_decays::get_partial_width<SM_cxx_diagrams::fields::hh,typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fu>::type,SM_cxx_diagrams::fields::Fu >(const SM_cxx_diagrams::context_base&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::hh >::type&, const typename SM_cxx_diagrams::field_indices<typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fu>::type >::type&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::Fu >::type&) const;
template <>
double SM_decays::get_partial_width<SM_cxx_diagrams::fields::hh,typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fd>::type,SM_cxx_diagrams::fields::Fd >(const SM_cxx_diagrams::context_base&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::hh >::type&, const typename SM_cxx_diagrams::field_indices<typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fd>::type >::type&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::Fd >::type&) const;
template <>
double SM_decays::get_partial_width<SM_cxx_diagrams::fields::hh,typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fe>::type,SM_cxx_diagrams::fields::Fe >(const SM_cxx_diagrams::context_base&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::hh >::type&, const typename SM_cxx_diagrams::field_indices<typename SM_cxx_diagrams::fields::bar<SM_cxx_diagrams::fields::Fe>::type >::type&, const typename SM_cxx_diagrams::field_indices<SM_cxx_diagrams::fields::Fe >::type&) const;

} // namespace flexiblesusy

#endif
