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
 * @file standard_model_decays.hpp
 *
 * @brief contains class for calculating particle decays
 *
 * This file was generated with FlexibleSUSY 2.7.1 and SARAH 4.14.5 .
 */

#ifndef STANDARDMODEL_DECAYS_H
#define STANDARDMODEL_DECAYS_H

#include "lowe.h"
#include "wrappers.hpp"
#include "error.hpp"
#include "physical_input.hpp"
#include "decays/flexibledecay_settings.hpp"
#include "decays/flexibledecay_problems.hpp"

#include "standard_model.hpp"
#include "cxx_qft/standard_model_qft.hpp"
#include "decays/decay_amplitudes.hpp"
#include "standard_model_decay_table.hpp"
#include "config.h"

#include <cmath>
#include <string>
#include <tuple>
#include <vector>

namespace flexiblesusy {

template<typename FieldIn, typename FieldOut1, typename FieldOut2>
struct standard_model_has_oneloop_amp : public std::false_type {};

template<>
struct standard_model_has_oneloop_amp<standard_model_cxx_diagrams::fields::hh, standard_model_cxx_diagrams::fields::VG, standard_model_cxx_diagrams::fields::VG> : public std::true_type {};
template<>
struct standard_model_has_oneloop_amp<standard_model_cxx_diagrams::fields::hh, standard_model_cxx_diagrams::fields::VP, standard_model_cxx_diagrams::fields::VP> : public std::true_type {};
template<>
struct standard_model_has_oneloop_amp<standard_model_cxx_diagrams::fields::hh, standard_model_cxx_diagrams::fields::VP, standard_model_cxx_diagrams::fields::VZ> : public std::true_type {};
template<>
struct standard_model_has_oneloop_amp<standard_model_cxx_diagrams::fields::hh, typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fv>::type, standard_model_cxx_diagrams::fields::Fv> : public std::true_type {};


template<typename FieldIn, typename FieldOut1, typename FieldOut2>
struct standard_model_has_tree_amp : public std::false_type {};

template<>
struct standard_model_has_tree_amp<standard_model_cxx_diagrams::fields::hh, standard_model_cxx_diagrams::fields::VZ, standard_model_cxx_diagrams::fields::VZ> : public std::true_type {};
template<>
struct standard_model_has_tree_amp<standard_model_cxx_diagrams::fields::hh, standard_model_cxx_diagrams::fields::VWp, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::VWp>::type> : public std::true_type {};
template<>
struct standard_model_has_tree_amp<standard_model_cxx_diagrams::fields::hh, typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fd>::type, standard_model_cxx_diagrams::fields::Fd> : public std::true_type {};
template<>
struct standard_model_has_tree_amp<standard_model_cxx_diagrams::fields::hh, typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fu>::type, standard_model_cxx_diagrams::fields::Fu> : public std::true_type {};
template<>
struct standard_model_has_tree_amp<standard_model_cxx_diagrams::fields::hh, typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fe>::type, standard_model_cxx_diagrams::fields::Fe> : public std::true_type {};

class Standard_model_decays {
public:
   Standard_model_decays() = default;
   Standard_model_decays(standard_model::Standard_model model_, softsusy::QedQcd const& qedqcd_,
         Physical_input const& physical_input_,
         FlexibleDecay_settings const& flexibledecay_settings_)
      : model(model_)
      , qedqcd(qedqcd_)
      , physical_input(physical_input_)
      , flexibledecay_settings(flexibledecay_settings_)
      {
         sm.initialise_from_input(qedqcd);
      }
   Standard_model_decays(const Standard_model_decays&) = default;
   Standard_model_decays(Standard_model_decays&&) = default;
   ~Standard_model_decays() = default;
   Standard_model_decays& operator=(const Standard_model_decays&) = default;
   Standard_model_decays& operator=(Standard_model_decays&&) = default;

   const standard_model_decay_table& get_decay_table() const;
   const FlexibleDecay_problems& get_problems() const;

   void clear();
   void clear_problems();
   void calculate_decays();
   EffectiveCoupling_list const& get_neutral_higgs_effc() const {return neutral_higgs_effc;};
   std::vector<std::tuple<int, int, int, double, std::string>> get_effhiggscouplings_block_input() const {
      return effhiggscouplings_block_input;
   }

   const Decays_list& get_hh_decays() const { return decay_table.get_hh_decays();
      }
   void calculate_hh_decays();

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   double get_partial_width(
      const standard_model_cxx_diagrams::context_base&,
      typename cxx_diagrams::field_indices<FieldIn>::type const&,
      typename cxx_diagrams::field_indices<FieldOut1>::type const&,
      typename cxx_diagrams::field_indices<FieldOut2>::type const&);
   double get_alphas(standard_model_cxx_diagrams::context_base const&) const;
   double get_alpha(standard_model_cxx_diagrams::context_base const&) const;

double partial_width_hh_to_VGVG(standard_model::Standard_model const&);
double partial_width_hh_to_VPVP(standard_model::Standard_model const&);
double partial_width_hh_to_VPVZ(standard_model::Standard_model const&);
double partial_width_hh_to_VZVZ(standard_model::Standard_model const&);
double partial_width_hh_to_conjVWpVWp(standard_model::Standard_model const&);
double partial_width_hh_to_barFvFv(standard_model::Standard_model const&, int, int);
double partial_width_hh_to_barFdFd(standard_model::Standard_model const&, int, int);
double partial_width_hh_to_barFuFu(standard_model::Standard_model const&, int, int);
double partial_width_hh_to_barFeFe(standard_model::Standard_model const&, int, int);

private:
   standard_model::Standard_model model{};
   standard_model::Standard_model sm {};
   softsusy::QedQcd qedqcd{};
   Physical_input physical_input;
   FlexibleDecay_settings flexibledecay_settings {};
   bool run_to_decay_particle_scale {true};
   standard_model_decay_table decay_table{};
   FlexibleDecay_problems problems{};
   EffectiveCoupling_list neutral_higgs_effc {};

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   typename Decay_amplitude_type<FieldIn, FieldOut1, FieldOut2>::type
   calculate_amplitude_1l(
      const standard_model_cxx_diagrams::context_base&,
      const typename cxx_diagrams::field_indices<FieldIn>::type&,
      const typename cxx_diagrams::field_indices<FieldOut1>::type&,
      const typename cxx_diagrams::field_indices<FieldOut2>::type&) const;
   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   typename Decay_amplitude_type<FieldIn, FieldOut1, FieldOut2>::type
   calculate_amplitude_tree(
      const standard_model_cxx_diagrams::context_base&,
      const typename cxx_diagrams::field_indices<FieldIn>::type&,
      const typename cxx_diagrams::field_indices<FieldOut1>::type&,
      const typename cxx_diagrams::field_indices<FieldOut2>::type&) const;

   template<typename FieldIn, typename FieldOut1, typename FieldOut2>
   double amplitude_squared(standard_model_cxx_diagrams::context_base const& context,
                  typename cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
                  typename cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
                  typename cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2) const;

   std::vector<std::tuple<int, int, int, double, std::string>> effhiggscouplings_block_input {};
};

template<>
Decay_amplitude_SVV Standard_model_decays::calculate_amplitude_1l<standard_model_cxx_diagrams::fields::hh, standard_model_cxx_diagrams::fields::VG, standard_model_cxx_diagrams::fields::VG>(const standard_model_cxx_diagrams::context_base&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::hh >::type&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::VG >::type&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::VG >::type&) const;

template<>
Decay_amplitude_SVV Standard_model_decays::calculate_amplitude_1l<standard_model_cxx_diagrams::fields::hh, standard_model_cxx_diagrams::fields::VP, standard_model_cxx_diagrams::fields::VP>(const standard_model_cxx_diagrams::context_base&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::hh >::type&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::VP >::type&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::VP >::type&) const;

template<>
Decay_amplitude_SVV Standard_model_decays::calculate_amplitude_1l<standard_model_cxx_diagrams::fields::hh, standard_model_cxx_diagrams::fields::VP, standard_model_cxx_diagrams::fields::VZ>(const standard_model_cxx_diagrams::context_base&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::hh >::type&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::VP >::type&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV Standard_model_decays::calculate_amplitude_1l<standard_model_cxx_diagrams::fields::hh, standard_model_cxx_diagrams::fields::VZ, standard_model_cxx_diagrams::fields::VZ>(const standard_model_cxx_diagrams::context_base&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::hh >::type&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::VZ >::type&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::VZ >::type&) const;

template<>
Decay_amplitude_SVV Standard_model_decays::calculate_amplitude_1l<standard_model_cxx_diagrams::fields::hh, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::VWp>::type, standard_model_cxx_diagrams::fields::VWp>(const standard_model_cxx_diagrams::context_base&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::hh >::type&, const typename cxx_diagrams::field_indices<typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::VWp>::type >::type&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::VWp >::type&) const;

template<>
Decay_amplitude_SFF Standard_model_decays::calculate_amplitude_1l<standard_model_cxx_diagrams::fields::hh, typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fv>::type, standard_model_cxx_diagrams::fields::Fv>(const standard_model_cxx_diagrams::context_base&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::hh >::type&, const typename cxx_diagrams::field_indices<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fv>::type >::type&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::Fv >::type&) const;

template<>
Decay_amplitude_SFF Standard_model_decays::calculate_amplitude_1l<standard_model_cxx_diagrams::fields::hh, typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fd>::type, standard_model_cxx_diagrams::fields::Fd>(const standard_model_cxx_diagrams::context_base&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::hh >::type&, const typename cxx_diagrams::field_indices<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fd>::type >::type&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::Fd >::type&) const;

template<>
Decay_amplitude_SFF Standard_model_decays::calculate_amplitude_1l<standard_model_cxx_diagrams::fields::hh, typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fu>::type, standard_model_cxx_diagrams::fields::Fu>(const standard_model_cxx_diagrams::context_base&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::hh >::type&, const typename cxx_diagrams::field_indices<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fu>::type >::type&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::Fu >::type&) const;

template<typename FieldIn, typename FieldOut1, typename FieldOut2>
typename Decay_amplitude_type<FieldIn, FieldOut1, FieldOut2>::type
Standard_model_decays::calculate_amplitude_1l(
   const standard_model_cxx_diagrams::context_base& context,
   typename cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
   typename cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
   typename cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2) const
{
      // amplitude type
   typename Decay_amplitude_type<FieldIn, FieldOut1, FieldOut2>::type result;

   // external particles' masses
   result.m_decay = context.physical_mass<FieldIn>(indexIn);
   if constexpr (cxx_diagrams::fields::is_fermion<FieldOut1>::value && cxx_diagrams::fields::is_fermion<FieldOut2>::value) {
      result.m_fermion_1 = context.physical_mass<FieldOut1>(indexOut1);
      result.m_fermion_2 = context.physical_mass<FieldOut2>(indexOut2);
   }

   return result;
}

// generic decay of FieldIn -> FieldOut1 FieldOut2
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
double Standard_model_decays::get_partial_width(
   const standard_model_cxx_diagrams::context_base& context,
   typename cxx_diagrams::field_indices<FieldIn>::type const& indexIn,
   typename cxx_diagrams::field_indices<FieldOut1>::type const& indexOut1,
   typename cxx_diagrams::field_indices<FieldOut2>::type const& indexOut2
   )
{

   // on-shell masses
   const double mIn = context.physical_mass<FieldIn>(indexIn);
   const double mOut1 = context.physical_mass<FieldOut1>(indexOut1);
   const double mOut2 = context.physical_mass<FieldOut2>(indexOut2);

   // is decay kinematically allowed?
   if(mIn < mOut1 + mOut2) {
      problems.add_warning(
         create_process_string<FieldIn,FieldOut1,FieldOut2>(indexIn, indexOut1, indexOut2)
            + ": Called kinematically forbidden decay."
      );
      return 0.;
   }

   // phase space without symmetry factor
   const double ps = 1./(8.*Pi) * std::sqrt(KallenLambda(1., Sqr(mOut1/mIn), Sqr(mOut2/mIn)));

   // phase space symmetry factor
   const double ps_symmetry =
      final_state_symmetry_factor<FieldOut1, FieldOut2>(indexOut1, indexOut2);

   // color factor
   constexpr double color_factor = squared_color_generator<FieldIn, FieldOut1, FieldOut2>();

   typename Decay_amplitude_type<FieldIn, FieldOut1, FieldOut2>::type amp;

   // matrix element squared
   if constexpr (standard_model_has_tree_amp<FieldIn, FieldOut1, FieldOut2>::value) {
      amp = calculate_amplitude_tree<FieldIn, FieldOut1, FieldOut2>(
         context, indexIn, indexOut1, indexOut2);
   }
   else if constexpr (standard_model_has_oneloop_amp<FieldIn, FieldOut1, FieldOut2>::value) {
      amp = calculate_amplitude_1l<FieldIn, FieldOut1, FieldOut2>(
         context, indexIn, indexOut1, indexOut2);
   }

   // flux * phase space factor * symmetry factor * color factor * |matrix element|^2
   const auto result = 0.5/mIn * ps * ps_symmetry * color_factor * amp.square();

   return result;
}

template <>
double Standard_model_decays::get_partial_width<standard_model_cxx_diagrams::fields::hh,standard_model_cxx_diagrams::fields::VZ,standard_model_cxx_diagrams::fields::VZ >(const standard_model_cxx_diagrams::context_base&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::hh >::type&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::VZ >::type&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::VZ >::type&);
template <>
double Standard_model_decays::get_partial_width<standard_model_cxx_diagrams::fields::hh,typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::VWp>::type,standard_model_cxx_diagrams::fields::VWp >(const standard_model_cxx_diagrams::context_base&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::hh >::type&, const typename cxx_diagrams::field_indices<typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::VWp>::type >::type&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::VWp >::type&);
template <>
double Standard_model_decays::get_partial_width<standard_model_cxx_diagrams::fields::hh,standard_model_cxx_diagrams::fields::VG,standard_model_cxx_diagrams::fields::VG >(const standard_model_cxx_diagrams::context_base&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::hh >::type&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::VG >::type&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::VG >::type&);
template <>
double Standard_model_decays::get_partial_width<standard_model_cxx_diagrams::fields::hh,standard_model_cxx_diagrams::fields::VP,standard_model_cxx_diagrams::fields::VP >(const standard_model_cxx_diagrams::context_base&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::hh >::type&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::VP >::type&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::VP >::type&);
template <>
double Standard_model_decays::get_partial_width<standard_model_cxx_diagrams::fields::hh,standard_model_cxx_diagrams::fields::VP,standard_model_cxx_diagrams::fields::VZ >(const standard_model_cxx_diagrams::context_base&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::hh >::type&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::VP >::type&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::VZ >::type&);
template <>
double Standard_model_decays::get_partial_width<standard_model_cxx_diagrams::fields::hh,typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fu>::type,standard_model_cxx_diagrams::fields::Fu >(const standard_model_cxx_diagrams::context_base&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::hh >::type&, const typename cxx_diagrams::field_indices<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fu>::type >::type&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::Fu >::type&);
template <>
double Standard_model_decays::get_partial_width<standard_model_cxx_diagrams::fields::hh,typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fd>::type,standard_model_cxx_diagrams::fields::Fd >(const standard_model_cxx_diagrams::context_base&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::hh >::type&, const typename cxx_diagrams::field_indices<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fd>::type >::type&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::Fd >::type&);
template <>
double Standard_model_decays::get_partial_width<standard_model_cxx_diagrams::fields::hh,typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fe>::type,standard_model_cxx_diagrams::fields::Fe >(const standard_model_cxx_diagrams::context_base&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::hh >::type&, const typename cxx_diagrams::field_indices<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fe>::type >::type&, const typename cxx_diagrams::field_indices<standard_model_cxx_diagrams::fields::Fe >::type&);

} // namespace flexiblesusy

#endif
