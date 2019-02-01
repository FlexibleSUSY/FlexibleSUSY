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

#ifndef DECAY_AMPLITUDES_H
#define DECAY_AMPLITUDES_H

#include "field_traits.hpp"

#include <complex>
#include <limits>

namespace flexiblesusy {

/**
 * @class Decay_amplitude_SSS
 * @brief generic amplitude for the decay of a scalar into two scalars
 */
struct Decay_amplitude_SSS {
   double m_decay{0.};
   double m_out_1{0.};
   double m_out_2{0.};
   std::complex<double> form_factor{};

   double square() const;
};

/**
 * @class Decay_amplitude_SSV
 * @brief generic amplitude for the decay of a scalar into a scalar and vector
 */
struct Decay_amplitude_SSV {
   double m_decay{0.};
   double m_scalar{0.};
   double m_vector{0.};
   double massless_vector_threshold{std::numeric_limits<double>::epsilon()};
   std::complex<double> form_factor{};

   double square() const;
};

/**
 * @class Decay_amplitude_SVV
 * @brief generic amplitude for the decay of a scalar into two vectors
 */
struct Decay_amplitude_SVV {
   double m_decay{0.};
   double m_out_1{0.};
   double m_out_2{0.};
   double massless_vector_threshold{std::numeric_limits<double>::epsilon()};
   std::complex<double> form_factor_g{};
   std::complex<double> form_factor_11{};
   std::complex<double> form_factor_12{};
   std::complex<double> form_factor_21{};
   std::complex<double> form_factor_22{};
   std::complex<double> form_factor_eps{};

   double square() const;
};

/**
 * @class Decay_amplitude_SFF
 * @brief generic amplitude for the decay of a scalar into two fermions
 */
struct Decay_amplitude_SFF {
   double m_decay{0.};
   double m_out_1{0.};
   double m_out_2{0.};
   std::complex<double> form_factor_left{};
   std::complex<double> form_factor_right{};

   double square() const;
};

/**
 * @class Decay_amplitude_FFS
 * @brief generic amplitude for the decay of a fermion into a fermion and scalar
 */
struct Decay_amplitude_FFS {
   double m_decay{0.};
   double m_fermion{0.};
   double m_scalar{0.};
   std::complex<double> form_factor_left{};
   std::complex<double> form_factor_right{};

   double square() const;
};

/**
 * @class Decay_amplitude_FFV
 * @brief generic amplitude for the decay of a fermion into a fermion and vector
 */
struct Decay_amplitude_FFV {
   double m_decay{0.};
   double m_fermion{0.};
   double m_vector{0.};
   double massless_vector_threshold{std::numeric_limits<double>::epsilon()};
   std::complex<double> form_factor_gam_left{};
   std::complex<double> form_factor_gam_right{};
   std::complex<double> form_factor_p_1{};
   std::complex<double> form_factor_p_2{};

   double square() const;
};

namespace detail {

template <class Field_in, class Field_out_1, class Field_out_2,
          class Amplitude_type = void>
struct Two_body_decay_amplitude_type { };

template <class Field_in, class Field_out_1, class Field_out_2>
struct Two_body_decay_amplitude_type<
   Field_in, Field_out_1, Field_out_2,
   typename std::enable_if<field_traits::is_scalar<Field_in>::value &&
                           field_traits::is_scalar<Field_out_1>::value &&
                           field_traits::is_scalar<Field_out_2>::value>::type > {
   using type = Decay_amplitude_SSS;
};

template <class Field_in, class Field_out_1, class Field_out_2>
struct Two_body_decay_amplitude_type<
   Field_in, Field_out_1, Field_out_2,
   typename std::enable_if<field_traits::is_scalar<Field_in>::value &&
                           field_traits::is_scalar<Field_out_1>::value &&
                           field_traits::is_vector<Field_out_2>::value>::type> {
   using type = Decay_amplitude_SSV;
};

template <class Field_in, class Field_out_1, class Field_out_2>
struct Two_body_decay_amplitude_type<
   Field_in, Field_out_1, Field_out_2,
   typename std::enable_if<field_traits::is_scalar<Field_in>::value &&
                           field_traits::is_vector<Field_out_1>::value &&
                           field_traits::is_scalar<Field_out_2>::value>::type> {
   using type = Decay_amplitude_SSV;
};

template <class Field_in, class Field_out_1, class Field_out_2>
struct Two_body_decay_amplitude_type<
   Field_in, Field_out_1, Field_out_2,
   typename std::enable_if<field_traits::is_scalar<Field_in>::value &&
                           field_traits::is_vector<Field_out_1>::value &&
                           field_traits::is_vector<Field_out_2>::value>::type> {
   using type = Decay_amplitude_SVV;
};

template <class Field_in, class Field_out_1, class Field_out_2>
struct Two_body_decay_amplitude_type<
   Field_in, Field_out_1, Field_out_2,
   typename std::enable_if<field_traits::is_scalar<Field_in>::value &&
                           field_traits::is_fermion<Field_out_1>::value &&
                           field_traits::is_fermion<Field_out_2>::value>::type> {
   using type = Decay_amplitude_SFF;
};

template <class Field_in, class Field_out_1, class Field_out_2>
struct Two_body_decay_amplitude_type<
   Field_in, Field_out_1, Field_out_2,
   typename std::enable_if<field_traits::is_fermion<Field_in>::value &&
                           field_traits::is_fermion<Field_out_1>::value &&
                           field_traits::is_scalar<Field_out_2>::value>::type> {
   using type = Decay_amplitude_FFS;
};

template <class Field_in, class Field_out_1, class Field_out_2>
struct Two_body_decay_amplitude_type<
   Field_in, Field_out_1, Field_out_2,
   typename std::enable_if<field_traits::is_fermion<Field_in>::value &&
                           field_traits::is_scalar<Field_out_1>::value &&
                           field_traits::is_fermion<Field_out_2>::value>::type> {
   using type = Decay_amplitude_FFS;
};

template <class Field_in, class Field_out_1, class Field_out_2>
struct Two_body_decay_amplitude_type<
   Field_in, Field_out_1, Field_out_2,
   typename std::enable_if<field_traits::is_fermion<Field_in>::value &&
                           field_traits::is_fermion<Field_out_1>::value &&
                           field_traits::is_vector<Field_out_2>::value>::type> {
   using type = Decay_amplitude_FFV;
};

template <class Field_in, class Field_out_1, class Field_out_2>
struct Two_body_decay_amplitude_type<
   Field_in, Field_out_1, Field_out_2,
   typename std::enable_if<field_traits::is_fermion<Field_in>::value &&
                           field_traits::is_vector<Field_out_1>::value &&
                           field_traits::is_fermion<Field_out_2>::value>::type> {
   using type = Decay_amplitude_FFV;
};

} // namespace detail

/**
 * @class Decay_amplitude_type
 * @brief helper class to determine amplitude type for a given set of fields
 */
template <class... Fields>
struct Decay_amplitude_type {
   using type =
      typename std::enable_if<
      sizeof...(Fields) == 3,
      typename detail::Two_body_decay_amplitude_type<Fields...>::type >::type;
};

template <typename Amplitude>
double square_amplitude(const Amplitude& a)
{
   return a.square();
}

} // namespace flexiblesusy

#endif