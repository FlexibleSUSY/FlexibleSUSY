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

#ifndef DECAY_AMPLITUDES_SRC_H
#define DECAY_AMPLITUDES_SRC_H

#include <complex>
#include <limits>

namespace flexiblesusy {

/**
 * @class Decay_amplitude_SSS
 * @brief generic amplitude for the decay of a scalar into two scalars
 */
struct Decay_amplitude_SSS {
   double m_decay {0.};
   double m_scalar_1 {0.};
   double m_scalar_2 {0.};
   std::complex<double> form_factor {0.};

   double square() const;
   Decay_amplitude_SSS operator+=(Decay_amplitude_SSS const& amp) {
      form_factor += amp.form_factor;
      return *this;
   }
};
Decay_amplitude_SSS operator*(std::complex<double>, Decay_amplitude_SSS const&);


/**
 * @class Decay_amplitude_SSV
 * @brief generic amplitude for the decay of a scalar into a scalar and vector
 */
struct Decay_amplitude_SSV {
   double m_decay {0.};
   double m_scalar {0.};
   double m_vector {0.};
   double massless_vector_threshold{std::numeric_limits<double>::epsilon()};
   std::complex<double> form_factor {0.};

   double square() const;
   Decay_amplitude_SSV operator+=(Decay_amplitude_SSV const& amp) {
      form_factor += amp.form_factor;
      return *this;
   }
};
Decay_amplitude_SSV operator*(std::complex<double>, Decay_amplitude_SSV const&);

/**
 * @class Decay_amplitude_SVV
 * @brief generic amplitude for the decay of a scalar into two vectors
 */
struct Decay_amplitude_SVV {
   double m_decay {0.};
   double m_vector_1 {0.};
   double m_vector_2 {0.};
   double massless_vector_threshold{std::numeric_limits<double>::epsilon()};
   std::complex<double> form_factor_g {0.};
   std::complex<double> form_factor_11 {0.};
   std::complex<double> form_factor_12 {0.};
   std::complex<double> form_factor_21 {0.};
   std::complex<double> form_factor_22 {0.};
   std::complex<double> form_factor_eps {0.};

   double square() const;
   Decay_amplitude_SVV operator+=(Decay_amplitude_SVV const& amp) {
      form_factor_g += amp.form_factor_g;
      form_factor_11 += amp.form_factor_11;
      form_factor_12 += amp.form_factor_12;
      form_factor_21 += amp.form_factor_21;
      form_factor_22 += amp.form_factor_22;
      form_factor_eps += amp.form_factor_eps;
      return *this;
   }
};
Decay_amplitude_SVV operator*(std::complex<double>, Decay_amplitude_SVV const&);

/**
 * @class Decay_amplitude_SFF
 * @brief generic amplitude for the decay of a scalar into two fermions
 */
struct Decay_amplitude_SFF {
   double m_decay {0.};
   double m_fermion_1 {0.};
   double m_fermion_2 {0.};
   std::complex<double> form_factor_left {0.};
   std::complex<double> form_factor_right {0.};

   double square() const;
   Decay_amplitude_SFF operator+=(Decay_amplitude_SFF const& amp) {
      form_factor_left += amp.form_factor_left;
      form_factor_right += amp.form_factor_right;
      return *this;
   }
};
Decay_amplitude_SFF operator*(std::complex<double>, Decay_amplitude_SFF const&);

/**
 * @class Decay_amplitude_FFS
 * @brief generic amplitude for the decay of a fermion into a fermion and scalar
 */
struct Decay_amplitude_FFS {
   double m_decay {0.};
   double m_fermion {0.};
   double m_scalar {0.};
   std::complex<double> form_factor_left {0.};
   std::complex<double> form_factor_right {0.};

   double square() const;
   Decay_amplitude_FFS operator+=(Decay_amplitude_FFS const& amp) {
      form_factor_left += amp.form_factor_left;
      form_factor_right += amp.form_factor_right;
      return *this;
   }
};
Decay_amplitude_FFS operator*(std::complex<double>, Decay_amplitude_FFS const&);

/**
 * @class Decay_amplitude_FFV
 * @brief generic amplitude for the decay of a fermion into a fermion and vector
 */
struct Decay_amplitude_FFV {
   double m_decay {0.};
   double m_fermion {0.};
   double m_vector {0.};
   double massless_vector_threshold{std::numeric_limits<double>::epsilon()};
   std::complex<double> form_factor_gam_left {0.};
   std::complex<double> form_factor_gam_right {0.};
   std::complex<double> form_factor_p_1 {0.};
   std::complex<double> form_factor_p_2 {0.};

   double square() const;
   Decay_amplitude_FFV operator+=(Decay_amplitude_FFV const& amp) {
      form_factor_gam_left += amp.form_factor_gam_left;
      form_factor_gam_right += amp.form_factor_gam_right;
      form_factor_p_1 += amp.form_factor_p_1;
      form_factor_p_2 += amp.form_factor_p_2;
      return *this;
   }
};
Decay_amplitude_FFV operator*(std::complex<double>, Decay_amplitude_FFV const&);

template <typename Amplitude>
double square_amplitude(const Amplitude& a)
{
   return a.square();
}

} // namespace flexiblesusy

#endif
