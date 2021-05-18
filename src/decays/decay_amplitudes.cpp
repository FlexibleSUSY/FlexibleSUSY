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
 * @todo implement formulas for most general forms of amplitudes
 */

#include <iomanip>

#include "decay_amplitudes.hpp"
#include "wrappers.hpp"
#include "src/numerics2.hpp"

namespace flexiblesusy {

double Decay_amplitude_SSS::square() const
{
   return AbsSqr(form_factor);
}

double Decay_amplitude_SSV::square() const
{
   if (m_vector <= massless_vector_threshold) {
      // eq. B36 of http://etheses.dur.ac.uk/2301 has incorrect oposite sign
      return 2. * AbsSqr(form_factor) * (Sqr(m_decay) + Sqr(m_scalar));
   }

   const double m_in_sq = Sqr(m_decay);
   const double m_s_sq = Sqr(m_scalar);
   const double m_v_sq = Sqr(m_vector);

   return (Sqr(m_in_sq) + Sqr(m_v_sq - m_s_sq)
           - 2. * m_in_sq * (m_s_sq + m_v_sq))
      * AbsSqr(form_factor) / m_v_sq;
}

// @todo handle massless vectors safely
// @todo check expressions correct
double Decay_amplitude_SVV::square() const
{
   if (m_decay < m_vector_1 + m_vector_2) {
      WARNING("Warning in Decay_amplitude_SVV::square(): decaying particle mass smaller than sum of product masses. Returning 0.");
      return 0.;
   }

   // eq. B.38 of http://etheses.dur.ac.uk/2301
   if (m_vector_1 <= massless_vector_threshold &&
       m_vector_2 <= massless_vector_threshold) {

      const double fgSqr = AbsSqr(form_factor_g);
      const double f21Sqr = AbsSqr(form_factor_21);
      const double fepsSqr = AbsSqr(form_factor_eps);

      if (!is_zero(form_factor_21) && !is_zero(form_factor_g)) {
         // use Ward identity to eliminate form_factor_g
         const double res1 =
            0.5*Power4(m_decay)*(f21Sqr + fepsSqr);
         // use Ward identity to eliminate form_factor_21
         const double res2 =
            2.*fgSqr + 0.5*Power4(m_decay)*fepsSqr;
         const double WI_violation = std::abs(1. - std::abs(res1/res2));
         if (WI_violation > 0.1) {
            std::stringstream ss;
            ss << std::setprecision(2) << 100.*WI_violation;
            WARNING("Warning: Ward identity violated in decay of scalar to massless vectors by " + ss.str() + "%");
         }
         // use res1 since form_factor_21 is not sensitive to the renormalization
         // scheme in which the Higgs mass is defined
         return res1;
      }
      // use full expression for tree-level decays where one of the form factors might be 0
      else {
         const double Refgf21 = Re(form_factor_g * Conj(form_factor_21));
         const double res3 = 4.*fgSqr + Sqr(m_decay)*Refgf21 + 0.5*Power4(m_decay)*fepsSqr;
         return res3;
      }
   } else if (m_vector_1 <= massless_vector_threshold || m_vector_2 <= massless_vector_threshold) {

      const double m_s_sq = Sqr(m_decay);
      const double m_vec_sq = m_vector_1 <= massless_vector_threshold ? Sqr(m_vector_2) : Sqr(m_vector_1);

      const double fgSqr = AbsSqr(form_factor_g);
      const double f21Sqr = AbsSqr(form_factor_21);
      const double fepsSqr = AbsSqr(form_factor_eps);
      const double prefactor = 0.25*Sqr(m_s_sq - m_vec_sq);

      if (!is_zero(form_factor_21) && !is_zero(form_factor_g)) {
         // use Ward identity to eliminate form_factor_g
         const double res1 =
            2.*prefactor*(f21Sqr + fepsSqr);
         // use Ward identity to eliminate form_factor_21
         const double res2 =
            2.*(fgSqr + prefactor*fepsSqr);
         // compare two results
         const double WI_violation = std::abs(1. - std::abs(res1/res2));
         if (WI_violation > 0.1) {
            std::stringstream ss;
            ss << std::setprecision(2) << 100.*WI_violation;
            WARNING("Warning: Ward identity violated in decay of scalar to massless and massive vector by " + ss.str() + "%");
         }
         // use res1 since form_factor_21 is not sensitive to the renormalization
         // scheme in which the Higgs mass is defined
         return res1;
      }
      // use full expression for tree-level decays where one of the form factors might be 0
      else {
         const double res3 =
            3.*fgSqr + prefactor*(2.*fepsSqr-f21Sqr);
         return res3;
      }
   }

   const double m_s_sq = Sqr(m_decay);
   const double m_s_4 = Power4(m_decay);
   const double m_1_sq = Sqr(m_vector_1);
   const double m_1_4 = Power4(m_vector_1);
   const double m_2_sq = Sqr(m_vector_2);
   const double m_2_4 = Power4(m_vector_2);

   const double mgg = 0.25 * AbsSqr(form_factor_g) * (
      m_s_4 + m_1_4 + m_2_4 + 10. * m_1_sq * m_2_sq
      - 2. * m_s_sq * (m_1_sq + m_2_sq)) / (m_1_sq * m_2_sq);

   const double m33 = 0.0625 * Sqr(m_s_4 + Sqr(m_1_sq - m_2_sq)
                                   - 2. * m_s_sq * (m_1_sq + m_2_sq))
      * AbsSqr(form_factor_21) / (m_1_sq * m_2_sq);

   const double mepseps = AbsSqr(form_factor_eps) * (
      0.5 * Sqr(m_s_sq - m_1_sq - m_2_sq) - 2. * m_1_sq * m_2_sq);

   const double mg3 = 0.25 * (m_s_4 * m_s_sq - 3. * m_s_4 * (m_1_sq + m_2_sq)
                              - Sqr(m_1_sq - m_2_sq) * (m_1_sq + m_2_sq)
                              + m_s_sq * (3. * m_1_4 + 2. * m_1_sq * m_2_sq
                                          + 3. * m_2_4))
      * Re(form_factor_g * Conj(form_factor_21)) / (m_1_sq * m_2_sq);

   return mgg + m33 + mepseps + mg3;
}

Decay_amplitude_SSS operator*(std::complex<double> factor, Decay_amplitude_SSS const& amp2) {
   Decay_amplitude_SSS amp;
   amp.m_decay = amp2.m_decay;
   amp.m_scalar_1 = amp2.m_scalar_1;
   amp.m_scalar_2 = amp2.m_scalar_2;
   amp.form_factor = factor * amp2.form_factor;
   return amp;
}
Decay_amplitude_SSS operator*(Decay_amplitude_SSS const& amp2, std::complex<double> factor) {
   return operator*(factor, amp2);
}

Decay_amplitude_SSV operator*(std::complex<double> factor, Decay_amplitude_SSV const& amp2) {
   Decay_amplitude_SSV amp;
   amp.m_decay = amp2.m_decay;
   amp.m_scalar = amp2.m_scalar;
   amp.m_vector = amp2.m_vector;
   amp.form_factor = factor * amp2.form_factor;
   return amp;
}
Decay_amplitude_SSV operator*(Decay_amplitude_SSV const& amp2, std::complex<double> factor) {
   return operator*(factor, amp2);
}

Decay_amplitude_SVV operator* (std::complex<double> factor, Decay_amplitude_SVV const& amp2) {
   Decay_amplitude_SVV amp;
   amp.m_decay = amp2.m_decay;
   amp.m_vector_1 = amp2.m_vector_1;
   amp.m_vector_2 = amp2.m_vector_2;
   amp.massless_vector_threshold = amp2.massless_vector_threshold;
   amp.form_factor_g = factor * amp2.form_factor_g;
   amp.form_factor_11 = factor * amp2.form_factor_11;
   amp.form_factor_12 = factor * amp2.form_factor_12;
   amp.form_factor_21 = factor * amp2.form_factor_21;
   amp.form_factor_22 = factor * amp2.form_factor_22;
   amp.form_factor_eps = factor * amp2.form_factor_eps;
   return amp;
}
Decay_amplitude_SVV operator*(Decay_amplitude_SVV const& amp2, std::complex<double> factors) {
   return operator*(factors, amp2);
}

Decay_amplitude_SFF operator*(std::complex<double> factor, Decay_amplitude_SFF const& amp2) {
   Decay_amplitude_SFF amp;
   amp.m_decay = amp2.m_decay;
   amp.m_fermion_1 = amp2.m_fermion_1;
   amp.m_fermion_2 = amp2.m_fermion_2;
   amp.form_factor_left = factor * amp2.form_factor_left;
   amp.form_factor_right = factor * amp2.form_factor_right;
   return amp;
}
Decay_amplitude_SFF operator*(Decay_amplitude_SFF const& amp, std::complex<double> factor) {
   return operator*(factor, amp);
}

double Decay_amplitude_SFF::square() const
{
   const double m_in_sq = Sqr(m_decay);
   const double m_1_sq = Sqr(m_fermion_1);
   const double m_2_sq = Sqr(m_fermion_2);

   return (m_in_sq - m_1_sq - m_2_sq) *
      (AbsSqr(form_factor_left) + AbsSqr(form_factor_right))
      - 4. * m_fermion_1 * m_fermion_2 * Re(
         form_factor_left * Conj(form_factor_right));
}

double Decay_amplitude_FFS::square() const
{
   const double m_in_sq = Sqr(m_decay);
   const double m_f_sq = Sqr(m_fermion);
   const double m_s_sq = Sqr(m_scalar);

   return 0.5 * (m_in_sq + m_f_sq - m_s_sq) *
      (AbsSqr(form_factor_left) + AbsSqr(form_factor_right))
      + 2. * m_fermion * m_scalar * Re(
         form_factor_left * Conj(form_factor_right));
}

// This routine is for decays that we currently do not support and is not used
// Currently in development for future versions
// @todo handle massless vectors safely
// @todo check these expressions, they appear not to agree with SARAH/SPheno
double Decay_amplitude_FFV::square() const
{
   const double m_in_sq = Sqr(m_decay);
   const double m_f_sq = Sqr(m_fermion);

   if (m_vector <= massless_vector_threshold) {
      const double c1 = m_in_sq + m_f_sq;
      const double c2 = -0.5 * m_in_sq * (m_in_sq + m_f_sq);
      const double c3 = -4. * m_decay * m_fermion;
      const double c4 = -m_in_sq * m_fermion;
      const double c5 = -0.5 * m_decay * (m_in_sq + m_f_sq);
      const double c6 = -0.5 * m_decay * m_fermion * (m_in_sq + m_f_sq);

      return c1 * (AbsSqr(form_factor_gam_left) + AbsSqr(form_factor_gam_right))
         + c2 * (AbsSqr(form_factor_p_1) + AbsSqr(form_factor_p_2))
         + 2. * c3 * Re(form_factor_gam_left * Conj(form_factor_gam_right))
         + 2. * c4 * (Re(form_factor_gam_left * Conj(form_factor_p_1))
                      + Re(form_factor_gam_right * Conj(form_factor_p_2)))
         + 2. * c5 * (Re(form_factor_gam_left * Conj(form_factor_p_2))
                      + Re(form_factor_gam_right * Conj(form_factor_p_1)))
         + 2. * c6 * Re(form_factor_p_1 * Conj(form_factor_p_2));
   }

   const double m_in_4 = Power4(m_decay);
   const double m_f_4 = Power4(m_fermion);
   const double m_v_sq = Sqr(m_vector);
   const double m_v_4 = Power4(m_vector);

   const double c1 = 0.5 * (m_in_4 + m_f_4 + m_f_sq * m_v_sq - 2. * m_v_4
                            + m_in_sq * (m_v_sq - 2. * m_f_sq)) / m_v_sq;

   const double c2 = 0.125 * ((m_in_sq + m_f_sq - m_v_sq) * (
                                 m_in_4  + Sqr(m_f_sq - m_v_sq)
                                 - 2. * m_in_sq * (m_f_sq + m_v_sq))) / m_v_sq;

   const double c3 = -3. * m_decay * m_fermion;

   const double c4 = 0.25 * (m_fermion * (m_in_4 + Sqr(m_f_sq - m_v_sq)
                                          - 2. * m_in_sq * (m_f_sq + m_v_sq)))
      / m_v_sq;

   const double c5 = 0.25 * (m_decay * (m_in_4 + Sqr(m_f_sq - m_v_sq)
                                        - 2. * m_in_sq * (m_f_sq + m_v_sq)))
      / m_v_sq;

   const double c6 = 0.25 * (m_decay * m_fermion * (
                                m_in_4 + Sqr(m_f_sq - m_v_sq) - 2. * m_in_sq * (
                                   m_f_sq + m_v_sq))) / m_v_sq;

   return c1 * (AbsSqr(form_factor_gam_left) + AbsSqr(form_factor_gam_right))
      + c2 * (AbsSqr(form_factor_p_1) + AbsSqr(form_factor_p_2))
      + 2. * c3 * Re(form_factor_gam_left * Conj(form_factor_gam_right))
      + 2. * c4 * Re(form_factor_gam_left * Conj(form_factor_p_1)
                     + form_factor_gam_right * Conj(form_factor_p_2))
      + 2. * c5 * Re(form_factor_gam_left * Conj(form_factor_p_2)
                     + form_factor_gam_right * Conj(form_factor_p_1))
      + 2. * c6 * Re(form_factor_p_1 * Conj(form_factor_p_2));
}

} // namespace flexiblesusy
