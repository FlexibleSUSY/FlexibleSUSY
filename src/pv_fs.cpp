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

#include "pv_fs.hpp"

#include <cmath>
#include <complex>

/**
 * @file PV.cpp
 *
 * @brief Implementation of real Passarino-Veltman loop functions with
 * squared arguments.
 */

namespace flexiblesusy {
namespace {

constexpr double EPSTOL = 1.0e-15; ///< underflow accuracy


constexpr double pow2(double x) noexcept { return x * x; }
constexpr double pow3(double x) noexcept { return x * x * x; }
constexpr double pow4(double x) noexcept { return pow2(pow2(x)); }
constexpr double pow5(double x) noexcept { return pow4(x) * x; }
constexpr double pow6(double x) noexcept { return pow2(pow2(x) * x); }


constexpr double sign(double x) noexcept
{
   return x >= 0.0 ? 1.0 : -1.0;
}


/// fast implementation of complex logarithm
template <class T>
std::complex<T> fast_log(const std::complex<T>& z) noexcept
{
   const T rz = std::real(z);
   const T iz = std::imag(z);

   return std::complex<T>(0.5*std::log(rz*rz + iz*iz), std::atan2(iz, rz));
}


double fB(const std::complex<double>& x) noexcept
{
   const double re = std::real(x);
   const double im = std::imag(x);

   if ((std::abs(re) == 0.0 || std::abs(re) == 1.0) && im == 0.0) {
      return -1.0;
   }

   if (std::abs(re - 1.0) < 1e-3) {
      const auto d = x - 1.0;
      const auto logd = fast_log(d);
      return std::real(-1.0 + d*(1.0 - logd + d*(0.5 - d/6.)));
   }

   if (std::abs(re) > 1e3) {
      return std::real(-0.5/x - 1.0/(6.0*x*x) + fast_log(-x));
   }

   return std::real(-1.0 + fast_log(1.0 - x) - x*fast_log(1.0 - 1.0/x));
}


} // anonymous namespace


double a0(double m2, double q2) noexcept
{
   m2 = std::abs(m2);

   if (m2 < EPSTOL)
      return 0;

   return m2 * (1 - std::log(m2 / q2));
}


/**
 * Re(B0(s,x,x,q2)), Eq.(2.4) from [hep-ph/0701051]
 *
 * @param p2 squared momentum
 * @param m2 squared mass
 * @param q2 squared renormalization scale
 *
 * @return Re(B0(s,x,x,q2))
 */
double b0xx(double p2, double m2, double q2) noexcept
{
   p2 = std::abs(p2);
   m2 = std::abs(m2);
   q2 = std::abs(q2);

   if (p2 < EPSTOL * q2 && m2 < EPSTOL * q2) {
      return 0; // IR divergence
   } else if (p2 < EPSTOL) {
      return -std::log(m2 / q2);
   } else if (p2 <= 4 * m2) {
      return 2 - std::log(m2 / q2) -
             2 * std::sqrt(4 * m2 / p2 - 1) *
                std::asin(std::sqrt(p2 / (4 * m2)));
   } else if (p2 <= 5e2 * m2) {
      const double sq = std::sqrt(1 - 4 * m2 / p2);
      return 2 - std::log(m2 / q2) +
         sq * std::log(p2 * (1 - sq) / (2 * m2) - 1);
   } else if (p2*EPSTOL < m2) {
      const double d = m2 / p2;
      const double logd = std::log(d);
      return 2 - std::log(p2 / q2)
         + d * (2 * (1 - logd)
         + d * (-1 - 2 * logd
         + d * (-10./3 - 4 * logd
         + d * (-59./6 - 10 * logd
         + d * (-449./15 - 28 * logd
         + d * (-1417./15 - 84 * logd
         + d * (-32254./105 - 264 * logd)))))));
   } else {
      return 2 - std::log(p2 / q2);
   }
}


/**
 * B0 function with squared arguments, from hep-ph/9606211.
 * @note returns only the real part of B0
 */
double b0(double p2, double m12, double m22, double q2) noexcept
{
   p2  = std::abs(p2);
   m12 = std::abs(m12);
   m22 = std::abs(m22);
   q2  = std::abs(q2);

   // protect against infrared divergence
   if (p2 < EPSTOL * q2 && m12 < EPSTOL * q2 && m22 < EPSTOL * q2) {
      return 0;
   }

   if (m12 > m22) {
      std::swap(m12, m22);
   }

   if (m12*(1 + EPSTOL) > m22) {
      return b0xx(p2, m12, q2);
   }

   if (p2 <= 1e-11 * m12) {
      if (m12 < EPSTOL * m22) {
         return 1 - std::log(m22 / q2);
      }

      return 1 - std::log(m22 / q2) + m12 * std::log(m22 / m12) / (m12 - m22);
   }

   if (m12 < EPSTOL * EPSTOL * m22) {
      const double del = std::hypot(m22 - p2, EPSTOL * m22);
      return 2 - std::log(m22 / q2) + (m22 - p2) / p2 * std::log(del / p2);
   }

   const double s = p2 - m22 + m12;
   const std::complex<double> imin(m12, -EPSTOL * m12);
   const std::complex<double> x = std::sqrt(pow2(s) - 4 * p2 * imin);
   const std::complex<double> xp = (s + sign(s) * x) / (2 * p2);
   const std::complex<double> xm = imin / (xp * p2);

   return -std::log(p2 / q2) - fB(xp) - fB(xm);
}

/**
 * Derivative of B0(p^2, m1^2, m2^2, Q^2) w.r.t. p^2, for p^2 = 0.
 *
 * @note Implemented only in the p^2 = 0 limit.
 *
 * @param m12 squared mass
 * @param m22 squared mass
 *
 * @return derivative of B0 w.r.t. p^2 at p^2 = 0
 */
double d1_b0(double m12, double m22) noexcept
{
   m12 = std::abs(m12);
   m22 = std::abs(m22);

   if ((m12 < 1e-14*m22) != (m22 < 1e-14*m12)) {
      return (m12 - m22) * (m12 + m22) / (2 * pow3(m12 - m22));
   } else if (m12 < EPSTOL && m22 < EPSTOL) {
      return 0;
   } else if (std::abs(m22 - m12) < 5e-2*m12) {
      return 1 / (6 * m12) +
             (m12 - m22)     / (12 * pow2(m12)) +
             pow2(m12 - m22) / (20 * pow3(m12)) +
             pow3(m12 - m22) / (30 * pow4(m12)) +
             pow4(m12 - m22) / (42 * pow5(m12)) +
             pow5(m12 - m22) / (56 * pow6(m12));
   }

   return ((m12 - m22) * (m12 + m22) + 2 * m12 * m22 * std::log(m22 / m12)) /
          (2 * pow3(m12 - m22));
}

} // namespace flexiblesusy
