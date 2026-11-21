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

#ifndef LOG_H
#define LOG_H

#include <cmath>
#include <complex>

namespace flexiblesusy {


/// returns log(1 + z) for complex z
template <typename T>
std::complex<T> log1p(const std::complex<T>& z) noexcept
{
   const std::complex<T> u = T(1) + z;

   if (std::real(u) == T(1) && std::imag(u) == T(0)) {
      return z;
   } else if (std::real(u) <= T(0)) {
      return std::log(u);
   }

   return std::log(u)*(z/(u - T(1)));
}


/**
 * returns log(z) for complex z
 * @param z complex argument
 * @note Points on the branch cut are treated differently from std::log(z):
 * Points with Re(z) < 0 and Im(z) == -0.0 are mapped to Im(z) == 0.0
 * @return log(z)
 */
template <typename T>
std::complex<T> pos_log(const std::complex<T>& z) noexcept
{
   if (std::imag(z) == T(0) && std::real(z) > T(0)) {
      return { std::log(std::real(z)), T(0) };
   } else if (std::imag(z) == T(0)) {
      return { std::log(-std::real(z)), 4*std::atan(T(1)) };
   }

   return std::log(z);
}


} // namespace flexiblesusy

#endif
