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

#ifndef HORNER_H
#define HORNER_H

#include <complex>

namespace flexiblesusy {

template <typename T, int N>
T horner(T x, const T (&c)[N]) noexcept
{
   T p = c[N - 1];
   for (int i = N - 2; i >= 0; --i) {
      p = p*x + c[i];
   }
   return p;
}


template <int Nstart, typename T, int N>
std::complex<T> horner(const std::complex<T>& z, const T (&coeffs)[N]) noexcept
{
   static_assert(0 <= Nstart && Nstart < N && N >= 2, "invalid array bounds");

   const T rz = std::real(z);
   const T iz = std::imag(z);
   const T r = rz + rz;
   const T s = std::norm(z);
   T a = coeffs[N - 1], b = coeffs[N - 2];

   for (int i = N - 3; i >= Nstart; --i) {
      const T t = a;
      a = b + r*a;
      b = coeffs[i] - s*t;
   }

   return { rz*a + b, iz*a };
}

} // namespace flexiblesusy

#endif
