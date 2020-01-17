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

#include "softsusy.hpp"

#define two_point_lt(NAME)\
   std::complex<double> Softsusy::NAME(\
      std::complex<double> p10_in,\
      std::complex<double> m02_in, std::complex<double> m12_in,\
      double scl2_in) noexcept\
{\
   return 0.; \
} // @ToDo
#define three_point_lt(NAME)\
   std::complex<double> Softsusy::NAME(\
      std::complex<double> p10_in, std::complex<double> p21_in, std::complex<double> p20_in,\
      std::complex<double> m02_in, std::complex<double> m12_in, std::complex<double> m22_in,\
      double scl2_in) noexcept\
{\
   return 0.; \
} // @ToDo
#define four_point_lt(NAME)\
   std::complex<double> Softsusy::NAME(\
      std::complex<double> p10_in, std::complex<double> p21_in, std::complex<double> p32_in,\
      std::complex<double> p30_in, std::complex<double> p20_in, std::complex<double> p31_in,\
      std::complex<double> m02_in, std::complex<double> m12_in, std::complex<double> m22_in, std::complex<double> m32_in,\
      double scl2_in) noexcept\
{\
   return 0.; \
} // @ToDo

namespace flexiblesusy {

two_point_lt(B0)
two_point_lt(B1)

three_point_lt(C0)
three_point_lt(C1)
three_point_lt(C2)
three_point_lt(C00)
three_point_lt(C11)
three_point_lt(C12)
three_point_lt(C22)

four_point_lt(D0)
four_point_lt(D00)
four_point_lt(D1)
four_point_lt(D11)
four_point_lt(D12)
four_point_lt(D13)
four_point_lt(D2)
four_point_lt(D22)
four_point_lt(D23)
four_point_lt(D3)
four_point_lt(D33)

void Softsusy::get_B(
   std::complex<double> (&b)[2],
   std::complex<double> p10_in,
   std::complex<double> m02_in, std::complex<double> m12_in,
   double scl2_in) noexcept
{
   // @ToDo
   for (int i = 0; i < 2; i++) {
      b[i] = 0.;
   }
}

void Softsusy::get_C(
   std::complex<double> (&c)[7],
   std::complex<double> p10_in, std::complex<double> p21_in, std::complex<double> p20_in,
   std::complex<double> m02_in, std::complex<double> m12_in, std::complex<double> m22_in,
   double scl2_in) noexcept
{
   // @ToDo
   for (int i = 0; i < 7; i++) {
      c[i] = 0.;
   }
}

void Softsusy::get_D(
   std::complex<double> (&d)[11],
   std::complex<double> p10_in, std::complex<double> p21_in, std::complex<double> p32_in,
   std::complex<double> p30_in, std::complex<double> p20_in, std::complex<double> p31_in,
   std::complex<double> m02_in, std::complex<double> m12_in, std::complex<double> m22_in, std::complex<double> m32_in,
   double scl2_in) noexcept
{
   // @ToDo
   for (int i = 0; i < 11; i++) {
      d[i] = 0.;
   }
}

}
