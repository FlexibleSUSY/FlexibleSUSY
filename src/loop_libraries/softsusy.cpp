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

#include <limits>
#include "softsusy.hpp"
#include "numerics.h"

#define UNDEFINED_C(NAME)\
   std::complex<double> Softsusy::NAME(\
      std::complex<double> p10_in, std::complex<double> p21_in, std::complex<double> p20_in,\
      std::complex<double> m02_in, std::complex<double> m12_in, std::complex<double> m22_in,\
      double scl2_in) noexcept\
{\
   return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()}; \
}
#define UNDEFINED_D(NAME)\
   std::complex<double> Softsusy::NAME(\
      std::complex<double> p10_in, std::complex<double> p21_in, std::complex<double> p32_in,\
      std::complex<double> p30_in, std::complex<double> p20_in, std::complex<double> p31_in,\
      std::complex<double> m02_in, std::complex<double> m12_in, std::complex<double> m22_in, std::complex<double> m32_in,\
      double scl2_in) noexcept\
{\
   return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()}; \
}

namespace looplibrary {

std::complex<double> Softsusy::A0(
   std::complex<double> m02_in,
   double scl2_in) noexcept
{
   double m = std::sqrt(m02_in.real());
   double q = std::sqrt(scl2_in);

   return {softsusy::a0(m, q), 0.0};
}

std::complex<double> Softsusy::B0(
   std::complex<double> p10_in,
   std::complex<double> m02_in, std::complex<double> m12_in,
   double scl2_in) noexcept
{
   double p = std::sqrt(p10_in.real());
   double m1 = std::sqrt(m02_in.real());
   double m2 = std::sqrt(m12_in.real());
   double q = std::sqrt(scl2_in);

   return {softsusy::b0(p, m1, m2, q), 0.0};
}

std::complex<double> Softsusy::B1(
         std::complex<double> p10_in,
         std::complex<double> m02_in, std::complex<double> m12_in,
   double scl2_in) noexcept
{
   double p = std::sqrt(p10_in.real());
   double m1 = std::sqrt(m02_in.real());
   double m2 = std::sqrt(m12_in.real());
   double q = std::sqrt(scl2_in);

   return {(-1)*softsusy::b1(p, m1, m2, q), 0.0};
}

std::complex<double> Softsusy::C0(
   std::complex<double> p10_in, std::complex<double> p21_in, std::complex<double> p20_in,
   std::complex<double> m02_in, std::complex<double> m12_in, std::complex<double> m22_in,
   double scl2_in) noexcept
{
   double m1 = std::sqrt(m02_in.real());
   double m2 = std::sqrt(m12_in.real());
   double m3 = std::sqrt(m22_in.real());

   return {softsusy::c0(m1, m2, m3), 0.0};
}

std::complex<double> Softsusy::C00(
   std::complex<double> p10_in, std::complex<double> p21_in, std::complex<double> p20_in,
   std::complex<double> m02_in, std::complex<double> m12_in, std::complex<double> m22_in,
   double scl2_in) noexcept
{
   double m1 = std::sqrt(m02_in.real());
   double m2 = std::sqrt(m12_in.real());
   double m3 = std::sqrt(m22_in.real());
   double q = std::sqrt(scl2_in);

   return {softsusy::c00(m1, m2, m3, q), 0.0};
}

UNDEFINED_C(C1)
UNDEFINED_C(C2)
UNDEFINED_C(C11)
UNDEFINED_C(C12)
UNDEFINED_C(C22)

std::complex<double> Softsusy::D0(
      std::complex<double> p10_in, std::complex<double> p21_in, std::complex<double> p32_in,
      std::complex<double> p30_in, std::complex<double> p20_in, std::complex<double> p31_in,
      std::complex<double> m02_in, std::complex<double> m12_in, std::complex<double> m22_in, std::complex<double> m32_in,
   double scl2_in) noexcept
{
   double m1 = std::sqrt(m02_in.real());
   double m2 = std::sqrt(m12_in.real());
   double m3 = std::sqrt(m22_in.real());
   double m4 = std::sqrt(m32_in.real());

   return {softsusy::d0(m1, m2, m3, m4), 0.0};
}

std::complex<double> Softsusy::D00(
      std::complex<double> p10_in, std::complex<double> p21_in, std::complex<double> p32_in,
      std::complex<double> p30_in, std::complex<double> p20_in, std::complex<double> p31_in,
      std::complex<double> m02_in, std::complex<double> m12_in, std::complex<double> m22_in, std::complex<double> m32_in,
      double scl2_in) noexcept
{
   double m1 = std::sqrt(m02_in.real());
   double m2 = std::sqrt(m12_in.real());
   double m3 = std::sqrt(m22_in.real());
   double m4 = std::sqrt(m32_in.real());

   return {softsusy::d27(m1, m2, m3, m4), 0.0};
}

UNDEFINED_D(D1)
UNDEFINED_D(D11)
UNDEFINED_D(D12)
UNDEFINED_D(D13)
UNDEFINED_D(D2)
UNDEFINED_D(D22)
UNDEFINED_D(D23)
UNDEFINED_D(D3)
UNDEFINED_D(D33)

void Softsusy::get_A(
   std::complex<double> (&a)[1],
   std::complex<double> m02_in,
   double scl2_in) noexcept
{
   double m = std::sqrt(m02_in.real());
   double q = std::sqrt(scl2_in);

   a[0] = {softsusy::a0(m, q), 0.0};
}

void Softsusy::get_B(
   std::complex<double> (&b)[2],
   std::complex<double> p10_in,
   std::complex<double> m02_in, std::complex<double> m12_in,
   double scl2_in) noexcept
{
   double p = std::sqrt(p10_in.real());
   double m1 = std::sqrt(m02_in.real());
   double m2 = std::sqrt(m12_in.real());
   double q = std::sqrt(scl2_in);

   b[0] = {softsusy::b0(p, m1, m2, q), 0.0};
   b[1] = {(-1)*softsusy::b1(p, m1, m2, q), 0.0};
}

void Softsusy::get_C(
   std::complex<double> (&c)[7],
   std::complex<double> p10_in, std::complex<double> p21_in, std::complex<double> p20_in,
   std::complex<double> m02_in, std::complex<double> m12_in, std::complex<double> m22_in,
   double scl2_in) noexcept
{
   double m1 = std::sqrt(m02_in.real());
   double m2 = std::sqrt(m12_in.real());
   double m3 = std::sqrt(m22_in.real());
   double q = std::sqrt(scl2_in);
   std::complex<double> undefined = {std::numeric_limits<double>::quiet_NaN(),
                                   std::numeric_limits<double>::quiet_NaN()};

   c[0] = {softsusy::c0(m1, m2, m3), 0.0};
   c[1] = undefined;
   c[2] = undefined;
   c[3] = {softsusy::c00(m1, m2, m3, q), 0.0};
   c[4] = undefined;
   c[5] = undefined;
   c[6] = undefined;
}

void Softsusy::get_D(
   std::complex<double> (&d)[11],
   std::complex<double> p10_in, std::complex<double> p21_in, std::complex<double> p32_in,
   std::complex<double> p30_in, std::complex<double> p20_in, std::complex<double> p31_in,
   std::complex<double> m02_in, std::complex<double> m12_in, std::complex<double> m22_in, std::complex<double> m32_in,
   double scl2_in) noexcept
{
   double m1 = std::sqrt(m02_in.real());
   double m2 = std::sqrt(m12_in.real());
   double m3 = std::sqrt(m22_in.real());
   double m4 = std::sqrt(m32_in.real());
   std::complex<double> undefined = {std::numeric_limits<double>::quiet_NaN(),
                                   std::numeric_limits<double>::quiet_NaN()};

   d[0] = {softsusy::d0(m1, m2, m3, m4), 0.0};
   d[1] = undefined;
   d[2] = undefined;
   d[3] = undefined;
   d[4] = {softsusy::d27(m1, m2, m3, m4), 0.0};
   d[5] = undefined;
   d[6] = undefined;
   d[7] = undefined;
   d[8] = undefined;
   d[9] = undefined;
   d[10] = undefined;
}

} // namespace looplibrary
