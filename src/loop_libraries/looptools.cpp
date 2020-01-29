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
#include "looptools.hpp"
#include "clooptools.h"

#define LOOPTOOLS_B(NAME,INDEX)\
   std::complex<double> Looptools::NAME(\
      std::complex<double> p10_in,\
      std::complex<double> m02_in, std::complex<double> m12_in,\
      double scl2_in) noexcept\
{\
   set_mu2_uv(scl2_in);\
   return B0i(INDEX, p10_in.real(), m02_in.real(), m12_in.real()); \
}
#define LOOPTOOLS_C(NAME,INDEX)\
   std::complex<double> Looptools::NAME(\
      std::complex<double> p10_in, std::complex<double> p21_in, std::complex<double> p20_in,\
      std::complex<double> m02_in, std::complex<double> m12_in, std::complex<double> m22_in,\
      double scl2_in) noexcept\
{\
   set_mu2_uv(scl2_in);\
   return C0i(INDEX, p10_in.real(), p21_in.real(), p20_in.real(), m02_in.real(), m12_in.real(), m22_in.real()); \
}
#define LOOPTOOLS_D(NAME,INDEX)\
   std::complex<double> Looptools::NAME(\
      std::complex<double> p10_in, std::complex<double> p21_in, std::complex<double> p32_in,\
      std::complex<double> p30_in, std::complex<double> p20_in, std::complex<double> p31_in,\
      std::complex<double> m02_in, std::complex<double> m12_in, std::complex<double> m22_in, std::complex<double> m32_in,\
      double scl2_in) noexcept\
{\
   set_mu2_uv(scl2_in);\
   return D0i(INDEX, p10_in.real(), p21_in.real(), p32_in.real(), p30_in.real(),\
                     p20_in.real(), p31_in.real(), m02_in.real(), m12_in.real(),\
                     m22_in.real(), m32_in.real());\
}

namespace looplibrary {

Looptools::Looptools() : current_mu2_uv(1.0)
{
   ltini();
}

void Looptools::set_mu2_uv(double scl2_in) noexcept
{
   if( std::abs(scl2_in - this->current_mu2_uv) > std::numeric_limits<double>::epsilon() )
   {
      setmudim(scl2_in);
      this->current_mu2_uv = scl2_in;
   }
}

std::complex<double> Looptools::A0(
   std::complex<double> m02_in,
   double scl2_in) noexcept
{
   set_mu2_uv(scl2_in);
   return A0i(aa0, m02_in.real());
}

LOOPTOOLS_B(B0,bb0)
LOOPTOOLS_B(B1,bb1)
LOOPTOOLS_B(B00,bb00)

LOOPTOOLS_C(C0,cc0)
LOOPTOOLS_C(C1,cc1)
LOOPTOOLS_C(C2,cc2)
LOOPTOOLS_C(C00,cc00)
LOOPTOOLS_C(C11,cc11)
LOOPTOOLS_C(C12,cc12)
LOOPTOOLS_C(C22,cc22)

LOOPTOOLS_D(D0,dd0)
LOOPTOOLS_D(D1,dd1)
LOOPTOOLS_D(D2,dd2)
LOOPTOOLS_D(D3,dd3)
LOOPTOOLS_D(D00,dd00)
LOOPTOOLS_D(D11,dd11)
LOOPTOOLS_D(D12,dd12)
LOOPTOOLS_D(D13,dd13)
LOOPTOOLS_D(D22,dd22)
LOOPTOOLS_D(D23,dd23)
LOOPTOOLS_D(D33,dd33)

void Looptools::A(
   std::complex<double> (&a)[1],
   std::complex<double> m02_in,
   double scl2_in) noexcept
{
   set_mu2_uv(scl2_in);
   a[0] = A0i(aa0, m02_in.real());
}

void Looptools::B(
   std::complex<double> (&b)[2],
   std::complex<double> p10_in,
   std::complex<double> m02_in, std::complex<double> m12_in,
   double scl2_in) noexcept
{
   double p10 = p10_in.real();
   double m02 = m02_in.real();
   double m12 = m12_in.real();
   set_mu2_uv(scl2_in);
   b[0] = B0i(bb0, p10, m02, m12);
   b[1] = B0i(bb1, p10, m02, m12);
}

void Looptools::C(
   std::complex<double> (&c)[7],
   std::complex<double> p10_in, std::complex<double> p21_in, std::complex<double> p20_in,
   std::complex<double> m02_in, std::complex<double> m12_in, std::complex<double> m22_in,
   double scl2_in) noexcept
{
   double p10 = p10_in.real();
   double p21 = p21_in.real();
   double p20 = p20_in.real();
   double m02 = m02_in.real();
   double m12 = m12_in.real();
   double m22 = m22_in.real();
   set_mu2_uv(scl2_in);
   c[0] = C0i(cc0, p10, p21, p20, m02, m12, m22);
   c[1] = C0i(cc1, p10, p21, p20, m02, m12, m22);
   c[2] = C0i(cc2, p10, p21, p20, m02, m12, m22);
   c[3] = C0i(cc00, p10, p21, p20, m02, m12, m22);
   c[4] = C0i(cc11, p10, p21, p20, m02, m12, m22);
   c[5] = C0i(cc12, p10, p21, p20, m02, m12, m22);
   c[6] = C0i(cc22, p10, p21, p20, m02, m12, m22);
}

void Looptools::D(
   std::complex<double> (&d)[11],
   std::complex<double> p10_in, std::complex<double> p21_in, std::complex<double> p32_in,
   std::complex<double> p30_in, std::complex<double> p20_in, std::complex<double> p31_in,
   std::complex<double> m02_in, std::complex<double> m12_in, std::complex<double> m22_in, std::complex<double> m32_in,
   double scl2_in) noexcept
{
   double p10 = p10_in.real();
   double p21 = p21_in.real();
   double p32 = p32_in.real();
   double p30 = p30_in.real();
   double p20 = p20_in.real();
   double p31 = p31_in.real();
   double m02 = m02_in.real();
   double m12 = m12_in.real();
   double m22 = m22_in.real();
   double m32 = m32_in.real();
   set_mu2_uv(scl2_in);
   d[0] = D0i(dd0, p10, p21, p32, p30, p20, p31, m02, m12, m22, m32);
   d[1] = D0i(dd1, p10, p21, p32, p30, p20, p31, m02, m12, m22, m32);
   d[2] = D0i(dd2, p10, p21, p32, p30, p20, p31, m02, m12, m22, m32);
   d[3] = D0i(dd3, p10, p21, p32, p30, p20, p31, m02, m12, m22, m32);
   d[4] = D0i(dd00, p10, p21, p32, p30, p20, p31, m02, m12, m22, m32);
   d[5] = D0i(dd11, p10, p21, p32, p30, p20, p31, m02, m12, m22, m32);
   d[6] = D0i(dd12, p10, p21, p32, p30, p20, p31, m02, m12, m22, m32);
   d[7] = D0i(dd13, p10, p21, p32, p30, p20, p31, m02, m12, m22, m32);
   d[8] = D0i(dd22, p10, p21, p32, p30, p20, p31, m02, m12, m22, m32);
   d[9] = D0i(dd23, p10, p21, p32, p30, p20, p31, m02, m12, m22, m32);
   d[10] = D0i(dd33, p10, p21, p32, p30, p20, p31, m02, m12, m22, m32);
}

} // namespace looplibrary
