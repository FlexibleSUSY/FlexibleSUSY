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

#include "collier.hpp"
#include <limits>

#define two_point_impl(NAME)\
   std::complex<double> NAME##_impl(\
      const std::complex<double>*,\
      const std::complex<double>*, const std::complex<double>*);
#define three_point_impl(NAME)\
   std::complex<double> NAME##_impl(\
      const std::complex<double>*, const std::complex<double>*, const std::complex<double>*,\
      const std::complex<double>*, const std::complex<double>*, const std::complex<double>*);
#define four_point_impl(NAME)\
   std::complex<double> NAME##_impl(\
      const std::complex<double>*, const std::complex<double>*, const std::complex<double>*,\
      const std::complex<double>*, const std::complex<double>*, const std::complex<double>*,\
      const std::complex<double>*, const std::complex<double>*, const std::complex<double>*, const std::complex<double>*);

/* Non-vanishing imaginary parts of momentum invariants are not yet
* suppoted by the current version (1.2.4) of COLLIER. */
#define two_point_collier(NAME)\
   std::complex<double> Collier::NAME(\
      std::complex<double> p10_in,\
      std::complex<double> m02_in, std::complex<double> m12_in,\
      double scl2_in) noexcept\
{\
   const std::complex<double> p10 (p10_in.real(), 0.);\
   const std::complex<double> m02 = m02_in;\
   const std::complex<double> m12 = m12_in;\
\
   set_mu2_uv(scl2_in);\
   return NAME##_impl(&p10, &m02, &m12);\
}
#define three_point_collier(NAME)\
   std::complex<double> Collier::NAME(\
      std::complex<double> p10_in, std::complex<double> p21_in, std::complex<double> p20_in,\
      std::complex<double> m02_in, std::complex<double> m12_in, std::complex<double> m22_in,\
      double scl2_in) noexcept\
{\
   const std::complex<double> p10 (p10_in.real(), 0.);\
   const std::complex<double> p21 (p21_in.real(), 0.);\
   const std::complex<double> p20 (p20_in.real(), 0.);\
   const std::complex<double> m02 = m02_in;\
   const std::complex<double> m12 = m12_in;\
   const std::complex<double> m22 = m22_in;\
\
   set_mu2_uv(scl2_in);\
   return NAME##_impl(&p10, &p21, &p20, &m02, &m12, &m22);\
}
#define four_point_collier(NAME)\
   std::complex<double> Collier::NAME(\
      std::complex<double> p10_in, std::complex<double> p21_in, std::complex<double> p32_in,\
      std::complex<double> p30_in, std::complex<double> p20_in, std::complex<double> p31_in,\
      std::complex<double> m02_in, std::complex<double> m12_in, std::complex<double> m22_in, std::complex<double> m32_in,\
      double scl2_in) noexcept\
{\
   const std::complex<double> p10 (p10_in.real(), 0.);\
   const std::complex<double> p21 (p21_in.real(), 0.);\
   const std::complex<double> p32 (p32_in.real(), 0.);\
   const std::complex<double> p30 (p30_in.real(), 0.);\
   const std::complex<double> p20 (p20_in.real(), 0.);\
   const std::complex<double> p31 (p31_in.real(), 0.);\
   const std::complex<double> m02 = m02_in;\
   const std::complex<double> m12 = m12_in;\
   const std::complex<double> m22 = m22_in;\
   const std::complex<double> m32 = m32_in;\
\
   set_mu2_uv(scl2_in);\
   return NAME##_impl(&p10, &p21, &p32, &p30, &p20, &p31, &m02, &m12, &m22, &m32);\
}

// Fortran wrapper routines
extern "C" {
   void initialize_collier_impl();
   void set_mu2_uv_impl(double*);

   two_point_impl(B0)
   two_point_impl(B1)

   three_point_impl(C0)
   three_point_impl(C1)
   three_point_impl(C2)
   three_point_impl(C00)
   three_point_impl(C11)
   three_point_impl(C12)
   three_point_impl(C22)

   four_point_impl(D0)
   four_point_impl(D00)
   four_point_impl(D1)
   four_point_impl(D11)
   four_point_impl(D12)
   four_point_impl(D13)
   four_point_impl(D2)
   four_point_impl(D22)
   four_point_impl(D23)
   four_point_impl(D3)
   four_point_impl(D33)
}

namespace flexiblesusy {

void Collier::initialize() noexcept
{
   initialize_collier_impl();
}

void Collier::set_mu2_uv(double scl2_in) noexcept
{
   double scl2 = scl2_in;
   if( std::abs(scl2 - this->current_mu2_uv) > std::numeric_limits<double>::epsilon() )
   {
      set_mu2_uv_impl(&scl2);
      this->current_mu2_uv = scl2;
   }
}

two_point_collier(B0)
two_point_collier(B1)

three_point_collier(C0)
three_point_collier(C1)
three_point_collier(C2)
three_point_collier(C00)
three_point_collier(C11)
three_point_collier(C12)
three_point_collier(C22)

four_point_collier(D0)
four_point_collier(D00)
four_point_collier(D1)
four_point_collier(D11)
four_point_collier(D12)
four_point_collier(D13)
four_point_collier(D2)
four_point_collier(D22)
four_point_collier(D23)
four_point_collier(D3)
four_point_collier(D33)

}
