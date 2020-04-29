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
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include "fortran_utils.hpp"
#include "library_collier.hpp"

#define COLLIER_TYPE(Z, N, TEXT) const std::complex<double>*,
#define COLLIER_ARGS(N) BOOST_PP_REPEAT(N,COLLIER_TYPE,) const std::complex<double>*

#define IMPL(R,ARGS,NAME) std::complex<double> BOOST_PP_CAT(NAME,_impl)(COLLIER_ARGS(ARGS));

/* Non-vanishing imaginary parts of momentum invariants are not yet
* suppoted by the current version (1.2.4) of COLLIER. */
#define COLLIER_B(R,ARGS,NAME) std::complex<double> Collier::NAME ARGS noexcept\
{\
   const std::complex<double> p10 (p10_in.real(), 0.);\
   const std::complex<double> m02 = m02_in;\
   const std::complex<double> m12 = m12_in;\
\
   set_mu2_uv(scl2_in);\
   return BOOST_PP_CAT(NAME,_impl)(&p10, &m02, &m12);\
}
#define COLLIER_C(R,ARGS,NAME) std::complex<double> Collier::NAME ARGS noexcept\
{\
   const std::complex<double> p10 (p10_in.real(), 0.);\
   const std::complex<double> p21 (p21_in.real(), 0.);\
   const std::complex<double> p20 (p20_in.real(), 0.);\
   const std::complex<double> m02 = m02_in;\
   const std::complex<double> m12 = m12_in;\
   const std::complex<double> m22 = m22_in;\
\
   set_mu2_uv(scl2_in);\
   return BOOST_PP_CAT(NAME,_impl)(&p10, &p21, &p20, &m02, &m12, &m22);\
}
#define COLLIER_D(R,ARGS,NAME) std::complex<double> Collier::NAME ARGS noexcept\
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
   return BOOST_PP_CAT(NAME,_impl)(&p10, &p21, &p32, &p30, &p20, &p31, &m02, &m12, &m22, &m32);\
}

// Fortran wrapper routines
extern "C" {
   void initialize_collier_impl();
   void set_mu2_uv_impl(double*);

   BOOST_PP_SEQ_FOR_EACH(IMPL,0,A_SEQ)
   BOOST_PP_SEQ_FOR_EACH(IMPL,2,B_SEQ)
   BOOST_PP_SEQ_FOR_EACH(IMPL,5,C_SEQ)
   BOOST_PP_SEQ_FOR_EACH(IMPL,9,D_SEQ)

   void get_A_impl(const std::complex<double> [1], COLLIER_ARGS(0));
   void get_B_impl(const std::complex<double> [2], COLLIER_ARGS(2));
   void get_C_impl(const std::complex<double> [7], COLLIER_ARGS(5));
   void get_D_impl(const std::complex<double> [11], COLLIER_ARGS(9));
}

namespace flexiblesusy
{
namespace looplibrary
{

void Collier::initialize() noexcept
{
   futils::swap();
   initialize_collier_impl();
   futils::flush();
   futils::swap();
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

std::complex<double> Collier::A0(A_ARGS) noexcept
{
   const std::complex<double> m02 = m02_in;

   set_mu2_uv(scl2_in);
   return A0_impl(&m02);
}

BOOST_PP_SEQ_FOR_EACH(COLLIER_B,(B_ARGS),B_SEQ)
BOOST_PP_SEQ_FOR_EACH(COLLIER_C,(C_ARGS),C_SEQ)
BOOST_PP_SEQ_FOR_EACH(COLLIER_D,(D_ARGS),D_SEQ)

void Collier::A(std::array<std::complex<double>, 1> &a, A_ARGS) noexcept
{
   const std::complex<double> m02 = m02_in;

   set_mu2_uv(scl2_in);
   get_A_impl(a.data(), &m02);
}

void Collier::B(std::array<std::complex<double>, 2> &b, B_ARGS) noexcept
{
   const std::complex<double> p10 (p10_in.real(), 0.);
   const std::complex<double> m02 = m02_in;
   const std::complex<double> m12 = m12_in;

   set_mu2_uv(scl2_in);
   get_B_impl(b.data(), &p10, &m02, &m12);
}

void Collier::C(std::array<std::complex<double>, 7> &c, C_ARGS) noexcept
{
   const std::complex<double> p10 (p10_in.real(), 0.);
   const std::complex<double> p21 (p21_in.real(), 0.);
   const std::complex<double> p20 (p20_in.real(), 0.);
   const std::complex<double> m02 = m02_in;
   const std::complex<double> m12 = m12_in;
   const std::complex<double> m22 = m22_in;

   set_mu2_uv(scl2_in);
   get_C_impl(c.data(), &p10, &p21, &p20, &m02, &m12, &m22);
}

void Collier::D(std::array<std::complex<double>, 11> &d, D_ARGS) noexcept
{
   const std::complex<double> p10 (p10_in.real(), 0.);
   const std::complex<double> p21 (p21_in.real(), 0.);
   const std::complex<double> p32 (p32_in.real(), 0.);
   const std::complex<double> p30 (p30_in.real(), 0.);
   const std::complex<double> p20 (p20_in.real(), 0.);
   const std::complex<double> p31 (p31_in.real(), 0.);
   const std::complex<double> m02 = m02_in;
   const std::complex<double> m12 = m12_in;
   const std::complex<double> m22 = m22_in;
   const std::complex<double> m32 = m32_in;

   set_mu2_uv(scl2_in);
   get_D_impl(d.data(), &p10, &p21, &p32, &p30, &p20, &p31, &m02, &m12, &m22, &m32);
}

} // namespace looplibrary
} // namespace flexiblesusy
