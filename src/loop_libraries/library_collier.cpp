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

#include "library_collier.hpp"
#include "fortran_utils.hpp"
#include <boost/preprocessor/arithmetic/sub.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <limits>

#define A_ARGS_N BOOST_PP_SUB(BOOST_PP_SEQ_SIZE(A_ARGS_SEQ), 1)
#define B_ARGS_N BOOST_PP_SUB(BOOST_PP_SEQ_SIZE(B_ARGS_SEQ), 1)
#define C_ARGS_N BOOST_PP_SUB(BOOST_PP_SEQ_SIZE(C_ARGS_SEQ), 1)
#define D_ARGS_N BOOST_PP_SUB(BOOST_PP_SEQ_SIZE(D_ARGS_SEQ), 1)

#define COLLIER_TYPE(Z, N, TEXT) const std::complex<double>*,
#define COLLIER_ARGS(N)                                                        \
   BOOST_PP_REPEAT(N, COLLIER_TYPE, ) const std::complex<double>*
#define IMPL(R, ARGS, NAME)                                                    \
   void BOOST_PP_CAT(NAME, _impl)(std::complex<double>*, COLLIER_ARGS(ARGS));

// Imaginary parts of momentum invariants are not suppoted by COLLIER 1.2.4
#define COLLIER_B(R, ARGS, NAME)                                               \
   std::complex<double> Collier::NAME ARGS noexcept                            \
   {                                                                           \
      const std::complex<double> p10(p10_in.real(), 0.);                       \
      const std::complex<double> m02 = m02_in;                                 \
      const std::complex<double> m12 = m12_in;                                 \
      std::complex<double> res = 0.0;                                          \
      set_mu2_uv(scl2_in);                                                     \
      BOOST_PP_CAT(NAME, _impl)(&res, &p10, &m02, &m12);                       \
      return res;                                                              \
   }
#define COLLIER_C(R, ARGS, NAME)                                               \
   std::complex<double> Collier::NAME ARGS noexcept                            \
   {                                                                           \
      const std::complex<double> p10(p10_in.real(), 0.);                       \
      const std::complex<double> p21(p21_in.real(), 0.);                       \
      const std::complex<double> p20(p20_in.real(), 0.);                       \
      const std::complex<double> m02 = m02_in;                                 \
      const std::complex<double> m12 = m12_in;                                 \
      const std::complex<double> m22 = m22_in;                                 \
      std::complex<double> res = 0.0;                                          \
      set_mu2_uv(scl2_in);                                                     \
      BOOST_PP_CAT(NAME, _impl)(&res, &p10, &p21, &p20, &m02, &m12, &m22);     \
      return res;                                                              \
   }
#define COLLIER_D(R, ARGS, NAME)                                               \
   std::complex<double> Collier::NAME ARGS noexcept                            \
   {                                                                           \
      const std::complex<double> p10(p10_in.real(), 0.);                       \
      const std::complex<double> p21(p21_in.real(), 0.);                       \
      const std::complex<double> p32(p32_in.real(), 0.);                       \
      const std::complex<double> p30(p30_in.real(), 0.);                       \
      const std::complex<double> p20(p20_in.real(), 0.);                       \
      const std::complex<double> p31(p31_in.real(), 0.);                       \
      const std::complex<double> m02 = m02_in;                                 \
      const std::complex<double> m12 = m12_in;                                 \
      const std::complex<double> m22 = m22_in;                                 \
      const std::complex<double> m32 = m32_in;                                 \
      std::complex<double> res = 0.0;                                          \
      set_mu2_uv(scl2_in);                                                     \
      BOOST_PP_CAT(NAME, _impl)                                                \
      (&res, &p10, &p21, &p32, &p30, &p20, &p31, &m02, &m12, &m22, &m32);      \
      return res;                                                              \
   }

// Fortran wrapper routines
extern "C" {
void initialize_collier_impl();
void set_mu2_uv_impl(double*);

BOOST_PP_SEQ_FOR_EACH(IMPL, A_ARGS_N, A_SEQ)
BOOST_PP_SEQ_FOR_EACH(IMPL, B_ARGS_N, B_SEQ)
BOOST_PP_SEQ_FOR_EACH(IMPL, C_ARGS_N, C_SEQ)
BOOST_PP_SEQ_FOR_EACH(IMPL, D_ARGS_N, D_SEQ)

void get_A_impl(const std::complex<double>[A_N], COLLIER_ARGS(A_ARGS_N));
void get_B_impl(const std::complex<double>[B_N], COLLIER_ARGS(B_ARGS_N));
void get_C_impl(const std::complex<double>[C_N], COLLIER_ARGS(C_ARGS_N));
void get_D_impl(const std::complex<double>[D_N], COLLIER_ARGS(D_ARGS_N));
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
   if (std::abs(scl2 - this->current_mu2_uv) >
       std::numeric_limits<double>::epsilon()) {
      set_mu2_uv_impl(&scl2);
      this->current_mu2_uv = scl2;
   }
}

std::complex<double> Collier::A0(A_ARGS) noexcept
{
   const std::complex<double> m02 = m02_in;
   std::complex<double> res = 0.0;
   set_mu2_uv(scl2_in);
   A0_impl(&res, &m02);
   return res;
}

BOOST_PP_SEQ_FOR_EACH(COLLIER_B, (B_ARGS), B_SEQ)
BOOST_PP_SEQ_FOR_EACH(COLLIER_C, (C_ARGS), C_SEQ)
BOOST_PP_SEQ_FOR_EACH(COLLIER_D, (D_ARGS), D_SEQ)

void Collier::A(Acoeff_t& a, A_ARGS) noexcept
{
   const std::complex<double> m02 = m02_in;

   set_mu2_uv(scl2_in);
   get_A_impl(a.data(), &m02);
}

void Collier::B(Bcoeff_t& b, B_ARGS) noexcept
{
   const std::complex<double> p10(p10_in.real(), 0.);
   const std::complex<double> m02 = m02_in;
   const std::complex<double> m12 = m12_in;

   set_mu2_uv(scl2_in);
   get_B_impl(b.data(), &p10, &m02, &m12);
}

void Collier::C(Ccoeff_t& c, C_ARGS) noexcept
{
   const std::complex<double> p10(p10_in.real(), 0.);
   const std::complex<double> p21(p21_in.real(), 0.);
   const std::complex<double> p20(p20_in.real(), 0.);
   const std::complex<double> m02 = m02_in;
   const std::complex<double> m12 = m12_in;
   const std::complex<double> m22 = m22_in;

   set_mu2_uv(scl2_in);
   get_C_impl(c.data(), &p10, &p21, &p20, &m02, &m12, &m22);
}

void Collier::D(Dcoeff_t& d, D_ARGS) noexcept
{
   const std::complex<double> p10(p10_in.real(), 0.);
   const std::complex<double> p21(p21_in.real(), 0.);
   const std::complex<double> p32(p32_in.real(), 0.);
   const std::complex<double> p30(p30_in.real(), 0.);
   const std::complex<double> p20(p20_in.real(), 0.);
   const std::complex<double> p31(p31_in.real(), 0.);
   const std::complex<double> m02 = m02_in;
   const std::complex<double> m12 = m12_in;
   const std::complex<double> m22 = m22_in;
   const std::complex<double> m32 = m32_in;

   set_mu2_uv(scl2_in);
   get_D_impl(d.data(), &p10, &p21, &p32, &p30, &p20, &p31, &m02, &m12, &m22,
              &m32);
}

} // namespace looplibrary
} // namespace flexiblesusy
