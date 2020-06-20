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

#include "library_fflite.hpp"
#include "fflite.hpp"
#include "fortran_utils.hpp"
#include <boost/preprocessor/repeat.hpp>
#include <limits>

#define NAN_Q std::numeric_limits<double>::quiet_NaN()
#define SET_TO_NAN(Z, INDEX, DATA) DATA[INDEX] = {NAN_Q, NAN_Q};
#define UNDEFINED(R, ARGS, NAME)                                               \
   std::complex<double> Fflite::NAME ARGS noexcept { return {NAN_Q, NAN_Q}; }

namespace flexiblesusy
{
namespace looplibrary
{

Fflite::Fflite()
{
   futils::swap();
   ltini_();
   futils::flush();
   futils::swap();
}

std::complex<double> Fflite::A0(A_ARGS) noexcept
{
   std::complex<double> ca0;
   int ier;
   ljffca0_(ca0, 0, scl2_in, m02_in, ier);
   return ca0;
}

std::complex<double> Fflite::B0(B_ARGS) noexcept
{
   std::complex<double> cb0;
   int ier;
   ljffcb0_(cb0, 0, scl2_in, p10_in, m02_in, m12_in, ier);
   return cb0;
}

std::complex<double> Fflite::B1(B_ARGS) noexcept
{
   std::complex<double> cb1;
   int ier;
   std::complex<double> cb0 = B0(p10_in, m02_in, m12_in, scl2_in);
   std::complex<double> ca0i[2] = {A0(m02_in, scl2_in), A0(m12_in, scl2_in)};
   std::complex<double> piDpj[9];
   ljffcot2_(piDpj, p10_in, m02_in, m12_in, m02_in - p10_in, m12_in - p10_in,
             m02_in - m12_in, ier);
   ljffcb1_(cb1, cb0, ca0i, p10_in, m02_in, m12_in, piDpj, ier);
   return cb1;
}

std::complex<double> Fflite::B00(B_ARGS) noexcept
{
   std::complex<double> cb2i[2];
   int ier;
   std::complex<double> cb1 = B1(p10_in, m02_in, m12_in, scl2_in);
   std::complex<double> cb0 = B0(p10_in, m02_in, m12_in, scl2_in);
   std::complex<double> ca0i[2] = {A0(m02_in, scl2_in), A0(m12_in, scl2_in)};
   std::complex<double> piDpj[9];
   ljffcot2_(piDpj, p10_in, m02_in, m12_in, m02_in - p10_in, m12_in - p10_in,
             m02_in - m12_in, ier);
   ljffcb2p_(cb2i, cb1, cb0, ca0i, p10_in, m02_in, m12_in, piDpj, ier);
   return cb2i[1];
}

BOOST_PP_SEQ_FOR_EACH(UNDEFINED, (C_ARGS), C_SEQ)
BOOST_PP_SEQ_FOR_EACH(UNDEFINED, (D_ARGS), D_SEQ)

void Fflite::A(Acoeff_t& a, A_ARGS) noexcept
{
   a.at(0) = A0(m02_in, scl2_in);
}

void Fflite::B(Bcoeff_t& b, B_ARGS) noexcept
{
   b.at(0) = B0(p10_in, m02_in, m12_in, scl2_in);
   b.at(1) = B1(p10_in, m02_in, m12_in, scl2_in);
   b.at(2) = B00(p10_in, m02_in, m12_in, scl2_in);
}

void Fflite::C(Ccoeff_t& c, C_ARGS) noexcept
{
   BOOST_PP_REPEAT(C_N, SET_TO_NAN, c)
}

void Fflite::D(Dcoeff_t& d, D_ARGS) noexcept
{
   BOOST_PP_REPEAT(D_N, SET_TO_NAN, d)
}

} // namespace looplibrary
} // namespace flexiblesusy
