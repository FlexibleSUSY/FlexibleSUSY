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

#include "library_softsusy.hpp"
#include "numerics.h"
#include <limits>

#define NAN_Q std::numeric_limits<double>::quiet_NaN()
#define UNDEFINED(R, ARGS, NAME)                                               \
   std::complex<double> Softsusy::NAME ARGS noexcept { return {NAN_Q, NAN_Q}; }

namespace flexiblesusy
{
namespace looplibrary
{

std::complex<double> Softsusy::A0(A_ARGS) noexcept
{
   double m = std::sqrt(m02_in.real());
   double q = std::sqrt(scl2_in);

   return {softsusy::a0(m, q), 0.0};
}

std::complex<double> Softsusy::B0(B_ARGS) noexcept
{
   double p = std::sqrt(p10_in.real());
   double m1 = std::sqrt(m02_in.real());
   double m2 = std::sqrt(m12_in.real());
   double q = std::sqrt(scl2_in);

   return {softsusy::b0(p, m1, m2, q), 0.0};
}

std::complex<double> Softsusy::B1(B_ARGS) noexcept
{
   double p = std::sqrt(p10_in.real());
   double m1 = std::sqrt(m02_in.real());
   double m2 = std::sqrt(m12_in.real());
   double q = std::sqrt(scl2_in);

   return {(-1) * softsusy::b1(p, m1, m2, q), 0.0};
}

std::complex<double> Softsusy::B00(B_ARGS) noexcept
{
   double p = std::sqrt(p10_in.real());
   double m1 = std::sqrt(m02_in.real());
   double m2 = std::sqrt(m12_in.real());
   double q = std::sqrt(scl2_in);

   return {softsusy::b22(p, m1, m2, q), 0.0};
}

std::complex<double> Softsusy::C0(C_ARGS) noexcept
{
   double m1 = std::sqrt(m02_in.real());
   double m2 = std::sqrt(m12_in.real());
   double m3 = std::sqrt(m22_in.real());

   return {softsusy::c0(m1, m2, m3), 0.0};
}

std::complex<double> Softsusy::C00(C_ARGS) noexcept
{
   double m1 = std::sqrt(m02_in.real());
   double m2 = std::sqrt(m12_in.real());
   double m3 = std::sqrt(m22_in.real());
   double q = std::sqrt(scl2_in);

   return {softsusy::c00(m1, m2, m3, q), 0.0};
}

BOOST_PP_SEQ_FOR_EACH(UNDEFINED, (C_ARGS), (C1)(C2)(C11)(C12)(C22))

std::complex<double> Softsusy::D0(D_ARGS) noexcept
{
   double m1 = std::sqrt(m02_in.real());
   double m2 = std::sqrt(m12_in.real());
   double m3 = std::sqrt(m22_in.real());
   double m4 = std::sqrt(m32_in.real());

   return {softsusy::d0(m1, m2, m3, m4), 0.0};
}

std::complex<double> Softsusy::D00(D_ARGS) noexcept
{
   double m1 = std::sqrt(m02_in.real());
   double m2 = std::sqrt(m12_in.real());
   double m3 = std::sqrt(m22_in.real());
   double m4 = std::sqrt(m32_in.real());

   return {softsusy::d27(m1, m2, m3, m4), 0.0};
}

BOOST_PP_SEQ_FOR_EACH(UNDEFINED, (D_ARGS),
                      (D1)(D11)(D12)(D13)(D2)(D22)(D23)(D3)(D33))

void Softsusy::A(Acoeff_t& a, A_ARGS) noexcept
{
   double m = std::sqrt(m02_in.real());
   double q = std::sqrt(scl2_in);

   a.at(0) = {softsusy::a0(m, q), 0.0};
}

void Softsusy::B(Bcoeff_t& b, B_ARGS) noexcept
{
   double p = std::sqrt(p10_in.real());
   double m1 = std::sqrt(m02_in.real());
   double m2 = std::sqrt(m12_in.real());
   double q = std::sqrt(scl2_in);

   b.at(0) = {softsusy::b0(p, m1, m2, q), 0.0};
   b.at(1) = {(-1) * softsusy::b1(p, m1, m2, q), 0.0};
   b.at(2) = {softsusy::b22(p, m1, m2, q), 0.0};
}

void Softsusy::C(Ccoeff_t& c, C_ARGS) noexcept
{
   double m1 = std::sqrt(m02_in.real());
   double m2 = std::sqrt(m12_in.real());
   double m3 = std::sqrt(m22_in.real());
   double q = std::sqrt(scl2_in);
   std::complex<double> undefined = {NAN_Q, NAN_Q};

   c.at(0) = {softsusy::c0(m1, m2, m3), 0.0};
   c.at(1) = undefined;
   c.at(2) = undefined;
   c.at(3) = {softsusy::c00(m1, m2, m3, q), 0.0};
   c.at(4) = undefined;
   c.at(5) = undefined;
   c.at(6) = undefined;
}

void Softsusy::D(Dcoeff_t& d, D_ARGS) noexcept
{
   double m1 = std::sqrt(m02_in.real());
   double m2 = std::sqrt(m12_in.real());
   double m3 = std::sqrt(m22_in.real());
   double m4 = std::sqrt(m32_in.real());
   std::complex<double> undefined = {NAN_Q, NAN_Q};

   d.at(0) = {softsusy::d0(m1, m2, m3, m4), 0.0};
   d.at(1) = undefined;
   d.at(2) = undefined;
   d.at(3) = undefined;
   d.at(4) = {softsusy::d27(m1, m2, m3, m4), 0.0};
   d.at(5) = undefined;
   d.at(6) = undefined;
   d.at(7) = undefined;
   d.at(8) = undefined;
   d.at(9) = undefined;
   d.at(10) = undefined;
}

} // namespace looplibrary
} // namespace flexiblesusy
