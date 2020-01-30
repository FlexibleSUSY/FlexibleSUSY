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

#ifndef LOOP_LIBRARY_INTERFACE
#define LOOP_LIBRARY_INTERFACE

#include <complex>
#include <boost/preprocessor/seq/for_each.hpp>

#define ARGS_TYPE(R,COMMA,ELEM) std::complex<double> ELEM,
#define VIRTUAL(R,ARGS,NAME) virtual std::complex<double> NAME ARGS = 0;

#define A_ARGS_SEQ (m02_in)
#define B_ARGS_SEQ (p10_in)(m02_in)(m12_in)
#define C_ARGS_SEQ (p10_in)(p21_in)(p20_in)(m02_in)(m12_in)(m22_in)
#define D_ARGS_SEQ (p10_in)(p21_in)(p32_in)(p30_in)(p20_in)(p31_in)(m02_in)(m12_in)(m22_in)(m32_in)

#define A_ARGS BOOST_PP_SEQ_FOR_EACH(ARGS_TYPE,,A_ARGS_SEQ) double scl2_in
#define B_ARGS BOOST_PP_SEQ_FOR_EACH(ARGS_TYPE,,B_ARGS_SEQ) double scl2_in
#define C_ARGS BOOST_PP_SEQ_FOR_EACH(ARGS_TYPE,,C_ARGS_SEQ) double scl2_in
#define D_ARGS BOOST_PP_SEQ_FOR_EACH(ARGS_TYPE,,D_ARGS_SEQ) double scl2_in

#define A_SEQ (A0)
#define B_SEQ (B0)(B1)(B00)
#define C_SEQ (C0)(C1)(C2)(C00)(C11)(C12)(C22)
#define D_SEQ (D0)(D00)(D1)(D11)(D12)(D13)(D2)(D22)(D23)(D3)(D33)

namespace looplibrary
{
class Loop_library_interface
{
   public:
      BOOST_PP_SEQ_FOR_EACH(VIRTUAL,(A_ARGS),A_SEQ)
      BOOST_PP_SEQ_FOR_EACH(VIRTUAL,(B_ARGS),B_SEQ)
      BOOST_PP_SEQ_FOR_EACH(VIRTUAL,(C_ARGS),C_SEQ)
      BOOST_PP_SEQ_FOR_EACH(VIRTUAL,(D_ARGS),D_SEQ)
      virtual void A(std::complex<double> (&)[1], A_ARGS) = 0;
      virtual void B(std::complex<double> (&)[2], B_ARGS) = 0;
      virtual void C(std::complex<double> (&)[7], C_ARGS) = 0;
      virtual void D(std::complex<double> (&)[11], D_ARGS) = 0;
};
} // namespace looplibrary

#endif // LOOP_LIBRARY_INTERFACE
