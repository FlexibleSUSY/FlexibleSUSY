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

#ifndef LOOP_LIBRARY_INTERFACE_H
#define LOOP_LIBRARY_INTERFACE_H

#include <boost/preprocessor/seq/enum.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/seq/size.hpp>
#include <boost/preprocessor/seq/transform.hpp>
#include <array>
#include <complex>

#define ARGS_TYPE(R, COMMA, ELEM) std::complex<double> ELEM,
#define VIRTUAL(R, ARGS, NAME) virtual std::complex<double> NAME ARGS = 0;

#define A_ARGS_SEQ (m02_in)
#define B_ARGS_SEQ (p10_in)(m02_in)(m12_in)
#define C_ARGS_SEQ (p10_in)(p21_in)(p20_in)(m02_in)(m12_in)(m22_in)
#define D_ARGS_SEQ                                                             \
   (p10_in)(p21_in)(p32_in)(p30_in)(p20_in)(p31_in)(m02_in)(m12_in)(m22_in)(   \
      m32_in)

#define A_ARGS BOOST_PP_SEQ_FOR_EACH(ARGS_TYPE, , A_ARGS_SEQ) double scl2_in
#define B_ARGS BOOST_PP_SEQ_FOR_EACH(ARGS_TYPE, , B_ARGS_SEQ) double scl2_in
#define C_ARGS BOOST_PP_SEQ_FOR_EACH(ARGS_TYPE, , C_ARGS_SEQ) double scl2_in
#define D_ARGS BOOST_PP_SEQ_FOR_EACH(ARGS_TYPE, , D_ARGS_SEQ) double scl2_in

#define A_CSEQ (0)
#define B_CSEQ (0)(1)(00)
#define C_CSEQ (0)(1)(2)(00)(11)(12)(22)
#define D_CSEQ (0)(1)(2)(3)(00)(11)(12)(13)(22)(23)(33)

#define A_N BOOST_PP_SEQ_SIZE(A_CSEQ)
#define B_N BOOST_PP_SEQ_SIZE(B_CSEQ)
#define C_N BOOST_PP_SEQ_SIZE(C_CSEQ)
#define D_N BOOST_PP_SEQ_SIZE(D_CSEQ)

#define APPEND(s, data, elem) data##elem
#define A_SEQ BOOST_PP_SEQ_TRANSFORM(APPEND, A, A_CSEQ)
#define B_SEQ BOOST_PP_SEQ_TRANSFORM(APPEND, B, B_CSEQ)
#define C_SEQ BOOST_PP_SEQ_TRANSFORM(APPEND, C, C_CSEQ)
#define D_SEQ BOOST_PP_SEQ_TRANSFORM(APPEND, D, D_CSEQ)

namespace flexiblesusy
{
namespace looplibrary
{

using Acoeff_t = std::array<std::complex<double>, A_N>;
using Bcoeff_t = std::array<std::complex<double>, B_N>;
using Ccoeff_t = std::array<std::complex<double>, C_N>;
using Dcoeff_t = std::array<std::complex<double>, D_N>;

enum Acoeffs : int {
   BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(APPEND, a, A_CSEQ))
};
enum Bcoeffs : int {
   BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(APPEND, b, B_CSEQ))
};
enum Ccoeffs : int {
   BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(APPEND, c, C_CSEQ))
};
enum Dcoeffs : int {
   BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(APPEND, d, D_CSEQ))
};

/**
 * @class Loop_library_interface
 * @brief interface for different one loop function libraries with
 * conventions of DE == [arXiv:0709.1075],
 * arguments order as in LT ==
 * [http://www.feynarts.de/looptools/LT215Guide.pdf], filling of given arrays as
 * in CO == [arXiv:1604.06792].
 *
 * Loop_library_interface is the abstract base class for one loop functions.
 * It defines the following set of loop functions names:
 *    one-loop ones:   A, A0
 *    two-loop ones:   B, B0, B1, B00
 *    three-loop ones: C, C0, C1, C2, C00, C11, C12, C22
 *    four-loop ones:  D, D0, D1, D2, D3, D00, D11, D12, D13, D22, D23, D33.
 *
 * For making the notation a little bit shorter the following abbreviations for
 * input momenta and masses are used:
 *  A_ARGS means a single std::complex<double> m02 - squared mass of particle
 *  in a loop (see sec. 1.3.1 of [LT]; name is different).
 *
 *  B_ARGS means a set of std::complex<double> p10, m02, m12 - squared momenta
 * and masses of particle in a loop (see sec. 1.3.2 of [LT]; names are
 * different, order is the same).
 *
 *  C_ARGS means a set of std::complex<double> p10, p21, p20, m02, m12, m22 -
 *  squared momenta and masses of particle in a loop (see sec. 1.3.4 of [LT];
 *  names are different, order is the same).
 *
 *  D_ARGS means a set of std::complex<double> p10, p21, p32, p30, p20, p31,
 * m02, m12, m22, m32 - squared momenta and masses of particle in a loop (see
 * sec. 1.3.5 of [LT]; names are different, order is the same).
 *
 * Functions with numeric indices return std::complex<double> of corresponding
 * Passarino-Veltman coefficient (see r.h.s. of eq. (4.4) in [DE]). Ti and Tij
 * functions accept T_ARGS, scl2  as arguments, where scl2 is squared scale
 * (squared mu of eq. (4.1) in [DE]).
 *
 * Example 1: B0(B_ARGS, scl2) means B0(p10, m02, m12, scl2) with:
 *            p10, m02, m12 being of std::complex<double> type;
 *            scl2 being of double type;
 *            p10 is p^2 from section 1.3.2 of [LT],
 *            m02 is m1^2 from section 1.3.2 of [LT],
 *            m12 is m2^2 from section 1.3.2 of [LT] - order is as in [LT];
 *            scl2 is squared scale (squared mu of eq. (4.1) in [DE]);
 *            returns T^1_0 from eq. (4.4) in [DE] of std::complex<double> type.
 *
 * A, B, C, D functions return void, their first arguments are
 * std::complex<double> arrays of fixed length (passed by a reference), which
 * equals to 1, 2, 7, 11. They fill given array with values of Passarino-Vertman
 * coefficients (inspect table 3 of [CO]). After the first argument goes T_ARGS
 * sequence, then scl2, which is described by the following example:
 *
 * Example 2: C(c, C_ARGS, scl2) means C(c, p10, p21, p20, m02, m12, m22, scl2)
 * with: c being array of std::complex<double> of fixed length 7 filled by C
 * function with C0, C1, C2, C00, C11, C12, C22 coefficients (order is defined
 * by table 3 of [CO], coefficients are defined by r.h.s. of eq. (4.4) in [DE]);
 *            p10 is p1^2 from section 1.3.4 of [LT],
 *            p21 is p2^2 from section 1.3.4 of [LT],
 *            p20 is (p1+p2)^2 from section 1.3.4 of [LT],
 *            m02 is m1^2 from section 1.3.4 of [LT],
 *            m12 is m2^2 from section 1.3.4 of [LT],
 *            m22 is m3^2 from section 1.3.4 of [LT] - order is as in [LT];
 *            scl2 is squared scale (squared mu of eq. (4.1) in [DE]);
 *            returns void.
 */
class Loop_library_interface
{
public:
   BOOST_PP_SEQ_FOR_EACH(VIRTUAL, (A_ARGS), A_SEQ)
   BOOST_PP_SEQ_FOR_EACH(VIRTUAL, (B_ARGS), B_SEQ)
   BOOST_PP_SEQ_FOR_EACH(VIRTUAL, (C_ARGS), C_SEQ)
   BOOST_PP_SEQ_FOR_EACH(VIRTUAL, (D_ARGS), D_SEQ)
   virtual void A(Acoeff_t&, A_ARGS) = 0;
   virtual void B(Bcoeff_t&, B_ARGS) = 0;
   virtual void C(Ccoeff_t&, C_ARGS) = 0;
   virtual void D(Dcoeff_t&, D_ARGS) = 0;
   virtual ~Loop_library_interface() {}
};
} // namespace looplibrary
} // namespace flexiblesusy

#endif // LOOP_LIBRARY_INTERFACE_H
