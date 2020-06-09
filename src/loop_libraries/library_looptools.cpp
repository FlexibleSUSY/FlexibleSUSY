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

#include "library_looptools.hpp"
#include "clooptools.h"
#include "fortran_utils.hpp"
#include <boost/preprocessor/seq/elem.hpp>
#include <boost/preprocessor/seq/for_each_i.hpp>
#include <limits>

#define A_PAIR (A)(aa)
#define B_PAIR (B)(bb)
#define C_PAIR (C)(cc)
#define D_PAIR (D)(dd)

#define REAL(R, DUMMY, ELEM) , ELEM.real()

#define CAT_(a, b) a##b
#define CAT(a, b) CAT_(a, b)

#define LIB(PAIR) BOOST_PP_SEQ_ELEM(0, PAIR)
#define LT(PAIR) BOOST_PP_SEQ_ELEM(1, PAIR)

#define LIB_NAME(PAIR, INDEX) CAT(LIB(PAIR), INDEX)
#define LIB_ARGS(PAIR) CAT(LIB(PAIR), _ARGS)

#define LT_ARGS(PAIR) CAT(LIB(PAIR), _ARGS_SEQ)
#define LT_NAME(PAIR) CAT(LIB(PAIR), 0i)
#define LT_IDX(PAIR, INDEX) CAT(LT(PAIR), INDEX)

#define LT_ONE(_, PAIR, I, INDEX)                                              \
   std::complex<double> Looptools::LIB_NAME(PAIR,                              \
                                            INDEX)(LIB_ARGS(PAIR)) noexcept    \
   {                                                                           \
      set_mu2_uv(scl2_in);                                                     \
      return LT_NAME(PAIR)(LT_IDX(PAIR, INDEX)                                 \
                              BOOST_PP_SEQ_FOR_EACH(REAL, , LT_ARGS(PAIR)));   \
   }

#define LT_ALL(PAIR)                                                           \
   void Looptools::LIB(PAIR)(CAT(LIB(PAIR), coeff_t) & \
                                arr,                                           \
                             CAT(LIB(PAIR), _ARGS)) noexcept                   \
   {                                                                           \
      const int coeffs[] = {BOOST_PP_SEQ_ENUM(                                 \
         BOOST_PP_SEQ_TRANSFORM(APPEND, LT(PAIR), CAT(LIB(PAIR), _CSEQ)))};    \
      ComplexType res[CAT(N, LT(PAIR))];                                       \
      set_mu2_uv(scl2_in);                                                     \
                                                                               \
      CAT(LIB(PAIR), put)                                                      \
      (res BOOST_PP_SEQ_FOR_EACH(REAL, , CAT(LIB(PAIR), _ARGS_SEQ)));          \
      for (int i = 0; i < CAT(LIB(PAIR), _N); ++i) {                           \
         arr.at(i) = res[coeffs[i]];                                           \
      }                                                                        \
   }

namespace flexiblesusy
{
namespace looplibrary
{

Looptools::Looptools() : current_mu2_uv(1.0)
{
   futils::swap();
   ltini();
   futils::flush();
   futils::swap();
}

void Looptools::set_mu2_uv(double scl2_in) noexcept
{
   if (std::abs(scl2_in - this->current_mu2_uv) >
       std::numeric_limits<double>::epsilon()) {
      setmudim(scl2_in);
      this->current_mu2_uv = scl2_in;
   }
}

BOOST_PP_SEQ_FOR_EACH_I(LT_ONE, A_PAIR, A_CSEQ)
BOOST_PP_SEQ_FOR_EACH_I(LT_ONE, B_PAIR, B_CSEQ)
BOOST_PP_SEQ_FOR_EACH_I(LT_ONE, C_PAIR, C_CSEQ)
BOOST_PP_SEQ_FOR_EACH_I(LT_ONE, D_PAIR, D_CSEQ)

LT_ALL(A_PAIR)
LT_ALL(B_PAIR)
LT_ALL(C_PAIR)
LT_ALL(D_PAIR)

} // namespace looplibrary
} // namespace flexiblesusy
