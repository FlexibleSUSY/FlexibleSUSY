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
#include "loop_library_interface.hpp"
#include "numerics.h"
#include <boost/preprocessor/seq/elem.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/seq/for_each_i.hpp>
#include <cmath>
#include <complex>
#include <limits>

#define A_PAIR (A)(aa)
#define B_PAIR (B)(bb)
#define C_PAIR (C)(cc)
#define D_PAIR (D)(dd)

#define REAL(R, DUMMY, ELEM) , ELEM.real()
#define PASS(R, DUMMY, ELEM) , ELEM

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
      return fix_zero_momentum(                                                \
         LT_IDX(PAIR, INDEX),                                                  \
         LT_NAME(PAIR)(LT_IDX(PAIR, INDEX)                                     \
                          BOOST_PP_SEQ_FOR_EACH(REAL, , LT_ARGS(PAIR)))        \
            BOOST_PP_SEQ_FOR_EACH(PASS, , LT_ARGS(PAIR)),                      \
         scl2_in);                                                             \
   }

#define LT_ONE_DERIVATIVE(_, PAIR, I, INDEX)                                              \
   std::complex<double> Looptools::CAT(D, LIB_NAME(PAIR, INDEX))(LIB_ARGS(PAIR)) noexcept    \
   {                                                                           \
      set_mu2_uv(scl2_in);                                                     \
      return LT_NAME(PAIR)(CAT(CAT(d, LT(PAIR)), INDEX)                                 \
                              BOOST_PP_SEQ_FOR_EACH(REAL, , LT_ARGS(PAIR)));   \
   } \


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
      fix_zero_momentum(                                                       \
         arr BOOST_PP_SEQ_FOR_EACH(PASS, , CAT(LIB(PAIR), _ARGS_SEQ)),         \
         scl2_in);                                                             \
   }

#define LT_ALL_DERIVATIVE(PAIR)                                                           \
   void Looptools::CAT(D, LIB(PAIR))(CAT(LIB(PAIR), coeff_t) & \
                                arr,                                           \
                             CAT(LIB(PAIR), _ARGS)) noexcept                   \
   {                                                                           \
      const int coeffs[] = {BOOST_PP_SEQ_ENUM(                                 \
         BOOST_PP_SEQ_TRANSFORM(APPEND, CAT(d, LT(PAIR)), CAT(LIB(PAIR), _CSEQ)))};    \
      ComplexType res[CAT(N, LT(PAIR))];                                       \
      set_mu2_uv(scl2_in);                                                     \
                                                                               \
      CAT(LIB(PAIR), put)                                                      \
      (res BOOST_PP_SEQ_FOR_EACH(REAL, , CAT(LIB(PAIR), _ARGS_SEQ)));          \
      for (int i = 0; i < CAT(CAT(D, LIB(PAIR)), _N); ++i) {                           \
         arr.at(i) = res[coeffs[i]];                                           \
      }                                                                        \
   }

namespace flexiblesusy
{
namespace looplibrary
{

namespace
{

using cd = std::complex<double>;

/// no correction needed for these functions
template <class T, class... Args>
T fix_zero_momentum(int /* index */, T value, Args... /* args */) noexcept
{
   return value;
}

/// no correction needed for these functions
template <class T, class... Args>
void fix_zero_momentum(T& /* coeffs */, Args... /* args */) noexcept
{
}

/// true if all external momenta vanish and the masses are real
bool is_zero_momentum_real_masses(cd p10, cd p21, cd p20, cd m02, cd m12,
                                  cd m22) noexcept
{
   return p10 == 0. && p21 == 0. && p20 == 0. && m02.imag() == 0. &&
          m12.imag() == 0. && m22.imag() == 0.;
}

/// C00 at vanishing external momenta
cd c00_zero_momentum(cd m02, cd m12, cd m22, double scl2) noexcept
{
   return {softsusy::c00(std::sqrt(m02.real()), std::sqrt(m12.real()),
                         std::sqrt(m22.real()), std::sqrt(scl2)),
           0.};
}

/**
 * LoopTools' tensor reduction breaks down at vanishing external momenta
 * (vanishing Gram determinant): C00 is NaN there.  It is replaced by the
 * analytic expression, which the loop-induced observables at zero
 * external momenta (l -> 3l, l -> l conversion) need.
 */
cd fix_zero_momentum(int index, cd value, cd p10, cd p21, cd p20, cd m02,
                     cd m12, cd m22, double scl2) noexcept
{
   if (index == cc00 &&
       is_zero_momentum_real_masses(p10, p21, p20, m02, m12, m22)) {
      return c00_zero_momentum(m02, m12, m22, scl2);
   }
   return value;
}

void fix_zero_momentum(Ccoeff_t& c, cd p10, cd p21, cd p20, cd m02, cd m12,
                       cd m22, double scl2) noexcept
{
   if (is_zero_momentum_real_masses(p10, p21, p20, m02, m12, m22)) {
      c.at(Ccoeffs::c00) = c00_zero_momentum(m02, m12, m22, scl2);
   }
}

} // anonymous namespace

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
BOOST_PP_SEQ_FOR_EACH_I(LT_ONE_DERIVATIVE, B_PAIR, DB_CSEQ)
BOOST_PP_SEQ_FOR_EACH_I(LT_ONE, C_PAIR, C_CSEQ)
BOOST_PP_SEQ_FOR_EACH_I(LT_ONE, D_PAIR, D_CSEQ)

LT_ALL(A_PAIR)
LT_ALL(B_PAIR)
LT_ALL_DERIVATIVE(B_PAIR)
LT_ALL(C_PAIR)
LT_ALL(D_PAIR)

} // namespace looplibrary
} // namespace flexiblesusy
