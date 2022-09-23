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

#include "sm_threeloophiggs.hpp"
#include <cmath>

namespace flexiblesusy {
namespace sm_threeloophiggs {

namespace {

constexpr double threeLoop = 2.539456721913701978e-07; // 1/(4 Pi)^6

double sqr(double x) noexcept { return x*x; }

} // anonymous namespace

/**
 * Standard Model Higgs self-energy 3-loop, \f$O(\alpha_t
 * \alpha_s^2)\f$.  Taken from arxiv:1407.4336, Eq. (3.3).
 *
 * @note The result contains the 3-loop tadpole diagrams.  It is
 * therefore not 1-particle irreducible (1PI).
 *
 * @param scale renormalization scale
 * @param mt MS-bar top mass
 * @param yt MS-bar Yukawa coupling
 * @param g3 MS-bar strong gauge coupling
 *
 * @return real part of 3-loop correction \f$O(\alpha_t \alpha_s^2)\f$
 */
double delta_mh_3loop_at_as_as_sm(
   double scale, double mt, double yt, double g3)
{
   const double yt2 = sqr(yt);
   const double mt2 = sqr(mt);
   const double g34 = sqr(sqr(g3));
   const double Q2 = sqr(scale);
   const double LogT = std::log(mt2 / Q2);
   const double LogT2 = sqr(LogT);
   const double LogT3 = LogT2*LogT;

   const double result =
      g34*yt2*mt2*(248.1215180432007 + 839.1966169377614*LogT
                   + 160*LogT2 - 736*LogT3);

   return result * threeLoop;
}

/**
 * Standard Model Higgs self-energy 3-loop,
 * \f$O(\alpha_t^2\alpha_s)\f$.  Taken from arxiv:1407.4336, Eq. (3.3).
 *
 * @note The result contains the 3-loop tadpole diagrams.  It is
 * therefore not 1-particle irreducible (1PI).
 *
 * @param scale renormalization scale
 * @param mt MS-bar top mass
 * @param yt MS-bar Yukawa coupling
 * @param g3 MS-bar strong gauge coupling
 *
 * @return real part of 3-loop correction \f$O(\alpha_t^2\alpha_s) \f$
 */
double delta_mh_3loop_at_at_as_sm(
   double scale, double mt, double yt, double g3)
{
   const double yt4 = sqr(sqr(yt));
   const double mt2 = sqr(mt);
   const double g32 = sqr(g3);
   const double Q2 = sqr(scale);
   const double LogT = std::log(mt2 / Q2);
   const double LogT2 = sqr(LogT);
   const double LogT3 = LogT2*LogT;

   const double result =
      g32*yt4*mt2*(2764.365124334015 + 1283.715638285500*LogT
                   - 360*LogT2 + 240*LogT3);

   return result * threeLoop;
}

/**
 * Standard Model Higgs self-energy 3-loop, \f$O(\alpha_t^3)\f$.
 * Taken from arxiv:1407.4336, Eq. (3.4).
 *
 * @note The result contains the 3-loop tadpole diagrams.  It is
 * therefore not 1-particle irreducible (1PI).
 *
 * @param scale renormalization scale
 * @param mt MS-bar top mass
 * @param yt MS-bar Yukawa coupling
 * @param mh MS-bar Higgs mass
 *
 * @return real part of 3-loop correction \f$O(\alpha_t^3) \f$
 */
double delta_mh_3loop_at_at_at_sm(
   double scale, double mt, double yt, double mh)
{
   const double yt6 = sqr(yt*yt*yt);
   const double mt2 = sqr(mt);
   const double mh2 = sqr(mh);
   const double Q2 = sqr(scale);
   const double LogH = std::log(mh2 / Q2);
   const double LogT = std::log(mt2 / Q2);
   const double LogT2 = sqr(LogT);
   const double LogT3 = LogT2*LogT;

   const double result =
      yt6*mt2*(-3199.016554815089 + 36*LogH - 2653.510765697467*LogT
               + 756*LogH*LogT + 27.*0.5*LogT2 + 324*LogH*LogT2 - 225*LogT3);

   return result * threeLoop;
}

} // namespace sm_threeloophiggs
} // namespace flexiblesusy
