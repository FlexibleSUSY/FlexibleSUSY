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

#include "sm_twoloophiggs.hpp"
#include "Li2.hpp"
#include "loop_libraries/loop_library.hpp"
#include "threshold_loop_functions.hpp"

#include <cmath>

namespace flexiblesusy {
namespace sm_twoloophiggs {

namespace {

constexpr double Pi      = 3.141592653589793;
constexpr double oneLoop = 6.332573977646110963e-03; // 1/(4 Pi)^2
constexpr double twoLoop = 4.010149318236068752e-05; // 1/(4 Pi)^4

double sqr(double x) noexcept { return x*x; }
double pow4(double x) noexcept { return sqr(sqr(x)); }

} // anonymous namespace

/**
 * Standard Model Higgs 1-loop contribution (Landau gauge).
 * Taken from arxiv:1205.6497, Eq. (16).
 *
 * @note The result contains the 1-loop tadpole diagrams.  It is
 * therefore not 1-particle irreducible (1PI).
 *
 * @warning The result is in Landau gauge (\f$\xi = 0\f$).
 *
 * @warning In the diagram that has a Higgs loop, the momentum has
 * been set to the tree-level Higgs mass, \f$p^2 = 2\lambda v^2\f$.
 *
 * @param p momentum
 * @param scale renormalization scale
 * @param mt MS-bar top mass
 * @param yt MS-bar Yukawa coupling
 * @param v MS-bar VEV
 * @param gY MS-bar hypercharge gauge coupling
 * @param g2 MS-bar SU(2)_L gauge coupling
 * @param lambda MS-bar 4-Higgs coupling
 *
 * @return real part of 1-loop correction
 */
double delta_mh_1loop_sm(
   double p, double scale, double mt, double yt,
   double v, double gY, double g2, double lambda)
{
   const double yt2 = sqr(yt);
   const double mt2 = sqr(mt);
   const double p2 = sqr(p);
   const double Q2 = sqr(scale);
   const double lambda2 = sqr(lambda);
   const double v2 = sqr(v);
   const double g22 = sqr(g2);
   const double g24 = sqr(g22);
   const double gp2 = sqr(gY);
   const double G2 = g22 + gp2;
   const double G4 = sqr(G2);
   const double mW2 = g22 * v2/4.;
   const double mZ2 = G2 * v2/4.;
   const double mH2 = 2*lambda*v2;
   const double LogW = std::log(mW2 / Q2);
   const double LogZ = std::log(mZ2 / Q2);
   const double LogH = std::log(mH2 / Q2);

   const double result =
      (+3*yt2*(4*mt2 - p2)*Loop_library::get().B0(p2,mt2,mt2,Q2).real()
       +6*lambda2*v2*(3*LogH-6+Pi*std::sqrt(3))
       -v2/4.*(3*g24-8*lambda*g22+16*lambda2)*Loop_library::get().B0(p2,mW2,mW2,Q2).real()
       -v2/8.*(3*G4-8*lambda*G2+16*lambda2)*Loop_library::get().B0(p2,mZ2,mZ2,Q2).real()
       +2*mW2*(g22-2*lambda*(LogW-1))
       +mZ2*(G2-2*lambda*(LogZ-1))
      );

   return result * oneLoop;
}

/**
 * Standard Model Higgs 1-loop contribution, \f$O(\alpha_t)\f$.
 * Taken from arxiv:1205.6497, Eq. (16).
 *
 * @note The result contains the 1-loop top quark tadpole diagram.  It
 * is therefore not 1-particle irreducible (1PI).
 *
 * @param p momentum
 * @param scale renormalization scale
 * @param mt MS-bar top mass
 * @param yt MS-bar Yukawa coupling
 *
 * @return real part of 1-loop correction O(alpha_t)
 */
double delta_mh_1loop_at_sm(
   double p, double scale, double mt, double yt)
{
   const double yt2 = sqr(yt);
   const double mt2 = sqr(mt);
   const double p2 = sqr(p);
   const double Q2 = sqr(scale);

   const double result =
      3*yt2*(4.*mt2 - p2)*Loop_library::get().B0(p2,mt2,mt2,Q2).real();

   return result * oneLoop;
}


double delta_mh_1loop_ab_sm(
   double p, double scale, double mb, double yb)
{
   const double yb2 = sqr(yb);
   const double mb2 = sqr(mb);
   const double p2 = sqr(p);
   const double Q2 = sqr(scale);

   const double result =
      3*yb2*(4.*mb2 - p2)*Loop_library::get().B0(p2,mb2,mb2,Q2).real();

   return result * oneLoop;
}

double delta_mh_1loop_atau_sm(
   double p, double scale, double mtau, double ytau)
{
   const double yt2 = sqr(ytau);
   const double mt2 = sqr(mtau);
   const double p2 = sqr(p);
   const double Q2 = sqr(scale);

   const double result =
      yt2*(4.*mt2 - p2)*Loop_library::get().B0(p2,mt2,mt2,Q2).real();

   return result * oneLoop;
}

/*
 * Derivatives w.r.t. to y_{t,b,tau}, v, p
 *
 * Note: Because 2L-Se are evaluated from the effective potential
 * all formal two-loop contr must be derived for p=0.
 *
 */

double delta_mh_1loop_at_sm_deriv_yt(
   double p, double scale, double mt, double yt)
{
   const double yt2 = sqr(yt);
   const double mt2 = sqr(mt);
   const double p2 = sqr(p);
   const double Q2 = sqr(scale);
   const double LogT = std::log(mt2 / Q2);


   const double result =
      -24. * yt * mt2 * (2. * LogT + 1.);

   return  result * oneLoop;
}

double delta_mh_1loop_ab_sm_deriv_yb(
   double p, double scale, double mb, double yb)
{

   const double yb2 = sqr(yb);
   const double mb2 = sqr(mb);
   const double p2 = sqr(p);
   const double Q2 = sqr(scale);
   const double LogB = std::log(mb2 / Q2);


   const double result =
      -24. * yb * mb2 *(2 * LogB + 1.);

   return  result * oneLoop;
}

double delta_mh_1loop_atau_sm_deriv_ytau(
   double p, double scale, double mtau, double ytau)
{

   const double ytau2 = sqr(ytau);
   const double mtau2 = sqr(mtau);
   const double p2 = sqr(p);
   const double Q2 = sqr(scale);
   const double LogT = std::log(mtau2 / Q2);


   const double result =
      -8. * ytau * mtau2 * (2. * LogT + 1.);

   return  result * oneLoop;
}

double delta_mh_1loop_at_sm_deriv_v(
   double p, double scale, double mt, double yt)
{

   const double yt2 = sqr(yt);
   const double mt2 = sqr(mt);
   const double p2 = sqr(p);
   const double Q2 = sqr(scale);
   const double LogT = std::log(mt2 / Q2);


   const double result =
      - 12. * 1.41421356237 * yt * yt2 * mt * ( LogT + 1.);

   return  result * oneLoop;
}

double delta_mh_1loop_ab_sm_deriv_v(
   double p, double scale, double mb, double yb)
{

   const double yb2 = sqr(yb);
   const double mb2 = sqr(mb);
   const double p2 = sqr(p);
   const double Q2 = sqr(scale);
   const double LogB = std::log(mb2 / Q2);


   const double result =
      - 12. * 1.41421356237 * yb * yb2 * mb * ( LogB + 1.);

   return  result * oneLoop;
}

double delta_mh_1loop_atau_sm_deriv_v(
   double p, double scale, double mtau, double ytau)
{

   const double ytau2 = sqr(ytau);
   const double mtau2 = sqr(mtau);
   const double p2 = sqr(p);
   const double Q2 = sqr(scale);
   const double LogTau = std::log(mtau2 / Q2);


   const double result =
      - 4. * 1.41421356237 * ytau * ytau2 * mtau * ( LogTau + 1.);

   return  result * oneLoop;
}


double delta_mh_1loop_at_sm_deriv_p2(
   double p, double scale, double mt, double yt)
{

   const double yt2 = sqr(yt);
   const double mt2 = sqr(mt);
   const double p2 = sqr(p);
   const double Q2 = sqr(scale);
   const double LogT = std::log(mt2 / Q2);


   const double result =
       yt2 * (3.* LogT + 2.);

   return  result * oneLoop;
}

double delta_mh_1loop_ab_sm_deriv_p2(
   double p, double scale, double mb, double yb)
{

   const double yb2 = sqr(yb);
   const double mb2 = sqr(mb);
   const double p2 = sqr(p);
   const double Q2 = sqr(scale);
   const double LogB = std::log(mb2 / Q2);


   const double result =
       yb2 * (3.* LogB + 2.);

   return  result * oneLoop;
}


double delta_mh_1loop_atau_sm_deriv_p2(
   double p, double scale, double mtau, double ytau)
{

   const double ytau2 = sqr(ytau);
   const double mtau2 = sqr(mtau);
   const double p2 = sqr(p);
   const double Q2 = sqr(scale);
   const double LogTau = std::log(mtau2 / Q2);


   const double result =
       ytau2 * ( LogTau + 2./3.);

   return  result * oneLoop;
}

/**
 * Standard Model Higgs self-energy 2-loop, \f$O(\alpha_t
 * \alpha_s)\f$, including momentum dependence.
 *
 * @warning The result is in Landau gauge (\f$\xi = 0\f$).
 *
 * @param p2    squared momentum
 * @param scale renormalization scale
 * @param mt MS-bar top mass
 * @param yt MS-bar Yukawa coupling
 * @param g3 MS-bar strong gauge coupling
 *
 * @return real part of 2-loop self-energy \f$O(\alpha_t \alpha_s)\f$
 */
double self_energy_higgs_2loop_at_as_sm(
   double p2, double scale, double mt, double yt, double g3)
{
   const double yt2 = sqr(yt);
   const double g32 = sqr(g3);
   const double t = sqr(mt);
   const double q = sqr(scale);
   const double lnt = std::log(t/q);

   const double result =
      1./135. * g32 * yt2 * (-10800*t - 1665*p2 + (122*sqr(p2))/t +
                             540*(12*t + 5*p2) * lnt -
                             1620*(12*t - p2) * sqr(lnt));

   return result * twoLoop;
}

/**
 * Standard Model Higgs self-energy 2-loop, \f$O(\alpha_b
 * \alpha_s)\f$, for zero momentum.
 *
 * @warning The result is in Landau gauge (\f$\xi = 0\f$).
 *
 * @param p2    squared momentum (not used)
 * @param scale renormalization scale
 * @param mb MS-bar bottom mass
 * @param yb MS-bar Yukawa coupling
 * @param g3 MS-bar strong gauge coupling
 *
 * @return real part of 2-loop self-energy \f$O(\alpha_b \alpha_s)\f$
 */
double self_energy_higgs_2loop_ab_as_sm(
   double /* p2 */, double scale, double mb, double yb, double g3)
{
   return self_energy_higgs_2loop_at_as_sm(0., scale, mb, yb, g3);
}

/**
 * Standard Model Higgs tadpole 2-loop, \f$O(\alpha_t \alpha_s)\f$.
 *
 * @warning The result is in Landau gauge (\f$\xi = 0\f$).
 *
 * @param scale renormalization scale
 * @param mt MS-bar top mass
 * @param yt MS-bar Yukawa coupling
 * @param g3 MS-bar strong gauge coupling
 *
 * @return real part of 2-loop self-energy \f$O(\alpha_t \alpha_s)\f$
 */
double tadpole_higgs_2loop_at_as_sm(
   double scale, double mt, double yt, double g3)
{
   const double yt2 = sqr(yt);
   const double g32 = sqr(g3);
   const double t = sqr(mt);
   const double q = sqr(scale);
   const double lnt = std::log(t/q);

   const double result =
      -16 * g32 * t * yt2 * (5 - 5*lnt + 3*sqr(lnt));

   return result * twoLoop;
}

/**
 * Standard Model Higgs tadpole 2-loop, \f$O(\alpha_b \alpha_s)\f$.
 *
 * @warning The result is in Landau gauge (\f$\xi = 0\f$).
 *
 * @param scale renormalization scale
 * @param mb MS-bar bottom mass
 * @param yb MS-bar Yukawa coupling
 * @param g3 MS-bar strong gauge coupling
 *
 * @return real part of 2-loop self-energy \f$O(\alpha_b \alpha_s)\f$
 */
double tadpole_higgs_2loop_ab_as_sm(
   double scale, double mb, double yb, double g3)
{
   return tadpole_higgs_2loop_at_as_sm(scale, mb, yb, g3);
}

double delta_mh_2loop_at_as_sm(
   double p2, double scale, double mt, double yt, double g3)
{
   return - self_energy_higgs_2loop_at_as_sm(p2, scale, mt, yt, g3)
      + tadpole_higgs_2loop_at_as_sm(scale, mt, yt, g3);
}

double delta_mh_2loop_ab_as_sm(
   double p2, double scale, double mb, double yb, double g3)
{
   return - self_energy_higgs_2loop_ab_as_sm(p2, scale, mb, yb, g3)
      + tadpole_higgs_2loop_ab_as_sm(scale, mb, yb, g3);
}

namespace {

/**
 * \f$m_h\f$-dependent part common to the tadpole and self-energy of
 * one quark (\f$N_c = 3\f$) with \f$h\f$ and \f$G^0\f$ exchange.
 * Contains \f$m_h^2 (m_h^2 - c m_f^2) \phi(m_f^2,m_f^2,m_h^2)\f$, with
 * \f$\phi\f$ of Davydychev and Tausk, Nucl. Phys. B397 (1993) 23,
 * evaluated as \f$m_h^2 (m_h^2 - c m_f^2)/m_f^2 \Phi(m_h^2/m_f^2, 1)\f$,
 * which is finite for \f$m_h \to 0\f$.
 *
 * @param mf fermion mass
 * @param mh Higgs mass, \f$m_h = 0\f$ gives the limit \f$m_h \to 0\f$
 * @param scale renormalization scale
 * @param c coefficient
 */
double higgs_terms(double mf, double mh, double scale, double c) noexcept
{
   // limit m_h -> 0 (avoids 0*log(0))
   if (mh == 0) {
      return 0;
   }

   const double mf2 = sqr(mf);
   const double mh2 = sqr(mh);
   const double Q2 = sqr(scale);
   const double L = std::log(mf2/Q2);
   const double Lh = std::log(mh2/Q2);

   return mh2*(13.5 - 4.5*sqr(L) + 9*L*Lh - 12*Lh)
      + 4.5*mh2*(mh2 - c*mf2)/mf2*threshold_loop_functions::phi_xyz(mh2, mf2, mf2);
}

/**
 * 2-loop tadpole \f$T\f$ of one quark (\f$N_c = 3\f$) with \f$h\f$ and
 * \f$G^0\f$ exchange (\f$m_G = 0\f$), without the factor
 * \f$1/(16\pi^2)^2\f$.
 *
 * @param yf Yukawa coupling of the fermion
 * @param mf fermion mass
 * @param mh Higgs mass
 * @param scale renormalization scale
 */
double tadpole_neutral(double yf, double mf, double mh, double scale) noexcept
{
   const double mf2 = sqr(mf);
   const double L = std::log(mf2/sqr(scale));

   return pow4(yf)*(36*mf2 - 30*mf2*L + 6*mf2*sqr(L)
                    + higgs_terms(mf, mh, scale, 4));
}

/**
 * 2-loop self-energy \f$\Sigma(0)\f$ of one quark (\f$N_c = 3\f$) with
 * \f$h\f$ and \f$G^0\f$ exchange (\f$m_G = 0\f$), without the factor
 * \f$1/(16\pi^2)^2\f$.
 *
 * @param yf Yukawa coupling of the fermion
 * @param mf fermion mass
 * @param mh Higgs mass
 * @param scale renormalization scale
 */
double self_energy_neutral(double yf, double mf, double mh, double scale) noexcept
{
   const double mf2 = sqr(mf);
   const double L = std::log(mf2/sqr(scale));

   return pow4(yf)*(48*mf2 - 66*mf2*L + 18*mf2*sqr(L)
                    + higgs_terms(mf, mh, scale, 8));
}

/**
 * 2-loop tadpole \f$T\f$ of the top-bottom doublet with \f$G^\pm\f$
 * exchange (\f$m_G = 0\f$), without the factor \f$1/(16\pi^2)^2\f$.
 *
 * @param yt top Yukawa coupling
 * @param yb bottom Yukawa coupling
 * @param mt top mass
 * @param mb bottom mass
 * @param scale renormalization scale
 */
double tadpole_charged(
   double yt, double yb, double mt, double mb, double scale) noexcept
{
   const double mt2 = sqr(mt);
   const double mb2 = sqr(mb);
   const double Q2 = sqr(scale);
   const double Lt = std::log(mt2/Q2);
   const double Lb = std::log(mb2/Q2);
   const double Lt2 = sqr(Lt);
   const double li = Li2(1 - mb2/mt2);
   const double yt2 = sqr(yt);
   const double yb2 = sqr(yb);

   return sqr(yt2)*(9*(mt2 + mb2) - 9*mt2*Lt + 3*mt2*Lt2 - 6*mb2*Lb + 6*mt2*li)
      + sqr(yb2)*(9*(mt2 + mb2) - 6*mt2*Lt - 9*mb2*Lb + 6*mb2*Lt*Lb
                  - 3*mb2*Lt2 - 6*mb2*li)
      + yt2*yb2*(-18*(mt2 + mb2) + 21*mt2*Lt + 21*mb2*Lb - 9*(mt2 - mb2)*Lt2
                 - 18*mb2*Lt*Lb - 18*(mt2 - mb2)*li);
}

/**
 * 2-loop self-energy \f$\Sigma(0)\f$ of the top-bottom doublet with
 * \f$G^\pm\f$ exchange (\f$m_G = 0\f$), without the factor
 * \f$1/(16\pi^2)^2\f$.
 *
 * @param yt top Yukawa coupling
 * @param yb bottom Yukawa coupling
 * @param mt top mass
 * @param mb bottom mass
 * @param scale renormalization scale
 */
double self_energy_charged(
   double yt, double yb, double mt, double mb, double scale) noexcept
{
   const double mt2 = sqr(mt);
   const double mb2 = sqr(mb);
   const double d = mt2 - mb2;
   const double Q2 = sqr(scale);
   const double Lt = std::log(mt2/Q2);
   const double Lb = std::log(mb2/Q2);
   const double Lt2 = sqr(Lt);
   const double li = Li2(1 - mb2/mt2);
   const double yt2 = sqr(yt);
   const double yb2 = sqr(yb);

   return sqr(yt2)*(9*(mt2 + mb2) + 9*mt2*Lt2 - 3*mt2*(5*mt2 - 9*mb2)/d*Lt
                    + 6*mb2*(mb2 - 3*mt2)/d*Lb + 18*mt2*li)
      + sqr(yb2)*(9*(mt2 + mb2) - 9*mb2*Lt2 + 18*mb2*Lt*Lb
                  - 6*mt2*(mt2 - 3*mb2)/d*Lt - 3*mb2*(9*mt2 - 5*mb2)/d*Lb
                  - 18*mb2*li)
      + yt2*yb2*(-24*(mt2 + mb2) - 9*d*Lt2 - 18*mb2*Lt*Lb
                 + 9*mt2*(3*mt2 + 5*mb2)/d*Lt - 9*mb2*(5*mt2 + 3*mb2)/d*Lb
                 - 18*d*li)
      + yt2*yt*yb*mt*mb*(18 - 18*Lt2 - 12*(mt2 + 3*mb2)/d*Lt + 48*mb2/d*Lb
                         - 36*li)
      + yt*yb2*yb*mt*mb*(18 + 18*Lt2 - 36*Lt*Lb - 48*mt2/d*Lt
                         + 12*(3*mt2 + mb2)/d*Lb + 36*li);
}

} // anonymous namespace

/**
 * Standard Model Higgs self-energy 2-loop, \f$O((\alpha_b + \alpha_t)^2)\f$.
 *
 * Closed form of the generic 2-loop diagrams [arxiv:1503.03098] with
 * one closed top or bottom loop and \f$h\f$, \f$G^0\f$ or \f$G^\pm\f$
 * exchange, with the full dependence on \f$m_t\f$, \f$m_b\f$ and
 * \f$m_h\f$ (MS-bar, Landau gauge, gaugeless limit, \f$m_G = 0\f$,
 * \f$p^2 = 0\f$).  The bottom Yukawa coupling is \f$y_b = y_t
 * m_b/m_t\f$.  Diagrams with the trilinear Higgs coupling are not
 * included.
 *
 * @param p2    squared momentum (not used)
 * @param scale renormalization scale
 * @param mt MS-bar top mass
 * @param yt MS-bar Yukawa coupling
 * @param mb MS-bar bottom mass
 * @param mh tree-level Higgs mass in the gaugeless limit,
 *           \f$\sqrt{\lambda} v\f$
 *
 * @return real part of 2-loop self-energy \f$O((\alpha_b + \alpha_t)^2)\f$
 */
double self_energy_higgs_2loop_at_at_sm(
   double /* p2 */, double scale, double mt, double yt, double mb, double mh)
{
   const double yb = yt*mb/mt;

   const double result =
      self_energy_neutral(yt, mt, mh, scale)
      + self_energy_neutral(yb, mb, mh, scale)
      + self_energy_charged(yt, yb, mt, mb, scale);

   return result * twoLoop;
}

/**
 * Standard Model Higgs tadpole 2-loop, \f$O((\alpha_b + \alpha_t)^2)\f$,
 * divided by the VEV.
 *
 * Closed form of the generic 2-loop diagrams, see
 * self_energy_higgs_2loop_at_at_sm().
 *
 * @param scale renormalization scale
 * @param mt MS-bar top mass
 * @param yt MS-bar Yukawa coupling
 * @param mb MS-bar bottom mass
 * @param mh tree-level Higgs mass in the gaugeless limit,
 *           \f$\sqrt{\lambda} v\f$
 *
 * @return real part of 2-loop tadpole \f$O((\alpha_b + \alpha_t)^2)\f$
 */
double tadpole_higgs_2loop_at_at_sm(
   double scale, double mt, double yt, double mb, double mh)
{
   const double yb = yt*mb/mt;

   const double result =
      tadpole_neutral(yt, mt, mh, scale)
      + tadpole_neutral(yb, mb, mh, scale)
      + tadpole_charged(yt, yb, mt, mb, scale);

   return result * twoLoop;
}

double delta_mh_2loop_at_at_sm(
   double p2, double scale, double mt, double yt, double mb, double mh)
{
   return - self_energy_higgs_2loop_at_at_sm(p2, scale, mt, yt, mb, mh)
      + tadpole_higgs_2loop_at_at_sm(scale, mt, yt, mb, mh);
}

/**
 * Standard Model Higgs self-energy 2-loop, \f$O(\alpha_\tau^2)\f$.
 *
 * Closed form of the generic 2-loop diagrams with one closed
 * \f$\tau\f$ loop, see self_energy_higgs_2loop_at_at_sm().  It is the
 * top-quark result with \f$N_c = 1\f$ and a massless neutrino.
 *
 * @param p2    squared momentum (not used)
 * @param scale renormalization scale
 * @param mtau MS-bar tau mass
 * @param ytau MS-bar Yukawa coupling
 * @param mh tree-level Higgs mass in the gaugeless limit,
 *           \f$\sqrt{\lambda} v\f$
 *
 * @return real part of 2-loop self-energy \f$O(\alpha_\tau^2) \f$
 */
double self_energy_higgs_2loop_atau_atau_sm(
   double /* p2 */, double scale, double mtau, double ytau, double mh)
{
   const double mtau2 = sqr(mtau);
   const double L = std::log(mtau2/sqr(scale));
   // G^+- exchange, limit m_b -> 0 of self_energy_charged()
   const double charged = pow4(ytau)*mtau2*(9 + 3*sqr(Pi) - 15*L + 9*sqr(L));

   const double result =
      (self_energy_neutral(ytau, mtau, mh, scale) + charged)/3;

   return result * twoLoop;
}

/**
 * Standard Model Higgs tadpole 2-loop, \f$O(\alpha_\tau^2)\f$, divided
 * by the VEV.
 *
 * Closed form of the generic 2-loop diagrams, see
 * self_energy_higgs_2loop_atau_atau_sm().
 *
 * @param scale renormalization scale
 * @param mtau MS-bar tau mass
 * @param ytau MS-bar Yukawa coupling
 * @param mh tree-level Higgs mass in the gaugeless limit,
 *           \f$\sqrt{\lambda} v\f$
 *
 * @return real part of 2-loop tadpole \f$O(\alpha_\tau^2) \f$
 */
double tadpole_higgs_2loop_atau_atau_sm(
   double scale, double mtau, double ytau, double mh)
{
   const double mtau2 = sqr(mtau);
   const double L = std::log(mtau2/sqr(scale));
   // G^+- exchange, limit m_b -> 0 of tadpole_charged()
   const double charged = pow4(ytau)*mtau2*(9 + sqr(Pi) - 9*L + 3*sqr(L));

   const double result = (tadpole_neutral(ytau, mtau, mh, scale) + charged)/3;

   return result * twoLoop;
}

double delta_mh_2loop_atau_atau_sm(
   double p2, double scale, double mtau, double ytau, double mh)
{
   return - self_energy_higgs_2loop_atau_atau_sm(p2, scale, mtau, ytau, mh)
      + tadpole_higgs_2loop_atau_atau_sm(scale, mtau, ytau, mh);
}

namespace {

double QA0(double m, double Q) {
   return flexiblesusy::Loop_library::get().A0(m*m, Q*Q).real();
}

double QB0(double p, double m1, double m2, double Q) {
   return flexiblesusy::Loop_library::get().B0(p*p, m1*m1, m2*m2, Q*Q).real();
}

} // anonymous namespace

/**
 * Standard Model Higgs 1-loop contribution as used in SUSYHD 1.0.2.
 *
 * @note The result contains the 1-loop tadpole diagrams.  It is
 * therefore not 1-particle irreducible (1PI).
 *
 * @param vev Higgs Vacuum expectation value (~ 246 GeV)
 * @param Mt top quark pole mass
 * @param mh MS-bar Higgs mass
 * @param MW W pole mass
 * @param MZ Z pole mass
 * @param Q renormalization scale
 *
 * @return real part of 1-loop self-energy
 */
double delta_mh_1loop_sm_SUSYHD(
   double vev, double Mt, double mh, double MW, double MZ, double Q)
{
   using namespace std;

   const double Pi = M_PI;

   const double delta_lambda =
      ((pow(mh,4) + pow(mh,2)*
         (-6*pow(Mt,2) + 2*pow(MW,2) + pow(MZ,2)) -
        8*(2*pow(MW,4) + pow(MZ,4)))/4. -
     (3*pow(mh,2)*pow(MW,2)*QA0(mh,Q))/
      (2.*(pow(mh,2) - pow(MW,2))) + 3*pow(mh,2)*QA0(Mt,Q) +
     (pow(mh,2)*(-11 + (3*pow(mh,2))/
           (pow(mh,2) - pow(MW,2)) -
          (3*pow(MW,2))/(-pow(MW,2) + pow(MZ,2)))*QA0(MW,Q))/
      2. + (pow(mh,2)*(7*pow(MW,2) - 4*pow(MZ,2))*QA0(MZ,Q))/
      (2.*(-pow(MW,2) + pow(MZ,2))) +
     (9*pow(mh,4)*QB0(mh,mh,mh,Q))/4. +
     3*pow(Mt,2)*(pow(mh,2) - 4*pow(Mt,2))*QB0(mh,Mt,Mt,Q) +
     ((pow(mh,4) - 4*pow(mh,2)*pow(MW,2) + 12*pow(MW,4))*
        QB0(mh,MW,MW,Q))/2. +
     ((pow(mh,4) - 4*pow(mh,2)*pow(MZ,2) + 12*pow(MZ,4))*
      QB0(mh,MZ,MZ,Q))/4.)/(16.*pow(Pi,2)*pow(vev,4));

   const double sigma = 2 * delta_lambda * vev * vev;

   return sigma;
}

/**
 * Standard Model Higgs 2-loop contribution as used in SUSYHD 1.0.2.
 *
 * @note The result contains the 2-loop tadpole diagrams.  It is
 * therefore not 1-particle irreducible (1PI).
 *
 * @param vev Higgs Vacuum expectation value (~ 246 GeV)
 * @param Mt top quark pole mass
 * @param Mh Higgs pole mass
 * @param g3 strong gauge coupling at the top quark pole mass scale
 *    \f$Q = M_t\f$
 *
 * @return real part of 2-loop self-energy
 */
double delta_mh_2loop_sm_SUSYHD(
   double vev, double Mt, double Mh, double g3)
{
   using namespace std;

   const double Pi = M_PI;

   const double delta_lambda_QCD = pow(g3,2)*(
      -23.88 + 0.12*(-125 + Mh) - 0.64*(-173 + Mt))/(256.*pow(Pi,4));
   const double delta_lambda_EW =
      (-9.45 - 0.12*(-125 + Mh) - 0.21*(-173 + Mt))/(256.*pow(Pi,4));
   const double delta_lambda = delta_lambda_QCD + delta_lambda_EW;

   const double sigma = 2 * delta_lambda * vev * vev;

   return sigma;
}

} // namespace sm_twoloophiggs
} // namespace flexiblesusy
