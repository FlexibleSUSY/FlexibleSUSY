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
#include "loop_libraries/loop_library.hpp"

#include <cmath>

namespace flexiblesusy {
namespace sm_twoloophiggs {

namespace {

constexpr double Pi      = 3.141592653589793;
constexpr double oneLoop = 6.332573977646110963e-03; // 1/(4 Pi)^2
constexpr double twoLoop = 4.010149318236068752e-05; // 1/(4 Pi)^4

double sqr(double x) noexcept { return x*x; }
double pow3(double x) noexcept { return x*x*x; }
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

/**
 * Standard Model Higgs self-energy 2-loop, \f$O((\alpha_b + \alpha_t)^2)\f$.
 *
 * @param p2    squared momentum (not used so far)
 * @param scale renormalization scale
 * @param mt MS-bar top mass
 * @param yt MS-bar Yukawa coupling
 * @param mb MS-bar bottom mass
 *
 * @return real part of 2-loop self-energy \f$O((\alpha_b + \alpha_t)^2)\f$
 */
double self_energy_higgs_2loop_at_at_sm(
   double /* p2 */, double scale, double mt, double yt, double mb)
{
   const double Pi2 = sqr(Pi);
   const double yt4 = pow4(yt);
   const double mt2 = sqr(mt);
   const double mb2 = sqr(mb);
   const double Q2 = sqr(scale);
   const double r = mb2/mt2;
   const double LogT = std::log(mt2 / Q2);
   const double LogT2 = sqr(LogT);
   const double LogTB = 0.5*std::log(r);
   const double LogTB2 = sqr(LogTB);

   const double result =
      3 * mt2 * yt4 * (19 + Pi2 - 27 * LogT + 9 * LogT2)
      - 3 * mt2 * yt4 * (5 + 3*Pi2 - 3*LogT + 9*LogT2)*r
      + 3./2. * mt2 * yt4 * (35 + 6*Pi2 + 6*LogT - 18*LogT2
                             - 24*LogTB*(1 + 3*LogT))*sqr(r)
      - 0.5 * mt2 * yt4 * (-29 + 6*Pi2 - 144*LogTB2 + 162*LogT - 54*LogT2
                           - 24*LogTB*(-8 + 9*LogT))*pow3(r);

   return result * twoLoop;
}

/**
 * Standard Model Higgs tadpole 2-loop, \f$O((\alpha_b + \alpha_t)^2)\f$.
 *
 * @param scale renormalization scale
 * @param mt MS-bar top mass
 * @param yt MS-bar Yukawa coupling
 * @param mb MS-bar bottom mass
 *
 * @return real part of 2-loop self-energy \f$O((\alpha_b + \alpha_t)^2)\f$
 */
double tadpole_higgs_2loop_at_at_sm(
   double scale, double mt, double yt, double mb)
{
   const double Pi2 = sqr(Pi);
   const double yt4 = pow4(yt);
   const double mt2 = sqr(mt);
   const double mb2 = sqr(mb);
   const double Q2 = sqr(scale);
   const double r = mb2/mt2;
   const double LogT = std::log(mt2 / Q2);
   const double LogT2 = sqr(LogT);
   const double LogTB = 0.5*std::log(r);
   const double LogTB2 = sqr(LogTB);

   const double result =
      mt2 * yt4 * (45 + Pi2 - 39*LogT + 9*LogT2)
      - 3 * mt2 * yt4 * (5 + Pi2 - 5*LogT + 3*LogT2)*r
      + 3./2. * mt2 * yt4 * (5 + 2*Pi2 + LogTB*(8 - 24*LogT) + 10*LogT
                             - 6*LogT2)*sqr(r)
      - 1./6. * mt2 * yt4 * (-185 + 6*Pi2 - 144*LogTB2 + 234*LogT - 54*LogT2
                             - 24*LogTB*(-14 + 9*LogT))*pow3(r);

   return result * twoLoop;
}

double delta_mh_2loop_at_at_sm(
   double p2, double scale, double mt, double yt, double mb)
{
   return - self_energy_higgs_2loop_at_at_sm(p2, scale, mt, yt, mb)
      + tadpole_higgs_2loop_at_at_sm(scale, mt, yt, mb);
}

/**
 * Standard Model Higgs self-energy 2-loop, \f$O(\alpha_tau^2)\f$.
 *
 * @param p2    squared momentum (not used so far)
 * @param scale renormalization scale
 * @param mtau MS-bar tau mass
 * @param ytau MS-bar Yukawa coupling
 *
 * @return real part of 2-loop self-energy \f$O(\alpha_tau^2) \f$
 */
double self_energy_higgs_2loop_atau_atau_sm(
   double /* p2 */, double scale, double mtau, double ytau)
{
   const double ytau4 = pow4(ytau);
   const double mtau2 = sqr(mtau);
   const double Q2 = sqr(scale);
   const double LogT = std::log(mtau2 / Q2);
   const double LogT2 = sqr(LogT);

   const double result =
      mtau2 * ytau4 * (19 + sqr(Pi) - 27 * LogT + 9 * LogT2);

   return result * twoLoop;
}

/**
 * Standard Model Higgs tadpole 2-loop, \f$O(\alpha_tau^2)\f$.
 *
 * @param scale renormalization scale
 * @param mtau MS-bar tau mass
 * @param ytau MS-bar Yukawa coupling
 *
 * @return real part of 2-loop self-energy \f$O(\alpha_tau^2) \f$
 */
double tadpole_higgs_2loop_atau_atau_sm(
   double scale, double mtau, double ytau)
{
   const double ytau4 = pow4(ytau);
   const double mtau2 = sqr(mtau);
   const double Q2 = sqr(scale);
   const double LogT = std::log(mtau2 / Q2);
   const double LogT2 = sqr(LogT);

   const double result =
      1./3. * mtau2 * ytau4 * (45 + sqr(Pi) - 39*LogT + 9*LogT2);

   return result * twoLoop;
}

double delta_mh_2loop_atau_atau_sm(
   double p2, double scale, double mtau, double ytau)
{
   return - self_energy_higgs_2loop_atau_atau_sm(p2, scale, mtau, ytau)
      + tadpole_higgs_2loop_atau_atau_sm(scale, mtau, ytau);
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
