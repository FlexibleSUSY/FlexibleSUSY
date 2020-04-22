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

// This file has been generated at Fri 17 Apr 2020 19:33:32
// with the script "as2_to_cpp.m".

#include "mssm_twoloop_as.hpp"
#include "dilog.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <ostream>

namespace flexiblesusy {
namespace mssm_twoloop_as {

namespace {
   const Real Pi = 3.1415926535897932384626433832795l;

   template <typename T> T power2(T x)  noexcept { return x*x; }
   template <typename T> T power3(T x)  noexcept { return x*x*x; }
   template <typename T> T power4(T x)  noexcept { return power2(power2(x)); }
   template <typename T> T power5(T x)  noexcept { return x*power4(x); }
   template <typename T> T power6(T x)  noexcept { return power2(power3(x)); }
   template <typename T> T power7(T x)  noexcept { return x*power6(x); }
   template <typename T> T power8(T x)  noexcept { return power2(power4(x)); }
   template <typename T> T power10(T x) noexcept { return power2(power5(x)); }
   template <typename T> T power12(T x) noexcept { return power3(power4(x)); }
   template <typename T> T power14(T x) noexcept { return power2(x)*power12(x); }

   const Real oneLoop = 1/power2(4*Pi);
   const Real twoLoop = power2(oneLoop);

   template <typename T>
   bool is_zero(T a, T prec = std::numeric_limits<T>::epsilon())
   {
      return std::abs(a) < prec;
   }

   template <typename T>
   bool is_equal(T a, T b, T prec = std::numeric_limits<T>::epsilon())
   {
      return is_zero(a - b, prec);
   }

   template <typename T>
   bool is_equal_rel(T a, T b, T prec = std::numeric_limits<T>::epsilon())
   {
      if (is_equal(a, b, std::numeric_limits<T>::epsilon())) {
         return true;
      }

      if (std::abs(a) < std::numeric_limits<T>::epsilon() ||
          std::abs(b) < std::numeric_limits<T>::epsilon()) {
         return false;
      }

      return std::abs((a - b)/a) < prec;
   }

   Real LambdaSquared(Real x, Real y) noexcept
   {
      return power2(1 - x - y) - 4*x*y;
   }

   /// ClausenCl[2,x]
   Real ClausenCl2(Real x) noexcept
   {
      const std::complex<Real> img(0.0l, 1.0l);

      return std::imag(dilog(std::exp(img*x)));
   }

   /// x < 1 && y < 1, LambdaSquared(x,y) > 0
   Real PhiPos(Real x, Real y) noexcept
   {
      const Real lambda = std::sqrt(LambdaSquared(x,y));

      return (-(std::log(x)*std::log(y))
              + 2*std::log((1 - lambda + x - y)/2)*std::log((1 - lambda - x + y)/2)
              - 2*dilog((1 - lambda + x - y)/2)
              - 2*dilog((1 - lambda - x + y)/2)
              + power2(Pi)/3)/lambda;
   }

   /// LambdaSquared(x,y) < 0
   Real PhiNeg(Real x, Real y) noexcept
   {
      const Real lambda = std::sqrt(-LambdaSquared(x,y));

      return 2*(+ ClausenCl2(2*std::acos((1 + x - y)/(2*std::sqrt(x))))
                + ClausenCl2(2*std::acos((1 - x + y)/(2*std::sqrt(y))))
                + ClausenCl2(2*std::acos((-1 + x + y)/(2*std::sqrt(x*y)))))/lambda;
   }

   Real Phi(Real x, Real y) noexcept
   {
      const Real lambda = LambdaSquared(x,y);

      if (lambda > 0) {
         return PhiPos(x,y);
      }

      return PhiNeg(x,y);
   }

   /**
    * Fin3[] function from twoloopbubble.m .
    *
    * @param mm1 squared mass \f$m_1^2\f$
    * @param mm2 squared mass \f$m_2^2\f$
    * @param mm3 squared mass \f$m_3^2\f$
    * @param mmu squared renormalization scale
    *
    * @return Fin3(m12, m22, m32, mmu)
    */
   Real Fin3(Real mm1, Real mm2, Real mm3, Real mmu) noexcept
   {
      std::array<Real,3> masses = { mm1, mm2, mm3 };
      std::sort(masses.begin(), masses.end());

      const Real mm = masses[2];
      const Real x = masses[0]/mm;
      const Real y = masses[1]/mm;
      const Real lambda = LambdaSquared(x,y);
      const Real logx = std::log(x);
      const Real logy = std::log(y);
      const Real logm = std::log(mm/mmu);

      if (is_zero(lambda, 1e-10l)) {
         return -(mm*(2*y*(-3 + 2*logm)*logy
                      + logx*(2*x*(-3 + 2*logm) + (-1 + x + y)*logy)
                      + (1 + x + y)*(7 - 6*logm + power2(Pi)/6 + 2*power2(logm))
                      + x*power2(logx) + y*power2(logy)))/2;
      }

      return mm*((-7 + 6*logm + logx*logy
                  - lambda*Phi(x,y) - power2(Pi)/6 - 2*power2(logm))/2
                 - (x*(7 - 6*logm + logx*(-6 + 4*logm + logy)
                       + power2(Pi)/6 + 2*power2(logm) + power2(logx)))/2
                 - (y*(7 - 6*logm + (
                     -6 + 4*logm + logx)*logy + power2(Pi)/6
                       + 2*power2(logm) + power2(logy)))/2);
   }

   /// Delta[m1,m2,m3,-1]
   Real DeltaInv(Real m1, Real m2, Real m3) noexcept
   {
      return 1/(power2(m1) + power2(m2) + power2(m3) - 2*(m1*m2 + m1*m3 + m2*m3));
   }

   /// calculates sin(theta)
   Real calc_sin_theta(Real mf, Real xf, Real msf12, Real msf22) noexcept
   {
      if (is_zero(mf, 1e-10l) || is_zero(xf, 1e-10l)) {
         return 0;
      }

      const Real sin_2theta = 2*mf*xf / (msf12 - msf22);
      const Real theta = 0.5l*std::asin(sin_2theta);

      return std::sin(theta);
   }

   /// calculates Higgs mixing angle from squarde Higgs masses and tan(beta)
   Real calc_alpha(Real mh2, Real mH2, Real tb) noexcept
   {
      const Real beta = std::atan(tb);
      const Real sin_2alpha = -(mH2 + mh2)/(mH2 - mh2) * std::sin(2*beta);

      return 0.5l*std::asin(sin_2alpha);
   }

   /// shift gluino mass away from mst1 and mst2 if too close
   Real shift_mg(Real mg, Real mst1, Real mst2) noexcept
   {
      if (is_equal_rel(std::min(mst1, mst2), mg, 0.01l)) {
         return mg * 0.98l;
      }

      if (is_equal_rel(std::max(mst1, mst2), mg, 0.01l)) {
         return mg * 1.02l;
      }

      return mg;
   }

} // anonymous namespace

/**
 * @brief 2-loop \f$O(\alpha_s^2)\f$ contributions to \f$\Delta g_3^{(2L)}\f$ [hep-ph/0509048, arXiv:0810.5101]
 *
 * The function returns \f$\Delta g_3^{(2L)}\f$ which appears in
 * conversion of the strong \f$\overline{\text{MS}}\f$ gauge coupling
 * of the Standard Model with 5 active quark flavours,
 * \f$g_3^{\text{SM(5)},\overline{\text{MS}}}\f$, to the strong
 * \f$\overline{\text{DR}}\f$ gauge coupling of the full MSSM
 * \f$g_3^{\text{MSSM},\overline{\text{DR}}}\f$,
 * \f{align*}{
  g_3^{\text{SM(5)},\overline{\text{MS}}} =
  g_3^{\text{MSSM},\overline{\text{DR}}} \left[
     1 + \Delta g_3^{(1L)} + \Delta g_3^{(2L)}
  \right]
 * \f}
 */
Real delta_alpha_s_2loop_as_as(const Parameters& pars)
{
   const Real g3        = pars.g3;
   const Real xt        = pars.xt;
   const Real xb        = pars.xb;
   const Real mt        = pars.mt;
   const Real mt2       = power2(pars.mt);
   const Real mb        = pars.mb;
   const Real mg        = shift_mg(pars.mg, pars.mst1, pars.mst2);
   const Real mg2       = power2(mg);
   const Real mg4       = power2(mg2);
   const Real mg6       = mg2*mg4;
   const Real mst12     = power2(pars.mst1);
   const Real mst14     = power2(mst12);
   const Real mst16     = mst12*mst14;
   const Real mst22     = power2(pars.mst2);
   const Real mst24     = power2(mst22);
   const Real mst26     = mst22*mst24;
   const Real msb12     = power2(pars.msb1);
   const Real msb22     = power2(pars.msb2);
   const Real msd12     = power2(pars.msd1);
   const Real msd14     = power2(msd12);
   const Real msd16     = msd12*msd14;
   const Real msd22     = power2(pars.msd2);
   const Real msd24     = power2(msd22);
   const Real msd26     = msd22*msd24;
   const Real Q2        = power2(pars.Q);
   const Real snt       = calc_sin_theta(mt, xt, mst12, mst22);
   const Real snb       = calc_sin_theta(mb, xb, msb12, msb22);
   const Real invdmst   = 1/(mst12 - mst22);
   const Real invdmsb1g = 1/(msb12 - mg2);
   const Real invdmsb2g = 1/(msb22 - mg2);
   const Real invdmsd1g = 1/(msd12 - mg2);
   const Real invdmsd2g = 1/(msd22 - mg2);
   const Real lmst12    = std::log(mst12/Q2);
   const Real lmst22    = std::log(mst22/Q2);
   const Real lmsb12    = std::log(msb12/Q2);
   const Real lmsb22    = std::log(msb22/Q2);
   const Real lmsd12    = std::log(msd12/Q2);
   const Real lmsd22    = std::log(msd22/Q2);
   const Real lmg2      = std::log(mg2/Q2);
   const Real lmt2      = std::log(mt2/Q2);

   const Real result =
   (invdmst*xt*(8*lmg2*mg2*mg4*(mst12 - mst22) + 36*lmst12*mst12*mst14*mst22 +
     54*lmg2*mg2*mst12*mst22*((2 + lmst12 - lmt2)*mst12 + (-2 - lmst22 + lmt2)*
     mst22) - 36*lmst22*mst12*mst22*mst24 + 8*lmg2*mg4*(mst12 - mst22)*mt2 - 54
     *mg2*mst12*(mst12 - mst22)*mst22*power2(lmt2) - 324*mg2*mst22*power2(mst12
     ) + 106*lmst12*mg2*mst22*power2(mst12) + 36*lmst12*mst22*mt2*power2(mst12)
     + 324*mg2*mst12*power2(mst22) - 106*lmst22*mg2*mst12*power2(mst22) - 36*
     lmst22*mst12*mt2*power2(mst22) + 4*lmt2*(2*mg4*mst22*mt2 - 9*mst22*mt2*
     power2(mst12) + mst12*(-9*mst14*mst22 + 9*mst22*mst24 - 2*mg4*mt2 + 9*mt2*
     power2(mst22))) + 8*mg2*mst12*power2(mt2) - 8*mg2*mst22*power2(mt2) - 2*
     lmt2*mg2*(4*mg4*(mst12 - mst22) + (-55 + 27*lmst12)*mst22*power2(mst12) +
     (55 - 27*lmst22)*mst12*power2(mst22) + 4*mst12*power2(mt2) - 4*mst22*
     power2(mt2)) - 9*mg2*mst22*power2(mst12)*power2(Pi) + 9*mg2*mst12*power2(
     mst22)*power2(Pi)))/(9.*mg*mst12*mst22*mt2) + DeltaInv(mg2,mt2,mst12)*((6*
     mst12*(44*mg4*mst12 + 9*mg2*(4*mst14 + mst12*mt2))*power2(lmt2) + 54*mg4*
     power2(lmg2)*power2(mst12) + mst12*(26*mg6 + 1880*mg2*mst14 - 666*lmst12*
     mg2*mst14 + 18*mst16 + 9*lmst12*mst16 - 42*mg4*mt2 + 262*mg2*mst12*mt2 -
     125*lmst12*mg2*mst12*mt2 - 18*mst14*mt2 - 27*lmst12*mst14*mt2 + 54*mg2*
     mst14*power2(lmst12) + 45*mg2*mst14*power2(Pi) + 9*mg2*mst12*mt2*power2(Pi
     ) + mg4*mst12*(2192 - 195*lmst12 + 53*power2(Pi))) - lmg2*(mg6*(23*mst12 -
     8*mt2) + mg4*mst12*((722 + 210*lmst12 - 318*lmt2)*mst12 + 43*mt2) + 9*mg2*
     mst12*((11 + 18*lmst12 - 18*lmt2)*mst14 + (13 + 6*lmst12 - 6*lmt2)*mst12*
     mt2) + 8*power4(mg2)) + lmt2*(mg6*(23*mst12 - 8*mt2) + mg2*mst12*(45*(-19
     + 6*lmst12)*mst14 + 2*(-41 + 27*lmst12)*mst12*mt2) + mst12*(-991*mg4*mst12
      + 210*lmst12*mg4*mst12 - 9*mst16 + 43*mg4*mt2 + 27*mst14*mt2) + 8*power4(
     mg2)))/(18.*mst12*mt2) - (invdmst*xt*(6*mg4*mst12*(35*mg2*mst12 + 27*mst14
      + 9*mst12*mt2)*power2(lmg2) + 6*mst12*(86*mg2*mg4*mst12 + 214*mg4*mst14 +
     36*mg2*mst16 + 53*mg4*mst12*mt2 + 45*mg2*mst14*mt2)*power2(lmt2) + mst12*(
     210*mg4*mst14*power2(lmst12) + mg2*(34*mg6 + 2700*mst16 - 1046*lmst12*
     mst16 - 50*mg4*mt2 + 2178*mst14*mt2 - 772*lmst12*mst14*mt2 + 162*mst16*
     power2(lmst12) + 54*mst14*mt2*power2(lmst12) + 63*mst16*power2(Pi) + 54*
     mst14*mt2*power2(Pi) + mg4*mst12*(5068 - 146*lmst12 + 121*power2(Pi))) + 2
     *mg4*(mst12*mt2*(1196 + 31*power2(Pi)) + mst14*(5759 + 138*power2(Pi))) +
     4*lmst12*(-533*mg4*mst14 - 87*mg4*mst12*mt2 - 9*mst16*mt2 + 9*power4(mst12
     ))) - 2*lmg2*(3*mg4*mst12*((263 + 152*lmst12 - 206*lmt2)*mst14 + 2*(72 +
     22*lmst12 - 31*lmt2)*mst12*mt2) + mg2*(-(mg6*(5*mst12 + 4*mt2)) + mst12*((
     910 + 153*lmst12 - 363*lmt2)*mg4*mst12 + 27*(-2 + lmst12 - lmt2)*mst16 +
     42*mg4*mt2 + 108*(1 + lmst12 - lmt2)*mst14*mt2)) + 4*power5(mg2)) + 2*lmt2
     *(mg2*(-(mg6*(5*mst12 + 4*mt2)) + mst12*((-1195 + 153*lmst12)*mg4*mst12 +
     7*(-95 + 27*lmst12)*mst16 + 42*mg4*mt2 + 2*(-239 + 81*lmst12)*mst14*mt2))
     + mst12*((-3113 + 666*lmst12)*mg4*mst14 + 6*(-85 + 22*lmst12)*mg4*mst12*
     mt2 + 18*mst16*mt2 - 18*power4(mst12)) + 4*power5(mg2))))/(9.*mg*mst12*mt2
     )) + DeltaInv(mg2,mst22,mt2)*((6*mst22*(44*mg4*mst22 + 9*mg2*(4*mst24 +
     mst22*mt2))*power2(lmt2) + 54*mg4*power2(lmg2)*power2(mst22) + mst22*(26*
     mg6 + 1880*mg2*mst24 - 666*lmst22*mg2*mst24 + 18*mst26 + 9*lmst22*mst26 -
     42*mg4*mt2 + 262*mg2*mst22*mt2 - 125*lmst22*mg2*mst22*mt2 - 18*mst24*mt2 -
     27*lmst22*mst24*mt2 + 54*mg2*mst24*power2(lmst22) + 45*mg2*mst24*power2(Pi
     ) + 9*mg2*mst22*mt2*power2(Pi) + mg4*mst22*(2192 - 195*lmst22 + 53*power2(
     Pi))) - lmg2*(mg6*(23*mst22 - 8*mt2) + mg4*mst22*((722 + 210*lmst22 - 318*
     lmt2)*mst22 + 43*mt2) + 9*mg2*mst22*((11 + 18*lmst22 - 18*lmt2)*mst24 + (
     13 + 6*lmst22 - 6*lmt2)*mst22*mt2) + 8*power4(mg2)) + lmt2*(mg6*(23*mst22
     - 8*mt2) + mg2*mst22*(45*(-19 + 6*lmst22)*mst24 + 2*(-41 + 27*lmst22)*
     mst22*mt2) + mst22*(-991*mg4*mst22 + 210*lmst22*mg4*mst22 - 9*mst26 + 43*
     mg4*mt2 + 27*mst24*mt2) + 8*power4(mg2)))/(18.*mst22*mt2) + (invdmst*xt*(6
     *mg4*mst22*(35*mg2*mst22 + 27*mst24 + 9*mst22*mt2)*power2(lmg2) + 6*mst22*
     (86*mg2*mg4*mst22 + 214*mg4*mst24 + 36*mg2*mst26 + 53*mg4*mst22*mt2 + 45*
     mg2*mst24*mt2)*power2(lmt2) + mst22*(210*mg4*mst24*power2(lmst22) + mg2*(
     34*mg6 + 2700*mst26 - 1046*lmst22*mst26 - 50*mg4*mt2 + 2178*mst24*mt2 -
     772*lmst22*mst24*mt2 + 162*mst26*power2(lmst22) + 54*mst24*mt2*power2(
     lmst22) + 63*mst26*power2(Pi) + 54*mst24*mt2*power2(Pi) + mg4*mst22*(5068
     - 146*lmst22 + 121*power2(Pi))) + 2*mg4*(mst22*mt2*(1196 + 31*power2(Pi))
     + mst24*(5759 + 138*power2(Pi))) + 4*lmst22*(-533*mg4*mst24 - 87*mg4*mst22
     *mt2 - 9*mst26*mt2 + 9*power4(mst22))) - 2*lmg2*(3*mg4*mst22*((263 + 152*
     lmst22 - 206*lmt2)*mst24 + 2*(72 + 22*lmst22 - 31*lmt2)*mst22*mt2) + mg2*(
     -(mg6*(5*mst22 + 4*mt2)) + mst22*((910 + 153*lmst22 - 363*lmt2)*mg4*mst22
     + 27*(-2 + lmst22 - lmt2)*mst26 + 42*mg4*mt2 + 108*(1 + lmst22 - lmt2)*
     mst24*mt2)) + 4*power5(mg2)) + 2*lmt2*(mg2*(-(mg6*(5*mst22 + 4*mt2)) +
     mst22*((-1195 + 153*lmst22)*mg4*mst22 + 7*(-95 + 27*lmst22)*mst26 + 42*mg4
     *mt2 + 2*(-239 + 81*lmst22)*mst24*mt2)) + mst22*((-3113 + 666*lmst22)*mg4*
     mst24 + 6*(-85 + 22*lmst22)*mg4*mst22*mt2 + 18*mst26*mt2 - 18*power4(mst22
     )) + 4*power5(mg2))))/(9.*mg*mst22*mt2)) + Fin3(mg2,mt2,mst12,Q2)*(((6*mg2
     *mst12)/mt2 - (4*invdmst*mg*(35*mg2*mst12 + 27*mst14 + 9*mst12*mt2)*xt)/(
     3.*mt2))*DeltaInv(mg2,mt2,mst12) + power2(DeltaInv(mg2,mt2,mst12))*((2*(33
     *mg6*mst12 + 2*mg4*(67*mst14 + 22*mst12*mt2) + 9*mg2*(mst16 + 4*mst14*mt2)
     ))/(3.*mt2) - (4*invdmst*mg*xt*(15*mg6*mst12 + 255*mg4*mst14 + 109*mg2*
     mst16 + 95*mg4*mst12*mt2 + 196*mg2*mst14*mt2 + 45*mst16*mt2 - 27*power4(
     mst12)))/(3.*mt2)) + power3(DeltaInv(mg2,mt2,mst12))*((4*(86*mg6*mst16 +
     175*mg6*mst14*mt2 + 139*mg4*mst16*mt2 + (-16*mst14 + 29*mst12*mt2)*power4(
     mg2) - 40*mg4*power4(mst12) - 9*mg2*(mst12 - mt2)*power4(mst12) - 21*mst12
     *power5(mg2)))/(3.*mt2) + (8*invdmst*mg*xt*(-142*mg6*mst16 - 269*mg6*mst14
     *mt2 - 339*mg4*mst16*mt2 + (56*mst14 - 33*mst12*mt2)*power4(mg2) + 16*mg4*
     power4(mst12) + 9*mg2*(5*mst12 - 7*mt2)*power4(mst12) + 25*mst12*power5(
     mg2)))/(3.*mt2))) + Fin3(mg2,mst22,mt2,Q2)*(((6*mg2*mst22)/mt2 + (4*
     invdmst*mg*(35*mg2*mst22 + 27*mst24 + 9*mst22*mt2)*xt)/(3.*mt2))*DeltaInv(
     mg2,mst22,mt2) + power2(DeltaInv(mg2,mst22,mt2))*((2*(33*mg6*mst22 + 2*mg4
     *(67*mst24 + 22*mst22*mt2) + 9*mg2*(mst26 + 4*mst24*mt2)))/(3.*mt2) + (4*
     invdmst*mg*xt*(15*mg6*mst22 + 255*mg4*mst24 + 109*mg2*mst26 + 95*mg4*mst22
     *mt2 + 196*mg2*mst24*mt2 + 45*mst26*mt2 - 27*power4(mst22)))/(3.*mt2)) +
     power3(DeltaInv(mg2,mst22,mt2))*((8*invdmst*mg*xt*(142*mg6*mst26 + 269*mg6
     *mst24*mt2 + 339*mg4*mst26*mt2 + (-56*mst24 + 33*mst22*mt2)*power4(mg2) -
     16*mg4*power4(mst22) - 9*mg2*(5*mst22 - 7*mt2)*power4(mst22) - 25*mst22*
     power5(mg2)))/(3.*mt2) + (4*(86*mg6*mst26 + 175*mg6*mst24*mt2 + 139*mg4*
     mst26*mt2 + (-16*mst24 + 29*mst22*mt2)*power4(mg2) - 40*mg4*power4(mst22)
     - 9*mg2*(mst22 - mt2)*power4(mst22) - 21*mst22*power5(mg2)))/(3.*mt2))) +
     power2(DeltaInv(mg2,mt2,mst12))*((23898*mg6*mst14 - 1272*lmst12*mg6*mst14
     + 18858*mg4*mst16 - 5268*lmst12*mg4*mst16 + 6930*mg6*mst12*mt2 - 300*
     lmst12*mg6*mst12*mt2 + 15708*mg4*mst14*mt2 - 2856*lmst12*mg4*mst14*mt2 +
     4914*mg2*mst16*mt2 - 1836*lmst12*mg2*mst16*mt2 + 198*mg6*mst14*power2(
     lmst12) + 804*mg4*mst16*power2(lmst12) + 264*mg4*mst14*mt2*power2(lmst12)
     + 216*mg2*mst16*mt2*power2(lmst12) + 569*mg6*mst14*power2(Pi) + 449*mg4*
     mst16*power2(Pi) + 165*mg6*mst12*mt2*power2(Pi) + 374*mg4*mst14*mt2*power2
     (Pi) + 117*mg2*mst16*mt2*power2(Pi) + 1974*mst12*power4(mg2) + 204*lmst12*
     mst12*power4(mg2) + 47*mst12*power2(Pi)*power4(mg2) + 6*power2(lmg2)*(2*
     mg6*(67*mst14 + 22*mst12*mt2) + 9*mg4*(mst16 + 4*mst14*mt2) + 33*mst12*
     power4(mg2)) - 378*mg2*power4(mst12) + 54*mg2*power2(lmst12)*power4(mst12)
     - 9*mg2*power2(Pi)*power4(mst12) + 6*power2(lmt2)*(402*mg6*mst14 + 306*mg4
     *mst16 + 121*mg6*mst12*mt2 + 294*mg4*mst14*mt2 + 81*mg2*mst16*mt2 + 14*
     mst12*power4(mg2) - 18*mg2*power4(mst12)) - 6*lmt2*((2300 - 301*lmst12)*
     mg6*mst14 + (1788 - 431*lmst12)*mg4*mst16 + (538 - 77*lmst12)*mg6*mst12*
     mt2 + 2*(682 - 151*lmst12)*mg4*mst14*mt2 + ((164 + 19*lmst12)*mst12 - 8*
     mt2)*power4(mg2) + 9*mg2*((42 - 13*lmst12)*mst16*mt2 + (-4 + lmst12)*
     power4(mst12)) + 8*power5(mg2)) + 6*lmg2*(-(mg6*((902 + 235*lmst12 - 503*
     lmt2)*mst14 + (402 + 77*lmst12 - 165*lmt2)*mst12*mt2)) - mg4*((28 + 163*
     lmst12 - 181*lmt2)*mst16 + 2*(202 + 107*lmst12 - 143*lmt2)*mst14*mt2) + ((
     -152 + 19*lmst12 + 47*lmt2)*mst12 - 8*mt2)*power4(mg2) + 9*mg2*((-2 - 5*
     lmst12 + 5*lmt2)*mst16*mt2 + (2 + 3*lmst12 - 3*lmt2)*power4(mst12)) + 8*
     power5(mg2)))/(18.*mt2) - (2*invdmst*mg*xt*(16842*mg6*mst14 - 174*lmst12*
     mg6*mst14 + 27174*mg4*mst16 - 5106*lmst12*mg4*mst16 + 6300*mg6*mst12*mt2 -
     174*lmst12*mg6*mst12*mt2 + 23688*mg4*mst14*mt2 - 2760*lmst12*mg4*mst14*mt2
      + 17472*mg2*mst16*mt2 - 4362*lmst12*mg2*mst16*mt2 + 45*mg6*mst14*power2(
     lmst12) + 765*mg4*mst16*power2(lmst12) + 285*mg4*mst14*mt2*power2(lmst12)
     + 588*mg2*mst16*mt2*power2(lmst12) + 401*mg6*mst14*power2(Pi) + 647*mg4*
     mst16*power2(Pi) + 150*mg6*mst12*mt2*power2(Pi) + 564*mg4*mst14*mt2*power2
     (Pi) + 416*mg2*mst16*mt2*power2(Pi) - 294*mst12*power4(mg2) + 126*lmst12*
     mst12*power4(mg2) - 7*mst12*power2(Pi)*power4(mg2) + 2142*mg2*power4(mst12
     ) - 1722*lmst12*mg2*power4(mst12) + 2268*mt2*power4(mst12) - 864*lmst12*
     mt2*power4(mst12) + 327*mg2*power2(lmst12)*power4(mst12) + 135*mt2*power2(
     lmst12)*power4(mst12) + 51*mg2*power2(Pi)*power4(mst12) + 54*mt2*power2(Pi
     )*power4(mst12) + 3*power2(lmg2)*(255*mg6*mst14 + 109*mg4*mst16 + 95*mg6*
     mst12*mt2 + 196*mg4*mst14*mt2 + 45*mg2*mst16*mt2 + 15*mst12*power4(mg2) -
     27*mg2*power4(mst12)) + 3*power2(lmt2)*(532*mg6*mst14 + 930*mg4*mst16 +
     205*mg6*mst12*mt2 + 837*mg4*mst14*mt2 + 591*mg2*mst16*mt2 - 29*mst12*
     power4(mg2) + 20*mg2*power4(mst12) + 63*mt2*power4(mst12) - 45*power5(
     mst12)) - 1512*power5(mst12) + 540*lmst12*power5(mst12) - 81*power2(lmst12
     )*power5(mst12) - 36*power2(Pi)*power5(mst12) - 6*lmt2*(1556*mg6*mst14 -
     146*lmst12*mg6*mst14 + 2680*mg4*mst16 - 538*lmst12*mg4*mst16 + (492 - 55*
     lmst12)*mg6*mst12*mt2 + 2136*mg4*mst14*mt2 - 368*lmst12*mg4*mst14*mt2 + ((
     -26 + 22*lmst12)*mst12 - 4*mt2)*power4(mg2) + 18*(10 - 3*lmst12)*mt2*
     power4(mst12) + mg2*((1580 - 371*lmst12)*mst16*mt2 + (136 - 78*lmst12)*
     power4(mst12)) + 4*power5(mg2) + 18*(-7 + 2*lmst12)*power5(mst12)) + 6*
     lmg2*(-821*mg6*mst14 + 386*lmt2*mg6*mst14 - 351*mg4*mst16 + 392*lmt2*mg4*
     mst16 - 379*mg6*mst12*mt2 + 150*lmt2*mg6*mst12*mt2 - 788*mg4*mst14*mt2 +
     469*lmt2*mg4*mst14*mt2 + ((-5 + 22*lmst12 - 7*lmt2)*mst12 - 4*mt2)*power4(
     mg2) + 9*lmt2*mt2*power4(mst12) + mg2*((-189 - 175*lmst12 + 220*lmt2)*
     mst16*mt2 + (117 + 31*lmst12 - 58*lmt2)*power4(mst12)) + 4*power5(mg2) - 9
     *lmt2*power5(mst12) + lmst12*(-131*mg6*mst14 - 283*mg4*mst16 - 55*mg6*
     mst12*mt2 - 273*mg4*mst14*mt2 - 9*mt2*power4(mst12) + 9*power5(mst12)))))/
     (9.*mt2)) + power2(DeltaInv(mg2,mst22,mt2))*((23898*mg6*mst24 - 1272*
     lmst22*mg6*mst24 + 18858*mg4*mst26 - 5268*lmst22*mg4*mst26 + 6930*mg6*
     mst22*mt2 - 300*lmst22*mg6*mst22*mt2 + 15708*mg4*mst24*mt2 - 2856*lmst22*
     mg4*mst24*mt2 + 4914*mg2*mst26*mt2 - 1836*lmst22*mg2*mst26*mt2 + 198*mg6*
     mst24*power2(lmst22) + 804*mg4*mst26*power2(lmst22) + 264*mg4*mst24*mt2*
     power2(lmst22) + 216*mg2*mst26*mt2*power2(lmst22) + 569*mg6*mst24*power2(
     Pi) + 449*mg4*mst26*power2(Pi) + 165*mg6*mst22*mt2*power2(Pi) + 374*mg4*
     mst24*mt2*power2(Pi) + 117*mg2*mst26*mt2*power2(Pi) + 1974*mst22*power4(
     mg2) + 204*lmst22*mst22*power4(mg2) + 47*mst22*power2(Pi)*power4(mg2) + 6*
     power2(lmg2)*(2*mg6*(67*mst24 + 22*mst22*mt2) + 9*mg4*(mst26 + 4*mst24*mt2
     ) + 33*mst22*power4(mg2)) - 378*mg2*power4(mst22) + 54*mg2*power2(lmst22)*
     power4(mst22) - 9*mg2*power2(Pi)*power4(mst22) + 6*power2(lmt2)*(402*mg6*
     mst24 + 306*mg4*mst26 + 121*mg6*mst22*mt2 + 294*mg4*mst24*mt2 + 81*mg2*
     mst26*mt2 + 14*mst22*power4(mg2) - 18*mg2*power4(mst22)) - 6*lmt2*((2300 -
     301*lmst22)*mg6*mst24 + (1788 - 431*lmst22)*mg4*mst26 + (538 - 77*lmst22)*
     mg6*mst22*mt2 + 2*(682 - 151*lmst22)*mg4*mst24*mt2 + ((164 + 19*lmst22)*
     mst22 - 8*mt2)*power4(mg2) + 9*mg2*((42 - 13*lmst22)*mst26*mt2 + (-4 +
     lmst22)*power4(mst22)) + 8*power5(mg2)) + 6*lmg2*(-(mg6*((902 + 235*lmst22
      - 503*lmt2)*mst24 + (402 + 77*lmst22 - 165*lmt2)*mst22*mt2)) - mg4*((28 +
     163*lmst22 - 181*lmt2)*mst26 + 2*(202 + 107*lmst22 - 143*lmt2)*mst24*mt2)
     + ((-152 + 19*lmst22 + 47*lmt2)*mst22 - 8*mt2)*power4(mg2) + 9*mg2*((-2 -
     5*lmst22 + 5*lmt2)*mst26*mt2 + (2 + 3*lmst22 - 3*lmt2)*power4(mst22)) + 8*
     power5(mg2)))/(18.*mt2) + (2*invdmst*mg*xt*(16842*mg6*mst24 - 174*lmst22*
     mg6*mst24 + 27174*mg4*mst26 - 5106*lmst22*mg4*mst26 + 6300*mg6*mst22*mt2 -
     174*lmst22*mg6*mst22*mt2 + 23688*mg4*mst24*mt2 - 2760*lmst22*mg4*mst24*mt2
      + 17472*mg2*mst26*mt2 - 4362*lmst22*mg2*mst26*mt2 + 45*mg6*mst24*power2(
     lmst22) + 765*mg4*mst26*power2(lmst22) + 285*mg4*mst24*mt2*power2(lmst22)
     + 588*mg2*mst26*mt2*power2(lmst22) + 401*mg6*mst24*power2(Pi) + 647*mg4*
     mst26*power2(Pi) + 150*mg6*mst22*mt2*power2(Pi) + 564*mg4*mst24*mt2*power2
     (Pi) + 416*mg2*mst26*mt2*power2(Pi) - 294*mst22*power4(mg2) + 126*lmst22*
     mst22*power4(mg2) - 7*mst22*power2(Pi)*power4(mg2) + 2142*mg2*power4(mst22
     ) - 1722*lmst22*mg2*power4(mst22) + 2268*mt2*power4(mst22) - 864*lmst22*
     mt2*power4(mst22) + 327*mg2*power2(lmst22)*power4(mst22) + 135*mt2*power2(
     lmst22)*power4(mst22) + 51*mg2*power2(Pi)*power4(mst22) + 54*mt2*power2(Pi
     )*power4(mst22) + 3*power2(lmg2)*(255*mg6*mst24 + 109*mg4*mst26 + 95*mg6*
     mst22*mt2 + 196*mg4*mst24*mt2 + 45*mg2*mst26*mt2 + 15*mst22*power4(mg2) -
     27*mg2*power4(mst22)) + 3*power2(lmt2)*(532*mg6*mst24 + 930*mg4*mst26 +
     205*mg6*mst22*mt2 + 837*mg4*mst24*mt2 + 591*mg2*mst26*mt2 - 29*mst22*
     power4(mg2) + 20*mg2*power4(mst22) + 63*mt2*power4(mst22) - 45*power5(
     mst22)) - 1512*power5(mst22) + 540*lmst22*power5(mst22) - 81*power2(lmst22
     )*power5(mst22) - 36*power2(Pi)*power5(mst22) - 6*lmt2*(1556*mg6*mst24 -
     146*lmst22*mg6*mst24 + 2680*mg4*mst26 - 538*lmst22*mg4*mst26 + (492 - 55*
     lmst22)*mg6*mst22*mt2 + 2136*mg4*mst24*mt2 - 368*lmst22*mg4*mst24*mt2 + ((
     -26 + 22*lmst22)*mst22 - 4*mt2)*power4(mg2) + 18*(10 - 3*lmst22)*mt2*
     power4(mst22) + mg2*((1580 - 371*lmst22)*mst26*mt2 + (136 - 78*lmst22)*
     power4(mst22)) + 4*power5(mg2) + 18*(-7 + 2*lmst22)*power5(mst22)) + 6*
     lmg2*(-821*mg6*mst24 + 386*lmt2*mg6*mst24 - 351*mg4*mst26 + 392*lmt2*mg4*
     mst26 - 379*mg6*mst22*mt2 + 150*lmt2*mg6*mst22*mt2 - 788*mg4*mst24*mt2 +
     469*lmt2*mg4*mst24*mt2 + ((-5 + 22*lmst22 - 7*lmt2)*mst22 - 4*mt2)*power4(
     mg2) + 9*lmt2*mt2*power4(mst22) + mg2*((-189 - 175*lmst22 + 220*lmt2)*
     mst26*mt2 + (117 + 31*lmst22 - 58*lmt2)*power4(mst22)) + 4*power5(mg2) - 9
     *lmt2*power5(mst22) + lmst22*(-131*mg6*mst24 - 283*mg4*mst26 - 55*mg6*
     mst22*mt2 - 273*mg4*mst24*mt2 - 9*mt2*power4(mst22) + 9*power5(mst22)))))/
     (9.*mt2)) + power3(DeltaInv(mg2,mt2,mst12))*((2*((6*mst16*(882 + 48*lmst12
      - 2*lmg2*(129 + 28*lmst12 - 71*lmt2) - 546*lmt2 + 40*lmst12*lmt2 + 43*
     power2(lmg2) - 8*power2(lmst12) + 91*power2(lmt2) + 21*power2(Pi)) + mst14
     *mt2*(12516 - 522*lmst12 - 6*lmg2*(525 + 94*lmst12 - 269*lmt2) - 7056*lmt2
      + 738*lmst12*lmt2 + 525*power2(lmg2) + 87*power2(lmst12) + 1176*power2(
     lmt2) + 298*power2(Pi)))*power4(mg2) + mg4*(-(mst12*(3570 - 54*lmg2*(3 + 4
     *lmst12 - 5*lmt2) - 2178*lmt2 + 24*lmst12*(-30 + 19*lmt2) + 27*power2(lmg2
     ) + 120*power2(lmst12) + 363*power2(lmt2) + 85*power2(Pi))) + mt2*(8484 -
     54*lmg2*(3 + 6*lmst12 - 7*lmt2) - 4608*lmt2 + 6*lmst12*(-417 + 193*lmt2) +
     27*power2(lmg2) + 417*power2(lmst12) + 768*power2(lmt2) + 202*power2(Pi)))
     *power4(mst12) + mg6*(mst16*mt2*(21588 - 6*lmg2*(417 + 200*lmst12 - 339*
     lmt2) - 12852*lmt2 + 450*lmst12*(-7 + 5*lmt2) + 417*power2(lmg2) + 525*
     power2(lmst12) + 2142*power2(lmt2) + 514*power2(Pi)) + 2*(1470 - 846*lmt2
     - 24*lmg2*(-15 + 3*lmst12 + 2*lmt2) + 6*lmst12*(-129 + 55*lmt2) - 60*
     power2(lmg2) + 129*power2(lmst12) + 141*power2(lmt2) + 35*power2(Pi))*
     power4(mst12)) + (-(mst14*(3234 - 48*lmg2*(6 + 5*lmst12 - 7*lmt2) - 2106*
     lmt2 + 6*lmst12*(-63 + 61*lmt2) + 48*power2(lmg2) + 63*power2(lmst12) +
     351*power2(lmt2) + 77*power2(Pi))) + 3*mst12*mt2*(-2*lmg2*(87 + 4*lmst12 -
     33*lmt2) + (-222 + 8*lmst12)*lmt2 + 29*power2(lmg2) + 37*power2(lmt2) + 11
     *(42 + power2(Pi))))*power5(mg2) - 9*mg2*(mst12 - mt2)*(42 + 6*lmst12*(-3
     + lmt2) - 18*lmt2 + 3*power2(lmst12) + 3*power2(lmt2) + power2(Pi))*power5
     (mst12) - mst12*(-6*lmg2*(63 + 4*lmst12 - 25*lmt2) + 6*(-87 + 4*lmst12)*
     lmt2 + 63*power2(lmg2) + 87*power2(lmt2) + 25*(42 + power2(Pi)))*power6(
     mg2)))/(9.*mt2) - (4*invdmst*mg*xt*((6*mst16*(882 + 168*lmst12 - 2*lmg2*(
     213 + 20*lmst12 - 91*lmt2) - 498*lmt2 - 16*lmst12*lmt2 + 71*power2(lmg2) -
     28*power2(lmst12) + 83*power2(lmt2) + 21*power2(Pi)) + mst14*mt2*(17850 -
     594*lmst12 - 6*lmg2*(807 + 123*lmst12 - 392*lmt2) - 9864*lmt2 + 936*lmst12
     *lmt2 + 807*power2(lmg2) + 99*power2(lmst12) + 1644*power2(lmt2) + 425*
     power2(Pi)))*power4(mg2) + mg4*(-(mst12*(5754 - 288*lmst12 - 6*lmg2*(135 +
     76*lmst12 - 121*lmt2) - 3834*lmt2 + 552*lmst12*lmt2 + 135*power2(lmg2) +
     48*power2(lmst12) + 639*power2(lmt2) + 137*power2(Pi))) + mt2*(24990 -
     6102*lmst12 - 6*lmg2*(189 + 193*lmst12 - 256*lmt2) - 14184*lmt2 + 3192*
     lmst12*lmt2 + 189*power2(lmg2) + 1017*power2(lmst12) + 2364*power2(lmt2) +
     595*power2(Pi)))*power4(mst12) + mg6*(mst16*mt2*(41286 - 18*lmg2*(339 +
     125*lmst12 - 238*lmt2) - 24444*lmt2 + 6*lmst12*(-807 + 644*lmt2) + 1017*
     power2(lmg2) + 807*power2(lmst12) + 4074*power2(lmt2) + 983*power2(Pi)) +
     2*(4956 - 6*lmg2*(-24 + 55*lmst12 - 47*lmt2) - 3114*lmt2 + 18*lmst12*(-71
     + 42*lmt2) - 24*power2(lmg2) + 213*power2(lmst12) + 519*power2(lmt2) + 118
     *power2(Pi))*power4(mst12)) + (-(mst14*(5964 - 6*lmg2*(168 + 61*lmst12 -
     117*lmt2) - 3654*lmt2 + 6*lmst12*(-75 + 86*lmt2) + 168*power2(lmg2) + 75*
     power2(lmst12) + 609*power2(lmt2) + 142*power2(Pi))) + mst12*mt2*(-6*lmg2*
     (99 + 4*lmst12 - 37*lmt2) + 6*(-123 + 4*lmst12)*lmt2 + 99*power2(lmg2) +
     123*power2(lmt2) + 37*(42 + power2(Pi))))*power5(mg2) - 9*mg2*(3*mst12*(-2
     *lmst12*(15 + lmg2 - 6*lmt2) + 2*(-21 + lmg2)*lmt2 + 5*power2(lmst12) + 7*
     power2(lmt2) + 2*(42 + power2(Pi))) - mt2*(-6*lmst12*(21 + lmg2 - 8*lmt2)
     + 6*(-27 + lmg2)*lmt2 + 21*power2(lmst12) + 27*power2(lmt2) + 8*(42 +
     power2(Pi))))*power5(mst12) - mst12*(-6*lmg2*(75 + 4*lmst12 - 29*lmt2) + 6
     *(-99 + 4*lmst12)*lmt2 + 75*power2(lmg2) + 99*power2(lmt2) + 29*(42 +
     power2(Pi)))*power6(mg2)))/(9.*mt2)) + power3(DeltaInv(mg2,mst22,mt2))*((2
     *((6*mst26*(882 + 48*lmst22 - 2*lmg2*(129 + 28*lmst22 - 71*lmt2) - 546*
     lmt2 + 40*lmst22*lmt2 + 43*power2(lmg2) - 8*power2(lmst22) + 91*power2(
     lmt2) + 21*power2(Pi)) + mst24*mt2*(12516 - 522*lmst22 - 6*lmg2*(525 + 94*
     lmst22 - 269*lmt2) - 7056*lmt2 + 738*lmst22*lmt2 + 525*power2(lmg2) + 87*
     power2(lmst22) + 1176*power2(lmt2) + 298*power2(Pi)))*power4(mg2) + mg4*(-
     (mst22*(3570 - 54*lmg2*(3 + 4*lmst22 - 5*lmt2) - 2178*lmt2 + 24*lmst22*(-
     30 + 19*lmt2) + 27*power2(lmg2) + 120*power2(lmst22) + 363*power2(lmt2) +
     85*power2(Pi))) + mt2*(8484 - 54*lmg2*(3 + 6*lmst22 - 7*lmt2) - 4608*lmt2
     + 6*lmst22*(-417 + 193*lmt2) + 27*power2(lmg2) + 417*power2(lmst22) + 768*
     power2(lmt2) + 202*power2(Pi)))*power4(mst22) + mg6*(mst26*mt2*(21588 - 6*
     lmg2*(417 + 200*lmst22 - 339*lmt2) - 12852*lmt2 + 450*lmst22*(-7 + 5*lmt2)
     + 417*power2(lmg2) + 525*power2(lmst22) + 2142*power2(lmt2) + 514*power2(
     Pi)) + 2*(1470 - 846*lmt2 - 24*lmg2*(-15 + 3*lmst22 + 2*lmt2) + 6*lmst22*(
     -129 + 55*lmt2) - 60*power2(lmg2) + 129*power2(lmst22) + 141*power2(lmt2)
     + 35*power2(Pi))*power4(mst22)) + (-(mst24*(3234 - 48*lmg2*(6 + 5*lmst22 -
     7*lmt2) - 2106*lmt2 + 6*lmst22*(-63 + 61*lmt2) + 48*power2(lmg2) + 63*
     power2(lmst22) + 351*power2(lmt2) + 77*power2(Pi))) + 3*mst22*mt2*(-2*lmg2
     *(87 + 4*lmst22 - 33*lmt2) + (-222 + 8*lmst22)*lmt2 + 29*power2(lmg2) + 37
     *power2(lmt2) + 11*(42 + power2(Pi))))*power5(mg2) - 9*mg2*(mst22 - mt2)*(
     42 + 6*lmst22*(-3 + lmt2) - 18*lmt2 + 3*power2(lmst22) + 3*power2(lmt2) +
     power2(Pi))*power5(mst22) - mst22*(-6*lmg2*(63 + 4*lmst22 - 25*lmt2) + 6*(
     -87 + 4*lmst22)*lmt2 + 63*power2(lmg2) + 87*power2(lmt2) + 25*(42 + power2
     (Pi)))*power6(mg2)))/(9.*mt2) + (4*invdmst*mg*xt*((6*mst26*(882 + 168*
     lmst22 - 2*lmg2*(213 + 20*lmst22 - 91*lmt2) - 498*lmt2 - 16*lmst22*lmt2 +
     71*power2(lmg2) - 28*power2(lmst22) + 83*power2(lmt2) + 21*power2(Pi)) +
     mst24*mt2*(17850 - 594*lmst22 - 6*lmg2*(807 + 123*lmst22 - 392*lmt2) -
     9864*lmt2 + 936*lmst22*lmt2 + 807*power2(lmg2) + 99*power2(lmst22) + 1644*
     power2(lmt2) + 425*power2(Pi)))*power4(mg2) + mg4*(-(mst22*(5754 - 288*
     lmst22 - 6*lmg2*(135 + 76*lmst22 - 121*lmt2) - 3834*lmt2 + 552*lmst22*lmt2
      + 135*power2(lmg2) + 48*power2(lmst22) + 639*power2(lmt2) + 137*power2(Pi
     ))) + mt2*(24990 - 6102*lmst22 - 6*lmg2*(189 + 193*lmst22 - 256*lmt2) -
     14184*lmt2 + 3192*lmst22*lmt2 + 189*power2(lmg2) + 1017*power2(lmst22) +
     2364*power2(lmt2) + 595*power2(Pi)))*power4(mst22) + mg6*(mst26*mt2*(41286
      - 18*lmg2*(339 + 125*lmst22 - 238*lmt2) - 24444*lmt2 + 6*lmst22*(-807 +
     644*lmt2) + 1017*power2(lmg2) + 807*power2(lmst22) + 4074*power2(lmt2) +
     983*power2(Pi)) + 2*(4956 - 6*lmg2*(-24 + 55*lmst22 - 47*lmt2) - 3114*lmt2
      + 18*lmst22*(-71 + 42*lmt2) - 24*power2(lmg2) + 213*power2(lmst22) + 519*
     power2(lmt2) + 118*power2(Pi))*power4(mst22)) + (-(mst24*(5964 - 6*lmg2*(
     168 + 61*lmst22 - 117*lmt2) - 3654*lmt2 + 6*lmst22*(-75 + 86*lmt2) + 168*
     power2(lmg2) + 75*power2(lmst22) + 609*power2(lmt2) + 142*power2(Pi))) +
     mst22*mt2*(-6*lmg2*(99 + 4*lmst22 - 37*lmt2) + 6*(-123 + 4*lmst22)*lmt2 +
     99*power2(lmg2) + 123*power2(lmt2) + 37*(42 + power2(Pi))))*power5(mg2) -
     9*mg2*(3*mst22*(-2*lmst22*(15 + lmg2 - 6*lmt2) + 2*(-21 + lmg2)*lmt2 + 5*
     power2(lmst22) + 7*power2(lmt2) + 2*(42 + power2(Pi))) - mt2*(-6*lmst22*(
     21 + lmg2 - 8*lmt2) + 6*(-27 + lmg2)*lmt2 + 21*power2(lmst22) + 27*power2(
     lmt2) + 8*(42 + power2(Pi))))*power5(mst22) - mst22*(-6*lmg2*(75 + 4*
     lmst22 - 29*lmt2) + 6*(-99 + 4*lmst22)*lmt2 + 75*power2(lmg2) + 99*power2(
     lmt2) + 29*(42 + power2(Pi)))*power6(mg2)))/(9.*mt2)) + (10332 + 1284*
     lmsb22 + 5136*lmsd12 + 72*lmsb22*lmsd12 + 5136*lmsd22 + 72*lmsb22*lmsd22 +
     288*lmsd12*lmsd22 + 1284*lmst12 + 18*lmsb22*lmst12 + 72*lmsd12*lmst12 + 72
     *lmsd22*lmst12 + 1284*lmst22 + 18*lmsb22*lmst22 + 72*lmsd12*lmst22 + 72*
     lmsd22*lmst22 + 18*lmst12*lmst22 + 1872*lmt2 + 72*lmsb22*lmt2 + 288*lmsd12
     *lmt2 + 288*lmsd22*lmt2 + 72*lmst12*lmt2 + 72*lmst22*lmt2 + (384*mg2)/
     msb12 + (384*mg2)/msb22 - 816*invdmsb2g*lmsb22*msb22 + (1536*mg2)/msd12 -
     2416*invdmsd1g*msd12 - 2592*invdmsd1g*lmsd12*msd12 - (2416*msd12)/mg2 + (
     672*lmsd12*msd12)/mg2 + (1536*mg2)/msd22 - 2416*invdmsd2g*msd22 - 2592*
     invdmsd2g*lmsd22*msd22 - (2416*msd22)/mg2 + (672*lmsd22*msd22)/mg2 + (384*
     mg2)/mst12 - (384*lmt2*mg2)/mst12 + (384*mg2)/mst22 - (384*lmt2*mg2)/mst22
      - (2496*mg2)/mt2 + (864*lmt2*mg2)/mt2 - (384*lmt2*mg4)/(mst12*mt2) - (864
     *mst12)/mt2 - (432*lmst12*mst12)/mt2 + (432*lmt2*mst12)/mt2 - (384*lmt2*
     mg4)/(mst22*mt2) - (864*mst22)/mt2 - (432*lmst22*mst22)/mt2 + (432*lmt2*
     mst22)/mt2 + (384*mt2)/mst12 - (384*lmt2*mt2)/mst12 + (384*mt2)/mst22 - (
     384*lmt2*mt2)/mst22 - 2416*msd14*power2(invdmsd1g) + 672*lmsd12*msd14*
     power2(invdmsd1g) - (1936*msd16*power2(invdmsd1g))/mg2 + (2976*lmsd12*
     msd16*power2(invdmsd1g))/mg2 - 2416*msd24*power2(invdmsd2g) + 672*lmsd22*
     msd24*power2(invdmsd2g) - (1936*msd26*power2(invdmsd2g))/mg2 + (2976*
     lmsd22*msd26*power2(invdmsd2g))/mg2 + 1296*power2(lmg2) + 9*power2(lmsb12)
     + 9*power2(lmsb22) + 144*power2(lmsd12) + 144*power2(lmsd22) + 9*power2(
     lmst12) + 9*power2(lmst22) + 144*power2(lmt2) + 768*power2(snb) - 384*
     lmsb22*power2(snb) - (384*msb12*power2(snb))/msb22 - (384*msb22*power2(snb
     ))/msb12 + (384*lmsb22*msb22*power2(snb))/msb12 - (384*lmsb12*msb12*(-1 +
     power2(snb))*power2(snb))/msb22 + 768*power2(snt) - 384*lmst12*power2(snt)
     - 384*lmst22*power2(snt) - (384*mst12*power2(snt))/mst22 + (384*lmst12*
     mst12*power2(snt))/mst22 - (384*mst22*power2(snt))/mst12 + (384*lmst22*
     mst22*power2(snt))/mst12 - 4352*msd16*power3(invdmsd1g) + 3648*lmsd12*
     msd16*power3(invdmsd1g) - 4352*msd26*power3(invdmsd2g) + 3648*lmsd22*msd26
     *power3(invdmsd2g) - 4352*power4(invdmsd1g)*power4(msd12) + 3648*lmsd12*
     power4(invdmsd1g)*power4(msd12) - 4352*power4(invdmsd2g)*power4(msd22) +
     3648*lmsd22*power4(invdmsd2g)*power4(msd22) - 768*power4(snb) + 384*lmsb22
     *power4(snb) + (384*msb12*power4(snb))/msb22 + (384*msb22*power4(snb))/
     msb12 - (384*lmsb22*msb22*power4(snb))/msb12 + 6*lmsb12*(214 + 3*lmsb22 +
     12*lmsd12 + 12*lmsd22 + 3*lmst12 + 3*lmst22 + 12*lmt2 - 136*invdmsb1g*
     msb12 - 64*power2(snb) + 64*power4(snb)) - 768*power4(snt) + 384*lmst12*
     power4(snt) + 384*lmst22*power4(snt) + (384*mst12*power4(snt))/mst22 - (
     384*lmst12*mst12*power4(snt))/mst22 + (384*mst22*power4(snt))/mst12 - (384
     *lmst22*mst22*power4(snt))/mst12 - (5664*power4(invdmsd1g)*power5(msd12))/
     mg2 - (448*lmsd12*power4(invdmsd1g)*power5(msd12))/mg2 - 10016*power5(
     invdmsd1g)*power5(msd12) + 3200*lmsd12*power5(invdmsd1g)*power5(msd12) - (
     5664*power4(invdmsd2g)*power5(msd22))/mg2 - (448*lmsd22*power4(invdmsd2g)*
     power5(msd22))/mg2 - 10016*power5(invdmsd2g)*power5(msd22) + 3200*lmsd22*
     power5(invdmsd2g)*power5(msd22) - 10016*power6(invdmsd1g)*power6(msd12) +
     3200*lmsd12*power6(invdmsd1g)*power6(msd12) - 10016*power6(invdmsd2g)*
     power6(msd22) + 3200*lmsd22*power6(invdmsd2g)*power6(msd22) + (10016*
     power6(invdmsd1g)*power7(msd12))/mg2 - (3200*lmsd12*power6(invdmsd1g)*
     power7(msd12))/mg2 + (10016*power6(invdmsd2g)*power7(msd22))/mg2 - (3200*
     lmsd22*power6(invdmsd2g)*power7(msd22))/mg2 + (8*lmg2*((48*mg4*(mst12 +
     mst22))/(mst12*mst22) + mt2*(-534 + 27*lmsb12 + 27*lmsb22 + 108*lmsd12 +
     108*lmsd22 + 27*lmst12 + 27*lmst22 + 108*lmt2 + 102*invdmsb1g*msb12 + 102*
     invdmsb2g*msb22 + 492*invdmsd1g*msd12 + 492*invdmsd2g*msd22 + 84*msd14*
     power2(invdmsd1g) + 84*msd24*power2(invdmsd2g) + 456*msd16*power3(
     invdmsd1g) + 456*msd26*power3(invdmsd2g) + 456*power4(invdmsd1g)*power4(
     msd12) + 456*power4(invdmsd2g)*power4(msd22) + 400*power5(invdmsd1g)*
     power5(msd12) + 400*power5(invdmsd2g)*power5(msd22) + 400*power6(invdmsd1g
     )*power6(msd12) + 400*power6(invdmsd2g)*power6(msd22)) + (4*(3*mg4*(-9 - (
     4*mt2)/msb12 - (4*mt2)/msb22 - (16*mt2)/msd12 - (16*mt2)/msd22) + mt2*(21*
     msd12 + 21*msd22 + 93*msd16*power2(invdmsd1g) + 93*msd26*power2(invdmsd2g)
     - 14*power4(invdmsd1g)*power5(msd12) - 14*power4(invdmsd2g)*power5(msd22)
     - 100*power6(invdmsd1g)*power7(msd12) - 100*power6(invdmsd2g)*power7(msd22
     ))))/mg2))/mt2)/864.;

   return power4(g3) * result * twoLoop;
}

/**
 * @brief 2-loop \f$O(\alpha_t\alpha_s)\f$ contributions to \f$\Delta g_3^{(2L)}\f$ [arXiv:1009.5455]
 *
 * The function returns \f$\Delta g_3^{(2L)}\f$ which appears in
 * conversion of the strong \f$\overline{\text{MS}}\f$ gauge coupling
 * of the Standard Model with 5 active quark flavours,
 * \f$g_3^{\text{SM(5)},\overline{\text{MS}}}\f$, to the strong
 * \f$\overline{\text{DR}}\f$ gauge coupling of the full MSSM
 * \f$g_3^{\text{MSSM},\overline{\text{DR}}}\f$,
 * \f{align*}{
  g_3^{\text{SM(5)},\overline{\text{MS}}} =
  g_3^{\text{MSSM},\overline{\text{DR}}} \left[
     1 + \Delta g_3^{(1L)} + \Delta g_3^{(2L)}
  \right]
 * \f}
 */
Real delta_alpha_s_2loop_at_as(const Parameters& pars)
{
   const Real g3        = pars.g3;
   const Real yt        = pars.yt;
   const Real xt        = pars.xt;
   const Real xb        = pars.xb;
   const Real mt        = pars.mt;
   const Real mt2       = power2(pars.mt);
   const Real mb        = pars.mb;
   const Real mst12     = power2(pars.mst1);
   const Real mst14     = power2(mst12);
   const Real mst16     = mst12*mst14;
   const Real mst22     = power2(pars.mst2);
   const Real mst24     = power2(mst22);
   const Real mst26     = mst22*mst24;
   const Real msb12     = power2(pars.msb1);
   const Real msb14     = power2(msb12);
   const Real msb16     = msb12*msb14;
   const Real msb22     = power2(pars.msb2);
   const Real msb24     = power2(msb22);
   const Real msb26     = msb22*msb24;
   const Real mw2       = power2(pars.mw);
   const Real mz2       = power2(pars.mz);
   const Real mh2       = power2(pars.mh);
   const Real mH2       = power2(pars.mH);
   const Real mC2       = power2(pars.mC);
   const Real mA2       = power2(pars.mA);
   const Real mu        = pars.mu;
   const Real mu2       = power2(pars.mu);
   const Real tb        = pars.tb;
   const Real sb        = tb / std::sqrt(1. + power2(tb));
   const Real cb        = 1. / std::sqrt(1. + power2(tb));
   const Real Q2        = power2(pars.Q);
   const Real snt       = calc_sin_theta(mt, xt, mst12, mst22);
   const Real snb       = calc_sin_theta(mb, xb, msb12, msb22);
   const Real alpha     = calc_alpha(mh2, mH2, tb);
   const Real sa        = std::sin(alpha);
   const Real ca        = std::cos(alpha);
   const Real At        = xt + mu*cb/sb;
   const Real invdmst   = 1/(mst12 - mst22);
   const Real invdct    = 1/(mC2 - mt2);
   const Real invdtw    = 1/(mt2 - mw2);
   const Real lmst12    = std::log(mst12/Q2);
   const Real lmst22    = std::log(mst22/Q2);
   const Real lmsb12    = std::log(msb12/Q2);
   const Real lmsb22    = std::log(msb22/Q2);
   const Real lmt2      = std::log(mt2/Q2);
   const Real lmw2      = std::log(mw2/Q2);
   const Real lmz2      = std::log(mz2/Q2);
   const Real lmu2      = std::log(mu2/Q2);
   const Real lmH2      = std::log(mH2/Q2);
   const Real lmA2      = std::log(mA2/Q2);
   const Real lmh2      = std::log(mh2/Q2);
   const Real lmC2      = std::log(mC2/Q2);

   const Real result =
   -(mt2*DeltaInv(mt2,mt2,mA2)*power2(cb)*(mA2*(38 - 26*lmt2 + 2*lmA2*(-5 + 6*
     lmt2) + 6*power2(lmA2) - 6*power2(lmt2) + power2(Pi)) + mt2*(294 - 96*lmt2
      + 12*lmA2*(-13 + 5*lmt2) + 30*power2(lmA2) - 6*power2(lmt2) + 7*power2(Pi
     ))))/12. - (mt2*DeltaInv(mt2,mt2,mz2)*(-1 + power2(cb))*(-(mz2*(38 - 10*
     lmz2 + 2*lmt2*(-13 + 6*lmz2) - 6*power2(lmt2) + 6*power2(lmz2) + power2(Pi
     ))) + mt2*(lmt2*(96 - 60*lmz2) + 156*lmz2 + 6*power2(lmt2) - 30*power2(
     lmz2) - 7*(42 + power2(Pi)))))/12. + (mu2*DeltaInv(mt2,mu2,mst12)*(mst12*(
     mst14 + mu2*(-3*mt2 + 2*lmt2*mt2 - 2*lmu2*mt2 + mu2 + 4*lmt2*mu2 - 4*lmu2*
     mu2)) - (lmt2 - lmu2)*(mt2 - mu2)*power2(mu2) + power2(mst12)*((-1 -
     lmst12 + lmt2)*mt2 + mu2*(40 - 17*lmt2 + 6*lmst12*(-2 + lmt2 - lmu2) - 7*
     lmu2 + 6*lmt2*lmu2 + 6*power2(lmt2) + power2(Pi)))))/(6.*mst12*mt2) + (mu2
     *DeltaInv(mst22,mt2,mu2)*(mst22*(mst24 + mu2*(-3*mt2 + 2*lmt2*mt2 - 2*lmu2
     *mt2 + mu2 + 4*lmt2*mu2 - 4*lmu2*mu2)) - (lmt2 - lmu2)*(mt2 - mu2)*power2(
     mu2) + power2(mst22)*((-1 - lmst22 + lmt2)*mt2 + mu2*(40 - 17*lmt2 + 6*
     lmst22*(-2 + lmt2 - lmu2) - 7*lmu2 + 6*lmt2*lmu2 + 6*power2(lmt2) + power2
     (Pi)))))/(6.*mst22*mt2) - (mt2*DeltaInv(mh2,mt2,mt2)*(mh2*(132 - 58*lmt2 +
     lmh2*(-58 + 36*lmt2) + 18*power2(lmh2) - 18*power2(lmt2) + 3*power2(Pi)) +
     mt2*(606 - 200*lmt2 + 4*lmh2*(-77 + 27*lmt2) + 54*power2(lmh2) + 18*power2
     (lmt2) + 15*power2(Pi)))*(-1 + power2(sa)))/36. + (mt2*DeltaInv(mH2,mt2,
     mt2)*(mH2*(132 - 58*lmt2 + lmH2*(-58 + 36*lmt2) + 18*power2(lmH2) - 18*
     power2(lmt2) + 3*power2(Pi)) + mt2*(606 - 200*lmt2 + 4*lmH2*(-77 + 27*lmt2
     ) + 54*power2(lmH2) + 18*power2(lmt2) + 15*power2(Pi)))*power2(sa))/36. +
     (mu2*DeltaInv(msb12,mt2,mu2)*(-(msb12*(msb14 + mu2*(-3*mt2 + 2*lmt2*mt2 -
     2*lmu2*mt2 + mu2 + 4*lmt2*mu2 - 4*lmu2*mu2))) + (lmt2 - lmu2)*(mt2 - mu2)*
     power2(mu2) + power2(msb12)*((1 + lmsb12 - lmt2)*mt2 - mu2*(40 - 17*lmt2 +
     6*lmsb12*(-2 + lmt2 - lmu2) - 7*lmu2 + 6*lmt2*lmu2 + 6*power2(lmt2) +
     power2(Pi))))*(-1 + power2(snb)))/(6.*msb12*mt2) + (mu2*DeltaInv(mt2,mu2,
     msb22)*(msb22*(msb24 + mu2*(-3*mt2 + 2*lmt2*mt2 - 2*lmu2*mt2 + mu2 + 4*
     lmt2*mu2 - 4*lmu2*mu2)) - (lmt2 - lmu2)*(mt2 - mu2)*power2(mu2) + power2(
     msb22)*((-1 - lmsb22 + lmt2)*mt2 + mu2*(40 - 17*lmt2 + 6*lmsb22*(-2 + lmt2
      - lmu2) - 7*lmu2 + 6*lmt2*lmu2 + 6*power2(lmt2) + power2(Pi))))*power2(
     snb))/(6.*msb22*mt2) + DeltaInv(mH2,mst12,mst12)*((invdmst*((1 + lmH2 -
     lmst12)*mH2 + 6*(-lmH2 + lmst12)*mst12)*mt2*sa*(At*sa*(-(cb*mu) + At*sb) +
     ca*(cb*mu2 - At*mu*sb)))/(3.*sb) + (((1 + lmH2 - lmst12)*mH2 + 6*(-lmH2 +
     lmst12)*mst12)*(mt2*power2(sa) + (At*sa*(2*ca*mu - At*sa) + mu2*(-1 +
     power2(sa)))*(-1 + power2(snt))*power2(snt)))/6.) + DeltaInv(mH2,mst22,
     mst22)*((invdmst*((1 + lmH2 - lmst22)*mH2 + 6*(-lmH2 + lmst22)*mst22)*mt2*
     sa*(At*sa*(cb*mu - At*sb) + ca*(-(cb*mu2) + At*mu*sb)))/(3.*sb) + (((1 +
     lmH2 - lmst22)*mH2 + 6*(-lmH2 + lmst22)*mst22)*(mt2*power2(sa) + (At*sa*(2
     *ca*mu - At*sa) + mu2*(-1 + power2(sa)))*(-1 + power2(snt))*power2(snt)))/
     6.) + DeltaInv(mh2,mst12,mst12)*(-(invdmst*((1 + lmh2 - lmst12)*mh2 + 6*(-
     lmh2 + lmst12)*mst12)*mt2*(ca*cb*mu2*sa - At*mu*(ca*sa*sb + cb*(-1 +
     power2(sa))) + sb*power2(At)*(-1 + power2(sa))))/(3.*sb) - (((1 + lmh2 -
     lmst12)*mh2 + 6*(-lmh2 + lmst12)*mst12)*(mt2*(-1 + power2(sa)) - (-2*At*ca
     *mu*sa + power2(At)*(-1 + power2(sa)) - mu2*power2(sa))*(-1 + power2(snt))
     *power2(snt)))/6.) + DeltaInv(mst22,mst22,mh2)*((invdmst*((1 + lmh2 -
     lmst22)*mh2 + 6*(-lmh2 + lmst22)*mst22)*mt2*(ca*cb*mu2*sa - At*mu*(ca*sa*
     sb + cb*(-1 + power2(sa))) + sb*power2(At)*(-1 + power2(sa))))/(3.*sb) - (
     ((1 + lmh2 - lmst22)*mh2 + 6*(-lmh2 + lmst22)*mst22)*(mt2*(-1 + power2(sa)
     ) - (-2*At*ca*mu*sa + power2(At)*(-1 + power2(sa)) - mu2*power2(sa))*(-1 +
     power2(snt))*power2(snt)))/6.) + ((invdmst*mt2*(6*mst16*(42 + 8*lmh2*(-3 +
     lmst12) - 12*lmst12 + 4*power2(lmh2) + power2(Pi)) + mh2*mst14*(42 + 12*
     lmh2*(-1 + lmst12) - 24*lmst12 + 6*power2(lmh2) - 6*power2(lmst12) +
     power2(Pi)))*(ca*cb*mu2*sa - At*mu*(ca*sa*sb + cb*(-1 + power2(sa))) + sb*
     power2(At)*(-1 + power2(sa))))/(3.*sb) + ((6*mst16*(42 + 8*lmh2*(-3 +
     lmst12) - 12*lmst12 + 4*power2(lmh2) + power2(Pi)) + mh2*mst14*(42 + 12*
     lmh2*(-1 + lmst12) - 24*lmst12 + 6*power2(lmh2) - 6*power2(lmst12) +
     power2(Pi)))*(mt2*(-1 + power2(sa)) - (-2*At*ca*mu*sa + power2(At)*(-1 +
     power2(sa)) - mu2*power2(sa))*(-1 + power2(snt))*power2(snt)))/6.)*power2(
     DeltaInv(mh2,mst12,mst12)) + Fin3(mh2,mt2,mt2,Q2)*(DeltaInv(mh2,mt2,mt2)*(
     mt2 - mt2*power2(sa)) - 3*mh2*power2(mt2)*(-1 + power2(sa))*power2(
     DeltaInv(mh2,mt2,mt2))) + ((invdmst*mt2*sa*(At*sa*(cb*mu - At*sb) + ca*(-(
     cb*mu2) + At*mu*sb))*(6*mst16*(42 + 8*lmH2*(-3 + lmst12) - 12*lmst12 + 4*
     power2(lmH2) + power2(Pi)) + mH2*mst14*(42 + 12*lmH2*(-1 + lmst12) - 24*
     lmst12 + 6*power2(lmH2) - 6*power2(lmst12) + power2(Pi))))/(3.*sb) - ((6*
     mst16*(42 + 8*lmH2*(-3 + lmst12) - 12*lmst12 + 4*power2(lmH2) + power2(Pi)
     ) + mH2*mst14*(42 + 12*lmH2*(-1 + lmst12) - 24*lmst12 + 6*power2(lmH2) - 6
     *power2(lmst12) + power2(Pi)))*(mt2*power2(sa) + (At*sa*(2*ca*mu - At*sa)
     + mu2*(-1 + power2(sa)))*(-1 + power2(snt))*power2(snt)))/6.)*power2(
     DeltaInv(mH2,mst12,mst12)) + ((invdmst*mt2*sa*(At*sa*(-(cb*mu) + At*sb) +
     ca*(cb*mu2 - At*mu*sb))*(6*mst26*(42 + 8*lmH2*(-3 + lmst22) - 12*lmst22 +
     4*power2(lmH2) + power2(Pi)) + mH2*mst24*(42 + 12*lmH2*(-1 + lmst22) - 24*
     lmst22 + 6*power2(lmH2) - 6*power2(lmst22) + power2(Pi))))/(3.*sb) - ((6*
     mst26*(42 + 8*lmH2*(-3 + lmst22) - 12*lmst22 + 4*power2(lmH2) + power2(Pi)
     ) + mH2*mst24*(42 + 12*lmH2*(-1 + lmst22) - 24*lmst22 + 6*power2(lmH2) - 6
     *power2(lmst22) + power2(Pi)))*(mt2*power2(sa) + (At*sa*(2*ca*mu - At*sa)
     + mu2*(-1 + power2(sa)))*(-1 + power2(snt))*power2(snt)))/6.)*power2(
     DeltaInv(mH2,mst22,mst22)) + Fin3(mH2,mt2,mt2,Q2)*(mt2*DeltaInv(mH2,mt2,
     mt2)*power2(sa) + 3*mH2*power2(mt2)*power2(sa)*power2(DeltaInv(mH2,mt2,mt2
     ))) - (power2(mu2)*(msb16*(42 + 6*lmu2 - 6*lmt2*(3 + lmu2) + 6*lmsb12*(-4
     + lmt2 + lmu2) + 6*power2(lmsb12) + power2(Pi)) + mu2*(6*(lmt2 - lmu2)*(
     mt2 - mu2)*mu2 + 6*msb12*mt2*(42 + 4*lmsb12*(-1 + lmt2 - lmu2) - 12*lmu2 +
     lmt2*(-20 + 6*lmu2) + 5*power2(lmt2) + power2(lmu2) + power2(Pi)) + msb12*
     mu2*(294 + 6*lmsb12*(2 + lmt2 - lmu2) - 114*lmu2 + 6*lmt2*(-25 + 7*lmu2) +
     24*power2(lmt2) + 18*power2(lmu2) + 7*power2(Pi))) + 2*msb14*(mt2*(84 - 39
     *lmt2 + 6*lmsb12*(-5 + 2*lmt2 - lmu2) - 3*lmu2 + 6*lmt2*lmu2 + 3*power2(
     lmsb12) + 9*power2(lmt2) + 2*power2(Pi)) + mu2*(336 - 201*lmt2 + 6*lmsb12*
     (-11 + 7*lmt2 - 4*lmu2) - 21*lmu2 + 30*lmt2*lmu2 + 9*power2(lmsb12) + 36*
     power2(lmt2) + 3*power2(lmu2) + 8*power2(Pi))))*(-1 + power2(snb))*power2(
     DeltaInv(msb12,mt2,mu2)))/(6.*mt2) + (-(invdmst*mt2*(6*mst26*(42 + 8*lmh2*
     (-3 + lmst22) - 12*lmst22 + 4*power2(lmh2) + power2(Pi)) + mh2*mst24*(42 +
     12*lmh2*(-1 + lmst22) - 24*lmst22 + 6*power2(lmh2) - 6*power2(lmst22) +
     power2(Pi)))*(ca*cb*mu2*sa - At*mu*(ca*sa*sb + cb*(-1 + power2(sa))) + sb*
     power2(At)*(-1 + power2(sa))))/(3.*sb) + ((6*mst26*(42 + 8*lmh2*(-3 +
     lmst22) - 12*lmst22 + 4*power2(lmh2) + power2(Pi)) + mh2*mst24*(42 + 12*
     lmh2*(-1 + lmst22) - 24*lmst22 + 6*power2(lmh2) - 6*power2(lmst22) +
     power2(Pi)))*(mt2*(-1 + power2(sa)) - (-2*At*ca*mu*sa + power2(At)*(-1 +
     power2(sa)) - mu2*power2(sa))*(-1 + power2(snt))*power2(snt)))/6.)*power2(
     DeltaInv(mst22,mst22,mh2)) + (power2(mu2)*(mst26*(42 + 6*lmu2 - 6*lmt2*(3
     + lmu2) + 6*lmst22*(-4 + lmt2 + lmu2) + 6*power2(lmst22) + power2(Pi)) +
     mu2*(6*(lmt2 - lmu2)*(mt2 - mu2)*mu2 + 6*mst22*mt2*(42 + 4*lmst22*(-1 +
     lmt2 - lmu2) - 12*lmu2 + lmt2*(-20 + 6*lmu2) + 5*power2(lmt2) + power2(
     lmu2) + power2(Pi)) + mst22*mu2*(294 + 6*lmst22*(2 + lmt2 - lmu2) - 114*
     lmu2 + 6*lmt2*(-25 + 7*lmu2) + 24*power2(lmt2) + 18*power2(lmu2) + 7*
     power2(Pi))) + 2*mst24*(mt2*(84 - 39*lmt2 + 6*lmst22*(-5 + 2*lmt2 - lmu2)
     - 3*lmu2 + 6*lmt2*lmu2 + 3*power2(lmst22) + 9*power2(lmt2) + 2*power2(Pi))
     + mu2*(336 - 201*lmt2 + 6*lmst22*(-11 + 7*lmt2 - 4*lmu2) - 21*lmu2 + 30*
     lmt2*lmu2 + 9*power2(lmst22) + 36*power2(lmt2) + 3*power2(lmu2) + 8*power2
     (Pi))))*power2(DeltaInv(mst22,mt2,mu2)))/(6.*mt2) + ((invdmst*mt2*mw2*(-3*
     lmst22*(-mst26 + mw2*(2*(7 + lmw2)*mst24 + (-1 + 2*lmw2)*mst22*mw2) +
     msb12*(mst24 + (3 + 2*lmw2)*mst22*mw2)) + 3*lmsb12*(-mst26 + mw2*((-23 + 6
     *lmst22 + 2*lmw2)*mst24 + ((-23 + 2*lmst22 + 6*lmw2)*mst22 - mw2)*mw2) +
     msb12*(mst24 + mw2*(2*(-3 + lmst22 + lmw2)*mst22 + mw2))) + 6*mw2*(msb12*
     mst22 + 2*(mst24 + mst22*mw2))*power2(lmsb12) + 6*mst24*mw2*power2(lmst22)
     + mw2*(-3*lmw2*msb12*mw2 + msb12*mst22*(42 - 9*lmw2 + power2(Pi)) + 3*
     mst24*(42 + lmw2 + power2(Pi)) + 3*mw2*(lmw2*mw2 + mst22*(42 - 14*lmw2 + 2
     *power2(lmw2) + power2(Pi)))))*(sb*power2(At)*(-1 + power2(cb)) - mu2*sb*
     power2(cb) + At*cb*mu*(1 - power2(cb) + power2(sb)))*(-1 + power2(snb)))/(
     3.*msb12*sb) + (mw2*(-3*lmst22*(-mst26 + mw2*(2*(7 + lmw2)*mst24 + (-1 + 2
     *lmw2)*mst22*mw2) + msb12*(mst24 + (3 + 2*lmw2)*mst22*mw2)) + 3*lmsb12*(-
     mst26 + mw2*((-23 + 6*lmst22 + 2*lmw2)*mst24 + ((-23 + 2*lmst22 + 6*lmw2)*
     mst22 - mw2)*mw2) + msb12*(mst24 + mw2*(2*(-3 + lmst22 + lmw2)*mst22 + mw2
     ))) + 6*mw2*(msb12*mst22 + 2*(mst24 + mst22*mw2))*power2(lmsb12) + 6*mst24
     *mw2*power2(lmst22) + mw2*(-3*lmw2*msb12*mw2 + msb12*mst22*(42 - 9*lmw2 +
     power2(Pi)) + 3*mst24*(42 + lmw2 + power2(Pi)) + 3*mw2*(lmw2*mw2 + mst22*(
     42 - 14*lmw2 + 2*power2(lmw2) + power2(Pi)))))*(-1 + power2(snb))*(2*At*cb
     *mu*sb*(-1 + power2(snt)) + power2(At)*(-1 + power2(cb))*(-1 + power2(snt)
     ) + mt2*power2(snt) + power2(cb)*(mu2 - mt2*power2(snt) - mu2*power2(snt))
     ))/(6.*msb12))*power2(DeltaInv(mst22,mw2,msb12)) + (-(invdmst*mt2*mw2*(-3*
     lmst22*(-mst26 + mw2*(2*(7 + lmw2)*mst24 + (-1 + 2*lmw2)*mst22*mw2) +
     msb22*(mst24 + (3 + 2*lmw2)*mst22*mw2)) + 3*lmsb22*(-mst26 + mw2*((-23 + 6
     *lmst22 + 2*lmw2)*mst24 + ((-23 + 2*lmst22 + 6*lmw2)*mst22 - mw2)*mw2) +
     msb22*(mst24 + mw2*(2*(-3 + lmst22 + lmw2)*mst22 + mw2))) + 6*mw2*(msb22*
     mst22 + 2*(mst24 + mst22*mw2))*power2(lmsb22) + 6*mst24*mw2*power2(lmst22)
     + mw2*(-3*lmw2*msb22*mw2 + msb22*mst22*(42 - 9*lmw2 + power2(Pi)) + 3*
     mst24*(42 + lmw2 + power2(Pi)) + 3*mw2*(lmw2*mw2 + mst22*(42 - 14*lmw2 + 2
     *power2(lmw2) + power2(Pi)))))*(sb*power2(At)*(-1 + power2(cb)) - mu2*sb*
     power2(cb) + At*cb*mu*(1 - power2(cb) + power2(sb)))*power2(snb))/(3.*
     msb22*sb) - (mw2*(-3*lmst22*(-mst26 + mw2*(2*(7 + lmw2)*mst24 + (-1 + 2*
     lmw2)*mst22*mw2) + msb22*(mst24 + (3 + 2*lmw2)*mst22*mw2)) + 3*lmsb22*(-
     mst26 + mw2*((-23 + 6*lmst22 + 2*lmw2)*mst24 + ((-23 + 2*lmst22 + 6*lmw2)*
     mst22 - mw2)*mw2) + msb22*(mst24 + mw2*(2*(-3 + lmst22 + lmw2)*mst22 + mw2
     ))) + 6*mw2*(msb22*mst22 + 2*(mst24 + mst22*mw2))*power2(lmsb22) + 6*mst24
     *mw2*power2(lmst22) + mw2*(-3*lmw2*msb22*mw2 + msb22*mst22*(42 - 9*lmw2 +
     power2(Pi)) + 3*mst24*(42 + lmw2 + power2(Pi)) + 3*mw2*(lmw2*mw2 + mst22*(
     42 - 14*lmw2 + 2*power2(lmw2) + power2(Pi)))))*power2(snb)*(2*At*cb*mu*sb*
     (-1 + power2(snt)) + power2(At)*(-1 + power2(cb))*(-1 + power2(snt)) + mt2
     *power2(snt) + power2(cb)*(mu2 - mt2*power2(snt) - mu2*power2(snt))))/(6.*
     msb22))*power2(DeltaInv(mst22,mw2,msb22)) + (mz2*(2*At*cb*mu*sb + power2(
     At)*(-1 + power2(cb)) - mu2*power2(cb))*(-3*lmst22*(-mst26 + mz2*(2*(7 +
     lmz2)*mst24 + (-1 + 2*lmz2)*mst22*mz2) + mst12*(mst24 + (3 + 2*lmz2)*mst22
     *mz2)) + 3*lmst12*(-mst26 + mz2*((-23 + 6*lmst22 + 2*lmz2)*mst24 + ((-23 +
     2*lmst22 + 6*lmz2)*mst22 - mz2)*mz2) + mst12*(mst24 + mz2*(2*(-3 + lmst22
     + lmz2)*mst22 + mz2))) + 6*mz2*(mst12*mst22 + 2*(mst24 + mst22*mz2))*
     power2(lmst12) + 6*mst24*mz2*power2(lmst22) + mz2*(-3*lmz2*mst12*mz2 +
     mst12*mst22*(42 - 9*lmz2 + power2(Pi)) + 3*mst24*(42 + lmz2 + power2(Pi))
     + 3*mz2*(lmz2*mz2 + mst22*(42 - 14*lmz2 + 2*power2(lmz2) + power2(Pi)))))*
     power2(DeltaInv(mst22,mz2,mst12)))/(12.*mst12) + (power2(mu2)*(msb26*(42 +
     6*lmu2 - 6*lmt2*(3 + lmu2) + 6*lmsb22*(-4 + lmt2 + lmu2) + 6*power2(lmsb22
     ) + power2(Pi)) + mu2*(6*(lmt2 - lmu2)*(mt2 - mu2)*mu2 + 6*msb22*mt2*(42 +
     4*lmsb22*(-1 + lmt2 - lmu2) - 12*lmu2 + lmt2*(-20 + 6*lmu2) + 5*power2(
     lmt2) + power2(lmu2) + power2(Pi)) + msb22*mu2*(294 + 6*lmsb22*(2 + lmt2 -
     lmu2) - 114*lmu2 + 6*lmt2*(-25 + 7*lmu2) + 24*power2(lmt2) + 18*power2(
     lmu2) + 7*power2(Pi))) + 2*msb24*(mt2*(84 - 39*lmt2 + 6*lmsb22*(-5 + 2*
     lmt2 - lmu2) - 3*lmu2 + 6*lmt2*lmu2 + 3*power2(lmsb22) + 9*power2(lmt2) +
     2*power2(Pi)) + mu2*(336 - 201*lmt2 + 6*lmsb22*(-11 + 7*lmt2 - 4*lmu2) -
     21*lmu2 + 30*lmt2*lmu2 + 9*power2(lmsb22) + 36*power2(lmt2) + 3*power2(
     lmu2) + 8*power2(Pi))))*power2(snb)*power2(DeltaInv(mt2,mu2,msb22)))/(6.*
     mt2) + (power2(mu2)*(mst16*(42 + 6*lmu2 - 6*lmt2*(3 + lmu2) + 6*lmst12*(-4
      + lmt2 + lmu2) + 6*power2(lmst12) + power2(Pi)) + mu2*(6*(lmt2 - lmu2)*(
     mt2 - mu2)*mu2 + 6*mst12*mt2*(42 + 4*lmst12*(-1 + lmt2 - lmu2) - 12*lmu2 +
     lmt2*(-20 + 6*lmu2) + 5*power2(lmt2) + power2(lmu2) + power2(Pi)) + mst12*
     mu2*(294 + 6*lmst12*(2 + lmt2 - lmu2) - 114*lmu2 + 6*lmt2*(-25 + 7*lmu2) +
     24*power2(lmt2) + 18*power2(lmu2) + 7*power2(Pi))) + 2*mst14*(mt2*(84 - 39
     *lmt2 + 6*lmst12*(-5 + 2*lmt2 - lmu2) - 3*lmu2 + 6*lmt2*lmu2 + 3*power2(
     lmst12) + 9*power2(lmt2) + 2*power2(Pi)) + mu2*(336 - 201*lmt2 + 6*lmst12*
     (-11 + 7*lmt2 - 4*lmu2) - 21*lmu2 + 30*lmt2*lmu2 + 9*power2(lmst12) + 36*
     power2(lmt2) + 3*power2(lmu2) + 8*power2(Pi))))*power2(DeltaInv(mt2,mu2,
     mst12)))/(6.*mt2) + (-(invdmst*mt2*mw2*(-3*lmst12*(-mst16 + mw2*(2*(7 +
     lmw2)*mst14 + (-1 + 2*lmw2)*mst12*mw2) + msb12*(mst14 + (3 + 2*lmw2)*mst12
     *mw2)) + 3*lmsb12*(-mst16 + mw2*((-23 + 6*lmst12 + 2*lmw2)*mst14 + ((-23 +
     2*lmst12 + 6*lmw2)*mst12 - mw2)*mw2) + msb12*(mst14 + mw2*(2*(-3 + lmst12
     + lmw2)*mst12 + mw2))) + 6*mw2*(msb12*mst12 + 2*(mst14 + mst12*mw2))*
     power2(lmsb12) + 6*mst14*mw2*power2(lmst12) + mw2*(-3*lmw2*msb12*mw2 +
     msb12*mst12*(42 - 9*lmw2 + power2(Pi)) + 3*mst14*(42 + lmw2 + power2(Pi))
     + 3*mw2*(lmw2*mw2 + mst12*(42 - 14*lmw2 + 2*power2(lmw2) + power2(Pi)))))*
     (sb*power2(At)*(-1 + power2(cb)) - mu2*sb*power2(cb) + At*cb*mu*(1 -
     power2(cb) + power2(sb)))*(-1 + power2(snb)))/(3.*msb12*sb) + (mw2*(-3*
     lmst12*(-mst16 + mw2*(2*(7 + lmw2)*mst14 + (-1 + 2*lmw2)*mst12*mw2) +
     msb12*(mst14 + (3 + 2*lmw2)*mst12*mw2)) + 3*lmsb12*(-mst16 + mw2*((-23 + 6
     *lmst12 + 2*lmw2)*mst14 + ((-23 + 2*lmst12 + 6*lmw2)*mst12 - mw2)*mw2) +
     msb12*(mst14 + mw2*(2*(-3 + lmst12 + lmw2)*mst12 + mw2))) + 6*mw2*(msb12*
     mst12 + 2*(mst14 + mst12*mw2))*power2(lmsb12) + 6*mst14*mw2*power2(lmst12)
     + mw2*(-3*lmw2*msb12*mw2 + msb12*mst12*(42 - 9*lmw2 + power2(Pi)) + 3*
     mst14*(42 + lmw2 + power2(Pi)) + 3*mw2*(lmw2*mw2 + mst12*(42 - 14*lmw2 + 2
     *power2(lmw2) + power2(Pi)))))*(-1 + power2(snb))*(mt2*(-1 + power2(cb))*(
     -1 + power2(snt)) + (-2*At*cb*mu*sb - power2(At)*(-1 + power2(cb)) + mu2*
     power2(cb))*power2(snt)))/(6.*msb12))*power2(DeltaInv(mw2,msb12,mst12)) +
     ((invdmst*mt2*mw2*(-3*lmst12*(-mst16 + mw2*(2*(7 + lmw2)*mst14 + (-1 + 2*
     lmw2)*mst12*mw2) + msb22*(mst14 + (3 + 2*lmw2)*mst12*mw2)) + 3*lmsb22*(-
     mst16 + mw2*((-23 + 6*lmst12 + 2*lmw2)*mst14 + ((-23 + 2*lmst12 + 6*lmw2)*
     mst12 - mw2)*mw2) + msb22*(mst14 + mw2*(2*(-3 + lmst12 + lmw2)*mst12 + mw2
     ))) + 6*mw2*(msb22*mst12 + 2*(mst14 + mst12*mw2))*power2(lmsb22) + 6*mst14
     *mw2*power2(lmst12) + mw2*(-3*lmw2*msb22*mw2 + msb22*mst12*(42 - 9*lmw2 +
     power2(Pi)) + 3*mst14*(42 + lmw2 + power2(Pi)) + 3*mw2*(lmw2*mw2 + mst12*(
     42 - 14*lmw2 + 2*power2(lmw2) + power2(Pi)))))*(sb*power2(At)*(-1 + power2
     (cb)) - mu2*sb*power2(cb) + At*cb*mu*(1 - power2(cb) + power2(sb)))*power2
     (snb))/(3.*msb22*sb) - (mw2*(-3*lmst12*(-mst16 + mw2*(2*(7 + lmw2)*mst14 +
     (-1 + 2*lmw2)*mst12*mw2) + msb22*(mst14 + (3 + 2*lmw2)*mst12*mw2)) + 3*
     lmsb22*(-mst16 + mw2*((-23 + 6*lmst12 + 2*lmw2)*mst14 + ((-23 + 2*lmst12 +
     6*lmw2)*mst12 - mw2)*mw2) + msb22*(mst14 + mw2*(2*(-3 + lmst12 + lmw2)*
     mst12 + mw2))) + 6*mw2*(msb22*mst12 + 2*(mst14 + mst12*mw2))*power2(lmsb22
     ) + 6*mst14*mw2*power2(lmst12) + mw2*(-3*lmw2*msb22*mw2 + msb22*mst12*(42
     - 9*lmw2 + power2(Pi)) + 3*mst14*(42 + lmw2 + power2(Pi)) + 3*mw2*(lmw2*
     mw2 + mst12*(42 - 14*lmw2 + 2*power2(lmw2) + power2(Pi)))))*power2(snb)*(
     mt2*(-1 + power2(cb))*(-1 + power2(snt)) + (-2*At*cb*mu*sb - power2(At)*(-
     1 + power2(cb)) + mu2*power2(cb))*power2(snt)))/(6.*msb22))*power2(
     DeltaInv(mw2,msb22,mst12)) - (invdmst*(mst12 - mst22)*mt2*(sb*power2(At)*(
     -((lmsb12 - lmw2)*msb22*mw2*(-1 + power2(snb))) + (lmsb22 - lmw2)*msb12*
     mw2*power2(snb) + msb12*msb22*(-3 + lmsb12 - 2*lmh2*(-1 + power2(sa)) + 2*
     lmH2*power2(sa) - lmsb12*power2(snb) + lmsb22*power2(snb)) + power2(cb)*(
     lmw2*msb22*mw2 - lmsb12*msb22*(mC2 - mw2)*(-1 + power2(snb)) + lmsb22*mC2*
     msb12*power2(snb) - lmsb22*msb12*mw2*power2(snb) + lmw2*msb12*mw2*power2(
     snb) - lmw2*msb22*mw2*power2(snb) - lmC2*mC2*(msb22 + msb12*power2(snb) -
     msb22*power2(snb)))) + cb*mu2*(2*ca*(-lmh2 + lmH2)*msb12*msb22*sa + cb*sb*
     (-(lmw2*msb22*mw2) + lmsb12*msb22*(mC2 - mw2)*(-1 + power2(snb)) - lmsb22*
     mC2*msb12*power2(snb) + lmsb22*msb12*mw2*power2(snb) - lmw2*msb12*mw2*
     power2(snb) + lmw2*msb22*mw2*power2(snb) + lmC2*mC2*(msb22 + msb12*power2(
     snb) - msb22*power2(snb)))) - At*mu*(2*ca*(-lmh2 + lmH2)*msb12*msb22*sa*sb
      + cb*(msb22*(-(lmC2*mC2*power2(sb)) + lmw2*mw2*(1 + power2(sb)) - lmsb12*
     (mw2 - mC2*power2(sb) + mw2*power2(sb)))*(-1 + power2(snb)) + msb12*((lmC2
     *mC2*power2(sb) - lmw2*mw2*(1 + power2(sb)) + lmsb22*(mw2 - mC2*power2(sb)
     + mw2*power2(sb)))*power2(snb) + msb22*(-3 + lmsb12 - 2*lmh2*(-1 + power2(
     sa)) + 2*lmH2*power2(sa) - lmsb12*power2(snb) + lmsb22*power2(snb)))) + (
     lmw2*msb22*mw2 - lmsb12*msb22*(mC2 - mw2)*(-1 + power2(snb)) + lmsb22*mC2*
     msb12*power2(snb) - lmsb22*msb12*mw2*power2(snb) + lmw2*msb12*mw2*power2(
     snb) - lmw2*msb22*mw2*power2(snb) - lmC2*mC2*(msb22 + msb12*power2(snb) -
     msb22*power2(snb)))*power3(cb))))/(6.*msb12*msb22*mst12*mst22*sb) - (mA2*(
     mu2 + 2*At*cb*mu*sb - mu2*power2(cb) + power2(At)*power2(cb))*power2(
     DeltaInv(mst22,mA2,mst12))*(42*mA2*mst12*mst22 - 9*lmst22*mA2*mst12*mst22
     + 126*mA2*mst24 - 42*lmst22*mA2*mst24 - 3*lmst22*mst12*mst24 + 3*lmst22*
     mst26 + 6*mA2*(2*mA2*mst22 + mst12*mst22 + 2*mst24)*power2(lmst12) + 6*mA2
     *mst24*power2(lmst22) + 126*mst22*power2(mA2) + 3*lmst22*mst22*power2(mA2)
     + 6*mst22*power2(lmA2)*power2(mA2) + 3*lmA2*mA2*((-3 + 2*lmst12 - 2*lmst22
     )*mst12*mst22 - mA2*(mst12 + 2*(7 - 3*lmst12 + lmst22)*mst22) + (1 + 2*
     lmst12 - 2*lmst22)*mst24 + power2(mA2)) + mA2*mst12*mst22*power2(Pi) + 3*
     mA2*mst24*power2(Pi) + 3*mst22*power2(mA2)*power2(Pi) - 3*lmst12*(-(mst12*
     mst24) + mA2*(-2*(-3 + lmst22)*mst12*mst22 + (23 - 6*lmst22)*mst24) +
     mst26 - (mst12 + (-23 + 2*lmst22)*mst22)*power2(mA2) + power3(mA2))))/(12.
     *mst12) + (DeltaInv(mst22,mA2,mst12)*(mu2 + 2*At*cb*mu*sb - mu2*power2(cb)
     + power2(At)*power2(cb))*(mst22*(2*mA2*mst12 + 5*lmst22*mA2*mst22 - lmst22
     *mst12*mst22 + lmst22*mst24) + lmA2*(mA2 - mst12 + 5*mst22)*power2(mA2) -
     lmst12*(mst22*(-(mst12*mst22) + mst24) - (mst12 - 5*mst22)*power2(mA2) + 5
     *mA2*power2(mst22) + power3(mA2))))/(24.*mst12*mst22) + power2(DeltaInv(
     msb12,mst12,mC2))*((cb*invdmst*mC2*mt2*(-(cb*mu2*sb) + cb*sb*power2(At) +
     At*mu*(-power2(cb) + power2(sb)))*(-1 + power2(snb))*(42*mC2*msb12*mst12 -
     9*lmst12*mC2*msb12*mst12 + 126*mC2*mst14 - 42*lmst12*mC2*mst14 - 3*lmst12*
     msb12*mst14 + 3*lmst12*mst16 + 6*mC2*(2*mC2*mst12 + msb12*mst12 + 2*mst14)
     *power2(lmsb12) + 6*mC2*mst14*power2(lmst12) + 126*mst12*power2(mC2) + 3*
     lmst12*mst12*power2(mC2) + 6*mst12*power2(lmC2)*power2(mC2) + 3*lmC2*mC2*(
     (-3 + 2*lmsb12 - 2*lmst12)*msb12*mst12 - mC2*(msb12 + 2*(7 - 3*lmsb12 +
     lmst12)*mst12) + (1 + 2*lmsb12 - 2*lmst12)*mst14 + power2(mC2)) + mC2*
     msb12*mst12*power2(Pi) + 3*mC2*mst14*power2(Pi) + 3*mst12*power2(mC2)*
     power2(Pi) - 3*lmsb12*(-(msb12*mst14) + mC2*(-2*(-3 + lmst12)*msb12*mst12
     + (23 - 6*lmst12)*mst14) + mst16 - (msb12 + (-23 + 2*lmst12)*mst12)*power2
     (mC2) + power3(mC2))))/(3.*msb12*sb) - (mC2*(-1 + power2(snb))*(-(mu2*
     power2(snt)) - 2*At*cb*mu*sb*power2(snt) + power2(cb)*(mt2*(-1 + power2(
     snt)) + (mu2 - power2(At))*power2(snt)))*(42*mC2*msb12*mst12 - 9*lmst12*
     mC2*msb12*mst12 + 126*mC2*mst14 - 42*lmst12*mC2*mst14 - 3*lmst12*msb12*
     mst14 + 3*lmst12*mst16 + 6*mC2*(2*mC2*mst12 + msb12*mst12 + 2*mst14)*
     power2(lmsb12) + 6*mC2*mst14*power2(lmst12) + 126*mst12*power2(mC2) + 3*
     lmst12*mst12*power2(mC2) + 6*mst12*power2(lmC2)*power2(mC2) + 3*lmC2*mC2*(
     (-3 + 2*lmsb12 - 2*lmst12)*msb12*mst12 - mC2*(msb12 + 2*(7 - 3*lmsb12 +
     lmst12)*mst12) + (1 + 2*lmsb12 - 2*lmst12)*mst14 + power2(mC2)) + mC2*
     msb12*mst12*power2(Pi) + 3*mC2*mst14*power2(Pi) + 3*mst12*power2(mC2)*
     power2(Pi) - 3*lmsb12*(-(msb12*mst14) + mC2*(-2*(-3 + lmst12)*msb12*mst12
     + (23 - 6*lmst12)*mst14) + mst16 - (msb12 + (-23 + 2*lmst12)*mst12)*power2
     (mC2) + power3(mC2))))/(6.*msb12)) + power2(DeltaInv(msb22,mst12,mC2))*((
     cb*invdmst*mC2*mt2*(cb*mu2*sb - cb*sb*power2(At) + At*mu*(power2(cb) -
     power2(sb)))*power2(snb)*(42*mC2*msb22*mst12 - 9*lmst12*mC2*msb22*mst12 +
     126*mC2*mst14 - 42*lmst12*mC2*mst14 - 3*lmst12*msb22*mst14 + 3*lmst12*
     mst16 + 6*mC2*(2*mC2*mst12 + msb22*mst12 + 2*mst14)*power2(lmsb22) + 6*mC2
     *mst14*power2(lmst12) + 126*mst12*power2(mC2) + 3*lmst12*mst12*power2(mC2)
     + 6*mst12*power2(lmC2)*power2(mC2) + 3*lmC2*mC2*((-3 + 2*lmsb22 - 2*lmst12
     )*msb22*mst12 - mC2*(msb22 + 2*(7 - 3*lmsb22 + lmst12)*mst12) + (1 + 2*
     lmsb22 - 2*lmst12)*mst14 + power2(mC2)) + mC2*msb22*mst12*power2(Pi) + 3*
     mC2*mst14*power2(Pi) + 3*mst12*power2(mC2)*power2(Pi) - 3*lmsb22*(-(msb22*
     mst14) + mC2*(-2*(-3 + lmst12)*msb22*mst12 + (23 - 6*lmst12)*mst14) +
     mst16 - (msb22 + (-23 + 2*lmst12)*mst12)*power2(mC2) + power3(mC2))))/(3.*
     msb22*sb) + (mC2*power2(snb)*(-(mu2*power2(snt)) - 2*At*cb*mu*sb*power2(
     snt) + power2(cb)*(mt2*(-1 + power2(snt)) + (mu2 - power2(At))*power2(snt)
     ))*(42*mC2*msb22*mst12 - 9*lmst12*mC2*msb22*mst12 + 126*mC2*mst14 - 42*
     lmst12*mC2*mst14 - 3*lmst12*msb22*mst14 + 3*lmst12*mst16 + 6*mC2*(2*mC2*
     mst12 + msb22*mst12 + 2*mst14)*power2(lmsb22) + 6*mC2*mst14*power2(lmst12)
     + 126*mst12*power2(mC2) + 3*lmst12*mst12*power2(mC2) + 6*mst12*power2(lmC2
     )*power2(mC2) + 3*lmC2*mC2*((-3 + 2*lmsb22 - 2*lmst12)*msb22*mst12 - mC2*(
     msb22 + 2*(7 - 3*lmsb22 + lmst12)*mst12) + (1 + 2*lmsb22 - 2*lmst12)*mst14
      + power2(mC2)) + mC2*msb22*mst12*power2(Pi) + 3*mC2*mst14*power2(Pi) + 3*
     mst12*power2(mC2)*power2(Pi) - 3*lmsb22*(-(msb22*mst14) + mC2*(-2*(-3 +
     lmst12)*msb22*mst12 + (23 - 6*lmst12)*mst14) + mst16 - (msb22 + (-23 + 2*
     lmst12)*mst12)*power2(mC2) + power3(mC2))))/(6.*msb22)) + power2(DeltaInv(
     mst22,msb12,mC2))*((cb*invdmst*mC2*mt2*(cb*mu2*sb - cb*sb*power2(At) + At*
     mu*(power2(cb) - power2(sb)))*(-1 + power2(snb))*(42*mC2*msb12*mst22 - 9*
     lmst22*mC2*msb12*mst22 + 126*mC2*mst24 - 42*lmst22*mC2*mst24 - 3*lmst22*
     msb12*mst24 + 3*lmst22*mst26 + 6*mC2*(2*mC2*mst22 + msb12*mst22 + 2*mst24)
     *power2(lmsb12) + 6*mC2*mst24*power2(lmst22) + 126*mst22*power2(mC2) + 3*
     lmst22*mst22*power2(mC2) + 6*mst22*power2(lmC2)*power2(mC2) + 3*lmC2*mC2*(
     (-3 + 2*lmsb12 - 2*lmst22)*msb12*mst22 - mC2*(msb12 + 2*(7 - 3*lmsb12 +
     lmst22)*mst22) + (1 + 2*lmsb12 - 2*lmst22)*mst24 + power2(mC2)) + mC2*
     msb12*mst22*power2(Pi) + 3*mC2*mst24*power2(Pi) + 3*mst22*power2(mC2)*
     power2(Pi) - 3*lmsb12*(-(msb12*mst24) + mC2*(-2*(-3 + lmst22)*msb12*mst22
     + (23 - 6*lmst22)*mst24) + mst26 - (msb12 + (-23 + 2*lmst22)*mst22)*power2
     (mC2) + power3(mC2))))/(3.*msb12*sb) + (mC2*(-1 + power2(snb))*(-2*At*cb*
     mu*sb*(-1 + power2(snt)) + mu2*(-1 + power2(cb))*(-1 + power2(snt)) -
     power2(At)*power2(cb)*(-1 + power2(snt)) + mt2*power2(cb)*power2(snt))*(42
     *mC2*msb12*mst22 - 9*lmst22*mC2*msb12*mst22 + 126*mC2*mst24 - 42*lmst22*
     mC2*mst24 - 3*lmst22*msb12*mst24 + 3*lmst22*mst26 + 6*mC2*(2*mC2*mst22 +
     msb12*mst22 + 2*mst24)*power2(lmsb12) + 6*mC2*mst24*power2(lmst22) + 126*
     mst22*power2(mC2) + 3*lmst22*mst22*power2(mC2) + 6*mst22*power2(lmC2)*
     power2(mC2) + 3*lmC2*mC2*((-3 + 2*lmsb12 - 2*lmst22)*msb12*mst22 - mC2*(
     msb12 + 2*(7 - 3*lmsb12 + lmst22)*mst22) + (1 + 2*lmsb12 - 2*lmst22)*mst24
      + power2(mC2)) + mC2*msb12*mst22*power2(Pi) + 3*mC2*mst24*power2(Pi) + 3*
     mst22*power2(mC2)*power2(Pi) - 3*lmsb12*(-(msb12*mst24) + mC2*(-2*(-3 +
     lmst22)*msb12*mst22 + (23 - 6*lmst22)*mst24) + mst26 - (msb12 + (-23 + 2*
     lmst22)*mst22)*power2(mC2) + power3(mC2))))/(6.*msb12)) + power2(DeltaInv(
     mst22,msb22,mC2))*((cb*invdmst*mC2*mt2*(-(cb*mu2*sb) + cb*sb*power2(At) +
     At*mu*(-power2(cb) + power2(sb)))*power2(snb)*(42*mC2*msb22*mst22 - 9*
     lmst22*mC2*msb22*mst22 + 126*mC2*mst24 - 42*lmst22*mC2*mst24 - 3*lmst22*
     msb22*mst24 + 3*lmst22*mst26 + 6*mC2*(2*mC2*mst22 + msb22*mst22 + 2*mst24)
     *power2(lmsb22) + 6*mC2*mst24*power2(lmst22) + 126*mst22*power2(mC2) + 3*
     lmst22*mst22*power2(mC2) + 6*mst22*power2(lmC2)*power2(mC2) + 3*lmC2*mC2*(
     (-3 + 2*lmsb22 - 2*lmst22)*msb22*mst22 - mC2*(msb22 + 2*(7 - 3*lmsb22 +
     lmst22)*mst22) + (1 + 2*lmsb22 - 2*lmst22)*mst24 + power2(mC2)) + mC2*
     msb22*mst22*power2(Pi) + 3*mC2*mst24*power2(Pi) + 3*mst22*power2(mC2)*
     power2(Pi) - 3*lmsb22*(-(msb22*mst24) + mC2*(-2*(-3 + lmst22)*msb22*mst22
     + (23 - 6*lmst22)*mst24) + mst26 - (msb22 + (-23 + 2*lmst22)*mst22)*power2
     (mC2) + power3(mC2))))/(3.*msb22*sb) + (mC2*power2(snb)*(2*At*cb*mu*sb*(-1
      + power2(snt)) - mu2*(-1 + power2(cb))*(-1 + power2(snt)) + power2(At)*
     power2(cb)*(-1 + power2(snt)) - mt2*power2(cb)*power2(snt))*(42*mC2*msb22*
     mst22 - 9*lmst22*mC2*msb22*mst22 + 126*mC2*mst24 - 42*lmst22*mC2*mst24 - 3
     *lmst22*msb22*mst24 + 3*lmst22*mst26 + 6*mC2*(2*mC2*mst22 + msb22*mst22 +
     2*mst24)*power2(lmsb22) + 6*mC2*mst24*power2(lmst22) + 126*mst22*power2(
     mC2) + 3*lmst22*mst22*power2(mC2) + 6*mst22*power2(lmC2)*power2(mC2) + 3*
     lmC2*mC2*((-3 + 2*lmsb22 - 2*lmst22)*msb22*mst22 - mC2*(msb22 + 2*(7 - 3*
     lmsb22 + lmst22)*mst22) + (1 + 2*lmsb22 - 2*lmst22)*mst24 + power2(mC2)) +
     mC2*msb22*mst22*power2(Pi) + 3*mC2*mst24*power2(Pi) + 3*mst22*power2(mC2)*
     power2(Pi) - 3*lmsb22*(-(msb22*mst24) + mC2*(-2*(-3 + lmst22)*msb22*mst22
     + (23 - 6*lmst22)*mst24) + mst26 - (msb22 + (-23 + 2*lmst22)*mst22)*power2
     (mC2) + power3(mC2))))/(6.*msb22)) + DeltaInv(msb12,mst12,mC2)*((cb*
     invdmst*mt2*(cb*mu2*sb - cb*sb*power2(At) + At*mu*(power2(cb) - power2(sb)
     ))*(-1 + power2(snb))*(mst12*(2*mC2*msb12 + 5*lmst12*mC2*mst12 - lmst12*
     msb12*mst12 + lmst12*mst14) + lmC2*(mC2 - msb12 + 5*mst12)*power2(mC2) -
     lmsb12*(mst12*(-(msb12*mst12) + mst14) - (msb12 - 5*mst12)*power2(mC2) + 5
     *mC2*power2(mst12) + power3(mC2))))/(6.*msb12*mst12*sb) + ((-1 + power2(
     snb))*(-(mu2*power2(snt)) - 2*At*cb*mu*sb*power2(snt) + power2(cb)*(mt2*(-
     1 + power2(snt)) + (mu2 - power2(At))*power2(snt)))*(mst12*(2*mC2*msb12 +
     5*lmst12*mC2*mst12 - lmst12*msb12*mst12 + lmst12*mst14) + lmC2*(mC2 -
     msb12 + 5*mst12)*power2(mC2) - lmsb12*(mst12*(-(msb12*mst12) + mst14) - (
     msb12 - 5*mst12)*power2(mC2) + 5*mC2*power2(mst12) + power3(mC2))))/(12.*
     msb12*mst12)) + DeltaInv(msb22,mst12,mC2)*((cb*invdmst*mt2*(-(cb*mu2*sb) +
     cb*sb*power2(At) + At*mu*(-power2(cb) + power2(sb)))*power2(snb)*(mst12*(2
     *mC2*msb22 + 5*lmst12*mC2*mst12 - lmst12*msb22*mst12 + lmst12*mst14) +
     lmC2*(mC2 - msb22 + 5*mst12)*power2(mC2) - lmsb22*(mst12*(-(msb22*mst12) +
     mst14) - (msb22 - 5*mst12)*power2(mC2) + 5*mC2*power2(mst12) + power3(mC2)
     )))/(6.*msb22*mst12*sb) - (power2(snb)*(-(mu2*power2(snt)) - 2*At*cb*mu*sb
     *power2(snt) + power2(cb)*(mt2*(-1 + power2(snt)) + (mu2 - power2(At))*
     power2(snt)))*(mst12*(2*mC2*msb22 + 5*lmst12*mC2*mst12 - lmst12*msb22*
     mst12 + lmst12*mst14) + lmC2*(mC2 - msb22 + 5*mst12)*power2(mC2) - lmsb22*
     (mst12*(-(msb22*mst12) + mst14) - (msb22 - 5*mst12)*power2(mC2) + 5*mC2*
     power2(mst12) + power3(mC2))))/(12.*msb22*mst12)) + DeltaInv(mst22,msb12,
     mC2)*((cb*invdmst*mt2*(-(cb*mu2*sb) + cb*sb*power2(At) + At*mu*(-power2(cb
     ) + power2(sb)))*(-1 + power2(snb))*(mst22*(2*mC2*msb12 + 5*lmst22*mC2*
     mst22 - lmst22*msb12*mst22 + lmst22*mst24) + lmC2*(mC2 - msb12 + 5*mst22)*
     power2(mC2) - lmsb12*(mst22*(-(msb12*mst22) + mst24) - (msb12 - 5*mst22)*
     power2(mC2) + 5*mC2*power2(mst22) + power3(mC2))))/(6.*msb12*mst22*sb) + (
     (-1 + power2(snb))*(2*At*cb*mu*sb*(-1 + power2(snt)) - mu2*(-1 + power2(cb
     ))*(-1 + power2(snt)) + power2(At)*power2(cb)*(-1 + power2(snt)) - mt2*
     power2(cb)*power2(snt))*(mst22*(2*mC2*msb12 + 5*lmst22*mC2*mst22 - lmst22*
     msb12*mst22 + lmst22*mst24) + lmC2*(mC2 - msb12 + 5*mst22)*power2(mC2) -
     lmsb12*(mst22*(-(msb12*mst22) + mst24) - (msb12 - 5*mst22)*power2(mC2) + 5
     *mC2*power2(mst22) + power3(mC2))))/(12.*msb12*mst22)) + DeltaInv(mst22,
     msb22,mC2)*((cb*invdmst*mt2*(cb*mu2*sb - cb*sb*power2(At) + At*mu*(power2(
     cb) - power2(sb)))*power2(snb)*(mst22*(2*mC2*msb22 + 5*lmst22*mC2*mst22 -
     lmst22*msb22*mst22 + lmst22*mst24) + lmC2*(mC2 - msb22 + 5*mst22)*power2(
     mC2) - lmsb22*(mst22*(-(msb22*mst22) + mst24) - (msb22 - 5*mst22)*power2(
     mC2) + 5*mC2*power2(mst22) + power3(mC2))))/(6.*msb22*mst22*sb) - (power2(
     snb)*(2*At*cb*mu*sb*(-1 + power2(snt)) - mu2*(-1 + power2(cb))*(-1 +
     power2(snt)) + power2(At)*power2(cb)*(-1 + power2(snt)) - mt2*power2(cb)*
     power2(snt))*(mst22*(2*mC2*msb22 + 5*lmst22*mC2*mst22 - lmst22*msb22*mst22
      + lmst22*mst24) + lmC2*(mC2 - msb22 + 5*mst22)*power2(mC2) - lmsb22*(
     mst22*(-(msb22*mst22) + mst24) - (msb22 - 5*mst22)*power2(mC2) + 5*mC2*
     power2(mst22) + power3(mC2))))/(12.*msb22*mst22)) + (mh2*(-2*At*ca*mu*sa +
     power2(At)*(-1 + power2(sa)) - mu2*power2(sa))*power2(DeltaInv(mst22,mh2,
     mst12))*power2(1 - 2*power2(snt))*(42*mh2*mst12*mst22 - 9*lmst22*mh2*mst12
     *mst22 + 126*mh2*mst24 - 42*lmst22*mh2*mst24 - 3*lmst22*mst12*mst24 + 3*
     lmst22*mst26 + 6*mh2*(2*mh2*mst22 + mst12*mst22 + 2*mst24)*power2(lmst12)
     + 6*mh2*mst24*power2(lmst22) + 126*mst22*power2(mh2) + 3*lmst22*mst22*
     power2(mh2) + 6*mst22*power2(lmh2)*power2(mh2) + 3*lmh2*mh2*((-3 + 2*
     lmst12 - 2*lmst22)*mst12*mst22 - mh2*(mst12 + 2*(7 - 3*lmst12 + lmst22)*
     mst22) + (1 + 2*lmst12 - 2*lmst22)*mst24 + power2(mh2)) + mh2*mst12*mst22*
     power2(Pi) + 3*mh2*mst24*power2(Pi) + 3*mst22*power2(mh2)*power2(Pi) - 3*
     lmst12*(-(mst12*mst24) + mh2*(-2*(-3 + lmst22)*mst12*mst22 + (23 - 6*
     lmst22)*mst24) + mst26 - (mst12 + (-23 + 2*lmst22)*mst22)*power2(mh2) +
     power3(mh2))))/(12.*mst12) - (DeltaInv(mst22,mh2,mst12)*(-2*At*ca*mu*sa +
     power2(At)*(-1 + power2(sa)) - mu2*power2(sa))*power2(1 - 2*power2(snt))*(
     mst22*(2*mh2*mst12 + 5*lmst22*mh2*mst22 - lmst22*mst12*mst22 + lmst22*
     mst24) + lmh2*(mh2 - mst12 + 5*mst22)*power2(mh2) - lmst12*(mst22*(-(mst12
     *mst22) + mst24) - (mst12 - 5*mst22)*power2(mh2) + 5*mh2*power2(mst22) +
     power3(mh2))))/(24.*mst12*mst22) + (mH2*(At*sa*(2*ca*mu - At*sa) + mu2*(-1
      + power2(sa)))*power2(DeltaInv(mH2,mst22,mst12))*power2(1 - 2*power2(snt)
     )*(42*mH2*mst12*mst22 - 9*lmst22*mH2*mst12*mst22 + 126*mH2*mst24 - 42*
     lmst22*mH2*mst24 - 3*lmst22*mst12*mst24 + 3*lmst22*mst26 + 6*mH2*(2*mH2*
     mst22 + mst12*mst22 + 2*mst24)*power2(lmst12) + 6*mH2*mst24*power2(lmst22)
     + 126*mst22*power2(mH2) + 3*lmst22*mst22*power2(mH2) + 6*mst22*power2(lmH2
     )*power2(mH2) + 3*lmH2*mH2*((-3 + 2*lmst12 - 2*lmst22)*mst12*mst22 - mH2*(
     mst12 + 2*(7 - 3*lmst12 + lmst22)*mst22) + (1 + 2*lmst12 - 2*lmst22)*mst24
      + power2(mH2)) + mH2*mst12*mst22*power2(Pi) + 3*mH2*mst24*power2(Pi) + 3*
     mst22*power2(mH2)*power2(Pi) - 3*lmst12*(-(mst12*mst24) + mH2*(-2*(-3 +
     lmst22)*mst12*mst22 + (23 - 6*lmst22)*mst24) + mst26 - (mst12 + (-23 + 2*
     lmst22)*mst22)*power2(mH2) + power3(mH2))))/(12.*mst12) - (DeltaInv(mH2,
     mst22,mst12)*(At*sa*(2*ca*mu - At*sa) + mu2*(-1 + power2(sa)))*power2(1 -
     2*power2(snt))*(mst22*(2*mH2*mst12 + 5*lmst22*mH2*mst22 - lmst22*mst12*
     mst22 + lmst22*mst24) + lmH2*(mH2 - mst12 + 5*mst22)*power2(mH2) - lmst12*
     (mst22*(-(mst12*mst22) + mst24) - (mst12 - 5*mst22)*power2(mH2) + 5*mH2*
     power2(mst22) + power3(mH2))))/(24.*mst12*mst22) - (3*mh2*(42 + 8*lmh2*(-3
      + lmt2) - 12*lmt2 + 4*power2(lmh2) + power2(Pi))*(-1 + power2(sa))*power2
     (DeltaInv(mh2,mt2,mt2))*power3(mt2))/2. + (3*mH2*(42 + 8*lmH2*(-3 + lmt2)
     - 12*lmt2 + 4*power2(lmH2) + power2(Pi))*power2(sa)*power2(DeltaInv(mH2,
     mt2,mt2))*power3(mt2))/2. - (power2(cb)*(12*mt2*(42 + 8*lmA2*(-3 + lmt2) -
     12*lmt2 + 4*power2(lmA2) + power2(Pi)) + mA2*(714 - 228*lmt2 + 48*lmA2*(-8
      + 3*lmt2) + 72*power2(lmA2) - 12*power2(lmt2) + 17*power2(Pi)))*power2(
     DeltaInv(mt2,mt2,mA2))*power3(mt2))/6. + ((-1 + power2(cb))*(12*mt2*(42 -
     24*lmz2 + 4*lmt2*(-3 + 2*lmz2) + 4*power2(lmz2) + power2(Pi)) + mz2*(714 -
     384*lmz2 + 12*lmt2*(-19 + 12*lmz2) - 12*power2(lmt2) + 72*power2(lmz2) +
     17*power2(Pi)))*power2(DeltaInv(mt2,mt2,mz2))*power3(mt2))/6. + DeltaInv(
     mw2,msb12,mst12)*((invdmst*mt2*(sb*power2(At)*(-1 + power2(cb)) - mu2*sb*
     power2(cb) + At*cb*mu*(1 - power2(cb) + power2(sb)))*(-1 + power2(snb))*(
     lmst12*mst12*(-(msb12*mst12) + mst14 + 5*mst12*mw2) + mw2*(lmw2*mw2*(5*
     mst12 + mw2) + msb12*(2*mst12 - lmw2*mw2)) - lmsb12*(mst12*mst14 + 5*mw2*
     power2(mst12) + 5*mst12*power2(mw2) - msb12*(power2(mst12) + power2(mw2))
     + power3(mw2))))/(6.*msb12*mst12*sb) - ((-1 + power2(snb))*(mt2*(-1 +
     power2(cb))*(-1 + power2(snt)) + (-2*At*cb*mu*sb - power2(At)*(-1 + power2
     (cb)) + mu2*power2(cb))*power2(snt))*(lmst12*mst12*(-(msb12*mst12) + mst14
      + 5*mst12*mw2) + mw2*(lmw2*mw2*(5*mst12 + mw2) + msb12*(2*mst12 - lmw2*
     mw2)) - lmsb12*(mst12*mst14 + 5*mw2*power2(mst12) + 5*mst12*power2(mw2) -
     msb12*(power2(mst12) + power2(mw2)) + power3(mw2))))/(12.*msb12*mst12)) +
     DeltaInv(mw2,msb22,mst12)*(-(invdmst*mt2*(sb*power2(At)*(-1 + power2(cb))
     - mu2*sb*power2(cb) + At*cb*mu*(1 - power2(cb) + power2(sb)))*power2(snb)*
     (lmst12*mst12*(-(msb22*mst12) + mst14 + 5*mst12*mw2) + mw2*(lmw2*mw2*(5*
     mst12 + mw2) + msb22*(2*mst12 - lmw2*mw2)) - lmsb22*(mst12*mst14 + 5*mw2*
     power2(mst12) + 5*mst12*power2(mw2) - msb22*(power2(mst12) + power2(mw2))
     + power3(mw2))))/(6.*msb22*mst12*sb) + (power2(snb)*(mt2*(-1 + power2(cb))
     *(-1 + power2(snt)) + (-2*At*cb*mu*sb - power2(At)*(-1 + power2(cb)) + mu2
     *power2(cb))*power2(snt))*(lmst12*mst12*(-(msb22*mst12) + mst14 + 5*mst12*
     mw2) + mw2*(lmw2*mw2*(5*mst12 + mw2) + msb22*(2*mst12 - lmw2*mw2)) -
     lmsb22*(mst12*mst14 + 5*mw2*power2(mst12) + 5*mst12*power2(mw2) - msb22*(
     power2(mst12) + power2(mw2)) + power3(mw2))))/(12.*msb22*mst12)) +
     DeltaInv(mst22,mw2,msb12)*(-(invdmst*mt2*(sb*power2(At)*(-1 + power2(cb))
     - mu2*sb*power2(cb) + At*cb*mu*(1 - power2(cb) + power2(sb)))*(-1 + power2
     (snb))*(lmst22*mst22*(-(msb12*mst22) + mst24 + 5*mst22*mw2) + mw2*(lmw2*
     mw2*(5*mst22 + mw2) + msb12*(2*mst22 - lmw2*mw2)) - lmsb12*(mst22*mst24 +
     5*mw2*power2(mst22) + 5*mst22*power2(mw2) - msb12*(power2(mst22) + power2(
     mw2)) + power3(mw2))))/(6.*msb12*mst22*sb) - ((-1 + power2(snb))*(2*At*cb*
     mu*sb*(-1 + power2(snt)) + power2(At)*(-1 + power2(cb))*(-1 + power2(snt))
     + mt2*power2(snt) + power2(cb)*(mu2 - mt2*power2(snt) - mu2*power2(snt)))*
     (lmst22*mst22*(-(msb12*mst22) + mst24 + 5*mst22*mw2) + mw2*(lmw2*mw2*(5*
     mst22 + mw2) + msb12*(2*mst22 - lmw2*mw2)) - lmsb12*(mst22*mst24 + 5*mw2*
     power2(mst22) + 5*mst22*power2(mw2) - msb12*(power2(mst22) + power2(mw2))
     + power3(mw2))))/(12.*msb12*mst22)) + DeltaInv(mst22,mw2,msb22)*((invdmst*
     mt2*(sb*power2(At)*(-1 + power2(cb)) - mu2*sb*power2(cb) + At*cb*mu*(1 -
     power2(cb) + power2(sb)))*power2(snb)*(lmst22*mst22*(-(msb22*mst22) +
     mst24 + 5*mst22*mw2) + mw2*(lmw2*mw2*(5*mst22 + mw2) + msb22*(2*mst22 -
     lmw2*mw2)) - lmsb22*(mst22*mst24 + 5*mw2*power2(mst22) + 5*mst22*power2(
     mw2) - msb22*(power2(mst22) + power2(mw2)) + power3(mw2))))/(6.*msb22*
     mst22*sb) + (power2(snb)*(2*At*cb*mu*sb*(-1 + power2(snt)) + power2(At)*(-
     1 + power2(cb))*(-1 + power2(snt)) + mt2*power2(snt) + power2(cb)*(mu2 -
     mt2*power2(snt) - mu2*power2(snt)))*(lmst22*mst22*(-(msb22*mst22) + mst24
     + 5*mst22*mw2) + mw2*(lmw2*mw2*(5*mst22 + mw2) + msb22*(2*mst22 - lmw2*mw2
     )) - lmsb22*(mst22*mst24 + 5*mw2*power2(mst22) + 5*mst22*power2(mw2) -
     msb22*(power2(mst22) + power2(mw2)) + power3(mw2))))/(12.*msb22*mst22)) -
     (DeltaInv(mst22,mz2,mst12)*(2*At*cb*mu*sb + power2(At)*(-1 + power2(cb)) -
     mu2*power2(cb))*(lmst22*mst22*(-(mst12*mst22) + mst24 + 5*mst22*mz2) + mz2
     *(lmz2*mz2*(5*mst22 + mz2) + mst12*(2*mst22 - lmz2*mz2)) - lmst12*(mst22*
     mst24 + 5*mz2*power2(mst22) + 5*mst22*power2(mz2) - mst12*(power2(mst22) +
     power2(mz2)) + power3(mz2))))/(24.*mst12*mst22) + Fin3(mh2,mst12,mst12,Q2)
     *(((4*invdmst*mst14*mt2*(ca*cb*mu2*sa - At*mu*(ca*sa*sb + cb*(-1 + power2(
     sa))) + sb*power2(At)*(-1 + power2(sa))))/sb + 2*mst14*(mt2*(-1 + power2(
     sa)) - (-2*At*ca*mu*sa + power2(At)*(-1 + power2(sa)) - mu2*power2(sa))*(-
     1 + power2(snt))*power2(snt)))*power2(DeltaInv(mh2,mst12,mst12)) + ((16*
     invdmst*mh2*mst16*mt2*(ca*cb*mu2*sa - At*mu*(ca*sa*sb + cb*(-1 + power2(sa
     ))) + sb*power2(At)*(-1 + power2(sa))))/sb + 8*mh2*mst16*(mt2*(-1 + power2
     (sa)) - (-2*At*ca*mu*sa + power2(At)*(-1 + power2(sa)) - mu2*power2(sa))*(
     -1 + power2(snt))*power2(snt)))*power3(DeltaInv(mh2,mst12,mst12))) + Fin3(
     mH2,mst12,mst12,Q2)*(((4*invdmst*mst14*mt2*sa*(At*sa*(cb*mu - At*sb) + ca*
     (-(cb*mu2) + At*mu*sb)))/sb - 2*mst14*(mt2*power2(sa) + (At*sa*(2*ca*mu -
     At*sa) + mu2*(-1 + power2(sa)))*(-1 + power2(snt))*power2(snt)))*power2(
     DeltaInv(mH2,mst12,mst12)) + ((16*invdmst*mH2*mst16*mt2*sa*(At*sa*(cb*mu -
     At*sb) + ca*(-(cb*mu2) + At*mu*sb)))/sb - 8*mH2*mst16*(mt2*power2(sa) + (
     At*sa*(2*ca*mu - At*sa) + mu2*(-1 + power2(sa)))*(-1 + power2(snt))*power2
     (snt)))*power3(DeltaInv(mH2,mst12,mst12))) + Fin3(mH2,mst22,mst12,Q2)*((
     mst22*(At*sa*(2*ca*mu - At*sa) + mu2*(-1 + power2(sa)))*power2(DeltaInv(
     mH2,mst22,mst12))*power2(mH2 - 2*mH2*power2(snt)))/mst12 - ((-2*mst12*
     mst24 - 2*mH2*(mst12*mst22 + mst24) + mst26 + mst22*power2(mH2))*(At*sa*(2
     *ca*mu - At*sa) + mu2*(-1 + power2(sa)))*power2(mH2 - 2*mH2*power2(snt))*
     power3(DeltaInv(mH2,mst22,mst12)))/mst12) + Fin3(mH2,mst22,mst22,Q2)*(((4*
     invdmst*mst24*mt2*sa*(At*sa*(-(cb*mu) + At*sb) + ca*(cb*mu2 - At*mu*sb)))/
     sb - 2*mst24*(mt2*power2(sa) + (At*sa*(2*ca*mu - At*sa) + mu2*(-1 + power2
     (sa)))*(-1 + power2(snt))*power2(snt)))*power2(DeltaInv(mH2,mst22,mst22))
     + ((16*invdmst*mH2*mst26*mt2*sa*(At*sa*(-(cb*mu) + At*sb) + ca*(cb*mu2 -
     At*mu*sb)))/sb - 8*mH2*mst26*(mt2*power2(sa) + (At*sa*(2*ca*mu - At*sa) +
     mu2*(-1 + power2(sa)))*(-1 + power2(snt))*power2(snt)))*power3(DeltaInv(
     mH2,mst22,mst22))) + Fin3(msb12,mst12,mC2,Q2)*(((4*cb*invdmst*mst12*mt2*
     power2(mC2)*(-(cb*mu2*sb) + cb*sb*power2(At) + At*mu*(-power2(cb) + power2
     (sb)))*(-1 + power2(snb)))/(msb12*sb) - (2*mst12*power2(mC2)*(-1 + power2(
     snb))*(-(mu2*power2(snt)) - 2*At*cb*mu*sb*power2(snt) + power2(cb)*(mt2*(-
     1 + power2(snt)) + (mu2 - power2(At))*power2(snt))))/msb12)*power2(
     DeltaInv(msb12,mst12,mC2)) + ((4*cb*invdmst*mt2*power2(mC2)*(-2*msb12*
     mst14 - 2*mC2*(msb12*mst12 + mst14) + mst16 + mst12*power2(mC2))*(cb*mu2*
     sb - cb*sb*power2(At) + At*mu*(power2(cb) - power2(sb)))*(-1 + power2(snb)
     ))/(msb12*sb) + (2*power2(mC2)*(-2*msb12*mst14 - 2*mC2*(msb12*mst12 +
     mst14) + mst16 + mst12*power2(mC2))*(-1 + power2(snb))*(-(mu2*power2(snt))
     - 2*At*cb*mu*sb*power2(snt) + power2(cb)*(mt2*(-1 + power2(snt)) + (mu2 -
     power2(At))*power2(snt))))/msb12)*power3(DeltaInv(msb12,mst12,mC2))) +
     Fin3(msb22,mst12,mC2,Q2)*(((4*cb*invdmst*mst12*mt2*power2(mC2)*(cb*mu2*sb
     - cb*sb*power2(At) + At*mu*(power2(cb) - power2(sb)))*power2(snb))/(msb22*
     sb) + (2*mst12*power2(mC2)*power2(snb)*(-(mu2*power2(snt)) - 2*At*cb*mu*sb
     *power2(snt) + power2(cb)*(mt2*(-1 + power2(snt)) + (mu2 - power2(At))*
     power2(snt))))/msb22)*power2(DeltaInv(msb22,mst12,mC2)) + ((4*cb*invdmst*
     mt2*power2(mC2)*(-2*msb22*mst14 - 2*mC2*(msb22*mst12 + mst14) + mst16 +
     mst12*power2(mC2))*(-(cb*mu2*sb) + cb*sb*power2(At) + At*mu*(-power2(cb) +
     power2(sb)))*power2(snb))/(msb22*sb) - (2*power2(mC2)*(-2*msb22*mst14 - 2*
     mC2*(msb22*mst12 + mst14) + mst16 + mst12*power2(mC2))*power2(snb)*(-(mu2*
     power2(snt)) - 2*At*cb*mu*sb*power2(snt) + power2(cb)*(mt2*(-1 + power2(
     snt)) + (mu2 - power2(At))*power2(snt))))/msb22)*power3(DeltaInv(msb22,
     mst12,mC2))) + Fin3(mst22,mA2,mst12,Q2)*(-((mst22*(mu2 + 2*At*cb*mu*sb -
     mu2*power2(cb) + power2(At)*power2(cb))*power2(mA2)*power2(DeltaInv(mst22,
     mA2,mst12)))/mst12) + ((mu2 + 2*At*cb*mu*sb - mu2*power2(cb) + power2(At)*
     power2(cb))*power2(mA2)*(-2*mst12*mst24 - 2*mA2*(mst12*mst22 + mst24) +
     mst26 + mst22*power2(mA2))*power3(DeltaInv(mst22,mA2,mst12)))/mst12) +
     Fin3(mst22,mh2,mst12,Q2)*((mst22*(-2*At*ca*mu*sa + power2(At)*(-1 + power2
     (sa)) - mu2*power2(sa))*power2(DeltaInv(mst22,mh2,mst12))*power2(mh2 - 2*
     mh2*power2(snt)))/mst12 - ((-2*mst12*mst24 - 2*mh2*(mst12*mst22 + mst24) +
     mst26 + mst22*power2(mh2))*(-2*At*ca*mu*sa + power2(At)*(-1 + power2(sa))
     - mu2*power2(sa))*power2(mh2 - 2*mh2*power2(snt))*power3(DeltaInv(mst22,
     mh2,mst12)))/mst12) + Fin3(mst22,msb12,mC2,Q2)*(((4*cb*invdmst*mst22*mt2*
     power2(mC2)*(cb*mu2*sb - cb*sb*power2(At) + At*mu*(power2(cb) - power2(sb)
     ))*(-1 + power2(snb)))/(msb12*sb) + (2*mst22*power2(mC2)*(-1 + power2(snb)
     )*(-2*At*cb*mu*sb*(-1 + power2(snt)) + mu2*(-1 + power2(cb))*(-1 + power2(
     snt)) - power2(At)*power2(cb)*(-1 + power2(snt)) + mt2*power2(cb)*power2(
     snt)))/msb12)*power2(DeltaInv(mst22,msb12,mC2)) + ((4*cb*invdmst*mt2*
     power2(mC2)*(-2*msb12*mst24 - 2*mC2*(msb12*mst22 + mst24) + mst26 + mst22*
     power2(mC2))*(-(cb*mu2*sb) + cb*sb*power2(At) + At*mu*(-power2(cb) +
     power2(sb)))*(-1 + power2(snb)))/(msb12*sb) + (2*power2(mC2)*(-2*msb12*
     mst24 - 2*mC2*(msb12*mst22 + mst24) + mst26 + mst22*power2(mC2))*(-1 +
     power2(snb))*(2*At*cb*mu*sb*(-1 + power2(snt)) - mu2*(-1 + power2(cb))*(-1
      + power2(snt)) + power2(At)*power2(cb)*(-1 + power2(snt)) - mt2*power2(cb
     )*power2(snt)))/msb12)*power3(DeltaInv(mst22,msb12,mC2))) + Fin3(mst22,
     msb22,mC2,Q2)*(((4*cb*invdmst*mst22*mt2*power2(mC2)*(-(cb*mu2*sb) + cb*sb*
     power2(At) + At*mu*(-power2(cb) + power2(sb)))*power2(snb))/(msb22*sb) + (
     2*mst22*power2(mC2)*power2(snb)*(2*At*cb*mu*sb*(-1 + power2(snt)) - mu2*(-
     1 + power2(cb))*(-1 + power2(snt)) + power2(At)*power2(cb)*(-1 + power2(
     snt)) - mt2*power2(cb)*power2(snt)))/msb22)*power2(DeltaInv(mst22,msb22,
     mC2)) + ((4*cb*invdmst*mt2*power2(mC2)*(-2*msb22*mst24 - 2*mC2*(msb22*
     mst22 + mst24) + mst26 + mst22*power2(mC2))*(cb*mu2*sb - cb*sb*power2(At)
     + At*mu*(power2(cb) - power2(sb)))*power2(snb))/(msb22*sb) + (2*power2(mC2
     )*(-2*msb22*mst24 - 2*mC2*(msb22*mst22 + mst24) + mst26 + mst22*power2(mC2
     ))*power2(snb)*(-2*At*cb*mu*sb*(-1 + power2(snt)) + mu2*(-1 + power2(cb))*
     (-1 + power2(snt)) - power2(At)*power2(cb)*(-1 + power2(snt)) + mt2*power2
     (cb)*power2(snt)))/msb22)*power3(DeltaInv(mst22,msb22,mC2))) + Fin3(mst22,
     mst22,mh2,Q2)*(((-4*invdmst*mst24*mt2*(ca*cb*mu2*sa - At*mu*(ca*sa*sb + cb
     *(-1 + power2(sa))) + sb*power2(At)*(-1 + power2(sa))))/sb + 2*mst24*(mt2*
     (-1 + power2(sa)) - (-2*At*ca*mu*sa + power2(At)*(-1 + power2(sa)) - mu2*
     power2(sa))*(-1 + power2(snt))*power2(snt)))*power2(DeltaInv(mst22,mst22,
     mh2)) + ((-16*invdmst*mh2*mst26*mt2*(ca*cb*mu2*sa - At*mu*(ca*sa*sb + cb*(
     -1 + power2(sa))) + sb*power2(At)*(-1 + power2(sa))))/sb + 8*mh2*mst26*(
     mt2*(-1 + power2(sa)) - (-2*At*ca*mu*sa + power2(At)*(-1 + power2(sa)) -
     mu2*power2(sa))*(-1 + power2(snt))*power2(snt)))*power3(DeltaInv(mst22,
     mst22,mh2))) + Fin3(mst22,mw2,msb12,Q2)*(((4*invdmst*mst22*mt2*power2(mw2)
     *(sb*power2(At)*(-1 + power2(cb)) - mu2*sb*power2(cb) + At*cb*mu*(1 -
     power2(cb) + power2(sb)))*(-1 + power2(snb)))/(msb12*sb) + (2*mst22*power2
     (mw2)*(-1 + power2(snb))*(2*At*cb*mu*sb*(-1 + power2(snt)) + power2(At)*(-
     1 + power2(cb))*(-1 + power2(snt)) + mt2*power2(snt) + power2(cb)*(mu2 -
     mt2*power2(snt) - mu2*power2(snt))))/msb12)*power2(DeltaInv(mst22,mw2,
     msb12)) + ((-4*invdmst*mt2*(mst26 + mw2*(-2*mst24 + mst22*mw2) - 2*msb12*(
     mst24 + mst22*mw2))*power2(mw2)*(sb*power2(At)*(-1 + power2(cb)) - mu2*sb*
     power2(cb) + At*cb*mu*(1 - power2(cb) + power2(sb)))*(-1 + power2(snb)))/(
     msb12*sb) - (2*(mst26 + mw2*(-2*mst24 + mst22*mw2) - 2*msb12*(mst24 +
     mst22*mw2))*power2(mw2)*(-1 + power2(snb))*(2*At*cb*mu*sb*(-1 + power2(snt
     )) + power2(At)*(-1 + power2(cb))*(-1 + power2(snt)) + mt2*power2(snt) +
     power2(cb)*(mu2 - mt2*power2(snt) - mu2*power2(snt))))/msb12)*power3(
     DeltaInv(mst22,mw2,msb12))) + Fin3(mst22,mw2,msb22,Q2)*(((-4*invdmst*mst22
     *mt2*power2(mw2)*(sb*power2(At)*(-1 + power2(cb)) - mu2*sb*power2(cb) + At
     *cb*mu*(1 - power2(cb) + power2(sb)))*power2(snb))/(msb22*sb) - (2*mst22*
     power2(mw2)*power2(snb)*(2*At*cb*mu*sb*(-1 + power2(snt)) + power2(At)*(-1
      + power2(cb))*(-1 + power2(snt)) + mt2*power2(snt) + power2(cb)*(mu2 -
     mt2*power2(snt) - mu2*power2(snt))))/msb22)*power2(DeltaInv(mst22,mw2,
     msb22)) + ((4*invdmst*mt2*(mst26 + mw2*(-2*mst24 + mst22*mw2) - 2*msb22*(
     mst24 + mst22*mw2))*power2(mw2)*(sb*power2(At)*(-1 + power2(cb)) - mu2*sb*
     power2(cb) + At*cb*mu*(1 - power2(cb) + power2(sb)))*power2(snb))/(msb22*
     sb) + (2*(mst26 + mw2*(-2*mst24 + mst22*mw2) - 2*msb22*(mst24 + mst22*mw2)
     )*power2(mw2)*power2(snb)*(2*At*cb*mu*sb*(-1 + power2(snt)) + power2(At)*(
     -1 + power2(cb))*(-1 + power2(snt)) + mt2*power2(snt) + power2(cb)*(mu2 -
     mt2*power2(snt) - mu2*power2(snt))))/msb22)*power3(DeltaInv(mst22,mw2,
     msb22))) + Fin3(mst22,mz2,mst12,Q2)*((mst22*(2*At*cb*mu*sb + power2(At)*(-
     1 + power2(cb)) - mu2*power2(cb))*power2(mz2)*power2(DeltaInv(mst22,mz2,
     mst12)))/mst12 - ((mst26 + mz2*(-2*mst24 + mst22*mz2) - 2*mst12*(mst24 +
     mst22*mz2))*(2*At*cb*mu*sb + power2(At)*(-1 + power2(cb)) - mu2*power2(cb)
     )*power2(mz2)*power3(DeltaInv(mst22,mz2,mst12)))/mst12) + Fin3(mw2,msb12,
     mst12,Q2)*(((-4*invdmst*mst12*mt2*power2(mw2)*(sb*power2(At)*(-1 + power2(
     cb)) - mu2*sb*power2(cb) + At*cb*mu*(1 - power2(cb) + power2(sb)))*(-1 +
     power2(snb)))/(msb12*sb) + (2*mst12*power2(mw2)*(-1 + power2(snb))*(mt2*(-
     1 + power2(cb))*(-1 + power2(snt)) + (-2*At*cb*mu*sb - power2(At)*(-1 +
     power2(cb)) + mu2*power2(cb))*power2(snt)))/msb12)*power2(DeltaInv(mw2,
     msb12,mst12)) + ((4*invdmst*mt2*(mst16 + mw2*(-2*mst14 + mst12*mw2) - 2*
     msb12*(mst14 + mst12*mw2))*power2(mw2)*(sb*power2(At)*(-1 + power2(cb)) -
     mu2*sb*power2(cb) + At*cb*mu*(1 - power2(cb) + power2(sb)))*(-1 + power2(
     snb)))/(msb12*sb) - (2*(mst16 + mw2*(-2*mst14 + mst12*mw2) - 2*msb12*(
     mst14 + mst12*mw2))*power2(mw2)*(-1 + power2(snb))*(mt2*(-1 + power2(cb))*
     (-1 + power2(snt)) + (-2*At*cb*mu*sb - power2(At)*(-1 + power2(cb)) + mu2*
     power2(cb))*power2(snt)))/msb12)*power3(DeltaInv(mw2,msb12,mst12))) + Fin3
     (mw2,msb22,mst12,Q2)*(((4*invdmst*mst12*mt2*power2(mw2)*(sb*power2(At)*(-1
      + power2(cb)) - mu2*sb*power2(cb) + At*cb*mu*(1 - power2(cb) + power2(sb)
     ))*power2(snb))/(msb22*sb) - (2*mst12*power2(mw2)*power2(snb)*(mt2*(-1 +
     power2(cb))*(-1 + power2(snt)) + (-2*At*cb*mu*sb - power2(At)*(-1 + power2
     (cb)) + mu2*power2(cb))*power2(snt)))/msb22)*power2(DeltaInv(mw2,msb22,
     mst12)) + ((-4*invdmst*mt2*(mst16 + mw2*(-2*mst14 + mst12*mw2) - 2*msb22*(
     mst14 + mst12*mw2))*power2(mw2)*(sb*power2(At)*(-1 + power2(cb)) - mu2*sb*
     power2(cb) + At*cb*mu*(1 - power2(cb) + power2(sb)))*power2(snb))/(msb22*
     sb) + (2*(mst16 + mw2*(-2*mst14 + mst12*mw2) - 2*msb22*(mst14 + mst12*mw2)
     )*power2(mw2)*power2(snb)*(mt2*(-1 + power2(cb))*(-1 + power2(snt)) + (-2*
     At*cb*mu*sb - power2(At)*(-1 + power2(cb)) + mu2*power2(cb))*power2(snt)))
     /msb22)*power3(DeltaInv(mw2,msb22,mst12))) + Fin3(msb12,mt2,mu2,Q2)*((-2*(
     msb14 + msb12*(mt2 + 3*mu2))*power2(mu2)*(-1 + power2(snb))*power2(
     DeltaInv(msb12,mt2,mu2)))/mt2 + (2*power2(mu2)*(-1 + power2(snb))*power3(
     DeltaInv(msb12,mt2,mu2))*(msb16*(-mt2 + mu2) - 5*msb14*mu2*(2*mt2 + mu2) +
     msb12*(-5*mt2 + 3*mu2)*power2(mu2) + power4(msb12)))/mt2) + Fin3(mt2,mu2,
     msb22,Q2)*((2*(msb24 + msb22*(mt2 + 3*mu2))*power2(mu2)*power2(snb)*power2
     (DeltaInv(mt2,mu2,msb22)))/mt2 - (2*power2(mu2)*power2(snb)*power3(
     DeltaInv(mt2,mu2,msb22))*(msb26*(-mt2 + mu2) - 5*msb24*mu2*(2*mt2 + mu2) +
     msb22*(-5*mt2 + 3*mu2)*power2(mu2) + power4(msb22)))/mt2) + power3(
     DeltaInv(mH2,mst12,mst12))*((8*invdmst*mH2*mt2*sa*(At*sa*(cb*mu - At*sb) +
     ca*(-(cb*mu2) + At*mu*sb))*(42 + 8*lmH2*(-3 + lmst12) - 12*lmst12 + 4*
     power2(lmH2) + power2(Pi))*power4(mst12))/sb - 4*mH2*(42 + 8*lmH2*(-3 +
     lmst12) - 12*lmst12 + 4*power2(lmH2) + power2(Pi))*(mt2*power2(sa) + (At*
     sa*(2*ca*mu - At*sa) + mu2*(-1 + power2(sa)))*(-1 + power2(snt))*power2(
     snt))*power4(mst12)) + power3(DeltaInv(mh2,mst12,mst12))*((8*invdmst*mh2*
     mt2*(42 + 8*lmh2*(-3 + lmst12) - 12*lmst12 + 4*power2(lmh2) + power2(Pi))*
     (ca*cb*mu2*sa - At*mu*(ca*sa*sb + cb*(-1 + power2(sa))) + sb*power2(At)*(-
     1 + power2(sa)))*power4(mst12))/sb + 4*mh2*(42 + 8*lmh2*(-3 + lmst12) - 12
     *lmst12 + 4*power2(lmh2) + power2(Pi))*(mt2*(-1 + power2(sa)) - (-2*At*ca*
     mu*sa + power2(At)*(-1 + power2(sa)) - mu2*power2(sa))*(-1 + power2(snt))*
     power2(snt))*power4(mst12)) + Fin3(mt2,mu2,mst12,Q2)*((2*(mst14 + mst12*(
     mt2 + 3*mu2))*power2(mu2)*power2(DeltaInv(mt2,mu2,mst12)))/mt2 - (2*power2
     (mu2)*power3(DeltaInv(mt2,mu2,mst12))*(mst16*(-mt2 + mu2) - 5*mst14*mu2*(2
     *mt2 + mu2) + mst12*(-5*mt2 + 3*mu2)*power2(mu2) + power4(mst12)))/mt2) +
     power3(DeltaInv(msb12,mst12,mC2))*((cb*invdmst*mt2*power2(mC2)*(cb*mu2*sb
     - cb*sb*power2(At) + At*mu*(power2(cb) - power2(sb)))*(-1 + power2(snb))*(
     -(msb12*mst16*(-6*(12 + lmC2)*lmst12 + 6*lmsb12*(-18 + lmC2 + 5*lmst12) +
     18*power2(lmsb12) + 12*power2(lmst12) + 5*(42 + power2(Pi)))) - mC2*(3*
     mst16*(42 - 24*lmst12 - 2*lmC2*(-6 + lmsb12 + lmst12) + 2*lmsb12*(-12 + 5*
     lmst12) - 2*power2(lmC2) + 4*power2(lmsb12) + 4*power2(lmst12) + power2(Pi
     )) + 2*msb12*mst14*(294 + 30*lmsb12*(-6 + lmst12) - 36*lmst12 + 6*lmC2*(5*
     lmsb12 - 3*(2 + lmst12)) + 6*power2(lmC2) + 30*power2(lmsb12) + 6*power2(
     lmst12) + 7*power2(Pi))) - power2(mC2)*(3*mst14*(42 + 2*lmC2*(-12 + 5*
     lmsb12 - lmst12) + 12*lmst12 - 2*lmsb12*(12 + lmst12) + 4*power2(lmC2) + 4
     *power2(lmsb12) - 2*power2(lmst12) + power2(Pi)) + msb12*mst12*(6*lmC2*(-
     12 + 5*lmsb12 - lmst12) + 6*lmsb12*(-18 + lmst12) + 12*power2(lmC2) + 18*
     power2(lmsb12) + 5*(42 + power2(Pi)))) + 3*mst12*(42 + 2*lmC2*(-6 + 3*
     lmsb12 - lmst12) + 2*lmsb12*(-12 + lmst12) + 2*power2(lmC2) + 4*power2(
     lmsb12) + power2(Pi))*power3(mC2) + 3*(42 + 2*lmsb12*(lmC2 + 3*(-4 +
     lmst12)) - 2*(6 + lmC2)*lmst12 + 4*power2(lmsb12) + 2*power2(lmst12) +
     power2(Pi))*power4(mst12)))/(3.*msb12*sb) + (power2(mC2)*(-1 + power2(snb)
     )*(-(mu2*power2(snt)) - 2*At*cb*mu*sb*power2(snt) + power2(cb)*(mt2*(-1 +
     power2(snt)) + (mu2 - power2(At))*power2(snt)))*(-(msb12*mst16*(-6*(12 +
     lmC2)*lmst12 + 6*lmsb12*(-18 + lmC2 + 5*lmst12) + 18*power2(lmsb12) + 12*
     power2(lmst12) + 5*(42 + power2(Pi)))) - mC2*(3*mst16*(42 - 24*lmst12 - 2*
     lmC2*(-6 + lmsb12 + lmst12) + 2*lmsb12*(-12 + 5*lmst12) - 2*power2(lmC2) +
     4*power2(lmsb12) + 4*power2(lmst12) + power2(Pi)) + 2*msb12*mst14*(294 +
     30*lmsb12*(-6 + lmst12) - 36*lmst12 + 6*lmC2*(5*lmsb12 - 3*(2 + lmst12)) +
     6*power2(lmC2) + 30*power2(lmsb12) + 6*power2(lmst12) + 7*power2(Pi))) -
     power2(mC2)*(3*mst14*(42 + 2*lmC2*(-12 + 5*lmsb12 - lmst12) + 12*lmst12 -
     2*lmsb12*(12 + lmst12) + 4*power2(lmC2) + 4*power2(lmsb12) - 2*power2(
     lmst12) + power2(Pi)) + msb12*mst12*(6*lmC2*(-12 + 5*lmsb12 - lmst12) + 6*
     lmsb12*(-18 + lmst12) + 12*power2(lmC2) + 18*power2(lmsb12) + 5*(42 +
     power2(Pi)))) + 3*mst12*(42 + 2*lmC2*(-6 + 3*lmsb12 - lmst12) + 2*lmsb12*(
     -12 + lmst12) + 2*power2(lmC2) + 4*power2(lmsb12) + power2(Pi))*power3(mC2
     ) + 3*(42 + 2*lmsb12*(lmC2 + 3*(-4 + lmst12)) - 2*(6 + lmC2)*lmst12 + 4*
     power2(lmsb12) + 2*power2(lmst12) + power2(Pi))*power4(mst12)))/(6.*msb12)
     ) + power3(DeltaInv(mw2,msb12,mst12))*(-(invdmst*mt2*power2(mw2)*(sb*
     power2(At)*(-1 + power2(cb)) - mu2*sb*power2(cb) + At*cb*mu*(1 - power2(cb
     ) + power2(sb)))*(-1 + power2(snb))*(3*mw2*(mst16*(42 + 2*lmsb12*(-12 + 5*
     lmst12 - lmw2) + 12*lmw2 - 2*lmst12*(12 + lmw2) + 4*power2(lmsb12) + 4*
     power2(lmst12) - 2*power2(lmw2) + power2(Pi)) + mst14*mw2*(42 - 2*lmsb12*(
     12 + lmst12 - 5*lmw2) - 2*lmst12*(-6 + lmw2) - 24*lmw2 + 4*power2(lmsb12)
     - 2*power2(lmst12) + 4*power2(lmw2) + power2(Pi))) + msb12*(2*mst14*mw2*(
     294 - 36*lmw2 - 18*lmst12*(2 + lmw2) + 30*lmsb12*(-6 + lmst12 + lmw2) + 30
     *power2(lmsb12) + 6*power2(lmst12) + 6*power2(lmw2) + 7*power2(Pi)) +
     mst16*(-6*lmst12*(12 + lmw2) + 6*lmsb12*(-18 + 5*lmst12 + lmw2) + 18*
     power2(lmsb12) + 12*power2(lmst12) + 5*(42 + power2(Pi)))) + mst12*power2(
     mw2)*(-3*mw2*(42 + 2*lmsb12*(lmst12 + 3*(-4 + lmw2)) - 2*(6 + lmst12)*lmw2
      + 4*power2(lmsb12) + 2*power2(lmw2) + power2(Pi)) + msb12*(-6*(12 +
     lmst12)*lmw2 + 6*lmsb12*(-18 + lmst12 + 5*lmw2) + 18*power2(lmsb12) + 12*
     power2(lmw2) + 5*(42 + power2(Pi)))) - 3*(42 - 2*lmst12*(6 + lmw2) + 2*
     lmsb12*(-12 + 3*lmst12 + lmw2) + 4*power2(lmsb12) + 2*power2(lmst12) +
     power2(Pi))*power4(mst12)))/(3.*msb12*sb) - (power2(mw2)*(-1 + power2(snb)
     )*(-(mt2*(-1 + power2(cb))*(-1 + power2(snt))) + (2*At*cb*mu*sb + power2(
     At)*(-1 + power2(cb)) - mu2*power2(cb))*power2(snt))*(3*mw2*(mst16*(42 + 2
     *lmsb12*(-12 + 5*lmst12 - lmw2) + 12*lmw2 - 2*lmst12*(12 + lmw2) + 4*
     power2(lmsb12) + 4*power2(lmst12) - 2*power2(lmw2) + power2(Pi)) + mst14*
     mw2*(42 - 2*lmsb12*(12 + lmst12 - 5*lmw2) - 2*lmst12*(-6 + lmw2) - 24*lmw2
      + 4*power2(lmsb12) - 2*power2(lmst12) + 4*power2(lmw2) + power2(Pi))) +
     msb12*(2*mst14*mw2*(294 - 36*lmw2 - 18*lmst12*(2 + lmw2) + 30*lmsb12*(-6 +
     lmst12 + lmw2) + 30*power2(lmsb12) + 6*power2(lmst12) + 6*power2(lmw2) + 7
     *power2(Pi)) + mst16*(-6*lmst12*(12 + lmw2) + 6*lmsb12*(-18 + 5*lmst12 +
     lmw2) + 18*power2(lmsb12) + 12*power2(lmst12) + 5*(42 + power2(Pi)))) +
     mst12*power2(mw2)*(-3*mw2*(42 + 2*lmsb12*(lmst12 + 3*(-4 + lmw2)) - 2*(6 +
     lmst12)*lmw2 + 4*power2(lmsb12) + 2*power2(lmw2) + power2(Pi)) + msb12*(-6
     *(12 + lmst12)*lmw2 + 6*lmsb12*(-18 + lmst12 + 5*lmw2) + 18*power2(lmsb12)
     + 12*power2(lmw2) + 5*(42 + power2(Pi)))) - 3*(42 - 2*lmst12*(6 + lmw2) +
     2*lmsb12*(-12 + 3*lmst12 + lmw2) + 4*power2(lmsb12) + 2*power2(lmst12) +
     power2(Pi))*power4(mst12)))/(6.*msb12)) + power3(DeltaInv(msb22,mst12,mC2)
     )*((cb*invdmst*mt2*power2(mC2)*(-(cb*mu2*sb) + cb*sb*power2(At) + At*mu*(-
     power2(cb) + power2(sb)))*power2(snb)*(-(msb22*mst16*(-6*(12 + lmC2)*
     lmst12 + 6*lmsb22*(-18 + lmC2 + 5*lmst12) + 18*power2(lmsb22) + 12*power2(
     lmst12) + 5*(42 + power2(Pi)))) - mC2*(3*mst16*(42 - 24*lmst12 - 2*lmC2*(-
     6 + lmsb22 + lmst12) + 2*lmsb22*(-12 + 5*lmst12) - 2*power2(lmC2) + 4*
     power2(lmsb22) + 4*power2(lmst12) + power2(Pi)) + 2*msb22*mst14*(294 + 30*
     lmsb22*(-6 + lmst12) - 36*lmst12 + 6*lmC2*(5*lmsb22 - 3*(2 + lmst12)) + 6*
     power2(lmC2) + 30*power2(lmsb22) + 6*power2(lmst12) + 7*power2(Pi))) -
     power2(mC2)*(3*mst14*(42 + 2*lmC2*(-12 + 5*lmsb22 - lmst12) + 12*lmst12 -
     2*lmsb22*(12 + lmst12) + 4*power2(lmC2) + 4*power2(lmsb22) - 2*power2(
     lmst12) + power2(Pi)) + msb22*mst12*(6*lmC2*(-12 + 5*lmsb22 - lmst12) + 6*
     lmsb22*(-18 + lmst12) + 12*power2(lmC2) + 18*power2(lmsb22) + 5*(42 +
     power2(Pi)))) + 3*mst12*(42 + 2*lmC2*(-6 + 3*lmsb22 - lmst12) + 2*lmsb22*(
     -12 + lmst12) + 2*power2(lmC2) + 4*power2(lmsb22) + power2(Pi))*power3(mC2
     ) + 3*(42 + 2*lmsb22*(lmC2 + 3*(-4 + lmst12)) - 2*(6 + lmC2)*lmst12 + 4*
     power2(lmsb22) + 2*power2(lmst12) + power2(Pi))*power4(mst12)))/(3.*msb22*
     sb) - (power2(mC2)*power2(snb)*(-(mu2*power2(snt)) - 2*At*cb*mu*sb*power2(
     snt) + power2(cb)*(mt2*(-1 + power2(snt)) + (mu2 - power2(At))*power2(snt)
     ))*(-(msb22*mst16*(-6*(12 + lmC2)*lmst12 + 6*lmsb22*(-18 + lmC2 + 5*lmst12
     ) + 18*power2(lmsb22) + 12*power2(lmst12) + 5*(42 + power2(Pi)))) - mC2*(3
     *mst16*(42 - 24*lmst12 - 2*lmC2*(-6 + lmsb22 + lmst12) + 2*lmsb22*(-12 + 5
     *lmst12) - 2*power2(lmC2) + 4*power2(lmsb22) + 4*power2(lmst12) + power2(
     Pi)) + 2*msb22*mst14*(294 + 30*lmsb22*(-6 + lmst12) - 36*lmst12 + 6*lmC2*(
     5*lmsb22 - 3*(2 + lmst12)) + 6*power2(lmC2) + 30*power2(lmsb22) + 6*power2
     (lmst12) + 7*power2(Pi))) - power2(mC2)*(3*mst14*(42 + 2*lmC2*(-12 + 5*
     lmsb22 - lmst12) + 12*lmst12 - 2*lmsb22*(12 + lmst12) + 4*power2(lmC2) + 4
     *power2(lmsb22) - 2*power2(lmst12) + power2(Pi)) + msb22*mst12*(6*lmC2*(-
     12 + 5*lmsb22 - lmst12) + 6*lmsb22*(-18 + lmst12) + 12*power2(lmC2) + 18*
     power2(lmsb22) + 5*(42 + power2(Pi)))) + 3*mst12*(42 + 2*lmC2*(-6 + 3*
     lmsb22 - lmst12) + 2*lmsb22*(-12 + lmst12) + 2*power2(lmC2) + 4*power2(
     lmsb22) + power2(Pi))*power3(mC2) + 3*(42 + 2*lmsb22*(lmC2 + 3*(-4 +
     lmst12)) - 2*(6 + lmC2)*lmst12 + 4*power2(lmsb22) + 2*power2(lmst12) +
     power2(Pi))*power4(mst12)))/(6.*msb22)) + power3(DeltaInv(mw2,msb22,mst12)
     )*((invdmst*mt2*power2(mw2)*(sb*power2(At)*(-1 + power2(cb)) - mu2*sb*
     power2(cb) + At*cb*mu*(1 - power2(cb) + power2(sb)))*power2(snb)*(3*mw2*(
     mst16*(42 + 2*lmsb22*(-12 + 5*lmst12 - lmw2) + 12*lmw2 - 2*lmst12*(12 +
     lmw2) + 4*power2(lmsb22) + 4*power2(lmst12) - 2*power2(lmw2) + power2(Pi))
     + mst14*mw2*(42 - 2*lmsb22*(12 + lmst12 - 5*lmw2) - 2*lmst12*(-6 + lmw2) -
     24*lmw2 + 4*power2(lmsb22) - 2*power2(lmst12) + 4*power2(lmw2) + power2(Pi
     ))) + msb22*(2*mst14*mw2*(294 - 36*lmw2 - 18*lmst12*(2 + lmw2) + 30*lmsb22
     *(-6 + lmst12 + lmw2) + 30*power2(lmsb22) + 6*power2(lmst12) + 6*power2(
     lmw2) + 7*power2(Pi)) + mst16*(-6*lmst12*(12 + lmw2) + 6*lmsb22*(-18 + 5*
     lmst12 + lmw2) + 18*power2(lmsb22) + 12*power2(lmst12) + 5*(42 + power2(Pi
     )))) + mst12*power2(mw2)*(-3*mw2*(42 + 2*lmsb22*(lmst12 + 3*(-4 + lmw2)) -
     2*(6 + lmst12)*lmw2 + 4*power2(lmsb22) + 2*power2(lmw2) + power2(Pi)) +
     msb22*(-6*(12 + lmst12)*lmw2 + 6*lmsb22*(-18 + lmst12 + 5*lmw2) + 18*
     power2(lmsb22) + 12*power2(lmw2) + 5*(42 + power2(Pi)))) - 3*(42 - 2*
     lmst12*(6 + lmw2) + 2*lmsb22*(-12 + 3*lmst12 + lmw2) + 4*power2(lmsb22) +
     2*power2(lmst12) + power2(Pi))*power4(mst12)))/(3.*msb22*sb) + (power2(mw2
     )*power2(snb)*(-(mt2*(-1 + power2(cb))*(-1 + power2(snt))) + (2*At*cb*mu*
     sb + power2(At)*(-1 + power2(cb)) - mu2*power2(cb))*power2(snt))*(3*mw2*(
     mst16*(42 + 2*lmsb22*(-12 + 5*lmst12 - lmw2) + 12*lmw2 - 2*lmst12*(12 +
     lmw2) + 4*power2(lmsb22) + 4*power2(lmst12) - 2*power2(lmw2) + power2(Pi))
     + mst14*mw2*(42 - 2*lmsb22*(12 + lmst12 - 5*lmw2) - 2*lmst12*(-6 + lmw2) -
     24*lmw2 + 4*power2(lmsb22) - 2*power2(lmst12) + 4*power2(lmw2) + power2(Pi
     ))) + msb22*(2*mst14*mw2*(294 - 36*lmw2 - 18*lmst12*(2 + lmw2) + 30*lmsb22
     *(-6 + lmst12 + lmw2) + 30*power2(lmsb22) + 6*power2(lmst12) + 6*power2(
     lmw2) + 7*power2(Pi)) + mst16*(-6*lmst12*(12 + lmw2) + 6*lmsb22*(-18 + 5*
     lmst12 + lmw2) + 18*power2(lmsb22) + 12*power2(lmst12) + 5*(42 + power2(Pi
     )))) + mst12*power2(mw2)*(-3*mw2*(42 + 2*lmsb22*(lmst12 + 3*(-4 + lmw2)) -
     2*(6 + lmst12)*lmw2 + 4*power2(lmsb22) + 2*power2(lmw2) + power2(Pi)) +
     msb22*(-6*(12 + lmst12)*lmw2 + 6*lmsb22*(-18 + lmst12 + 5*lmw2) + 18*
     power2(lmsb22) + 12*power2(lmw2) + 5*(42 + power2(Pi)))) - 3*(42 - 2*
     lmst12*(6 + lmw2) + 2*lmsb22*(-12 + 3*lmst12 + lmw2) + 4*power2(lmsb22) +
     2*power2(lmst12) + power2(Pi))*power4(mst12)))/(6.*msb22)) + ((mu2 + 2*At*
     cb*mu*sb - mu2*power2(cb) + power2(At)*power2(cb))*power2(mA2)*power3(
     DeltaInv(mst22,mA2,mst12))*(-(mst12*mst26*(-6*(12 + lmA2)*lmst22 + 6*
     lmst12*(-18 + lmA2 + 5*lmst22) + 18*power2(lmst12) + 12*power2(lmst22) + 5
     *(42 + power2(Pi)))) - mA2*(3*mst26*(42 - 24*lmst22 - 2*lmA2*(-6 + lmst12
     + lmst22) + 2*lmst12*(-12 + 5*lmst22) - 2*power2(lmA2) + 4*power2(lmst12)
     + 4*power2(lmst22) + power2(Pi)) + 2*mst12*mst24*(294 + 30*lmst12*(-6 +
     lmst22) - 36*lmst22 + 6*lmA2*(5*lmst12 - 3*(2 + lmst22)) + 6*power2(lmA2)
     + 30*power2(lmst12) + 6*power2(lmst22) + 7*power2(Pi))) - power2(mA2)*(3*
     mst24*(42 + 2*lmA2*(-12 + 5*lmst12 - lmst22) + 12*lmst22 - 2*lmst12*(12 +
     lmst22) + 4*power2(lmA2) + 4*power2(lmst12) - 2*power2(lmst22) + power2(Pi
     )) + mst12*mst22*(6*lmA2*(-12 + 5*lmst12 - lmst22) + 6*lmst12*(-18 +
     lmst22) + 12*power2(lmA2) + 18*power2(lmst12) + 5*(42 + power2(Pi)))) + 3*
     mst22*(42 + 2*lmA2*(-6 + 3*lmst12 - lmst22) + 2*lmst12*(-12 + lmst22) + 2*
     power2(lmA2) + 4*power2(lmst12) + power2(Pi))*power3(mA2) + 3*(42 + 2*
     lmst12*(lmA2 + 3*(-4 + lmst22)) - 2*(6 + lmA2)*lmst22 + 4*power2(lmst12) +
     2*power2(lmst22) + power2(Pi))*power4(mst22)))/(12.*mst12) + ((-2*At*ca*mu
     *sa + power2(At)*(-1 + power2(sa)) - mu2*power2(sa))*power2(mh2 - 2*mh2*
     power2(snt))*power3(DeltaInv(mst22,mh2,mst12))*(mst12*mst26*(-6*(12 + lmh2
     )*lmst22 + 6*lmst12*(-18 + lmh2 + 5*lmst22) + 18*power2(lmst12) + 12*
     power2(lmst22) + 5*(42 + power2(Pi))) + mh2*(3*mst26*(42 - 24*lmst22 - 2*
     lmh2*(-6 + lmst12 + lmst22) + 2*lmst12*(-12 + 5*lmst22) - 2*power2(lmh2) +
     4*power2(lmst12) + 4*power2(lmst22) + power2(Pi)) + 2*mst12*mst24*(294 +
     30*lmst12*(-6 + lmst22) - 36*lmst22 + 6*lmh2*(5*lmst12 - 3*(2 + lmst22)) +
     6*power2(lmh2) + 30*power2(lmst12) + 6*power2(lmst22) + 7*power2(Pi))) +
     power2(mh2)*(3*mst24*(42 + 2*lmh2*(-12 + 5*lmst12 - lmst22) + 12*lmst22 -
     2*lmst12*(12 + lmst22) + 4*power2(lmh2) + 4*power2(lmst12) - 2*power2(
     lmst22) + power2(Pi)) + mst12*mst22*(6*lmh2*(-12 + 5*lmst12 - lmst22) + 6*
     lmst12*(-18 + lmst22) + 12*power2(lmh2) + 18*power2(lmst12) + 5*(42 +
     power2(Pi)))) - 3*mst22*(42 + 2*lmh2*(-6 + 3*lmst12 - lmst22) + 2*lmst12*(
     -12 + lmst22) + 2*power2(lmh2) + 4*power2(lmst12) + power2(Pi))*power3(mh2
     ) - 3*(42 + 2*lmst12*(lmh2 + 3*(-4 + lmst22)) - 2*(6 + lmh2)*lmst22 + 4*
     power2(lmst12) + 2*power2(lmst22) + power2(Pi))*power4(mst22)))/(12.*mst12
     ) - ((mu2 + At*sa*(-2*ca*mu + At*sa) - mu2*power2(sa))*power2(mH2 - 2*mH2*
     power2(snt))*power3(DeltaInv(mH2,mst22,mst12))*(mst12*mst26*(-6*(12 + lmH2
     )*lmst22 + 6*lmst12*(-18 + lmH2 + 5*lmst22) + 18*power2(lmst12) + 12*
     power2(lmst22) + 5*(42 + power2(Pi))) + mH2*(3*mst26*(42 - 24*lmst22 - 2*
     lmH2*(-6 + lmst12 + lmst22) + 2*lmst12*(-12 + 5*lmst22) - 2*power2(lmH2) +
     4*power2(lmst12) + 4*power2(lmst22) + power2(Pi)) + 2*mst12*mst24*(294 +
     30*lmst12*(-6 + lmst22) - 36*lmst22 + 6*lmH2*(5*lmst12 - 3*(2 + lmst22)) +
     6*power2(lmH2) + 30*power2(lmst12) + 6*power2(lmst22) + 7*power2(Pi))) +
     power2(mH2)*(3*mst24*(42 + 2*lmH2*(-12 + 5*lmst12 - lmst22) + 12*lmst22 -
     2*lmst12*(12 + lmst22) + 4*power2(lmH2) + 4*power2(lmst12) - 2*power2(
     lmst22) + power2(Pi)) + mst12*mst22*(6*lmH2*(-12 + 5*lmst12 - lmst22) + 6*
     lmst12*(-18 + lmst22) + 12*power2(lmH2) + 18*power2(lmst12) + 5*(42 +
     power2(Pi)))) - 3*mst22*(42 + 2*lmH2*(-6 + 3*lmst12 - lmst22) + 2*lmst12*(
     -12 + lmst22) + 2*power2(lmH2) + 4*power2(lmst12) + power2(Pi))*power3(mH2
     ) - 3*(42 + 2*lmst12*(lmH2 + 3*(-4 + lmst22)) - 2*(6 + lmH2)*lmst22 + 4*
     power2(lmst12) + 2*power2(lmst22) + power2(Pi))*power4(mst22)))/(12.*mst12
     ) + ((2*At*cb*mu*sb + power2(At)*(-1 + power2(cb)) - mu2*power2(cb))*
     power2(mz2)*power3(DeltaInv(mst22,mz2,mst12))*(3*mz2*(mst26*(42 + 2*lmst12
     *(-12 + 5*lmst22 - lmz2) + 12*lmz2 - 2*lmst22*(12 + lmz2) + 4*power2(
     lmst12) + 4*power2(lmst22) - 2*power2(lmz2) + power2(Pi)) + mst24*mz2*(42
     - 2*lmst12*(12 + lmst22 - 5*lmz2) - 2*lmst22*(-6 + lmz2) - 24*lmz2 + 4*
     power2(lmst12) - 2*power2(lmst22) + 4*power2(lmz2) + power2(Pi))) + mst12*
     (2*mst24*mz2*(294 - 36*lmz2 - 18*lmst22*(2 + lmz2) + 30*lmst12*(-6 +
     lmst22 + lmz2) + 30*power2(lmst12) + 6*power2(lmst22) + 6*power2(lmz2) + 7
     *power2(Pi)) + mst26*(-6*lmst22*(12 + lmz2) + 6*lmst12*(-18 + 5*lmst22 +
     lmz2) + 18*power2(lmst12) + 12*power2(lmst22) + 5*(42 + power2(Pi)))) +
     mst22*power2(mz2)*(-3*mz2*(42 + 2*lmst12*(lmst22 + 3*(-4 + lmz2)) - 2*(6 +
     lmst22)*lmz2 + 4*power2(lmst12) + 2*power2(lmz2) + power2(Pi)) + mst12*(-6
     *(12 + lmst22)*lmz2 + 6*lmst12*(-18 + lmst22 + 5*lmz2) + 18*power2(lmst12)
     + 12*power2(lmz2) + 5*(42 + power2(Pi)))) - 3*(42 - 2*lmst22*(6 + lmz2) +
     2*lmst12*(-12 + 3*lmst22 + lmz2) + 4*power2(lmst12) + 2*power2(lmst22) +
     power2(Pi))*power4(mst22)))/(12.*mst12) + power3(DeltaInv(mH2,mst22,mst22)
     )*((8*invdmst*mH2*mt2*sa*(At*sa*(-(cb*mu) + At*sb) + ca*(cb*mu2 - At*mu*sb
     ))*(42 + 8*lmH2*(-3 + lmst22) - 12*lmst22 + 4*power2(lmH2) + power2(Pi))*
     power4(mst22))/sb - 4*mH2*(42 + 8*lmH2*(-3 + lmst22) - 12*lmst22 + 4*
     power2(lmH2) + power2(Pi))*(mt2*power2(sa) + (At*sa*(2*ca*mu - At*sa) +
     mu2*(-1 + power2(sa)))*(-1 + power2(snt))*power2(snt))*power4(mst22)) +
     power3(DeltaInv(mst22,mst22,mh2))*((-8*invdmst*mh2*mt2*(42 + 8*lmh2*(-3 +
     lmst22) - 12*lmst22 + 4*power2(lmh2) + power2(Pi))*(ca*cb*mu2*sa - At*mu*(
     ca*sa*sb + cb*(-1 + power2(sa))) + sb*power2(At)*(-1 + power2(sa)))*power4
     (mst22))/sb + 4*mh2*(42 + 8*lmh2*(-3 + lmst22) - 12*lmst22 + 4*power2(lmh2
     ) + power2(Pi))*(mt2*(-1 + power2(sa)) - (-2*At*ca*mu*sa + power2(At)*(-1
     + power2(sa)) - mu2*power2(sa))*(-1 + power2(snt))*power2(snt))*power4(
     mst22)) + Fin3(mst22,mt2,mu2,Q2)*((2*(mst24 + mst22*(mt2 + 3*mu2))*power2(
     mu2)*power2(DeltaInv(mst22,mt2,mu2)))/mt2 - (2*power2(mu2)*power3(DeltaInv
     (mst22,mt2,mu2))*(mst26*(-mt2 + mu2) - 5*mst24*mu2*(2*mt2 + mu2) + mst22*(
     -5*mt2 + 3*mu2)*power2(mu2) + power4(mst22)))/mt2) + power3(DeltaInv(mst22
     ,msb12,mC2))*((cb*invdmst*mt2*power2(mC2)*(-(cb*mu2*sb) + cb*sb*power2(At)
     + At*mu*(-power2(cb) + power2(sb)))*(-1 + power2(snb))*(-(msb12*mst26*(-6*
     (12 + lmC2)*lmst22 + 6*lmsb12*(-18 + lmC2 + 5*lmst22) + 18*power2(lmsb12)
     + 12*power2(lmst22) + 5*(42 + power2(Pi)))) - mC2*(3*mst26*(42 - 24*lmst22
      - 2*lmC2*(-6 + lmsb12 + lmst22) + 2*lmsb12*(-12 + 5*lmst22) - 2*power2(
     lmC2) + 4*power2(lmsb12) + 4*power2(lmst22) + power2(Pi)) + 2*msb12*mst24*
     (294 + 30*lmsb12*(-6 + lmst22) - 36*lmst22 + 6*lmC2*(5*lmsb12 - 3*(2 +
     lmst22)) + 6*power2(lmC2) + 30*power2(lmsb12) + 6*power2(lmst22) + 7*
     power2(Pi))) - power2(mC2)*(3*mst24*(42 + 2*lmC2*(-12 + 5*lmsb12 - lmst22)
     + 12*lmst22 - 2*lmsb12*(12 + lmst22) + 4*power2(lmC2) + 4*power2(lmsb12) -
     2*power2(lmst22) + power2(Pi)) + msb12*mst22*(6*lmC2*(-12 + 5*lmsb12 -
     lmst22) + 6*lmsb12*(-18 + lmst22) + 12*power2(lmC2) + 18*power2(lmsb12) +
     5*(42 + power2(Pi)))) + 3*mst22*(42 + 2*lmC2*(-6 + 3*lmsb12 - lmst22) + 2*
     lmsb12*(-12 + lmst22) + 2*power2(lmC2) + 4*power2(lmsb12) + power2(Pi))*
     power3(mC2) + 3*(42 + 2*lmsb12*(lmC2 + 3*(-4 + lmst22)) - 2*(6 + lmC2)*
     lmst22 + 4*power2(lmsb12) + 2*power2(lmst22) + power2(Pi))*power4(mst22)))
     /(3.*msb12*sb) + (power2(mC2)*(-1 + power2(snb))*(2*At*cb*mu*sb*(-1 +
     power2(snt)) - mu2*(-1 + power2(cb))*(-1 + power2(snt)) + power2(At)*
     power2(cb)*(-1 + power2(snt)) - mt2*power2(cb)*power2(snt))*(-(msb12*mst26
     *(-6*(12 + lmC2)*lmst22 + 6*lmsb12*(-18 + lmC2 + 5*lmst22) + 18*power2(
     lmsb12) + 12*power2(lmst22) + 5*(42 + power2(Pi)))) - mC2*(3*mst26*(42 -
     24*lmst22 - 2*lmC2*(-6 + lmsb12 + lmst22) + 2*lmsb12*(-12 + 5*lmst22) - 2*
     power2(lmC2) + 4*power2(lmsb12) + 4*power2(lmst22) + power2(Pi)) + 2*msb12
     *mst24*(294 + 30*lmsb12*(-6 + lmst22) - 36*lmst22 + 6*lmC2*(5*lmsb12 - 3*(
     2 + lmst22)) + 6*power2(lmC2) + 30*power2(lmsb12) + 6*power2(lmst22) + 7*
     power2(Pi))) - power2(mC2)*(3*mst24*(42 + 2*lmC2*(-12 + 5*lmsb12 - lmst22)
     + 12*lmst22 - 2*lmsb12*(12 + lmst22) + 4*power2(lmC2) + 4*power2(lmsb12) -
     2*power2(lmst22) + power2(Pi)) + msb12*mst22*(6*lmC2*(-12 + 5*lmsb12 -
     lmst22) + 6*lmsb12*(-18 + lmst22) + 12*power2(lmC2) + 18*power2(lmsb12) +
     5*(42 + power2(Pi)))) + 3*mst22*(42 + 2*lmC2*(-6 + 3*lmsb12 - lmst22) + 2*
     lmsb12*(-12 + lmst22) + 2*power2(lmC2) + 4*power2(lmsb12) + power2(Pi))*
     power3(mC2) + 3*(42 + 2*lmsb12*(lmC2 + 3*(-4 + lmst22)) - 2*(6 + lmC2)*
     lmst22 + 4*power2(lmsb12) + 2*power2(lmst22) + power2(Pi))*power4(mst22)))
     /(6.*msb12)) + power3(DeltaInv(mst22,mw2,msb12))*((invdmst*mt2*power2(mw2)
     *(sb*power2(At)*(-1 + power2(cb)) - mu2*sb*power2(cb) + At*cb*mu*(1 -
     power2(cb) + power2(sb)))*(-1 + power2(snb))*(3*mw2*(mst26*(42 + 2*lmsb12*
     (-12 + 5*lmst22 - lmw2) + 12*lmw2 - 2*lmst22*(12 + lmw2) + 4*power2(lmsb12
     ) + 4*power2(lmst22) - 2*power2(lmw2) + power2(Pi)) + mst24*mw2*(42 - 2*
     lmsb12*(12 + lmst22 - 5*lmw2) - 2*lmst22*(-6 + lmw2) - 24*lmw2 + 4*power2(
     lmsb12) - 2*power2(lmst22) + 4*power2(lmw2) + power2(Pi))) + msb12*(2*
     mst24*mw2*(294 - 36*lmw2 - 18*lmst22*(2 + lmw2) + 30*lmsb12*(-6 + lmst22 +
     lmw2) + 30*power2(lmsb12) + 6*power2(lmst22) + 6*power2(lmw2) + 7*power2(
     Pi)) + mst26*(-6*lmst22*(12 + lmw2) + 6*lmsb12*(-18 + 5*lmst22 + lmw2) +
     18*power2(lmsb12) + 12*power2(lmst22) + 5*(42 + power2(Pi)))) + mst22*
     power2(mw2)*(-3*mw2*(42 + 2*lmsb12*(lmst22 + 3*(-4 + lmw2)) - 2*(6 +
     lmst22)*lmw2 + 4*power2(lmsb12) + 2*power2(lmw2) + power2(Pi)) + msb12*(-6
     *(12 + lmst22)*lmw2 + 6*lmsb12*(-18 + lmst22 + 5*lmw2) + 18*power2(lmsb12)
     + 12*power2(lmw2) + 5*(42 + power2(Pi)))) - 3*(42 - 2*lmst22*(6 + lmw2) +
     2*lmsb12*(-12 + 3*lmst22 + lmw2) + 4*power2(lmsb12) + 2*power2(lmst22) +
     power2(Pi))*power4(mst22)))/(3.*msb12*sb) - (power2(mw2)*(-1 + power2(snb)
     )*(2*At*cb*mu*sb*(-1 + power2(snt)) + power2(At)*(-1 + power2(cb))*(-1 +
     power2(snt)) + mt2*power2(snt) + power2(cb)*(mu2 - mt2*power2(snt) - mu2*
     power2(snt)))*(-3*mw2*(mst26*(42 + 2*lmsb12*(-12 + 5*lmst22 - lmw2) + 12*
     lmw2 - 2*lmst22*(12 + lmw2) + 4*power2(lmsb12) + 4*power2(lmst22) - 2*
     power2(lmw2) + power2(Pi)) + mst24*mw2*(42 - 2*lmsb12*(12 + lmst22 - 5*
     lmw2) - 2*lmst22*(-6 + lmw2) - 24*lmw2 + 4*power2(lmsb12) - 2*power2(
     lmst22) + 4*power2(lmw2) + power2(Pi))) - msb12*(2*mst24*mw2*(294 - 36*
     lmw2 - 18*lmst22*(2 + lmw2) + 30*lmsb12*(-6 + lmst22 + lmw2) + 30*power2(
     lmsb12) + 6*power2(lmst22) + 6*power2(lmw2) + 7*power2(Pi)) + mst26*(-6*
     lmst22*(12 + lmw2) + 6*lmsb12*(-18 + 5*lmst22 + lmw2) + 18*power2(lmsb12)
     + 12*power2(lmst22) + 5*(42 + power2(Pi)))) + mst22*power2(mw2)*(3*mw2*(42
      + 2*lmsb12*(lmst22 + 3*(-4 + lmw2)) - 2*(6 + lmst22)*lmw2 + 4*power2(
     lmsb12) + 2*power2(lmw2) + power2(Pi)) - msb12*(-6*(12 + lmst22)*lmw2 + 6*
     lmsb12*(-18 + lmst22 + 5*lmw2) + 18*power2(lmsb12) + 12*power2(lmw2) + 5*(
     42 + power2(Pi)))) + 3*(42 - 2*lmst22*(6 + lmw2) + 2*lmsb12*(-12 + 3*
     lmst22 + lmw2) + 4*power2(lmsb12) + 2*power2(lmst22) + power2(Pi))*power4(
     mst22)))/(6.*msb12)) + power3(DeltaInv(mst22,msb22,mC2))*((cb*invdmst*mt2*
     power2(mC2)*(cb*mu2*sb - cb*sb*power2(At) + At*mu*(power2(cb) - power2(sb)
     ))*power2(snb)*(-(msb22*mst26*(-6*(12 + lmC2)*lmst22 + 6*lmsb22*(-18 +
     lmC2 + 5*lmst22) + 18*power2(lmsb22) + 12*power2(lmst22) + 5*(42 + power2(
     Pi)))) - mC2*(3*mst26*(42 - 24*lmst22 - 2*lmC2*(-6 + lmsb22 + lmst22) + 2*
     lmsb22*(-12 + 5*lmst22) - 2*power2(lmC2) + 4*power2(lmsb22) + 4*power2(
     lmst22) + power2(Pi)) + 2*msb22*mst24*(294 + 30*lmsb22*(-6 + lmst22) - 36*
     lmst22 + 6*lmC2*(5*lmsb22 - 3*(2 + lmst22)) + 6*power2(lmC2) + 30*power2(
     lmsb22) + 6*power2(lmst22) + 7*power2(Pi))) - power2(mC2)*(3*mst24*(42 + 2
     *lmC2*(-12 + 5*lmsb22 - lmst22) + 12*lmst22 - 2*lmsb22*(12 + lmst22) + 4*
     power2(lmC2) + 4*power2(lmsb22) - 2*power2(lmst22) + power2(Pi)) + msb22*
     mst22*(6*lmC2*(-12 + 5*lmsb22 - lmst22) + 6*lmsb22*(-18 + lmst22) + 12*
     power2(lmC2) + 18*power2(lmsb22) + 5*(42 + power2(Pi)))) + 3*mst22*(42 + 2
     *lmC2*(-6 + 3*lmsb22 - lmst22) + 2*lmsb22*(-12 + lmst22) + 2*power2(lmC2)
     + 4*power2(lmsb22) + power2(Pi))*power3(mC2) + 3*(42 + 2*lmsb22*(lmC2 + 3*
     (-4 + lmst22)) - 2*(6 + lmC2)*lmst22 + 4*power2(lmsb22) + 2*power2(lmst22)
     + power2(Pi))*power4(mst22)))/(3.*msb22*sb) + (power2(mC2)*power2(snb)*(-2
     *At*cb*mu*sb*(-1 + power2(snt)) + mu2*(-1 + power2(cb))*(-1 + power2(snt))
     - power2(At)*power2(cb)*(-1 + power2(snt)) + mt2*power2(cb)*power2(snt))*(
     -(msb22*mst26*(-6*(12 + lmC2)*lmst22 + 6*lmsb22*(-18 + lmC2 + 5*lmst22) +
     18*power2(lmsb22) + 12*power2(lmst22) + 5*(42 + power2(Pi)))) - mC2*(3*
     mst26*(42 - 24*lmst22 - 2*lmC2*(-6 + lmsb22 + lmst22) + 2*lmsb22*(-12 + 5*
     lmst22) - 2*power2(lmC2) + 4*power2(lmsb22) + 4*power2(lmst22) + power2(Pi
     )) + 2*msb22*mst24*(294 + 30*lmsb22*(-6 + lmst22) - 36*lmst22 + 6*lmC2*(5*
     lmsb22 - 3*(2 + lmst22)) + 6*power2(lmC2) + 30*power2(lmsb22) + 6*power2(
     lmst22) + 7*power2(Pi))) - power2(mC2)*(3*mst24*(42 + 2*lmC2*(-12 + 5*
     lmsb22 - lmst22) + 12*lmst22 - 2*lmsb22*(12 + lmst22) + 4*power2(lmC2) + 4
     *power2(lmsb22) - 2*power2(lmst22) + power2(Pi)) + msb22*mst22*(6*lmC2*(-
     12 + 5*lmsb22 - lmst22) + 6*lmsb22*(-18 + lmst22) + 12*power2(lmC2) + 18*
     power2(lmsb22) + 5*(42 + power2(Pi)))) + 3*mst22*(42 + 2*lmC2*(-6 + 3*
     lmsb22 - lmst22) + 2*lmsb22*(-12 + lmst22) + 2*power2(lmC2) + 4*power2(
     lmsb22) + power2(Pi))*power3(mC2) + 3*(42 + 2*lmsb22*(lmC2 + 3*(-4 +
     lmst22)) - 2*(6 + lmC2)*lmst22 + 4*power2(lmsb22) + 2*power2(lmst22) +
     power2(Pi))*power4(mst22)))/(6.*msb22)) + power3(DeltaInv(mst22,mw2,msb22)
     )*(-(invdmst*mt2*power2(mw2)*(sb*power2(At)*(-1 + power2(cb)) - mu2*sb*
     power2(cb) + At*cb*mu*(1 - power2(cb) + power2(sb)))*power2(snb)*(3*mw2*(
     mst26*(42 + 2*lmsb22*(-12 + 5*lmst22 - lmw2) + 12*lmw2 - 2*lmst22*(12 +
     lmw2) + 4*power2(lmsb22) + 4*power2(lmst22) - 2*power2(lmw2) + power2(Pi))
     + mst24*mw2*(42 - 2*lmsb22*(12 + lmst22 - 5*lmw2) - 2*lmst22*(-6 + lmw2) -
     24*lmw2 + 4*power2(lmsb22) - 2*power2(lmst22) + 4*power2(lmw2) + power2(Pi
     ))) + msb22*(2*mst24*mw2*(294 - 36*lmw2 - 18*lmst22*(2 + lmw2) + 30*lmsb22
     *(-6 + lmst22 + lmw2) + 30*power2(lmsb22) + 6*power2(lmst22) + 6*power2(
     lmw2) + 7*power2(Pi)) + mst26*(-6*lmst22*(12 + lmw2) + 6*lmsb22*(-18 + 5*
     lmst22 + lmw2) + 18*power2(lmsb22) + 12*power2(lmst22) + 5*(42 + power2(Pi
     )))) + mst22*power2(mw2)*(-3*mw2*(42 + 2*lmsb22*(lmst22 + 3*(-4 + lmw2)) -
     2*(6 + lmst22)*lmw2 + 4*power2(lmsb22) + 2*power2(lmw2) + power2(Pi)) +
     msb22*(-6*(12 + lmst22)*lmw2 + 6*lmsb22*(-18 + lmst22 + 5*lmw2) + 18*
     power2(lmsb22) + 12*power2(lmw2) + 5*(42 + power2(Pi)))) - 3*(42 - 2*
     lmst22*(6 + lmw2) + 2*lmsb22*(-12 + 3*lmst22 + lmw2) + 4*power2(lmsb22) +
     2*power2(lmst22) + power2(Pi))*power4(mst22)))/(3.*msb22*sb) + (power2(mw2
     )*power2(snb)*(2*At*cb*mu*sb*(-1 + power2(snt)) + power2(At)*(-1 + power2(
     cb))*(-1 + power2(snt)) + mt2*power2(snt) + power2(cb)*(mu2 - mt2*power2(
     snt) - mu2*power2(snt)))*(-3*mw2*(mst26*(42 + 2*lmsb22*(-12 + 5*lmst22 -
     lmw2) + 12*lmw2 - 2*lmst22*(12 + lmw2) + 4*power2(lmsb22) + 4*power2(
     lmst22) - 2*power2(lmw2) + power2(Pi)) + mst24*mw2*(42 - 2*lmsb22*(12 +
     lmst22 - 5*lmw2) - 2*lmst22*(-6 + lmw2) - 24*lmw2 + 4*power2(lmsb22) - 2*
     power2(lmst22) + 4*power2(lmw2) + power2(Pi))) - msb22*(2*mst24*mw2*(294 -
     36*lmw2 - 18*lmst22*(2 + lmw2) + 30*lmsb22*(-6 + lmst22 + lmw2) + 30*
     power2(lmsb22) + 6*power2(lmst22) + 6*power2(lmw2) + 7*power2(Pi)) + mst26
     *(-6*lmst22*(12 + lmw2) + 6*lmsb22*(-18 + 5*lmst22 + lmw2) + 18*power2(
     lmsb22) + 12*power2(lmst22) + 5*(42 + power2(Pi)))) + mst22*power2(mw2)*(3
     *mw2*(42 + 2*lmsb22*(lmst22 + 3*(-4 + lmw2)) - 2*(6 + lmst22)*lmw2 + 4*
     power2(lmsb22) + 2*power2(lmw2) + power2(Pi)) - msb22*(-6*(12 + lmst22)*
     lmw2 + 6*lmsb22*(-18 + lmst22 + 5*lmw2) + 18*power2(lmsb22) + 12*power2(
     lmw2) + 5*(42 + power2(Pi)))) + 3*(42 - 2*lmst22*(6 + lmw2) + 2*lmsb22*(-
     12 + 3*lmst22 + lmw2) + 4*power2(lmsb22) + 2*power2(lmst22) + power2(Pi))*
     power4(mst22)))/(6.*msb22)) + Fin3(mt2,mt2,mA2,Q2)*(-(mt2*DeltaInv(mt2,mt2
     ,mA2)*power2(cb)) - (5*mA2 + 4*mt2)*power2(cb)*power2(mt2)*power2(DeltaInv
     (mt2,mt2,mA2)) - 16*mA2*power2(cb)*power3(DeltaInv(mt2,mt2,mA2))*power4(
     mt2)) + Fin3(mt2,mt2,mz2,Q2)*(mt2*DeltaInv(mt2,mt2,mz2)*(-1 + power2(cb))
     + (4*mt2 + 5*mz2)*(-1 + power2(cb))*power2(mt2)*power2(DeltaInv(mt2,mt2,
     mz2)) + 16*mz2*(-1 + power2(cb))*power3(DeltaInv(mt2,mt2,mz2))*power4(mt2)
     ) + (power2(mu2)*(-1 + power2(snb))*power3(DeltaInv(msb12,mt2,mu2))*(mu2*(
     -(msb14*mu2*(mu2*(42 + 6*(-15 + 4*lmt2)*lmu2 + 6*lmsb12*(9 - 4*lmt2 + lmu2
     ) - 9*power2(lmsb12) + 15*power2(lmu2) + power2(Pi)) + 5*mt2*(210 + 6*
     lmsb12*(-3 + 3*lmt2 - 2*lmu2) - 36*lmu2 + 6*lmt2*(-21 + 4*lmu2) + 3*power2
     (lmsb12) + 21*power2(lmt2) + 6*power2(lmu2) + 5*power2(Pi)))) - msb16*(3*
     mu2*(126 - 84*lmt2 + 10*lmsb12*(-3 + 2*lmt2 - lmu2) + 6*lmu2 + 8*lmt2*lmu2
      + 5*power2(lmsb12) + 14*power2(lmt2) - power2(lmu2) + 3*power2(Pi)) + mt2
     *(672 - 378*lmt2 + 30*lmsb12*(-6 + 3*lmt2 - lmu2) - 18*lmu2 + 36*lmt2*lmu2
      + 30*power2(lmsb12) + 63*power2(lmt2) + 3*power2(lmu2) + 16*power2(Pi))))
     + msb12*(-3*mt2*(-2*(15 + lmsb12)*lmu2 + 2*lmt2*(-21 + lmsb12 + 6*lmu2) +
     7*power2(lmt2) + 5*power2(lmu2) + 2*(42 + power2(Pi))) + mu2*(-6*(9 +
     lmsb12)*lmu2 + 6*lmt2*(-15 + lmsb12 + 4*lmu2) + 15*power2(lmt2) + 9*power2
     (lmu2) + 4*(42 + power2(Pi))))*power3(mu2) + (-(mt2*(42 + 6*lmsb12*(-3 +
     lmt2) - 18*lmt2 + 3*power2(lmsb12) + 3*power2(lmt2) + power2(Pi))) + mu2*(
     210 + 24*lmt2*(-6 + lmu2) - 18*lmu2 + 6*lmsb12*(4*lmt2 - 3*(1 + lmu2)) + 3
     *power2(lmsb12) + 24*power2(lmt2) + 3*power2(lmu2) + 5*power2(Pi)))*power4
     (msb12) + (42 + 6*lmsb12*(-3 + lmt2) - 18*lmt2 + 3*power2(lmsb12) + 3*
     power2(lmt2) + power2(Pi))*power5(msb12)))/(3.*mt2) + (power2(mu2)*power2(
     snb)*power3(DeltaInv(mt2,mu2,msb22))*(mu2*(msb24*mu2*(mu2*(42 + 6*(-15 + 4
     *lmt2)*lmu2 + 6*lmsb22*(9 - 4*lmt2 + lmu2) - 9*power2(lmsb22) + 15*power2(
     lmu2) + power2(Pi)) + 5*mt2*(210 + 6*lmsb22*(-3 + 3*lmt2 - 2*lmu2) - 36*
     lmu2 + 6*lmt2*(-21 + 4*lmu2) + 3*power2(lmsb22) + 21*power2(lmt2) + 6*
     power2(lmu2) + 5*power2(Pi))) + msb26*(3*mu2*(126 - 84*lmt2 + 10*lmsb22*(-
     3 + 2*lmt2 - lmu2) + 6*lmu2 + 8*lmt2*lmu2 + 5*power2(lmsb22) + 14*power2(
     lmt2) - power2(lmu2) + 3*power2(Pi)) + mt2*(672 - 378*lmt2 + 30*lmsb22*(-6
      + 3*lmt2 - lmu2) - 18*lmu2 + 36*lmt2*lmu2 + 30*power2(lmsb22) + 63*power2
     (lmt2) + 3*power2(lmu2) + 16*power2(Pi)))) + msb22*(3*mt2*(-2*(15 + lmsb22
     )*lmu2 + 2*lmt2*(-21 + lmsb22 + 6*lmu2) + 7*power2(lmt2) + 5*power2(lmu2)
     + 2*(42 + power2(Pi))) - mu2*(-6*(9 + lmsb22)*lmu2 + 6*lmt2*(-15 + lmsb22
     + 4*lmu2) + 15*power2(lmt2) + 9*power2(lmu2) + 4*(42 + power2(Pi))))*
     power3(mu2) + (mt2*(42 + 6*lmsb22*(-3 + lmt2) - 18*lmt2 + 3*power2(lmsb22)
     + 3*power2(lmt2) + power2(Pi)) - mu2*(210 + 24*lmt2*(-6 + lmu2) - 18*lmu2
     + 6*lmsb22*(4*lmt2 - 3*(1 + lmu2)) + 3*power2(lmsb22) + 24*power2(lmt2) +
     3*power2(lmu2) + 5*power2(Pi)))*power4(msb22) - (42 + 6*lmsb22*(-3 + lmt2)
     - 18*lmt2 + 3*power2(lmsb22) + 3*power2(lmt2) + power2(Pi))*power5(msb22))
     )/(3.*mt2) + (power2(mu2)*power3(DeltaInv(mt2,mu2,mst12))*(mu2*(mst14*mu2*
     (mu2*(42 + 6*(-15 + 4*lmt2)*lmu2 + 6*lmst12*(9 - 4*lmt2 + lmu2) - 9*power2
     (lmst12) + 15*power2(lmu2) + power2(Pi)) + 5*mt2*(210 + 6*lmst12*(-3 + 3*
     lmt2 - 2*lmu2) - 36*lmu2 + 6*lmt2*(-21 + 4*lmu2) + 3*power2(lmst12) + 21*
     power2(lmt2) + 6*power2(lmu2) + 5*power2(Pi))) + mst16*(3*mu2*(126 - 84*
     lmt2 + 10*lmst12*(-3 + 2*lmt2 - lmu2) + 6*lmu2 + 8*lmt2*lmu2 + 5*power2(
     lmst12) + 14*power2(lmt2) - power2(lmu2) + 3*power2(Pi)) + mt2*(672 - 378*
     lmt2 + 30*lmst12*(-6 + 3*lmt2 - lmu2) - 18*lmu2 + 36*lmt2*lmu2 + 30*power2
     (lmst12) + 63*power2(lmt2) + 3*power2(lmu2) + 16*power2(Pi)))) + mst12*(3*
     mt2*(-2*(15 + lmst12)*lmu2 + 2*lmt2*(-21 + lmst12 + 6*lmu2) + 7*power2(
     lmt2) + 5*power2(lmu2) + 2*(42 + power2(Pi))) - mu2*(-6*(9 + lmst12)*lmu2
     + 6*lmt2*(-15 + lmst12 + 4*lmu2) + 15*power2(lmt2) + 9*power2(lmu2) + 4*(
     42 + power2(Pi))))*power3(mu2) + (mt2*(42 + 6*lmst12*(-3 + lmt2) - 18*lmt2
      + 3*power2(lmst12) + 3*power2(lmt2) + power2(Pi)) - mu2*(210 + 24*lmt2*(-
     6 + lmu2) - 18*lmu2 + 6*lmst12*(4*lmt2 - 3*(1 + lmu2)) + 3*power2(lmst12)
     + 24*power2(lmt2) + 3*power2(lmu2) + 5*power2(Pi)))*power4(mst12) - (42 +
     6*lmst12*(-3 + lmt2) - 18*lmt2 + 3*power2(lmst12) + 3*power2(lmt2) +
     power2(Pi))*power5(mst12)))/(3.*mt2) + (power2(mu2)*power3(DeltaInv(mst22,
     mt2,mu2))*(mu2*(mst24*mu2*(mu2*(42 + 6*(-15 + 4*lmt2)*lmu2 + 6*lmst22*(9 -
     4*lmt2 + lmu2) - 9*power2(lmst22) + 15*power2(lmu2) + power2(Pi)) + 5*mt2*
     (210 + 6*lmst22*(-3 + 3*lmt2 - 2*lmu2) - 36*lmu2 + 6*lmt2*(-21 + 4*lmu2) +
     3*power2(lmst22) + 21*power2(lmt2) + 6*power2(lmu2) + 5*power2(Pi))) +
     mst26*(3*mu2*(126 - 84*lmt2 + 10*lmst22*(-3 + 2*lmt2 - lmu2) + 6*lmu2 + 8*
     lmt2*lmu2 + 5*power2(lmst22) + 14*power2(lmt2) - power2(lmu2) + 3*power2(
     Pi)) + mt2*(672 - 378*lmt2 + 30*lmst22*(-6 + 3*lmt2 - lmu2) - 18*lmu2 + 36
     *lmt2*lmu2 + 30*power2(lmst22) + 63*power2(lmt2) + 3*power2(lmu2) + 16*
     power2(Pi)))) + mst22*(3*mt2*(-2*(15 + lmst22)*lmu2 + 2*lmt2*(-21 + lmst22
      + 6*lmu2) + 7*power2(lmt2) + 5*power2(lmu2) + 2*(42 + power2(Pi))) - mu2*
     (-6*(9 + lmst22)*lmu2 + 6*lmt2*(-15 + lmst22 + 4*lmu2) + 15*power2(lmt2) +
     9*power2(lmu2) + 4*(42 + power2(Pi))))*power3(mu2) + (mt2*(42 + 6*lmst22*(
     -3 + lmt2) - 18*lmt2 + 3*power2(lmst22) + 3*power2(lmt2) + power2(Pi)) -
     mu2*(210 + 24*lmt2*(-6 + lmu2) - 18*lmu2 + 6*lmst22*(4*lmt2 - 3*(1 + lmu2)
     ) + 3*power2(lmst22) + 24*power2(lmt2) + 3*power2(lmu2) + 5*power2(Pi)))*
     power4(mst22) - (42 + 6*lmst22*(-3 + lmt2) - 18*lmt2 + 3*power2(lmst22) +
     3*power2(lmt2) + power2(Pi))*power5(mst22)))/(3.*mt2) - 8*mA2*power2(cb)*(
     42 + 8*lmA2*(-3 + lmt2) - 12*lmt2 + 4*power2(lmA2) + power2(Pi))*power3(
     DeltaInv(mt2,mt2,mA2))*power5(mt2) + 8*mz2*(-1 + power2(cb))*(42 - 24*lmz2
      + 4*lmt2*(-3 + 2*lmz2) + 4*power2(lmz2) + power2(Pi))*power3(DeltaInv(mt2
     ,mt2,mz2))*power5(mt2) + (-36 - 72*lmt2 - 36*lmu2 - 90*lmz2 - (9*mh2)/
     mst12 - (9*mh2)/mst22 - (18*msb12)/mst22 + (18*lmsb12*msb12)/mst22 - (18*
     mst12)/mst22 + (18*lmst12*mst12)/mst22 - (18*mst22)/msb12 + (18*lmst22*
     mst22)/msb12 - (18*mst22)/mst12 + (18*lmst22*mst22)/mst12 - (36*mt2)/mh2 +
     (24*lmt2*mt2)/mh2 + (18*mt2)/msb12 + (18*lmsb12*mt2)/msb12 - (36*lmt2*mt2)
     /msb12 - (18*mt2)/mst12 + (18*lmsb12*mt2)/mst12 - (36*lmt2*mt2)/mst12 - (
     36*lmt2*mt2)/mst22 + (36*lmu2*mst22)/(mst22 - mu2) + (18*mu2)/msb12 + (18*
     lmsb12*mu2)/msb12 - (36*lmt2*mu2)/msb12 + (18*mu2)/mst12 + (18*lmst12*mu2)
     /mst12 - (36*lmt2*mu2)/mst12 + (36*mu2)/mst22 + (18*lmsb12*mu2)/mst22 + (
     18*lmst12*mu2)/mst22 - (36*lmt2*mu2)/mst22 - (36*lmu2*mu2)/mst22 - (18*
     lmC2*mC2*mu2)/(msb12*mst22) + (18*lmsb12*mC2*mu2)/(msb12*mst22) - (9*lmA2*
     mA2*mu2)/(mst12*mst22) + (9*lmst12*mA2*mu2)/(mst12*mst22) - (9*lmH2*mH2*
     mu2)/(mst12*mst22) + (9*lmst12*mH2*mu2)/(mst12*mst22) - (108*mu2)/mt2 + (
     36*lmst22*mst22)/(-mst22 + mu2) + 35*invdtw*mw2 - 84*invdtw*lmt2*mw2 + 60*
     invdtw*lmw2*mw2 - (18*mw2)/msb12 + (18*lmw2*mw2)/msb12 - (18*mw2)/mst22 +
     (18*lmw2*mw2)/mst22 - (35*mw2)/mt2 + (12*lmt2*mw2)/mt2 + (12*lmw2*mw2)/mt2
      + (18*lmsb12*mt2*mw2)/(msb12*mst12) - (18*lmw2*mt2*mw2)/(msb12*mst12) - (
     9*mz2)/mst12 + (9*lmz2*mz2)/mst12 - (9*mz2)/mst22 + (9*lmz2*mz2)/mst22 + (
     18*At*ca*lmst12*mh2*mu*sa)/(mst12*mst22) + (18*At*ca*lmH2*mH2*mu*sa)/(
     mst12*mst22) - (18*At*ca*lmst12*mH2*mu*sa)/(mst12*mst22) - (18*power2(At))
     /msb12 + (18*lmsb12*power2(At))/msb12 - (18*power2(At))/mst12 + (18*lmst12
     *power2(At))/mst12 - (36*power2(At))/mst22 + (18*lmsb12*power2(At))/mst22
     + (18*lmst12*power2(At))/mst22 + (9*lmst12*mh2*power2(At))/(mst12*mst22) +
     (18*lmsb12*mw2*power2(At))/(msb12*mst22) - (18*lmw2*mw2*power2(At))/(msb12
     *mst22) + (9*lmst12*mz2*power2(At))/(mst12*mst22) - (9*lmz2*mz2*power2(At)
     )/(mst12*mst22) - (36*lmt2*power2(mu2))/(msb12*mt2) + (36*lmu2*power2(mu2)
     )/(msb12*mt2) - (36*lmt2*power2(mu2))/(mst12*mt2) + (36*lmu2*power2(mu2))/
     (mst12*mt2) - (36*lmt2*power2(mu2))/(mst22*mt2) + (36*lmu2*power2(mu2))/(
     mst22*mt2) - 35*power2(invdtw)*power2(mw2) + 12*lmt2*power2(invdtw)*power2
     (mw2) + 12*lmw2*power2(invdtw)*power2(mw2) + 54*lmH2*power2(sa) + (9*mh2*
     power2(sa))/mst12 - (9*mH2*power2(sa))/mst12 + (9*lmH2*mH2*power2(sa))/
     mst12 + (9*mh2*power2(sa))/mst22 - (9*mH2*power2(sa))/mst22 + (9*lmH2*mH2*
     power2(sa))/mst22 + (36*mt2*power2(sa))/mh2 - (24*lmt2*mt2*power2(sa))/mh2
      - (36*mt2*power2(sa))/mH2 + (24*lmH2*mt2*power2(sa))/mH2 + (24*lmt2*mt2*
     power2(sa))/mH2 + (36*lmH2*mt2*power2(sa))/mst12 + (36*lmH2*mt2*power2(sa)
     )/mst22 + (9*lmst12*mh2*mu2*power2(sa))/(mst12*mst22) + (9*lmH2*mH2*mu2*
     power2(sa))/(mst12*mst22) - (9*lmst12*mH2*mu2*power2(sa))/(mst12*mst22) -
     (9*lmst12*mh2*power2(At)*power2(sa))/(mst12*mst22) - (9*lmH2*mH2*power2(At
     )*power2(sa))/(mst12*mst22) + (9*lmst12*mH2*power2(At)*power2(sa))/(mst12*
     mst22) + (18*msb12*power2(snb))/mst22 - (18*lmsb12*msb12*power2(snb))/
     mst22 - (18*msb22*power2(snb))/mst22 + (18*lmsb22*msb22*power2(snb))/mst22
      + (18*mst22*power2(snb))/msb12 - (18*lmst22*mst22*power2(snb))/msb12 - (
     18*mst22*power2(snb))/msb22 + (18*lmst22*mst22*power2(snb))/msb22 - (18*
     mt2*power2(snb))/msb12 - (18*lmsb12*mt2*power2(snb))/msb12 + (36*lmt2*mt2*
     power2(snb))/msb12 + (18*mt2*power2(snb))/msb22 + (18*lmsb22*mt2*power2(
     snb))/msb22 - (36*lmt2*mt2*power2(snb))/msb22 - (18*lmsb12*mt2*power2(snb)
     )/mst12 + (18*lmsb22*mt2*power2(snb))/mst12 - (18*mu2*power2(snb))/msb12 -
     (18*lmsb12*mu2*power2(snb))/msb12 + (36*lmt2*mu2*power2(snb))/msb12 + (18*
     mu2*power2(snb))/msb22 + (18*lmsb22*mu2*power2(snb))/msb22 - (36*lmt2*mu2*
     power2(snb))/msb22 - (18*lmsb12*mu2*power2(snb))/mst22 + (18*lmsb22*mu2*
     power2(snb))/mst22 + (18*lmC2*mC2*mu2*power2(snb))/(msb12*mst22) - (18*
     lmsb12*mC2*mu2*power2(snb))/(msb12*mst22) - (18*lmC2*mC2*mu2*power2(snb))/
     (msb22*mst22) + (18*lmsb22*mC2*mu2*power2(snb))/(msb22*mst22) + (18*mw2*
     power2(snb))/msb12 - (18*lmw2*mw2*power2(snb))/msb12 - (18*mw2*power2(snb)
     )/msb22 + (18*lmw2*mw2*power2(snb))/msb22 - (18*lmsb12*mt2*mw2*power2(snb)
     )/(msb12*mst12) + (18*lmw2*mt2*mw2*power2(snb))/(msb12*mst12) + (18*lmsb22
     *mt2*mw2*power2(snb))/(msb22*mst12) - (18*lmw2*mt2*mw2*power2(snb))/(msb22
     *mst12) + (18*power2(At)*power2(snb))/msb12 - (18*lmsb12*power2(At)*power2
     (snb))/msb12 - (18*power2(At)*power2(snb))/msb22 + (18*lmsb22*power2(At)*
     power2(snb))/msb22 - (18*lmsb12*power2(At)*power2(snb))/mst22 + (18*lmsb22
     *power2(At)*power2(snb))/mst22 - (18*lmsb12*mw2*power2(At)*power2(snb))/(
     msb12*mst22) + (18*lmw2*mw2*power2(At)*power2(snb))/(msb12*mst22) + (18*
     lmsb22*mw2*power2(At)*power2(snb))/(msb22*mst22) - (18*lmw2*mw2*power2(At)
     *power2(snb))/(msb22*mst22) + (36*lmt2*power2(mu2)*power2(snb))/(msb12*mt2
     ) - (36*lmu2*power2(mu2)*power2(snb))/(msb12*mt2) - (36*lmt2*power2(mu2)*
     power2(snb))/(msb22*mt2) + (36*lmu2*power2(mu2)*power2(snb))/(msb22*mt2) -
     288*power2(snt) + 144*lmst12*power2(snt) + 144*lmst22*power2(snt) - (18*
     msb12*power2(snt))/mst12 + (18*lmsb12*msb12*power2(snt))/mst12 - (18*mst12
     *power2(snt))/msb12 + (18*lmst12*mst12*power2(snt))/msb12 + (18*msb12*
     power2(snt))/mst22 - (18*lmsb12*msb12*power2(snt))/mst22 + (144*mst12*
     power2(snt))/mst22 - (144*lmst12*mst12*power2(snt))/mst22 + (18*mst22*
     power2(snt))/msb12 - (18*lmst22*mst22*power2(snt))/msb12 + (144*mst22*
     power2(snt))/mst12 - (144*lmst22*mst22*power2(snt))/mst12 + (18*mt2*power2
     (snt))/mst12 - (18*lmsb12*mt2*power2(snt))/mst12 - (18*mt2*power2(snt))/
     mst22 + (18*lmsb12*mt2*power2(snt))/mst22 - (36*lmst12*mst12*power2(snt))/
     (mst12 - mu2) + (36*lmu2*mst12*power2(snt))/(mst12 - mu2) + (36*lmst22*
     mst22*power2(snt))/(mst22 - mu2) + (18*mu2*power2(snt))/mst12 + (36*lmH2*
     mu2*power2(snt))/mst12 + (18*lmsb12*mu2*power2(snt))/mst12 - (36*lmst12*
     mu2*power2(snt))/mst12 - (36*lmu2*mu2*power2(snt))/mst12 - (18*lmC2*mC2*
     mu2*power2(snt))/(msb12*mst12) + (18*lmsb12*mC2*mu2*power2(snt))/(msb12*
     mst12) - (18*mu2*power2(snt))/mst22 + (36*lmH2*mu2*power2(snt))/mst22 - (
     18*lmsb12*mu2*power2(snt))/mst22 - (36*lmst12*mu2*power2(snt))/mst22 + (36
     *lmu2*mu2*power2(snt))/mst22 + (18*lmC2*mC2*mu2*power2(snt))/(msb12*mst22)
     - (18*lmsb12*mC2*mu2*power2(snt))/(msb12*mst22) + (36*lmH2*mH2*mu2*power2(
     snt))/(mst12*mst22) - (36*lmst12*mH2*mu2*power2(snt))/(mst12*mst22) + (36*
     lmu2*mst22*power2(snt))/(-mst22 + mu2) - (18*mw2*power2(snt))/mst12 + (18*
     lmw2*mw2*power2(snt))/mst12 + (18*mw2*power2(snt))/mst22 - (18*lmw2*mw2*
     power2(snt))/mst22 - (18*lmsb12*mt2*mw2*power2(snt))/(msb12*mst12) + (18*
     lmw2*mt2*mw2*power2(snt))/(msb12*mst12) + (18*lmsb12*mt2*mw2*power2(snt))/
     (msb12*mst22) - (18*lmw2*mt2*mw2*power2(snt))/(msb12*mst22) - (72*At*ca*
     lmH2*mu*sa*power2(snt))/mst12 - (72*At*ca*lmH2*mu*sa*power2(snt))/mst22 -
     (72*At*ca*lmst12*mh2*mu*sa*power2(snt))/(mst12*mst22) - (72*At*ca*lmH2*mH2
     *mu*sa*power2(snt))/(mst12*mst22) + (72*At*ca*lmst12*mH2*mu*sa*power2(snt)
     )/(mst12*mst22) - (18*power2(At)*power2(snt))/mst12 + (18*lmsb12*power2(At
     )*power2(snt))/mst12 - (36*lmst12*power2(At)*power2(snt))/mst12 + (18*
     power2(At)*power2(snt))/mst22 - (18*lmsb12*power2(At)*power2(snt))/mst22 -
     (36*lmst12*power2(At)*power2(snt))/mst22 - (36*lmst12*mh2*power2(At)*
     power2(snt))/(mst12*mst22) + (18*lmsb12*mw2*power2(At)*power2(snt))/(msb12
     *mst12) - (18*lmw2*mw2*power2(At)*power2(snt))/(msb12*mst12) - (18*lmsb12*
     mw2*power2(At)*power2(snt))/(msb12*mst22) + (18*lmw2*mw2*power2(At)*power2
     (snt))/(msb12*mst22) - (36*lmH2*mu2*power2(sa)*power2(snt))/mst12 - (36*
     lmH2*mu2*power2(sa)*power2(snt))/mst22 - (36*lmst12*mh2*mu2*power2(sa)*
     power2(snt))/(mst12*mst22) - (36*lmH2*mH2*mu2*power2(sa)*power2(snt))/(
     mst12*mst22) + (36*lmst12*mH2*mu2*power2(sa)*power2(snt))/(mst12*mst22) +
     (36*lmH2*power2(At)*power2(sa)*power2(snt))/mst12 + (36*lmH2*power2(At)*
     power2(sa)*power2(snt))/mst22 + (36*lmst12*mh2*power2(At)*power2(sa)*
     power2(snt))/(mst12*mst22) + (36*lmH2*mH2*power2(At)*power2(sa)*power2(snt
     ))/(mst12*mst22) - (36*lmst12*mH2*power2(At)*power2(sa)*power2(snt))/(
     mst12*mst22) + (18*msb12*power2(snb)*power2(snt))/mst12 - (18*lmsb12*msb12
     *power2(snb)*power2(snt))/mst12 - (18*msb22*power2(snb)*power2(snt))/mst12
      + (18*lmsb22*msb22*power2(snb)*power2(snt))/mst12 + (18*mst12*power2(snb)
     *power2(snt))/msb12 - (18*lmst12*mst12*power2(snb)*power2(snt))/msb12 - (
     18*mst12*power2(snb)*power2(snt))/msb22 + (18*lmst12*mst12*power2(snb)*
     power2(snt))/msb22 - (18*msb12*power2(snb)*power2(snt))/mst22 + (18*lmsb12
     *msb12*power2(snb)*power2(snt))/mst22 + (18*msb22*power2(snb)*power2(snt))
     /mst22 - (18*lmsb22*msb22*power2(snb)*power2(snt))/mst22 - (18*mst22*
     power2(snb)*power2(snt))/msb12 + (18*lmst22*mst22*power2(snb)*power2(snt))
     /msb12 + (18*mst22*power2(snb)*power2(snt))/msb22 - (18*lmst22*mst22*
     power2(snb)*power2(snt))/msb22 + (18*lmsb12*mt2*power2(snb)*power2(snt))/
     mst12 - (18*lmsb22*mt2*power2(snb)*power2(snt))/mst12 - (18*lmsb12*mt2*
     power2(snb)*power2(snt))/mst22 + (18*lmsb22*mt2*power2(snb)*power2(snt))/
     mst22 - (18*lmsb12*mu2*power2(snb)*power2(snt))/mst12 + (18*lmsb22*mu2*
     power2(snb)*power2(snt))/mst12 + (18*lmC2*mC2*mu2*power2(snb)*power2(snt))
     /(msb12*mst12) - (18*lmsb12*mC2*mu2*power2(snb)*power2(snt))/(msb12*mst12)
     - (18*lmC2*mC2*mu2*power2(snb)*power2(snt))/(msb22*mst12) + (18*lmsb22*mC2
     *mu2*power2(snb)*power2(snt))/(msb22*mst12) + (18*lmsb12*mu2*power2(snb)*
     power2(snt))/mst22 - (18*lmsb22*mu2*power2(snb)*power2(snt))/mst22 - (18*
     lmC2*mC2*mu2*power2(snb)*power2(snt))/(msb12*mst22) + (18*lmsb12*mC2*mu2*
     power2(snb)*power2(snt))/(msb12*mst22) + (18*lmC2*mC2*mu2*power2(snb)*
     power2(snt))/(msb22*mst22) - (18*lmsb22*mC2*mu2*power2(snb)*power2(snt))/(
     msb22*mst22) + (18*lmsb12*mt2*mw2*power2(snb)*power2(snt))/(msb12*mst12) -
     (18*lmw2*mt2*mw2*power2(snb)*power2(snt))/(msb12*mst12) - (18*lmsb22*mt2*
     mw2*power2(snb)*power2(snt))/(msb22*mst12) + (18*lmw2*mt2*mw2*power2(snb)*
     power2(snt))/(msb22*mst12) - (18*lmsb12*mt2*mw2*power2(snb)*power2(snt))/(
     msb12*mst22) + (18*lmw2*mt2*mw2*power2(snb)*power2(snt))/(msb12*mst22) + (
     18*lmsb22*mt2*mw2*power2(snb)*power2(snt))/(msb22*mst22) - (18*lmw2*mt2*
     mw2*power2(snb)*power2(snt))/(msb22*mst22) - (18*lmsb12*power2(At)*power2(
     snb)*power2(snt))/mst12 + (18*lmsb22*power2(At)*power2(snb)*power2(snt))/
     mst12 + (18*lmsb12*power2(At)*power2(snb)*power2(snt))/mst22 - (18*lmsb22*
     power2(At)*power2(snb)*power2(snt))/mst22 - (18*lmsb12*mw2*power2(At)*
     power2(snb)*power2(snt))/(msb12*mst12) + (18*lmw2*mw2*power2(At)*power2(
     snb)*power2(snt))/(msb12*mst12) + (18*lmsb22*mw2*power2(At)*power2(snb)*
     power2(snt))/(msb22*mst12) - (18*lmw2*mw2*power2(At)*power2(snb)*power2(
     snt))/(msb22*mst12) + (18*lmsb12*mw2*power2(At)*power2(snb)*power2(snt))/(
     msb12*mst22) - (18*lmw2*mw2*power2(At)*power2(snb)*power2(snt))/(msb12*
     mst22) - (18*lmsb22*mw2*power2(At)*power2(snb)*power2(snt))/(msb22*mst22)
     + (18*lmw2*mw2*power2(At)*power2(snb)*power2(snt))/(msb22*mst22) - (18*At*
     cb*mu*sb*(lmA2*mA2*msb12*msb22 + 2*lmC2*mC2*msb22*mst12 - 2*lmsb12*mC2*
     msb22*mst12 + 2*lmsb12*msb22*mst12*mw2 - 2*lmw2*msb22*mst12*mw2 - lmz2*
     msb12*msb22*mz2 + lmst12*msb12*msb22*(-mA2 + mz2) + 2*lmC2*mC2*msb12*mst12
     *power2(snb) - 2*lmsb22*mC2*msb12*mst12*power2(snb) - 2*lmC2*mC2*msb22*
     mst12*power2(snb) + 2*lmsb12*mC2*msb22*mst12*power2(snb) + 2*lmsb22*msb12*
     mst12*mw2*power2(snb) - 2*lmw2*msb12*mst12*mw2*power2(snb) - 2*lmsb12*
     msb22*mst12*mw2*power2(snb) + 2*lmw2*msb22*mst12*mw2*power2(snb) - 2*lmC2*
     mC2*msb22*mst12*power2(snt) + 2*lmsb12*mC2*msb22*mst12*power2(snt) + 2*
     lmC2*mC2*msb22*mst22*power2(snt) - 2*lmsb12*mC2*msb22*mst22*power2(snt) -
     2*lmsb12*msb22*mst12*mw2*power2(snt) + 2*lmw2*msb22*mst12*mw2*power2(snt)
     + 2*lmsb12*msb22*mst22*mw2*power2(snt) - 2*lmw2*msb22*mst22*mw2*power2(snt
     ) - 2*lmC2*mC2*msb12*mst12*power2(snb)*power2(snt) + 2*lmsb22*mC2*msb12*
     mst12*power2(snb)*power2(snt) + 2*lmC2*mC2*msb22*mst12*power2(snb)*power2(
     snt) - 2*lmsb12*mC2*msb22*mst12*power2(snb)*power2(snt) + 2*lmC2*mC2*msb12
     *mst22*power2(snb)*power2(snt) - 2*lmsb22*mC2*msb12*mst22*power2(snb)*
     power2(snt) - 2*lmC2*mC2*msb22*mst22*power2(snb)*power2(snt) + 2*lmsb12*
     mC2*msb22*mst22*power2(snb)*power2(snt) - 2*lmsb22*msb12*mst12*mw2*power2(
     snb)*power2(snt) + 2*lmw2*msb12*mst12*mw2*power2(snb)*power2(snt) + 2*
     lmsb12*msb22*mst12*mw2*power2(snb)*power2(snt) - 2*lmw2*msb22*mst12*mw2*
     power2(snb)*power2(snt) + 2*lmsb22*msb12*mst22*mw2*power2(snb)*power2(snt)
     - 2*lmw2*msb12*mst22*mw2*power2(snb)*power2(snt) - 2*lmsb12*msb22*mst22*
     mw2*power2(snb)*power2(snt) + 2*lmw2*msb22*mst22*mw2*power2(snb)*power2(
     snt)))/(msb12*msb22*mst12*mst22) - (3*lmh2*(8*mst12*mst22*mt2*(-1 + power2
     (sa)) + 6*mh2*(3*mst12*mst22*(-1 + power2(sa)) + 2*mst22*mt2*(-1 + power2(
     sa)) + 2*mst22*mu2*power2(sa)*(-1 + power2(snt))*power2(snt) + 2*mst12*(
     mt2*(-1 + power2(sa)) + mu2*power2(sa)*(-1 + power2(snt))*power2(snt))) +
     6*At*ca*mh2*mu*sa*(4*(mst12 + mst22)*(-1 + power2(snt))*power2(snt) + mh2*
     power2(1 - 2*power2(snt))) - 3*mh2*power2(At)*(-1 + power2(sa))*(4*(mst12
     + mst22)*(-1 + power2(snt))*power2(snt) + mh2*power2(1 - 2*power2(snt))) +
     3*power2(mh2)*(mst12*(-1 + power2(sa)) + mst22*(-1 + power2(sa)) + mu2*
     power2(sa)*power2(1 - 2*power2(snt)))))/(mh2*mst12*mst22) + (45*power2(
     invdtw)*power3(mw2))/mt2 + (12*lmt2*power2(invdtw)*power3(mw2))/mt2 + (12*
     lmw2*power2(invdtw)*power3(mw2))/mt2 - 10*power3(invdtw)*power3(mw2) - 24*
     lmt2*power3(invdtw)*power3(mw2) - 24*lmw2*power3(invdtw)*power3(mw2) + 10*
     power4(invdtw)*power4(mw2) + 24*lmt2*power4(invdtw)*power4(mw2) + 24*lmw2*
     power4(invdtw)*power4(mw2) + 288*power4(snt) - 144*lmst12*power4(snt) -
     144*lmst22*power4(snt) - (144*mst12*power4(snt))/mst22 + (144*lmst12*mst12
     *power4(snt))/mst22 - (144*mst22*power4(snt))/mst12 + (144*lmst22*mst22*
     power4(snt))/mst12 - (36*lmH2*mu2*power4(snt))/mst12 + (36*lmst12*mu2*
     power4(snt))/mst12 - (36*lmH2*mu2*power4(snt))/mst22 + (36*lmst12*mu2*
     power4(snt))/mst22 - (36*lmH2*mH2*mu2*power4(snt))/(mst12*mst22) + (36*
     lmst12*mH2*mu2*power4(snt))/(mst12*mst22) + (72*At*ca*lmH2*mu*sa*power4(
     snt))/mst12 + (72*At*ca*lmH2*mu*sa*power4(snt))/mst22 + (72*At*ca*lmst12*
     mh2*mu*sa*power4(snt))/(mst12*mst22) + (72*At*ca*lmH2*mH2*mu*sa*power4(snt
     ))/(mst12*mst22) - (72*At*ca*lmst12*mH2*mu*sa*power4(snt))/(mst12*mst22) +
     (36*lmst12*power2(At)*power4(snt))/mst12 + (36*lmst12*power2(At)*power4(
     snt))/mst22 + (36*lmst12*mh2*power2(At)*power4(snt))/(mst12*mst22) + (36*
     lmH2*mu2*power2(sa)*power4(snt))/mst12 + (36*lmH2*mu2*power2(sa)*power4(
     snt))/mst22 + (36*lmst12*mh2*mu2*power2(sa)*power4(snt))/(mst12*mst22) + (
     36*lmH2*mH2*mu2*power2(sa)*power4(snt))/(mst12*mst22) - (36*lmst12*mH2*mu2
     *power2(sa)*power4(snt))/(mst12*mst22) - (36*lmH2*power2(At)*power2(sa)*
     power4(snt))/mst12 - (36*lmH2*power2(At)*power2(sa)*power4(snt))/mst22 - (
     36*lmst12*mh2*power2(At)*power2(sa)*power4(snt))/(mst12*mst22) - (36*lmH2*
     mH2*power2(At)*power2(sa)*power4(snt))/(mst12*mst22) + (36*lmst12*mH2*
     power2(At)*power2(sa)*power4(snt))/(mst12*mst22) - (223*power4(invdtw)*
     power5(mw2))/mt2 + (42*lmt2*power4(invdtw)*power5(mw2))/mt2 + (42*lmw2*
     power4(invdtw)*power5(mw2))/mt2 + 213*power5(invdtw)*power5(mw2) - 66*lmt2
     *power5(invdtw)*power5(mw2) - 66*lmw2*power5(invdtw)*power5(mw2) - 213*
     power6(invdtw)*power6(mw2) + 66*lmt2*power6(invdtw)*power6(mw2) + 66*lmw2*
     power6(invdtw)*power6(mw2) + (213*power6(invdtw)*power7(mw2))/mt2 - (66*
     lmt2*power6(invdtw)*power7(mw2))/mt2 - (66*lmw2*power6(invdtw)*power7(mw2)
     )/mt2 + (power2(cb)*(-9*mA2*msb12*msb22*mst12*mt2 + 9*lmA2*mA2*msb12*msb22
     *mst12*mt2 - 9*mA2*msb12*msb22*mst22*mt2 + 9*lmA2*mA2*msb12*msb22*mst22*
     mt2 - 90*lmA2*msb12*msb22*mst12*mst22*mt2 + 90*lmz2*msb12*msb22*mst12*
     mst22*mt2 + 9*lmA2*mA2*msb12*msb22*mt2*mu2 - 9*lmst12*mA2*msb12*msb22*mt2*
     mu2 + 35*msb12*msb22*mst12*mst22*mw2 - 12*lmt2*msb12*msb22*mst12*mst22*mw2
      - 12*lmw2*msb12*msb22*mst12*mst22*mw2 + 18*msb12*msb22*mst12*mt2*mw2 - 18
     *lmw2*msb12*msb22*mst12*mt2*mw2 + 18*msb22*mst12*mst22*mt2*mw2 - 18*lmw2*
     msb22*mst12*mst22*mt2*mw2 - 35*invdtw*msb12*msb22*mst12*mst22*mt2*mw2 + 84
     *invdtw*lmt2*msb12*msb22*mst12*mst22*mt2*mw2 - 60*invdtw*lmw2*msb12*msb22*
     mst12*mst22*mt2*mw2 + 18*lmsb12*msb22*mst12*mt2*mu2*mw2 - 18*lmw2*msb22*
     mst12*mt2*mu2*mw2 + 9*msb12*msb22*mst12*mt2*mz2 - 9*lmz2*msb12*msb22*mst12
     *mt2*mz2 + 9*msb12*msb22*mst22*mt2*mz2 - 9*lmz2*msb12*msb22*mst22*mt2*mz2
     + 9*lmst12*msb12*msb22*mt2*mu2*mz2 - 9*lmz2*msb12*msb22*mt2*mu2*mz2 + (-35
      + 12*lmC2 + 12*lmt2)*msb12*msb22*mst12*mst22*mt2*power2(invdct)*power2(
     mC2) - 18*lmsb12*msb22*mst22*mw2*power2(mt2) + 18*lmw2*msb22*mst22*mw2*
     power2(mt2) + 35*msb12*msb22*mst12*mst22*mt2*power2(invdtw)*power2(mw2) -
     12*lmt2*msb12*msb22*mst12*mst22*mt2*power2(invdtw)*power2(mw2) - 12*lmw2*
     msb12*msb22*mst12*mst22*mt2*power2(invdtw)*power2(mw2) + 18*msb12*mst12*
     mst22*mt2*mw2*power2(snb) - 18*lmw2*msb12*mst12*mst22*mt2*mw2*power2(snb)
     - 18*msb22*mst12*mst22*mt2*mw2*power2(snb) + 18*lmw2*msb22*mst12*mst22*mt2
     *mw2*power2(snb) + 18*lmsb22*msb12*mst12*mt2*mu2*mw2*power2(snb) - 18*lmw2
     *msb12*mst12*mt2*mu2*mw2*power2(snb) - 18*lmsb12*msb22*mst12*mt2*mu2*mw2*
     power2(snb) + 18*lmw2*msb22*mst12*mt2*mu2*mw2*power2(snb) - 18*lmsb22*
     msb12*mst22*mw2*power2(mt2)*power2(snb) + 18*lmw2*msb12*mst22*mw2*power2(
     mt2)*power2(snb) + 18*lmsb12*msb22*mst22*mw2*power2(mt2)*power2(snb) - 18*
     lmw2*msb22*mst22*mw2*power2(mt2)*power2(snb) - 18*msb12*msb22*mst12*mt2*
     mw2*power2(snt) + 18*lmw2*msb12*msb22*mst12*mt2*mw2*power2(snt) + 18*msb12
     *msb22*mst22*mt2*mw2*power2(snt) - 18*lmw2*msb12*msb22*mst22*mt2*mw2*
     power2(snt) - 18*lmsb12*msb22*mst12*mt2*mu2*mw2*power2(snt) + 18*lmw2*
     msb22*mst12*mt2*mu2*mw2*power2(snt) + 18*lmsb12*msb22*mst22*mt2*mu2*mw2*
     power2(snt) - 18*lmw2*msb22*mst22*mt2*mu2*mw2*power2(snt) - 18*lmsb12*
     msb22*mst12*mw2*power2(mt2)*power2(snt) + 18*lmw2*msb22*mst12*mw2*power2(
     mt2)*power2(snt) + 18*lmsb12*msb22*mst22*mw2*power2(mt2)*power2(snt) - 18*
     lmw2*msb22*mst22*mw2*power2(mt2)*power2(snt) - 18*lmsb22*msb12*mst12*mt2*
     mu2*mw2*power2(snb)*power2(snt) + 18*lmw2*msb12*mst12*mt2*mu2*mw2*power2(
     snb)*power2(snt) + 18*lmsb12*msb22*mst12*mt2*mu2*mw2*power2(snb)*power2(
     snt) - 18*lmw2*msb22*mst12*mt2*mu2*mw2*power2(snb)*power2(snt) + 18*lmsb22
     *msb12*mst22*mt2*mu2*mw2*power2(snb)*power2(snt) - 18*lmw2*msb12*mst22*mt2
     *mu2*mw2*power2(snb)*power2(snt) - 18*lmsb12*msb22*mst22*mt2*mu2*mw2*
     power2(snb)*power2(snt) + 18*lmw2*msb22*mst22*mt2*mu2*mw2*power2(snb)*
     power2(snt) - 18*lmsb22*msb12*mst12*mw2*power2(mt2)*power2(snb)*power2(snt
     ) + 18*lmw2*msb12*mst12*mw2*power2(mt2)*power2(snb)*power2(snt) + 18*
     lmsb12*msb22*mst12*mw2*power2(mt2)*power2(snb)*power2(snt) - 18*lmw2*msb22
     *mst12*mw2*power2(mt2)*power2(snb)*power2(snt) + 18*lmsb22*msb12*mst22*mw2
     *power2(mt2)*power2(snb)*power2(snt) - 18*lmw2*msb12*mst22*mw2*power2(mt2)
     *power2(snb)*power2(snt) - 18*lmsb12*msb22*mst22*mw2*power2(mt2)*power2(
     snb)*power2(snt) + 18*lmw2*msb22*mst22*mw2*power2(mt2)*power2(snb)*power2(
     snt) - 9*mt2*power2(At)*(lmA2*mA2*msb12*msb22 + 2*lmsb12*msb22*mst12*mw2 -
     2*lmw2*msb22*mst12*mw2 - lmz2*msb12*msb22*mz2 + lmst12*msb12*msb22*(-mA2 +
     mz2) + 2*lmsb22*msb12*mst12*mw2*power2(snb) - 2*lmw2*msb12*mst12*mw2*
     power2(snb) - 2*lmsb12*msb22*mst12*mw2*power2(snb) + 2*lmw2*msb22*mst12*
     mw2*power2(snb) - 2*lmsb12*msb22*mst12*mw2*power2(snt) + 2*lmw2*msb22*
     mst12*mw2*power2(snt) + 2*lmsb12*msb22*mst22*mw2*power2(snt) - 2*lmw2*
     msb22*mst22*mw2*power2(snt) - 2*lmsb22*msb12*mst12*mw2*power2(snb)*power2(
     snt) + 2*lmw2*msb12*mst12*mw2*power2(snb)*power2(snt) + 2*lmsb12*msb22*
     mst12*mw2*power2(snb)*power2(snt) - 2*lmw2*msb22*mst12*mw2*power2(snb)*
     power2(snt) + 2*lmsb22*msb12*mst22*mw2*power2(snb)*power2(snt) - 2*lmw2*
     msb12*mst22*mw2*power2(snb)*power2(snt) - 2*lmsb12*msb22*mst22*mw2*power2(
     snb)*power2(snt) + 2*lmw2*msb22*mst22*mw2*power2(snb)*power2(snt)) + mC2*(
     18*msb22*mt2*(-1 + power2(snb))*((lmC2 - lmsb12)*power2(At)*(mst12 - mst12
     *power2(snt) + mst22*power2(snt)) - (lmC2 - lmsb12)*mst22*(mt2*(-1 +
     power2(snt)) + mu2*power2(snt)) + mst12*(mst22 - lmC2*mst22 + (lmC2 -
     lmsb12)*(mu2*(-1 + power2(snt)) + mt2*power2(snt)))) + msb12*(msb22*(mst12
     *mst22*(lmC2*(12 - 60*invdct*mt2) - 35*(1 + invdct*mt2) + 12*lmt2*(1 + 7*
     invdct*mt2)) - 18*(-1 + lmC2)*mst12*mt2*(-1 + power2(snt)) + 18*(-1 + lmC2
     )*mst22*mt2*power2(snt)) + 18*mt2*power2(snb)*(-((lmC2 - lmsb22)*power2(At
     )*(mst12 - mst12*power2(snt) + mst22*power2(snt))) + (lmC2 - lmsb22)*mst22
     *(mt2*(-1 + power2(snt)) + mu2*power2(snt)) + mst12*((-1 + lmC2)*mst22 - (
     lmC2 - lmsb22)*(mu2*(-1 + power2(snt)) + mt2*power2(snt)))))) + msb12*
     msb22*mst12*mst22*(45 + 12*lmt2 + 10*invdct*mt2 + 24*invdct*lmt2*mt2 + 12*
     lmC2*(1 + 2*invdct*mt2))*power2(invdct)*power3(mC2) - 45*msb12*msb22*mst12
     *mst22*power2(invdtw)*power3(mw2) - 12*lmt2*msb12*msb22*mst12*mst22*power2
     (invdtw)*power3(mw2) - 12*lmw2*msb12*msb22*mst12*mst22*power2(invdtw)*
     power3(mw2) + 10*msb12*msb22*mst12*mst22*mt2*power3(invdtw)*power3(mw2) +
     24*lmt2*msb12*msb22*mst12*mst22*mt2*power3(invdtw)*power3(mw2) + 24*lmw2*
     msb12*msb22*mst12*mst22*mt2*power3(invdtw)*power3(mw2) + 2*(5 + 12*lmC2 +
     12*lmt2)*msb12*msb22*mst12*mst22*mt2*power4(invdct)*power4(mC2) - 10*msb12
     *msb22*mst12*mst22*mt2*power4(invdtw)*power4(mw2) - 24*lmt2*msb12*msb22*
     mst12*mst22*mt2*power4(invdtw)*power4(mw2) - 24*lmw2*msb12*msb22*mst12*
     mst22*mt2*power4(invdtw)*power4(mw2) + msb12*msb22*mst12*mst22*(-223 - 213
     *invdct*mt2 + 6*lmC2*(7 + 11*invdct*mt2) + 6*lmt2*(7 + 11*invdct*mt2))*
     power4(invdct)*power5(mC2) + 223*msb12*msb22*mst12*mst22*power4(invdtw)*
     power5(mw2) - 42*lmt2*msb12*msb22*mst12*mst22*power4(invdtw)*power5(mw2) -
     42*lmw2*msb12*msb22*mst12*mst22*power4(invdtw)*power5(mw2) - 213*msb12*
     msb22*mst12*mst22*mt2*power5(invdtw)*power5(mw2) + 66*lmt2*msb12*msb22*
     mst12*mst22*mt2*power5(invdtw)*power5(mw2) + 66*lmw2*msb12*msb22*mst12*
     mst22*mt2*power5(invdtw)*power5(mw2) + 3*(-71 + 22*lmC2 + 22*lmt2)*msb12*
     msb22*mst12*mst22*mt2*power6(invdct)*power6(mC2) + 213*msb12*msb22*mst12*
     mst22*mt2*power6(invdtw)*power6(mw2) - 66*lmt2*msb12*msb22*mst12*mst22*mt2
     *power6(invdtw)*power6(mw2) - 66*lmw2*msb12*msb22*mst12*mst22*mt2*power6(
     invdtw)*power6(mw2) - 3*(-71 + 22*lmC2 + 22*lmt2)*msb12*msb22*mst12*mst22*
     power6(invdct)*power7(mC2) - 213*msb12*msb22*mst12*mst22*power6(invdtw)*
     power7(mw2) + 66*lmt2*msb12*msb22*mst12*mst22*power6(invdtw)*power7(mw2) +
     66*lmw2*msb12*msb22*mst12*mst22*power6(invdtw)*power7(mw2)))/(msb12*msb22*
     mst12*mst22*mt2))/216.;

   return power2(g3) * power2(yt) * result * twoLoop;
}

/**
 * @brief 2-loop \f$O(\alpha_b\alpha_s)\f$ contributions to \f$\Delta g_3^{(2L)}\f$ [arXiv:1009.5455]
 *
 * The function returns \f$\Delta g_3^{(2L)}\f$ which appears in
 * conversion of the strong \f$\overline{\text{MS}}\f$ gauge coupling
 * of the Standard Model with 5 active quark flavours,
 * \f$g_3^{\text{SM(5)},\overline{\text{MS}}}\f$, to the strong
 * \f$\overline{\text{DR}}\f$ gauge coupling of the full MSSM
 * \f$g_3^{\text{MSSM},\overline{\text{DR}}}\f$,
 * \f{align*}{
  g_3^{\text{SM(5)},\overline{\text{MS}}} =
  g_3^{\text{MSSM},\overline{\text{DR}}} \left[
     1 + \Delta g_3^{(1L)} + \Delta g_3^{(2L)}
  \right]
 * \f}
 */
Real delta_alpha_s_2loop_ab_as(const Parameters& pars)
{
   const Real g3        = pars.g3;
   const Real yb        = pars.yb;
   const Real xt        = pars.xt;
   const Real xb        = pars.xb;
   const Real mt        = pars.mt;
   const Real mt2       = power2(pars.mt);
   const Real mb        = pars.mb;
   const Real mst12     = power2(pars.mst1);
   const Real mst14     = power2(mst12);
   const Real mst16     = mst12*mst14;
   const Real mst22     = power2(pars.mst2);
   const Real mst24     = power2(mst22);
   const Real mst26     = mst22*mst24;
   const Real msb12     = power2(pars.msb1);
   const Real msb14     = power2(msb12);
   const Real msb16     = msb12*msb14;
   const Real msb22     = power2(pars.msb2);
   const Real msb24     = power2(msb22);
   const Real msb26     = msb22*msb24;
   const Real mw2       = power2(pars.mw);
   const Real mz2       = power2(pars.mz);
   const Real mh2       = power2(pars.mh);
   const Real mH2       = power2(pars.mH);
   const Real mC2       = power2(pars.mC);
   const Real mA2       = power2(pars.mA);
   const Real mu        = pars.mu;
   const Real mu2       = power2(pars.mu);
   const Real tb        = pars.tb;
   const Real sb        = tb / std::sqrt(1. + power2(tb));
   const Real cb        = 1. / std::sqrt(1. + power2(tb));
   const Real Q2        = power2(pars.Q);
   const Real snt       = calc_sin_theta(mt, xt, mst12, mst22);
   const Real snb       = calc_sin_theta(mb, xb, msb12, msb22);
   const Real alpha     = calc_alpha(mh2, mH2, tb);
   const Real sa        = std::sin(alpha);
   const Real ca        = std::cos(alpha);
   const Real Ab        = xb + mu*sb/cb;
   const Real invdmst   = 1/(mst12 - mst22);
   const Real invdct    = 1/(mC2 - mt2);
   const Real invdtw    = 1/(mt2 - mw2);
   const Real lmst12    = std::log(mst12/Q2);
   const Real lmst22    = std::log(mst22/Q2);
   const Real lmsb12    = std::log(msb12/Q2);
   const Real lmsb22    = std::log(msb22/Q2);
   const Real lmt2      = std::log(mt2/Q2);
   const Real lmw2      = std::log(mw2/Q2);
   const Real lmz2      = std::log(mz2/Q2);
   const Real lmu2      = std::log(mu2/Q2);
   const Real lmH2      = std::log(mH2/Q2);
   const Real lmA2      = std::log(mA2/Q2);
   const Real lmh2      = std::log(mh2/Q2);
   const Real lmC2      = std::log(mC2/Q2);

   const Real result =
   (mu2*DeltaInv(mt2,mu2,msb22)*(-(msb22*(msb24 + mu2*(-3*mt2 + 2*lmt2*mt2 - 2*
     lmu2*mt2 + mu2 + 4*lmt2*mu2 - 4*lmu2*mu2))) + (lmt2 - lmu2)*(mt2 - mu2)*
     power2(mu2) + power2(msb22)*((1 + lmsb22 - lmt2)*mt2 - mu2*(40 - 17*lmt2 +
     6*lmsb22*(-2 + lmt2 - lmu2) - 7*lmu2 + 6*lmt2*lmu2 + 6*power2(lmt2) +
     power2(Pi))))*(-1 + power2(snb)))/(6.*msb22*mt2) + (mu2*DeltaInv(msb12,mt2
     ,mu2)*(msb12*(msb14 + mu2*(-3*mt2 + 2*lmt2*mt2 - 2*lmu2*mt2 + mu2 + 4*lmt2
     *mu2 - 4*lmu2*mu2)) - (lmt2 - lmu2)*(mt2 - mu2)*power2(mu2) + power2(msb12
     )*((-1 - lmsb12 + lmt2)*mt2 + mu2*(40 - 17*lmt2 + 6*lmsb12*(-2 + lmt2 -
     lmu2) - 7*lmu2 + 6*lmt2*lmu2 + 6*power2(lmt2) + power2(Pi))))*power2(snb))
     /(6.*msb12*mt2) - (((1 + lmh2 - lmsb12)*mh2 + 6*(-lmh2 + lmsb12)*msb12)*
     DeltaInv(mh2,msb12,msb12)*(mu2 + Ab*sa*(2*ca*mu + Ab*sa) - mu2*power2(sa))
     *(-1 + power2(snb))*power2(snb))/6. - (((1 + lmh2 - lmsb22)*mh2 + 6*(-lmh2
      + lmsb22)*msb22)*DeltaInv(mh2,msb22,msb22)*(mu2 + Ab*sa*(2*ca*mu + Ab*sa)
     - mu2*power2(sa))*(-1 + power2(snb))*power2(snb))/6. + (((1 + lmH2 -
     lmsb12)*mH2 + 6*(-lmH2 + lmsb12)*msb12)*DeltaInv(mH2,msb12,msb12)*(2*Ab*ca
     *mu*sa + power2(Ab)*(-1 + power2(sa)) - mu2*power2(sa))*(-1 + power2(snb))
     *power2(snb))/6. + (((1 + lmH2 - lmsb22)*mH2 + 6*(-lmH2 + lmsb22)*msb22)*
     DeltaInv(mH2,msb22,msb22)*(2*Ab*ca*mu*sa + power2(Ab)*(-1 + power2(sa)) -
     mu2*power2(sa))*(-1 + power2(snb))*power2(snb))/6. + (invdmst*(mst12 -
     mst22)*mt2*xt*(-(cb*mu*sb) + Ab*(-1 + power2(cb)))*(-(lmsb22*msb12*(mC2 +
     msb22)*(-1 + power2(snb))) + msb22*(-msb12 + lmsb12*mC2*power2(snb) +
     lmsb12*msb12*power2(snb)) - lmC2*mC2*(msb12 - msb12*power2(snb) + msb22*
     power2(snb))))/(6.*msb12*msb22*mst12*mst22) + ((6*msb16*(42 + 8*lmh2*(-3 +
     lmsb12) - 12*lmsb12 + 4*power2(lmh2) + power2(Pi)) + mh2*msb14*(42 + 12*
     lmh2*(-1 + lmsb12) - 24*lmsb12 + 6*power2(lmh2) - 6*power2(lmsb12) +
     power2(Pi)))*(mu2 + Ab*sa*(2*ca*mu + Ab*sa) - mu2*power2(sa))*(-1 + power2
     (snb))*power2(snb)*power2(DeltaInv(mh2,msb12,msb12)))/6. + ((6*msb26*(42 +
     8*lmh2*(-3 + lmsb22) - 12*lmsb22 + 4*power2(lmh2) + power2(Pi)) + mh2*
     msb24*(42 + 12*lmh2*(-1 + lmsb22) - 24*lmsb22 + 6*power2(lmh2) - 6*power2(
     lmsb22) + power2(Pi)))*(mu2 + Ab*sa*(2*ca*mu + Ab*sa) - mu2*power2(sa))*(-
     1 + power2(snb))*power2(snb)*power2(DeltaInv(mh2,msb22,msb22)))/6. - ((6*
     msb16*(42 + 8*lmH2*(-3 + lmsb12) - 12*lmsb12 + 4*power2(lmH2) + power2(Pi)
     ) + mH2*msb14*(42 + 12*lmH2*(-1 + lmsb12) - 24*lmsb12 + 6*power2(lmH2) - 6
     *power2(lmsb12) + power2(Pi)))*(2*Ab*ca*mu*sa + power2(Ab)*(-1 + power2(sa
     )) - mu2*power2(sa))*(-1 + power2(snb))*power2(snb)*power2(DeltaInv(mH2,
     msb12,msb12)))/6. - ((6*msb26*(42 + 8*lmH2*(-3 + lmsb22) - 12*lmsb22 + 4*
     power2(lmH2) + power2(Pi)) + mH2*msb24*(42 + 12*lmH2*(-1 + lmsb22) - 24*
     lmsb22 + 6*power2(lmH2) - 6*power2(lmsb22) + power2(Pi)))*(2*Ab*ca*mu*sa +
     power2(Ab)*(-1 + power2(sa)) - mu2*power2(sa))*(-1 + power2(snb))*power2(
     snb)*power2(DeltaInv(mH2,msb22,msb22)))/6. + (power2(mu2)*(msb16*(42 + 6*
     lmu2 - 6*lmt2*(3 + lmu2) + 6*lmsb12*(-4 + lmt2 + lmu2) + 6*power2(lmsb12)
     + power2(Pi)) + mu2*(6*(lmt2 - lmu2)*(mt2 - mu2)*mu2 + 6*msb12*mt2*(42 + 4
     *lmsb12*(-1 + lmt2 - lmu2) - 12*lmu2 + lmt2*(-20 + 6*lmu2) + 5*power2(lmt2
     ) + power2(lmu2) + power2(Pi)) + msb12*mu2*(294 + 6*lmsb12*(2 + lmt2 -
     lmu2) - 114*lmu2 + 6*lmt2*(-25 + 7*lmu2) + 24*power2(lmt2) + 18*power2(
     lmu2) + 7*power2(Pi))) + 2*msb14*(mt2*(84 - 39*lmt2 + 6*lmsb12*(-5 + 2*
     lmt2 - lmu2) - 3*lmu2 + 6*lmt2*lmu2 + 3*power2(lmsb12) + 9*power2(lmt2) +
     2*power2(Pi)) + mu2*(336 - 201*lmt2 + 6*lmsb12*(-11 + 7*lmt2 - 4*lmu2) -
     21*lmu2 + 30*lmt2*lmu2 + 9*power2(lmsb12) + 36*power2(lmt2) + 3*power2(
     lmu2) + 8*power2(Pi))))*power2(snb)*power2(DeltaInv(msb12,mt2,mu2)))/(6.*
     mt2) - (mz2*(mu2 - 2*Ab*cb*mu*sb - mu2*power2(cb) + power2(Ab)*power2(cb))
     *(-3*lmsb22*(-msb26 + mz2*(2*(7 + lmz2)*msb24 + (-1 + 2*lmz2)*msb22*mz2) +
     msb12*(msb24 + (3 + 2*lmz2)*msb22*mz2)) + 3*lmsb12*(-msb26 + mz2*((-23 + 6
     *lmsb22 + 2*lmz2)*msb24 + ((-23 + 2*lmsb22 + 6*lmz2)*msb22 - mz2)*mz2) +
     msb12*(msb24 + mz2*(2*(-3 + lmsb22 + lmz2)*msb22 + mz2))) + 6*mz2*(msb12*
     msb22 + 2*(msb24 + msb22*mz2))*power2(lmsb12) + 6*msb24*mz2*power2(lmsb22)
     + mz2*(-3*lmz2*msb12*mz2 + msb12*msb22*(42 - 9*lmz2 + power2(Pi)) + 3*
     msb24*(42 + lmz2 + power2(Pi)) + 3*mz2*(lmz2*mz2 + msb22*(42 - 14*lmz2 + 2
     *power2(lmz2) + power2(Pi)))))*power2(DeltaInv(msb12,mz2,msb22)))/(12.*
     msb12) - (mw2*(mu2 - 2*Ab*cb*mu*sb - mu2*power2(cb) + power2(Ab)*power2(cb
     ))*(-3*lmst22*(-mst26 + mw2*(2*(7 + lmw2)*mst24 + (-1 + 2*lmw2)*mst22*mw2)
     + msb12*(mst24 + (3 + 2*lmw2)*mst22*mw2)) + 3*lmsb12*(-mst26 + mw2*((-23 +
     6*lmst22 + 2*lmw2)*mst24 + ((-23 + 2*lmst22 + 6*lmw2)*mst22 - mw2)*mw2) +
     msb12*(mst24 + mw2*(2*(-3 + lmst22 + lmw2)*mst22 + mw2))) + 6*mw2*(msb12*
     mst22 + 2*(mst24 + mst22*mw2))*power2(lmsb12) + 6*mst24*mw2*power2(lmst22)
     + mw2*(-3*lmw2*msb12*mw2 + msb12*mst22*(42 - 9*lmw2 + power2(Pi)) + 3*
     mst24*(42 + lmw2 + power2(Pi)) + 3*mw2*(lmw2*mw2 + mst22*(42 - 14*lmw2 + 2
     *power2(lmw2) + power2(Pi)))))*power2(snb)*power2(snt)*power2(DeltaInv(
     mst22,mw2,msb12)))/(6.*msb12) + (mw2*(mu2 - 2*Ab*cb*mu*sb - mu2*power2(cb)
     + power2(Ab)*power2(cb))*(-3*lmst22*(-mst26 + mw2*(2*(7 + lmw2)*mst24 + (-
     1 + 2*lmw2)*mst22*mw2) + msb22*(mst24 + (3 + 2*lmw2)*mst22*mw2)) + 3*
     lmsb22*(-mst26 + mw2*((-23 + 6*lmst22 + 2*lmw2)*mst24 + ((-23 + 2*lmst22 +
     6*lmw2)*mst22 - mw2)*mw2) + msb22*(mst24 + mw2*(2*(-3 + lmst22 + lmw2)*
     mst22 + mw2))) + 6*mw2*(msb22*mst22 + 2*(mst24 + mst22*mw2))*power2(lmsb22
     ) + 6*mst24*mw2*power2(lmst22) + mw2*(-3*lmw2*msb22*mw2 + msb22*mst22*(42
     - 9*lmw2 + power2(Pi)) + 3*mst24*(42 + lmw2 + power2(Pi)) + 3*mw2*(lmw2*
     mw2 + mst22*(42 - 14*lmw2 + 2*power2(lmw2) + power2(Pi)))))*(-1 + power2(
     snb))*power2(snt)*power2(DeltaInv(mst22,mw2,msb22)))/(6.*msb22) - (power2(
     mu2)*(msb26*(42 + 6*lmu2 - 6*lmt2*(3 + lmu2) + 6*lmsb22*(-4 + lmt2 + lmu2)
     + 6*power2(lmsb22) + power2(Pi)) + mu2*(6*(lmt2 - lmu2)*(mt2 - mu2)*mu2 +
     6*msb22*mt2*(42 + 4*lmsb22*(-1 + lmt2 - lmu2) - 12*lmu2 + lmt2*(-20 + 6*
     lmu2) + 5*power2(lmt2) + power2(lmu2) + power2(Pi)) + msb22*mu2*(294 + 6*
     lmsb22*(2 + lmt2 - lmu2) - 114*lmu2 + 6*lmt2*(-25 + 7*lmu2) + 24*power2(
     lmt2) + 18*power2(lmu2) + 7*power2(Pi))) + 2*msb24*(mt2*(84 - 39*lmt2 + 6*
     lmsb22*(-5 + 2*lmt2 - lmu2) - 3*lmu2 + 6*lmt2*lmu2 + 3*power2(lmsb22) + 9*
     power2(lmt2) + 2*power2(Pi)) + mu2*(336 - 201*lmt2 + 6*lmsb22*(-11 + 7*
     lmt2 - 4*lmu2) - 21*lmu2 + 30*lmt2*lmu2 + 9*power2(lmsb22) + 36*power2(
     lmt2) + 3*power2(lmu2) + 8*power2(Pi))))*(-1 + power2(snb))*power2(
     DeltaInv(mt2,mu2,msb22)))/(6.*mt2) + (mw2*(mu2 - 2*Ab*cb*mu*sb - mu2*
     power2(cb) + power2(Ab)*power2(cb))*(-3*lmst12*(-mst16 + mw2*(2*(7 + lmw2)
     *mst14 + (-1 + 2*lmw2)*mst12*mw2) + msb12*(mst14 + (3 + 2*lmw2)*mst12*mw2)
     ) + 3*lmsb12*(-mst16 + mw2*((-23 + 6*lmst12 + 2*lmw2)*mst14 + ((-23 + 2*
     lmst12 + 6*lmw2)*mst12 - mw2)*mw2) + msb12*(mst14 + mw2*(2*(-3 + lmst12 +
     lmw2)*mst12 + mw2))) + 6*mw2*(msb12*mst12 + 2*(mst14 + mst12*mw2))*power2(
     lmsb12) + 6*mst14*mw2*power2(lmst12) + mw2*(-3*lmw2*msb12*mw2 + msb12*
     mst12*(42 - 9*lmw2 + power2(Pi)) + 3*mst14*(42 + lmw2 + power2(Pi)) + 3*
     mw2*(lmw2*mw2 + mst12*(42 - 14*lmw2 + 2*power2(lmw2) + power2(Pi)))))*
     power2(snb)*(-1 + power2(snt))*power2(DeltaInv(mw2,msb12,mst12)))/(6.*
     msb12) - (mw2*(mu2 - 2*Ab*cb*mu*sb - mu2*power2(cb) + power2(Ab)*power2(cb
     ))*(-3*lmst12*(-mst16 + mw2*(2*(7 + lmw2)*mst14 + (-1 + 2*lmw2)*mst12*mw2)
     + msb22*(mst14 + (3 + 2*lmw2)*mst12*mw2)) + 3*lmsb22*(-mst16 + mw2*((-23 +
     6*lmst12 + 2*lmw2)*mst14 + ((-23 + 2*lmst12 + 6*lmw2)*mst12 - mw2)*mw2) +
     msb22*(mst14 + mw2*(2*(-3 + lmst12 + lmw2)*mst12 + mw2))) + 6*mw2*(msb22*
     mst12 + 2*(mst14 + mst12*mw2))*power2(lmsb22) + 6*mst14*mw2*power2(lmst12)
     + mw2*(-3*lmw2*msb22*mw2 + msb22*mst12*(42 - 9*lmw2 + power2(Pi)) + 3*
     mst14*(42 + lmw2 + power2(Pi)) + 3*mw2*(lmw2*mw2 + mst12*(42 - 14*lmw2 + 2
     *power2(lmw2) + power2(Pi)))))*(-1 + power2(snb))*(-1 + power2(snt))*
     power2(DeltaInv(mw2,msb22,mst12)))/(6.*msb22) + (mA2*(-2*Ab*cb*mu*sb +
     power2(Ab)*(-1 + power2(cb)) - mu2*power2(cb))*power2(DeltaInv(msb12,mA2,
     msb22))*(42*mA2*msb12*msb22 - 9*lmsb22*mA2*msb12*msb22 + 126*mA2*msb24 -
     42*lmsb22*mA2*msb24 - 3*lmsb22*msb12*msb24 + 3*lmsb22*msb26 + 6*mA2*(2*mA2
     *msb22 + msb12*msb22 + 2*msb24)*power2(lmsb12) + 6*mA2*msb24*power2(lmsb22
     ) + 126*msb22*power2(mA2) + 3*lmsb22*msb22*power2(mA2) + 6*msb22*power2(
     lmA2)*power2(mA2) + 3*lmA2*mA2*((-3 + 2*lmsb12 - 2*lmsb22)*msb12*msb22 -
     mA2*(msb12 + 2*(7 - 3*lmsb12 + lmsb22)*msb22) + (1 + 2*lmsb12 - 2*lmsb22)*
     msb24 + power2(mA2)) + mA2*msb12*msb22*power2(Pi) + 3*mA2*msb24*power2(Pi)
     + 3*msb22*power2(mA2)*power2(Pi) - 3*lmsb12*(-(msb12*msb24) + mA2*(-2*(-3
     + lmsb22)*msb12*msb22 + (23 - 6*lmsb22)*msb24) + msb26 - (msb12 + (-23 + 2
     *lmsb22)*msb22)*power2(mA2) + power3(mA2))))/(12.*msb12) - (DeltaInv(msb12
     ,mA2,msb22)*(-2*Ab*cb*mu*sb + power2(Ab)*(-1 + power2(cb)) - mu2*power2(cb
     ))*(msb22*(2*mA2*msb12 + 5*lmsb22*mA2*msb22 - lmsb22*msb12*msb22 + lmsb22*
     msb24) + lmA2*(mA2 - msb12 + 5*msb22)*power2(mA2) - lmsb12*(msb22*(-(msb12
     *msb22) + msb24) - (msb12 - 5*msb22)*power2(mA2) + 5*mA2*power2(msb22) +
     power3(mA2))))/(24.*msb12*msb22) + power2(DeltaInv(msb12,mst12,mC2))*((
     invdmst*mC2*mt2*xt*(-(cb*mu*sb) + Ab*(-1 + power2(cb)))*power2(snb)*(42*
     mC2*msb12*mst12 - 9*lmst12*mC2*msb12*mst12 + 126*mC2*mst14 - 42*lmst12*mC2
     *mst14 - 3*lmst12*msb12*mst14 + 3*lmst12*mst16 + 6*mC2*(2*mC2*mst12 +
     msb12*mst12 + 2*mst14)*power2(lmsb12) + 6*mC2*mst14*power2(lmst12) + 126*
     mst12*power2(mC2) + 3*lmst12*mst12*power2(mC2) + 6*mst12*power2(lmC2)*
     power2(mC2) + 3*lmC2*mC2*((-3 + 2*lmsb12 - 2*lmst12)*msb12*mst12 - mC2*(
     msb12 + 2*(7 - 3*lmsb12 + lmst12)*mst12) + (1 + 2*lmsb12 - 2*lmst12)*mst14
      + power2(mC2)) + mC2*msb12*mst12*power2(Pi) + 3*mC2*mst14*power2(Pi) + 3*
     mst12*power2(mC2)*power2(Pi) - 3*lmsb12*(-(msb12*mst14) + mC2*(-2*(-3 +
     lmst12)*msb12*mst12 + (23 - 6*lmst12)*mst14) + mst16 - (msb12 + (-23 + 2*
     lmst12)*mst12)*power2(mC2) + power3(mC2))))/(3.*msb12) - (mC2*power2(snb)*
     (-2*Ab*cb*mu*sb*(-1 + power2(snt)) + power2(Ab)*(-1 + power2(cb))*(-1 +
     power2(snt)) + mt2*power2(snt) + power2(cb)*(mu2 - mt2*power2(snt) - mu2*
     power2(snt)))*(42*mC2*msb12*mst12 - 9*lmst12*mC2*msb12*mst12 + 126*mC2*
     mst14 - 42*lmst12*mC2*mst14 - 3*lmst12*msb12*mst14 + 3*lmst12*mst16 + 6*
     mC2*(2*mC2*mst12 + msb12*mst12 + 2*mst14)*power2(lmsb12) + 6*mC2*mst14*
     power2(lmst12) + 126*mst12*power2(mC2) + 3*lmst12*mst12*power2(mC2) + 6*
     mst12*power2(lmC2)*power2(mC2) + 3*lmC2*mC2*((-3 + 2*lmsb12 - 2*lmst12)*
     msb12*mst12 - mC2*(msb12 + 2*(7 - 3*lmsb12 + lmst12)*mst12) + (1 + 2*
     lmsb12 - 2*lmst12)*mst14 + power2(mC2)) + mC2*msb12*mst12*power2(Pi) + 3*
     mC2*mst14*power2(Pi) + 3*mst12*power2(mC2)*power2(Pi) - 3*lmsb12*(-(msb12*
     mst14) + mC2*(-2*(-3 + lmst12)*msb12*mst12 + (23 - 6*lmst12)*mst14) +
     mst16 - (msb12 + (-23 + 2*lmst12)*mst12)*power2(mC2) + power3(mC2))))/(6.*
     msb12)) + power2(DeltaInv(msb22,mst12,mC2))*(-(invdmst*mC2*mt2*xt*(-(cb*mu
     *sb) + Ab*(-1 + power2(cb)))*(-1 + power2(snb))*(42*mC2*msb22*mst12 - 9*
     lmst12*mC2*msb22*mst12 + 126*mC2*mst14 - 42*lmst12*mC2*mst14 - 3*lmst12*
     msb22*mst14 + 3*lmst12*mst16 + 6*mC2*(2*mC2*mst12 + msb22*mst12 + 2*mst14)
     *power2(lmsb22) + 6*mC2*mst14*power2(lmst12) + 126*mst12*power2(mC2) + 3*
     lmst12*mst12*power2(mC2) + 6*mst12*power2(lmC2)*power2(mC2) + 3*lmC2*mC2*(
     (-3 + 2*lmsb22 - 2*lmst12)*msb22*mst12 - mC2*(msb22 + 2*(7 - 3*lmsb22 +
     lmst12)*mst12) + (1 + 2*lmsb22 - 2*lmst12)*mst14 + power2(mC2)) + mC2*
     msb22*mst12*power2(Pi) + 3*mC2*mst14*power2(Pi) + 3*mst12*power2(mC2)*
     power2(Pi) - 3*lmsb22*(-(msb22*mst14) + mC2*(-2*(-3 + lmst12)*msb22*mst12
     + (23 - 6*lmst12)*mst14) + mst16 - (msb22 + (-23 + 2*lmst12)*mst12)*power2
     (mC2) + power3(mC2))))/(3.*msb22) + (mC2*(-1 + power2(snb))*(-2*Ab*cb*mu*
     sb*(-1 + power2(snt)) + power2(Ab)*(-1 + power2(cb))*(-1 + power2(snt)) +
     mt2*power2(snt) + power2(cb)*(mu2 - mt2*power2(snt) - mu2*power2(snt)))*(
     42*mC2*msb22*mst12 - 9*lmst12*mC2*msb22*mst12 + 126*mC2*mst14 - 42*lmst12*
     mC2*mst14 - 3*lmst12*msb22*mst14 + 3*lmst12*mst16 + 6*mC2*(2*mC2*mst12 +
     msb22*mst12 + 2*mst14)*power2(lmsb22) + 6*mC2*mst14*power2(lmst12) + 126*
     mst12*power2(mC2) + 3*lmst12*mst12*power2(mC2) + 6*mst12*power2(lmC2)*
     power2(mC2) + 3*lmC2*mC2*((-3 + 2*lmsb22 - 2*lmst12)*msb22*mst12 - mC2*(
     msb22 + 2*(7 - 3*lmsb22 + lmst12)*mst12) + (1 + 2*lmsb22 - 2*lmst12)*mst14
      + power2(mC2)) + mC2*msb22*mst12*power2(Pi) + 3*mC2*mst14*power2(Pi) + 3*
     mst12*power2(mC2)*power2(Pi) - 3*lmsb22*(-(msb22*mst14) + mC2*(-2*(-3 +
     lmst12)*msb22*mst12 + (23 - 6*lmst12)*mst14) + mst16 - (msb22 + (-23 + 2*
     lmst12)*mst12)*power2(mC2) + power3(mC2))))/(6.*msb22)) + power2(DeltaInv(
     mst22,msb12,mC2))*(-(invdmst*mC2*mt2*xt*(-(cb*mu*sb) + Ab*(-1 + power2(cb)
     ))*power2(snb)*(42*mC2*msb12*mst22 - 9*lmst22*mC2*msb12*mst22 + 126*mC2*
     mst24 - 42*lmst22*mC2*mst24 - 3*lmst22*msb12*mst24 + 3*lmst22*mst26 + 6*
     mC2*(2*mC2*mst22 + msb12*mst22 + 2*mst24)*power2(lmsb12) + 6*mC2*mst24*
     power2(lmst22) + 126*mst22*power2(mC2) + 3*lmst22*mst22*power2(mC2) + 6*
     mst22*power2(lmC2)*power2(mC2) + 3*lmC2*mC2*((-3 + 2*lmsb12 - 2*lmst22)*
     msb12*mst22 - mC2*(msb12 + 2*(7 - 3*lmsb12 + lmst22)*mst22) + (1 + 2*
     lmsb12 - 2*lmst22)*mst24 + power2(mC2)) + mC2*msb12*mst22*power2(Pi) + 3*
     mC2*mst24*power2(Pi) + 3*mst22*power2(mC2)*power2(Pi) - 3*lmsb12*(-(msb12*
     mst24) + mC2*(-2*(-3 + lmst22)*msb12*mst22 + (23 - 6*lmst22)*mst24) +
     mst26 - (msb12 + (-23 + 2*lmst22)*mst22)*power2(mC2) + power3(mC2))))/(3.*
     msb12) - (mC2*power2(snb)*(mt2*(-1 + power2(cb))*(-1 + power2(snt)) + (2*
     Ab*cb*mu*sb - power2(Ab)*(-1 + power2(cb)) + mu2*power2(cb))*power2(snt))*
     (42*mC2*msb12*mst22 - 9*lmst22*mC2*msb12*mst22 + 126*mC2*mst24 - 42*lmst22
     *mC2*mst24 - 3*lmst22*msb12*mst24 + 3*lmst22*mst26 + 6*mC2*(2*mC2*mst22 +
     msb12*mst22 + 2*mst24)*power2(lmsb12) + 6*mC2*mst24*power2(lmst22) + 126*
     mst22*power2(mC2) + 3*lmst22*mst22*power2(mC2) + 6*mst22*power2(lmC2)*
     power2(mC2) + 3*lmC2*mC2*((-3 + 2*lmsb12 - 2*lmst22)*msb12*mst22 - mC2*(
     msb12 + 2*(7 - 3*lmsb12 + lmst22)*mst22) + (1 + 2*lmsb12 - 2*lmst22)*mst24
      + power2(mC2)) + mC2*msb12*mst22*power2(Pi) + 3*mC2*mst24*power2(Pi) + 3*
     mst22*power2(mC2)*power2(Pi) - 3*lmsb12*(-(msb12*mst24) + mC2*(-2*(-3 +
     lmst22)*msb12*mst22 + (23 - 6*lmst22)*mst24) + mst26 - (msb12 + (-23 + 2*
     lmst22)*mst22)*power2(mC2) + power3(mC2))))/(6.*msb12)) + power2(DeltaInv(
     mst22,msb22,mC2))*((invdmst*mC2*mt2*xt*(-(cb*mu*sb) + Ab*(-1 + power2(cb))
     )*(-1 + power2(snb))*(42*mC2*msb22*mst22 - 9*lmst22*mC2*msb22*mst22 + 126*
     mC2*mst24 - 42*lmst22*mC2*mst24 - 3*lmst22*msb22*mst24 + 3*lmst22*mst26 +
     6*mC2*(2*mC2*mst22 + msb22*mst22 + 2*mst24)*power2(lmsb22) + 6*mC2*mst24*
     power2(lmst22) + 126*mst22*power2(mC2) + 3*lmst22*mst22*power2(mC2) + 6*
     mst22*power2(lmC2)*power2(mC2) + 3*lmC2*mC2*((-3 + 2*lmsb22 - 2*lmst22)*
     msb22*mst22 - mC2*(msb22 + 2*(7 - 3*lmsb22 + lmst22)*mst22) + (1 + 2*
     lmsb22 - 2*lmst22)*mst24 + power2(mC2)) + mC2*msb22*mst22*power2(Pi) + 3*
     mC2*mst24*power2(Pi) + 3*mst22*power2(mC2)*power2(Pi) - 3*lmsb22*(-(msb22*
     mst24) + mC2*(-2*(-3 + lmst22)*msb22*mst22 + (23 - 6*lmst22)*mst24) +
     mst26 - (msb22 + (-23 + 2*lmst22)*mst22)*power2(mC2) + power3(mC2))))/(3.*
     msb22) + (mC2*(-1 + power2(snb))*(mt2*(-1 + power2(cb))*(-1 + power2(snt))
     + (2*Ab*cb*mu*sb - power2(Ab)*(-1 + power2(cb)) + mu2*power2(cb))*power2(
     snt))*(42*mC2*msb22*mst22 - 9*lmst22*mC2*msb22*mst22 + 126*mC2*mst24 - 42*
     lmst22*mC2*mst24 - 3*lmst22*msb22*mst24 + 3*lmst22*mst26 + 6*mC2*(2*mC2*
     mst22 + msb22*mst22 + 2*mst24)*power2(lmsb22) + 6*mC2*mst24*power2(lmst22)
     + 126*mst22*power2(mC2) + 3*lmst22*mst22*power2(mC2) + 6*mst22*power2(lmC2
     )*power2(mC2) + 3*lmC2*mC2*((-3 + 2*lmsb22 - 2*lmst22)*msb22*mst22 - mC2*(
     msb22 + 2*(7 - 3*lmsb22 + lmst22)*mst22) + (1 + 2*lmsb22 - 2*lmst22)*mst24
      + power2(mC2)) + mC2*msb22*mst22*power2(Pi) + 3*mC2*mst24*power2(Pi) + 3*
     mst22*power2(mC2)*power2(Pi) - 3*lmsb22*(-(msb22*mst24) + mC2*(-2*(-3 +
     lmst22)*msb22*mst22 + (23 - 6*lmst22)*mst24) + mst26 - (msb22 + (-23 + 2*
     lmst22)*mst22)*power2(mC2) + power3(mC2))))/(6.*msb22)) + DeltaInv(msb12,
     mst12,mC2)*(-(invdmst*mt2*xt*(-(cb*mu*sb) + Ab*(-1 + power2(cb)))*power2(
     snb)*(mst12*(2*mC2*msb12 + 5*lmst12*mC2*mst12 - lmst12*msb12*mst12 +
     lmst12*mst14) + lmC2*(mC2 - msb12 + 5*mst12)*power2(mC2) - lmsb12*(mst12*(
     -(msb12*mst12) + mst14) - (msb12 - 5*mst12)*power2(mC2) + 5*mC2*power2(
     mst12) + power3(mC2))))/(6.*msb12*mst12) + (power2(snb)*(-2*Ab*cb*mu*sb*(-
     1 + power2(snt)) + power2(Ab)*(-1 + power2(cb))*(-1 + power2(snt)) + mt2*
     power2(snt) + power2(cb)*(mu2 - mt2*power2(snt) - mu2*power2(snt)))*(mst12
     *(2*mC2*msb12 + 5*lmst12*mC2*mst12 - lmst12*msb12*mst12 + lmst12*mst14) +
     lmC2*(mC2 - msb12 + 5*mst12)*power2(mC2) - lmsb12*(mst12*(-(msb12*mst12) +
     mst14) - (msb12 - 5*mst12)*power2(mC2) + 5*mC2*power2(mst12) + power3(mC2)
     )))/(12.*msb12*mst12)) + DeltaInv(msb22,mst12,mC2)*((invdmst*mt2*xt*(-(cb*
     mu*sb) + Ab*(-1 + power2(cb)))*(-1 + power2(snb))*(mst12*(2*mC2*msb22 + 5*
     lmst12*mC2*mst12 - lmst12*msb22*mst12 + lmst12*mst14) + lmC2*(mC2 - msb22
     + 5*mst12)*power2(mC2) - lmsb22*(mst12*(-(msb22*mst12) + mst14) - (msb22 -
     5*mst12)*power2(mC2) + 5*mC2*power2(mst12) + power3(mC2))))/(6.*msb22*
     mst12) - ((-1 + power2(snb))*(-2*Ab*cb*mu*sb*(-1 + power2(snt)) + power2(
     Ab)*(-1 + power2(cb))*(-1 + power2(snt)) + mt2*power2(snt) + power2(cb)*(
     mu2 - mt2*power2(snt) - mu2*power2(snt)))*(mst12*(2*mC2*msb22 + 5*lmst12*
     mC2*mst12 - lmst12*msb22*mst12 + lmst12*mst14) + lmC2*(mC2 - msb22 + 5*
     mst12)*power2(mC2) - lmsb22*(mst12*(-(msb22*mst12) + mst14) - (msb22 - 5*
     mst12)*power2(mC2) + 5*mC2*power2(mst12) + power3(mC2))))/(12.*msb22*mst12
     )) + DeltaInv(mst22,msb12,mC2)*((invdmst*mt2*xt*(-(cb*mu*sb) + Ab*(-1 +
     power2(cb)))*power2(snb)*(mst22*(2*mC2*msb12 + 5*lmst22*mC2*mst22 - lmst22
     *msb12*mst22 + lmst22*mst24) + lmC2*(mC2 - msb12 + 5*mst22)*power2(mC2) -
     lmsb12*(mst22*(-(msb12*mst22) + mst24) - (msb12 - 5*mst22)*power2(mC2) + 5
     *mC2*power2(mst22) + power3(mC2))))/(6.*msb12*mst22) + (power2(snb)*(mt2*(
     -1 + power2(cb))*(-1 + power2(snt)) + (2*Ab*cb*mu*sb - power2(Ab)*(-1 +
     power2(cb)) + mu2*power2(cb))*power2(snt))*(mst22*(2*mC2*msb12 + 5*lmst22*
     mC2*mst22 - lmst22*msb12*mst22 + lmst22*mst24) + lmC2*(mC2 - msb12 + 5*
     mst22)*power2(mC2) - lmsb12*(mst22*(-(msb12*mst22) + mst24) - (msb12 - 5*
     mst22)*power2(mC2) + 5*mC2*power2(mst22) + power3(mC2))))/(12.*msb12*mst22
     )) + DeltaInv(mst22,msb22,mC2)*(-(invdmst*mt2*xt*(-(cb*mu*sb) + Ab*(-1 +
     power2(cb)))*(-1 + power2(snb))*(mst22*(2*mC2*msb22 + 5*lmst22*mC2*mst22 -
     lmst22*msb22*mst22 + lmst22*mst24) + lmC2*(mC2 - msb22 + 5*mst22)*power2(
     mC2) - lmsb22*(mst22*(-(msb22*mst22) + mst24) - (msb22 - 5*mst22)*power2(
     mC2) + 5*mC2*power2(mst22) + power3(mC2))))/(6.*msb22*mst22) - ((-1 +
     power2(snb))*(mt2*(-1 + power2(cb))*(-1 + power2(snt)) + (2*Ab*cb*mu*sb -
     power2(Ab)*(-1 + power2(cb)) + mu2*power2(cb))*power2(snt))*(mst22*(2*mC2*
     msb22 + 5*lmst22*mC2*mst22 - lmst22*msb22*mst22 + lmst22*mst24) + lmC2*(
     mC2 - msb22 + 5*mst22)*power2(mC2) - lmsb22*(mst22*(-(msb22*mst22) + mst24
     ) - (msb22 - 5*mst22)*power2(mC2) + 5*mC2*power2(mst22) + power3(mC2))))/(
     12.*msb22*mst22)) - (mh2*(mu2 + Ab*sa*(2*ca*mu + Ab*sa) - mu2*power2(sa))*
     power2(DeltaInv(mh2,msb12,msb22))*power2(1 - 2*power2(snb))*(42*mh2*msb12*
     msb22 - 9*lmsb22*mh2*msb12*msb22 + 126*mh2*msb24 - 42*lmsb22*mh2*msb24 - 3
     *lmsb22*msb12*msb24 + 3*lmsb22*msb26 + 6*mh2*(2*mh2*msb22 + msb12*msb22 +
     2*msb24)*power2(lmsb12) + 6*mh2*msb24*power2(lmsb22) + 126*msb22*power2(
     mh2) + 3*lmsb22*msb22*power2(mh2) + 6*msb22*power2(lmh2)*power2(mh2) + 3*
     lmh2*mh2*((-3 + 2*lmsb12 - 2*lmsb22)*msb12*msb22 - mh2*(msb12 + 2*(7 - 3*
     lmsb12 + lmsb22)*msb22) + (1 + 2*lmsb12 - 2*lmsb22)*msb24 + power2(mh2)) +
     mh2*msb12*msb22*power2(Pi) + 3*mh2*msb24*power2(Pi) + 3*msb22*power2(mh2)*
     power2(Pi) - 3*lmsb12*(-(msb12*msb24) + mh2*(-2*(-3 + lmsb22)*msb12*msb22
     + (23 - 6*lmsb22)*msb24) + msb26 - (msb12 + (-23 + 2*lmsb22)*msb22)*power2
     (mh2) + power3(mh2))))/(12.*msb12) - (DeltaInv(mh2,msb12,msb22)*(-(Ab*sa*(
     2*ca*mu + Ab*sa)) + mu2*(-1 + power2(sa)))*power2(1 - 2*power2(snb))*(
     msb22*(2*mh2*msb12 + 5*lmsb22*mh2*msb22 - lmsb22*msb12*msb22 + lmsb22*
     msb24) + lmh2*(mh2 - msb12 + 5*msb22)*power2(mh2) - lmsb12*(msb22*(-(msb12
     *msb22) + msb24) - (msb12 - 5*msb22)*power2(mh2) + 5*mh2*power2(msb22) +
     power3(mh2))))/(24.*msb12*msb22) + (mH2*(2*Ab*ca*mu*sa + power2(Ab)*(-1 +
     power2(sa)) - mu2*power2(sa))*power2(DeltaInv(mH2,msb12,msb22))*power2(1 -
     2*power2(snb))*(42*mH2*msb12*msb22 - 9*lmsb22*mH2*msb12*msb22 + 126*mH2*
     msb24 - 42*lmsb22*mH2*msb24 - 3*lmsb22*msb12*msb24 + 3*lmsb22*msb26 + 6*
     mH2*(2*mH2*msb22 + msb12*msb22 + 2*msb24)*power2(lmsb12) + 6*mH2*msb24*
     power2(lmsb22) + 126*msb22*power2(mH2) + 3*lmsb22*msb22*power2(mH2) + 6*
     msb22*power2(lmH2)*power2(mH2) + 3*lmH2*mH2*((-3 + 2*lmsb12 - 2*lmsb22)*
     msb12*msb22 - mH2*(msb12 + 2*(7 - 3*lmsb12 + lmsb22)*msb22) + (1 + 2*
     lmsb12 - 2*lmsb22)*msb24 + power2(mH2)) + mH2*msb12*msb22*power2(Pi) + 3*
     mH2*msb24*power2(Pi) + 3*msb22*power2(mH2)*power2(Pi) - 3*lmsb12*(-(msb12*
     msb24) + mH2*(-2*(-3 + lmsb22)*msb12*msb22 + (23 - 6*lmsb22)*msb24) +
     msb26 - (msb12 + (-23 + 2*lmsb22)*msb22)*power2(mH2) + power3(mH2))))/(12.
     *msb12) - (DeltaInv(mH2,msb12,msb22)*(2*Ab*ca*mu*sa + power2(Ab)*(-1 +
     power2(sa)) - mu2*power2(sa))*power2(1 - 2*power2(snb))*(msb22*(2*mH2*
     msb12 + 5*lmsb22*mH2*msb22 - lmsb22*msb12*msb22 + lmsb22*msb24) + lmH2*(
     mH2 - msb12 + 5*msb22)*power2(mH2) - lmsb12*(msb22*(-(msb12*msb22) + msb24
     ) - (msb12 - 5*msb22)*power2(mH2) + 5*mH2*power2(msb22) + power3(mH2))))/(
     24.*msb12*msb22) - (DeltaInv(mw2,msb12,mst12)*(mu2 - 2*Ab*cb*mu*sb - mu2*
     power2(cb) + power2(Ab)*power2(cb))*power2(snb)*(-1 + power2(snt))*(lmst12
     *mst12*(-(msb12*mst12) + mst14 + 5*mst12*mw2) + mw2*(lmw2*mw2*(5*mst12 +
     mw2) + msb12*(2*mst12 - lmw2*mw2)) - lmsb12*(mst12*mst14 + 5*mw2*power2(
     mst12) + 5*mst12*power2(mw2) - msb12*(power2(mst12) + power2(mw2)) +
     power3(mw2))))/(12.*msb12*mst12) + (DeltaInv(mw2,msb22,mst12)*(mu2 - 2*Ab*
     cb*mu*sb - mu2*power2(cb) + power2(Ab)*power2(cb))*(-1 + power2(snb))*(-1
     + power2(snt))*(lmst12*mst12*(-(msb22*mst12) + mst14 + 5*mst12*mw2) + mw2*
     (lmw2*mw2*(5*mst12 + mw2) + msb22*(2*mst12 - lmw2*mw2)) - lmsb22*(mst12*
     mst14 + 5*mw2*power2(mst12) + 5*mst12*power2(mw2) - msb22*(power2(mst12) +
     power2(mw2)) + power3(mw2))))/(12.*msb22*mst12) + (DeltaInv(mst22,mw2,
     msb12)*(mu2 - 2*Ab*cb*mu*sb - mu2*power2(cb) + power2(Ab)*power2(cb))*
     power2(snb)*power2(snt)*(lmst22*mst22*(-(msb12*mst22) + mst24 + 5*mst22*
     mw2) + mw2*(lmw2*mw2*(5*mst22 + mw2) + msb12*(2*mst22 - lmw2*mw2)) -
     lmsb12*(mst22*mst24 + 5*mw2*power2(mst22) + 5*mst22*power2(mw2) - msb12*(
     power2(mst22) + power2(mw2)) + power3(mw2))))/(12.*msb12*mst22) - (
     DeltaInv(mst22,mw2,msb22)*(mu2 - 2*Ab*cb*mu*sb - mu2*power2(cb) + power2(
     Ab)*power2(cb))*(-1 + power2(snb))*power2(snt)*(lmst22*mst22*(-(msb22*
     mst22) + mst24 + 5*mst22*mw2) + mw2*(lmw2*mw2*(5*mst22 + mw2) + msb22*(2*
     mst22 - lmw2*mw2)) - lmsb22*(mst22*mst24 + 5*mw2*power2(mst22) + 5*mst22*
     power2(mw2) - msb22*(power2(mst22) + power2(mw2)) + power3(mw2))))/(12.*
     msb22*mst22) + (DeltaInv(msb12,mz2,msb22)*(mu2 - 2*Ab*cb*mu*sb - mu2*
     power2(cb) + power2(Ab)*power2(cb))*(lmsb22*msb22*(-(msb12*msb22) + msb24
     + 5*msb22*mz2) + mz2*(lmz2*mz2*(5*msb22 + mz2) + msb12*(2*msb22 - lmz2*mz2
     )) - lmsb12*(msb22*msb24 + 5*mz2*power2(msb22) + 5*msb22*power2(mz2) -
     msb12*(power2(msb22) + power2(mz2)) + power3(mz2))))/(24.*msb12*msb22) +
     Fin3(mh2,msb12,msb12,Q2)*(2*msb14*(mu2 + Ab*sa*(2*ca*mu + Ab*sa) - mu2*
     power2(sa))*(-1 + power2(snb))*power2(snb)*power2(DeltaInv(mh2,msb12,msb12
     )) + 8*mh2*msb16*(mu2 + Ab*sa*(2*ca*mu + Ab*sa) - mu2*power2(sa))*(-1 +
     power2(snb))*power2(snb)*power3(DeltaInv(mh2,msb12,msb12))) + Fin3(mh2,
     msb12,msb22,Q2)*((msb22*(-(Ab*sa*(2*ca*mu + Ab*sa)) + mu2*(-1 + power2(sa)
     ))*power2(DeltaInv(mh2,msb12,msb22))*power2(mh2 - 2*mh2*power2(snb)))/
     msb12 + ((-2*msb12*msb24 - 2*mh2*(msb12*msb22 + msb24) + msb26 + msb22*
     power2(mh2))*(mu2 + Ab*sa*(2*ca*mu + Ab*sa) - mu2*power2(sa))*power2(mh2 -
     2*mh2*power2(snb))*power3(DeltaInv(mh2,msb12,msb22)))/msb12) + Fin3(mh2,
     msb22,msb22,Q2)*(2*msb24*(mu2 + Ab*sa*(2*ca*mu + Ab*sa) - mu2*power2(sa))*
     (-1 + power2(snb))*power2(snb)*power2(DeltaInv(mh2,msb22,msb22)) + 8*mh2*
     msb26*(mu2 + Ab*sa*(2*ca*mu + Ab*sa) - mu2*power2(sa))*(-1 + power2(snb))*
     power2(snb)*power3(DeltaInv(mh2,msb22,msb22))) + Fin3(mH2,msb12,msb12,Q2)*
     (-2*msb14*(2*Ab*ca*mu*sa + power2(Ab)*(-1 + power2(sa)) - mu2*power2(sa))*
     (-1 + power2(snb))*power2(snb)*power2(DeltaInv(mH2,msb12,msb12)) - 8*mH2*
     msb16*(2*Ab*ca*mu*sa + power2(Ab)*(-1 + power2(sa)) - mu2*power2(sa))*(-1
     + power2(snb))*power2(snb)*power3(DeltaInv(mH2,msb12,msb12))) + Fin3(mH2,
     msb12,msb22,Q2)*((msb22*(2*Ab*ca*mu*sa + power2(Ab)*(-1 + power2(sa)) -
     mu2*power2(sa))*power2(DeltaInv(mH2,msb12,msb22))*power2(mH2 - 2*mH2*
     power2(snb)))/msb12 - ((-2*msb12*msb24 - 2*mH2*(msb12*msb22 + msb24) +
     msb26 + msb22*power2(mH2))*(2*Ab*ca*mu*sa + power2(Ab)*(-1 + power2(sa)) -
     mu2*power2(sa))*power2(mH2 - 2*mH2*power2(snb))*power3(DeltaInv(mH2,msb12,
     msb22)))/msb12) + Fin3(mH2,msb22,msb22,Q2)*(-2*msb24*(2*Ab*ca*mu*sa +
     power2(Ab)*(-1 + power2(sa)) - mu2*power2(sa))*(-1 + power2(snb))*power2(
     snb)*power2(DeltaInv(mH2,msb22,msb22)) - 8*mH2*msb26*(2*Ab*ca*mu*sa +
     power2(Ab)*(-1 + power2(sa)) - mu2*power2(sa))*(-1 + power2(snb))*power2(
     snb)*power3(DeltaInv(mH2,msb22,msb22))) + Fin3(msb12,mA2,msb22,Q2)*((msb22
     *(-2*Ab*cb*mu*sb + power2(Ab)*(-1 + power2(cb)) - mu2*power2(cb))*power2(
     mA2)*power2(DeltaInv(msb12,mA2,msb22)))/msb12 - ((-2*Ab*cb*mu*sb + power2(
     Ab)*(-1 + power2(cb)) - mu2*power2(cb))*power2(mA2)*(-2*msb12*msb24 - 2*
     mA2*(msb12*msb22 + msb24) + msb26 + msb22*power2(mA2))*power3(DeltaInv(
     msb12,mA2,msb22)))/msb12) + Fin3(msb12,mst12,mC2,Q2)*(((4*invdmst*mst12*
     mt2*xt*(-(cb*mu*sb) + Ab*(-1 + power2(cb)))*power2(mC2)*power2(snb))/msb12
      - (2*mst12*power2(mC2)*power2(snb)*(-2*Ab*cb*mu*sb*(-1 + power2(snt)) +
     power2(Ab)*(-1 + power2(cb))*(-1 + power2(snt)) + mt2*power2(snt) + power2
     (cb)*(mu2 - mt2*power2(snt) - mu2*power2(snt))))/msb12)*power2(DeltaInv(
     msb12,mst12,mC2)) + ((-4*invdmst*mt2*xt*(-(cb*mu*sb) + Ab*(-1 + power2(cb)
     ))*power2(mC2)*(-2*msb12*mst14 - 2*mC2*(msb12*mst12 + mst14) + mst16 +
     mst12*power2(mC2))*power2(snb))/msb12 + (2*power2(mC2)*(-2*msb12*mst14 - 2
     *mC2*(msb12*mst12 + mst14) + mst16 + mst12*power2(mC2))*power2(snb)*(-2*Ab
     *cb*mu*sb*(-1 + power2(snt)) + power2(Ab)*(-1 + power2(cb))*(-1 + power2(
     snt)) + mt2*power2(snt) + power2(cb)*(mu2 - mt2*power2(snt) - mu2*power2(
     snt))))/msb12)*power3(DeltaInv(msb12,mst12,mC2))) + Fin3(msb12,mz2,msb22,
     Q2)*(-((msb22*(mu2 - 2*Ab*cb*mu*sb - mu2*power2(cb) + power2(Ab)*power2(cb
     ))*power2(mz2)*power2(DeltaInv(msb12,mz2,msb22)))/msb12) + ((msb26 + mz2*(
     -2*msb24 + msb22*mz2) - 2*msb12*(msb24 + msb22*mz2))*(mu2 - 2*Ab*cb*mu*sb
     - mu2*power2(cb) + power2(Ab)*power2(cb))*power2(mz2)*power3(DeltaInv(
     msb12,mz2,msb22)))/msb12) + Fin3(msb22,mst12,mC2,Q2)*(((-4*invdmst*mst12*
     mt2*xt*(-(cb*mu*sb) + Ab*(-1 + power2(cb)))*power2(mC2)*(-1 + power2(snb))
     )/msb22 + (2*mst12*power2(mC2)*(-1 + power2(snb))*(-2*Ab*cb*mu*sb*(-1 +
     power2(snt)) + power2(Ab)*(-1 + power2(cb))*(-1 + power2(snt)) + mt2*
     power2(snt) + power2(cb)*(mu2 - mt2*power2(snt) - mu2*power2(snt))))/msb22
     )*power2(DeltaInv(msb22,mst12,mC2)) + ((4*invdmst*mt2*xt*(-(cb*mu*sb) + Ab
     *(-1 + power2(cb)))*power2(mC2)*(-2*msb22*mst14 - 2*mC2*(msb22*mst12 +
     mst14) + mst16 + mst12*power2(mC2))*(-1 + power2(snb)))/msb22 - (2*power2(
     mC2)*(-2*msb22*mst14 - 2*mC2*(msb22*mst12 + mst14) + mst16 + mst12*power2(
     mC2))*(-1 + power2(snb))*(-2*Ab*cb*mu*sb*(-1 + power2(snt)) + power2(Ab)*(
     -1 + power2(cb))*(-1 + power2(snt)) + mt2*power2(snt) + power2(cb)*(mu2 -
     mt2*power2(snt) - mu2*power2(snt))))/msb22)*power3(DeltaInv(msb22,mst12,
     mC2))) + Fin3(mst22,msb12,mC2,Q2)*(((4*invdmst*mst22*mt2*xt*(Ab + cb*mu*sb
      - Ab*power2(cb))*power2(mC2)*power2(snb))/msb12 - (2*mst22*power2(mC2)*
     power2(snb)*(mt2*(-1 + power2(cb))*(-1 + power2(snt)) + (2*Ab*cb*mu*sb -
     power2(Ab)*(-1 + power2(cb)) + mu2*power2(cb))*power2(snt)))/msb12)*power2
     (DeltaInv(mst22,msb12,mC2)) + ((4*invdmst*mt2*xt*(-(cb*mu*sb) + Ab*(-1 +
     power2(cb)))*power2(mC2)*(-2*msb12*mst24 - 2*mC2*(msb12*mst22 + mst24) +
     mst26 + mst22*power2(mC2))*power2(snb))/msb12 + (2*power2(mC2)*(-2*msb12*
     mst24 - 2*mC2*(msb12*mst22 + mst24) + mst26 + mst22*power2(mC2))*power2(
     snb)*(mt2*(-1 + power2(cb))*(-1 + power2(snt)) + (2*Ab*cb*mu*sb - power2(
     Ab)*(-1 + power2(cb)) + mu2*power2(cb))*power2(snt)))/msb12)*power3(
     DeltaInv(mst22,msb12,mC2))) + Fin3(mst22,msb22,mC2,Q2)*(((4*invdmst*mst22*
     mt2*xt*(-(cb*mu*sb) + Ab*(-1 + power2(cb)))*power2(mC2)*(-1 + power2(snb))
     )/msb22 + (2*mst22*power2(mC2)*(-1 + power2(snb))*(mt2*(-1 + power2(cb))*(
     -1 + power2(snt)) + (2*Ab*cb*mu*sb - power2(Ab)*(-1 + power2(cb)) + mu2*
     power2(cb))*power2(snt)))/msb22)*power2(DeltaInv(mst22,msb22,mC2)) + ((-4*
     invdmst*mt2*xt*(-(cb*mu*sb) + Ab*(-1 + power2(cb)))*power2(mC2)*(-2*msb22*
     mst24 - 2*mC2*(msb22*mst22 + mst24) + mst26 + mst22*power2(mC2))*(-1 +
     power2(snb)))/msb22 - (2*power2(mC2)*(-2*msb22*mst24 - 2*mC2*(msb22*mst22
     + mst24) + mst26 + mst22*power2(mC2))*(-1 + power2(snb))*(mt2*(-1 + power2
     (cb))*(-1 + power2(snt)) + (2*Ab*cb*mu*sb - power2(Ab)*(-1 + power2(cb)) +
     mu2*power2(cb))*power2(snt)))/msb22)*power3(DeltaInv(mst22,msb22,mC2))) +
     Fin3(mst22,mw2,msb12,Q2)*((-2*mst22*(mu2 - 2*Ab*cb*mu*sb - mu2*power2(cb)
     + power2(Ab)*power2(cb))*power2(mw2)*power2(snb)*power2(snt)*power2(
     DeltaInv(mst22,mw2,msb12)))/msb12 + (2*(mst26 + mw2*(-2*mst24 + mst22*mw2)
     - 2*msb12*(mst24 + mst22*mw2))*(mu2 - 2*Ab*cb*mu*sb - mu2*power2(cb) +
     power2(Ab)*power2(cb))*power2(mw2)*power2(snb)*power2(snt)*power3(DeltaInv
     (mst22,mw2,msb12)))/msb12) + Fin3(mst22,mw2,msb22,Q2)*((2*mst22*(mu2 - 2*
     Ab*cb*mu*sb - mu2*power2(cb) + power2(Ab)*power2(cb))*power2(mw2)*(-1 +
     power2(snb))*power2(snt)*power2(DeltaInv(mst22,mw2,msb22)))/msb22 - (2*(
     mst26 + mw2*(-2*mst24 + mst22*mw2) - 2*msb22*(mst24 + mst22*mw2))*(mu2 - 2
     *Ab*cb*mu*sb - mu2*power2(cb) + power2(Ab)*power2(cb))*power2(mw2)*(-1 +
     power2(snb))*power2(snt)*power3(DeltaInv(mst22,mw2,msb22)))/msb22) + Fin3(
     mw2,msb12,mst12,Q2)*((2*mst12*(mu2 - 2*Ab*cb*mu*sb - mu2*power2(cb) +
     power2(Ab)*power2(cb))*power2(mw2)*power2(snb)*(-1 + power2(snt))*power2(
     DeltaInv(mw2,msb12,mst12)))/msb12 - (2*(mst16 + mw2*(-2*mst14 + mst12*mw2)
     - 2*msb12*(mst14 + mst12*mw2))*(mu2 - 2*Ab*cb*mu*sb - mu2*power2(cb) +
     power2(Ab)*power2(cb))*power2(mw2)*power2(snb)*(-1 + power2(snt))*power3(
     DeltaInv(mw2,msb12,mst12)))/msb12) + Fin3(mw2,msb22,mst12,Q2)*((-2*mst12*(
     mu2 - 2*Ab*cb*mu*sb - mu2*power2(cb) + power2(Ab)*power2(cb))*power2(mw2)*
     (-1 + power2(snb))*(-1 + power2(snt))*power2(DeltaInv(mw2,msb22,mst12)))/
     msb22 + (2*(mst16 + mw2*(-2*mst14 + mst12*mw2) - 2*msb22*(mst14 + mst12*
     mw2))*(mu2 - 2*Ab*cb*mu*sb - mu2*power2(cb) + power2(Ab)*power2(cb))*
     power2(mw2)*(-1 + power2(snb))*(-1 + power2(snt))*power3(DeltaInv(mw2,
     msb22,mst12)))/msb22) - 4*mh2*(42 + 8*lmh2*(-3 + lmsb12) - 12*lmsb12 + 4*
     power2(lmh2) + power2(Pi))*(-(Ab*sa*(2*ca*mu + Ab*sa)) + mu2*(-1 + power2(
     sa)))*(-1 + power2(snb))*power2(snb)*power3(DeltaInv(mh2,msb12,msb12))*
     power4(msb12) - 4*mH2*(42 + 8*lmH2*(-3 + lmsb12) - 12*lmsb12 + 4*power2(
     lmH2) + power2(Pi))*(2*Ab*ca*mu*sa + power2(Ab)*(-1 + power2(sa)) - mu2*
     power2(sa))*(-1 + power2(snb))*power2(snb)*power3(DeltaInv(mH2,msb12,msb12
     ))*power4(msb12) + Fin3(msb12,mt2,mu2,Q2)*((2*(msb14 + msb12*(mt2 + 3*mu2)
     )*power2(mu2)*power2(snb)*power2(DeltaInv(msb12,mt2,mu2)))/mt2 - (2*power2
     (mu2)*power2(snb)*power3(DeltaInv(msb12,mt2,mu2))*(msb16*(-mt2 + mu2) - 5*
     msb14*mu2*(2*mt2 + mu2) + msb12*(-5*mt2 + 3*mu2)*power2(mu2) + power4(
     msb12)))/mt2) - 4*mh2*(42 + 8*lmh2*(-3 + lmsb22) - 12*lmsb22 + 4*power2(
     lmh2) + power2(Pi))*(-(Ab*sa*(2*ca*mu + Ab*sa)) + mu2*(-1 + power2(sa)))*(
     -1 + power2(snb))*power2(snb)*power3(DeltaInv(mh2,msb22,msb22))*power4(
     msb22) - 4*mH2*(42 + 8*lmH2*(-3 + lmsb22) - 12*lmsb22 + 4*power2(lmH2) +
     power2(Pi))*(2*Ab*ca*mu*sa + power2(Ab)*(-1 + power2(sa)) - mu2*power2(sa)
     )*(-1 + power2(snb))*power2(snb)*power3(DeltaInv(mH2,msb22,msb22))*power4(
     msb22) - ((-2*Ab*cb*mu*sb + power2(Ab)*(-1 + power2(cb)) - mu2*power2(cb))
     *power2(mA2)*power3(DeltaInv(msb12,mA2,msb22))*(-(msb12*msb26*(-6*(12 +
     lmA2)*lmsb22 + 6*lmsb12*(-18 + lmA2 + 5*lmsb22) + 18*power2(lmsb12) + 12*
     power2(lmsb22) + 5*(42 + power2(Pi)))) - mA2*(3*msb26*(42 - 24*lmsb22 - 2*
     lmA2*(-6 + lmsb12 + lmsb22) + 2*lmsb12*(-12 + 5*lmsb22) - 2*power2(lmA2) +
     4*power2(lmsb12) + 4*power2(lmsb22) + power2(Pi)) + 2*msb12*msb24*(294 +
     30*lmsb12*(-6 + lmsb22) - 36*lmsb22 + 6*lmA2*(5*lmsb12 - 3*(2 + lmsb22)) +
     6*power2(lmA2) + 30*power2(lmsb12) + 6*power2(lmsb22) + 7*power2(Pi))) -
     power2(mA2)*(3*msb24*(42 + 2*lmA2*(-12 + 5*lmsb12 - lmsb22) + 12*lmsb22 -
     2*lmsb12*(12 + lmsb22) + 4*power2(lmA2) + 4*power2(lmsb12) - 2*power2(
     lmsb22) + power2(Pi)) + msb12*msb22*(6*lmA2*(-12 + 5*lmsb12 - lmsb22) + 6*
     lmsb12*(-18 + lmsb22) + 12*power2(lmA2) + 18*power2(lmsb12) + 5*(42 +
     power2(Pi)))) + 3*msb22*(42 + 2*lmA2*(-6 + 3*lmsb12 - lmsb22) + 2*lmsb12*(
     -12 + lmsb22) + 2*power2(lmA2) + 4*power2(lmsb12) + power2(Pi))*power3(mA2
     ) + 3*(42 + 2*lmsb12*(lmA2 + 3*(-4 + lmsb22)) - 2*(6 + lmA2)*lmsb22 + 4*
     power2(lmsb12) + 2*power2(lmsb22) + power2(Pi))*power4(msb22)))/(12.*msb12
     ) - ((mu2 + Ab*sa*(2*ca*mu + Ab*sa) - mu2*power2(sa))*power2(mh2 - 2*mh2*
     power2(snb))*power3(DeltaInv(mh2,msb12,msb22))*(msb12*msb26*(-6*(12 + lmh2
     )*lmsb22 + 6*lmsb12*(-18 + lmh2 + 5*lmsb22) + 18*power2(lmsb12) + 12*
     power2(lmsb22) + 5*(42 + power2(Pi))) + mh2*(3*msb26*(42 - 24*lmsb22 - 2*
     lmh2*(-6 + lmsb12 + lmsb22) + 2*lmsb12*(-12 + 5*lmsb22) - 2*power2(lmh2) +
     4*power2(lmsb12) + 4*power2(lmsb22) + power2(Pi)) + 2*msb12*msb24*(294 +
     30*lmsb12*(-6 + lmsb22) - 36*lmsb22 + 6*lmh2*(5*lmsb12 - 3*(2 + lmsb22)) +
     6*power2(lmh2) + 30*power2(lmsb12) + 6*power2(lmsb22) + 7*power2(Pi))) +
     power2(mh2)*(3*msb24*(42 + 2*lmh2*(-12 + 5*lmsb12 - lmsb22) + 12*lmsb22 -
     2*lmsb12*(12 + lmsb22) + 4*power2(lmh2) + 4*power2(lmsb12) - 2*power2(
     lmsb22) + power2(Pi)) + msb12*msb22*(6*lmh2*(-12 + 5*lmsb12 - lmsb22) + 6*
     lmsb12*(-18 + lmsb22) + 12*power2(lmh2) + 18*power2(lmsb12) + 5*(42 +
     power2(Pi)))) - 3*msb22*(42 + 2*lmh2*(-6 + 3*lmsb12 - lmsb22) + 2*lmsb12*(
     -12 + lmsb22) + 2*power2(lmh2) + 4*power2(lmsb12) + power2(Pi))*power3(mh2
     ) - 3*(42 + 2*lmsb12*(lmh2 + 3*(-4 + lmsb22)) - 2*(6 + lmh2)*lmsb22 + 4*
     power2(lmsb12) + 2*power2(lmsb22) + power2(Pi))*power4(msb22)))/(12.*msb12
     ) + ((2*Ab*ca*mu*sa + power2(Ab)*(-1 + power2(sa)) - mu2*power2(sa))*
     power2(mH2 - 2*mH2*power2(snb))*power3(DeltaInv(mH2,msb12,msb22))*(msb12*
     msb26*(-6*(12 + lmH2)*lmsb22 + 6*lmsb12*(-18 + lmH2 + 5*lmsb22) + 18*
     power2(lmsb12) + 12*power2(lmsb22) + 5*(42 + power2(Pi))) + mH2*(3*msb26*(
     42 - 24*lmsb22 - 2*lmH2*(-6 + lmsb12 + lmsb22) + 2*lmsb12*(-12 + 5*lmsb22)
     - 2*power2(lmH2) + 4*power2(lmsb12) + 4*power2(lmsb22) + power2(Pi)) + 2*
     msb12*msb24*(294 + 30*lmsb12*(-6 + lmsb22) - 36*lmsb22 + 6*lmH2*(5*lmsb12
     - 3*(2 + lmsb22)) + 6*power2(lmH2) + 30*power2(lmsb12) + 6*power2(lmsb22)
     + 7*power2(Pi))) + power2(mH2)*(3*msb24*(42 + 2*lmH2*(-12 + 5*lmsb12 -
     lmsb22) + 12*lmsb22 - 2*lmsb12*(12 + lmsb22) + 4*power2(lmH2) + 4*power2(
     lmsb12) - 2*power2(lmsb22) + power2(Pi)) + msb12*msb22*(6*lmH2*(-12 + 5*
     lmsb12 - lmsb22) + 6*lmsb12*(-18 + lmsb22) + 12*power2(lmH2) + 18*power2(
     lmsb12) + 5*(42 + power2(Pi)))) - 3*msb22*(42 + 2*lmH2*(-6 + 3*lmsb12 -
     lmsb22) + 2*lmsb12*(-12 + lmsb22) + 2*power2(lmH2) + 4*power2(lmsb12) +
     power2(Pi))*power3(mH2) - 3*(42 + 2*lmsb12*(lmH2 + 3*(-4 + lmsb22)) - 2*(6
      + lmH2)*lmsb22 + 4*power2(lmsb12) + 2*power2(lmsb22) + power2(Pi))*power4
     (msb22)))/(12.*msb12) - ((mu2 - 2*Ab*cb*mu*sb - mu2*power2(cb) + power2(Ab
     )*power2(cb))*power2(mz2)*power3(DeltaInv(msb12,mz2,msb22))*(3*mz2*(msb26*
     (42 + 2*lmsb12*(-12 + 5*lmsb22 - lmz2) + 12*lmz2 - 2*lmsb22*(12 + lmz2) +
     4*power2(lmsb12) + 4*power2(lmsb22) - 2*power2(lmz2) + power2(Pi)) + msb24
     *mz2*(42 - 2*lmsb12*(12 + lmsb22 - 5*lmz2) - 2*lmsb22*(-6 + lmz2) - 24*
     lmz2 + 4*power2(lmsb12) - 2*power2(lmsb22) + 4*power2(lmz2) + power2(Pi)))
     + msb12*(2*msb24*mz2*(294 - 36*lmz2 - 18*lmsb22*(2 + lmz2) + 30*lmsb12*(-6
      + lmsb22 + lmz2) + 30*power2(lmsb12) + 6*power2(lmsb22) + 6*power2(lmz2)
     + 7*power2(Pi)) + msb26*(-6*lmsb22*(12 + lmz2) + 6*lmsb12*(-18 + 5*lmsb22
     + lmz2) + 18*power2(lmsb12) + 12*power2(lmsb22) + 5*(42 + power2(Pi)))) +
     msb22*power2(mz2)*(-3*mz2*(42 + 2*lmsb12*(lmsb22 + 3*(-4 + lmz2)) - 2*(6 +
     lmsb22)*lmz2 + 4*power2(lmsb12) + 2*power2(lmz2) + power2(Pi)) + msb12*(-6
     *(12 + lmsb22)*lmz2 + 6*lmsb12*(-18 + lmsb22 + 5*lmz2) + 18*power2(lmsb12)
     + 12*power2(lmz2) + 5*(42 + power2(Pi)))) - 3*(42 - 2*lmsb22*(6 + lmz2) +
     2*lmsb12*(-12 + 3*lmsb22 + lmz2) + 4*power2(lmsb12) + 2*power2(lmsb22) +
     power2(Pi))*power4(msb22)))/(12.*msb12) + Fin3(mt2,mu2,msb22,Q2)*((-2*(
     msb24 + msb22*(mt2 + 3*mu2))*power2(mu2)*(-1 + power2(snb))*power2(
     DeltaInv(mt2,mu2,msb22)))/mt2 + (2*power2(mu2)*(-1 + power2(snb))*power3(
     DeltaInv(mt2,mu2,msb22))*(msb26*(-mt2 + mu2) - 5*msb24*mu2*(2*mt2 + mu2) +
     msb22*(-5*mt2 + 3*mu2)*power2(mu2) + power4(msb22)))/mt2) + ((mu2 - 2*Ab*
     cb*mu*sb - mu2*power2(cb) + power2(Ab)*power2(cb))*power2(mw2)*power2(snb)
     *(-1 + power2(snt))*power3(DeltaInv(mw2,msb12,mst12))*(3*mw2*(mst16*(42 +
     2*lmsb12*(-12 + 5*lmst12 - lmw2) + 12*lmw2 - 2*lmst12*(12 + lmw2) + 4*
     power2(lmsb12) + 4*power2(lmst12) - 2*power2(lmw2) + power2(Pi)) + mst14*
     mw2*(42 - 2*lmsb12*(12 + lmst12 - 5*lmw2) - 2*lmst12*(-6 + lmw2) - 24*lmw2
      + 4*power2(lmsb12) - 2*power2(lmst12) + 4*power2(lmw2) + power2(Pi))) +
     msb12*(2*mst14*mw2*(294 - 36*lmw2 - 18*lmst12*(2 + lmw2) + 30*lmsb12*(-6 +
     lmst12 + lmw2) + 30*power2(lmsb12) + 6*power2(lmst12) + 6*power2(lmw2) + 7
     *power2(Pi)) + mst16*(-6*lmst12*(12 + lmw2) + 6*lmsb12*(-18 + 5*lmst12 +
     lmw2) + 18*power2(lmsb12) + 12*power2(lmst12) + 5*(42 + power2(Pi)))) +
     mst12*power2(mw2)*(-3*mw2*(42 + 2*lmsb12*(lmst12 + 3*(-4 + lmw2)) - 2*(6 +
     lmst12)*lmw2 + 4*power2(lmsb12) + 2*power2(lmw2) + power2(Pi)) + msb12*(-6
     *(12 + lmst12)*lmw2 + 6*lmsb12*(-18 + lmst12 + 5*lmw2) + 18*power2(lmsb12)
     + 12*power2(lmw2) + 5*(42 + power2(Pi)))) - 3*(42 - 2*lmst12*(6 + lmw2) +
     2*lmsb12*(-12 + 3*lmst12 + lmw2) + 4*power2(lmsb12) + 2*power2(lmst12) +
     power2(Pi))*power4(mst12)))/(6.*msb12) - ((mu2 - 2*Ab*cb*mu*sb - mu2*
     power2(cb) + power2(Ab)*power2(cb))*power2(mw2)*(-1 + power2(snb))*(-1 +
     power2(snt))*power3(DeltaInv(mw2,msb22,mst12))*(3*mw2*(mst16*(42 + 2*
     lmsb22*(-12 + 5*lmst12 - lmw2) + 12*lmw2 - 2*lmst12*(12 + lmw2) + 4*power2
     (lmsb22) + 4*power2(lmst12) - 2*power2(lmw2) + power2(Pi)) + mst14*mw2*(42
      - 2*lmsb22*(12 + lmst12 - 5*lmw2) - 2*lmst12*(-6 + lmw2) - 24*lmw2 + 4*
     power2(lmsb22) - 2*power2(lmst12) + 4*power2(lmw2) + power2(Pi))) + msb22*
     (2*mst14*mw2*(294 - 36*lmw2 - 18*lmst12*(2 + lmw2) + 30*lmsb22*(-6 +
     lmst12 + lmw2) + 30*power2(lmsb22) + 6*power2(lmst12) + 6*power2(lmw2) + 7
     *power2(Pi)) + mst16*(-6*lmst12*(12 + lmw2) + 6*lmsb22*(-18 + 5*lmst12 +
     lmw2) + 18*power2(lmsb22) + 12*power2(lmst12) + 5*(42 + power2(Pi)))) +
     mst12*power2(mw2)*(-3*mw2*(42 + 2*lmsb22*(lmst12 + 3*(-4 + lmw2)) - 2*(6 +
     lmst12)*lmw2 + 4*power2(lmsb22) + 2*power2(lmw2) + power2(Pi)) + msb22*(-6
     *(12 + lmst12)*lmw2 + 6*lmsb22*(-18 + lmst12 + 5*lmw2) + 18*power2(lmsb22)
     + 12*power2(lmw2) + 5*(42 + power2(Pi)))) - 3*(42 - 2*lmst12*(6 + lmw2) +
     2*lmsb22*(-12 + 3*lmst12 + lmw2) + 4*power2(lmsb22) + 2*power2(lmst12) +
     power2(Pi))*power4(mst12)))/(6.*msb22) + power3(DeltaInv(msb12,mst12,mC2))
     *(-(invdmst*mt2*xt*(-(cb*mu*sb) + Ab*(-1 + power2(cb)))*power2(mC2)*power2
     (snb)*(-(msb12*mst16*(-6*(12 + lmC2)*lmst12 + 6*lmsb12*(-18 + lmC2 + 5*
     lmst12) + 18*power2(lmsb12) + 12*power2(lmst12) + 5*(42 + power2(Pi)))) -
     mC2*(3*mst16*(42 - 24*lmst12 - 2*lmC2*(-6 + lmsb12 + lmst12) + 2*lmsb12*(-
     12 + 5*lmst12) - 2*power2(lmC2) + 4*power2(lmsb12) + 4*power2(lmst12) +
     power2(Pi)) + 2*msb12*mst14*(294 + 30*lmsb12*(-6 + lmst12) - 36*lmst12 + 6
     *lmC2*(5*lmsb12 - 3*(2 + lmst12)) + 6*power2(lmC2) + 30*power2(lmsb12) + 6
     *power2(lmst12) + 7*power2(Pi))) - power2(mC2)*(3*mst14*(42 + 2*lmC2*(-12
     + 5*lmsb12 - lmst12) + 12*lmst12 - 2*lmsb12*(12 + lmst12) + 4*power2(lmC2)
     + 4*power2(lmsb12) - 2*power2(lmst12) + power2(Pi)) + msb12*mst12*(6*lmC2*
     (-12 + 5*lmsb12 - lmst12) + 6*lmsb12*(-18 + lmst12) + 12*power2(lmC2) + 18
     *power2(lmsb12) + 5*(42 + power2(Pi)))) + 3*mst12*(42 + 2*lmC2*(-6 + 3*
     lmsb12 - lmst12) + 2*lmsb12*(-12 + lmst12) + 2*power2(lmC2) + 4*power2(
     lmsb12) + power2(Pi))*power3(mC2) + 3*(42 + 2*lmsb12*(lmC2 + 3*(-4 +
     lmst12)) - 2*(6 + lmC2)*lmst12 + 4*power2(lmsb12) + 2*power2(lmst12) +
     power2(Pi))*power4(mst12)))/(3.*msb12) + (power2(mC2)*power2(snb)*(-2*Ab*
     cb*mu*sb*(-1 + power2(snt)) + power2(Ab)*(-1 + power2(cb))*(-1 + power2(
     snt)) + mt2*power2(snt) + power2(cb)*(mu2 - mt2*power2(snt) - mu2*power2(
     snt)))*(-(msb12*mst16*(-6*(12 + lmC2)*lmst12 + 6*lmsb12*(-18 + lmC2 + 5*
     lmst12) + 18*power2(lmsb12) + 12*power2(lmst12) + 5*(42 + power2(Pi)))) -
     mC2*(3*mst16*(42 - 24*lmst12 - 2*lmC2*(-6 + lmsb12 + lmst12) + 2*lmsb12*(-
     12 + 5*lmst12) - 2*power2(lmC2) + 4*power2(lmsb12) + 4*power2(lmst12) +
     power2(Pi)) + 2*msb12*mst14*(294 + 30*lmsb12*(-6 + lmst12) - 36*lmst12 + 6
     *lmC2*(5*lmsb12 - 3*(2 + lmst12)) + 6*power2(lmC2) + 30*power2(lmsb12) + 6
     *power2(lmst12) + 7*power2(Pi))) - power2(mC2)*(3*mst14*(42 + 2*lmC2*(-12
     + 5*lmsb12 - lmst12) + 12*lmst12 - 2*lmsb12*(12 + lmst12) + 4*power2(lmC2)
     + 4*power2(lmsb12) - 2*power2(lmst12) + power2(Pi)) + msb12*mst12*(6*lmC2*
     (-12 + 5*lmsb12 - lmst12) + 6*lmsb12*(-18 + lmst12) + 12*power2(lmC2) + 18
     *power2(lmsb12) + 5*(42 + power2(Pi)))) + 3*mst12*(42 + 2*lmC2*(-6 + 3*
     lmsb12 - lmst12) + 2*lmsb12*(-12 + lmst12) + 2*power2(lmC2) + 4*power2(
     lmsb12) + power2(Pi))*power3(mC2) + 3*(42 + 2*lmsb12*(lmC2 + 3*(-4 +
     lmst12)) - 2*(6 + lmC2)*lmst12 + 4*power2(lmsb12) + 2*power2(lmst12) +
     power2(Pi))*power4(mst12)))/(6.*msb12)) + power3(DeltaInv(msb22,mst12,mC2)
     )*((invdmst*mt2*xt*(-(cb*mu*sb) + Ab*(-1 + power2(cb)))*power2(mC2)*(-1 +
     power2(snb))*(-(msb22*mst16*(-6*(12 + lmC2)*lmst12 + 6*lmsb22*(-18 + lmC2
     + 5*lmst12) + 18*power2(lmsb22) + 12*power2(lmst12) + 5*(42 + power2(Pi)))
     ) - mC2*(3*mst16*(42 - 24*lmst12 - 2*lmC2*(-6 + lmsb22 + lmst12) + 2*
     lmsb22*(-12 + 5*lmst12) - 2*power2(lmC2) + 4*power2(lmsb22) + 4*power2(
     lmst12) + power2(Pi)) + 2*msb22*mst14*(294 + 30*lmsb22*(-6 + lmst12) - 36*
     lmst12 + 6*lmC2*(5*lmsb22 - 3*(2 + lmst12)) + 6*power2(lmC2) + 30*power2(
     lmsb22) + 6*power2(lmst12) + 7*power2(Pi))) - power2(mC2)*(3*mst14*(42 + 2
     *lmC2*(-12 + 5*lmsb22 - lmst12) + 12*lmst12 - 2*lmsb22*(12 + lmst12) + 4*
     power2(lmC2) + 4*power2(lmsb22) - 2*power2(lmst12) + power2(Pi)) + msb22*
     mst12*(6*lmC2*(-12 + 5*lmsb22 - lmst12) + 6*lmsb22*(-18 + lmst12) + 12*
     power2(lmC2) + 18*power2(lmsb22) + 5*(42 + power2(Pi)))) + 3*mst12*(42 + 2
     *lmC2*(-6 + 3*lmsb22 - lmst12) + 2*lmsb22*(-12 + lmst12) + 2*power2(lmC2)
     + 4*power2(lmsb22) + power2(Pi))*power3(mC2) + 3*(42 + 2*lmsb22*(lmC2 + 3*
     (-4 + lmst12)) - 2*(6 + lmC2)*lmst12 + 4*power2(lmsb22) + 2*power2(lmst12)
     + power2(Pi))*power4(mst12)))/(3.*msb22) - (power2(mC2)*(-1 + power2(snb))
     *(-2*Ab*cb*mu*sb*(-1 + power2(snt)) + power2(Ab)*(-1 + power2(cb))*(-1 +
     power2(snt)) + mt2*power2(snt) + power2(cb)*(mu2 - mt2*power2(snt) - mu2*
     power2(snt)))*(-(msb22*mst16*(-6*(12 + lmC2)*lmst12 + 6*lmsb22*(-18 + lmC2
      + 5*lmst12) + 18*power2(lmsb22) + 12*power2(lmst12) + 5*(42 + power2(Pi))
     )) - mC2*(3*mst16*(42 - 24*lmst12 - 2*lmC2*(-6 + lmsb22 + lmst12) + 2*
     lmsb22*(-12 + 5*lmst12) - 2*power2(lmC2) + 4*power2(lmsb22) + 4*power2(
     lmst12) + power2(Pi)) + 2*msb22*mst14*(294 + 30*lmsb22*(-6 + lmst12) - 36*
     lmst12 + 6*lmC2*(5*lmsb22 - 3*(2 + lmst12)) + 6*power2(lmC2) + 30*power2(
     lmsb22) + 6*power2(lmst12) + 7*power2(Pi))) - power2(mC2)*(3*mst14*(42 + 2
     *lmC2*(-12 + 5*lmsb22 - lmst12) + 12*lmst12 - 2*lmsb22*(12 + lmst12) + 4*
     power2(lmC2) + 4*power2(lmsb22) - 2*power2(lmst12) + power2(Pi)) + msb22*
     mst12*(6*lmC2*(-12 + 5*lmsb22 - lmst12) + 6*lmsb22*(-18 + lmst12) + 12*
     power2(lmC2) + 18*power2(lmsb22) + 5*(42 + power2(Pi)))) + 3*mst12*(42 + 2
     *lmC2*(-6 + 3*lmsb22 - lmst12) + 2*lmsb22*(-12 + lmst12) + 2*power2(lmC2)
     + 4*power2(lmsb22) + power2(Pi))*power3(mC2) + 3*(42 + 2*lmsb22*(lmC2 + 3*
     (-4 + lmst12)) - 2*(6 + lmC2)*lmst12 + 4*power2(lmsb22) + 2*power2(lmst12)
     + power2(Pi))*power4(mst12)))/(6.*msb22)) - ((mu2 - 2*Ab*cb*mu*sb - mu2*
     power2(cb) + power2(Ab)*power2(cb))*power2(mw2)*power2(snb)*power2(snt)*
     power3(DeltaInv(mst22,mw2,msb12))*(3*mw2*(mst26*(42 + 2*lmsb12*(-12 + 5*
     lmst22 - lmw2) + 12*lmw2 - 2*lmst22*(12 + lmw2) + 4*power2(lmsb12) + 4*
     power2(lmst22) - 2*power2(lmw2) + power2(Pi)) + mst24*mw2*(42 - 2*lmsb12*(
     12 + lmst22 - 5*lmw2) - 2*lmst22*(-6 + lmw2) - 24*lmw2 + 4*power2(lmsb12)
     - 2*power2(lmst22) + 4*power2(lmw2) + power2(Pi))) + msb12*(2*mst24*mw2*(
     294 - 36*lmw2 - 18*lmst22*(2 + lmw2) + 30*lmsb12*(-6 + lmst22 + lmw2) + 30
     *power2(lmsb12) + 6*power2(lmst22) + 6*power2(lmw2) + 7*power2(Pi)) +
     mst26*(-6*lmst22*(12 + lmw2) + 6*lmsb12*(-18 + 5*lmst22 + lmw2) + 18*
     power2(lmsb12) + 12*power2(lmst22) + 5*(42 + power2(Pi)))) + mst22*power2(
     mw2)*(-3*mw2*(42 + 2*lmsb12*(lmst22 + 3*(-4 + lmw2)) - 2*(6 + lmst22)*lmw2
      + 4*power2(lmsb12) + 2*power2(lmw2) + power2(Pi)) + msb12*(-6*(12 +
     lmst22)*lmw2 + 6*lmsb12*(-18 + lmst22 + 5*lmw2) + 18*power2(lmsb12) + 12*
     power2(lmw2) + 5*(42 + power2(Pi)))) - 3*(42 - 2*lmst22*(6 + lmw2) + 2*
     lmsb12*(-12 + 3*lmst22 + lmw2) + 4*power2(lmsb12) + 2*power2(lmst22) +
     power2(Pi))*power4(mst22)))/(6.*msb12) + ((mu2 - 2*Ab*cb*mu*sb - mu2*
     power2(cb) + power2(Ab)*power2(cb))*power2(mw2)*(-1 + power2(snb))*power2(
     snt)*power3(DeltaInv(mst22,mw2,msb22))*(3*mw2*(mst26*(42 + 2*lmsb22*(-12 +
     5*lmst22 - lmw2) + 12*lmw2 - 2*lmst22*(12 + lmw2) + 4*power2(lmsb22) + 4*
     power2(lmst22) - 2*power2(lmw2) + power2(Pi)) + mst24*mw2*(42 - 2*lmsb22*(
     12 + lmst22 - 5*lmw2) - 2*lmst22*(-6 + lmw2) - 24*lmw2 + 4*power2(lmsb22)
     - 2*power2(lmst22) + 4*power2(lmw2) + power2(Pi))) + msb22*(2*mst24*mw2*(
     294 - 36*lmw2 - 18*lmst22*(2 + lmw2) + 30*lmsb22*(-6 + lmst22 + lmw2) + 30
     *power2(lmsb22) + 6*power2(lmst22) + 6*power2(lmw2) + 7*power2(Pi)) +
     mst26*(-6*lmst22*(12 + lmw2) + 6*lmsb22*(-18 + 5*lmst22 + lmw2) + 18*
     power2(lmsb22) + 12*power2(lmst22) + 5*(42 + power2(Pi)))) + mst22*power2(
     mw2)*(-3*mw2*(42 + 2*lmsb22*(lmst22 + 3*(-4 + lmw2)) - 2*(6 + lmst22)*lmw2
      + 4*power2(lmsb22) + 2*power2(lmw2) + power2(Pi)) + msb22*(-6*(12 +
     lmst22)*lmw2 + 6*lmsb22*(-18 + lmst22 + 5*lmw2) + 18*power2(lmsb22) + 12*
     power2(lmw2) + 5*(42 + power2(Pi)))) - 3*(42 - 2*lmst22*(6 + lmw2) + 2*
     lmsb22*(-12 + 3*lmst22 + lmw2) + 4*power2(lmsb22) + 2*power2(lmst22) +
     power2(Pi))*power4(mst22)))/(6.*msb22) + power3(DeltaInv(mst22,msb12,mC2))
     *((invdmst*mt2*xt*(-(cb*mu*sb) + Ab*(-1 + power2(cb)))*power2(mC2)*power2(
     snb)*(-(msb12*mst26*(-6*(12 + lmC2)*lmst22 + 6*lmsb12*(-18 + lmC2 + 5*
     lmst22) + 18*power2(lmsb12) + 12*power2(lmst22) + 5*(42 + power2(Pi)))) -
     mC2*(3*mst26*(42 - 24*lmst22 - 2*lmC2*(-6 + lmsb12 + lmst22) + 2*lmsb12*(-
     12 + 5*lmst22) - 2*power2(lmC2) + 4*power2(lmsb12) + 4*power2(lmst22) +
     power2(Pi)) + 2*msb12*mst24*(294 + 30*lmsb12*(-6 + lmst22) - 36*lmst22 + 6
     *lmC2*(5*lmsb12 - 3*(2 + lmst22)) + 6*power2(lmC2) + 30*power2(lmsb12) + 6
     *power2(lmst22) + 7*power2(Pi))) - power2(mC2)*(3*mst24*(42 + 2*lmC2*(-12
     + 5*lmsb12 - lmst22) + 12*lmst22 - 2*lmsb12*(12 + lmst22) + 4*power2(lmC2)
     + 4*power2(lmsb12) - 2*power2(lmst22) + power2(Pi)) + msb12*mst22*(6*lmC2*
     (-12 + 5*lmsb12 - lmst22) + 6*lmsb12*(-18 + lmst22) + 12*power2(lmC2) + 18
     *power2(lmsb12) + 5*(42 + power2(Pi)))) + 3*mst22*(42 + 2*lmC2*(-6 + 3*
     lmsb12 - lmst22) + 2*lmsb12*(-12 + lmst22) + 2*power2(lmC2) + 4*power2(
     lmsb12) + power2(Pi))*power3(mC2) + 3*(42 + 2*lmsb12*(lmC2 + 3*(-4 +
     lmst22)) - 2*(6 + lmC2)*lmst22 + 4*power2(lmsb12) + 2*power2(lmst22) +
     power2(Pi))*power4(mst22)))/(3.*msb12) + (power2(mC2)*power2(snb)*(mt2*(-1
      + power2(cb))*(-1 + power2(snt)) + (2*Ab*cb*mu*sb - power2(Ab)*(-1 +
     power2(cb)) + mu2*power2(cb))*power2(snt))*(-(msb12*mst26*(-6*(12 + lmC2)*
     lmst22 + 6*lmsb12*(-18 + lmC2 + 5*lmst22) + 18*power2(lmsb12) + 12*power2(
     lmst22) + 5*(42 + power2(Pi)))) - mC2*(3*mst26*(42 - 24*lmst22 - 2*lmC2*(-
     6 + lmsb12 + lmst22) + 2*lmsb12*(-12 + 5*lmst22) - 2*power2(lmC2) + 4*
     power2(lmsb12) + 4*power2(lmst22) + power2(Pi)) + 2*msb12*mst24*(294 + 30*
     lmsb12*(-6 + lmst22) - 36*lmst22 + 6*lmC2*(5*lmsb12 - 3*(2 + lmst22)) + 6*
     power2(lmC2) + 30*power2(lmsb12) + 6*power2(lmst22) + 7*power2(Pi))) -
     power2(mC2)*(3*mst24*(42 + 2*lmC2*(-12 + 5*lmsb12 - lmst22) + 12*lmst22 -
     2*lmsb12*(12 + lmst22) + 4*power2(lmC2) + 4*power2(lmsb12) - 2*power2(
     lmst22) + power2(Pi)) + msb12*mst22*(6*lmC2*(-12 + 5*lmsb12 - lmst22) + 6*
     lmsb12*(-18 + lmst22) + 12*power2(lmC2) + 18*power2(lmsb12) + 5*(42 +
     power2(Pi)))) + 3*mst22*(42 + 2*lmC2*(-6 + 3*lmsb12 - lmst22) + 2*lmsb12*(
     -12 + lmst22) + 2*power2(lmC2) + 4*power2(lmsb12) + power2(Pi))*power3(mC2
     ) + 3*(42 + 2*lmsb12*(lmC2 + 3*(-4 + lmst22)) - 2*(6 + lmC2)*lmst22 + 4*
     power2(lmsb12) + 2*power2(lmst22) + power2(Pi))*power4(mst22)))/(6.*msb12)
     ) + power3(DeltaInv(mst22,msb22,mC2))*(-(invdmst*mt2*xt*(-(cb*mu*sb) + Ab*
     (-1 + power2(cb)))*power2(mC2)*(-1 + power2(snb))*(-(msb22*mst26*(-6*(12 +
     lmC2)*lmst22 + 6*lmsb22*(-18 + lmC2 + 5*lmst22) + 18*power2(lmsb22) + 12*
     power2(lmst22) + 5*(42 + power2(Pi)))) - mC2*(3*mst26*(42 - 24*lmst22 - 2*
     lmC2*(-6 + lmsb22 + lmst22) + 2*lmsb22*(-12 + 5*lmst22) - 2*power2(lmC2) +
     4*power2(lmsb22) + 4*power2(lmst22) + power2(Pi)) + 2*msb22*mst24*(294 +
     30*lmsb22*(-6 + lmst22) - 36*lmst22 + 6*lmC2*(5*lmsb22 - 3*(2 + lmst22)) +
     6*power2(lmC2) + 30*power2(lmsb22) + 6*power2(lmst22) + 7*power2(Pi))) -
     power2(mC2)*(3*mst24*(42 + 2*lmC2*(-12 + 5*lmsb22 - lmst22) + 12*lmst22 -
     2*lmsb22*(12 + lmst22) + 4*power2(lmC2) + 4*power2(lmsb22) - 2*power2(
     lmst22) + power2(Pi)) + msb22*mst22*(6*lmC2*(-12 + 5*lmsb22 - lmst22) + 6*
     lmsb22*(-18 + lmst22) + 12*power2(lmC2) + 18*power2(lmsb22) + 5*(42 +
     power2(Pi)))) + 3*mst22*(42 + 2*lmC2*(-6 + 3*lmsb22 - lmst22) + 2*lmsb22*(
     -12 + lmst22) + 2*power2(lmC2) + 4*power2(lmsb22) + power2(Pi))*power3(mC2
     ) + 3*(42 + 2*lmsb22*(lmC2 + 3*(-4 + lmst22)) - 2*(6 + lmC2)*lmst22 + 4*
     power2(lmsb22) + 2*power2(lmst22) + power2(Pi))*power4(mst22)))/(3.*msb22)
     - (power2(mC2)*(-1 + power2(snb))*(mt2*(-1 + power2(cb))*(-1 + power2(snt)
     ) + (2*Ab*cb*mu*sb - power2(Ab)*(-1 + power2(cb)) + mu2*power2(cb))*power2
     (snt))*(-(msb22*mst26*(-6*(12 + lmC2)*lmst22 + 6*lmsb22*(-18 + lmC2 + 5*
     lmst22) + 18*power2(lmsb22) + 12*power2(lmst22) + 5*(42 + power2(Pi)))) -
     mC2*(3*mst26*(42 - 24*lmst22 - 2*lmC2*(-6 + lmsb22 + lmst22) + 2*lmsb22*(-
     12 + 5*lmst22) - 2*power2(lmC2) + 4*power2(lmsb22) + 4*power2(lmst22) +
     power2(Pi)) + 2*msb22*mst24*(294 + 30*lmsb22*(-6 + lmst22) - 36*lmst22 + 6
     *lmC2*(5*lmsb22 - 3*(2 + lmst22)) + 6*power2(lmC2) + 30*power2(lmsb22) + 6
     *power2(lmst22) + 7*power2(Pi))) - power2(mC2)*(3*mst24*(42 + 2*lmC2*(-12
     + 5*lmsb22 - lmst22) + 12*lmst22 - 2*lmsb22*(12 + lmst22) + 4*power2(lmC2)
     + 4*power2(lmsb22) - 2*power2(lmst22) + power2(Pi)) + msb22*mst22*(6*lmC2*
     (-12 + 5*lmsb22 - lmst22) + 6*lmsb22*(-18 + lmst22) + 12*power2(lmC2) + 18
     *power2(lmsb22) + 5*(42 + power2(Pi)))) + 3*mst22*(42 + 2*lmC2*(-6 + 3*
     lmsb22 - lmst22) + 2*lmsb22*(-12 + lmst22) + 2*power2(lmC2) + 4*power2(
     lmsb22) + power2(Pi))*power3(mC2) + 3*(42 + 2*lmsb22*(lmC2 + 3*(-4 +
     lmst22)) - 2*(6 + lmC2)*lmst22 + 4*power2(lmsb22) + 2*power2(lmst22) +
     power2(Pi))*power4(mst22)))/(6.*msb22)) + (power2(mu2)*power2(snb)*power3(
     DeltaInv(msb12,mt2,mu2))*(mu2*(msb14*mu2*(mu2*(42 + 6*(-15 + 4*lmt2)*lmu2
     + 6*lmsb12*(9 - 4*lmt2 + lmu2) - 9*power2(lmsb12) + 15*power2(lmu2) +
     power2(Pi)) + 5*mt2*(210 + 6*lmsb12*(-3 + 3*lmt2 - 2*lmu2) - 36*lmu2 + 6*
     lmt2*(-21 + 4*lmu2) + 3*power2(lmsb12) + 21*power2(lmt2) + 6*power2(lmu2)
     + 5*power2(Pi))) + msb16*(3*mu2*(126 - 84*lmt2 + 10*lmsb12*(-3 + 2*lmt2 -
     lmu2) + 6*lmu2 + 8*lmt2*lmu2 + 5*power2(lmsb12) + 14*power2(lmt2) - power2
     (lmu2) + 3*power2(Pi)) + mt2*(672 - 378*lmt2 + 30*lmsb12*(-6 + 3*lmt2 -
     lmu2) - 18*lmu2 + 36*lmt2*lmu2 + 30*power2(lmsb12) + 63*power2(lmt2) + 3*
     power2(lmu2) + 16*power2(Pi)))) + msb12*(3*mt2*(-2*(15 + lmsb12)*lmu2 + 2*
     lmt2*(-21 + lmsb12 + 6*lmu2) + 7*power2(lmt2) + 5*power2(lmu2) + 2*(42 +
     power2(Pi))) - mu2*(-6*(9 + lmsb12)*lmu2 + 6*lmt2*(-15 + lmsb12 + 4*lmu2)
     + 15*power2(lmt2) + 9*power2(lmu2) + 4*(42 + power2(Pi))))*power3(mu2) + (
     mt2*(42 + 6*lmsb12*(-3 + lmt2) - 18*lmt2 + 3*power2(lmsb12) + 3*power2(
     lmt2) + power2(Pi)) - mu2*(210 + 24*lmt2*(-6 + lmu2) - 18*lmu2 + 6*lmsb12*
     (4*lmt2 - 3*(1 + lmu2)) + 3*power2(lmsb12) + 24*power2(lmt2) + 3*power2(
     lmu2) + 5*power2(Pi)))*power4(msb12) - (42 + 6*lmsb12*(-3 + lmt2) - 18*
     lmt2 + 3*power2(lmsb12) + 3*power2(lmt2) + power2(Pi))*power5(msb12)))/(3.
     *mt2) + (power2(mu2)*(-1 + power2(snb))*power3(DeltaInv(mt2,mu2,msb22))*(
     mu2*(-(msb24*mu2*(mu2*(42 + 6*(-15 + 4*lmt2)*lmu2 + 6*lmsb22*(9 - 4*lmt2 +
     lmu2) - 9*power2(lmsb22) + 15*power2(lmu2) + power2(Pi)) + 5*mt2*(210 + 6*
     lmsb22*(-3 + 3*lmt2 - 2*lmu2) - 36*lmu2 + 6*lmt2*(-21 + 4*lmu2) + 3*power2
     (lmsb22) + 21*power2(lmt2) + 6*power2(lmu2) + 5*power2(Pi)))) - msb26*(3*
     mu2*(126 - 84*lmt2 + 10*lmsb22*(-3 + 2*lmt2 - lmu2) + 6*lmu2 + 8*lmt2*lmu2
      + 5*power2(lmsb22) + 14*power2(lmt2) - power2(lmu2) + 3*power2(Pi)) + mt2
     *(672 - 378*lmt2 + 30*lmsb22*(-6 + 3*lmt2 - lmu2) - 18*lmu2 + 36*lmt2*lmu2
      + 30*power2(lmsb22) + 63*power2(lmt2) + 3*power2(lmu2) + 16*power2(Pi))))
     + msb22*(-3*mt2*(-2*(15 + lmsb22)*lmu2 + 2*lmt2*(-21 + lmsb22 + 6*lmu2) +
     7*power2(lmt2) + 5*power2(lmu2) + 2*(42 + power2(Pi))) + mu2*(-6*(9 +
     lmsb22)*lmu2 + 6*lmt2*(-15 + lmsb22 + 4*lmu2) + 15*power2(lmt2) + 9*power2
     (lmu2) + 4*(42 + power2(Pi))))*power3(mu2) + (-(mt2*(42 + 6*lmsb22*(-3 +
     lmt2) - 18*lmt2 + 3*power2(lmsb22) + 3*power2(lmt2) + power2(Pi))) + mu2*(
     210 + 24*lmt2*(-6 + lmu2) - 18*lmu2 + 6*lmsb22*(4*lmt2 - 3*(1 + lmu2)) + 3
     *power2(lmsb22) + 24*power2(lmt2) + 3*power2(lmu2) + 5*power2(Pi)))*power4
     (msb22) + (42 + 6*lmsb22*(-3 + lmt2) - 18*lmt2 + 3*power2(lmsb22) + 3*
     power2(lmt2) + power2(Pi))*power5(msb22)))/(3.*mt2) + (36 - 72*lmt2 - 108*
     lmu2 - 35*invdct*mC2 - 60*invdct*lmC2*mC2 + 84*invdct*lmt2*mC2 - (9*mA2)/
     msb12 - (9*mH2)/msb12 - (9*mA2)/msb22 - (18*mC2)/msb22 + (18*lmC2*mC2)/
     msb22 - (9*mH2)/msb22 - (18*msb12)/msb22 + (18*lmsb12*msb12)/msb22 - (18*
     msb22)/msb12 + (18*lmsb22*msb22)/msb12 - (18*mC2)/mst12 + (18*lmC2*mC2)/
     mst12 - (18*msb22)/mst12 + (18*lmsb22*msb22)/mst12 - (18*mst12)/msb22 + (
     18*lmst12*mst12)/msb22 - (35*mC2)/mt2 + (12*lmC2*mC2)/mt2 + (12*lmt2*mC2)/
     mt2 + (18*mt2)/msb22 + (18*lmsb22*mt2)/msb22 - (36*lmt2*mt2)/msb22 - (18*
     mt2)/mst22 + (18*lmsb22*mt2)/mst22 - (18*lmC2*mC2*mt2)/(msb22*mst22) + (18
     *lmsb22*mC2*mt2)/(msb22*mst22) - (36*lmsb12*msb12)/(msb12 - mu2) + (36*
     lmu2*msb12)/(msb12 - mu2) - (36*lmsb22*msb22)/(msb22 - mu2) + (36*lmu2*
     msb22)/(msb22 - mu2) - (36*lmst12*mst12)/(mst12 - mu2) + (36*lmu2*mst12)/(
     mst12 - mu2) + (18*mu2)/msb12 + (18*lmsb12*mu2)/msb12 - (36*lmu2*mu2)/
     msb12 + (36*mu2)/msb22 + (18*lmsb12*mu2)/msb22 + (18*lmsb22*mu2)/msb22 - (
     36*lmt2*mu2)/msb22 - (36*lmu2*mu2)/msb22 - (9*lmh2*mh2*mu2)/(msb12*msb22)
     + (9*lmsb12*mh2*mu2)/(msb12*msb22) + (18*mu2)/mst12 + (18*lmsb22*mu2)/
     mst12 - (36*lmu2*mu2)/mst12 - (36*mu2)/mt2 + (18*lmsb22*mu2*mw2)/(msb22*
     mst12) - (18*lmw2*mu2*mw2)/(msb22*mst12) + (9*lmsb12*mu2*mz2)/(msb12*msb22
     ) - (9*lmz2*mu2*mz2)/(msb12*msb22) - (18*Ab*ca*lmh2*mh2*mu*sa)/(msb12*
     msb22) + (18*Ab*ca*lmsb12*mh2*mu*sa)/(msb12*msb22) - (18*Ab*ca*lmsb12*mH2*
     mu*sa)/(msb12*msb22) + (18*Ab*cb*lmsb12*mA2*mu*sb)/(msb12*msb22) - (36*Ab*
     cb*lmC2*mC2*mu*sb)/(msb22*mst12) + (36*Ab*cb*lmsb22*mC2*mu*sb)/(msb22*
     mst12) - (36*Ab*cb*lmsb22*mu*mw2*sb)/(msb22*mst12) + (36*Ab*cb*lmw2*mu*mw2
     *sb)/(msb22*mst12) - (18*Ab*cb*lmsb12*mu*mz2*sb)/(msb12*msb22) + (18*Ab*cb
     *lmz2*mu*mz2*sb)/(msb12*msb22) - (18*power2(Ab))/msb12 + (18*lmsb12*power2
     (Ab))/msb12 - (36*power2(Ab))/msb22 + (18*lmsb12*power2(Ab))/msb22 + (18*
     lmsb22*power2(Ab))/msb22 + (9*lmsb12*mA2*power2(Ab))/(msb12*msb22) + (9*
     lmsb12*mH2*power2(Ab))/(msb12*msb22) - (18*power2(Ab))/mst12 + (18*lmsb22*
     power2(Ab))/mst12 - (18*lmC2*mC2*power2(Ab))/(msb22*mst12) + (18*lmsb22*
     mC2*power2(Ab))/(msb22*mst12) - 54*lmz2*power2(cb) + 35*invdct*mC2*power2(
     cb) + 60*invdct*lmC2*mC2*power2(cb) - 84*invdct*lmt2*mC2*power2(cb) + (9*
     mA2*power2(cb))/msb12 + (9*mA2*power2(cb))/msb22 + (18*mC2*power2(cb))/
     msb22 - (18*lmC2*mC2*power2(cb))/msb22 + (18*mC2*power2(cb))/mst12 - (18*
     lmC2*mC2*power2(cb))/mst12 + (35*mC2*power2(cb))/mt2 - (12*lmC2*mC2*power2
     (cb))/mt2 - (12*lmt2*mC2*power2(cb))/mt2 + (18*mt2*power2(cb))/msb22 - (18
     *lmsb22*mt2*power2(cb))/msb22 + (18*mt2*power2(cb))/mst22 - (18*lmsb22*mt2
     *power2(cb))/mst22 + (18*lmC2*mC2*mt2*power2(cb))/(msb22*mst22) - (18*
     lmsb22*mC2*mt2*power2(cb))/(msb22*mst22) + (9*lmsb12*mA2*mu2*power2(cb))/(
     msb12*msb22) - (18*lmC2*mC2*mu2*power2(cb))/(msb22*mst12) + (18*lmsb22*mC2
     *mu2*power2(cb))/(msb22*mst12) + 35*invdtw*mw2*power2(cb) - 84*invdtw*lmt2
     *mw2*power2(cb) + 60*invdtw*lmw2*mw2*power2(cb) - (18*mw2*power2(cb))/
     msb22 + (18*lmw2*mw2*power2(cb))/msb22 - (18*mw2*power2(cb))/mst12 + (18*
     lmw2*mw2*power2(cb))/mst12 - (35*mw2*power2(cb))/mt2 + (12*lmt2*mw2*power2
     (cb))/mt2 + (12*lmw2*mw2*power2(cb))/mt2 - (18*lmsb22*mu2*mw2*power2(cb))/
     (msb22*mst12) + (18*lmw2*mu2*mw2*power2(cb))/(msb22*mst12) - (9*mz2*power2
     (cb))/msb12 + (9*lmz2*mz2*power2(cb))/msb12 - (9*mz2*power2(cb))/msb22 + (
     9*lmz2*mz2*power2(cb))/msb22 - (9*lmsb12*mu2*mz2*power2(cb))/(msb12*msb22)
     + (9*lmz2*mu2*mz2*power2(cb))/(msb12*msb22) - (9*lmsb12*mA2*power2(Ab)*
     power2(cb))/(msb12*msb22) + (18*lmC2*mC2*power2(Ab)*power2(cb))/(msb22*
     mst12) - (18*lmsb22*mC2*power2(Ab)*power2(cb))/(msb22*mst12) + (18*lmsb22*
     mw2*power2(Ab)*power2(cb))/(msb22*mst12) - (18*lmw2*mw2*power2(Ab)*power2(
     cb))/(msb22*mst12) + (9*lmsb12*mz2*power2(Ab)*power2(cb))/(msb12*msb22) -
     (9*lmz2*mz2*power2(Ab)*power2(cb))/(msb12*msb22) + (9*lmA2*(-2*Ab*cb*mA2*
     mu*sb + 6*msb12*msb22*(-1 + power2(cb)) + mA2*power2(Ab)*(-1 + power2(cb))
     + mA2*(msb12 + msb22 - msb12*power2(cb) - msb22*power2(cb) - mu2*power2(cb
     ))))/(msb12*msb22) - 35*power2(invdct)*power2(mC2) + 12*lmC2*power2(invdct
     )*power2(mC2) + 12*lmt2*power2(invdct)*power2(mC2) + 35*power2(cb)*power2(
     invdct)*power2(mC2) - 12*lmC2*power2(cb)*power2(invdct)*power2(mC2) - 12*
     lmt2*power2(cb)*power2(invdct)*power2(mC2) - (36*lmt2*power2(mu2))/(msb22*
     mt2) + (36*lmu2*power2(mu2))/(msb22*mt2) - 35*power2(cb)*power2(invdtw)*
     power2(mw2) + 12*lmt2*power2(cb)*power2(invdtw)*power2(mw2) + 12*lmw2*
     power2(cb)*power2(invdtw)*power2(mw2) - 54*lmh2*power2(sa) - (9*mh2*power2
     (sa))/msb12 + (9*lmh2*mh2*power2(sa))/msb12 + (9*mH2*power2(sa))/msb12 - (
     9*mh2*power2(sa))/msb22 + (9*lmh2*mh2*power2(sa))/msb22 + (9*mH2*power2(sa
     ))/msb22 + (9*lmh2*mh2*mu2*power2(sa))/(msb12*msb22) - (9*lmsb12*mh2*mu2*
     power2(sa))/(msb12*msb22) + (9*lmsb12*mH2*mu2*power2(sa))/(msb12*msb22) -
     (9*lmh2*mh2*power2(Ab)*power2(sa))/(msb12*msb22) + (9*lmsb12*mh2*power2(Ab
     )*power2(sa))/(msb12*msb22) - (9*lmsb12*mH2*power2(Ab)*power2(sa))/(msb12*
     msb22) - 288*power2(snb) + 144*lmsb12*power2(snb) + 144*lmsb22*power2(snb)
     - (18*mC2*power2(snb))/msb12 + (18*lmC2*mC2*power2(snb))/msb12 + (18*mC2*
     power2(snb))/msb22 - (18*lmC2*mC2*power2(snb))/msb22 + (144*msb12*power2(
     snb))/msb22 - (144*lmsb12*msb12*power2(snb))/msb22 + (144*msb22*power2(snb
     ))/msb12 - (144*lmsb22*msb22*power2(snb))/msb12 - (18*msb12*power2(snb))/
     mst12 + (18*lmsb12*msb12*power2(snb))/mst12 + (18*msb22*power2(snb))/mst12
      - (18*lmsb22*msb22*power2(snb))/mst12 - (18*mst12*power2(snb))/msb12 + (
     18*lmst12*mst12*power2(snb))/msb12 + (18*mst12*power2(snb))/msb22 - (18*
     lmst12*mst12*power2(snb))/msb22 + (18*mt2*power2(snb))/msb12 + (18*lmsb12*
     mt2*power2(snb))/msb12 - (36*lmt2*mt2*power2(snb))/msb12 - (18*mt2*power2(
     snb))/msb22 - (18*lmsb22*mt2*power2(snb))/msb22 + (36*lmt2*mt2*power2(snb)
     )/msb22 + (18*lmsb12*mt2*power2(snb))/mst22 - (18*lmsb22*mt2*power2(snb))/
     mst22 - (18*lmC2*mC2*mt2*power2(snb))/(msb12*mst22) + (18*lmsb12*mC2*mt2*
     power2(snb))/(msb12*mst22) + (18*lmC2*mC2*mt2*power2(snb))/(msb22*mst22) -
     (18*lmsb22*mC2*mt2*power2(snb))/(msb22*mst22) + (18*mu2*power2(snb))/msb12
      + (36*lmh2*mu2*power2(snb))/msb12 - (18*lmsb12*mu2*power2(snb))/msb12 - (
     36*lmt2*mu2*power2(snb))/msb12 - (18*mu2*power2(snb))/msb22 + (36*lmh2*mu2
     *power2(snb))/msb22 - (36*lmsb12*mu2*power2(snb))/msb22 - (18*lmsb22*mu2*
     power2(snb))/msb22 + (36*lmt2*mu2*power2(snb))/msb22 + (36*lmh2*mh2*mu2*
     power2(snb))/(msb12*msb22) - (36*lmsb12*mh2*mu2*power2(snb))/(msb12*msb22)
     + (18*lmsb12*mu2*power2(snb))/mst12 - (18*lmsb22*mu2*power2(snb))/mst12 +
     (18*lmsb12*mu2*mw2*power2(snb))/(msb12*mst12) - (18*lmw2*mu2*mw2*power2(
     snb))/(msb12*mst12) - (18*lmsb22*mu2*mw2*power2(snb))/(msb22*mst12) + (18*
     lmw2*mu2*mw2*power2(snb))/(msb22*mst12) + (72*Ab*ca*lmh2*mu*sa*power2(snb)
     )/msb12 + (72*Ab*ca*lmh2*mu*sa*power2(snb))/msb22 + (72*Ab*ca*lmh2*mh2*mu*
     sa*power2(snb))/(msb12*msb22) - (72*Ab*ca*lmsb12*mh2*mu*sa*power2(snb))/(
     msb12*msb22) + (72*Ab*ca*lmsb12*mH2*mu*sa*power2(snb))/(msb12*msb22) - (36
     *Ab*cb*lmC2*mC2*mu*sb*power2(snb))/(msb12*mst12) + (36*Ab*cb*lmsb12*mC2*mu
     *sb*power2(snb))/(msb12*mst12) + (36*Ab*cb*lmC2*mC2*mu*sb*power2(snb))/(
     msb22*mst12) - (36*Ab*cb*lmsb22*mC2*mu*sb*power2(snb))/(msb22*mst12) - (36
     *Ab*cb*lmsb12*mu*mw2*sb*power2(snb))/(msb12*mst12) + (36*Ab*cb*lmw2*mu*mw2
     *sb*power2(snb))/(msb12*mst12) + (36*Ab*cb*lmsb22*mu*mw2*sb*power2(snb))/(
     msb22*mst12) - (36*Ab*cb*lmw2*mu*mw2*sb*power2(snb))/(msb22*mst12) - (18*
     power2(Ab)*power2(snb))/msb12 - (18*lmsb12*power2(Ab)*power2(snb))/msb12 +
     (18*power2(Ab)*power2(snb))/msb22 - (36*lmsb12*power2(Ab)*power2(snb))/
     msb22 - (18*lmsb22*power2(Ab)*power2(snb))/msb22 - (36*lmsb12*mH2*power2(
     Ab)*power2(snb))/(msb12*msb22) + (18*lmsb12*power2(Ab)*power2(snb))/mst12
     - (18*lmsb22*power2(Ab)*power2(snb))/mst12 - (18*lmC2*mC2*power2(Ab)*
     power2(snb))/(msb12*mst12) + (18*lmsb12*mC2*power2(Ab)*power2(snb))/(msb12
     *mst12) + (18*lmC2*mC2*power2(Ab)*power2(snb))/(msb22*mst12) - (18*lmsb22*
     mC2*power2(Ab)*power2(snb))/(msb22*mst12) + (18*mC2*power2(cb)*power2(snb)
     )/msb12 - (18*lmC2*mC2*power2(cb)*power2(snb))/msb12 - (18*mC2*power2(cb)*
     power2(snb))/msb22 + (18*lmC2*mC2*power2(cb)*power2(snb))/msb22 + (18*mt2*
     power2(cb)*power2(snb))/msb12 - (18*lmsb12*mt2*power2(cb)*power2(snb))/
     msb12 - (18*mt2*power2(cb)*power2(snb))/msb22 + (18*lmsb22*mt2*power2(cb)*
     power2(snb))/msb22 - (18*lmsb12*mt2*power2(cb)*power2(snb))/mst22 + (18*
     lmsb22*mt2*power2(cb)*power2(snb))/mst22 + (18*lmC2*mC2*mt2*power2(cb)*
     power2(snb))/(msb12*mst22) - (18*lmsb12*mC2*mt2*power2(cb)*power2(snb))/(
     msb12*mst22) - (18*lmC2*mC2*mt2*power2(cb)*power2(snb))/(msb22*mst22) + (
     18*lmsb22*mC2*mt2*power2(cb)*power2(snb))/(msb22*mst22) - (18*lmC2*mC2*mu2
     *power2(cb)*power2(snb))/(msb12*mst12) + (18*lmsb12*mC2*mu2*power2(cb)*
     power2(snb))/(msb12*mst12) + (18*lmC2*mC2*mu2*power2(cb)*power2(snb))/(
     msb22*mst12) - (18*lmsb22*mC2*mu2*power2(cb)*power2(snb))/(msb22*mst12) -
     (18*mw2*power2(cb)*power2(snb))/msb12 + (18*lmw2*mw2*power2(cb)*power2(snb
     ))/msb12 + (18*mw2*power2(cb)*power2(snb))/msb22 - (18*lmw2*mw2*power2(cb)
     *power2(snb))/msb22 - (18*lmsb12*mu2*mw2*power2(cb)*power2(snb))/(msb12*
     mst12) + (18*lmw2*mu2*mw2*power2(cb)*power2(snb))/(msb12*mst12) + (18*
     lmsb22*mu2*mw2*power2(cb)*power2(snb))/(msb22*mst12) - (18*lmw2*mu2*mw2*
     power2(cb)*power2(snb))/(msb22*mst12) + (18*lmC2*mC2*power2(Ab)*power2(cb)
     *power2(snb))/(msb12*mst12) - (18*lmsb12*mC2*power2(Ab)*power2(cb)*power2(
     snb))/(msb12*mst12) - (18*lmC2*mC2*power2(Ab)*power2(cb)*power2(snb))/(
     msb22*mst12) + (18*lmsb22*mC2*power2(Ab)*power2(cb)*power2(snb))/(msb22*
     mst12) + (18*lmsb12*mw2*power2(Ab)*power2(cb)*power2(snb))/(msb12*mst12) -
     (18*lmw2*mw2*power2(Ab)*power2(cb)*power2(snb))/(msb12*mst12) - (18*lmsb22
     *mw2*power2(Ab)*power2(cb)*power2(snb))/(msb22*mst12) + (18*lmw2*mw2*
     power2(Ab)*power2(cb)*power2(snb))/(msb22*mst12) - (36*lmt2*power2(mu2)*
     power2(snb))/(msb12*mt2) + (36*lmu2*power2(mu2)*power2(snb))/(msb12*mt2) +
     (36*lmt2*power2(mu2)*power2(snb))/(msb22*mt2) - (36*lmu2*power2(mu2)*
     power2(snb))/(msb22*mt2) - (36*lmh2*mu2*power2(sa)*power2(snb))/msb12 - (
     36*lmh2*mu2*power2(sa)*power2(snb))/msb22 - (36*lmh2*mh2*mu2*power2(sa)*
     power2(snb))/(msb12*msb22) + (36*lmsb12*mh2*mu2*power2(sa)*power2(snb))/(
     msb12*msb22) - (36*lmsb12*mH2*mu2*power2(sa)*power2(snb))/(msb12*msb22) +
     (36*lmh2*power2(Ab)*power2(sa)*power2(snb))/msb12 + (36*lmh2*power2(Ab)*
     power2(sa)*power2(snb))/msb22 + (36*lmh2*mh2*power2(Ab)*power2(sa)*power2(
     snb))/(msb12*msb22) - (36*lmsb12*mh2*power2(Ab)*power2(sa)*power2(snb))/(
     msb12*msb22) + (36*lmsb12*mH2*power2(Ab)*power2(sa)*power2(snb))/(msb12*
     msb22) + (18*mC2*power2(snt))/mst12 - (18*lmC2*mC2*power2(snt))/mst12 + (
     18*msb22*power2(snt))/mst12 - (18*lmsb22*msb22*power2(snt))/mst12 + (18*
     mst12*power2(snt))/msb22 - (18*lmst12*mst12*power2(snt))/msb22 - (18*mC2*
     power2(snt))/mst22 + (18*lmC2*mC2*power2(snt))/mst22 - (18*msb22*power2(
     snt))/mst22 + (18*lmsb22*msb22*power2(snt))/mst22 - (18*mst22*power2(snt))
     /msb22 + (18*lmst22*mst22*power2(snt))/msb22 - (18*mt2*power2(snt))/mst12
     + (18*lmsb22*mt2*power2(snt))/mst12 - (18*lmC2*mC2*mt2*power2(snt))/(msb22
     *mst12) + (18*lmsb22*mC2*mt2*power2(snt))/(msb22*mst12) + (18*mt2*power2(
     snt))/mst22 - (18*lmsb22*mt2*power2(snt))/mst22 + (18*lmC2*mC2*mt2*power2(
     snt))/(msb22*mst22) - (18*lmsb22*mC2*mt2*power2(snt))/(msb22*mst22) + (36*
     lmst12*mst12*power2(snt))/(mst12 - mu2) - (36*lmu2*mst12*power2(snt))/(
     mst12 - mu2) + (36*lmu2*mst22*power2(snt))/(mst22 - mu2) - (18*mu2*power2(
     snt))/mst12 - (18*lmsb22*mu2*power2(snt))/mst12 + (36*lmu2*mu2*power2(snt)
     )/mst12 + (18*mu2*power2(snt))/mst22 + (18*lmsb22*mu2*power2(snt))/mst22 -
     (36*lmu2*mu2*power2(snt))/mst22 + (36*lmst22*mst22*power2(snt))/(-mst22 +
     mu2) - (18*lmsb22*mu2*mw2*power2(snt))/(msb22*mst12) + (18*lmw2*mu2*mw2*
     power2(snt))/(msb22*mst12) + (18*lmsb22*mu2*mw2*power2(snt))/(msb22*mst22)
     - (18*lmw2*mu2*mw2*power2(snt))/(msb22*mst22) + (36*Ab*cb*lmC2*mC2*mu*sb*
     power2(snt))/(msb22*mst12) - (36*Ab*cb*lmsb22*mC2*mu*sb*power2(snt))/(
     msb22*mst12) - (36*Ab*cb*lmC2*mC2*mu*sb*power2(snt))/(msb22*mst22) + (36*
     Ab*cb*lmsb22*mC2*mu*sb*power2(snt))/(msb22*mst22) + (36*Ab*cb*lmsb22*mu*
     mw2*sb*power2(snt))/(msb22*mst12) - (36*Ab*cb*lmw2*mu*mw2*sb*power2(snt))/
     (msb22*mst12) - (36*Ab*cb*lmsb22*mu*mw2*sb*power2(snt))/(msb22*mst22) + (
     36*Ab*cb*lmw2*mu*mw2*sb*power2(snt))/(msb22*mst22) + (18*power2(Ab)*power2
     (snt))/mst12 - (18*lmsb22*power2(Ab)*power2(snt))/mst12 + (18*lmC2*mC2*
     power2(Ab)*power2(snt))/(msb22*mst12) - (18*lmsb22*mC2*power2(Ab)*power2(
     snt))/(msb22*mst12) - (18*power2(Ab)*power2(snt))/mst22 + (18*lmsb22*
     power2(Ab)*power2(snt))/mst22 - (18*lmC2*mC2*power2(Ab)*power2(snt))/(
     msb22*mst22) + (18*lmsb22*mC2*power2(Ab)*power2(snt))/(msb22*mst22) - (18*
     mC2*power2(cb)*power2(snt))/mst12 + (18*lmC2*mC2*power2(cb)*power2(snt))/
     mst12 + (18*mC2*power2(cb)*power2(snt))/mst22 - (18*lmC2*mC2*power2(cb)*
     power2(snt))/mst22 + (18*mt2*power2(cb)*power2(snt))/mst12 - (18*lmsb22*
     mt2*power2(cb)*power2(snt))/mst12 + (18*lmC2*mC2*mt2*power2(cb)*power2(snt
     ))/(msb22*mst12) - (18*lmsb22*mC2*mt2*power2(cb)*power2(snt))/(msb22*mst12
     ) - (18*mt2*power2(cb)*power2(snt))/mst22 + (18*lmsb22*mt2*power2(cb)*
     power2(snt))/mst22 - (18*lmC2*mC2*mt2*power2(cb)*power2(snt))/(msb22*mst22
     ) + (18*lmsb22*mC2*mt2*power2(cb)*power2(snt))/(msb22*mst22) + (18*lmC2*
     mC2*mu2*power2(cb)*power2(snt))/(msb22*mst12) - (18*lmsb22*mC2*mu2*power2(
     cb)*power2(snt))/(msb22*mst12) - (18*lmC2*mC2*mu2*power2(cb)*power2(snt))/
     (msb22*mst22) + (18*lmsb22*mC2*mu2*power2(cb)*power2(snt))/(msb22*mst22) +
     (18*mw2*power2(cb)*power2(snt))/mst12 - (18*lmw2*mw2*power2(cb)*power2(snt
     ))/mst12 - (18*mw2*power2(cb)*power2(snt))/mst22 + (18*lmw2*mw2*power2(cb)
     *power2(snt))/mst22 + (18*lmsb22*mu2*mw2*power2(cb)*power2(snt))/(msb22*
     mst12) - (18*lmw2*mu2*mw2*power2(cb)*power2(snt))/(msb22*mst12) - (18*
     lmsb22*mu2*mw2*power2(cb)*power2(snt))/(msb22*mst22) + (18*lmw2*mu2*mw2*
     power2(cb)*power2(snt))/(msb22*mst22) - (18*lmC2*mC2*power2(Ab)*power2(cb)
     *power2(snt))/(msb22*mst12) + (18*lmsb22*mC2*power2(Ab)*power2(cb)*power2(
     snt))/(msb22*mst12) + (18*lmC2*mC2*power2(Ab)*power2(cb)*power2(snt))/(
     msb22*mst22) - (18*lmsb22*mC2*power2(Ab)*power2(cb)*power2(snt))/(msb22*
     mst22) - (18*lmsb22*mw2*power2(Ab)*power2(cb)*power2(snt))/(msb22*mst12) +
     (18*lmw2*mw2*power2(Ab)*power2(cb)*power2(snt))/(msb22*mst12) + (18*lmsb22
     *mw2*power2(Ab)*power2(cb)*power2(snt))/(msb22*mst22) - (18*lmw2*mw2*
     power2(Ab)*power2(cb)*power2(snt))/(msb22*mst22) + (18*msb12*power2(snb)*
     power2(snt))/mst12 - (18*lmsb12*msb12*power2(snb)*power2(snt))/mst12 - (18
     *msb22*power2(snb)*power2(snt))/mst12 + (18*lmsb22*msb22*power2(snb)*
     power2(snt))/mst12 + (18*mst12*power2(snb)*power2(snt))/msb12 - (18*lmst12
     *mst12*power2(snb)*power2(snt))/msb12 - (18*mst12*power2(snb)*power2(snt))
     /msb22 + (18*lmst12*mst12*power2(snb)*power2(snt))/msb22 - (18*msb12*
     power2(snb)*power2(snt))/mst22 + (18*lmsb12*msb12*power2(snb)*power2(snt))
     /mst22 + (18*msb22*power2(snb)*power2(snt))/mst22 - (18*lmsb22*msb22*
     power2(snb)*power2(snt))/mst22 - (18*mst22*power2(snb)*power2(snt))/msb12
     + (18*lmst22*mst22*power2(snb)*power2(snt))/msb12 + (18*mst22*power2(snb)*
     power2(snt))/msb22 - (18*lmst22*mst22*power2(snb)*power2(snt))/msb22 + (18
     *lmsb12*mt2*power2(snb)*power2(snt))/mst12 - (18*lmsb22*mt2*power2(snb)*
     power2(snt))/mst12 - (18*lmC2*mC2*mt2*power2(snb)*power2(snt))/(msb12*
     mst12) + (18*lmsb12*mC2*mt2*power2(snb)*power2(snt))/(msb12*mst12) + (18*
     lmC2*mC2*mt2*power2(snb)*power2(snt))/(msb22*mst12) - (18*lmsb22*mC2*mt2*
     power2(snb)*power2(snt))/(msb22*mst12) - (18*lmsb12*mt2*power2(snb)*power2
     (snt))/mst22 + (18*lmsb22*mt2*power2(snb)*power2(snt))/mst22 + (18*lmC2*
     mC2*mt2*power2(snb)*power2(snt))/(msb12*mst22) - (18*lmsb12*mC2*mt2*power2
     (snb)*power2(snt))/(msb12*mst22) - (18*lmC2*mC2*mt2*power2(snb)*power2(snt
     ))/(msb22*mst22) + (18*lmsb22*mC2*mt2*power2(snb)*power2(snt))/(msb22*
     mst22) - (18*lmsb12*mu2*power2(snb)*power2(snt))/mst12 + (18*lmsb22*mu2*
     power2(snb)*power2(snt))/mst12 + (18*lmsb12*mu2*power2(snb)*power2(snt))/
     mst22 - (18*lmsb22*mu2*power2(snb)*power2(snt))/mst22 - (18*lmsb12*mu2*mw2
     *power2(snb)*power2(snt))/(msb12*mst12) + (18*lmw2*mu2*mw2*power2(snb)*
     power2(snt))/(msb12*mst12) + (18*lmsb22*mu2*mw2*power2(snb)*power2(snt))/(
     msb22*mst12) - (18*lmw2*mu2*mw2*power2(snb)*power2(snt))/(msb22*mst12) + (
     18*lmsb12*mu2*mw2*power2(snb)*power2(snt))/(msb12*mst22) - (18*lmw2*mu2*
     mw2*power2(snb)*power2(snt))/(msb12*mst22) - (18*lmsb22*mu2*mw2*power2(snb
     )*power2(snt))/(msb22*mst22) + (18*lmw2*mu2*mw2*power2(snb)*power2(snt))/(
     msb22*mst22) + (36*Ab*cb*lmC2*mC2*mu*sb*power2(snb)*power2(snt))/(msb12*
     mst12) - (36*Ab*cb*lmsb12*mC2*mu*sb*power2(snb)*power2(snt))/(msb12*mst12)
     - (36*Ab*cb*lmC2*mC2*mu*sb*power2(snb)*power2(snt))/(msb22*mst12) + (36*Ab
     *cb*lmsb22*mC2*mu*sb*power2(snb)*power2(snt))/(msb22*mst12) - (36*Ab*cb*
     lmC2*mC2*mu*sb*power2(snb)*power2(snt))/(msb12*mst22) + (36*Ab*cb*lmsb12*
     mC2*mu*sb*power2(snb)*power2(snt))/(msb12*mst22) + (36*Ab*cb*lmC2*mC2*mu*
     sb*power2(snb)*power2(snt))/(msb22*mst22) - (36*Ab*cb*lmsb22*mC2*mu*sb*
     power2(snb)*power2(snt))/(msb22*mst22) + (36*Ab*cb*lmsb12*mu*mw2*sb*power2
     (snb)*power2(snt))/(msb12*mst12) - (36*Ab*cb*lmw2*mu*mw2*sb*power2(snb)*
     power2(snt))/(msb12*mst12) - (36*Ab*cb*lmsb22*mu*mw2*sb*power2(snb)*power2
     (snt))/(msb22*mst12) + (36*Ab*cb*lmw2*mu*mw2*sb*power2(snb)*power2(snt))/(
     msb22*mst12) - (36*Ab*cb*lmsb12*mu*mw2*sb*power2(snb)*power2(snt))/(msb12*
     mst22) + (36*Ab*cb*lmw2*mu*mw2*sb*power2(snb)*power2(snt))/(msb12*mst22) +
     (36*Ab*cb*lmsb22*mu*mw2*sb*power2(snb)*power2(snt))/(msb22*mst22) - (36*Ab
     *cb*lmw2*mu*mw2*sb*power2(snb)*power2(snt))/(msb22*mst22) - (18*lmsb12*
     power2(Ab)*power2(snb)*power2(snt))/mst12 + (18*lmsb22*power2(Ab)*power2(
     snb)*power2(snt))/mst12 + (18*lmC2*mC2*power2(Ab)*power2(snb)*power2(snt))
     /(msb12*mst12) - (18*lmsb12*mC2*power2(Ab)*power2(snb)*power2(snt))/(msb12
     *mst12) - (18*lmC2*mC2*power2(Ab)*power2(snb)*power2(snt))/(msb22*mst12) +
     (18*lmsb22*mC2*power2(Ab)*power2(snb)*power2(snt))/(msb22*mst12) + (18*
     lmsb12*power2(Ab)*power2(snb)*power2(snt))/mst22 - (18*lmsb22*power2(Ab)*
     power2(snb)*power2(snt))/mst22 - (18*lmC2*mC2*power2(Ab)*power2(snb)*
     power2(snt))/(msb12*mst22) + (18*lmsb12*mC2*power2(Ab)*power2(snb)*power2(
     snt))/(msb12*mst22) + (18*lmC2*mC2*power2(Ab)*power2(snb)*power2(snt))/(
     msb22*mst22) - (18*lmsb22*mC2*power2(Ab)*power2(snb)*power2(snt))/(msb22*
     mst22) - (18*lmsb12*mt2*power2(cb)*power2(snb)*power2(snt))/mst12 + (18*
     lmsb22*mt2*power2(cb)*power2(snb)*power2(snt))/mst12 + (18*lmC2*mC2*mt2*
     power2(cb)*power2(snb)*power2(snt))/(msb12*mst12) - (18*lmsb12*mC2*mt2*
     power2(cb)*power2(snb)*power2(snt))/(msb12*mst12) - (18*lmC2*mC2*mt2*
     power2(cb)*power2(snb)*power2(snt))/(msb22*mst12) + (18*lmsb22*mC2*mt2*
     power2(cb)*power2(snb)*power2(snt))/(msb22*mst12) + (18*lmsb12*mt2*power2(
     cb)*power2(snb)*power2(snt))/mst22 - (18*lmsb22*mt2*power2(cb)*power2(snb)
     *power2(snt))/mst22 - (18*lmC2*mC2*mt2*power2(cb)*power2(snb)*power2(snt))
     /(msb12*mst22) + (18*lmsb12*mC2*mt2*power2(cb)*power2(snb)*power2(snt))/(
     msb12*mst22) + (18*lmC2*mC2*mt2*power2(cb)*power2(snb)*power2(snt))/(msb22
     *mst22) - (18*lmsb22*mC2*mt2*power2(cb)*power2(snb)*power2(snt))/(msb22*
     mst22) + (18*lmC2*mC2*mu2*power2(cb)*power2(snb)*power2(snt))/(msb12*mst12
     ) - (18*lmsb12*mC2*mu2*power2(cb)*power2(snb)*power2(snt))/(msb12*mst12) -
     (18*lmC2*mC2*mu2*power2(cb)*power2(snb)*power2(snt))/(msb22*mst12) + (18*
     lmsb22*mC2*mu2*power2(cb)*power2(snb)*power2(snt))/(msb22*mst12) - (18*
     lmC2*mC2*mu2*power2(cb)*power2(snb)*power2(snt))/(msb12*mst22) + (18*
     lmsb12*mC2*mu2*power2(cb)*power2(snb)*power2(snt))/(msb12*mst22) + (18*
     lmC2*mC2*mu2*power2(cb)*power2(snb)*power2(snt))/(msb22*mst22) - (18*
     lmsb22*mC2*mu2*power2(cb)*power2(snb)*power2(snt))/(msb22*mst22) + (18*
     lmsb12*mu2*mw2*power2(cb)*power2(snb)*power2(snt))/(msb12*mst12) - (18*
     lmw2*mu2*mw2*power2(cb)*power2(snb)*power2(snt))/(msb12*mst12) - (18*
     lmsb22*mu2*mw2*power2(cb)*power2(snb)*power2(snt))/(msb22*mst12) + (18*
     lmw2*mu2*mw2*power2(cb)*power2(snb)*power2(snt))/(msb22*mst12) - (18*
     lmsb12*mu2*mw2*power2(cb)*power2(snb)*power2(snt))/(msb12*mst22) + (18*
     lmw2*mu2*mw2*power2(cb)*power2(snb)*power2(snt))/(msb12*mst22) + (18*
     lmsb22*mu2*mw2*power2(cb)*power2(snb)*power2(snt))/(msb22*mst22) - (18*
     lmw2*mu2*mw2*power2(cb)*power2(snb)*power2(snt))/(msb22*mst22) - (18*lmC2*
     mC2*power2(Ab)*power2(cb)*power2(snb)*power2(snt))/(msb12*mst12) + (18*
     lmsb12*mC2*power2(Ab)*power2(cb)*power2(snb)*power2(snt))/(msb12*mst12) +
     (18*lmC2*mC2*power2(Ab)*power2(cb)*power2(snb)*power2(snt))/(msb22*mst12)
     - (18*lmsb22*mC2*power2(Ab)*power2(cb)*power2(snb)*power2(snt))/(msb22*
     mst12) + (18*lmC2*mC2*power2(Ab)*power2(cb)*power2(snb)*power2(snt))/(
     msb12*mst22) - (18*lmsb12*mC2*power2(Ab)*power2(cb)*power2(snb)*power2(snt
     ))/(msb12*mst22) - (18*lmC2*mC2*power2(Ab)*power2(cb)*power2(snb)*power2(
     snt))/(msb22*mst22) + (18*lmsb22*mC2*power2(Ab)*power2(cb)*power2(snb)*
     power2(snt))/(msb22*mst22) - (18*lmsb12*mw2*power2(Ab)*power2(cb)*power2(
     snb)*power2(snt))/(msb12*mst12) + (18*lmw2*mw2*power2(Ab)*power2(cb)*
     power2(snb)*power2(snt))/(msb12*mst12) + (18*lmsb22*mw2*power2(Ab)*power2(
     cb)*power2(snb)*power2(snt))/(msb22*mst12) - (18*lmw2*mw2*power2(Ab)*
     power2(cb)*power2(snb)*power2(snt))/(msb22*mst12) + (18*lmsb12*mw2*power2(
     Ab)*power2(cb)*power2(snb)*power2(snt))/(msb12*mst22) - (18*lmw2*mw2*
     power2(Ab)*power2(cb)*power2(snb)*power2(snt))/(msb12*mst22) - (18*lmsb22*
     mw2*power2(Ab)*power2(cb)*power2(snb)*power2(snt))/(msb22*mst22) + (18*
     lmw2*mw2*power2(Ab)*power2(cb)*power2(snb)*power2(snt))/(msb22*mst22) + (
     45*power2(invdct)*power3(mC2))/mt2 + (12*lmC2*power2(invdct)*power3(mC2))/
     mt2 + (12*lmt2*power2(invdct)*power3(mC2))/mt2 - (45*power2(cb)*power2(
     invdct)*power3(mC2))/mt2 - (12*lmC2*power2(cb)*power2(invdct)*power3(mC2))
     /mt2 - (12*lmt2*power2(cb)*power2(invdct)*power3(mC2))/mt2 + 10*power3(
     invdct)*power3(mC2) + 24*lmC2*power3(invdct)*power3(mC2) + 24*lmt2*power3(
     invdct)*power3(mC2) - 10*power2(cb)*power3(invdct)*power3(mC2) - 24*lmC2*
     power2(cb)*power3(invdct)*power3(mC2) - 24*lmt2*power2(cb)*power3(invdct)*
     power3(mC2) + (45*power2(cb)*power2(invdtw)*power3(mw2))/mt2 + (12*lmt2*
     power2(cb)*power2(invdtw)*power3(mw2))/mt2 + (12*lmw2*power2(cb)*power2(
     invdtw)*power3(mw2))/mt2 - 10*power2(cb)*power3(invdtw)*power3(mw2) - 24*
     lmt2*power2(cb)*power3(invdtw)*power3(mw2) - 24*lmw2*power2(cb)*power3(
     invdtw)*power3(mw2) + 10*power4(invdct)*power4(mC2) + 24*lmC2*power4(
     invdct)*power4(mC2) + 24*lmt2*power4(invdct)*power4(mC2) - 10*power2(cb)*
     power4(invdct)*power4(mC2) - 24*lmC2*power2(cb)*power4(invdct)*power4(mC2)
     - 24*lmt2*power2(cb)*power4(invdct)*power4(mC2) + 10*power2(cb)*power4(
     invdtw)*power4(mw2) + 24*lmt2*power2(cb)*power4(invdtw)*power4(mw2) + 24*
     lmw2*power2(cb)*power4(invdtw)*power4(mw2) + 288*power4(snb) - 144*lmsb12*
     power4(snb) - 144*lmsb22*power4(snb) - (144*msb12*power4(snb))/msb22 + (
     144*lmsb12*msb12*power4(snb))/msb22 - (144*msb22*power4(snb))/msb12 + (144
     *lmsb22*msb22*power4(snb))/msb12 - (36*lmh2*mu2*power4(snb))/msb12 + (36*
     lmsb12*mu2*power4(snb))/msb12 - (36*lmh2*mu2*power4(snb))/msb22 + (36*
     lmsb12*mu2*power4(snb))/msb22 - (36*lmh2*mh2*mu2*power4(snb))/(msb12*msb22
     ) + (36*lmsb12*mh2*mu2*power4(snb))/(msb12*msb22) - (72*Ab*ca*lmh2*mu*sa*
     power4(snb))/msb12 - (72*Ab*ca*lmh2*mu*sa*power4(snb))/msb22 - (72*Ab*ca*
     lmh2*mh2*mu*sa*power4(snb))/(msb12*msb22) + (72*Ab*ca*lmsb12*mh2*mu*sa*
     power4(snb))/(msb12*msb22) - (72*Ab*ca*lmsb12*mH2*mu*sa*power4(snb))/(
     msb12*msb22) + (36*lmsb12*power2(Ab)*power4(snb))/msb12 + (36*lmsb12*
     power2(Ab)*power4(snb))/msb22 + (36*lmsb12*mH2*power2(Ab)*power4(snb))/(
     msb12*msb22) + (36*lmh2*mu2*power2(sa)*power4(snb))/msb12 + (36*lmh2*mu2*
     power2(sa)*power4(snb))/msb22 + (36*lmh2*mh2*mu2*power2(sa)*power4(snb))/(
     msb12*msb22) - (36*lmsb12*mh2*mu2*power2(sa)*power4(snb))/(msb12*msb22) +
     (36*lmsb12*mH2*mu2*power2(sa)*power4(snb))/(msb12*msb22) - (36*lmh2*power2
     (Ab)*power2(sa)*power4(snb))/msb12 - (36*lmh2*power2(Ab)*power2(sa)*power4
     (snb))/msb22 - (36*lmh2*mh2*power2(Ab)*power2(sa)*power4(snb))/(msb12*
     msb22) + (36*lmsb12*mh2*power2(Ab)*power2(sa)*power4(snb))/(msb12*msb22) -
     (36*lmsb12*mH2*power2(Ab)*power2(sa)*power4(snb))/(msb12*msb22) + (9*lmH2*
     (-6*msb12*msb22 + 6*msb12*msb22*power2(sa) + 4*msb12*mu2*power2(sa)*power2
     (snb) + 4*msb22*mu2*power2(sa)*power2(snb) + 2*Ab*ca*mu*sa*(4*(msb12 +
     msb22)*(-1 + power2(snb))*power2(snb) + mH2*power2(1 - 2*power2(snb))) +
     power2(Ab)*(-1 + power2(sa))*(4*(msb12 + msb22)*(-1 + power2(snb))*power2(
     snb) + mH2*power2(1 - 2*power2(snb))) + mH2*(msb12 + msb22 - msb12*power2(
     sa) - msb22*power2(sa) - mu2*power2(sa)*power2(1 - 2*power2(snb))) - 4*
     msb12*mu2*power2(sa)*power4(snb) - 4*msb22*mu2*power2(sa)*power4(snb)))/(
     msb12*msb22) - (223*power4(invdct)*power5(mC2))/mt2 + (42*lmC2*power4(
     invdct)*power5(mC2))/mt2 + (42*lmt2*power4(invdct)*power5(mC2))/mt2 + (223
     *power2(cb)*power4(invdct)*power5(mC2))/mt2 - (42*lmC2*power2(cb)*power4(
     invdct)*power5(mC2))/mt2 - (42*lmt2*power2(cb)*power4(invdct)*power5(mC2))
     /mt2 - 213*power5(invdct)*power5(mC2) + 66*lmC2*power5(invdct)*power5(mC2)
     + 66*lmt2*power5(invdct)*power5(mC2) + 213*power2(cb)*power5(invdct)*
     power5(mC2) - 66*lmC2*power2(cb)*power5(invdct)*power5(mC2) - 66*lmt2*
     power2(cb)*power5(invdct)*power5(mC2) - (223*power2(cb)*power4(invdtw)*
     power5(mw2))/mt2 + (42*lmt2*power2(cb)*power4(invdtw)*power5(mw2))/mt2 + (
     42*lmw2*power2(cb)*power4(invdtw)*power5(mw2))/mt2 + 213*power2(cb)*power5
     (invdtw)*power5(mw2) - 66*lmt2*power2(cb)*power5(invdtw)*power5(mw2) - 66*
     lmw2*power2(cb)*power5(invdtw)*power5(mw2) - 213*power6(invdct)*power6(mC2
     ) + 66*lmC2*power6(invdct)*power6(mC2) + 66*lmt2*power6(invdct)*power6(mC2
     ) + 213*power2(cb)*power6(invdct)*power6(mC2) - 66*lmC2*power2(cb)*power6(
     invdct)*power6(mC2) - 66*lmt2*power2(cb)*power6(invdct)*power6(mC2) - 213*
     power2(cb)*power6(invdtw)*power6(mw2) + 66*lmt2*power2(cb)*power6(invdtw)*
     power6(mw2) + 66*lmw2*power2(cb)*power6(invdtw)*power6(mw2) + (213*power6(
     invdct)*power7(mC2))/mt2 - (66*lmC2*power6(invdct)*power7(mC2))/mt2 - (66*
     lmt2*power6(invdct)*power7(mC2))/mt2 - (213*power2(cb)*power6(invdct)*
     power7(mC2))/mt2 + (66*lmC2*power2(cb)*power6(invdct)*power7(mC2))/mt2 + (
     66*lmt2*power2(cb)*power6(invdct)*power7(mC2))/mt2 + (213*power2(cb)*
     power6(invdtw)*power7(mw2))/mt2 - (66*lmt2*power2(cb)*power6(invdtw)*
     power7(mw2))/mt2 - (66*lmw2*power2(cb)*power6(invdtw)*power7(mw2))/mt2)/
     216.;

   return power2(g3) * power2(yb) * result * twoLoop;
}

std::ostream& operator<<(std::ostream& out, const Parameters& pars)
{
   out <<
      "Delta alpha_s 2L parameters:\n"
      "g3   = " <<  pars.g3   << '\n' <<
      "yt   = " <<  pars.yt   << '\n' <<
      "yb   = " <<  pars.yb   << '\n' <<
      "mt   = " <<  pars.mt   << '\n' <<
      "mb   = " <<  pars.mb   << '\n' <<
      "mg   = " <<  pars.mg   << '\n' <<
      "mst1 = " <<  pars.mst1 << '\n' <<
      "mst2 = " <<  pars.mst2 << '\n' <<
      "msb1 = " <<  pars.msb1 << '\n' <<
      "msb2 = " <<  pars.msb2 << '\n' <<
      "msd1 = " <<  pars.msd1 << '\n' <<
      "msd2 = " <<  pars.msd2 << '\n' <<
      "xt   = " <<  pars.xt   << '\n' <<
      "xb   = " <<  pars.xb   << '\n' <<
      "mw   = " <<  pars.mw   << '\n' <<
      "mz   = " <<  pars.mz   << '\n' <<
      "mh   = " <<  pars.mh   << '\n' <<
      "mH   = " <<  pars.mH   << '\n' <<
      "mC   = " <<  pars.mC   << '\n' <<
      "mA   = " <<  pars.mA   << '\n' <<
      "mu   = " <<  pars.mu   << '\n' <<
      "tb   = " <<  pars.tb   << '\n' <<
      "Q    = " <<  pars.Q    << '\n';

   return out;
}

} // namespace mssm_twoloop_as
} // namespace flexiblesusy
