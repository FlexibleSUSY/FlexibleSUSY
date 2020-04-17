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

// This file has been generated at Sat 18 Apr 2020 01:20:19
// with the script "bquark_to_cpp.m".

#include "mssm_twoloop_mb.hpp"
#include "dilog.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>

namespace flexiblesusy {
namespace mssm_twoloop_mb {

namespace {
   const double Pi  = 3.1415926535897932384626433832795;
   const double zt2 = 1.6449340668482264364724151666460;

   template <typename T> T pow2(T x) noexcept { return x*x; }
   template <typename T> T pow3(T x) noexcept { return x*x*x; }
   template <typename T> T pow4(T x) noexcept { return pow2(pow2(x)); }
   template <typename T> T pow5(T x) noexcept { return x*pow4(x); }

   const double oneLoop = 1/pow2(4*Pi);
   const double twoLoop = pow2(oneLoop);

   template <typename T>
   bool is_zero(T a, T prec = std::numeric_limits<T>::epsilon()) noexcept
   {
      return std::abs(a) < prec;
   }

   template <typename T>
   bool is_equal(T a, T b, T prec = std::numeric_limits<T>::epsilon()) noexcept
   {
      return is_zero(a - b, prec);
   }

   /**
    * Fin20[] function from twoloopbubble.m .
    *
    * @param mm1 squared mass \f$m_1^2\f$
    * @param mm2 squared mass \f$m_2^2\f$
    * @param mmu squared renormalization scale
    *
    * @return Fin20(m12, m22, mmu)
    */
   double Fin20(double mm1, double mm2, double mmu) noexcept
   {
      const double log12 = std::log(mm1/mm2);
      const double log1u = std::log(mm1/mmu);
      const double log2u = std::log(mm2/mmu);

      return (6*(mm1*log1u + mm2*log2u) +
         (-mm1 - mm2)*(7 + pow2(Pi)/6) +
         (mm1 - mm2)*(2*dilog(1 - mm1/mm2) +
            pow2(log12)/2) +
         ((mm1 + mm2)*pow2(log12))/2 -
         2*(mm1*pow2(log1u) + mm2*pow2(log2u)))/2;
   }

   double LambdaSquared(double x, double y) noexcept
   {
      return pow2(1 - x - y) - 4*x*y;
   }

   /// ClausenCl[2,x]
   double ClausenCl2(double x) noexcept
   {
      const std::complex<double> img(0.0, 1.0);

      return std::imag(dilog(std::exp(img*x)));
   }

   /// x < 1 && y < 1, LambdaSquared(x,y) > 0
   double PhiPos(double x, double y) noexcept
   {
      const double lambda = std::sqrt(LambdaSquared(x,y));

      return (-(std::log(x)*std::log(y))
              + 2*std::log((1 - lambda + x - y)/2)*std::log((1 - lambda - x + y)/2)
              - 2*dilog((1 - lambda + x - y)/2)
              - 2*dilog((1 - lambda - x + y)/2)
              + pow2(Pi)/3)/lambda;
   }

   /// LambdaSquared(x,y) < 0
   double PhiNeg(double x, double y) noexcept
   {
      const double lambda = std::sqrt(-LambdaSquared(x,y));

      return 2*(+ ClausenCl2(2*std::acos((1 + x - y)/(2*std::sqrt(x))))
                + ClausenCl2(2*std::acos((1 - x + y)/(2*std::sqrt(y))))
                + ClausenCl2(2*std::acos((-1 + x + y)/(2*std::sqrt(x*y)))))/lambda;
   }

   double Phi(double x, double y) noexcept
   {
      const double lambda = LambdaSquared(x,y);

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
   double Fin3(double mm1, double mm2, double mm3, double mmu) noexcept
   {
      std::array<double,3> masses = { mm1, mm2, mm3 };
      std::sort(masses.begin(), masses.end());

      const double mm = masses[2];
      const double x = masses[0]/mm;
      const double y = masses[1]/mm;
      const double lambda = LambdaSquared(x,y);
      const double logx = std::log(x);
      const double logy = std::log(y);
      const double logm = std::log(mm/mmu);

      if (is_zero(lambda, 1e-10)) {
         return -(mm*(2*y*(-3 + 2*logm)*logy
                      + logx*(2*x*(-3 + 2*logm) + (-1 + x + y)*logy)
                      + (1 + x + y)*(7 - 6*logm + pow2(Pi)/6 + 2*pow2(logm))
                      + x*pow2(logx) + y*pow2(logy)))/2;
      }

      return mm*((-7 + 6*logm + logx*logy
                  - lambda*Phi(x,y) - pow2(Pi)/6 - 2*pow2(logm))/2
                 - (x*(7 - 6*logm + logx*(-6 + 4*logm + logy)
                       + pow2(Pi)/6 + 2*pow2(logm) + pow2(logx)))/2
                 - (y*(7 - 6*logm + (
                     -6 + 4*logm + logx)*logy + pow2(Pi)/6
                       + 2*pow2(logm) + pow2(logy)))/2);
   }

   /// Delta[m1,m2,m3,-1]
   double DeltaInv(double m1, double m2, double m3) noexcept
   {
      return 1/(pow2(m1) + pow2(m2) + pow2(m3) - 2*(m1*m2 + m1*m3 + m2*m3));
   }

} // anonymous namespace

/**
 * The function returns the 2-loop SQCD (QCD + SUSY) relation between
 * the Standard Model MS-bar bottom mass
 * \f$m_b^{\text{SM},\overline{\text{MS}}}\f$ and the MSSM DR-bar
 * bottom mass \f$m_b^{\text{MSSM},\overline{\text{DR}}}\f$,
 * Eq. (61) of [arxiv:0707.0650].  The relation has the form
 *
 * \f[
    m_b^{\text{SM},\overline{\text{MS}}} =
    m_b^{\text{MSSM},\overline{\text{DR}}} \left[
       1 + \left(\frac{\Delta m_b}{m_b}\right)_{1L}
         + \left(\frac{\Delta m_b}{m_b}\right)_{2L}
    \right]
   \f]
 *
 * The function returns \f$(\Delta m_b/m_b)_{2L}\f$.
 */
double delta_mb_2loop(const Parameters& pars)
{
   const double g3      = pars.g3;
   const double Xt      = pars.xt;
   const double Xb      = pars.xb;
   const double mgl     = pars.mg;
   const double mt      = pars.mt;
   const double mb      = pars.mb;
   const double mmt     = pow2(pars.mt);
   const double mmgl    = pow2(pars.mg);
   const double mmgl2   = pow2(mmgl);
   const double mmu     = pow2(pars.Q);
   const double mmst1   = pow2(pars.mst1);
   const double mmst2   = pow2(pars.mst2);
   const double mmsb1   = pow2(pars.msb1);
   const double mmsb2   = pow2(pars.msb2);
   const double mmsb12  = pow2(mmsb1);
   const double mmsb22  = pow2(mmsb2);
   const double mmsusy  = pow2(pars.msusy);
   const double lgu     = std::log(mmgl/mmu);
   const double lt1u    = std::log(mmst1/mmu);
   const double lt2u    = std::log(mmst2/mmu);
   const double lb1u    = std::log(mmsb1/mmu);
   const double lb2u    = std::log(mmsb2/mmu);
   const double lsu     = std::log(mmsusy/mmu);
   const double ltu     = std::log(mmt/mmu);
   const double s2t     = 2*mt*Xt / (mmst1 - mmst2);
   const double s2b     = 2*mb*Xb / (mmsb1 - mmsb2);
   const double invdgb1 = 1/(-mmgl + mmsb1);
   const double invdgb2 = 1/(-mmgl + mmsb2);
   const double invdb12 = 1/(mmsb1 - mmsb2);

   const double result =
   Fin20(mmsb1,mmsusy,mmu)*((32*mmgl*(mmsb1 - mmsusy)*s2b*pow2(invdgb1))/(3.*mb
     *mgl) + (8*(2*invdgb1 + (-7*mmsb1 + 3*mmsusy)*pow2(invdgb1) + 4*(mmsb12 -
     mmsb1*mmsusy)*pow3(invdgb1)))/3.) + Fin20(mmsb2,mmgl,mmu)*((-4*mmgl*s2b*(4
     *invdgb1 - 22*invdgb2 + 3*(mmsb1 - mmsb2)*pow2(invdgb1)))/(9.*mb*mgl) + (
     128*invdb12 - invdgb1 - 37*invdgb2 - 8*invdb12*invdgb1*mmsb2 - 120*invdb12
     *invdgb2*mmsb2 + 108*mmsb2*pow2(invdgb2) - 128*invdb12*pow2(s2b) - 15*
     invdgb1*pow2(s2b) + 15*invdgb2*pow2(s2b) - 26*invdb12*invdgb1*mmsb2*pow2(
     s2b) + 154*invdb12*invdgb2*mmsb2*pow2(s2b) + pow2(invdgb1)*(9*mmsb1 - 9*
     mmsb2 - 26*mmsb1*pow2(s2b)) - 12*mmsb12*pow3(invdgb1) + 12*mmsb1*mmsb2*
     pow3(invdgb1) - 48*mmsb22*pow3(invdgb2))/9.) + Fin20(mmsb1,mmgl,mmu)*((-4*
     mmgl*s2b*(22*invdgb1 - 4*invdgb2 + 3*(mmsb1 - mmsb2)*pow2(invdgb2)))/(9.*
     mb*mgl) + (-128*invdb12 + 83*invdgb1 + 7*invdgb2 + 120*invdb12*invdgb1*
     mmsb2 + 8*invdb12*invdgb2*mmsb2 + 108*mmsb1*pow2(invdgb1) + 128*invdb12*
     pow2(s2b) - 139*invdgb1*pow2(s2b) + 11*invdgb2*pow2(s2b) - 154*invdb12*
     invdgb1*mmsb2*pow2(s2b) + 26*invdb12*invdgb2*mmsb2*pow2(s2b) + pow2(
     invdgb2)*(-9*mmsb1 + 9*mmsb2 - 26*mmsb2*pow2(s2b)) - 48*mmsb12*pow3(
     invdgb1) + 12*mmsb1*mmsb2*pow3(invdgb2) - 12*mmsb22*pow3(invdgb2))/9.) +
     Fin20(mmsb1,mmsb2,mmu)*((4*mmgl*s2b*(invdgb1 - invdgb2 + 3*(mmsb1 - mmsb2)
     *pow2(invdgb1) + 3*(mmsb1 - mmsb2)*pow2(invdgb2)))/(9.*mb*mgl) + (10*
     invdgb1 + 2*invdgb2 + 8*invdb12*invdgb1*mmsb2 - 8*invdb12*invdgb2*mmsb2 -
     11*invdgb1*pow2(s2b) - 11*invdgb2*pow2(s2b) + pow2(invdgb1)*(-21*mmsb1 + 9
     *mmsb2 + 26*mmsb1*pow2(s2b)) + pow2(invdgb2)*(9*mmsb1 - 21*mmsb2 + 26*
     mmsb2*pow2(s2b)) + 12*mmsb12*pow3(invdgb1) - 12*mmsb1*mmsb2*pow3(invdgb1)
     - 12*mmsb1*mmsb2*pow3(invdgb2) + 12*mmsb22*pow3(invdgb2))/9.) + Fin20(mmgl
     ,mmsusy,mmu)*((-32*mmgl*s2b*(invdgb1 - invdgb2 + (mmsb1 - mmsusy)*pow2(
     invdgb1) + (-mmsb2 + mmsusy)*pow2(invdgb2)))/(3.*mb*mgl) + (8*(invdgb1 +
     invdgb2 + 3*(mmsb1 - mmsusy)*pow2(invdgb1) + 3*(mmsb2 - mmsusy)*pow2(
     invdgb2) - 4*mmsb12*pow3(invdgb1) + 4*mmsb1*mmsusy*pow3(invdgb1) - 4*
     mmsb22*pow3(invdgb2) + 4*mmsb2*mmsusy*pow3(invdgb2)))/3.) + Fin20(mmsb2,
     mmsusy,mmu)*((32*mmgl*(-mmsb2 + mmsusy)*s2b*pow2(invdgb2))/(3.*mb*mgl) + (
     8*(2*invdgb2 + (-7*mmsb2 + 3*mmsusy)*pow2(invdgb2) + 4*(mmsb22 - mmsb2*
     mmsusy)*pow3(invdgb2)))/3.) + Fin3(mmt,mmsb1,mmst1,mmu)*(-((mmsb1 - mmst1
     + mmt)*pow2(invdgb1)) + (4*mmgl*(mmsb1 - mmst1 + mmt)*s2b*pow2(invdgb1))/(
     3.*mb*mgl) + (4*(mmsb12 + mmsb1*(-mmst1 + mmt))*pow3(invdgb1))/3. + s2t*((
     4*mmt*s2b*(invdgb1 - 2*mmsb1*pow2(invdgb1)))/(3.*mb*mt) + (4*mmgl*mmt*(
     pow2(invdgb1) - 2*mmsb1*pow3(invdgb1)))/(3.*mgl*mt)) + DeltaInv(mmt,mmsb1,
     mmst1)*((4*mmgl*s2t*pow2(invdgb1)*((mmsb12 - mmsb1*mmst1)*mmt - mmsb1*pow2
     (mmt)))/(3.*mgl*mt) + (2*(invdgb1*(mmsb12 - 2*mmsb1*mmst1 - mmsb1*mmt -
     mmst1*mmt + pow2(mmst1)) + 2*pow2(invdgb1)*(2*mmsb12*mmst1 + mmsb12*mmt +
     mmsb1*mmst1*mmt - mmsb1*pow2(mmst1) - pow3(mmsb1))))/3.)) + Fin3(mmt,mmsb1
     ,mmst2,mmu)*(-((mmsb1 - mmst2 + mmt)*pow2(invdgb1)) + (4*mmgl*(mmsb1 -
     mmst2 + mmt)*s2b*pow2(invdgb1))/(3.*mb*mgl) + (4*(mmsb12 + mmsb1*(-mmst2 +
     mmt))*pow3(invdgb1))/3. + s2t*((-4*mmt*s2b*(invdgb1 - 2*mmsb1*pow2(invdgb1
     )))/(3.*mb*mt) - (4*mmgl*mmt*(pow2(invdgb1) - 2*mmsb1*pow3(invdgb1)))/(3.*
     mgl*mt)) + DeltaInv(mmt,mmsb1,mmst2)*((4*mmgl*s2t*pow2(invdgb1)*(-(mmsb12*
     mmt) + mmsb1*mmst2*mmt + mmsb1*pow2(mmt)))/(3.*mgl*mt) + (2*(invdgb1*(
     mmsb12 - 2*mmsb1*mmst2 - mmsb1*mmt - mmst2*mmt + pow2(mmst2)) + 2*pow2(
     invdgb1)*(2*mmsb12*mmst2 + mmsb12*mmt + mmsb1*mmst2*mmt - mmsb1*pow2(mmst2
     ) - pow3(mmsb1))))/3.)) + Fin3(mmt,mmsb2,mmst1,mmu)*(-((mmsb2 - mmst1 +
     mmt)*pow2(invdgb2)) - (4*mmgl*(mmsb2 - mmst1 + mmt)*s2b*pow2(invdgb2))/(3.
     *mb*mgl) + (4*(mmsb22 + mmsb2*(-mmst1 + mmt))*pow3(invdgb2))/3. + s2t*((-4
     *mmt*s2b*(invdgb2 - 2*mmsb2*pow2(invdgb2)))/(3.*mb*mt) + (4*mmgl*mmt*(pow2
     (invdgb2) - 2*mmsb2*pow3(invdgb2)))/(3.*mgl*mt)) + DeltaInv(mmt,mmsb2,
     mmst1)*((4*mmgl*s2t*pow2(invdgb2)*((mmsb22 - mmsb2*mmst1)*mmt - mmsb2*pow2
     (mmt)))/(3.*mgl*mt) + (2*(invdgb2*(mmsb22 - 2*mmsb2*mmst1 - mmsb2*mmt -
     mmst1*mmt + pow2(mmst1)) + 2*pow2(invdgb2)*(2*mmsb22*mmst1 + mmsb22*mmt +
     mmsb2*mmst1*mmt - mmsb2*pow2(mmst1) - pow3(mmsb2))))/3.)) + Fin3(mmt,mmsb2
     ,mmst2,mmu)*(-((mmsb2 - mmst2 + mmt)*pow2(invdgb2)) - (4*mmgl*(mmsb2 -
     mmst2 + mmt)*s2b*pow2(invdgb2))/(3.*mb*mgl) + (4*(mmsb22 + mmsb2*(-mmst2 +
     mmt))*pow3(invdgb2))/3. + s2t*((4*mmt*s2b*(invdgb2 - 2*mmsb2*pow2(invdgb2)
     ))/(3.*mb*mt) - (4*mmgl*mmt*(pow2(invdgb2) - 2*mmsb2*pow3(invdgb2)))/(3.*
     mgl*mt)) + DeltaInv(mmt,mmsb2,mmst2)*((4*mmgl*s2t*pow2(invdgb2)*(-(mmsb22*
     mmt) + mmsb2*mmst2*mmt + mmsb2*pow2(mmt)))/(3.*mgl*mt) + (2*(invdgb2*(
     mmsb22 - 2*mmsb2*mmst2 - mmsb2*mmt - mmst2*mmt + pow2(mmst2)) + 2*pow2(
     invdgb2)*(2*mmsb22*mmst2 + mmsb22*mmt + mmsb2*mmst2*mmt - mmsb2*pow2(mmst2
     ) - pow3(mmsb2))))/3.)) - (2*mmgl*s2b*(-32*lb1u + 32*lb2u + 32*lb1u*lgu -
     32*lb2u*lgu + 162*invdgb1*mmsb1 - 36*invdgb2*mmsb1 + 28*invdgb1*lb1u*mmsb1
      + 24*invdgb2*lb1u*mmsb1 - 216*invdgb1*lgu*mmsb1 + 12*invdgb2*lgu*mmsb1 +
     50*invdgb1*lb1u*lgu*mmsb1 - 6*invdgb2*lb1u*lgu*mmsb1 + 72*invdgb1*ltu*
     mmsb1 - 24*invdgb1*lgu*ltu*mmsb1 + 36*invdgb1*mmsb2 - 162*invdgb2*mmsb2 -
     4*invdgb1*lb1u*mmsb2 + 4*invdgb2*lb1u*mmsb2 - 20*invdgb1*lb2u*mmsb2 - 32*
     invdgb2*lb2u*mmsb2 + 12*invdgb1*lb1u*lb2u*mmsb2 + 20*invdgb2*lb1u*lb2u*
     mmsb2 - 12*invdgb1*lgu*mmsb2 + 216*invdgb2*lgu*mmsb2 + 20*invdgb1*lb1u*lgu
     *mmsb2 - 20*invdgb2*lb1u*lgu*mmsb2 - 6*invdgb1*lb2u*lgu*mmsb2 - 38*invdgb2
     *lb2u*lgu*mmsb2 - 72*invdgb2*ltu*mmsb2 + 24*invdgb2*lgu*ltu*mmsb2 - 4*
     invdb12*invdgb1*lb1u*mmsb22 + 4*invdb12*invdgb2*lb1u*mmsb22 + 4*invdb12*
     invdgb1*lb2u*mmsb22 - 4*invdb12*invdgb2*lb2u*mmsb22 + 4*invdb12*invdgb1*
     lb1u*lb2u*mmsb22 + 28*invdb12*invdgb2*lb1u*lb2u*mmsb22 + 28*invdb12*
     invdgb1*lb1u*lgu*mmsb22 - 28*invdb12*invdgb2*lb1u*lgu*mmsb22 - 4*invdb12*
     invdgb1*lb2u*lgu*mmsb22 + 4*invdb12*invdgb2*lb2u*lgu*mmsb22 - 48*invdgb1*
     mmst1 + 48*invdgb2*mmst1 - 12*invdgb1*lgu*mmst1 + 12*invdgb2*lgu*mmst1 +
     12*invdgb1*lt1u*mmst1 - 12*invdgb2*lt1u*mmst1 + 6*invdgb1*lgu*lt1u*mmst1 -
     6*invdgb2*lgu*lt1u*mmst1 + 36*invdgb1*ltu*mmst1 - 36*invdgb2*ltu*mmst1 -
     12*invdgb1*lt1u*ltu*mmst1 + 12*invdgb2*lt1u*ltu*mmst1 - 48*invdgb1*mmst2 +
     48*invdgb2*mmst2 - 12*invdgb1*lgu*mmst2 + 12*invdgb2*lgu*mmst2 + 12*
     invdgb1*lt2u*mmst2 - 12*invdgb2*lt2u*mmst2 + 6*invdgb1*lgu*lt2u*mmst2 - 6*
     invdgb2*lgu*lt2u*mmst2 + 36*invdgb1*ltu*mmst2 - 36*invdgb2*ltu*mmst2 - 12*
     invdgb1*lt2u*ltu*mmst2 + 12*invdgb2*lt2u*ltu*mmst2 + 288*invdgb1*mmsusy -
     288*invdgb2*mmsusy - 96*invdgb1*lgu*mmsusy + 96*invdgb2*lgu*mmsusy - 192*
     invdgb1*lsu*mmsusy + 192*invdgb2*lsu*mmsusy + 48*invdgb1*lgu*lsu*mmsusy -
     48*invdgb2*lgu*lsu*mmsusy - 72*invdgb1*mmt + 72*invdgb2*mmt + 24*invdgb1*
     lgu*mmt - 24*invdgb2*lgu*mmt + 48*invdgb1*ltu*mmt - 48*invdgb2*ltu*mmt -
     12*invdgb1*lgu*ltu*mmt + 12*invdgb2*lgu*ltu*mmt + 20*invdgb1*mmsb1*zt2 - 6
     *invdgb2*mmsb1*zt2 + 6*invdgb1*mmsb2*zt2 - 20*invdgb2*mmsb2*zt2 - 6*
     invdgb1*mmst1*zt2 + 6*invdgb2*mmst1*zt2 - 6*invdgb1*mmst2*zt2 + 6*invdgb2*
     mmst2*zt2 + 48*invdgb1*mmsusy*zt2 - 48*invdgb2*mmsusy*zt2 - 12*invdgb1*mmt
     *zt2 + 12*invdgb2*mmt*zt2 - 11*invdgb1*mmsb1*pow2(lb1u) - 3*invdgb2*mmsb1*
     pow2(lb1u) - 16*invdgb1*mmsb2*pow2(lb1u) - 16*invdb12*invdgb1*mmsb22*pow2(
     lb1u) + 3*invdgb1*mmsb2*pow2(lb2u) - 5*invdgb2*mmsb2*pow2(lb2u) - 16*
     invdb12*invdgb2*mmsb22*pow2(lb2u) + 37*invdgb1*mmsb1*pow2(lgu) - 3*invdgb2
     *mmsb1*pow2(lgu) - invdgb1*mmsb2*pow2(lgu) - 33*invdgb2*mmsb2*pow2(lgu) -
     12*invdb12*invdgb1*mmsb22*pow2(lgu) + 12*invdb12*invdgb2*mmsb22*pow2(lgu)
     + 3*invdgb1*mmst1*pow2(lgu) - 3*invdgb2*mmst1*pow2(lgu) + 3*invdgb1*mmst2*
     pow2(lgu) - 3*invdgb2*mmst2*pow2(lgu) + 24*invdgb1*mmsusy*pow2(lgu) - 24*
     invdgb2*mmsusy*pow2(lgu) - 6*invdgb1*mmt*pow2(lgu) + 6*invdgb2*mmt*pow2(
     lgu) + 24*invdgb1*mmsusy*pow2(lsu) - 24*invdgb2*mmsusy*pow2(lsu) - 3*
     invdgb1*mmst1*pow2(lt1u) + 3*invdgb2*mmst1*pow2(lt1u) - 3*invdgb1*mmst2*
     pow2(lt2u) + 3*invdgb2*mmst2*pow2(lt2u) - 12*invdgb1*mmsb1*pow2(ltu) + 12*
     invdgb2*mmsb2*pow2(ltu) - 6*invdgb1*mmst1*pow2(ltu) + 6*invdgb2*mmst1*pow2
     (ltu) - 6*invdgb1*mmst2*pow2(ltu) + 6*invdgb2*mmst2*pow2(ltu) - 6*invdgb1*
     mmt*pow2(ltu) + 6*invdgb2*mmt*pow2(ltu) + 8*lb1u*lgu*pow3(invdgb1)*pow3(
     mmsb1) - 4*pow2(lb1u)*pow3(invdgb1)*pow3(mmsb1) - 4*pow2(lgu)*pow3(invdgb1
     )*pow3(mmsb1) - 8*invdgb1*lb1u*lb2u*pow2(invdb12)*pow3(mmsb2) + 8*invdgb2*
     lb1u*lb2u*pow2(invdb12)*pow3(mmsb2) + 8*invdgb1*lb1u*lgu*pow2(invdb12)*
     pow3(mmsb2) - 8*invdgb2*lb1u*lgu*pow2(invdb12)*pow3(mmsb2) + 8*invdgb1*
     lb2u*lgu*pow2(invdb12)*pow3(mmsb2) - 8*invdgb2*lb2u*lgu*pow2(invdb12)*pow3
     (mmsb2) - 8*invdgb1*pow2(invdb12)*pow2(lgu)*pow3(mmsb2) + 8*invdgb2*pow2(
     invdb12)*pow2(lgu)*pow3(mmsb2) - 8*lb2u*lgu*pow3(invdgb2)*pow3(mmsb2) + 4*
     pow2(lb2u)*pow3(invdgb2)*pow3(mmsb2) + 4*pow2(lgu)*pow3(invdgb2)*pow3(
     mmsb2) + pow2(invdgb2)*(2*lgu*(3*(-2 + lb1u)*mmsb1*mmsb2 + 2*(45 + lb1u)*
     mmsb22 + 3*mmsb2*((-2 + lt1u)*mmst1 + (-2 + lt2u)*mmst2 - 16*mmsusy + 8*
     lsu*mmsusy + 4*mmt - 2*ltu*mmt)) - 2*lb2u*(3*(-2 + lb1u)*mmsb1*mmsb2 + (90
      + 2*lb1u - 39*lgu)*mmsb22 + 3*mmsb2*((-2 + lt1u)*mmst1 + (-2 + lt2u)*
     mmst2 - 16*mmsusy + 8*lsu*mmsusy + 4*mmt - 2*ltu*mmt)) - (3*mmsb1*mmsb2 -
     31*mmsb22 + 3*mmsb2*(mmst1 + mmst2 + 8*mmsusy - 2*mmt))*pow2(lb2u) - 4*
     invdb12*(lb1u*(lb2u - lgu) - lb2u*lgu)*pow3(mmsb2) + pow2(lgu)*(3*mmsb1*
     mmsb2 - 109*mmsb22 + 3*mmsb2*mmst1 + 3*mmsb2*mmst2 + 24*mmsb2*mmsusy - 6*
     mmsb2*mmt - 4*invdb12*pow3(mmsb2))) + pow2(invdgb1)*((-31*mmsb12 + 3*mmsb1
     *(mmsb2 + mmst1 + mmst2 + 8*mmsusy - 2*mmt))*pow2(lb1u) + pow2(lgu)*(105*
     mmsb12 - 7*mmsb1*mmsb2 - 4*mmsb22 - 3*mmsb1*mmst1 - 3*mmsb1*mmst2 - 24*
     mmsb1*mmsusy + 6*mmsb1*mmt - 4*invdb12*pow3(mmsb2)) - 2*(lgu*(90*mmsb12 +
     (-6 + lb2u)*mmsb1*mmsb2 - 2*lb2u*mmsb22 + 3*mmsb1*((-2 + lt1u)*mmst1 + (-2
      + lt2u)*mmst2 - 16*mmsusy + 8*lsu*mmsusy + 4*mmt - 2*ltu*mmt)) + lb1u*((-
     90 + 37*lgu)*mmsb12 + 2*(lb2u - lgu)*mmsb22 - mmsb1*((-6 + lb2u + 2*lgu)*
     mmsb2 + 3*((-2 + lt1u)*mmst1 + (-2 + lt2u)*mmst2 - 16*mmsusy + 8*lsu*
     mmsusy + 4*mmt - 2*ltu*mmt))) + 2*invdb12*(lb1u*(lb2u - lgu) - lb2u*lgu)*
     pow3(mmsb2)))))/(9.*mb*mgl) + Fin3(mmt,mmst1,mmgl,mmu)*((4*mmgl*s2b*(
     invdgb1 - invdgb2 - (mmsb1 - mmst1 + mmt)*pow2(invdgb1) + (mmsb2 - mmst1 +
     mmt)*pow2(invdgb2)))/(3.*mb*mgl) + (-3*invdgb1 - 3*invdgb2 + (7*mmsb1 - 3*
     mmst1 + 3*mmt)*pow2(invdgb1) + (7*mmsb2 - 3*mmst1 + 3*mmt)*pow2(invdgb2) -
     4*mmsb12*pow3(invdgb1) + 4*mmsb1*mmst1*pow3(invdgb1) - 4*mmsb1*mmt*pow3(
     invdgb1) - 4*mmsb22*pow3(invdgb2) + 4*mmsb2*mmst1*pow3(invdgb2) - 4*mmsb2*
     mmt*pow3(invdgb2))/3. + s2t*((4*mmt*s2b*(-invdgb1 + invdgb2 + 2*mmsb1*pow2
     (invdgb1) - 2*mmsb2*pow2(invdgb2)))/(3.*mb*mt) - (4*mmgl*mmt*(pow2(invdgb1
     ) + pow2(invdgb2) - 2*(mmsb1*pow3(invdgb1) + mmsb2*pow3(invdgb2))))/(3.*
     mgl*mt)) + DeltaInv(mmt,mmst1,mmgl)*((8*mmgl*s2b*(mmsb1 - invdgb1*mmsb12 -
     mmsb2 + invdgb2*mmsb22 + 2*invdgb1*mmsb1*mmst1 - 2*invdgb2*mmsb2*mmst1 +
     invdgb1*mmsb1*mmt - invdgb2*mmsb2*mmt + invdgb1*mmst1*mmt - invdgb2*mmst1*
     mmt + (-invdgb1 + invdgb2)*pow2(mmst1)))/(3.*mb*mgl) + s2t*((-8*s2b*((
     mmsb1 - invdgb1*mmsb12 - mmsb2 + invdgb2*mmsb22 + invdgb1*mmsb1*mmst1 -
     invdgb2*mmsb2*mmst1)*mmt + (invdgb1*mmsb1 - invdgb2*mmsb2)*pow2(mmt)))/(3.
     *mb*mt) - (4*mmgl*(-2*mmt + 2*invdgb1*mmsb1*mmt + 2*invdgb2*mmsb2*mmt -
     invdgb1*mmst1*mmt - invdgb2*mmst1*mmt - invdgb1*pow2(mmt) - invdgb2*pow2(
     mmt) + pow2(invdgb1)*(-(mmsb12*mmt) + mmsb1*mmst1*mmt + mmsb1*pow2(mmt)) +
     pow2(invdgb2)*(-(mmsb22*mmt) + mmsb2*mmst1*mmt + mmsb2*pow2(mmt))))/(3.*
     mgl*mt)) - (4*(2*mmgl + 2*mmsb1 - 3*invdgb1*mmsb12 + 2*mmsb2 - 3*invdgb2*
     mmsb22 - 4*mmst1 + 4*invdgb1*mmsb1*mmst1 + 4*invdgb2*mmsb2*mmst1 - 2*mmt +
     2*invdgb1*mmsb1*mmt + 2*invdgb2*mmsb2*mmt + invdgb1*mmst1*mmt + invdgb2*
     mmst1*mmt - invdgb1*pow2(mmst1) - invdgb2*pow2(mmst1) + pow2(invdgb1)*(-2*
     mmsb12*mmst1 - mmsb12*mmt - mmsb1*mmst1*mmt + mmsb1*pow2(mmst1) + pow3(
     mmsb1)) + pow2(invdgb2)*(-2*mmsb22*mmst1 - mmsb22*mmt - mmsb2*mmst1*mmt +
     mmsb2*pow2(mmst1) + pow3(mmsb2))))/3.)) + Fin3(mmt,mmst2,mmgl,mmu)*((4*
     mmgl*s2b*(invdgb1 - invdgb2 - (mmsb1 - mmst2 + mmt)*pow2(invdgb1) + (mmsb2
      - mmst2 + mmt)*pow2(invdgb2)))/(3.*mb*mgl) + (-3*invdgb1 - 3*invdgb2 + (7
     *mmsb1 - 3*mmst2 + 3*mmt)*pow2(invdgb1) + (7*mmsb2 - 3*mmst2 + 3*mmt)*pow2
     (invdgb2) - 4*mmsb12*pow3(invdgb1) + 4*mmsb1*mmst2*pow3(invdgb1) - 4*mmsb1
     *mmt*pow3(invdgb1) - 4*mmsb22*pow3(invdgb2) + 4*mmsb2*mmst2*pow3(invdgb2)
     - 4*mmsb2*mmt*pow3(invdgb2))/3. + s2t*((4*mmt*s2b*(invdgb1 - invdgb2 - 2*
     mmsb1*pow2(invdgb1) + 2*mmsb2*pow2(invdgb2)))/(3.*mb*mt) + (4*mmgl*mmt*(
     pow2(invdgb1) + pow2(invdgb2) - 2*(mmsb1*pow3(invdgb1) + mmsb2*pow3(
     invdgb2))))/(3.*mgl*mt)) + DeltaInv(mmt,mmst2,mmgl)*((8*mmgl*s2b*(mmsb1 -
     invdgb1*mmsb12 - mmsb2 + invdgb2*mmsb22 + 2*invdgb1*mmsb1*mmst2 - 2*
     invdgb2*mmsb2*mmst2 + invdgb1*mmsb1*mmt - invdgb2*mmsb2*mmt + invdgb1*
     mmst2*mmt - invdgb2*mmst2*mmt + (-invdgb1 + invdgb2)*pow2(mmst2)))/(3.*mb*
     mgl) + s2t*((8*s2b*((mmsb1 - invdgb1*mmsb12 - mmsb2 + invdgb2*mmsb22 +
     invdgb1*mmsb1*mmst2 - invdgb2*mmsb2*mmst2)*mmt + (invdgb1*mmsb1 - invdgb2*
     mmsb2)*pow2(mmt)))/(3.*mb*mt) + (4*mmgl*(-2*mmt + 2*invdgb1*mmsb1*mmt + 2*
     invdgb2*mmsb2*mmt - invdgb1*mmst2*mmt - invdgb2*mmst2*mmt - invdgb1*pow2(
     mmt) - invdgb2*pow2(mmt) + pow2(invdgb1)*(-(mmsb12*mmt) + mmsb1*mmst2*mmt
     + mmsb1*pow2(mmt)) + pow2(invdgb2)*(-(mmsb22*mmt) + mmsb2*mmst2*mmt +
     mmsb2*pow2(mmt))))/(3.*mgl*mt)) - (4*(2*mmgl + 2*mmsb1 - 3*invdgb1*mmsb12
     + 2*mmsb2 - 3*invdgb2*mmsb22 - 4*mmst2 + 4*invdgb1*mmsb1*mmst2 + 4*invdgb2
     *mmsb2*mmst2 - 2*mmt + 2*invdgb1*mmsb1*mmt + 2*invdgb2*mmsb2*mmt + invdgb1
     *mmst2*mmt + invdgb2*mmst2*mmt - invdgb1*pow2(mmst2) - invdgb2*pow2(mmst2)
     + pow2(invdgb1)*(-2*mmsb12*mmst2 - mmsb12*mmt - mmsb1*mmst2*mmt + mmsb1*
     pow2(mmst2) + pow3(mmsb1)) + pow2(invdgb2)*(-2*mmsb22*mmst2 - mmsb22*mmt -
     mmsb2*mmst2*mmt + mmsb2*pow2(mmst2) + pow3(mmsb2))))/3.)) + DeltaInv(mmt,
     mmsb1,mmst1)*((2*mmgl*s2t*pow2(invdgb1)*(-14*mmsb12*mmst1*mmt + 6*lb1u*
     mmsb12*mmst1*mmt - 6*lt1u*mmsb12*mmst1*mmt + 2*lb1u*lt1u*mmsb12*mmst1*mmt
     + 12*ltu*mmsb12*mmst1*mmt - 4*lb1u*ltu*mmsb12*mmst1*mmt - 2*mmsb12*mmst1*
     mmt*zt2 - 2*mmsb12*mmst1*mmt*pow2(ltu) + 6*lt1u*mmsb1*mmt*pow2(mmst1) - 2*
     lb1u*lt1u*mmsb1*mmt*pow2(mmst1) - 6*ltu*mmsb1*mmt*pow2(mmst1) + 2*lb1u*ltu
     *mmsb1*mmt*pow2(mmst1) + mmsb1*mmt*pow2(ltu)*pow2(mmst1) - 14*mmsb12*pow2(
     mmt) + 6*lb1u*mmsb12*pow2(mmt) + 6*ltu*mmsb12*pow2(mmt) - 2*lb1u*ltu*
     mmsb12*pow2(mmt) - 28*mmsb1*mmst1*pow2(mmt) + 6*lt1u*mmsb1*mmst1*pow2(mmt)
     + 2*lb1u*lt1u*mmsb1*mmst1*pow2(mmt) + 18*ltu*mmsb1*mmst1*pow2(mmt) - 2*
     lb1u*ltu*mmsb1*mmst1*pow2(mmt) - 4*lt1u*ltu*mmsb1*mmst1*pow2(mmt) - 2*
     mmsb12*zt2*pow2(mmt) - 4*mmsb1*mmst1*zt2*pow2(mmt) - mmsb12*pow2(ltu)*pow2
     (mmt) - 3*mmsb1*mmst1*pow2(ltu)*pow2(mmt) + pow2(lt1u)*(mmsb12*mmst1*mmt -
     mmsb1*mmt*pow2(mmst1) - mmsb1*mmst1*pow2(mmt)) + 14*mmt*pow3(mmsb1) - 6*
     lb1u*mmt*pow3(mmsb1) - 6*ltu*mmt*pow3(mmsb1) + 2*lb1u*ltu*mmt*pow3(mmsb1)
     + 2*mmt*zt2*pow3(mmsb1) + mmt*pow2(ltu)*pow3(mmsb1) - pow2(lb1u)*(mmsb12*
     mmst1*mmt + mmsb12*pow2(mmt) - mmt*pow3(mmsb1))))/(3.*mgl*mt) + (-(invdgb1
     *(14*mmsb12*mmst1 - 12*lb1u*mmsb12*mmst1 + 6*lt1u*mmsb12*mmst1 - 6*ltu*
     mmsb12*mmst1 + 4*lb1u*ltu*mmsb12*mmst1 - 2*lt1u*ltu*mmsb12*mmst1 + 14*
     mmsb12*mmt - 6*lb1u*mmsb12*mmt - 6*ltu*mmsb12*mmt + 2*lb1u*ltu*mmsb12*mmt
     + 56*mmsb1*mmst1*mmt - 6*lb1u*mmsb1*mmst1*mmt - 6*lt1u*mmsb1*mmst1*mmt - 4
     *lb1u*lt1u*mmsb1*mmst1*mmt - 36*ltu*mmsb1*mmst1*mmt + 6*lb1u*ltu*mmsb1*
     mmst1*mmt + 6*lt1u*ltu*mmsb1*mmst1*mmt + 2*mmsb12*mmst1*zt2 + 2*mmsb12*mmt
     *zt2 + 8*mmsb1*mmst1*mmt*zt2 + mmsb12*mmst1*pow2(ltu) + mmsb12*mmt*pow2(
     ltu) + 6*mmsb1*mmst1*mmt*pow2(ltu) + 14*mmsb1*pow2(mmst1) + 6*lb1u*mmsb1*
     pow2(mmst1) - 12*lt1u*mmsb1*pow2(mmst1) - 6*ltu*mmsb1*pow2(mmst1) - 2*lb1u
     *ltu*mmsb1*pow2(mmst1) + 4*lt1u*ltu*mmsb1*pow2(mmst1) + 14*mmt*pow2(mmst1)
     - 6*lt1u*mmt*pow2(mmst1) - 6*ltu*mmt*pow2(mmst1) + 2*lt1u*ltu*mmt*pow2(
     mmst1) + 2*mmsb1*zt2*pow2(mmst1) + 2*mmt*zt2*pow2(mmst1) + mmsb1*pow2(ltu)
     *pow2(mmst1) + mmt*pow2(ltu)*pow2(mmst1) + pow2(lb1u)*(2*mmsb12*mmst1 +
     mmsb12*mmt + mmsb1*mmst1*mmt - mmsb1*pow2(mmst1) - pow3(mmsb1)) - 14*pow3(
     mmsb1) + 6*lb1u*pow3(mmsb1) + 6*ltu*pow3(mmsb1) - 2*lb1u*ltu*pow3(mmsb1) -
     2*zt2*pow3(mmsb1) - pow2(ltu)*pow3(mmsb1) + pow2(lt1u)*(-(mmsb12*mmst1) +
     mmsb1*mmst1*mmt + (2*mmsb1 + mmt)*pow2(mmst1) - pow3(mmst1)) - 14*pow3(
     mmst1) + 6*lt1u*pow3(mmst1) + 6*ltu*pow3(mmst1) - 2*lt1u*ltu*pow3(mmst1) -
     2*zt2*pow3(mmst1) - pow2(ltu)*pow3(mmst1))) - 2*pow2(invdgb1)*(-56*mmsb12*
     mmst1*mmt + 6*lb1u*mmsb12*mmst1*mmt + 6*lt1u*mmsb12*mmst1*mmt + 4*lb1u*
     lt1u*mmsb12*mmst1*mmt + 36*ltu*mmsb12*mmst1*mmt - 6*lb1u*ltu*mmsb12*mmst1*
     mmt - 6*lt1u*ltu*mmsb12*mmst1*mmt - 8*mmsb12*mmst1*mmt*zt2 - 6*mmsb12*
     mmst1*mmt*pow2(ltu) - 14*mmsb12*pow2(mmst1) - 6*lb1u*mmsb12*pow2(mmst1) +
     12*lt1u*mmsb12*pow2(mmst1) + 6*ltu*mmsb12*pow2(mmst1) + 2*lb1u*ltu*mmsb12*
     pow2(mmst1) - 4*lt1u*ltu*mmsb12*pow2(mmst1) - 14*mmsb1*mmt*pow2(mmst1) + 6
     *lt1u*mmsb1*mmt*pow2(mmst1) + 6*ltu*mmsb1*mmt*pow2(mmst1) - 2*lt1u*ltu*
     mmsb1*mmt*pow2(mmst1) - 2*mmsb12*zt2*pow2(mmst1) - 2*mmsb1*mmt*zt2*pow2(
     mmst1) - mmsb12*pow2(ltu)*pow2(mmst1) - mmsb1*mmt*pow2(ltu)*pow2(mmst1) -
     14*mmst1*pow3(mmsb1) + 12*lb1u*mmst1*pow3(mmsb1) - 6*lt1u*mmst1*pow3(mmsb1
     ) + 6*ltu*mmst1*pow3(mmsb1) - 4*lb1u*ltu*mmst1*pow3(mmsb1) + 2*lt1u*ltu*
     mmst1*pow3(mmsb1) - 14*mmt*pow3(mmsb1) + 6*lb1u*mmt*pow3(mmsb1) + 6*ltu*
     mmt*pow3(mmsb1) - 2*lb1u*ltu*mmt*pow3(mmsb1) - 2*mmst1*zt2*pow3(mmsb1) - 2
     *mmt*zt2*pow3(mmsb1) - mmst1*pow2(ltu)*pow3(mmsb1) - mmt*pow2(ltu)*pow3(
     mmsb1) + 14*mmsb1*pow3(mmst1) - 6*lt1u*mmsb1*pow3(mmst1) - 6*ltu*mmsb1*
     pow3(mmst1) + 2*lt1u*ltu*mmsb1*pow3(mmst1) + 2*mmsb1*zt2*pow3(mmst1) +
     mmsb1*pow2(ltu)*pow3(mmst1) + pow2(lt1u)*(-(mmsb12*mmst1*mmt) - (2*mmsb12
     + mmsb1*mmt)*pow2(mmst1) + mmst1*pow3(mmsb1) + mmsb1*pow3(mmst1)) + 14*
     pow4(mmsb1) - 6*lb1u*pow4(mmsb1) - 6*ltu*pow4(mmsb1) + 2*lb1u*ltu*pow4(
     mmsb1) + 2*zt2*pow4(mmsb1) + pow2(ltu)*pow4(mmsb1) + pow2(lb1u)*(-(mmsb12*
     mmst1*mmt) + mmsb12*pow2(mmst1) - (2*mmst1 + mmt)*pow3(mmsb1) + pow4(mmsb1
     ))))/3.) + DeltaInv(mmt,mmsb1,mmst2)*((-2*mmgl*s2t*pow2(invdgb1)*(-14*
     mmsb12*mmst2*mmt + 6*lb1u*mmsb12*mmst2*mmt - 6*lt2u*mmsb12*mmst2*mmt + 2*
     lb1u*lt2u*mmsb12*mmst2*mmt + 12*ltu*mmsb12*mmst2*mmt - 4*lb1u*ltu*mmsb12*
     mmst2*mmt - 2*mmsb12*mmst2*mmt*zt2 - 2*mmsb12*mmst2*mmt*pow2(ltu) + 6*lt2u
     *mmsb1*mmt*pow2(mmst2) - 2*lb1u*lt2u*mmsb1*mmt*pow2(mmst2) - 6*ltu*mmsb1*
     mmt*pow2(mmst2) + 2*lb1u*ltu*mmsb1*mmt*pow2(mmst2) + mmsb1*mmt*pow2(ltu)*
     pow2(mmst2) - 14*mmsb12*pow2(mmt) + 6*lb1u*mmsb12*pow2(mmt) + 6*ltu*mmsb12
     *pow2(mmt) - 2*lb1u*ltu*mmsb12*pow2(mmt) - 28*mmsb1*mmst2*pow2(mmt) + 6*
     lt2u*mmsb1*mmst2*pow2(mmt) + 2*lb1u*lt2u*mmsb1*mmst2*pow2(mmt) + 18*ltu*
     mmsb1*mmst2*pow2(mmt) - 2*lb1u*ltu*mmsb1*mmst2*pow2(mmt) - 4*lt2u*ltu*
     mmsb1*mmst2*pow2(mmt) - 2*mmsb12*zt2*pow2(mmt) - 4*mmsb1*mmst2*zt2*pow2(
     mmt) - mmsb12*pow2(ltu)*pow2(mmt) - 3*mmsb1*mmst2*pow2(ltu)*pow2(mmt) +
     pow2(lt2u)*(mmsb12*mmst2*mmt - mmsb1*mmt*pow2(mmst2) - mmsb1*mmst2*pow2(
     mmt)) + 14*mmt*pow3(mmsb1) - 6*lb1u*mmt*pow3(mmsb1) - 6*ltu*mmt*pow3(mmsb1
     ) + 2*lb1u*ltu*mmt*pow3(mmsb1) + 2*mmt*zt2*pow3(mmsb1) + mmt*pow2(ltu)*
     pow3(mmsb1) - pow2(lb1u)*(mmsb12*mmst2*mmt + mmsb12*pow2(mmt) - mmt*pow3(
     mmsb1))))/(3.*mgl*mt) + (-(invdgb1*(14*mmsb12*mmst2 - 12*lb1u*mmsb12*mmst2
      + 6*lt2u*mmsb12*mmst2 - 6*ltu*mmsb12*mmst2 + 4*lb1u*ltu*mmsb12*mmst2 - 2*
     lt2u*ltu*mmsb12*mmst2 + 14*mmsb12*mmt - 6*lb1u*mmsb12*mmt - 6*ltu*mmsb12*
     mmt + 2*lb1u*ltu*mmsb12*mmt + 56*mmsb1*mmst2*mmt - 6*lb1u*mmsb1*mmst2*mmt
     - 6*lt2u*mmsb1*mmst2*mmt - 4*lb1u*lt2u*mmsb1*mmst2*mmt - 36*ltu*mmsb1*
     mmst2*mmt + 6*lb1u*ltu*mmsb1*mmst2*mmt + 6*lt2u*ltu*mmsb1*mmst2*mmt + 2*
     mmsb12*mmst2*zt2 + 2*mmsb12*mmt*zt2 + 8*mmsb1*mmst2*mmt*zt2 + mmsb12*mmst2
     *pow2(ltu) + mmsb12*mmt*pow2(ltu) + 6*mmsb1*mmst2*mmt*pow2(ltu) + 14*mmsb1
     *pow2(mmst2) + 6*lb1u*mmsb1*pow2(mmst2) - 12*lt2u*mmsb1*pow2(mmst2) - 6*
     ltu*mmsb1*pow2(mmst2) - 2*lb1u*ltu*mmsb1*pow2(mmst2) + 4*lt2u*ltu*mmsb1*
     pow2(mmst2) + 14*mmt*pow2(mmst2) - 6*lt2u*mmt*pow2(mmst2) - 6*ltu*mmt*pow2
     (mmst2) + 2*lt2u*ltu*mmt*pow2(mmst2) + 2*mmsb1*zt2*pow2(mmst2) + 2*mmt*zt2
     *pow2(mmst2) + mmsb1*pow2(ltu)*pow2(mmst2) + mmt*pow2(ltu)*pow2(mmst2) +
     pow2(lb1u)*(2*mmsb12*mmst2 + mmsb12*mmt + mmsb1*mmst2*mmt - mmsb1*pow2(
     mmst2) - pow3(mmsb1)) - 14*pow3(mmsb1) + 6*lb1u*pow3(mmsb1) + 6*ltu*pow3(
     mmsb1) - 2*lb1u*ltu*pow3(mmsb1) - 2*zt2*pow3(mmsb1) - pow2(ltu)*pow3(mmsb1
     ) + pow2(lt2u)*(-(mmsb12*mmst2) + mmsb1*mmst2*mmt + (2*mmsb1 + mmt)*pow2(
     mmst2) - pow3(mmst2)) - 14*pow3(mmst2) + 6*lt2u*pow3(mmst2) + 6*ltu*pow3(
     mmst2) - 2*lt2u*ltu*pow3(mmst2) - 2*zt2*pow3(mmst2) - pow2(ltu)*pow3(mmst2
     ))) - 2*pow2(invdgb1)*(-56*mmsb12*mmst2*mmt + 6*lb1u*mmsb12*mmst2*mmt + 6*
     lt2u*mmsb12*mmst2*mmt + 4*lb1u*lt2u*mmsb12*mmst2*mmt + 36*ltu*mmsb12*mmst2
     *mmt - 6*lb1u*ltu*mmsb12*mmst2*mmt - 6*lt2u*ltu*mmsb12*mmst2*mmt - 8*
     mmsb12*mmst2*mmt*zt2 - 6*mmsb12*mmst2*mmt*pow2(ltu) - 14*mmsb12*pow2(mmst2
     ) - 6*lb1u*mmsb12*pow2(mmst2) + 12*lt2u*mmsb12*pow2(mmst2) + 6*ltu*mmsb12*
     pow2(mmst2) + 2*lb1u*ltu*mmsb12*pow2(mmst2) - 4*lt2u*ltu*mmsb12*pow2(mmst2
     ) - 14*mmsb1*mmt*pow2(mmst2) + 6*lt2u*mmsb1*mmt*pow2(mmst2) + 6*ltu*mmsb1*
     mmt*pow2(mmst2) - 2*lt2u*ltu*mmsb1*mmt*pow2(mmst2) - 2*mmsb12*zt2*pow2(
     mmst2) - 2*mmsb1*mmt*zt2*pow2(mmst2) - mmsb12*pow2(ltu)*pow2(mmst2) -
     mmsb1*mmt*pow2(ltu)*pow2(mmst2) - 14*mmst2*pow3(mmsb1) + 12*lb1u*mmst2*
     pow3(mmsb1) - 6*lt2u*mmst2*pow3(mmsb1) + 6*ltu*mmst2*pow3(mmsb1) - 4*lb1u*
     ltu*mmst2*pow3(mmsb1) + 2*lt2u*ltu*mmst2*pow3(mmsb1) - 14*mmt*pow3(mmsb1)
     + 6*lb1u*mmt*pow3(mmsb1) + 6*ltu*mmt*pow3(mmsb1) - 2*lb1u*ltu*mmt*pow3(
     mmsb1) - 2*mmst2*zt2*pow3(mmsb1) - 2*mmt*zt2*pow3(mmsb1) - mmst2*pow2(ltu)
     *pow3(mmsb1) - mmt*pow2(ltu)*pow3(mmsb1) + 14*mmsb1*pow3(mmst2) - 6*lt2u*
     mmsb1*pow3(mmst2) - 6*ltu*mmsb1*pow3(mmst2) + 2*lt2u*ltu*mmsb1*pow3(mmst2)
     + 2*mmsb1*zt2*pow3(mmst2) + mmsb1*pow2(ltu)*pow3(mmst2) + pow2(lt2u)*(-(
     mmsb12*mmst2*mmt) - (2*mmsb12 + mmsb1*mmt)*pow2(mmst2) + mmst2*pow3(mmsb1)
     + mmsb1*pow3(mmst2)) + 14*pow4(mmsb1) - 6*lb1u*pow4(mmsb1) - 6*ltu*pow4(
     mmsb1) + 2*lb1u*ltu*pow4(mmsb1) + 2*zt2*pow4(mmsb1) + pow2(ltu)*pow4(mmsb1
     ) + pow2(lb1u)*(-(mmsb12*mmst2*mmt) + mmsb12*pow2(mmst2) - (2*mmst2 + mmt)
     *pow3(mmsb1) + pow4(mmsb1))))/3.) + (-3659 + 450*lb1u - 318*lb2u + 6096*
     lgu + 384*lb1u*lgu + 384*lb2u*lgu + 528*lsu + 66*lt1u + 66*lt2u - 744*ltu
     + 288*lgu*ltu - 432*invdgb1*mmsb1 + 324*invdgb2*mmsb1 + 708*invdgb1*lb1u*
     mmsb1 - 216*invdgb2*lb1u*mmsb1 - 1044*invdgb1*lgu*mmsb1 - 108*invdgb2*lgu*
     mmsb1 - 330*invdgb1*lb1u*lgu*mmsb1 + 54*invdgb2*lb1u*lgu*mmsb1 + 1080*
     invdgb1*ltu*mmsb1 - 72*invdgb1*lb1u*ltu*mmsb1 - 288*invdgb1*lgu*ltu*mmsb1
     + 5340*invdgb1*mmsb2 - 5448*invdgb2*mmsb2 + 768*invdb12*lb1u*mmsb2 - 792*
     invdgb1*lb1u*mmsb2 + 24*invdgb2*lb1u*mmsb2 - 768*invdb12*lb2u*mmsb2 - 264*
     invdgb1*lb2u*mmsb2 + 1524*invdgb2*lb2u*mmsb2 - 24*invdgb1*lb1u*lb2u*mmsb2
     + 24*invdgb2*lb1u*lb2u*mmsb2 - 3588*invdgb1*lgu*mmsb2 + 2436*invdgb2*lgu*
     mmsb2 + 24*invdgb1*lb1u*lgu*mmsb2 - 24*invdgb2*lb1u*lgu*mmsb2 + 78*invdgb1
     *lb2u*lgu*mmsb2 - 354*invdgb2*lb2u*lgu*mmsb2 + 1080*invdgb2*ltu*mmsb2 - 72
     *invdgb2*lb2u*ltu*mmsb2 - 288*invdgb2*lgu*ltu*mmsb2 + 5016*invdb12*invdgb1
     *mmsb22 - 5016*invdb12*invdgb2*mmsb22 - 768*invdb12*invdgb1*lb1u*mmsb22 -
     72*invdb12*invdgb1*lb2u*mmsb22 + 840*invdb12*invdgb2*lb2u*mmsb22 + 24*
     invdb12*invdgb1*lb1u*lb2u*mmsb22 - 24*invdb12*invdgb2*lb1u*lb2u*mmsb22 -
     3480*invdb12*invdgb1*lgu*mmsb22 + 3480*invdb12*invdgb2*lgu*mmsb22 - 24*
     invdb12*invdgb1*lb1u*lgu*mmsb22 + 24*invdb12*invdgb2*lb1u*lgu*mmsb22 - 24*
     invdb12*invdgb1*lb2u*lgu*mmsb22 + 24*invdb12*invdgb2*lb2u*lgu*mmsb22 - 432
     *invdgb1*mmst1 - 432*invdgb2*mmst1 - 108*invdgb1*lgu*mmst1 - 108*invdgb2*
     lgu*mmst1 + 108*invdgb1*lt1u*mmst1 + 108*invdgb2*lt1u*mmst1 + 54*invdgb1*
     lgu*lt1u*mmst1 + 54*invdgb2*lgu*lt1u*mmst1 + 324*invdgb1*ltu*mmst1 + 324*
     invdgb2*ltu*mmst1 - 108*invdgb1*lt1u*ltu*mmst1 - 108*invdgb2*lt1u*ltu*
     mmst1 - 432*invdgb1*mmst2 - 432*invdgb2*mmst2 - 108*invdgb1*lgu*mmst2 -
     108*invdgb2*lgu*mmst2 + 108*invdgb1*lt2u*mmst2 + 108*invdgb2*lt2u*mmst2 +
     54*invdgb1*lgu*lt2u*mmst2 + 54*invdgb2*lgu*lt2u*mmst2 + 324*invdgb1*ltu*
     mmst2 + 324*invdgb2*ltu*mmst2 - 108*invdgb1*lt2u*ltu*mmst2 - 108*invdgb2*
     lt2u*ltu*mmst2 + 2592*invdgb1*mmsusy + 2592*invdgb2*mmsusy - 864*invdgb1*
     lgu*mmsusy - 864*invdgb2*lgu*mmsusy - 1728*invdgb1*lsu*mmsusy - 1728*
     invdgb2*lsu*mmsusy + 432*invdgb1*lgu*lsu*mmsusy + 432*invdgb2*lgu*lsu*
     mmsusy - 648*invdgb1*mmt - 648*invdgb2*mmt + 216*invdgb1*lgu*mmt + 216*
     invdgb2*lgu*mmt + 432*invdgb1*ltu*mmt + 432*invdgb2*ltu*mmt - 108*invdgb1*
     lgu*ltu*mmt - 108*invdgb2*lgu*ltu*mmt - 432*zt2 - 186*invdgb1*mmsb1*zt2 +
     54*invdgb2*mmsb1*zt2 + 774*invdgb1*mmsb2*zt2 - 906*invdgb2*mmsb2*zt2 + 720
     *invdb12*invdgb1*mmsb22*zt2 - 720*invdb12*invdgb2*mmsb22*zt2 - 54*invdgb1*
     mmst1*zt2 - 54*invdgb2*mmst1*zt2 - 54*invdgb1*mmst2*zt2 - 54*invdgb2*mmst2
     *zt2 + 432*invdgb1*mmsusy*zt2 + 432*invdgb2*mmsusy*zt2 - 108*invdgb1*mmt*
     zt2 - 108*invdgb2*mmt*zt2 - 366*pow2(lb1u) + 243*invdgb1*mmsb1*pow2(lb1u)
     + 27*invdgb2*mmsb1*pow2(lb1u) - 384*invdb12*mmsb2*pow2(lb1u) + 384*invdgb1
     *mmsb2*pow2(lb1u) + 384*invdb12*invdgb1*mmsb22*pow2(lb1u) + 18*pow2(lb2u)
     + 384*invdb12*mmsb2*pow2(lb2u) + 27*invdgb1*mmsb2*pow2(lb2u) - 141*invdgb2
     *mmsb2*pow2(lb2u) - 384*invdb12*invdgb2*mmsb22*pow2(lb2u) - 1524*pow2(lgu)
     + 255*invdgb1*mmsb1*pow2(lgu) + 27*invdgb2*mmsb1*pow2(lgu) + 1059*invdgb1*
     mmsb2*pow2(lgu) - 777*invdgb2*mmsb2*pow2(lgu) + 1080*invdb12*invdgb1*
     mmsb22*pow2(lgu) - 1080*invdb12*invdgb2*mmsb22*pow2(lgu) + 27*invdgb1*
     mmst1*pow2(lgu) + 27*invdgb2*mmst1*pow2(lgu) + 27*invdgb1*mmst2*pow2(lgu)
     + 27*invdgb2*mmst2*pow2(lgu) + 216*invdgb1*mmsusy*pow2(lgu) + 216*invdgb2*
     mmsusy*pow2(lgu) - 54*invdgb1*mmt*pow2(lgu) - 54*invdgb2*mmt*pow2(lgu) +
     144*pow2(lsu) + 216*invdgb1*mmsusy*pow2(lsu) + 216*invdgb2*mmsusy*pow2(lsu
     ) + 18*pow2(lt1u) - 27*invdgb1*mmst1*pow2(lt1u) - 27*invdgb2*mmst1*pow2(
     lt1u) + 18*pow2(lt2u) - 27*invdgb1*mmst2*pow2(lt2u) - 27*invdgb2*mmst2*
     pow2(lt2u) + 216*pow2(ltu) - 180*invdgb1*mmsb1*pow2(ltu) - 180*invdgb2*
     mmsb2*pow2(ltu) - 54*invdgb1*mmst1*pow2(ltu) - 54*invdgb2*mmst1*pow2(ltu)
     - 54*invdgb1*mmst2*pow2(ltu) - 54*invdgb2*mmst2*pow2(ltu) - 54*invdgb1*mmt
     *pow2(ltu) - 54*invdgb2*mmt*pow2(ltu) + 3840*pow2(s2b) - 1152*lb1u*pow2(
     s2b) - 384*lb2u*pow2(s2b) - 2304*lgu*pow2(s2b) + 384*lb1u*lgu*pow2(s2b) +
     384*lb2u*lgu*pow2(s2b) - 4926*invdgb1*mmsb1*pow2(s2b) - 144*invdgb2*mmsb1*
     pow2(s2b) + 2268*invdgb1*lb1u*mmsb1*pow2(s2b) + 144*invdgb2*lb1u*mmsb1*
     pow2(s2b) - 96*invdgb2*lb2u*mmsb1*pow2(s2b) + 96*invdgb2*lb1u*lb2u*mmsb1*
     pow2(s2b) + 2208*invdgb1*lgu*mmsb1*pow2(s2b) + 96*invdgb2*lgu*mmsb1*pow2(
     s2b) - 672*invdgb1*lb1u*lgu*mmsb1*pow2(s2b) - 96*invdgb2*lb1u*lgu*mmsb1*
     pow2(s2b) - 7728*invdgb1*mmsb2*pow2(s2b) + 2658*invdgb2*mmsb2*pow2(s2b) -
     768*invdb12*lb1u*mmsb2*pow2(s2b) + 1116*invdgb1*lb1u*mmsb2*pow2(s2b) - 444
     *invdgb2*lb1u*mmsb2*pow2(s2b) + 768*invdb12*lb2u*mmsb2*pow2(s2b) + 564*
     invdgb1*lb2u*mmsb2*pow2(s2b) + 1080*invdgb2*lb2u*mmsb2*pow2(s2b) + 18*
     invdgb1*lb1u*lb2u*mmsb2*pow2(s2b) + 78*invdgb2*lb1u*lb2u*mmsb2*pow2(s2b) +
     4944*invdgb1*lgu*mmsb2*pow2(s2b) - 2640*invdgb2*lgu*mmsb2*pow2(s2b) - 78*
     invdgb1*lb1u*lgu*mmsb2*pow2(s2b) + 78*invdgb2*lb1u*lgu*mmsb2*pow2(s2b) -
     174*invdgb1*lb2u*lgu*mmsb2*pow2(s2b) - 594*invdgb2*lb2u*lgu*mmsb2*pow2(s2b
     ) - 7584*invdb12*invdgb1*mmsb22*pow2(s2b) + 7584*invdb12*invdgb2*mmsb22*
     pow2(s2b) + 1236*invdb12*invdgb1*lb1u*mmsb22*pow2(s2b) - 468*invdb12*
     invdgb2*lb1u*mmsb22*pow2(s2b) + 396*invdb12*invdgb1*lb2u*mmsb22*pow2(s2b)
     - 1164*invdb12*invdgb2*lb2u*mmsb22*pow2(s2b) - 30*invdb12*invdgb1*lb1u*
     lb2u*mmsb22*pow2(s2b) + 30*invdb12*invdgb2*lb1u*lb2u*mmsb22*pow2(s2b) +
     4848*invdb12*invdgb1*lgu*mmsb22*pow2(s2b) - 4848*invdb12*invdgb2*lgu*
     mmsb22*pow2(s2b) - 126*invdb12*invdgb1*lb1u*lgu*mmsb22*pow2(s2b) + 126*
     invdb12*invdgb2*lb1u*lgu*mmsb22*pow2(s2b) - 126*invdb12*invdgb1*lb2u*lgu*
     mmsb22*pow2(s2b) + 126*invdb12*invdgb2*lb2u*lgu*mmsb22*pow2(s2b) + 768*zt2
     *pow2(s2b) - 834*invdgb1*mmsb1*zt2*pow2(s2b) - 1080*invdgb1*mmsb2*zt2*pow2
     (s2b) + 246*invdgb2*mmsb2*zt2*pow2(s2b) - 1080*invdb12*invdgb1*mmsb22*zt2*
     pow2(s2b) + 1080*invdb12*invdgb2*mmsb22*zt2*pow2(s2b) + 384*pow2(lb1u)*
     pow2(s2b) - 612*invdgb1*mmsb1*pow2(lb1u)*pow2(s2b) + 384*invdb12*mmsb2*
     pow2(lb1u)*pow2(s2b) - 462*invdgb1*mmsb2*pow2(lb1u)*pow2(s2b) + 78*invdgb2
     *mmsb2*pow2(lb1u)*pow2(s2b) - 462*invdb12*invdgb1*mmsb22*pow2(lb1u)*pow2(
     s2b) + 78*invdb12*invdgb2*mmsb22*pow2(lb1u)*pow2(s2b) - 384*invdb12*mmsb2*
     pow2(lb2u)*pow2(s2b) - 78*invdgb1*mmsb2*pow2(lb2u)*pow2(s2b) - 150*invdgb2
     *mmsb2*pow2(lb2u)*pow2(s2b) - 78*invdb12*invdgb1*mmsb22*pow2(lb2u)*pow2(
     s2b) + 462*invdb12*invdgb2*mmsb22*pow2(lb2u)*pow2(s2b) + 384*pow2(lgu)*
     pow2(s2b) - 384*invdgb1*mmsb1*pow2(lgu)*pow2(s2b) - 1386*invdgb1*mmsb2*
     pow2(lgu)*pow2(s2b) + 1002*invdgb2*mmsb2*pow2(lgu)*pow2(s2b) - 1338*
     invdb12*invdgb1*mmsb22*pow2(lgu)*pow2(s2b) + 1338*invdb12*invdgb2*mmsb22*
     pow2(lgu)*pow2(s2b) + 144*lb1u*mmsb12*mmsb2*pow3(invdgb1) - 72*lb1u*lb2u*
     mmsb12*mmsb2*pow3(invdgb1) - 144*lgu*mmsb12*mmsb2*pow3(invdgb1) + 72*lb2u*
     lgu*mmsb12*mmsb2*pow3(invdgb1) + 144*lb1u*mmsb12*mmst1*pow3(invdgb1) - 144
     *lgu*mmsb12*mmst1*pow3(invdgb1) - 72*lb1u*lt1u*mmsb12*mmst1*pow3(invdgb1)
     + 72*lgu*lt1u*mmsb12*mmst1*pow3(invdgb1) + 144*lb1u*mmsb12*mmst2*pow3(
     invdgb1) - 144*lgu*mmsb12*mmst2*pow3(invdgb1) - 72*lb1u*lt2u*mmsb12*mmst2*
     pow3(invdgb1) + 72*lgu*lt2u*mmsb12*mmst2*pow3(invdgb1) + 1152*lb1u*mmsb12*
     mmsusy*pow3(invdgb1) - 1152*lgu*mmsb12*mmsusy*pow3(invdgb1) - 576*lb1u*lsu
     *mmsb12*mmsusy*pow3(invdgb1) + 576*lgu*lsu*mmsb12*mmsusy*pow3(invdgb1) -
     288*lb1u*mmsb12*mmt*pow3(invdgb1) + 288*lgu*mmsb12*mmt*pow3(invdgb1) + 144
     *lb1u*ltu*mmsb12*mmt*pow3(invdgb1) - 144*lgu*ltu*mmsb12*mmt*pow3(invdgb1)
     - 36*mmsb12*mmsb2*pow2(lb1u)*pow3(invdgb1) - 36*mmsb12*mmst1*pow2(lb1u)*
     pow3(invdgb1) - 36*mmsb12*mmst2*pow2(lb1u)*pow3(invdgb1) - 288*mmsb12*
     mmsusy*pow2(lb1u)*pow3(invdgb1) + 72*mmsb12*mmt*pow2(lb1u)*pow3(invdgb1) +
     36*mmsb12*mmsb2*pow2(lgu)*pow3(invdgb1) + 36*mmsb12*mmst1*pow2(lgu)*pow3(
     invdgb1) + 36*mmsb12*mmst2*pow2(lgu)*pow3(invdgb1) + 288*mmsb12*mmsusy*
     pow2(lgu)*pow3(invdgb1) - 72*mmsb12*mmt*pow2(lgu)*pow3(invdgb1) - 96*lb1u*
     mmsb12*mmsb2*pow2(s2b)*pow3(invdgb1) + 96*lb1u*lb2u*mmsb12*mmsb2*pow2(s2b)
     *pow3(invdgb1) + 96*lgu*mmsb12*mmsb2*pow2(s2b)*pow3(invdgb1) - 96*lb2u*lgu
     *mmsb12*mmsb2*pow2(s2b)*pow3(invdgb1) + 144*lb2u*mmsb1*mmsb22*pow3(invdgb2
     ) - 72*lb1u*lb2u*mmsb1*mmsb22*pow3(invdgb2) - 144*lgu*mmsb1*mmsb22*pow3(
     invdgb2) + 72*lb1u*lgu*mmsb1*mmsb22*pow3(invdgb2) + 144*lb2u*mmsb22*mmst1*
     pow3(invdgb2) - 144*lgu*mmsb22*mmst1*pow3(invdgb2) - 72*lb2u*lt1u*mmsb22*
     mmst1*pow3(invdgb2) + 72*lgu*lt1u*mmsb22*mmst1*pow3(invdgb2) + 144*lb2u*
     mmsb22*mmst2*pow3(invdgb2) - 144*lgu*mmsb22*mmst2*pow3(invdgb2) - 72*lb2u*
     lt2u*mmsb22*mmst2*pow3(invdgb2) + 72*lgu*lt2u*mmsb22*mmst2*pow3(invdgb2) +
     1152*lb2u*mmsb22*mmsusy*pow3(invdgb2) - 1152*lgu*mmsb22*mmsusy*pow3(
     invdgb2) - 576*lb2u*lsu*mmsb22*mmsusy*pow3(invdgb2) + 576*lgu*lsu*mmsb22*
     mmsusy*pow3(invdgb2) - 288*lb2u*mmsb22*mmt*pow3(invdgb2) + 288*lgu*mmsb22*
     mmt*pow3(invdgb2) + 144*lb2u*ltu*mmsb22*mmt*pow3(invdgb2) - 144*lgu*ltu*
     mmsb22*mmt*pow3(invdgb2) - 36*mmsb1*mmsb22*pow2(lb2u)*pow3(invdgb2) - 36*
     mmsb22*mmst1*pow2(lb2u)*pow3(invdgb2) - 36*mmsb22*mmst2*pow2(lb2u)*pow3(
     invdgb2) - 288*mmsb22*mmsusy*pow2(lb2u)*pow3(invdgb2) + 72*mmsb22*mmt*pow2
     (lb2u)*pow3(invdgb2) + 36*mmsb1*mmsb22*pow2(lgu)*pow3(invdgb2) + 36*mmsb22
     *mmst1*pow2(lgu)*pow3(invdgb2) + 36*mmsb22*mmst2*pow2(lgu)*pow3(invdgb2) +
     288*mmsb22*mmsusy*pow2(lgu)*pow3(invdgb2) - 72*mmsb22*mmt*pow2(lgu)*pow3(
     invdgb2) - 96*lb2u*mmsb1*mmsb22*pow2(s2b)*pow3(invdgb2) + 96*lb1u*lb2u*
     mmsb1*mmsb22*pow2(s2b)*pow3(invdgb2) + 96*lgu*mmsb1*mmsb22*pow2(s2b)*pow3(
     invdgb2) - 96*lb1u*lgu*mmsb1*mmsb22*pow2(s2b)*pow3(invdgb2) - 2016*pow3(
     invdgb1)*pow3(mmsb1) - 1320*lb1u*pow3(invdgb1)*pow3(mmsb1) + 3048*lgu*pow3
     (invdgb1)*pow3(mmsb1) + 696*lb1u*lgu*pow3(invdgb1)*pow3(mmsb1) - 288*zt2*
     pow3(invdgb1)*pow3(mmsb1) + 180*pow2(lb1u)*pow3(invdgb1)*pow3(mmsb1) -
     1452*pow2(lgu)*pow3(invdgb1)*pow3(mmsb1) + 120*lb1u*pow2(s2b)*pow3(invdgb1
     )*pow3(mmsb1) - 120*lgu*pow2(s2b)*pow3(invdgb1)*pow3(mmsb1) - 48*pow2(lb1u
     )*pow2(s2b)*pow3(invdgb1)*pow3(mmsb1) + 48*pow2(lgu)*pow2(s2b)*pow3(
     invdgb1)*pow3(mmsb1) + 24*invdgb1*lb1u*pow2(invdb12)*pow3(mmsb2) - 24*
     invdgb2*lb1u*pow2(invdb12)*pow3(mmsb2) - 24*invdgb1*lb2u*pow2(invdb12)*
     pow3(mmsb2) + 24*invdgb2*lb2u*pow2(invdb12)*pow3(mmsb2) + 96*invdgb1*lb1u*
     lb2u*pow2(invdb12)*pow3(mmsb2) - 96*invdgb2*lb1u*lb2u*pow2(invdb12)*pow3(
     mmsb2) - 96*invdgb1*lb1u*lgu*pow2(invdb12)*pow3(mmsb2) + 96*invdgb2*lb1u*
     lgu*pow2(invdb12)*pow3(mmsb2) - 96*invdgb1*lb2u*lgu*pow2(invdb12)*pow3(
     mmsb2) + 96*invdgb2*lb2u*lgu*pow2(invdb12)*pow3(mmsb2) + 96*invdgb1*pow2(
     invdb12)*pow2(lgu)*pow3(mmsb2) - 96*invdgb2*pow2(invdb12)*pow2(lgu)*pow3(
     mmsb2) + 24*invdgb1*lb1u*pow2(invdb12)*pow2(s2b)*pow3(mmsb2) - 24*invdgb2*
     lb1u*pow2(invdb12)*pow2(s2b)*pow3(mmsb2) - 24*invdgb1*lb2u*pow2(invdb12)*
     pow2(s2b)*pow3(mmsb2) + 24*invdgb2*lb2u*pow2(invdb12)*pow2(s2b)*pow3(mmsb2
     ) + 96*invdgb1*lb1u*lb2u*pow2(invdb12)*pow2(s2b)*pow3(mmsb2) - 96*invdgb2*
     lb1u*lb2u*pow2(invdb12)*pow2(s2b)*pow3(mmsb2) - 96*invdgb1*lb1u*lgu*pow2(
     invdb12)*pow2(s2b)*pow3(mmsb2) + 96*invdgb2*lb1u*lgu*pow2(invdb12)*pow2(
     s2b)*pow3(mmsb2) - 96*invdgb1*lb2u*lgu*pow2(invdb12)*pow2(s2b)*pow3(mmsb2)
     + 96*invdgb2*lb2u*lgu*pow2(invdb12)*pow2(s2b)*pow3(mmsb2) + 96*invdgb1*
     pow2(invdb12)*pow2(lgu)*pow2(s2b)*pow3(mmsb2) - 96*invdgb2*pow2(invdb12)*
     pow2(lgu)*pow2(s2b)*pow3(mmsb2) - 2016*pow3(invdgb2)*pow3(mmsb2) - 1320*
     lb2u*pow3(invdgb2)*pow3(mmsb2) + 3048*lgu*pow3(invdgb2)*pow3(mmsb2) + 696*
     lb2u*lgu*pow3(invdgb2)*pow3(mmsb2) - 288*zt2*pow3(invdgb2)*pow3(mmsb2) +
     180*pow2(lb2u)*pow3(invdgb2)*pow3(mmsb2) - 1452*pow2(lgu)*pow3(invdgb2)*
     pow3(mmsb2) + 120*lb2u*pow2(s2b)*pow3(invdgb2)*pow3(mmsb2) - 120*lgu*pow2(
     s2b)*pow3(invdgb2)*pow3(mmsb2) - 48*pow2(lb2u)*pow2(s2b)*pow3(invdgb2)*
     pow3(mmsb2) + 48*pow2(lgu)*pow2(s2b)*pow3(invdgb2)*pow3(mmsb2) - (96*
     invdgb1*mmgl*mmsb1*pow3(s2b))/(mb*mgl) - (96*invdgb2*mmgl*mmsb1*pow3(s2b))
     /(mb*mgl) + (192*invdgb1*lb1u*mmgl*mmsb1*pow3(s2b))/(mb*mgl) + (96*invdgb2
     *lb1u*mmgl*mmsb1*pow3(s2b))/(mb*mgl) - (96*invdgb2*lb2u*mmgl*mmsb1*pow3(
     s2b))/(mb*mgl) + (96*invdgb2*lb1u*lb2u*mmgl*mmsb1*pow3(s2b))/(mb*mgl) - (
     96*invdgb1*lgu*mmgl*mmsb1*pow3(s2b))/(mb*mgl) + (96*invdgb2*lgu*mmgl*mmsb1
     *pow3(s2b))/(mb*mgl) + (96*invdgb1*lb1u*lgu*mmgl*mmsb1*pow3(s2b))/(mb*mgl)
     - (96*invdgb2*lb1u*lgu*mmgl*mmsb1*pow3(s2b))/(mb*mgl) + (96*invdgb1*mmgl*
     mmsb2*pow3(s2b))/(mb*mgl) + (96*invdgb2*mmgl*mmsb2*pow3(s2b))/(mb*mgl) + (
     96*invdgb1*lb1u*mmgl*mmsb2*pow3(s2b))/(mb*mgl) - (96*invdgb1*lb2u*mmgl*
     mmsb2*pow3(s2b))/(mb*mgl) - (192*invdgb2*lb2u*mmgl*mmsb2*pow3(s2b))/(mb*
     mgl) + (96*invdgb1*lb1u*lb2u*mmgl*mmsb2*pow3(s2b))/(mb*mgl) + (192*invdgb2
     *lb1u*lb2u*mmgl*mmsb2*pow3(s2b))/(mb*mgl) - (96*invdgb1*lgu*mmgl*mmsb2*
     pow3(s2b))/(mb*mgl) + (96*invdgb2*lgu*mmgl*mmsb2*pow3(s2b))/(mb*mgl) + (
     192*invdgb1*lb1u*lgu*mmgl*mmsb2*pow3(s2b))/(mb*mgl) - (192*invdgb2*lb1u*
     lgu*mmgl*mmsb2*pow3(s2b))/(mb*mgl) - (96*invdgb1*lb2u*lgu*mmgl*mmsb2*pow3(
     s2b))/(mb*mgl) + (96*invdgb2*lb2u*lgu*mmgl*mmsb2*pow3(s2b))/(mb*mgl) + (
     192*invdb12*invdgb1*lb1u*lb2u*mmgl*mmsb22*pow3(s2b))/(mb*mgl) + (192*
     invdb12*invdgb2*lb1u*lb2u*mmgl*mmsb22*pow3(s2b))/(mb*mgl) + (192*invdb12*
     invdgb1*lb1u*lgu*mmgl*mmsb22*pow3(s2b))/(mb*mgl) - (192*invdb12*invdgb2*
     lb1u*lgu*mmgl*mmsb22*pow3(s2b))/(mb*mgl) - (192*invdb12*invdgb1*lb2u*lgu*
     mmgl*mmsb22*pow3(s2b))/(mb*mgl) + (192*invdb12*invdgb2*lb2u*lgu*mmgl*
     mmsb22*pow3(s2b))/(mb*mgl) - (96*invdgb1*mmgl*mmsb1*pow2(lb1u)*pow3(s2b))/
     (mb*mgl) - (192*invdgb1*mmgl*mmsb2*pow2(lb1u)*pow3(s2b))/(mb*mgl) - (192*
     invdb12*invdgb1*mmgl*mmsb22*pow2(lb1u)*pow3(s2b))/(mb*mgl) - (96*invdgb2*
     mmgl*mmsb2*pow2(lb2u)*pow3(s2b))/(mb*mgl) - (192*invdb12*invdgb2*mmgl*
     mmsb22*pow2(lb2u)*pow3(s2b))/(mb*mgl) - 72*lb1u*lgu*pow4(invdgb1)*pow4(
     mmsb1) + 36*pow2(lb1u)*pow4(invdgb1)*pow4(mmsb1) + 36*pow2(lgu)*pow4(
     invdgb1)*pow4(mmsb1) + 24*lb1u*lgu*pow2(s2b)*pow4(invdgb1)*pow4(mmsb1) -
     12*pow2(lb1u)*pow2(s2b)*pow4(invdgb1)*pow4(mmsb1) - 12*pow2(lgu)*pow2(s2b)
     *pow4(invdgb1)*pow4(mmsb1) + 48*invdgb1*lb1u*lb2u*pow3(invdb12)*pow4(mmsb2
     ) - 48*invdgb2*lb1u*lb2u*pow3(invdb12)*pow4(mmsb2) - 48*invdgb1*lb1u*lgu*
     pow3(invdb12)*pow4(mmsb2) + 48*invdgb2*lb1u*lgu*pow3(invdb12)*pow4(mmsb2)
     - 48*invdgb1*lb2u*lgu*pow3(invdb12)*pow4(mmsb2) + 48*invdgb2*lb2u*lgu*pow3
     (invdb12)*pow4(mmsb2) + 48*invdgb1*pow2(lgu)*pow3(invdb12)*pow4(mmsb2) -
     48*invdgb2*pow2(lgu)*pow3(invdb12)*pow4(mmsb2) + 48*invdgb1*lb1u*lb2u*pow2
     (s2b)*pow3(invdb12)*pow4(mmsb2) - 48*invdgb2*lb1u*lb2u*pow2(s2b)*pow3(
     invdb12)*pow4(mmsb2) - 48*invdgb1*lb1u*lgu*pow2(s2b)*pow3(invdb12)*pow4(
     mmsb2) + 48*invdgb2*lb1u*lgu*pow2(s2b)*pow3(invdb12)*pow4(mmsb2) - 48*
     invdgb1*lb2u*lgu*pow2(s2b)*pow3(invdb12)*pow4(mmsb2) + 48*invdgb2*lb2u*lgu
     *pow2(s2b)*pow3(invdb12)*pow4(mmsb2) + 48*invdgb1*pow2(lgu)*pow2(s2b)*pow3
     (invdb12)*pow4(mmsb2) - 48*invdgb2*pow2(lgu)*pow2(s2b)*pow3(invdb12)*pow4(
     mmsb2) - 72*lb2u*lgu*pow4(invdgb2)*pow4(mmsb2) + 36*pow2(lb2u)*pow4(
     invdgb2)*pow4(mmsb2) + 36*pow2(lgu)*pow4(invdgb2)*pow4(mmsb2) + 24*lb2u*
     lgu*pow2(s2b)*pow4(invdgb2)*pow4(mmsb2) - 12*pow2(lb2u)*pow2(s2b)*pow4(
     invdgb2)*pow4(mmsb2) - 12*pow2(lgu)*pow2(s2b)*pow4(invdgb2)*pow4(mmsb2) +
     3*pow2(invdgb1)*(2444*mmsb12 + 320*lb1u*mmsb12 - 2208*lgu*mmsb12 - 226*
     lb1u*lgu*mmsb12 - 288*ltu*mmsb12 + 48*lb1u*ltu*mmsb12 + 48*lgu*ltu*mmsb12
     - 144*mmsb1*mmsb2 - 28*lb1u*mmsb1*mmsb2 + 96*lb2u*mmsb1*mmsb2 + 34*lb1u*
     lb2u*mmsb1*mmsb2 + 76*lgu*mmsb1*mmsb2 - 16*lb1u*lgu*mmsb1*mmsb2 - 58*lb2u*
     lgu*mmsb1*mmsb2 + 8*lb1u*mmsb22 + 24*lb1u*lb2u*mmsb22 - 8*lgu*mmsb22 - 24*
     lb1u*lgu*mmsb22 - 24*lb2u*lgu*mmsb22 + 192*mmsb1*mmst1 - 36*lb1u*mmsb1*
     mmst1 + 84*lgu*mmsb1*mmst1 - 48*lt1u*mmsb1*mmst1 + 18*lb1u*lt1u*mmsb1*
     mmst1 - 42*lgu*lt1u*mmsb1*mmst1 - 144*ltu*mmsb1*mmst1 + 48*lt1u*ltu*mmsb1*
     mmst1 + 192*mmsb1*mmst2 - 36*lb1u*mmsb1*mmst2 + 84*lgu*mmsb1*mmst2 - 48*
     lt2u*mmsb1*mmst2 + 18*lb1u*lt2u*mmsb1*mmst2 - 42*lgu*lt2u*mmsb1*mmst2 -
     144*ltu*mmsb1*mmst2 + 48*lt2u*ltu*mmsb1*mmst2 - 1152*mmsb1*mmsusy - 288*
     lb1u*mmsb1*mmsusy + 672*lgu*mmsb1*mmsusy + 768*lsu*mmsb1*mmsusy + 144*lb1u
     *lsu*mmsb1*mmsusy - 336*lgu*lsu*mmsb1*mmsusy + 288*mmsb1*mmt + 72*lb1u*
     mmsb1*mmt - 168*lgu*mmsb1*mmt - 192*ltu*mmsb1*mmt - 36*lb1u*ltu*mmsb1*mmt
     + 84*lgu*ltu*mmsb1*mmt + 360*mmsb12*zt2 - 24*mmsb1*mmsb2*zt2 + 24*mmsb1*
     mmst1*zt2 + 24*mmsb1*mmst2*zt2 - 192*mmsb1*mmsusy*zt2 + 48*mmsb1*mmt*zt2 -
     12*mmsb1*mmsb2*pow2(lb2u) + 853*mmsb12*pow2(lgu) - 5*mmsb1*mmsb2*pow2(lgu)
     + 24*mmsb22*pow2(lgu) - 21*mmsb1*mmst1*pow2(lgu) - 21*mmsb1*mmst2*pow2(lgu
     ) - 168*mmsb1*mmsusy*pow2(lgu) + 42*mmsb1*mmt*pow2(lgu) - 96*mmsb1*mmsusy*
     pow2(lsu) + 12*mmsb1*mmst1*pow2(lt1u) + 12*mmsb1*mmst2*pow2(lt2u) + 48*
     mmsb12*pow2(ltu) + 24*mmsb1*mmst1*pow2(ltu) + 24*mmsb1*mmst2*pow2(ltu) +
     24*mmsb1*mmt*pow2(ltu) - 36*mmsb12*pow2(s2b) - 460*lb1u*mmsb12*pow2(s2b) +
     492*lgu*mmsb12*pow2(s2b) + 54*lb1u*lgu*mmsb12*pow2(s2b) + 32*mmsb1*mmsb2*
     pow2(s2b) + 72*lb1u*mmsb1*mmsb2*pow2(s2b) - 32*lb2u*mmsb1*mmsb2*pow2(s2b)
     - 48*lb1u*lb2u*mmsb1*mmsb2*pow2(s2b) - 72*lgu*mmsb1*mmsb2*pow2(s2b) - 16*
     lb1u*lgu*mmsb1*mmsb2*pow2(s2b) + 48*lb2u*lgu*mmsb1*mmsb2*pow2(s2b) + 8*
     lb1u*mmsb22*pow2(s2b) + 24*lb1u*lb2u*mmsb22*pow2(s2b) - 8*lgu*mmsb22*pow2(
     s2b) - 24*lb1u*lgu*mmsb22*pow2(s2b) - 24*lb2u*lgu*mmsb22*pow2(s2b) - 175*
     mmsb12*pow2(lgu)*pow2(s2b) + 16*mmsb1*mmsb2*pow2(lgu)*pow2(s2b) + 24*
     mmsb22*pow2(lgu)*pow2(s2b) + 8*invdb12*lb1u*pow3(mmsb2) + 32*invdb12*lb1u*
     lb2u*pow3(mmsb2) - 8*invdb12*lgu*pow3(mmsb2) - 32*invdb12*lb1u*lgu*pow3(
     mmsb2) - 32*invdb12*lb2u*lgu*pow3(mmsb2) + 32*invdb12*pow2(lgu)*pow3(mmsb2
     ) + 8*invdb12*lb1u*pow2(s2b)*pow3(mmsb2) + 32*invdb12*lb1u*lb2u*pow2(s2b)*
     pow3(mmsb2) - 8*invdb12*lgu*pow2(s2b)*pow3(mmsb2) - 32*invdb12*lb1u*lgu*
     pow2(s2b)*pow3(mmsb2) - 32*invdb12*lb2u*lgu*pow2(s2b)*pow3(mmsb2) + 32*
     invdb12*pow2(lgu)*pow2(s2b)*pow3(mmsb2) + (32*lb1u*mmgl*mmsb12*pow3(s2b))/
     (mb*mgl) - (32*lgu*mmgl*mmsb12*pow3(s2b))/(mb*mgl) + (32*lb1u*lgu*mmgl*
     mmsb12*pow3(s2b))/(mb*mgl) - (32*lb1u*mmgl*mmsb1*mmsb2*pow3(s2b))/(mb*mgl)
     + (32*lb1u*lb2u*mmgl*mmsb1*mmsb2*pow3(s2b))/(mb*mgl) + (32*lgu*mmgl*mmsb1*
     mmsb2*pow3(s2b))/(mb*mgl) - (32*lb2u*lgu*mmgl*mmsb1*mmsb2*pow3(s2b))/(mb*
     mgl) + pow2(lb1u)*(-51*mmsb12 + 9*mmsb1*(mmsb2 + mmst1 + mmst2 + 8*mmsusy
     - 2*mmt) + 121*mmsb12*pow2(s2b) - (32*mmgl*mmsb12*pow3(s2b))/(mb*mgl)) + 8
     *lb1u*lb2u*pow2(invdb12)*pow4(mmsb2) - 8*lb1u*lgu*pow2(invdb12)*pow4(mmsb2
     ) - 8*lb2u*lgu*pow2(invdb12)*pow4(mmsb2) + 8*pow2(invdb12)*pow2(lgu)*pow4(
     mmsb2) + 8*lb1u*lb2u*pow2(invdb12)*pow2(s2b)*pow4(mmsb2) - 8*lb1u*lgu*pow2
     (invdb12)*pow2(s2b)*pow4(mmsb2) - 8*lb2u*lgu*pow2(invdb12)*pow2(s2b)*pow4(
     mmsb2) + 8*pow2(invdb12)*pow2(lgu)*pow2(s2b)*pow4(mmsb2)) + 3*pow2(invdgb2
     )*(-144*mmsb1*mmsb2 + 96*lb1u*mmsb1*mmsb2 - 36*lb2u*mmsb1*mmsb2 + 18*lb1u*
     lb2u*mmsb1*mmsb2 + 84*lgu*mmsb1*mmsb2 - 42*lb1u*lgu*mmsb1*mmsb2 + 2444*
     mmsb22 + 312*lb2u*mmsb22 - 8*lb1u*lb2u*mmsb22 - 2200*lgu*mmsb22 + 8*lb1u*
     lgu*mmsb22 - 218*lb2u*lgu*mmsb22 - 288*ltu*mmsb22 + 48*lb2u*ltu*mmsb22 +
     48*lgu*ltu*mmsb22 + 192*mmsb2*mmst1 - 36*lb2u*mmsb2*mmst1 + 84*lgu*mmsb2*
     mmst1 - 48*lt1u*mmsb2*mmst1 + 18*lb2u*lt1u*mmsb2*mmst1 - 42*lgu*lt1u*mmsb2
     *mmst1 - 144*ltu*mmsb2*mmst1 + 48*lt1u*ltu*mmsb2*mmst1 + 192*mmsb2*mmst2 -
     36*lb2u*mmsb2*mmst2 + 84*lgu*mmsb2*mmst2 - 48*lt2u*mmsb2*mmst2 + 18*lb2u*
     lt2u*mmsb2*mmst2 - 42*lgu*lt2u*mmsb2*mmst2 - 144*ltu*mmsb2*mmst2 + 48*lt2u
     *ltu*mmsb2*mmst2 - 1152*mmsb2*mmsusy - 288*lb2u*mmsb2*mmsusy + 672*lgu*
     mmsb2*mmsusy + 768*lsu*mmsb2*mmsusy + 144*lb2u*lsu*mmsb2*mmsusy - 336*lgu*
     lsu*mmsb2*mmsusy + 288*mmsb2*mmt + 72*lb2u*mmsb2*mmt - 168*lgu*mmsb2*mmt -
     192*ltu*mmsb2*mmt - 36*lb2u*ltu*mmsb2*mmt + 84*lgu*ltu*mmsb2*mmt - 24*
     mmsb1*mmsb2*zt2 + 360*mmsb22*zt2 + 24*mmsb2*mmst1*zt2 + 24*mmsb2*mmst2*zt2
      - 192*mmsb2*mmsusy*zt2 + 48*mmsb2*mmt*zt2 - 12*mmsb1*mmsb2*pow2(lb1u) -
     21*mmsb1*mmsb2*pow2(lgu) + 845*mmsb22*pow2(lgu) - 21*mmsb2*mmst1*pow2(lgu)
     - 21*mmsb2*mmst2*pow2(lgu) - 168*mmsb2*mmsusy*pow2(lgu) + 42*mmsb2*mmt*
     pow2(lgu) - 96*mmsb2*mmsusy*pow2(lsu) + 12*mmsb2*mmst1*pow2(lt1u) + 12*
     mmsb2*mmst2*pow2(lt2u) + 48*mmsb22*pow2(ltu) + 24*mmsb2*mmst1*pow2(ltu) +
     24*mmsb2*mmst2*pow2(ltu) + 24*mmsb2*mmt*pow2(ltu) + 32*mmsb1*mmsb2*pow2(
     s2b) - 32*lb1u*mmsb1*mmsb2*pow2(s2b) + 64*lb2u*mmsb1*mmsb2*pow2(s2b) - 64*
     lb1u*lb2u*mmsb1*mmsb2*pow2(s2b) - 64*lgu*mmsb1*mmsb2*pow2(s2b) + 64*lb1u*
     lgu*mmsb1*mmsb2*pow2(s2b) - 36*mmsb22*pow2(s2b) - 468*lb2u*mmsb22*pow2(s2b
     ) - 8*lb1u*lb2u*mmsb22*pow2(s2b) + 500*lgu*mmsb22*pow2(s2b) + 8*lb1u*lgu*
     mmsb22*pow2(s2b) + 62*lb2u*lgu*mmsb22*pow2(s2b) - 183*mmsb22*pow2(lgu)*
     pow2(s2b) - 8*invdb12*lb2u*pow3(mmsb2) + 8*invdb12*lgu*pow3(mmsb2) - 8*
     invdb12*lb2u*pow2(s2b)*pow3(mmsb2) + 8*invdb12*lgu*pow2(s2b)*pow3(mmsb2) +
     (32*lb2u*mmgl*mmsb1*mmsb2*pow3(s2b))/(mb*mgl) - (32*lb1u*lb2u*mmgl*mmsb1*
     mmsb2*pow3(s2b))/(mb*mgl) - (32*lgu*mmgl*mmsb1*mmsb2*pow3(s2b))/(mb*mgl) +
     (32*lb1u*lgu*mmgl*mmsb1*mmsb2*pow3(s2b))/(mb*mgl) - (32*lb2u*mmgl*mmsb22*
     pow3(s2b))/(mb*mgl) + (32*lgu*mmgl*mmsb22*pow3(s2b))/(mb*mgl) - (32*lb2u*
     lgu*mmgl*mmsb22*pow3(s2b))/(mb*mgl) + pow2(lb2u)*(9*mmsb1*mmsb2 - 51*
     mmsb22 + 9*mmsb2*mmst1 + 9*mmsb2*mmst2 + 72*mmsb2*mmsusy - 18*mmsb2*mmt +
     121*mmsb22*pow2(s2b) + (32*mmgl*mmsb22*pow3(s2b))/(mb*mgl)) + 8*lb1u*lb2u*
     pow2(invdb12)*pow4(mmsb2) - 8*lb1u*lgu*pow2(invdb12)*pow4(mmsb2) - 8*lb2u*
     lgu*pow2(invdb12)*pow4(mmsb2) + 8*pow2(invdb12)*pow2(lgu)*pow4(mmsb2) + 8*
     lb1u*lb2u*pow2(invdb12)*pow2(s2b)*pow4(mmsb2) - 8*lb1u*lgu*pow2(invdb12)*
     pow2(s2b)*pow4(mmsb2) - 8*lb2u*lgu*pow2(invdb12)*pow2(s2b)*pow4(mmsb2) + 8
     *pow2(invdb12)*pow2(lgu)*pow2(s2b)*pow4(mmsb2)))/54. + DeltaInv(mmt,mmsb2,
     mmst1)*((2*mmgl*s2t*pow2(invdgb2)*(-14*mmsb22*mmst1*mmt + 6*lb2u*mmsb22*
     mmst1*mmt - 6*lt1u*mmsb22*mmst1*mmt + 2*lb2u*lt1u*mmsb22*mmst1*mmt + 12*
     ltu*mmsb22*mmst1*mmt - 4*lb2u*ltu*mmsb22*mmst1*mmt - 2*mmsb22*mmst1*mmt*
     zt2 - 2*mmsb22*mmst1*mmt*pow2(ltu) + 6*lt1u*mmsb2*mmt*pow2(mmst1) - 2*lb2u
     *lt1u*mmsb2*mmt*pow2(mmst1) - 6*ltu*mmsb2*mmt*pow2(mmst1) + 2*lb2u*ltu*
     mmsb2*mmt*pow2(mmst1) + mmsb2*mmt*pow2(ltu)*pow2(mmst1) - 14*mmsb22*pow2(
     mmt) + 6*lb2u*mmsb22*pow2(mmt) + 6*ltu*mmsb22*pow2(mmt) - 2*lb2u*ltu*
     mmsb22*pow2(mmt) - 28*mmsb2*mmst1*pow2(mmt) + 6*lt1u*mmsb2*mmst1*pow2(mmt)
     + 2*lb2u*lt1u*mmsb2*mmst1*pow2(mmt) + 18*ltu*mmsb2*mmst1*pow2(mmt) - 2*
     lb2u*ltu*mmsb2*mmst1*pow2(mmt) - 4*lt1u*ltu*mmsb2*mmst1*pow2(mmt) - 2*
     mmsb22*zt2*pow2(mmt) - 4*mmsb2*mmst1*zt2*pow2(mmt) - mmsb22*pow2(ltu)*pow2
     (mmt) - 3*mmsb2*mmst1*pow2(ltu)*pow2(mmt) + pow2(lt1u)*(mmsb22*mmst1*mmt -
     mmsb2*mmt*pow2(mmst1) - mmsb2*mmst1*pow2(mmt)) + 14*mmt*pow3(mmsb2) - 6*
     lb2u*mmt*pow3(mmsb2) - 6*ltu*mmt*pow3(mmsb2) + 2*lb2u*ltu*mmt*pow3(mmsb2)
     + 2*mmt*zt2*pow3(mmsb2) + mmt*pow2(ltu)*pow3(mmsb2) - pow2(lb2u)*(mmsb22*
     mmst1*mmt + mmsb22*pow2(mmt) - mmt*pow3(mmsb2))))/(3.*mgl*mt) + (-(invdgb2
     *(14*mmsb22*mmst1 - 12*lb2u*mmsb22*mmst1 + 6*lt1u*mmsb22*mmst1 - 6*ltu*
     mmsb22*mmst1 + 4*lb2u*ltu*mmsb22*mmst1 - 2*lt1u*ltu*mmsb22*mmst1 + 14*
     mmsb22*mmt - 6*lb2u*mmsb22*mmt - 6*ltu*mmsb22*mmt + 2*lb2u*ltu*mmsb22*mmt
     + 56*mmsb2*mmst1*mmt - 6*lb2u*mmsb2*mmst1*mmt - 6*lt1u*mmsb2*mmst1*mmt - 4
     *lb2u*lt1u*mmsb2*mmst1*mmt - 36*ltu*mmsb2*mmst1*mmt + 6*lb2u*ltu*mmsb2*
     mmst1*mmt + 6*lt1u*ltu*mmsb2*mmst1*mmt + 2*mmsb22*mmst1*zt2 + 2*mmsb22*mmt
     *zt2 + 8*mmsb2*mmst1*mmt*zt2 + mmsb22*mmst1*pow2(ltu) + mmsb22*mmt*pow2(
     ltu) + 6*mmsb2*mmst1*mmt*pow2(ltu) + 14*mmsb2*pow2(mmst1) + 6*lb2u*mmsb2*
     pow2(mmst1) - 12*lt1u*mmsb2*pow2(mmst1) - 6*ltu*mmsb2*pow2(mmst1) - 2*lb2u
     *ltu*mmsb2*pow2(mmst1) + 4*lt1u*ltu*mmsb2*pow2(mmst1) + 14*mmt*pow2(mmst1)
     - 6*lt1u*mmt*pow2(mmst1) - 6*ltu*mmt*pow2(mmst1) + 2*lt1u*ltu*mmt*pow2(
     mmst1) + 2*mmsb2*zt2*pow2(mmst1) + 2*mmt*zt2*pow2(mmst1) + mmsb2*pow2(ltu)
     *pow2(mmst1) + mmt*pow2(ltu)*pow2(mmst1) + pow2(lb2u)*(2*mmsb22*mmst1 +
     mmsb22*mmt + mmsb2*mmst1*mmt - mmsb2*pow2(mmst1) - pow3(mmsb2)) - 14*pow3(
     mmsb2) + 6*lb2u*pow3(mmsb2) + 6*ltu*pow3(mmsb2) - 2*lb2u*ltu*pow3(mmsb2) -
     2*zt2*pow3(mmsb2) - pow2(ltu)*pow3(mmsb2) + pow2(lt1u)*(-(mmsb22*mmst1) +
     mmsb2*mmst1*mmt + (2*mmsb2 + mmt)*pow2(mmst1) - pow3(mmst1)) - 14*pow3(
     mmst1) + 6*lt1u*pow3(mmst1) + 6*ltu*pow3(mmst1) - 2*lt1u*ltu*pow3(mmst1) -
     2*zt2*pow3(mmst1) - pow2(ltu)*pow3(mmst1))) - 2*pow2(invdgb2)*(-56*mmsb22*
     mmst1*mmt + 6*lb2u*mmsb22*mmst1*mmt + 6*lt1u*mmsb22*mmst1*mmt + 4*lb2u*
     lt1u*mmsb22*mmst1*mmt + 36*ltu*mmsb22*mmst1*mmt - 6*lb2u*ltu*mmsb22*mmst1*
     mmt - 6*lt1u*ltu*mmsb22*mmst1*mmt - 8*mmsb22*mmst1*mmt*zt2 - 6*mmsb22*
     mmst1*mmt*pow2(ltu) - 14*mmsb22*pow2(mmst1) - 6*lb2u*mmsb22*pow2(mmst1) +
     12*lt1u*mmsb22*pow2(mmst1) + 6*ltu*mmsb22*pow2(mmst1) + 2*lb2u*ltu*mmsb22*
     pow2(mmst1) - 4*lt1u*ltu*mmsb22*pow2(mmst1) - 14*mmsb2*mmt*pow2(mmst1) + 6
     *lt1u*mmsb2*mmt*pow2(mmst1) + 6*ltu*mmsb2*mmt*pow2(mmst1) - 2*lt1u*ltu*
     mmsb2*mmt*pow2(mmst1) - 2*mmsb22*zt2*pow2(mmst1) - 2*mmsb2*mmt*zt2*pow2(
     mmst1) - mmsb22*pow2(ltu)*pow2(mmst1) - mmsb2*mmt*pow2(ltu)*pow2(mmst1) -
     14*mmst1*pow3(mmsb2) + 12*lb2u*mmst1*pow3(mmsb2) - 6*lt1u*mmst1*pow3(mmsb2
     ) + 6*ltu*mmst1*pow3(mmsb2) - 4*lb2u*ltu*mmst1*pow3(mmsb2) + 2*lt1u*ltu*
     mmst1*pow3(mmsb2) - 14*mmt*pow3(mmsb2) + 6*lb2u*mmt*pow3(mmsb2) + 6*ltu*
     mmt*pow3(mmsb2) - 2*lb2u*ltu*mmt*pow3(mmsb2) - 2*mmst1*zt2*pow3(mmsb2) - 2
     *mmt*zt2*pow3(mmsb2) - mmst1*pow2(ltu)*pow3(mmsb2) - mmt*pow2(ltu)*pow3(
     mmsb2) + 14*mmsb2*pow3(mmst1) - 6*lt1u*mmsb2*pow3(mmst1) - 6*ltu*mmsb2*
     pow3(mmst1) + 2*lt1u*ltu*mmsb2*pow3(mmst1) + 2*mmsb2*zt2*pow3(mmst1) +
     mmsb2*pow2(ltu)*pow3(mmst1) + pow2(lt1u)*(-(mmsb22*mmst1*mmt) - (2*mmsb22
     + mmsb2*mmt)*pow2(mmst1) + mmst1*pow3(mmsb2) + mmsb2*pow3(mmst1)) + 14*
     pow4(mmsb2) - 6*lb2u*pow4(mmsb2) - 6*ltu*pow4(mmsb2) + 2*lb2u*ltu*pow4(
     mmsb2) + 2*zt2*pow4(mmsb2) + pow2(ltu)*pow4(mmsb2) + pow2(lb2u)*(-(mmsb22*
     mmst1*mmt) + mmsb22*pow2(mmst1) - (2*mmst1 + mmt)*pow3(mmsb2) + pow4(mmsb2
     ))))/3.) + DeltaInv(mmt,mmst1,mmgl)*(s2t*((-4*s2b*(14*mmgl*mmsb1*mmt - 6*
     lgu*mmgl*mmsb1*mmt - 6*ltu*mmgl*mmsb1*mmt + 2*lgu*ltu*mmgl*mmsb1*mmt + 14*
     mmsb12*mmt - 6*lgu*mmsb12*mmt - 6*ltu*mmsb12*mmt + 2*lgu*ltu*mmsb12*mmt -
     14*mmgl*mmsb2*mmt + 6*lgu*mmgl*mmsb2*mmt + 6*ltu*mmgl*mmsb2*mmt - 2*lgu*
     ltu*mmgl*mmsb2*mmt - 14*mmsb22*mmt + 6*lgu*mmsb22*mmt + 6*ltu*mmsb22*mmt -
     2*lgu*ltu*mmsb22*mmt - 14*mmsb1*mmst1*mmt + 6*lgu*mmsb1*mmst1*mmt - 6*lt1u
     *mmsb1*mmst1*mmt + 2*lgu*lt1u*mmsb1*mmst1*mmt + 12*ltu*mmsb1*mmst1*mmt - 4
     *lgu*ltu*mmsb1*mmst1*mmt + 14*invdgb1*mmsb12*mmst1*mmt - 6*invdgb1*lgu*
     mmsb12*mmst1*mmt + 6*invdgb1*lt1u*mmsb12*mmst1*mmt - 2*invdgb1*lgu*lt1u*
     mmsb12*mmst1*mmt - 12*invdgb1*ltu*mmsb12*mmst1*mmt + 4*invdgb1*lgu*ltu*
     mmsb12*mmst1*mmt + 14*mmsb2*mmst1*mmt - 6*lgu*mmsb2*mmst1*mmt + 6*lt1u*
     mmsb2*mmst1*mmt - 2*lgu*lt1u*mmsb2*mmst1*mmt - 12*ltu*mmsb2*mmst1*mmt + 4*
     lgu*ltu*mmsb2*mmst1*mmt - 14*invdgb2*mmsb22*mmst1*mmt + 6*invdgb2*lgu*
     mmsb22*mmst1*mmt - 6*invdgb2*lt1u*mmsb22*mmst1*mmt + 2*invdgb2*lgu*lt1u*
     mmsb22*mmst1*mmt + 12*invdgb2*ltu*mmsb22*mmst1*mmt - 4*invdgb2*lgu*ltu*
     mmsb22*mmst1*mmt + 2*mmgl*mmsb1*mmt*zt2 + 2*mmsb12*mmt*zt2 - 2*mmgl*mmsb2*
     mmt*zt2 - 2*mmsb22*mmt*zt2 - 2*mmsb1*mmst1*mmt*zt2 + 2*invdgb1*mmsb12*
     mmst1*mmt*zt2 + 2*mmsb2*mmst1*mmt*zt2 - 2*invdgb2*mmsb22*mmst1*mmt*zt2 +
     mmgl*mmsb1*mmt*pow2(ltu) + mmsb12*mmt*pow2(ltu) - mmgl*mmsb2*mmt*pow2(ltu)
     - mmsb22*mmt*pow2(ltu) - 2*mmsb1*mmst1*mmt*pow2(ltu) + 2*invdgb1*mmsb12*
     mmst1*mmt*pow2(ltu) + 2*mmsb2*mmst1*mmt*pow2(ltu) - 2*invdgb2*mmsb22*mmst1
     *mmt*pow2(ltu) - 6*invdgb1*lt1u*mmsb1*mmt*pow2(mmst1) + 2*invdgb1*lgu*lt1u
     *mmsb1*mmt*pow2(mmst1) + 6*invdgb1*ltu*mmsb1*mmt*pow2(mmst1) - 2*invdgb1*
     lgu*ltu*mmsb1*mmt*pow2(mmst1) + 6*invdgb2*lt1u*mmsb2*mmt*pow2(mmst1) - 2*
     invdgb2*lgu*lt1u*mmsb2*mmt*pow2(mmst1) - 6*invdgb2*ltu*mmsb2*mmt*pow2(
     mmst1) + 2*invdgb2*lgu*ltu*mmsb2*mmt*pow2(mmst1) - invdgb1*mmsb1*mmt*pow2(
     ltu)*pow2(mmst1) + invdgb2*mmsb2*mmt*pow2(ltu)*pow2(mmst1) - 14*mmsb1*pow2
     (mmt) + 6*lgu*mmsb1*pow2(mmt) + 6*ltu*mmsb1*pow2(mmt) - 2*lgu*ltu*mmsb1*
     pow2(mmt) + 14*invdgb1*mmsb12*pow2(mmt) - 6*invdgb1*lgu*mmsb12*pow2(mmt) -
     6*invdgb1*ltu*mmsb12*pow2(mmt) + 2*invdgb1*lgu*ltu*mmsb12*pow2(mmt) + 14*
     mmsb2*pow2(mmt) - 6*lgu*mmsb2*pow2(mmt) - 6*ltu*mmsb2*pow2(mmt) + 2*lgu*
     ltu*mmsb2*pow2(mmt) - 14*invdgb2*mmsb22*pow2(mmt) + 6*invdgb2*lgu*mmsb22*
     pow2(mmt) + 6*invdgb2*ltu*mmsb22*pow2(mmt) - 2*invdgb2*lgu*ltu*mmsb22*pow2
     (mmt) + 28*invdgb1*mmsb1*mmst1*pow2(mmt) - 6*invdgb1*lt1u*mmsb1*mmst1*pow2
     (mmt) - 2*invdgb1*lgu*lt1u*mmsb1*mmst1*pow2(mmt) - 18*invdgb1*ltu*mmsb1*
     mmst1*pow2(mmt) + 2*invdgb1*lgu*ltu*mmsb1*mmst1*pow2(mmt) + 4*invdgb1*lt1u
     *ltu*mmsb1*mmst1*pow2(mmt) - 28*invdgb2*mmsb2*mmst1*pow2(mmt) + 6*invdgb2*
     lt1u*mmsb2*mmst1*pow2(mmt) + 2*invdgb2*lgu*lt1u*mmsb2*mmst1*pow2(mmt) + 18
     *invdgb2*ltu*mmsb2*mmst1*pow2(mmt) - 2*invdgb2*lgu*ltu*mmsb2*mmst1*pow2(
     mmt) - 4*invdgb2*lt1u*ltu*mmsb2*mmst1*pow2(mmt) - 2*mmsb1*zt2*pow2(mmt) +
     2*invdgb1*mmsb12*zt2*pow2(mmt) + 2*mmsb2*zt2*pow2(mmt) - 2*invdgb2*mmsb22*
     zt2*pow2(mmt) + 4*invdgb1*mmsb1*mmst1*zt2*pow2(mmt) - 4*invdgb2*mmsb2*
     mmst1*zt2*pow2(mmt) - mmsb1*pow2(ltu)*pow2(mmt) + invdgb1*mmsb12*pow2(ltu)
     *pow2(mmt) + mmsb2*pow2(ltu)*pow2(mmt) - invdgb2*mmsb22*pow2(ltu)*pow2(mmt
     ) + 3*invdgb1*mmsb1*mmst1*pow2(ltu)*pow2(mmt) - 3*invdgb2*mmsb2*mmst1*pow2
     (ltu)*pow2(mmt) + pow2(lt1u)*((mmsb1 - invdgb1*mmsb12 - mmsb2 + invdgb2*
     mmsb22)*mmst1*mmt + (invdgb1*mmsb1 - invdgb2*mmsb2)*mmt*pow2(mmst1) + (
     invdgb1*mmsb1 - invdgb2*mmsb2)*mmst1*pow2(mmt)) - 14*invdgb1*mmt*pow3(
     mmsb1) + 6*invdgb1*lgu*mmt*pow3(mmsb1) + 6*invdgb1*ltu*mmt*pow3(mmsb1) - 2
     *invdgb1*lgu*ltu*mmt*pow3(mmsb1) - 2*invdgb1*mmt*zt2*pow3(mmsb1) - invdgb1
     *mmt*pow2(ltu)*pow3(mmsb1) + 14*invdgb2*mmt*pow3(mmsb2) - 6*invdgb2*lgu*
     mmt*pow3(mmsb2) - 6*invdgb2*ltu*mmt*pow3(mmsb2) + 2*invdgb2*lgu*ltu*mmt*
     pow3(mmsb2) + 2*invdgb2*mmt*zt2*pow3(mmsb2) + invdgb2*mmt*pow2(ltu)*pow3(
     mmsb2) + pow2(lgu)*((-mmsb1 + invdgb1*mmsb12 + mmsb2 - invdgb2*mmsb22)*
     pow2(mmt) + mmt*(mmgl*mmsb1 + mmsb12 - mmgl*mmsb2 - mmsb22 - mmsb1*mmst1 +
     invdgb1*mmsb12*mmst1 + mmsb2*mmst1 - invdgb2*mmsb22*mmst1 - invdgb1*pow3(
     mmsb1) + invdgb2*pow3(mmsb2)))))/(3.*mb*mt) + (2*(28*mmgl2*mmt - 12*lgu*
     mmgl2*mmt - 12*ltu*mmgl2*mmt + 4*lgu*ltu*mmgl2*mmt + 28*mmgl*mmsb1*mmt -
     12*lgu*mmgl*mmsb1*mmt - 12*ltu*mmgl*mmsb1*mmt + 4*lgu*ltu*mmgl*mmsb1*mmt -
     42*invdgb1*mmgl*mmsb12*mmt + 18*invdgb1*lgu*mmgl*mmsb12*mmt + 18*invdgb1*
     ltu*mmgl*mmsb12*mmt - 6*invdgb1*lgu*ltu*mmgl*mmsb12*mmt + 28*mmgl*mmsb2*
     mmt - 12*lgu*mmgl*mmsb2*mmt - 12*ltu*mmgl*mmsb2*mmt + 4*lgu*ltu*mmgl*mmsb2
     *mmt - 42*invdgb2*mmgl*mmsb22*mmt + 18*invdgb2*lgu*mmgl*mmsb22*mmt + 18*
     invdgb2*ltu*mmgl*mmsb22*mmt - 6*invdgb2*lgu*ltu*mmgl*mmsb22*mmt - 28*mmgl*
     mmst1*mmt + 12*lgu*mmgl*mmst1*mmt - 12*lt1u*mmgl*mmst1*mmt + 4*lgu*lt1u*
     mmgl*mmst1*mmt + 24*ltu*mmgl*mmst1*mmt - 8*lgu*ltu*mmgl*mmst1*mmt + 28*
     invdgb1*mmgl*mmsb1*mmst1*mmt - 12*invdgb1*lgu*mmgl*mmsb1*mmst1*mmt + 12*
     invdgb1*lt1u*mmgl*mmsb1*mmst1*mmt - 4*invdgb1*lgu*lt1u*mmgl*mmsb1*mmst1*
     mmt - 24*invdgb1*ltu*mmgl*mmsb1*mmst1*mmt + 8*invdgb1*lgu*ltu*mmgl*mmsb1*
     mmst1*mmt + 28*invdgb2*mmgl*mmsb2*mmst1*mmt - 12*invdgb2*lgu*mmgl*mmsb2*
     mmst1*mmt + 12*invdgb2*lt1u*mmgl*mmsb2*mmst1*mmt - 4*invdgb2*lgu*lt1u*mmgl
     *mmsb2*mmst1*mmt - 24*invdgb2*ltu*mmgl*mmsb2*mmst1*mmt + 8*invdgb2*lgu*ltu
     *mmgl*mmsb2*mmst1*mmt + 4*mmgl2*mmt*zt2 + 4*mmgl*mmsb1*mmt*zt2 - 6*invdgb1
     *mmgl*mmsb12*mmt*zt2 + 4*mmgl*mmsb2*mmt*zt2 - 6*invdgb2*mmgl*mmsb22*mmt*
     zt2 - 4*mmgl*mmst1*mmt*zt2 + 4*invdgb1*mmgl*mmsb1*mmst1*mmt*zt2 + 4*
     invdgb2*mmgl*mmsb2*mmst1*mmt*zt2 + 2*mmgl2*mmt*pow2(lgu) + 2*mmgl*mmsb1*
     mmt*pow2(lgu) - 3*invdgb1*mmgl*mmsb12*mmt*pow2(lgu) + 2*mmgl*mmsb2*mmt*
     pow2(lgu) - 3*invdgb2*mmgl*mmsb22*mmt*pow2(lgu) - 2*mmgl*mmst1*mmt*pow2(
     lgu) + 2*invdgb1*mmgl*mmsb1*mmst1*mmt*pow2(lgu) + 2*invdgb2*mmgl*mmsb2*
     mmst1*mmt*pow2(lgu) + 2*mmgl*mmst1*mmt*pow2(lt1u) - 2*invdgb1*mmgl*mmsb1*
     mmst1*mmt*pow2(lt1u) - 2*invdgb2*mmgl*mmsb2*mmst1*mmt*pow2(lt1u) + 2*mmgl2
     *mmt*pow2(ltu) + 2*mmgl*mmsb1*mmt*pow2(ltu) - 3*invdgb1*mmgl*mmsb12*mmt*
     pow2(ltu) + 2*mmgl*mmsb2*mmt*pow2(ltu) - 3*invdgb2*mmgl*mmsb22*mmt*pow2(
     ltu) - 4*mmgl*mmst1*mmt*pow2(ltu) + 4*invdgb1*mmgl*mmsb1*mmst1*mmt*pow2(
     ltu) + 4*invdgb2*mmgl*mmsb2*mmst1*mmt*pow2(ltu) - 6*invdgb1*lt1u*mmgl*mmt*
     pow2(mmst1) - 6*invdgb2*lt1u*mmgl*mmt*pow2(mmst1) + 2*invdgb1*lgu*lt1u*
     mmgl*mmt*pow2(mmst1) + 2*invdgb2*lgu*lt1u*mmgl*mmt*pow2(mmst1) + 6*invdgb1
     *ltu*mmgl*mmt*pow2(mmst1) + 6*invdgb2*ltu*mmgl*mmt*pow2(mmst1) - 2*invdgb1
     *lgu*ltu*mmgl*mmt*pow2(mmst1) - 2*invdgb2*lgu*ltu*mmgl*mmt*pow2(mmst1) +
     invdgb1*mmgl*mmt*pow2(lt1u)*pow2(mmst1) + invdgb2*mmgl*mmt*pow2(lt1u)*pow2
     (mmst1) - invdgb1*mmgl*mmt*pow2(ltu)*pow2(mmst1) - invdgb2*mmgl*mmt*pow2(
     ltu)*pow2(mmst1) - 28*mmgl*pow2(mmt) + 12*lgu*mmgl*pow2(mmt) + 12*ltu*mmgl
     *pow2(mmt) - 4*lgu*ltu*mmgl*pow2(mmt) + 28*invdgb1*mmgl*mmsb1*pow2(mmt) -
     12*invdgb1*lgu*mmgl*mmsb1*pow2(mmt) - 12*invdgb1*ltu*mmgl*mmsb1*pow2(mmt)
     + 4*invdgb1*lgu*ltu*mmgl*mmsb1*pow2(mmt) + 28*invdgb2*mmgl*mmsb2*pow2(mmt)
     - 12*invdgb2*lgu*mmgl*mmsb2*pow2(mmt) - 12*invdgb2*ltu*mmgl*mmsb2*pow2(mmt
     ) + 4*invdgb2*lgu*ltu*mmgl*mmsb2*pow2(mmt) + 28*invdgb1*mmgl*mmst1*pow2(
     mmt) + 28*invdgb2*mmgl*mmst1*pow2(mmt) - 6*invdgb1*lt1u*mmgl*mmst1*pow2(
     mmt) - 6*invdgb2*lt1u*mmgl*mmst1*pow2(mmt) - 2*invdgb1*lgu*lt1u*mmgl*mmst1
     *pow2(mmt) - 2*invdgb2*lgu*lt1u*mmgl*mmst1*pow2(mmt) - 18*invdgb1*ltu*mmgl
     *mmst1*pow2(mmt) - 18*invdgb2*ltu*mmgl*mmst1*pow2(mmt) + 2*invdgb1*lgu*ltu
     *mmgl*mmst1*pow2(mmt) + 2*invdgb2*lgu*ltu*mmgl*mmst1*pow2(mmt) + 4*invdgb1
     *lt1u*ltu*mmgl*mmst1*pow2(mmt) + 4*invdgb2*lt1u*ltu*mmgl*mmst1*pow2(mmt) -
     4*mmgl*zt2*pow2(mmt) + 4*invdgb1*mmgl*mmsb1*zt2*pow2(mmt) + 4*invdgb2*mmgl
     *mmsb2*zt2*pow2(mmt) + 4*invdgb1*mmgl*mmst1*zt2*pow2(mmt) + 4*invdgb2*mmgl
     *mmst1*zt2*pow2(mmt) - 2*mmgl*pow2(lgu)*pow2(mmt) + 2*invdgb1*mmgl*mmsb1*
     pow2(lgu)*pow2(mmt) + 2*invdgb2*mmgl*mmsb2*pow2(lgu)*pow2(mmt) + invdgb1*
     mmgl*mmst1*pow2(lt1u)*pow2(mmt) + invdgb2*mmgl*mmst1*pow2(lt1u)*pow2(mmt)
     - 2*mmgl*pow2(ltu)*pow2(mmt) + 2*invdgb1*mmgl*mmsb1*pow2(ltu)*pow2(mmt) +
     2*invdgb2*mmgl*mmsb2*pow2(ltu)*pow2(mmt) + 3*invdgb1*mmgl*mmst1*pow2(ltu)*
     pow2(mmt) + 3*invdgb2*mmgl*mmst1*pow2(ltu)*pow2(mmt) + mmgl*pow2(invdgb1)*
     (-14*mmsb12*mmst1*mmt + 6*lgu*mmsb12*mmst1*mmt - 6*lt1u*mmsb12*mmst1*mmt +
     2*lgu*lt1u*mmsb12*mmst1*mmt + 12*ltu*mmsb12*mmst1*mmt - 4*lgu*ltu*mmsb12*
     mmst1*mmt - 2*mmsb12*mmst1*mmt*zt2 - 2*mmsb12*mmst1*mmt*pow2(ltu) + 6*lt1u
     *mmsb1*mmt*pow2(mmst1) - 2*lgu*lt1u*mmsb1*mmt*pow2(mmst1) - 6*ltu*mmsb1*
     mmt*pow2(mmst1) + 2*lgu*ltu*mmsb1*mmt*pow2(mmst1) + mmsb1*mmt*pow2(ltu)*
     pow2(mmst1) - 14*mmsb12*pow2(mmt) + 6*lgu*mmsb12*pow2(mmt) + 6*ltu*mmsb12*
     pow2(mmt) - 2*lgu*ltu*mmsb12*pow2(mmt) - 28*mmsb1*mmst1*pow2(mmt) + 6*lt1u
     *mmsb1*mmst1*pow2(mmt) + 2*lgu*lt1u*mmsb1*mmst1*pow2(mmt) + 18*ltu*mmsb1*
     mmst1*pow2(mmt) - 2*lgu*ltu*mmsb1*mmst1*pow2(mmt) - 4*lt1u*ltu*mmsb1*mmst1
     *pow2(mmt) - 2*mmsb12*zt2*pow2(mmt) - 4*mmsb1*mmst1*zt2*pow2(mmt) - mmsb12
     *pow2(ltu)*pow2(mmt) - 3*mmsb1*mmst1*pow2(ltu)*pow2(mmt) + pow2(lt1u)*(
     mmsb12*mmst1*mmt - mmsb1*mmt*pow2(mmst1) - mmsb1*mmst1*pow2(mmt)) + 14*mmt
     *pow3(mmsb1) - 6*lgu*mmt*pow3(mmsb1) - 6*ltu*mmt*pow3(mmsb1) + 2*lgu*ltu*
     mmt*pow3(mmsb1) + 2*mmt*zt2*pow3(mmsb1) + mmt*pow2(ltu)*pow3(mmsb1) - pow2
     (lgu)*(mmsb12*mmst1*mmt + mmsb12*pow2(mmt) - mmt*pow3(mmsb1))) + mmgl*pow2
     (invdgb2)*(-14*mmsb22*mmst1*mmt + 6*lgu*mmsb22*mmst1*mmt - 6*lt1u*mmsb22*
     mmst1*mmt + 2*lgu*lt1u*mmsb22*mmst1*mmt + 12*ltu*mmsb22*mmst1*mmt - 4*lgu*
     ltu*mmsb22*mmst1*mmt - 2*mmsb22*mmst1*mmt*zt2 - 2*mmsb22*mmst1*mmt*pow2(
     ltu) + 6*lt1u*mmsb2*mmt*pow2(mmst1) - 2*lgu*lt1u*mmsb2*mmt*pow2(mmst1) - 6
     *ltu*mmsb2*mmt*pow2(mmst1) + 2*lgu*ltu*mmsb2*mmt*pow2(mmst1) + mmsb2*mmt*
     pow2(ltu)*pow2(mmst1) - 14*mmsb22*pow2(mmt) + 6*lgu*mmsb22*pow2(mmt) + 6*
     ltu*mmsb22*pow2(mmt) - 2*lgu*ltu*mmsb22*pow2(mmt) - 28*mmsb2*mmst1*pow2(
     mmt) + 6*lt1u*mmsb2*mmst1*pow2(mmt) + 2*lgu*lt1u*mmsb2*mmst1*pow2(mmt) +
     18*ltu*mmsb2*mmst1*pow2(mmt) - 2*lgu*ltu*mmsb2*mmst1*pow2(mmt) - 4*lt1u*
     ltu*mmsb2*mmst1*pow2(mmt) - 2*mmsb22*zt2*pow2(mmt) - 4*mmsb2*mmst1*zt2*
     pow2(mmt) - mmsb22*pow2(ltu)*pow2(mmt) - 3*mmsb2*mmst1*pow2(ltu)*pow2(mmt)
     + pow2(lt1u)*(mmsb22*mmst1*mmt - mmsb2*mmt*pow2(mmst1) - mmsb2*mmst1*pow2(
     mmt)) + 14*mmt*pow3(mmsb2) - 6*lgu*mmt*pow3(mmsb2) - 6*ltu*mmt*pow3(mmsb2)
     + 2*lgu*ltu*mmt*pow3(mmsb2) + 2*mmt*zt2*pow3(mmsb2) + mmt*pow2(ltu)*pow3(
     mmsb2) - pow2(lgu)*(mmsb22*mmst1*mmt + mmsb22*pow2(mmt) - mmt*pow3(mmsb2))
     )))/(3.*mgl*mt)) + (4*s2b*(14*mmgl2*mmsb1 - 6*lgu*mmgl2*mmsb1 - 6*ltu*
     mmgl2*mmsb1 + 2*lgu*ltu*mmgl2*mmsb1 + 14*mmgl*mmsb12 - 6*lgu*mmgl*mmsb12 -
     6*ltu*mmgl*mmsb12 + 2*lgu*ltu*mmgl*mmsb12 - 14*mmgl2*mmsb2 + 6*lgu*mmgl2*
     mmsb2 + 6*ltu*mmgl2*mmsb2 - 2*lgu*ltu*mmgl2*mmsb2 - 14*mmgl*mmsb22 + 6*lgu
     *mmgl*mmsb22 + 6*ltu*mmgl*mmsb22 - 2*lgu*ltu*mmgl*mmsb22 - 14*mmgl*mmsb1*
     mmst1 + 12*lgu*mmgl*mmsb1*mmst1 - 6*lt1u*mmgl*mmsb1*mmst1 + 6*ltu*mmgl*
     mmsb1*mmst1 - 4*lgu*ltu*mmgl*mmsb1*mmst1 + 2*lt1u*ltu*mmgl*mmsb1*mmst1 +
     14*invdgb1*mmgl*mmsb12*mmst1 - 12*invdgb1*lgu*mmgl*mmsb12*mmst1 + 6*
     invdgb1*lt1u*mmgl*mmsb12*mmst1 - 6*invdgb1*ltu*mmgl*mmsb12*mmst1 + 4*
     invdgb1*lgu*ltu*mmgl*mmsb12*mmst1 - 2*invdgb1*lt1u*ltu*mmgl*mmsb12*mmst1 +
     14*mmgl*mmsb2*mmst1 - 12*lgu*mmgl*mmsb2*mmst1 + 6*lt1u*mmgl*mmsb2*mmst1 -
     6*ltu*mmgl*mmsb2*mmst1 + 4*lgu*ltu*mmgl*mmsb2*mmst1 - 2*lt1u*ltu*mmgl*
     mmsb2*mmst1 - 14*invdgb2*mmgl*mmsb22*mmst1 + 12*invdgb2*lgu*mmgl*mmsb22*
     mmst1 - 6*invdgb2*lt1u*mmgl*mmsb22*mmst1 + 6*invdgb2*ltu*mmgl*mmsb22*mmst1
      - 4*invdgb2*lgu*ltu*mmgl*mmsb22*mmst1 + 2*invdgb2*lt1u*ltu*mmgl*mmsb22*
     mmst1 - 14*mmgl*mmsb1*mmt + 6*lgu*mmgl*mmsb1*mmt + 6*ltu*mmgl*mmsb1*mmt -
     2*lgu*ltu*mmgl*mmsb1*mmt + 14*invdgb1*mmgl*mmsb12*mmt - 6*invdgb1*lgu*mmgl
     *mmsb12*mmt - 6*invdgb1*ltu*mmgl*mmsb12*mmt + 2*invdgb1*lgu*ltu*mmgl*
     mmsb12*mmt + 14*mmgl*mmsb2*mmt - 6*lgu*mmgl*mmsb2*mmt - 6*ltu*mmgl*mmsb2*
     mmt + 2*lgu*ltu*mmgl*mmsb2*mmt - 14*invdgb2*mmgl*mmsb22*mmt + 6*invdgb2*
     lgu*mmgl*mmsb22*mmt + 6*invdgb2*ltu*mmgl*mmsb22*mmt - 2*invdgb2*lgu*ltu*
     mmgl*mmsb22*mmt + 56*invdgb1*mmgl*mmsb1*mmst1*mmt - 6*invdgb1*lgu*mmgl*
     mmsb1*mmst1*mmt - 6*invdgb1*lt1u*mmgl*mmsb1*mmst1*mmt - 4*invdgb1*lgu*lt1u
     *mmgl*mmsb1*mmst1*mmt - 36*invdgb1*ltu*mmgl*mmsb1*mmst1*mmt + 6*invdgb1*
     lgu*ltu*mmgl*mmsb1*mmst1*mmt + 6*invdgb1*lt1u*ltu*mmgl*mmsb1*mmst1*mmt -
     56*invdgb2*mmgl*mmsb2*mmst1*mmt + 6*invdgb2*lgu*mmgl*mmsb2*mmst1*mmt + 6*
     invdgb2*lt1u*mmgl*mmsb2*mmst1*mmt + 4*invdgb2*lgu*lt1u*mmgl*mmsb2*mmst1*
     mmt + 36*invdgb2*ltu*mmgl*mmsb2*mmst1*mmt - 6*invdgb2*lgu*ltu*mmgl*mmsb2*
     mmst1*mmt - 6*invdgb2*lt1u*ltu*mmgl*mmsb2*mmst1*mmt + 2*mmgl2*mmsb1*zt2 +
     2*mmgl*mmsb12*zt2 - 2*mmgl2*mmsb2*zt2 - 2*mmgl*mmsb22*zt2 - 2*mmgl*mmsb1*
     mmst1*zt2 + 2*invdgb1*mmgl*mmsb12*mmst1*zt2 + 2*mmgl*mmsb2*mmst1*zt2 - 2*
     invdgb2*mmgl*mmsb22*mmst1*zt2 - 2*mmgl*mmsb1*mmt*zt2 + 2*invdgb1*mmgl*
     mmsb12*mmt*zt2 + 2*mmgl*mmsb2*mmt*zt2 - 2*invdgb2*mmgl*mmsb22*mmt*zt2 + 8*
     invdgb1*mmgl*mmsb1*mmst1*mmt*zt2 - 8*invdgb2*mmgl*mmsb2*mmst1*mmt*zt2 +
     mmgl2*mmsb1*pow2(ltu) + mmgl*mmsb12*pow2(ltu) - mmgl2*mmsb2*pow2(ltu) -
     mmgl*mmsb22*pow2(ltu) - mmgl*mmsb1*mmst1*pow2(ltu) + invdgb1*mmgl*mmsb12*
     mmst1*pow2(ltu) + mmgl*mmsb2*mmst1*pow2(ltu) - invdgb2*mmgl*mmsb22*mmst1*
     pow2(ltu) - mmgl*mmsb1*mmt*pow2(ltu) + invdgb1*mmgl*mmsb12*mmt*pow2(ltu) +
     mmgl*mmsb2*mmt*pow2(ltu) - invdgb2*mmgl*mmsb22*mmt*pow2(ltu) + 6*invdgb1*
     mmgl*mmsb1*mmst1*mmt*pow2(ltu) - 6*invdgb2*mmgl*mmsb2*mmst1*mmt*pow2(ltu)
     + 14*invdgb1*mmgl*mmsb1*pow2(mmst1) + 6*invdgb1*lgu*mmgl*mmsb1*pow2(mmst1)
     - 12*invdgb1*lt1u*mmgl*mmsb1*pow2(mmst1) - 6*invdgb1*ltu*mmgl*mmsb1*pow2(
     mmst1) - 2*invdgb1*lgu*ltu*mmgl*mmsb1*pow2(mmst1) + 4*invdgb1*lt1u*ltu*
     mmgl*mmsb1*pow2(mmst1) - 14*invdgb2*mmgl*mmsb2*pow2(mmst1) - 6*invdgb2*lgu
     *mmgl*mmsb2*pow2(mmst1) + 12*invdgb2*lt1u*mmgl*mmsb2*pow2(mmst1) + 6*
     invdgb2*ltu*mmgl*mmsb2*pow2(mmst1) + 2*invdgb2*lgu*ltu*mmgl*mmsb2*pow2(
     mmst1) - 4*invdgb2*lt1u*ltu*mmgl*mmsb2*pow2(mmst1) + 14*invdgb1*mmgl*mmt*
     pow2(mmst1) - 14*invdgb2*mmgl*mmt*pow2(mmst1) - 6*invdgb1*lt1u*mmgl*mmt*
     pow2(mmst1) + 6*invdgb2*lt1u*mmgl*mmt*pow2(mmst1) - 6*invdgb1*ltu*mmgl*mmt
     *pow2(mmst1) + 6*invdgb2*ltu*mmgl*mmt*pow2(mmst1) + 2*invdgb1*lt1u*ltu*
     mmgl*mmt*pow2(mmst1) - 2*invdgb2*lt1u*ltu*mmgl*mmt*pow2(mmst1) + 2*invdgb1
     *mmgl*mmsb1*zt2*pow2(mmst1) - 2*invdgb2*mmgl*mmsb2*zt2*pow2(mmst1) + 2*
     invdgb1*mmgl*mmt*zt2*pow2(mmst1) - 2*invdgb2*mmgl*mmt*zt2*pow2(mmst1) +
     invdgb1*mmgl*mmsb1*pow2(ltu)*pow2(mmst1) - invdgb2*mmgl*mmsb2*pow2(ltu)*
     pow2(mmst1) + invdgb1*mmgl*mmt*pow2(ltu)*pow2(mmst1) - invdgb2*mmgl*mmt*
     pow2(ltu)*pow2(mmst1) - 14*invdgb1*mmgl*pow3(mmsb1) + 6*invdgb1*lgu*mmgl*
     pow3(mmsb1) + 6*invdgb1*ltu*mmgl*pow3(mmsb1) - 2*invdgb1*lgu*ltu*mmgl*pow3
     (mmsb1) - 2*invdgb1*mmgl*zt2*pow3(mmsb1) - invdgb1*mmgl*pow2(ltu)*pow3(
     mmsb1) + 14*invdgb2*mmgl*pow3(mmsb2) - 6*invdgb2*lgu*mmgl*pow3(mmsb2) - 6*
     invdgb2*ltu*mmgl*pow3(mmsb2) + 2*invdgb2*lgu*ltu*mmgl*pow3(mmsb2) + 2*
     invdgb2*mmgl*zt2*pow3(mmsb2) + invdgb2*mmgl*pow2(ltu)*pow3(mmsb2) + pow2(
     lgu)*(mmgl2*mmsb1 + mmgl*mmsb12 - mmgl2*mmsb2 - mmgl*mmsb22 - 2*mmgl*mmsb1
     *mmst1 + 2*invdgb1*mmgl*mmsb12*mmst1 + 2*mmgl*mmsb2*mmst1 - 2*invdgb2*mmgl
     *mmsb22*mmst1 - mmgl*mmsb1*mmt + invdgb1*mmgl*mmsb12*mmt + mmgl*mmsb2*mmt
     - invdgb2*mmgl*mmsb22*mmt + invdgb1*mmgl*mmsb1*mmst1*mmt - invdgb2*mmgl*
     mmsb2*mmst1*mmt + (-(invdgb1*mmgl*mmsb1) + invdgb2*mmgl*mmsb2)*pow2(mmst1)
     - invdgb1*mmgl*pow3(mmsb1) + invdgb2*mmgl*pow3(mmsb2)) - 14*invdgb1*mmgl*
     pow3(mmst1) + 14*invdgb2*mmgl*pow3(mmst1) + 6*invdgb1*lt1u*mmgl*pow3(mmst1
     ) - 6*invdgb2*lt1u*mmgl*pow3(mmst1) + 6*invdgb1*ltu*mmgl*pow3(mmst1) - 6*
     invdgb2*ltu*mmgl*pow3(mmst1) - 2*invdgb1*lt1u*ltu*mmgl*pow3(mmst1) + 2*
     invdgb2*lt1u*ltu*mmgl*pow3(mmst1) - 2*invdgb1*mmgl*zt2*pow3(mmst1) + 2*
     invdgb2*mmgl*zt2*pow3(mmst1) - invdgb1*mmgl*pow2(ltu)*pow3(mmst1) +
     invdgb2*mmgl*pow2(ltu)*pow3(mmst1) + mmgl*pow2(lt1u)*(mmst1*(mmsb1 -
     invdgb1*mmsb12 - mmsb2 + invdgb2*mmsb22 + invdgb1*mmsb1*mmt - invdgb2*
     mmsb2*mmt) + (invdgb1*(2*mmsb1 + mmt) - invdgb2*(2*mmsb2 + mmt))*pow2(
     mmst1) + (-invdgb1 + invdgb2)*pow3(mmst1))))/(3.*mb*mgl) - (2*(28*mmgl2 -
     12*lgu*mmgl2 - 12*ltu*mmgl2 + 4*lgu*ltu*mmgl2 + 28*mmgl*mmsb1 - 12*lgu*
     mmgl*mmsb1 - 12*ltu*mmgl*mmsb1 + 4*lgu*ltu*mmgl*mmsb1 + 42*mmsb12 - 18*lgu
     *mmsb12 - 18*ltu*mmsb12 + 6*lgu*ltu*mmsb12 + 28*mmgl*mmsb2 - 12*lgu*mmgl*
     mmsb2 - 12*ltu*mmgl*mmsb2 + 4*lgu*ltu*mmgl*mmsb2 + 42*mmsb22 - 18*lgu*
     mmsb22 - 18*ltu*mmsb22 + 6*lgu*ltu*mmsb22 - 28*mmgl*mmst1 + 24*lgu*mmgl*
     mmst1 - 12*lt1u*mmgl*mmst1 + 12*ltu*mmgl*mmst1 - 8*lgu*ltu*mmgl*mmst1 + 4*
     lt1u*ltu*mmgl*mmst1 - 28*mmsb1*mmst1 + 24*lgu*mmsb1*mmst1 - 12*lt1u*mmsb1*
     mmst1 + 12*ltu*mmsb1*mmst1 - 8*lgu*ltu*mmsb1*mmst1 + 4*lt1u*ltu*mmsb1*
     mmst1 + 42*invdgb1*mmsb12*mmst1 - 36*invdgb1*lgu*mmsb12*mmst1 + 18*invdgb1
     *lt1u*mmsb12*mmst1 - 18*invdgb1*ltu*mmsb12*mmst1 + 12*invdgb1*lgu*ltu*
     mmsb12*mmst1 - 6*invdgb1*lt1u*ltu*mmsb12*mmst1 - 28*mmsb2*mmst1 + 24*lgu*
     mmsb2*mmst1 - 12*lt1u*mmsb2*mmst1 + 12*ltu*mmsb2*mmst1 - 8*lgu*ltu*mmsb2*
     mmst1 + 4*lt1u*ltu*mmsb2*mmst1 + 42*invdgb2*mmsb22*mmst1 - 36*invdgb2*lgu*
     mmsb22*mmst1 + 18*invdgb2*lt1u*mmsb22*mmst1 - 18*invdgb2*ltu*mmsb22*mmst1
     + 12*invdgb2*lgu*ltu*mmsb22*mmst1 - 6*invdgb2*lt1u*ltu*mmsb22*mmst1 - 28*
     mmgl*mmt + 12*lgu*mmgl*mmt + 12*ltu*mmgl*mmt - 4*lgu*ltu*mmgl*mmt - 28*
     mmsb1*mmt + 12*lgu*mmsb1*mmt + 12*ltu*mmsb1*mmt - 4*lgu*ltu*mmsb1*mmt + 42
     *invdgb1*mmsb12*mmt - 18*invdgb1*lgu*mmsb12*mmt - 18*invdgb1*ltu*mmsb12*
     mmt + 6*invdgb1*lgu*ltu*mmsb12*mmt - 28*mmsb2*mmt + 12*lgu*mmsb2*mmt + 12*
     ltu*mmsb2*mmt - 4*lgu*ltu*mmsb2*mmt + 42*invdgb2*mmsb22*mmt - 18*invdgb2*
     lgu*mmsb22*mmt - 18*invdgb2*ltu*mmsb22*mmt + 6*invdgb2*lgu*ltu*mmsb22*mmt
     - 112*mmst1*mmt + 12*lgu*mmst1*mmt + 12*lt1u*mmst1*mmt + 8*lgu*lt1u*mmst1*
     mmt + 72*ltu*mmst1*mmt - 12*lgu*ltu*mmst1*mmt - 12*lt1u*ltu*mmst1*mmt +
     112*invdgb1*mmsb1*mmst1*mmt - 12*invdgb1*lgu*mmsb1*mmst1*mmt - 12*invdgb1*
     lt1u*mmsb1*mmst1*mmt - 8*invdgb1*lgu*lt1u*mmsb1*mmst1*mmt - 72*invdgb1*ltu
     *mmsb1*mmst1*mmt + 12*invdgb1*lgu*ltu*mmsb1*mmst1*mmt + 12*invdgb1*lt1u*
     ltu*mmsb1*mmst1*mmt + 112*invdgb2*mmsb2*mmst1*mmt - 12*invdgb2*lgu*mmsb2*
     mmst1*mmt - 12*invdgb2*lt1u*mmsb2*mmst1*mmt - 8*invdgb2*lgu*lt1u*mmsb2*
     mmst1*mmt - 72*invdgb2*ltu*mmsb2*mmst1*mmt + 12*invdgb2*lgu*ltu*mmsb2*
     mmst1*mmt + 12*invdgb2*lt1u*ltu*mmsb2*mmst1*mmt + 4*mmgl2*zt2 + 4*mmgl*
     mmsb1*zt2 + 6*mmsb12*zt2 + 4*mmgl*mmsb2*zt2 + 6*mmsb22*zt2 - 4*mmgl*mmst1*
     zt2 - 4*mmsb1*mmst1*zt2 + 6*invdgb1*mmsb12*mmst1*zt2 - 4*mmsb2*mmst1*zt2 +
     6*invdgb2*mmsb22*mmst1*zt2 - 4*mmgl*mmt*zt2 - 4*mmsb1*mmt*zt2 + 6*invdgb1*
     mmsb12*mmt*zt2 - 4*mmsb2*mmt*zt2 + 6*invdgb2*mmsb22*mmt*zt2 - 16*mmst1*mmt
     *zt2 + 16*invdgb1*mmsb1*mmst1*mmt*zt2 + 16*invdgb2*mmsb2*mmst1*mmt*zt2 + 2
     *mmgl2*pow2(lgu) + 2*mmgl*mmsb1*pow2(lgu) + 3*mmsb12*pow2(lgu) + 2*mmgl*
     mmsb2*pow2(lgu) + 3*mmsb22*pow2(lgu) - 4*mmgl*mmst1*pow2(lgu) - 4*mmsb1*
     mmst1*pow2(lgu) + 6*invdgb1*mmsb12*mmst1*pow2(lgu) - 4*mmsb2*mmst1*pow2(
     lgu) + 6*invdgb2*mmsb22*mmst1*pow2(lgu) - 2*mmgl*mmt*pow2(lgu) - 2*mmsb1*
     mmt*pow2(lgu) + 3*invdgb1*mmsb12*mmt*pow2(lgu) - 2*mmsb2*mmt*pow2(lgu) + 3
     *invdgb2*mmsb22*mmt*pow2(lgu) - 2*mmst1*mmt*pow2(lgu) + 2*invdgb1*mmsb1*
     mmst1*mmt*pow2(lgu) + 2*invdgb2*mmsb2*mmst1*mmt*pow2(lgu) + 2*mmgl*mmst1*
     pow2(lt1u) + 2*mmsb1*mmst1*pow2(lt1u) - 3*invdgb1*mmsb12*mmst1*pow2(lt1u)
     + 2*mmsb2*mmst1*pow2(lt1u) - 3*invdgb2*mmsb22*mmst1*pow2(lt1u) - 2*mmst1*
     mmt*pow2(lt1u) + 2*invdgb1*mmsb1*mmst1*mmt*pow2(lt1u) + 2*invdgb2*mmsb2*
     mmst1*mmt*pow2(lt1u) + 2*mmgl2*pow2(ltu) + 2*mmgl*mmsb1*pow2(ltu) + 3*
     mmsb12*pow2(ltu) + 2*mmgl*mmsb2*pow2(ltu) + 3*mmsb22*pow2(ltu) - 2*mmgl*
     mmst1*pow2(ltu) - 2*mmsb1*mmst1*pow2(ltu) + 3*invdgb1*mmsb12*mmst1*pow2(
     ltu) - 2*mmsb2*mmst1*pow2(ltu) + 3*invdgb2*mmsb22*mmst1*pow2(ltu) - 2*mmgl
     *mmt*pow2(ltu) - 2*mmsb1*mmt*pow2(ltu) + 3*invdgb1*mmsb12*mmt*pow2(ltu) -
     2*mmsb2*mmt*pow2(ltu) + 3*invdgb2*mmsb22*mmt*pow2(ltu) - 12*mmst1*mmt*pow2
     (ltu) + 12*invdgb1*mmsb1*mmst1*mmt*pow2(ltu) + 12*invdgb2*mmsb2*mmst1*mmt*
     pow2(ltu) - 28*pow2(mmst1) - 12*lgu*pow2(mmst1) + 24*lt1u*pow2(mmst1) + 12
     *ltu*pow2(mmst1) + 4*lgu*ltu*pow2(mmst1) - 8*lt1u*ltu*pow2(mmst1) + 28*
     invdgb1*mmsb1*pow2(mmst1) + 12*invdgb1*lgu*mmsb1*pow2(mmst1) - 24*invdgb1*
     lt1u*mmsb1*pow2(mmst1) - 12*invdgb1*ltu*mmsb1*pow2(mmst1) - 4*invdgb1*lgu*
     ltu*mmsb1*pow2(mmst1) + 8*invdgb1*lt1u*ltu*mmsb1*pow2(mmst1) + 28*invdgb2*
     mmsb2*pow2(mmst1) + 12*invdgb2*lgu*mmsb2*pow2(mmst1) - 24*invdgb2*lt1u*
     mmsb2*pow2(mmst1) - 12*invdgb2*ltu*mmsb2*pow2(mmst1) - 4*invdgb2*lgu*ltu*
     mmsb2*pow2(mmst1) + 8*invdgb2*lt1u*ltu*mmsb2*pow2(mmst1) + 14*invdgb1*mmt*
     pow2(mmst1) + 14*invdgb2*mmt*pow2(mmst1) - 6*invdgb1*lt1u*mmt*pow2(mmst1)
     - 6*invdgb2*lt1u*mmt*pow2(mmst1) - 6*invdgb1*ltu*mmt*pow2(mmst1) - 6*
     invdgb2*ltu*mmt*pow2(mmst1) + 2*invdgb1*lt1u*ltu*mmt*pow2(mmst1) + 2*
     invdgb2*lt1u*ltu*mmt*pow2(mmst1) - 4*zt2*pow2(mmst1) + 4*invdgb1*mmsb1*zt2
     *pow2(mmst1) + 4*invdgb2*mmsb2*zt2*pow2(mmst1) + 2*invdgb1*mmt*zt2*pow2(
     mmst1) + 2*invdgb2*mmt*zt2*pow2(mmst1) + 2*pow2(lgu)*pow2(mmst1) - 2*
     invdgb1*mmsb1*pow2(lgu)*pow2(mmst1) - 2*invdgb2*mmsb2*pow2(lgu)*pow2(mmst1
     ) - 4*pow2(lt1u)*pow2(mmst1) + 4*invdgb1*mmsb1*pow2(lt1u)*pow2(mmst1) + 4*
     invdgb2*mmsb2*pow2(lt1u)*pow2(mmst1) + invdgb1*mmt*pow2(lt1u)*pow2(mmst1)
     + invdgb2*mmt*pow2(lt1u)*pow2(mmst1) - 2*pow2(ltu)*pow2(mmst1) + 2*invdgb1
     *mmsb1*pow2(ltu)*pow2(mmst1) + 2*invdgb2*mmsb2*pow2(ltu)*pow2(mmst1) +
     invdgb1*mmt*pow2(ltu)*pow2(mmst1) + invdgb2*mmt*pow2(ltu)*pow2(mmst1) - 56
     *invdgb1*pow3(mmsb1) + 24*invdgb1*lgu*pow3(mmsb1) + 24*invdgb1*ltu*pow3(
     mmsb1) - 8*invdgb1*lgu*ltu*pow3(mmsb1) - 8*invdgb1*zt2*pow3(mmsb1) - 4*
     invdgb1*pow2(lgu)*pow3(mmsb1) - 4*invdgb1*pow2(ltu)*pow3(mmsb1) - 56*
     invdgb2*pow3(mmsb2) + 24*invdgb2*lgu*pow3(mmsb2) + 24*invdgb2*ltu*pow3(
     mmsb2) - 8*invdgb2*lgu*ltu*pow3(mmsb2) - 8*invdgb2*zt2*pow3(mmsb2) - 4*
     invdgb2*pow2(lgu)*pow3(mmsb2) - 4*invdgb2*pow2(ltu)*pow3(mmsb2) - 14*
     invdgb1*pow3(mmst1) - 14*invdgb2*pow3(mmst1) + 6*invdgb1*lt1u*pow3(mmst1)
     + 6*invdgb2*lt1u*pow3(mmst1) + 6*invdgb1*ltu*pow3(mmst1) + 6*invdgb2*ltu*
     pow3(mmst1) - 2*invdgb1*lt1u*ltu*pow3(mmst1) - 2*invdgb2*lt1u*ltu*pow3(
     mmst1) - 2*invdgb1*zt2*pow3(mmst1) - 2*invdgb2*zt2*pow3(mmst1) - invdgb1*
     pow2(lt1u)*pow3(mmst1) - invdgb2*pow2(lt1u)*pow3(mmst1) - invdgb1*pow2(ltu
     )*pow3(mmst1) - invdgb2*pow2(ltu)*pow3(mmst1) + pow2(invdgb1)*(-56*mmsb12*
     mmst1*mmt + 6*lgu*mmsb12*mmst1*mmt + 6*lt1u*mmsb12*mmst1*mmt + 4*lgu*lt1u*
     mmsb12*mmst1*mmt + 36*ltu*mmsb12*mmst1*mmt - 6*lgu*ltu*mmsb12*mmst1*mmt -
     6*lt1u*ltu*mmsb12*mmst1*mmt - 8*mmsb12*mmst1*mmt*zt2 - 6*mmsb12*mmst1*mmt*
     pow2(ltu) - 14*mmsb12*pow2(mmst1) - 6*lgu*mmsb12*pow2(mmst1) + 12*lt1u*
     mmsb12*pow2(mmst1) + 6*ltu*mmsb12*pow2(mmst1) + 2*lgu*ltu*mmsb12*pow2(
     mmst1) - 4*lt1u*ltu*mmsb12*pow2(mmst1) - 14*mmsb1*mmt*pow2(mmst1) + 6*lt1u
     *mmsb1*mmt*pow2(mmst1) + 6*ltu*mmsb1*mmt*pow2(mmst1) - 2*lt1u*ltu*mmsb1*
     mmt*pow2(mmst1) - 2*mmsb12*zt2*pow2(mmst1) - 2*mmsb1*mmt*zt2*pow2(mmst1) -
     mmsb12*pow2(ltu)*pow2(mmst1) - mmsb1*mmt*pow2(ltu)*pow2(mmst1) - 14*mmst1*
     pow3(mmsb1) + 12*lgu*mmst1*pow3(mmsb1) - 6*lt1u*mmst1*pow3(mmsb1) + 6*ltu*
     mmst1*pow3(mmsb1) - 4*lgu*ltu*mmst1*pow3(mmsb1) + 2*lt1u*ltu*mmst1*pow3(
     mmsb1) - 14*mmt*pow3(mmsb1) + 6*lgu*mmt*pow3(mmsb1) + 6*ltu*mmt*pow3(mmsb1
     ) - 2*lgu*ltu*mmt*pow3(mmsb1) - 2*mmst1*zt2*pow3(mmsb1) - 2*mmt*zt2*pow3(
     mmsb1) - mmst1*pow2(ltu)*pow3(mmsb1) - mmt*pow2(ltu)*pow3(mmsb1) + 14*
     mmsb1*pow3(mmst1) - 6*lt1u*mmsb1*pow3(mmst1) - 6*ltu*mmsb1*pow3(mmst1) + 2
     *lt1u*ltu*mmsb1*pow3(mmst1) + 2*mmsb1*zt2*pow3(mmst1) + mmsb1*pow2(ltu)*
     pow3(mmst1) + pow2(lt1u)*(-(mmsb12*mmst1*mmt) - (2*mmsb12 + mmsb1*mmt)*
     pow2(mmst1) + mmst1*pow3(mmsb1) + mmsb1*pow3(mmst1)) + 14*pow4(mmsb1) - 6*
     lgu*pow4(mmsb1) - 6*ltu*pow4(mmsb1) + 2*lgu*ltu*pow4(mmsb1) + 2*zt2*pow4(
     mmsb1) + pow2(ltu)*pow4(mmsb1) + pow2(lgu)*(-(mmsb12*mmst1*mmt) + mmsb12*
     pow2(mmst1) - (2*mmst1 + mmt)*pow3(mmsb1) + pow4(mmsb1))) + pow2(invdgb2)*
     (-56*mmsb22*mmst1*mmt + 6*lgu*mmsb22*mmst1*mmt + 6*lt1u*mmsb22*mmst1*mmt +
     4*lgu*lt1u*mmsb22*mmst1*mmt + 36*ltu*mmsb22*mmst1*mmt - 6*lgu*ltu*mmsb22*
     mmst1*mmt - 6*lt1u*ltu*mmsb22*mmst1*mmt - 8*mmsb22*mmst1*mmt*zt2 - 6*
     mmsb22*mmst1*mmt*pow2(ltu) - 14*mmsb22*pow2(mmst1) - 6*lgu*mmsb22*pow2(
     mmst1) + 12*lt1u*mmsb22*pow2(mmst1) + 6*ltu*mmsb22*pow2(mmst1) + 2*lgu*ltu
     *mmsb22*pow2(mmst1) - 4*lt1u*ltu*mmsb22*pow2(mmst1) - 14*mmsb2*mmt*pow2(
     mmst1) + 6*lt1u*mmsb2*mmt*pow2(mmst1) + 6*ltu*mmsb2*mmt*pow2(mmst1) - 2*
     lt1u*ltu*mmsb2*mmt*pow2(mmst1) - 2*mmsb22*zt2*pow2(mmst1) - 2*mmsb2*mmt*
     zt2*pow2(mmst1) - mmsb22*pow2(ltu)*pow2(mmst1) - mmsb2*mmt*pow2(ltu)*pow2(
     mmst1) - 14*mmst1*pow3(mmsb2) + 12*lgu*mmst1*pow3(mmsb2) - 6*lt1u*mmst1*
     pow3(mmsb2) + 6*ltu*mmst1*pow3(mmsb2) - 4*lgu*ltu*mmst1*pow3(mmsb2) + 2*
     lt1u*ltu*mmst1*pow3(mmsb2) - 14*mmt*pow3(mmsb2) + 6*lgu*mmt*pow3(mmsb2) +
     6*ltu*mmt*pow3(mmsb2) - 2*lgu*ltu*mmt*pow3(mmsb2) - 2*mmst1*zt2*pow3(mmsb2
     ) - 2*mmt*zt2*pow3(mmsb2) - mmst1*pow2(ltu)*pow3(mmsb2) - mmt*pow2(ltu)*
     pow3(mmsb2) + 14*mmsb2*pow3(mmst1) - 6*lt1u*mmsb2*pow3(mmst1) - 6*ltu*
     mmsb2*pow3(mmst1) + 2*lt1u*ltu*mmsb2*pow3(mmst1) + 2*mmsb2*zt2*pow3(mmst1)
     + mmsb2*pow2(ltu)*pow3(mmst1) + pow2(lt1u)*(-(mmsb22*mmst1*mmt) - (2*
     mmsb22 + mmsb2*mmt)*pow2(mmst1) + mmst1*pow3(mmsb2) + mmsb2*pow3(mmst1)) +
     14*pow4(mmsb2) - 6*lgu*pow4(mmsb2) - 6*ltu*pow4(mmsb2) + 2*lgu*ltu*pow4(
     mmsb2) + 2*zt2*pow4(mmsb2) + pow2(ltu)*pow4(mmsb2) + pow2(lgu)*(-(mmsb22*
     mmst1*mmt) + mmsb22*pow2(mmst1) - (2*mmst1 + mmt)*pow3(mmsb2) + pow4(mmsb2
     )))))/3.) + DeltaInv(mmt,mmsb2,mmst2)*((-2*mmgl*s2t*pow2(invdgb2)*(-14*
     mmsb22*mmst2*mmt + 6*lb2u*mmsb22*mmst2*mmt - 6*lt2u*mmsb22*mmst2*mmt + 2*
     lb2u*lt2u*mmsb22*mmst2*mmt + 12*ltu*mmsb22*mmst2*mmt - 4*lb2u*ltu*mmsb22*
     mmst2*mmt - 2*mmsb22*mmst2*mmt*zt2 - 2*mmsb22*mmst2*mmt*pow2(ltu) + 6*lt2u
     *mmsb2*mmt*pow2(mmst2) - 2*lb2u*lt2u*mmsb2*mmt*pow2(mmst2) - 6*ltu*mmsb2*
     mmt*pow2(mmst2) + 2*lb2u*ltu*mmsb2*mmt*pow2(mmst2) + mmsb2*mmt*pow2(ltu)*
     pow2(mmst2) - 14*mmsb22*pow2(mmt) + 6*lb2u*mmsb22*pow2(mmt) + 6*ltu*mmsb22
     *pow2(mmt) - 2*lb2u*ltu*mmsb22*pow2(mmt) - 28*mmsb2*mmst2*pow2(mmt) + 6*
     lt2u*mmsb2*mmst2*pow2(mmt) + 2*lb2u*lt2u*mmsb2*mmst2*pow2(mmt) + 18*ltu*
     mmsb2*mmst2*pow2(mmt) - 2*lb2u*ltu*mmsb2*mmst2*pow2(mmt) - 4*lt2u*ltu*
     mmsb2*mmst2*pow2(mmt) - 2*mmsb22*zt2*pow2(mmt) - 4*mmsb2*mmst2*zt2*pow2(
     mmt) - mmsb22*pow2(ltu)*pow2(mmt) - 3*mmsb2*mmst2*pow2(ltu)*pow2(mmt) +
     pow2(lt2u)*(mmsb22*mmst2*mmt - mmsb2*mmt*pow2(mmst2) - mmsb2*mmst2*pow2(
     mmt)) + 14*mmt*pow3(mmsb2) - 6*lb2u*mmt*pow3(mmsb2) - 6*ltu*mmt*pow3(mmsb2
     ) + 2*lb2u*ltu*mmt*pow3(mmsb2) + 2*mmt*zt2*pow3(mmsb2) + mmt*pow2(ltu)*
     pow3(mmsb2) - pow2(lb2u)*(mmsb22*mmst2*mmt + mmsb22*pow2(mmt) - mmt*pow3(
     mmsb2))))/(3.*mgl*mt) + (-(invdgb2*(14*mmsb22*mmst2 - 12*lb2u*mmsb22*mmst2
      + 6*lt2u*mmsb22*mmst2 - 6*ltu*mmsb22*mmst2 + 4*lb2u*ltu*mmsb22*mmst2 - 2*
     lt2u*ltu*mmsb22*mmst2 + 14*mmsb22*mmt - 6*lb2u*mmsb22*mmt - 6*ltu*mmsb22*
     mmt + 2*lb2u*ltu*mmsb22*mmt + 56*mmsb2*mmst2*mmt - 6*lb2u*mmsb2*mmst2*mmt
     - 6*lt2u*mmsb2*mmst2*mmt - 4*lb2u*lt2u*mmsb2*mmst2*mmt - 36*ltu*mmsb2*
     mmst2*mmt + 6*lb2u*ltu*mmsb2*mmst2*mmt + 6*lt2u*ltu*mmsb2*mmst2*mmt + 2*
     mmsb22*mmst2*zt2 + 2*mmsb22*mmt*zt2 + 8*mmsb2*mmst2*mmt*zt2 + mmsb22*mmst2
     *pow2(ltu) + mmsb22*mmt*pow2(ltu) + 6*mmsb2*mmst2*mmt*pow2(ltu) + 14*mmsb2
     *pow2(mmst2) + 6*lb2u*mmsb2*pow2(mmst2) - 12*lt2u*mmsb2*pow2(mmst2) - 6*
     ltu*mmsb2*pow2(mmst2) - 2*lb2u*ltu*mmsb2*pow2(mmst2) + 4*lt2u*ltu*mmsb2*
     pow2(mmst2) + 14*mmt*pow2(mmst2) - 6*lt2u*mmt*pow2(mmst2) - 6*ltu*mmt*pow2
     (mmst2) + 2*lt2u*ltu*mmt*pow2(mmst2) + 2*mmsb2*zt2*pow2(mmst2) + 2*mmt*zt2
     *pow2(mmst2) + mmsb2*pow2(ltu)*pow2(mmst2) + mmt*pow2(ltu)*pow2(mmst2) +
     pow2(lb2u)*(2*mmsb22*mmst2 + mmsb22*mmt + mmsb2*mmst2*mmt - mmsb2*pow2(
     mmst2) - pow3(mmsb2)) - 14*pow3(mmsb2) + 6*lb2u*pow3(mmsb2) + 6*ltu*pow3(
     mmsb2) - 2*lb2u*ltu*pow3(mmsb2) - 2*zt2*pow3(mmsb2) - pow2(ltu)*pow3(mmsb2
     ) + pow2(lt2u)*(-(mmsb22*mmst2) + mmsb2*mmst2*mmt + (2*mmsb2 + mmt)*pow2(
     mmst2) - pow3(mmst2)) - 14*pow3(mmst2) + 6*lt2u*pow3(mmst2) + 6*ltu*pow3(
     mmst2) - 2*lt2u*ltu*pow3(mmst2) - 2*zt2*pow3(mmst2) - pow2(ltu)*pow3(mmst2
     ))) - 2*pow2(invdgb2)*(-56*mmsb22*mmst2*mmt + 6*lb2u*mmsb22*mmst2*mmt + 6*
     lt2u*mmsb22*mmst2*mmt + 4*lb2u*lt2u*mmsb22*mmst2*mmt + 36*ltu*mmsb22*mmst2
     *mmt - 6*lb2u*ltu*mmsb22*mmst2*mmt - 6*lt2u*ltu*mmsb22*mmst2*mmt - 8*
     mmsb22*mmst2*mmt*zt2 - 6*mmsb22*mmst2*mmt*pow2(ltu) - 14*mmsb22*pow2(mmst2
     ) - 6*lb2u*mmsb22*pow2(mmst2) + 12*lt2u*mmsb22*pow2(mmst2) + 6*ltu*mmsb22*
     pow2(mmst2) + 2*lb2u*ltu*mmsb22*pow2(mmst2) - 4*lt2u*ltu*mmsb22*pow2(mmst2
     ) - 14*mmsb2*mmt*pow2(mmst2) + 6*lt2u*mmsb2*mmt*pow2(mmst2) + 6*ltu*mmsb2*
     mmt*pow2(mmst2) - 2*lt2u*ltu*mmsb2*mmt*pow2(mmst2) - 2*mmsb22*zt2*pow2(
     mmst2) - 2*mmsb2*mmt*zt2*pow2(mmst2) - mmsb22*pow2(ltu)*pow2(mmst2) -
     mmsb2*mmt*pow2(ltu)*pow2(mmst2) - 14*mmst2*pow3(mmsb2) + 12*lb2u*mmst2*
     pow3(mmsb2) - 6*lt2u*mmst2*pow3(mmsb2) + 6*ltu*mmst2*pow3(mmsb2) - 4*lb2u*
     ltu*mmst2*pow3(mmsb2) + 2*lt2u*ltu*mmst2*pow3(mmsb2) - 14*mmt*pow3(mmsb2)
     + 6*lb2u*mmt*pow3(mmsb2) + 6*ltu*mmt*pow3(mmsb2) - 2*lb2u*ltu*mmt*pow3(
     mmsb2) - 2*mmst2*zt2*pow3(mmsb2) - 2*mmt*zt2*pow3(mmsb2) - mmst2*pow2(ltu)
     *pow3(mmsb2) - mmt*pow2(ltu)*pow3(mmsb2) + 14*mmsb2*pow3(mmst2) - 6*lt2u*
     mmsb2*pow3(mmst2) - 6*ltu*mmsb2*pow3(mmst2) + 2*lt2u*ltu*mmsb2*pow3(mmst2)
     + 2*mmsb2*zt2*pow3(mmst2) + mmsb2*pow2(ltu)*pow3(mmst2) + pow2(lt2u)*(-(
     mmsb22*mmst2*mmt) - (2*mmsb22 + mmsb2*mmt)*pow2(mmst2) + mmst2*pow3(mmsb2)
     + mmsb2*pow3(mmst2)) + 14*pow4(mmsb2) - 6*lb2u*pow4(mmsb2) - 6*ltu*pow4(
     mmsb2) + 2*lb2u*ltu*pow4(mmsb2) + 2*zt2*pow4(mmsb2) + pow2(ltu)*pow4(mmsb2
     ) + pow2(lb2u)*(-(mmsb22*mmst2*mmt) + mmsb22*pow2(mmst2) - (2*mmst2 + mmt)
     *pow3(mmsb2) + pow4(mmsb2))))/3.) + DeltaInv(mmt,mmst2,mmgl)*(s2t*((4*s2b*
     (14*mmgl*mmsb1*mmt - 6*lgu*mmgl*mmsb1*mmt - 6*ltu*mmgl*mmsb1*mmt + 2*lgu*
     ltu*mmgl*mmsb1*mmt + 14*mmsb12*mmt - 6*lgu*mmsb12*mmt - 6*ltu*mmsb12*mmt +
     2*lgu*ltu*mmsb12*mmt - 14*mmgl*mmsb2*mmt + 6*lgu*mmgl*mmsb2*mmt + 6*ltu*
     mmgl*mmsb2*mmt - 2*lgu*ltu*mmgl*mmsb2*mmt - 14*mmsb22*mmt + 6*lgu*mmsb22*
     mmt + 6*ltu*mmsb22*mmt - 2*lgu*ltu*mmsb22*mmt - 14*mmsb1*mmst2*mmt + 6*lgu
     *mmsb1*mmst2*mmt - 6*lt2u*mmsb1*mmst2*mmt + 2*lgu*lt2u*mmsb1*mmst2*mmt +
     12*ltu*mmsb1*mmst2*mmt - 4*lgu*ltu*mmsb1*mmst2*mmt + 14*invdgb1*mmsb12*
     mmst2*mmt - 6*invdgb1*lgu*mmsb12*mmst2*mmt + 6*invdgb1*lt2u*mmsb12*mmst2*
     mmt - 2*invdgb1*lgu*lt2u*mmsb12*mmst2*mmt - 12*invdgb1*ltu*mmsb12*mmst2*
     mmt + 4*invdgb1*lgu*ltu*mmsb12*mmst2*mmt + 14*mmsb2*mmst2*mmt - 6*lgu*
     mmsb2*mmst2*mmt + 6*lt2u*mmsb2*mmst2*mmt - 2*lgu*lt2u*mmsb2*mmst2*mmt - 12
     *ltu*mmsb2*mmst2*mmt + 4*lgu*ltu*mmsb2*mmst2*mmt - 14*invdgb2*mmsb22*mmst2
     *mmt + 6*invdgb2*lgu*mmsb22*mmst2*mmt - 6*invdgb2*lt2u*mmsb22*mmst2*mmt +
     2*invdgb2*lgu*lt2u*mmsb22*mmst2*mmt + 12*invdgb2*ltu*mmsb22*mmst2*mmt - 4*
     invdgb2*lgu*ltu*mmsb22*mmst2*mmt + 2*mmgl*mmsb1*mmt*zt2 + 2*mmsb12*mmt*zt2
      - 2*mmgl*mmsb2*mmt*zt2 - 2*mmsb22*mmt*zt2 - 2*mmsb1*mmst2*mmt*zt2 + 2*
     invdgb1*mmsb12*mmst2*mmt*zt2 + 2*mmsb2*mmst2*mmt*zt2 - 2*invdgb2*mmsb22*
     mmst2*mmt*zt2 + mmgl*mmsb1*mmt*pow2(ltu) + mmsb12*mmt*pow2(ltu) - mmgl*
     mmsb2*mmt*pow2(ltu) - mmsb22*mmt*pow2(ltu) - 2*mmsb1*mmst2*mmt*pow2(ltu) +
     2*invdgb1*mmsb12*mmst2*mmt*pow2(ltu) + 2*mmsb2*mmst2*mmt*pow2(ltu) - 2*
     invdgb2*mmsb22*mmst2*mmt*pow2(ltu) - 6*invdgb1*lt2u*mmsb1*mmt*pow2(mmst2)
     + 2*invdgb1*lgu*lt2u*mmsb1*mmt*pow2(mmst2) + 6*invdgb1*ltu*mmsb1*mmt*pow2(
     mmst2) - 2*invdgb1*lgu*ltu*mmsb1*mmt*pow2(mmst2) + 6*invdgb2*lt2u*mmsb2*
     mmt*pow2(mmst2) - 2*invdgb2*lgu*lt2u*mmsb2*mmt*pow2(mmst2) - 6*invdgb2*ltu
     *mmsb2*mmt*pow2(mmst2) + 2*invdgb2*lgu*ltu*mmsb2*mmt*pow2(mmst2) - invdgb1
     *mmsb1*mmt*pow2(ltu)*pow2(mmst2) + invdgb2*mmsb2*mmt*pow2(ltu)*pow2(mmst2)
     - 14*mmsb1*pow2(mmt) + 6*lgu*mmsb1*pow2(mmt) + 6*ltu*mmsb1*pow2(mmt) - 2*
     lgu*ltu*mmsb1*pow2(mmt) + 14*invdgb1*mmsb12*pow2(mmt) - 6*invdgb1*lgu*
     mmsb12*pow2(mmt) - 6*invdgb1*ltu*mmsb12*pow2(mmt) + 2*invdgb1*lgu*ltu*
     mmsb12*pow2(mmt) + 14*mmsb2*pow2(mmt) - 6*lgu*mmsb2*pow2(mmt) - 6*ltu*
     mmsb2*pow2(mmt) + 2*lgu*ltu*mmsb2*pow2(mmt) - 14*invdgb2*mmsb22*pow2(mmt)
     + 6*invdgb2*lgu*mmsb22*pow2(mmt) + 6*invdgb2*ltu*mmsb22*pow2(mmt) - 2*
     invdgb2*lgu*ltu*mmsb22*pow2(mmt) + 28*invdgb1*mmsb1*mmst2*pow2(mmt) - 6*
     invdgb1*lt2u*mmsb1*mmst2*pow2(mmt) - 2*invdgb1*lgu*lt2u*mmsb1*mmst2*pow2(
     mmt) - 18*invdgb1*ltu*mmsb1*mmst2*pow2(mmt) + 2*invdgb1*lgu*ltu*mmsb1*
     mmst2*pow2(mmt) + 4*invdgb1*lt2u*ltu*mmsb1*mmst2*pow2(mmt) - 28*invdgb2*
     mmsb2*mmst2*pow2(mmt) + 6*invdgb2*lt2u*mmsb2*mmst2*pow2(mmt) + 2*invdgb2*
     lgu*lt2u*mmsb2*mmst2*pow2(mmt) + 18*invdgb2*ltu*mmsb2*mmst2*pow2(mmt) - 2*
     invdgb2*lgu*ltu*mmsb2*mmst2*pow2(mmt) - 4*invdgb2*lt2u*ltu*mmsb2*mmst2*
     pow2(mmt) - 2*mmsb1*zt2*pow2(mmt) + 2*invdgb1*mmsb12*zt2*pow2(mmt) + 2*
     mmsb2*zt2*pow2(mmt) - 2*invdgb2*mmsb22*zt2*pow2(mmt) + 4*invdgb1*mmsb1*
     mmst2*zt2*pow2(mmt) - 4*invdgb2*mmsb2*mmst2*zt2*pow2(mmt) - mmsb1*pow2(ltu
     )*pow2(mmt) + invdgb1*mmsb12*pow2(ltu)*pow2(mmt) + mmsb2*pow2(ltu)*pow2(
     mmt) - invdgb2*mmsb22*pow2(ltu)*pow2(mmt) + 3*invdgb1*mmsb1*mmst2*pow2(ltu
     )*pow2(mmt) - 3*invdgb2*mmsb2*mmst2*pow2(ltu)*pow2(mmt) + pow2(lt2u)*((
     mmsb1 - invdgb1*mmsb12 - mmsb2 + invdgb2*mmsb22)*mmst2*mmt + (invdgb1*
     mmsb1 - invdgb2*mmsb2)*mmt*pow2(mmst2) + (invdgb1*mmsb1 - invdgb2*mmsb2)*
     mmst2*pow2(mmt)) - 14*invdgb1*mmt*pow3(mmsb1) + 6*invdgb1*lgu*mmt*pow3(
     mmsb1) + 6*invdgb1*ltu*mmt*pow3(mmsb1) - 2*invdgb1*lgu*ltu*mmt*pow3(mmsb1)
     - 2*invdgb1*mmt*zt2*pow3(mmsb1) - invdgb1*mmt*pow2(ltu)*pow3(mmsb1) + 14*
     invdgb2*mmt*pow3(mmsb2) - 6*invdgb2*lgu*mmt*pow3(mmsb2) - 6*invdgb2*ltu*
     mmt*pow3(mmsb2) + 2*invdgb2*lgu*ltu*mmt*pow3(mmsb2) + 2*invdgb2*mmt*zt2*
     pow3(mmsb2) + invdgb2*mmt*pow2(ltu)*pow3(mmsb2) + pow2(lgu)*((-mmsb1 +
     invdgb1*mmsb12 + mmsb2 - invdgb2*mmsb22)*pow2(mmt) + mmt*(mmgl*mmsb1 +
     mmsb12 - mmgl*mmsb2 - mmsb22 - mmsb1*mmst2 + invdgb1*mmsb12*mmst2 + mmsb2*
     mmst2 - invdgb2*mmsb22*mmst2 - invdgb1*pow3(mmsb1) + invdgb2*pow3(mmsb2)))
     ))/(3.*mb*mt) - (2*(28*mmgl2*mmt - 12*lgu*mmgl2*mmt - 12*ltu*mmgl2*mmt + 4
     *lgu*ltu*mmgl2*mmt + 28*mmgl*mmsb1*mmt - 12*lgu*mmgl*mmsb1*mmt - 12*ltu*
     mmgl*mmsb1*mmt + 4*lgu*ltu*mmgl*mmsb1*mmt - 42*invdgb1*mmgl*mmsb12*mmt +
     18*invdgb1*lgu*mmgl*mmsb12*mmt + 18*invdgb1*ltu*mmgl*mmsb12*mmt - 6*
     invdgb1*lgu*ltu*mmgl*mmsb12*mmt + 28*mmgl*mmsb2*mmt - 12*lgu*mmgl*mmsb2*
     mmt - 12*ltu*mmgl*mmsb2*mmt + 4*lgu*ltu*mmgl*mmsb2*mmt - 42*invdgb2*mmgl*
     mmsb22*mmt + 18*invdgb2*lgu*mmgl*mmsb22*mmt + 18*invdgb2*ltu*mmgl*mmsb22*
     mmt - 6*invdgb2*lgu*ltu*mmgl*mmsb22*mmt - 28*mmgl*mmst2*mmt + 12*lgu*mmgl*
     mmst2*mmt - 12*lt2u*mmgl*mmst2*mmt + 4*lgu*lt2u*mmgl*mmst2*mmt + 24*ltu*
     mmgl*mmst2*mmt - 8*lgu*ltu*mmgl*mmst2*mmt + 28*invdgb1*mmgl*mmsb1*mmst2*
     mmt - 12*invdgb1*lgu*mmgl*mmsb1*mmst2*mmt + 12*invdgb1*lt2u*mmgl*mmsb1*
     mmst2*mmt - 4*invdgb1*lgu*lt2u*mmgl*mmsb1*mmst2*mmt - 24*invdgb1*ltu*mmgl*
     mmsb1*mmst2*mmt + 8*invdgb1*lgu*ltu*mmgl*mmsb1*mmst2*mmt + 28*invdgb2*mmgl
     *mmsb2*mmst2*mmt - 12*invdgb2*lgu*mmgl*mmsb2*mmst2*mmt + 12*invdgb2*lt2u*
     mmgl*mmsb2*mmst2*mmt - 4*invdgb2*lgu*lt2u*mmgl*mmsb2*mmst2*mmt - 24*
     invdgb2*ltu*mmgl*mmsb2*mmst2*mmt + 8*invdgb2*lgu*ltu*mmgl*mmsb2*mmst2*mmt
     + 4*mmgl2*mmt*zt2 + 4*mmgl*mmsb1*mmt*zt2 - 6*invdgb1*mmgl*mmsb12*mmt*zt2 +
     4*mmgl*mmsb2*mmt*zt2 - 6*invdgb2*mmgl*mmsb22*mmt*zt2 - 4*mmgl*mmst2*mmt*
     zt2 + 4*invdgb1*mmgl*mmsb1*mmst2*mmt*zt2 + 4*invdgb2*mmgl*mmsb2*mmst2*mmt*
     zt2 + 2*mmgl2*mmt*pow2(lgu) + 2*mmgl*mmsb1*mmt*pow2(lgu) - 3*invdgb1*mmgl*
     mmsb12*mmt*pow2(lgu) + 2*mmgl*mmsb2*mmt*pow2(lgu) - 3*invdgb2*mmgl*mmsb22*
     mmt*pow2(lgu) - 2*mmgl*mmst2*mmt*pow2(lgu) + 2*invdgb1*mmgl*mmsb1*mmst2*
     mmt*pow2(lgu) + 2*invdgb2*mmgl*mmsb2*mmst2*mmt*pow2(lgu) + 2*mmgl*mmst2*
     mmt*pow2(lt2u) - 2*invdgb1*mmgl*mmsb1*mmst2*mmt*pow2(lt2u) - 2*invdgb2*
     mmgl*mmsb2*mmst2*mmt*pow2(lt2u) + 2*mmgl2*mmt*pow2(ltu) + 2*mmgl*mmsb1*mmt
     *pow2(ltu) - 3*invdgb1*mmgl*mmsb12*mmt*pow2(ltu) + 2*mmgl*mmsb2*mmt*pow2(
     ltu) - 3*invdgb2*mmgl*mmsb22*mmt*pow2(ltu) - 4*mmgl*mmst2*mmt*pow2(ltu) +
     4*invdgb1*mmgl*mmsb1*mmst2*mmt*pow2(ltu) + 4*invdgb2*mmgl*mmsb2*mmst2*mmt*
     pow2(ltu) - 6*invdgb1*lt2u*mmgl*mmt*pow2(mmst2) - 6*invdgb2*lt2u*mmgl*mmt*
     pow2(mmst2) + 2*invdgb1*lgu*lt2u*mmgl*mmt*pow2(mmst2) + 2*invdgb2*lgu*lt2u
     *mmgl*mmt*pow2(mmst2) + 6*invdgb1*ltu*mmgl*mmt*pow2(mmst2) + 6*invdgb2*ltu
     *mmgl*mmt*pow2(mmst2) - 2*invdgb1*lgu*ltu*mmgl*mmt*pow2(mmst2) - 2*invdgb2
     *lgu*ltu*mmgl*mmt*pow2(mmst2) + invdgb1*mmgl*mmt*pow2(lt2u)*pow2(mmst2) +
     invdgb2*mmgl*mmt*pow2(lt2u)*pow2(mmst2) - invdgb1*mmgl*mmt*pow2(ltu)*pow2(
     mmst2) - invdgb2*mmgl*mmt*pow2(ltu)*pow2(mmst2) - 28*mmgl*pow2(mmt) + 12*
     lgu*mmgl*pow2(mmt) + 12*ltu*mmgl*pow2(mmt) - 4*lgu*ltu*mmgl*pow2(mmt) + 28
     *invdgb1*mmgl*mmsb1*pow2(mmt) - 12*invdgb1*lgu*mmgl*mmsb1*pow2(mmt) - 12*
     invdgb1*ltu*mmgl*mmsb1*pow2(mmt) + 4*invdgb1*lgu*ltu*mmgl*mmsb1*pow2(mmt)
     + 28*invdgb2*mmgl*mmsb2*pow2(mmt) - 12*invdgb2*lgu*mmgl*mmsb2*pow2(mmt) -
     12*invdgb2*ltu*mmgl*mmsb2*pow2(mmt) + 4*invdgb2*lgu*ltu*mmgl*mmsb2*pow2(
     mmt) + 28*invdgb1*mmgl*mmst2*pow2(mmt) + 28*invdgb2*mmgl*mmst2*pow2(mmt) -
     6*invdgb1*lt2u*mmgl*mmst2*pow2(mmt) - 6*invdgb2*lt2u*mmgl*mmst2*pow2(mmt)
     - 2*invdgb1*lgu*lt2u*mmgl*mmst2*pow2(mmt) - 2*invdgb2*lgu*lt2u*mmgl*mmst2*
     pow2(mmt) - 18*invdgb1*ltu*mmgl*mmst2*pow2(mmt) - 18*invdgb2*ltu*mmgl*
     mmst2*pow2(mmt) + 2*invdgb1*lgu*ltu*mmgl*mmst2*pow2(mmt) + 2*invdgb2*lgu*
     ltu*mmgl*mmst2*pow2(mmt) + 4*invdgb1*lt2u*ltu*mmgl*mmst2*pow2(mmt) + 4*
     invdgb2*lt2u*ltu*mmgl*mmst2*pow2(mmt) - 4*mmgl*zt2*pow2(mmt) + 4*invdgb1*
     mmgl*mmsb1*zt2*pow2(mmt) + 4*invdgb2*mmgl*mmsb2*zt2*pow2(mmt) + 4*invdgb1*
     mmgl*mmst2*zt2*pow2(mmt) + 4*invdgb2*mmgl*mmst2*zt2*pow2(mmt) - 2*mmgl*
     pow2(lgu)*pow2(mmt) + 2*invdgb1*mmgl*mmsb1*pow2(lgu)*pow2(mmt) + 2*invdgb2
     *mmgl*mmsb2*pow2(lgu)*pow2(mmt) + invdgb1*mmgl*mmst2*pow2(lt2u)*pow2(mmt)
     + invdgb2*mmgl*mmst2*pow2(lt2u)*pow2(mmt) - 2*mmgl*pow2(ltu)*pow2(mmt) + 2
     *invdgb1*mmgl*mmsb1*pow2(ltu)*pow2(mmt) + 2*invdgb2*mmgl*mmsb2*pow2(ltu)*
     pow2(mmt) + 3*invdgb1*mmgl*mmst2*pow2(ltu)*pow2(mmt) + 3*invdgb2*mmgl*
     mmst2*pow2(ltu)*pow2(mmt) + mmgl*pow2(invdgb1)*(-14*mmsb12*mmst2*mmt + 6*
     lgu*mmsb12*mmst2*mmt - 6*lt2u*mmsb12*mmst2*mmt + 2*lgu*lt2u*mmsb12*mmst2*
     mmt + 12*ltu*mmsb12*mmst2*mmt - 4*lgu*ltu*mmsb12*mmst2*mmt - 2*mmsb12*
     mmst2*mmt*zt2 - 2*mmsb12*mmst2*mmt*pow2(ltu) + 6*lt2u*mmsb1*mmt*pow2(mmst2
     ) - 2*lgu*lt2u*mmsb1*mmt*pow2(mmst2) - 6*ltu*mmsb1*mmt*pow2(mmst2) + 2*lgu
     *ltu*mmsb1*mmt*pow2(mmst2) + mmsb1*mmt*pow2(ltu)*pow2(mmst2) - 14*mmsb12*
     pow2(mmt) + 6*lgu*mmsb12*pow2(mmt) + 6*ltu*mmsb12*pow2(mmt) - 2*lgu*ltu*
     mmsb12*pow2(mmt) - 28*mmsb1*mmst2*pow2(mmt) + 6*lt2u*mmsb1*mmst2*pow2(mmt)
     + 2*lgu*lt2u*mmsb1*mmst2*pow2(mmt) + 18*ltu*mmsb1*mmst2*pow2(mmt) - 2*lgu*
     ltu*mmsb1*mmst2*pow2(mmt) - 4*lt2u*ltu*mmsb1*mmst2*pow2(mmt) - 2*mmsb12*
     zt2*pow2(mmt) - 4*mmsb1*mmst2*zt2*pow2(mmt) - mmsb12*pow2(ltu)*pow2(mmt) -
     3*mmsb1*mmst2*pow2(ltu)*pow2(mmt) + pow2(lt2u)*(mmsb12*mmst2*mmt - mmsb1*
     mmt*pow2(mmst2) - mmsb1*mmst2*pow2(mmt)) + 14*mmt*pow3(mmsb1) - 6*lgu*mmt*
     pow3(mmsb1) - 6*ltu*mmt*pow3(mmsb1) + 2*lgu*ltu*mmt*pow3(mmsb1) + 2*mmt*
     zt2*pow3(mmsb1) + mmt*pow2(ltu)*pow3(mmsb1) - pow2(lgu)*(mmsb12*mmst2*mmt
     + mmsb12*pow2(mmt) - mmt*pow3(mmsb1))) + mmgl*pow2(invdgb2)*(-14*mmsb22*
     mmst2*mmt + 6*lgu*mmsb22*mmst2*mmt - 6*lt2u*mmsb22*mmst2*mmt + 2*lgu*lt2u*
     mmsb22*mmst2*mmt + 12*ltu*mmsb22*mmst2*mmt - 4*lgu*ltu*mmsb22*mmst2*mmt -
     2*mmsb22*mmst2*mmt*zt2 - 2*mmsb22*mmst2*mmt*pow2(ltu) + 6*lt2u*mmsb2*mmt*
     pow2(mmst2) - 2*lgu*lt2u*mmsb2*mmt*pow2(mmst2) - 6*ltu*mmsb2*mmt*pow2(
     mmst2) + 2*lgu*ltu*mmsb2*mmt*pow2(mmst2) + mmsb2*mmt*pow2(ltu)*pow2(mmst2)
     - 14*mmsb22*pow2(mmt) + 6*lgu*mmsb22*pow2(mmt) + 6*ltu*mmsb22*pow2(mmt) -
     2*lgu*ltu*mmsb22*pow2(mmt) - 28*mmsb2*mmst2*pow2(mmt) + 6*lt2u*mmsb2*mmst2
     *pow2(mmt) + 2*lgu*lt2u*mmsb2*mmst2*pow2(mmt) + 18*ltu*mmsb2*mmst2*pow2(
     mmt) - 2*lgu*ltu*mmsb2*mmst2*pow2(mmt) - 4*lt2u*ltu*mmsb2*mmst2*pow2(mmt)
     - 2*mmsb22*zt2*pow2(mmt) - 4*mmsb2*mmst2*zt2*pow2(mmt) - mmsb22*pow2(ltu)*
     pow2(mmt) - 3*mmsb2*mmst2*pow2(ltu)*pow2(mmt) + pow2(lt2u)*(mmsb22*mmst2*
     mmt - mmsb2*mmt*pow2(mmst2) - mmsb2*mmst2*pow2(mmt)) + 14*mmt*pow3(mmsb2)
     - 6*lgu*mmt*pow3(mmsb2) - 6*ltu*mmt*pow3(mmsb2) + 2*lgu*ltu*mmt*pow3(mmsb2
     ) + 2*mmt*zt2*pow3(mmsb2) + mmt*pow2(ltu)*pow3(mmsb2) - pow2(lgu)*(mmsb22*
     mmst2*mmt + mmsb22*pow2(mmt) - mmt*pow3(mmsb2)))))/(3.*mgl*mt)) + (4*s2b*(
     14*mmgl2*mmsb1 - 6*lgu*mmgl2*mmsb1 - 6*ltu*mmgl2*mmsb1 + 2*lgu*ltu*mmgl2*
     mmsb1 + 14*mmgl*mmsb12 - 6*lgu*mmgl*mmsb12 - 6*ltu*mmgl*mmsb12 + 2*lgu*ltu
     *mmgl*mmsb12 - 14*mmgl2*mmsb2 + 6*lgu*mmgl2*mmsb2 + 6*ltu*mmgl2*mmsb2 - 2*
     lgu*ltu*mmgl2*mmsb2 - 14*mmgl*mmsb22 + 6*lgu*mmgl*mmsb22 + 6*ltu*mmgl*
     mmsb22 - 2*lgu*ltu*mmgl*mmsb22 - 14*mmgl*mmsb1*mmst2 + 12*lgu*mmgl*mmsb1*
     mmst2 - 6*lt2u*mmgl*mmsb1*mmst2 + 6*ltu*mmgl*mmsb1*mmst2 - 4*lgu*ltu*mmgl*
     mmsb1*mmst2 + 2*lt2u*ltu*mmgl*mmsb1*mmst2 + 14*invdgb1*mmgl*mmsb12*mmst2 -
     12*invdgb1*lgu*mmgl*mmsb12*mmst2 + 6*invdgb1*lt2u*mmgl*mmsb12*mmst2 - 6*
     invdgb1*ltu*mmgl*mmsb12*mmst2 + 4*invdgb1*lgu*ltu*mmgl*mmsb12*mmst2 - 2*
     invdgb1*lt2u*ltu*mmgl*mmsb12*mmst2 + 14*mmgl*mmsb2*mmst2 - 12*lgu*mmgl*
     mmsb2*mmst2 + 6*lt2u*mmgl*mmsb2*mmst2 - 6*ltu*mmgl*mmsb2*mmst2 + 4*lgu*ltu
     *mmgl*mmsb2*mmst2 - 2*lt2u*ltu*mmgl*mmsb2*mmst2 - 14*invdgb2*mmgl*mmsb22*
     mmst2 + 12*invdgb2*lgu*mmgl*mmsb22*mmst2 - 6*invdgb2*lt2u*mmgl*mmsb22*
     mmst2 + 6*invdgb2*ltu*mmgl*mmsb22*mmst2 - 4*invdgb2*lgu*ltu*mmgl*mmsb22*
     mmst2 + 2*invdgb2*lt2u*ltu*mmgl*mmsb22*mmst2 - 14*mmgl*mmsb1*mmt + 6*lgu*
     mmgl*mmsb1*mmt + 6*ltu*mmgl*mmsb1*mmt - 2*lgu*ltu*mmgl*mmsb1*mmt + 14*
     invdgb1*mmgl*mmsb12*mmt - 6*invdgb1*lgu*mmgl*mmsb12*mmt - 6*invdgb1*ltu*
     mmgl*mmsb12*mmt + 2*invdgb1*lgu*ltu*mmgl*mmsb12*mmt + 14*mmgl*mmsb2*mmt -
     6*lgu*mmgl*mmsb2*mmt - 6*ltu*mmgl*mmsb2*mmt + 2*lgu*ltu*mmgl*mmsb2*mmt -
     14*invdgb2*mmgl*mmsb22*mmt + 6*invdgb2*lgu*mmgl*mmsb22*mmt + 6*invdgb2*ltu
     *mmgl*mmsb22*mmt - 2*invdgb2*lgu*ltu*mmgl*mmsb22*mmt + 56*invdgb1*mmgl*
     mmsb1*mmst2*mmt - 6*invdgb1*lgu*mmgl*mmsb1*mmst2*mmt - 6*invdgb1*lt2u*mmgl
     *mmsb1*mmst2*mmt - 4*invdgb1*lgu*lt2u*mmgl*mmsb1*mmst2*mmt - 36*invdgb1*
     ltu*mmgl*mmsb1*mmst2*mmt + 6*invdgb1*lgu*ltu*mmgl*mmsb1*mmst2*mmt + 6*
     invdgb1*lt2u*ltu*mmgl*mmsb1*mmst2*mmt - 56*invdgb2*mmgl*mmsb2*mmst2*mmt +
     6*invdgb2*lgu*mmgl*mmsb2*mmst2*mmt + 6*invdgb2*lt2u*mmgl*mmsb2*mmst2*mmt +
     4*invdgb2*lgu*lt2u*mmgl*mmsb2*mmst2*mmt + 36*invdgb2*ltu*mmgl*mmsb2*mmst2*
     mmt - 6*invdgb2*lgu*ltu*mmgl*mmsb2*mmst2*mmt - 6*invdgb2*lt2u*ltu*mmgl*
     mmsb2*mmst2*mmt + 2*mmgl2*mmsb1*zt2 + 2*mmgl*mmsb12*zt2 - 2*mmgl2*mmsb2*
     zt2 - 2*mmgl*mmsb22*zt2 - 2*mmgl*mmsb1*mmst2*zt2 + 2*invdgb1*mmgl*mmsb12*
     mmst2*zt2 + 2*mmgl*mmsb2*mmst2*zt2 - 2*invdgb2*mmgl*mmsb22*mmst2*zt2 - 2*
     mmgl*mmsb1*mmt*zt2 + 2*invdgb1*mmgl*mmsb12*mmt*zt2 + 2*mmgl*mmsb2*mmt*zt2
     - 2*invdgb2*mmgl*mmsb22*mmt*zt2 + 8*invdgb1*mmgl*mmsb1*mmst2*mmt*zt2 - 8*
     invdgb2*mmgl*mmsb2*mmst2*mmt*zt2 + mmgl2*mmsb1*pow2(ltu) + mmgl*mmsb12*
     pow2(ltu) - mmgl2*mmsb2*pow2(ltu) - mmgl*mmsb22*pow2(ltu) - mmgl*mmsb1*
     mmst2*pow2(ltu) + invdgb1*mmgl*mmsb12*mmst2*pow2(ltu) + mmgl*mmsb2*mmst2*
     pow2(ltu) - invdgb2*mmgl*mmsb22*mmst2*pow2(ltu) - mmgl*mmsb1*mmt*pow2(ltu)
     + invdgb1*mmgl*mmsb12*mmt*pow2(ltu) + mmgl*mmsb2*mmt*pow2(ltu) - invdgb2*
     mmgl*mmsb22*mmt*pow2(ltu) + 6*invdgb1*mmgl*mmsb1*mmst2*mmt*pow2(ltu) - 6*
     invdgb2*mmgl*mmsb2*mmst2*mmt*pow2(ltu) + 14*invdgb1*mmgl*mmsb1*pow2(mmst2)
     + 6*invdgb1*lgu*mmgl*mmsb1*pow2(mmst2) - 12*invdgb1*lt2u*mmgl*mmsb1*pow2(
     mmst2) - 6*invdgb1*ltu*mmgl*mmsb1*pow2(mmst2) - 2*invdgb1*lgu*ltu*mmgl*
     mmsb1*pow2(mmst2) + 4*invdgb1*lt2u*ltu*mmgl*mmsb1*pow2(mmst2) - 14*invdgb2
     *mmgl*mmsb2*pow2(mmst2) - 6*invdgb2*lgu*mmgl*mmsb2*pow2(mmst2) + 12*
     invdgb2*lt2u*mmgl*mmsb2*pow2(mmst2) + 6*invdgb2*ltu*mmgl*mmsb2*pow2(mmst2)
     + 2*invdgb2*lgu*ltu*mmgl*mmsb2*pow2(mmst2) - 4*invdgb2*lt2u*ltu*mmgl*mmsb2
     *pow2(mmst2) + 14*invdgb1*mmgl*mmt*pow2(mmst2) - 14*invdgb2*mmgl*mmt*pow2(
     mmst2) - 6*invdgb1*lt2u*mmgl*mmt*pow2(mmst2) + 6*invdgb2*lt2u*mmgl*mmt*
     pow2(mmst2) - 6*invdgb1*ltu*mmgl*mmt*pow2(mmst2) + 6*invdgb2*ltu*mmgl*mmt*
     pow2(mmst2) + 2*invdgb1*lt2u*ltu*mmgl*mmt*pow2(mmst2) - 2*invdgb2*lt2u*ltu
     *mmgl*mmt*pow2(mmst2) + 2*invdgb1*mmgl*mmsb1*zt2*pow2(mmst2) - 2*invdgb2*
     mmgl*mmsb2*zt2*pow2(mmst2) + 2*invdgb1*mmgl*mmt*zt2*pow2(mmst2) - 2*
     invdgb2*mmgl*mmt*zt2*pow2(mmst2) + invdgb1*mmgl*mmsb1*pow2(ltu)*pow2(mmst2
     ) - invdgb2*mmgl*mmsb2*pow2(ltu)*pow2(mmst2) + invdgb1*mmgl*mmt*pow2(ltu)*
     pow2(mmst2) - invdgb2*mmgl*mmt*pow2(ltu)*pow2(mmst2) - 14*invdgb1*mmgl*
     pow3(mmsb1) + 6*invdgb1*lgu*mmgl*pow3(mmsb1) + 6*invdgb1*ltu*mmgl*pow3(
     mmsb1) - 2*invdgb1*lgu*ltu*mmgl*pow3(mmsb1) - 2*invdgb1*mmgl*zt2*pow3(
     mmsb1) - invdgb1*mmgl*pow2(ltu)*pow3(mmsb1) + 14*invdgb2*mmgl*pow3(mmsb2)
     - 6*invdgb2*lgu*mmgl*pow3(mmsb2) - 6*invdgb2*ltu*mmgl*pow3(mmsb2) + 2*
     invdgb2*lgu*ltu*mmgl*pow3(mmsb2) + 2*invdgb2*mmgl*zt2*pow3(mmsb2) +
     invdgb2*mmgl*pow2(ltu)*pow3(mmsb2) + pow2(lgu)*(mmgl2*mmsb1 + mmgl*mmsb12
     - mmgl2*mmsb2 - mmgl*mmsb22 - 2*mmgl*mmsb1*mmst2 + 2*invdgb1*mmgl*mmsb12*
     mmst2 + 2*mmgl*mmsb2*mmst2 - 2*invdgb2*mmgl*mmsb22*mmst2 - mmgl*mmsb1*mmt
     + invdgb1*mmgl*mmsb12*mmt + mmgl*mmsb2*mmt - invdgb2*mmgl*mmsb22*mmt +
     invdgb1*mmgl*mmsb1*mmst2*mmt - invdgb2*mmgl*mmsb2*mmst2*mmt + (-(invdgb1*
     mmgl*mmsb1) + invdgb2*mmgl*mmsb2)*pow2(mmst2) - invdgb1*mmgl*pow3(mmsb1) +
     invdgb2*mmgl*pow3(mmsb2)) - 14*invdgb1*mmgl*pow3(mmst2) + 14*invdgb2*mmgl*
     pow3(mmst2) + 6*invdgb1*lt2u*mmgl*pow3(mmst2) - 6*invdgb2*lt2u*mmgl*pow3(
     mmst2) + 6*invdgb1*ltu*mmgl*pow3(mmst2) - 6*invdgb2*ltu*mmgl*pow3(mmst2) -
     2*invdgb1*lt2u*ltu*mmgl*pow3(mmst2) + 2*invdgb2*lt2u*ltu*mmgl*pow3(mmst2)
     - 2*invdgb1*mmgl*zt2*pow3(mmst2) + 2*invdgb2*mmgl*zt2*pow3(mmst2) -
     invdgb1*mmgl*pow2(ltu)*pow3(mmst2) + invdgb2*mmgl*pow2(ltu)*pow3(mmst2) +
     mmgl*pow2(lt2u)*(mmst2*(mmsb1 - invdgb1*mmsb12 - mmsb2 + invdgb2*mmsb22 +
     invdgb1*mmsb1*mmt - invdgb2*mmsb2*mmt) + (invdgb1*(2*mmsb1 + mmt) -
     invdgb2*(2*mmsb2 + mmt))*pow2(mmst2) + (-invdgb1 + invdgb2)*pow3(mmst2))))
     /(3.*mb*mgl) - (2*(28*mmgl2 - 12*lgu*mmgl2 - 12*ltu*mmgl2 + 4*lgu*ltu*
     mmgl2 + 28*mmgl*mmsb1 - 12*lgu*mmgl*mmsb1 - 12*ltu*mmgl*mmsb1 + 4*lgu*ltu*
     mmgl*mmsb1 + 42*mmsb12 - 18*lgu*mmsb12 - 18*ltu*mmsb12 + 6*lgu*ltu*mmsb12
     + 28*mmgl*mmsb2 - 12*lgu*mmgl*mmsb2 - 12*ltu*mmgl*mmsb2 + 4*lgu*ltu*mmgl*
     mmsb2 + 42*mmsb22 - 18*lgu*mmsb22 - 18*ltu*mmsb22 + 6*lgu*ltu*mmsb22 - 28*
     mmgl*mmst2 + 24*lgu*mmgl*mmst2 - 12*lt2u*mmgl*mmst2 + 12*ltu*mmgl*mmst2 -
     8*lgu*ltu*mmgl*mmst2 + 4*lt2u*ltu*mmgl*mmst2 - 28*mmsb1*mmst2 + 24*lgu*
     mmsb1*mmst2 - 12*lt2u*mmsb1*mmst2 + 12*ltu*mmsb1*mmst2 - 8*lgu*ltu*mmsb1*
     mmst2 + 4*lt2u*ltu*mmsb1*mmst2 + 42*invdgb1*mmsb12*mmst2 - 36*invdgb1*lgu*
     mmsb12*mmst2 + 18*invdgb1*lt2u*mmsb12*mmst2 - 18*invdgb1*ltu*mmsb12*mmst2
     + 12*invdgb1*lgu*ltu*mmsb12*mmst2 - 6*invdgb1*lt2u*ltu*mmsb12*mmst2 - 28*
     mmsb2*mmst2 + 24*lgu*mmsb2*mmst2 - 12*lt2u*mmsb2*mmst2 + 12*ltu*mmsb2*
     mmst2 - 8*lgu*ltu*mmsb2*mmst2 + 4*lt2u*ltu*mmsb2*mmst2 + 42*invdgb2*mmsb22
     *mmst2 - 36*invdgb2*lgu*mmsb22*mmst2 + 18*invdgb2*lt2u*mmsb22*mmst2 - 18*
     invdgb2*ltu*mmsb22*mmst2 + 12*invdgb2*lgu*ltu*mmsb22*mmst2 - 6*invdgb2*
     lt2u*ltu*mmsb22*mmst2 - 28*mmgl*mmt + 12*lgu*mmgl*mmt + 12*ltu*mmgl*mmt -
     4*lgu*ltu*mmgl*mmt - 28*mmsb1*mmt + 12*lgu*mmsb1*mmt + 12*ltu*mmsb1*mmt -
     4*lgu*ltu*mmsb1*mmt + 42*invdgb1*mmsb12*mmt - 18*invdgb1*lgu*mmsb12*mmt -
     18*invdgb1*ltu*mmsb12*mmt + 6*invdgb1*lgu*ltu*mmsb12*mmt - 28*mmsb2*mmt +
     12*lgu*mmsb2*mmt + 12*ltu*mmsb2*mmt - 4*lgu*ltu*mmsb2*mmt + 42*invdgb2*
     mmsb22*mmt - 18*invdgb2*lgu*mmsb22*mmt - 18*invdgb2*ltu*mmsb22*mmt + 6*
     invdgb2*lgu*ltu*mmsb22*mmt - 112*mmst2*mmt + 12*lgu*mmst2*mmt + 12*lt2u*
     mmst2*mmt + 8*lgu*lt2u*mmst2*mmt + 72*ltu*mmst2*mmt - 12*lgu*ltu*mmst2*mmt
      - 12*lt2u*ltu*mmst2*mmt + 112*invdgb1*mmsb1*mmst2*mmt - 12*invdgb1*lgu*
     mmsb1*mmst2*mmt - 12*invdgb1*lt2u*mmsb1*mmst2*mmt - 8*invdgb1*lgu*lt2u*
     mmsb1*mmst2*mmt - 72*invdgb1*ltu*mmsb1*mmst2*mmt + 12*invdgb1*lgu*ltu*
     mmsb1*mmst2*mmt + 12*invdgb1*lt2u*ltu*mmsb1*mmst2*mmt + 112*invdgb2*mmsb2*
     mmst2*mmt - 12*invdgb2*lgu*mmsb2*mmst2*mmt - 12*invdgb2*lt2u*mmsb2*mmst2*
     mmt - 8*invdgb2*lgu*lt2u*mmsb2*mmst2*mmt - 72*invdgb2*ltu*mmsb2*mmst2*mmt
     + 12*invdgb2*lgu*ltu*mmsb2*mmst2*mmt + 12*invdgb2*lt2u*ltu*mmsb2*mmst2*mmt
      + 4*mmgl2*zt2 + 4*mmgl*mmsb1*zt2 + 6*mmsb12*zt2 + 4*mmgl*mmsb2*zt2 + 6*
     mmsb22*zt2 - 4*mmgl*mmst2*zt2 - 4*mmsb1*mmst2*zt2 + 6*invdgb1*mmsb12*mmst2
     *zt2 - 4*mmsb2*mmst2*zt2 + 6*invdgb2*mmsb22*mmst2*zt2 - 4*mmgl*mmt*zt2 - 4
     *mmsb1*mmt*zt2 + 6*invdgb1*mmsb12*mmt*zt2 - 4*mmsb2*mmt*zt2 + 6*invdgb2*
     mmsb22*mmt*zt2 - 16*mmst2*mmt*zt2 + 16*invdgb1*mmsb1*mmst2*mmt*zt2 + 16*
     invdgb2*mmsb2*mmst2*mmt*zt2 + 2*mmgl2*pow2(lgu) + 2*mmgl*mmsb1*pow2(lgu) +
     3*mmsb12*pow2(lgu) + 2*mmgl*mmsb2*pow2(lgu) + 3*mmsb22*pow2(lgu) - 4*mmgl*
     mmst2*pow2(lgu) - 4*mmsb1*mmst2*pow2(lgu) + 6*invdgb1*mmsb12*mmst2*pow2(
     lgu) - 4*mmsb2*mmst2*pow2(lgu) + 6*invdgb2*mmsb22*mmst2*pow2(lgu) - 2*mmgl
     *mmt*pow2(lgu) - 2*mmsb1*mmt*pow2(lgu) + 3*invdgb1*mmsb12*mmt*pow2(lgu) -
     2*mmsb2*mmt*pow2(lgu) + 3*invdgb2*mmsb22*mmt*pow2(lgu) - 2*mmst2*mmt*pow2(
     lgu) + 2*invdgb1*mmsb1*mmst2*mmt*pow2(lgu) + 2*invdgb2*mmsb2*mmst2*mmt*
     pow2(lgu) + 2*mmgl*mmst2*pow2(lt2u) + 2*mmsb1*mmst2*pow2(lt2u) - 3*invdgb1
     *mmsb12*mmst2*pow2(lt2u) + 2*mmsb2*mmst2*pow2(lt2u) - 3*invdgb2*mmsb22*
     mmst2*pow2(lt2u) - 2*mmst2*mmt*pow2(lt2u) + 2*invdgb1*mmsb1*mmst2*mmt*pow2
     (lt2u) + 2*invdgb2*mmsb2*mmst2*mmt*pow2(lt2u) + 2*mmgl2*pow2(ltu) + 2*mmgl
     *mmsb1*pow2(ltu) + 3*mmsb12*pow2(ltu) + 2*mmgl*mmsb2*pow2(ltu) + 3*mmsb22*
     pow2(ltu) - 2*mmgl*mmst2*pow2(ltu) - 2*mmsb1*mmst2*pow2(ltu) + 3*invdgb1*
     mmsb12*mmst2*pow2(ltu) - 2*mmsb2*mmst2*pow2(ltu) + 3*invdgb2*mmsb22*mmst2*
     pow2(ltu) - 2*mmgl*mmt*pow2(ltu) - 2*mmsb1*mmt*pow2(ltu) + 3*invdgb1*
     mmsb12*mmt*pow2(ltu) - 2*mmsb2*mmt*pow2(ltu) + 3*invdgb2*mmsb22*mmt*pow2(
     ltu) - 12*mmst2*mmt*pow2(ltu) + 12*invdgb1*mmsb1*mmst2*mmt*pow2(ltu) + 12*
     invdgb2*mmsb2*mmst2*mmt*pow2(ltu) - 28*pow2(mmst2) - 12*lgu*pow2(mmst2) +
     24*lt2u*pow2(mmst2) + 12*ltu*pow2(mmst2) + 4*lgu*ltu*pow2(mmst2) - 8*lt2u*
     ltu*pow2(mmst2) + 28*invdgb1*mmsb1*pow2(mmst2) + 12*invdgb1*lgu*mmsb1*pow2
     (mmst2) - 24*invdgb1*lt2u*mmsb1*pow2(mmst2) - 12*invdgb1*ltu*mmsb1*pow2(
     mmst2) - 4*invdgb1*lgu*ltu*mmsb1*pow2(mmst2) + 8*invdgb1*lt2u*ltu*mmsb1*
     pow2(mmst2) + 28*invdgb2*mmsb2*pow2(mmst2) + 12*invdgb2*lgu*mmsb2*pow2(
     mmst2) - 24*invdgb2*lt2u*mmsb2*pow2(mmst2) - 12*invdgb2*ltu*mmsb2*pow2(
     mmst2) - 4*invdgb2*lgu*ltu*mmsb2*pow2(mmst2) + 8*invdgb2*lt2u*ltu*mmsb2*
     pow2(mmst2) + 14*invdgb1*mmt*pow2(mmst2) + 14*invdgb2*mmt*pow2(mmst2) - 6*
     invdgb1*lt2u*mmt*pow2(mmst2) - 6*invdgb2*lt2u*mmt*pow2(mmst2) - 6*invdgb1*
     ltu*mmt*pow2(mmst2) - 6*invdgb2*ltu*mmt*pow2(mmst2) + 2*invdgb1*lt2u*ltu*
     mmt*pow2(mmst2) + 2*invdgb2*lt2u*ltu*mmt*pow2(mmst2) - 4*zt2*pow2(mmst2) +
     4*invdgb1*mmsb1*zt2*pow2(mmst2) + 4*invdgb2*mmsb2*zt2*pow2(mmst2) + 2*
     invdgb1*mmt*zt2*pow2(mmst2) + 2*invdgb2*mmt*zt2*pow2(mmst2) + 2*pow2(lgu)*
     pow2(mmst2) - 2*invdgb1*mmsb1*pow2(lgu)*pow2(mmst2) - 2*invdgb2*mmsb2*pow2
     (lgu)*pow2(mmst2) - 4*pow2(lt2u)*pow2(mmst2) + 4*invdgb1*mmsb1*pow2(lt2u)*
     pow2(mmst2) + 4*invdgb2*mmsb2*pow2(lt2u)*pow2(mmst2) + invdgb1*mmt*pow2(
     lt2u)*pow2(mmst2) + invdgb2*mmt*pow2(lt2u)*pow2(mmst2) - 2*pow2(ltu)*pow2(
     mmst2) + 2*invdgb1*mmsb1*pow2(ltu)*pow2(mmst2) + 2*invdgb2*mmsb2*pow2(ltu)
     *pow2(mmst2) + invdgb1*mmt*pow2(ltu)*pow2(mmst2) + invdgb2*mmt*pow2(ltu)*
     pow2(mmst2) - 56*invdgb1*pow3(mmsb1) + 24*invdgb1*lgu*pow3(mmsb1) + 24*
     invdgb1*ltu*pow3(mmsb1) - 8*invdgb1*lgu*ltu*pow3(mmsb1) - 8*invdgb1*zt2*
     pow3(mmsb1) - 4*invdgb1*pow2(lgu)*pow3(mmsb1) - 4*invdgb1*pow2(ltu)*pow3(
     mmsb1) - 56*invdgb2*pow3(mmsb2) + 24*invdgb2*lgu*pow3(mmsb2) + 24*invdgb2*
     ltu*pow3(mmsb2) - 8*invdgb2*lgu*ltu*pow3(mmsb2) - 8*invdgb2*zt2*pow3(mmsb2
     ) - 4*invdgb2*pow2(lgu)*pow3(mmsb2) - 4*invdgb2*pow2(ltu)*pow3(mmsb2) - 14
     *invdgb1*pow3(mmst2) - 14*invdgb2*pow3(mmst2) + 6*invdgb1*lt2u*pow3(mmst2)
     + 6*invdgb2*lt2u*pow3(mmst2) + 6*invdgb1*ltu*pow3(mmst2) + 6*invdgb2*ltu*
     pow3(mmst2) - 2*invdgb1*lt2u*ltu*pow3(mmst2) - 2*invdgb2*lt2u*ltu*pow3(
     mmst2) - 2*invdgb1*zt2*pow3(mmst2) - 2*invdgb2*zt2*pow3(mmst2) - invdgb1*
     pow2(lt2u)*pow3(mmst2) - invdgb2*pow2(lt2u)*pow3(mmst2) - invdgb1*pow2(ltu
     )*pow3(mmst2) - invdgb2*pow2(ltu)*pow3(mmst2) + pow2(invdgb1)*(-56*mmsb12*
     mmst2*mmt + 6*lgu*mmsb12*mmst2*mmt + 6*lt2u*mmsb12*mmst2*mmt + 4*lgu*lt2u*
     mmsb12*mmst2*mmt + 36*ltu*mmsb12*mmst2*mmt - 6*lgu*ltu*mmsb12*mmst2*mmt -
     6*lt2u*ltu*mmsb12*mmst2*mmt - 8*mmsb12*mmst2*mmt*zt2 - 6*mmsb12*mmst2*mmt*
     pow2(ltu) - 14*mmsb12*pow2(mmst2) - 6*lgu*mmsb12*pow2(mmst2) + 12*lt2u*
     mmsb12*pow2(mmst2) + 6*ltu*mmsb12*pow2(mmst2) + 2*lgu*ltu*mmsb12*pow2(
     mmst2) - 4*lt2u*ltu*mmsb12*pow2(mmst2) - 14*mmsb1*mmt*pow2(mmst2) + 6*lt2u
     *mmsb1*mmt*pow2(mmst2) + 6*ltu*mmsb1*mmt*pow2(mmst2) - 2*lt2u*ltu*mmsb1*
     mmt*pow2(mmst2) - 2*mmsb12*zt2*pow2(mmst2) - 2*mmsb1*mmt*zt2*pow2(mmst2) -
     mmsb12*pow2(ltu)*pow2(mmst2) - mmsb1*mmt*pow2(ltu)*pow2(mmst2) - 14*mmst2*
     pow3(mmsb1) + 12*lgu*mmst2*pow3(mmsb1) - 6*lt2u*mmst2*pow3(mmsb1) + 6*ltu*
     mmst2*pow3(mmsb1) - 4*lgu*ltu*mmst2*pow3(mmsb1) + 2*lt2u*ltu*mmst2*pow3(
     mmsb1) - 14*mmt*pow3(mmsb1) + 6*lgu*mmt*pow3(mmsb1) + 6*ltu*mmt*pow3(mmsb1
     ) - 2*lgu*ltu*mmt*pow3(mmsb1) - 2*mmst2*zt2*pow3(mmsb1) - 2*mmt*zt2*pow3(
     mmsb1) - mmst2*pow2(ltu)*pow3(mmsb1) - mmt*pow2(ltu)*pow3(mmsb1) + 14*
     mmsb1*pow3(mmst2) - 6*lt2u*mmsb1*pow3(mmst2) - 6*ltu*mmsb1*pow3(mmst2) + 2
     *lt2u*ltu*mmsb1*pow3(mmst2) + 2*mmsb1*zt2*pow3(mmst2) + mmsb1*pow2(ltu)*
     pow3(mmst2) + pow2(lt2u)*(-(mmsb12*mmst2*mmt) - (2*mmsb12 + mmsb1*mmt)*
     pow2(mmst2) + mmst2*pow3(mmsb1) + mmsb1*pow3(mmst2)) + 14*pow4(mmsb1) - 6*
     lgu*pow4(mmsb1) - 6*ltu*pow4(mmsb1) + 2*lgu*ltu*pow4(mmsb1) + 2*zt2*pow4(
     mmsb1) + pow2(ltu)*pow4(mmsb1) + pow2(lgu)*(-(mmsb12*mmst2*mmt) + mmsb12*
     pow2(mmst2) - (2*mmst2 + mmt)*pow3(mmsb1) + pow4(mmsb1))) + pow2(invdgb2)*
     (-56*mmsb22*mmst2*mmt + 6*lgu*mmsb22*mmst2*mmt + 6*lt2u*mmsb22*mmst2*mmt +
     4*lgu*lt2u*mmsb22*mmst2*mmt + 36*ltu*mmsb22*mmst2*mmt - 6*lgu*ltu*mmsb22*
     mmst2*mmt - 6*lt2u*ltu*mmsb22*mmst2*mmt - 8*mmsb22*mmst2*mmt*zt2 - 6*
     mmsb22*mmst2*mmt*pow2(ltu) - 14*mmsb22*pow2(mmst2) - 6*lgu*mmsb22*pow2(
     mmst2) + 12*lt2u*mmsb22*pow2(mmst2) + 6*ltu*mmsb22*pow2(mmst2) + 2*lgu*ltu
     *mmsb22*pow2(mmst2) - 4*lt2u*ltu*mmsb22*pow2(mmst2) - 14*mmsb2*mmt*pow2(
     mmst2) + 6*lt2u*mmsb2*mmt*pow2(mmst2) + 6*ltu*mmsb2*mmt*pow2(mmst2) - 2*
     lt2u*ltu*mmsb2*mmt*pow2(mmst2) - 2*mmsb22*zt2*pow2(mmst2) - 2*mmsb2*mmt*
     zt2*pow2(mmst2) - mmsb22*pow2(ltu)*pow2(mmst2) - mmsb2*mmt*pow2(ltu)*pow2(
     mmst2) - 14*mmst2*pow3(mmsb2) + 12*lgu*mmst2*pow3(mmsb2) - 6*lt2u*mmst2*
     pow3(mmsb2) + 6*ltu*mmst2*pow3(mmsb2) - 4*lgu*ltu*mmst2*pow3(mmsb2) + 2*
     lt2u*ltu*mmst2*pow3(mmsb2) - 14*mmt*pow3(mmsb2) + 6*lgu*mmt*pow3(mmsb2) +
     6*ltu*mmt*pow3(mmsb2) - 2*lgu*ltu*mmt*pow3(mmsb2) - 2*mmst2*zt2*pow3(mmsb2
     ) - 2*mmt*zt2*pow3(mmsb2) - mmst2*pow2(ltu)*pow3(mmsb2) - mmt*pow2(ltu)*
     pow3(mmsb2) + 14*mmsb2*pow3(mmst2) - 6*lt2u*mmsb2*pow3(mmst2) - 6*ltu*
     mmsb2*pow3(mmst2) + 2*lt2u*ltu*mmsb2*pow3(mmst2) + 2*mmsb2*zt2*pow3(mmst2)
     + mmsb2*pow2(ltu)*pow3(mmst2) + pow2(lt2u)*(-(mmsb22*mmst2*mmt) - (2*
     mmsb22 + mmsb2*mmt)*pow2(mmst2) + mmst2*pow3(mmsb2) + mmsb2*pow3(mmst2)) +
     14*pow4(mmsb2) - 6*lgu*pow4(mmsb2) - 6*ltu*pow4(mmsb2) + 2*lgu*ltu*pow4(
     mmsb2) + 2*zt2*pow4(mmsb2) + pow2(ltu)*pow4(mmsb2) + pow2(lgu)*(-(mmsb22*
     mmst2*mmt) + mmsb22*pow2(mmst2) - (2*mmst2 + mmt)*pow3(mmsb2) + pow4(mmsb2
     )))))/3.);

   return pow4(g3) * result * twoLoop;
}

} // namespace mssm_twoloop_mb
} // namespace flexiblesusy
