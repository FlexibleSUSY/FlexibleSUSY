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

// This file has been generated at Sat 18 Apr 2020 01:11:28
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
   const double g3     = pars.g3;
   const double Xt     = pars.xt;
   const double Xb     = pars.xb;
   const double mgl    = pars.mg;
   const double mt     = pars.mt;
   const double mb     = pars.mb;
   const double mmt    = pow2(pars.mt);
   const double mmb    = pow2(pars.mb);
   const double mmgl   = pow2(pars.mg);
   const double mmgl2  = pow2(mmgl);
   const double mmgl3  = mmgl*mmgl2;
   const double mmu    = pow2(pars.Q);
   const double mmst1  = pow2(pars.mst1);
   const double mmst2  = pow2(pars.mst2);
   const double mmsb1  = pow2(pars.msb1);
   const double mmsb2  = pow2(pars.msb2);
   const double mmsb12 = pow2(mmsb1);
   const double mmsb22 = pow2(mmsb2);
   const double mmsusy = pow2(pars.msusy);
   const double lgu    = std::log(mmgl/mmu);
   const double lt1u   = std::log(mmst1/mmu);
   const double lt2u   = std::log(mmst2/mmu);
   const double lb1u   = std::log(mmsb1/mmu);
   const double lb2u   = std::log(mmsb2/mmu);
   const double lsu    = std::log(mmsusy/mmu);
   const double ltu    = std::log(mmt/mmu);
   const double s2t    = 2*mt*Xt / (mmst1 - mmst2);
   const double s2b    = 2*mb*Xb / (mmsb1 - mmsb2);
   const double invdgb1 = 1/(-mmgl + mmsb1);
   const double invdgb2 = 1/(-mmgl + mmsb2);
   const double invdb12 = 1/(mmsb1 - mmsb2);

   const double result =
   16.85185185185185 + (11*lb1u)/9. + (11*lb2u)/9. + 8*lgu + (11*lt1u)/9. + (11
     *lt2u)/9. - (124*ltu)/9. + (16*lgu*ltu)/3. - (16*invdgb1*mmsb1)/3. + 6*
     invdgb2*mmsb1 - 8*invdgb1*lb1u*mmsb1 - 4*invdgb2*lb1u*mmsb1 + 6*invdgb1*
     lgu*mmsb1 - 2*invdgb2*lgu*mmsb1 + 20*invdgb1*ltu*mmsb1 - (4*invdgb1*lb1u*
     ltu*mmsb1)/3. - (16*invdgb1*lgu*ltu*mmsb1)/3. + 6*invdgb1*mmsb2 - (16*
     invdgb2*mmsb2)/3. - 4*invdgb1*lb2u*mmsb2 - 8*invdgb2*lb2u*mmsb2 - 2*
     invdgb1*lgu*mmsb2 + 6*invdgb2*lgu*mmsb2 + 20*invdgb2*ltu*mmsb2 - (4*
     invdgb2*lb2u*ltu*mmsb2)/3. - (16*invdgb2*lgu*ltu*mmsb2)/3. - 8*invdgb1*
     mmst1 - 8*invdgb2*mmst1 - 2*invdgb1*lgu*mmst1 - 2*invdgb2*lgu*mmst1 + 2*
     invdgb1*lt1u*mmst1 + 2*invdgb2*lt1u*mmst1 + invdgb1*lgu*lt1u*mmst1 +
     invdgb2*lgu*lt1u*mmst1 + 6*invdgb1*ltu*mmst1 + 6*invdgb2*ltu*mmst1 - 2*
     invdgb1*lt1u*ltu*mmst1 - 2*invdgb2*lt1u*ltu*mmst1 - 8*invdgb1*mmst2 - 8*
     invdgb2*mmst2 - 2*invdgb1*lgu*mmst2 - 2*invdgb2*lgu*mmst2 + 2*invdgb1*lt2u
     *mmst2 + 2*invdgb2*lt2u*mmst2 + invdgb1*lgu*lt2u*mmst2 + invdgb2*lgu*lt2u*
     mmst2 + 6*invdgb1*ltu*mmst2 + 6*invdgb2*ltu*mmst2 - 2*invdgb1*lt2u*ltu*
     mmst2 - 2*invdgb2*lt2u*ltu*mmst2 + 48*invdgb1*mmsusy + 48*invdgb2*mmsusy -
     16*invdgb1*lgu*mmsusy - 16*invdgb2*lgu*mmsusy - 12*invdgb1*mmt - 12*
     invdgb2*mmt + 4*invdgb1*lgu*mmt + 4*invdgb2*lgu*mmt + 8*invdgb1*ltu*mmt +
     8*invdgb2*ltu*mmt - 2*invdgb1*lgu*ltu*mmt - 2*invdgb2*lgu*ltu*mmt + (4*zt2
     )/3. - (19*invdgb1*mmsb1*zt2)/3. + invdgb2*mmsb1*zt2 + invdgb1*mmsb2*zt2 -
     (19*invdgb2*mmsb2*zt2)/3. - invdgb1*mmst1*zt2 - invdgb2*mmst1*zt2 -
     invdgb1*mmst2*zt2 - invdgb2*mmst2*zt2 + 8*invdgb1*mmsusy*zt2 + 8*invdgb2*
     mmsusy*zt2 - 2*invdgb1*mmt*zt2 - 2*invdgb2*mmt*zt2 + (20*mmsb12*pow2(
     invdgb1))/3. + (70*lb1u*mmsb12*pow2(invdgb1))/3. - (70*lgu*mmsb12*pow2(
     invdgb1))/3. - 16*ltu*mmsb12*pow2(invdgb1) + (8*lb1u*ltu*mmsb12*pow2(
     invdgb1))/3. + (8*lgu*ltu*mmsb12*pow2(invdgb1))/3. - 8*mmsb1*mmsb2*pow2(
     invdgb1) - 2*lb1u*mmsb1*mmsb2*pow2(invdgb1) + (16*lb2u*mmsb1*mmsb2*pow2(
     invdgb1))/3. + (14*lgu*mmsb1*mmsb2*pow2(invdgb1))/3. + (32*mmsb1*mmst1*
     pow2(invdgb1))/3. - 2*lb1u*mmsb1*mmst1*pow2(invdgb1) + (14*lgu*mmsb1*mmst1
     *pow2(invdgb1))/3. - (8*lt1u*mmsb1*mmst1*pow2(invdgb1))/3. + lb1u*lt1u*
     mmsb1*mmst1*pow2(invdgb1) - (7*lgu*lt1u*mmsb1*mmst1*pow2(invdgb1))/3. - 8*
     ltu*mmsb1*mmst1*pow2(invdgb1) + (8*lt1u*ltu*mmsb1*mmst1*pow2(invdgb1))/3.
     + (32*mmsb1*mmst2*pow2(invdgb1))/3. - 2*lb1u*mmsb1*mmst2*pow2(invdgb1) + (
     14*lgu*mmsb1*mmst2*pow2(invdgb1))/3. - (8*lt2u*mmsb1*mmst2*pow2(invdgb1))/
     3. + lb1u*lt2u*mmsb1*mmst2*pow2(invdgb1) - (7*lgu*lt2u*mmsb1*mmst2*pow2(
     invdgb1))/3. - 8*ltu*mmsb1*mmst2*pow2(invdgb1) + (8*lt2u*ltu*mmsb1*mmst2*
     pow2(invdgb1))/3. - 64*mmsb1*mmsusy*pow2(invdgb1) - 16*lb1u*mmsb1*mmsusy*
     pow2(invdgb1) + (112*lgu*mmsb1*mmsusy*pow2(invdgb1))/3. + 16*mmsb1*mmt*
     pow2(invdgb1) + 4*lb1u*mmsb1*mmt*pow2(invdgb1) - (28*lgu*mmsb1*mmt*pow2(
     invdgb1))/3. - (32*ltu*mmsb1*mmt*pow2(invdgb1))/3. - 2*lb1u*ltu*mmsb1*mmt*
     pow2(invdgb1) + (14*lgu*ltu*mmsb1*mmt*pow2(invdgb1))/3. + (16*mmsb12*zt2*
     pow2(invdgb1))/3. - (4*mmsb1*mmsb2*zt2*pow2(invdgb1))/3. + (4*mmsb1*mmst1*
     zt2*pow2(invdgb1))/3. + (4*mmsb1*mmst2*zt2*pow2(invdgb1))/3. - (32*mmsb1*
     mmsusy*zt2*pow2(invdgb1))/3. + (8*mmsb1*mmt*zt2*pow2(invdgb1))/3. - 8*
     mmsb1*mmsb2*pow2(invdgb2) + (16*lb1u*mmsb1*mmsb2*pow2(invdgb2))/3. - 2*
     lb2u*mmsb1*mmsb2*pow2(invdgb2) + (14*lgu*mmsb1*mmsb2*pow2(invdgb2))/3. + (
     20*mmsb22*pow2(invdgb2))/3. + (70*lb2u*mmsb22*pow2(invdgb2))/3. - (70*lgu*
     mmsb22*pow2(invdgb2))/3. - 16*ltu*mmsb22*pow2(invdgb2) + (8*lb2u*ltu*
     mmsb22*pow2(invdgb2))/3. + (8*lgu*ltu*mmsb22*pow2(invdgb2))/3. + (32*mmsb2
     *mmst1*pow2(invdgb2))/3. - 2*lb2u*mmsb2*mmst1*pow2(invdgb2) + (14*lgu*
     mmsb2*mmst1*pow2(invdgb2))/3. - (8*lt1u*mmsb2*mmst1*pow2(invdgb2))/3. +
     lb2u*lt1u*mmsb2*mmst1*pow2(invdgb2) - (7*lgu*lt1u*mmsb2*mmst1*pow2(invdgb2
     ))/3. - 8*ltu*mmsb2*mmst1*pow2(invdgb2) + (8*lt1u*ltu*mmsb2*mmst1*pow2(
     invdgb2))/3. + (32*mmsb2*mmst2*pow2(invdgb2))/3. - 2*lb2u*mmsb2*mmst2*pow2
     (invdgb2) + (14*lgu*mmsb2*mmst2*pow2(invdgb2))/3. - (8*lt2u*mmsb2*mmst2*
     pow2(invdgb2))/3. + lb2u*lt2u*mmsb2*mmst2*pow2(invdgb2) - (7*lgu*lt2u*
     mmsb2*mmst2*pow2(invdgb2))/3. - 8*ltu*mmsb2*mmst2*pow2(invdgb2) + (8*lt2u*
     ltu*mmsb2*mmst2*pow2(invdgb2))/3. - 64*mmsb2*mmsusy*pow2(invdgb2) - 16*
     lb2u*mmsb2*mmsusy*pow2(invdgb2) + (112*lgu*mmsb2*mmsusy*pow2(invdgb2))/3.
     + 16*mmsb2*mmt*pow2(invdgb2) + 4*lb2u*mmsb2*mmt*pow2(invdgb2) - (28*lgu*
     mmsb2*mmt*pow2(invdgb2))/3. - (32*ltu*mmsb2*mmt*pow2(invdgb2))/3. - 2*lb2u
     *ltu*mmsb2*mmt*pow2(invdgb2) + (14*lgu*ltu*mmsb2*mmt*pow2(invdgb2))/3. - (
     4*mmsb1*mmsb2*zt2*pow2(invdgb2))/3. + (16*mmsb22*zt2*pow2(invdgb2))/3. + (
     4*mmsb2*mmst1*zt2*pow2(invdgb2))/3. + (4*mmsb2*mmst2*zt2*pow2(invdgb2))/3.
      - (32*mmsb2*mmsusy*zt2*pow2(invdgb2))/3. + (8*mmsb2*mmt*zt2*pow2(invdgb2)
     )/3. + (4*lsu*(7.333333333333333 - 24*invdgb1*mmsusy - 24*invdgb2*mmsusy +
     32*mmsb1*mmsusy*pow2(invdgb1) + 32*mmsb2*mmsusy*pow2(invdgb2)))/3. + pow2(
     lb1u)/3. + (5*invdgb1*mmsb1*pow2(lb1u))/2. + (invdgb2*mmsb1*pow2(lb1u))/2.
      - (55*mmsb12*pow2(invdgb1)*pow2(lb1u))/6. + (mmsb1*mmsb2*pow2(invdgb1)*
     pow2(lb1u))/2. + (mmsb1*mmst1*pow2(invdgb1)*pow2(lb1u))/2. + (mmsb1*mmst2*
     pow2(invdgb1)*pow2(lb1u))/2. + 4*mmsb1*mmsusy*pow2(invdgb1)*pow2(lb1u) -
     mmsb1*mmt*pow2(invdgb1)*pow2(lb1u) - (2*mmsb1*mmsb2*pow2(invdgb2)*pow2(
     lb1u))/3. + pow2(lb2u)/3. + (invdgb1*mmsb2*pow2(lb2u))/2. + (5*invdgb2*
     mmsb2*pow2(lb2u))/2. - (2*mmsb1*mmsb2*pow2(invdgb1)*pow2(lb2u))/3. + (
     mmsb1*mmsb2*pow2(invdgb2)*pow2(lb2u))/2. - (55*mmsb22*pow2(invdgb2)*pow2(
     lb2u))/6. + (mmsb2*mmst1*pow2(invdgb2)*pow2(lb2u))/2. + (mmsb2*mmst2*pow2(
     invdgb2)*pow2(lb2u))/2. + 4*mmsb2*mmsusy*pow2(invdgb2)*pow2(lb2u) - mmsb2*
     mmt*pow2(invdgb2)*pow2(lb2u) - (4*pow2(lgu))/3. - (37*invdgb1*mmsb1*pow2(
     lgu))/6. + (invdgb2*mmsb1*pow2(lgu))/2. + (invdgb1*mmsb2*pow2(lgu))/2. - (
     37*invdgb2*mmsb2*pow2(lgu))/6. + (invdgb1*mmst1*pow2(lgu))/2. + (invdgb2*
     mmst1*pow2(lgu))/2. + (invdgb1*mmst2*pow2(lgu))/2. + (invdgb2*mmst2*pow2(
     lgu))/2. + 4*invdgb1*mmsusy*pow2(lgu) + 4*invdgb2*mmsusy*pow2(lgu) -
     invdgb1*mmt*pow2(lgu) - invdgb2*mmt*pow2(lgu) + (85*mmsb12*pow2(invdgb1)*
     pow2(lgu))/6. - (7*mmsb1*mmsb2*pow2(invdgb1)*pow2(lgu))/6. - (7*mmsb1*
     mmst1*pow2(invdgb1)*pow2(lgu))/6. - (7*mmsb1*mmst2*pow2(invdgb1)*pow2(lgu)
     )/6. - (28*mmsb1*mmsusy*pow2(invdgb1)*pow2(lgu))/3. + (7*mmsb1*mmt*pow2(
     invdgb1)*pow2(lgu))/3. - (7*mmsb1*mmsb2*pow2(invdgb2)*pow2(lgu))/6. + (85*
     mmsb22*pow2(invdgb2)*pow2(lgu))/6. - (7*mmsb2*mmst1*pow2(invdgb2)*pow2(lgu
     ))/6. - (7*mmsb2*mmst2*pow2(invdgb2)*pow2(lgu))/6. - (28*mmsb2*mmsusy*pow2
     (invdgb2)*pow2(lgu))/3. + (7*mmsb2*mmt*pow2(invdgb2)*pow2(lgu))/3. + (4*(2
      + 3*invdgb1*mmsusy + 3*invdgb2*mmsusy - 4*mmsb1*mmsusy*pow2(invdgb1) - 4*
     mmsb2*mmsusy*pow2(invdgb2))*pow2(lsu))/3. + pow2(lt1u)/3. - (invdgb1*mmst1
     *pow2(lt1u))/2. - (invdgb2*mmst1*pow2(lt1u))/2. + (2*mmsb1*mmst1*pow2(
     invdgb1)*pow2(lt1u))/3. + (2*mmsb2*mmst1*pow2(invdgb2)*pow2(lt1u))/3. +
     pow2(lt2u)/3. - (invdgb1*mmst2*pow2(lt2u))/2. - (invdgb2*mmst2*pow2(lt2u))
     /2. + (2*mmsb1*mmst2*pow2(invdgb1)*pow2(lt2u))/3. + (2*mmsb2*mmst2*pow2(
     invdgb2)*pow2(lt2u))/3. + 4*pow2(ltu) - (10*invdgb1*mmsb1*pow2(ltu))/3. -
     (10*invdgb2*mmsb2*pow2(ltu))/3. - invdgb1*mmst1*pow2(ltu) - invdgb2*mmst1*
     pow2(ltu) - invdgb1*mmst2*pow2(ltu) - invdgb2*mmst2*pow2(ltu) - invdgb1*
     mmt*pow2(ltu) - invdgb2*mmt*pow2(ltu) + (8*mmsb12*pow2(invdgb1)*pow2(ltu))
     /3. + (4*mmsb1*mmst1*pow2(invdgb1)*pow2(ltu))/3. + (4*mmsb1*mmst2*pow2(
     invdgb1)*pow2(ltu))/3. + (4*mmsb1*mmt*pow2(invdgb1)*pow2(ltu))/3. + (8*
     mmsb22*pow2(invdgb2)*pow2(ltu))/3. + (4*mmsb2*mmst1*pow2(invdgb2)*pow2(ltu
     ))/3. + (4*mmsb2*mmst2*pow2(invdgb2)*pow2(ltu))/3. + (4*mmsb2*mmt*pow2(
     invdgb2)*pow2(ltu))/3. - (28*invdgb1*mmsb1*pow2(s2b))/3. + 8*invdgb1*lb1u*
     mmsb1*pow2(s2b) - (112*invdgb1*mmsb2*pow2(s2b))/3. + 28*invdgb2*mmsb2*pow2
     (s2b) + 8*invdgb1*lb1u*mmsb2*pow2(s2b) - 8*invdgb2*lb1u*mmsb2*pow2(s2b) +
     8*invdgb1*lb2u*mmsb2*pow2(s2b) + 16*invdgb1*lgu*mmsb2*pow2(s2b) - 16*
     invdgb2*lgu*mmsb2*pow2(s2b) - (112*invdb12*invdgb1*mmsb22*pow2(s2b))/3. +
     (112*invdb12*invdgb2*mmsb22*pow2(s2b))/3. + 8*invdb12*invdgb1*lb1u*mmsb22*
     pow2(s2b) - 8*invdb12*invdgb2*lb1u*mmsb22*pow2(s2b) + 8*invdb12*invdgb1*
     lb2u*mmsb22*pow2(s2b) - 8*invdb12*invdgb2*lb2u*mmsb22*pow2(s2b) + 16*
     invdb12*invdgb1*lgu*mmsb22*pow2(s2b) - 16*invdb12*invdgb2*lgu*mmsb22*pow2(
     s2b) - (4*invdgb1*mmsb1*zt2*pow2(s2b))/3. - (16*invdgb1*mmsb2*zt2*pow2(s2b
     ))/3. + 4*invdgb2*mmsb2*zt2*pow2(s2b) - (16*invdb12*invdgb1*mmsb22*zt2*
     pow2(s2b))/3. + (16*invdb12*invdgb2*mmsb22*zt2*pow2(s2b))/3. - 8*lb1u*
     mmsb12*pow2(invdgb1)*pow2(s2b) + 8*lgu*mmsb12*pow2(invdgb1)*pow2(s2b) - 8*
     lb2u*mmsb22*pow2(invdgb2)*pow2(s2b) + 8*lgu*mmsb22*pow2(invdgb2)*pow2(s2b)
     - (8*invdgb1*mmsb1*pow2(lb1u)*pow2(s2b))/3. - (4*invdgb1*mmsb2*pow2(lb1u)*
     pow2(s2b))/3. + (4*invdgb2*mmsb2*pow2(lb1u)*pow2(s2b))/3. - (4*invdb12*
     invdgb1*mmsb22*pow2(lb1u)*pow2(s2b))/3. + (4*invdb12*invdgb2*mmsb22*pow2(
     lb1u)*pow2(s2b))/3. + 4*mmsb12*pow2(invdgb1)*pow2(lb1u)*pow2(s2b) - (4*
     invdgb1*mmsb2*pow2(lb2u)*pow2(s2b))/3. - (4*invdgb2*mmsb2*pow2(lb2u)*pow2(
     s2b))/3. - (4*invdb12*invdgb1*mmsb22*pow2(lb2u)*pow2(s2b))/3. + (4*invdb12
     *invdgb2*mmsb22*pow2(lb2u)*pow2(s2b))/3. + 4*mmsb22*pow2(invdgb2)*pow2(
     lb2u)*pow2(s2b) - (8*invdgb1*mmsb2*pow2(lgu)*pow2(s2b))/3. + (8*invdgb2*
     mmsb2*pow2(lgu)*pow2(s2b))/3. - (8*invdb12*invdgb1*mmsb22*pow2(lgu)*pow2(
     s2b))/3. + (8*invdb12*invdgb2*mmsb22*pow2(lgu)*pow2(s2b))/3. - (4*mmsb12*
     pow2(invdgb1)*pow2(lgu)*pow2(s2b))/3. - (4*mmsb22*pow2(invdgb2)*pow2(lgu)*
     pow2(s2b))/3. + 4*zt2*(1 - (7*invdgb1*mmsb1)/2. - 2*invdgb1*mmsb2 - (3*
     invdgb2*mmsb2)/2. - 2*invdb12*invdgb1*mmsb22 + 2*invdb12*invdgb2*mmsb22 +
     3*mmsb12*pow2(invdgb1) + 3*mmsb22*pow2(invdgb2) + (invdgb1*mmsb1*pow2(s2b)
     )/4. - invdgb1*mmsb2*pow2(s2b) + (5*invdgb2*mmsb2*pow2(s2b))/4. - invdb12*
     invdgb1*mmsb22*pow2(s2b) + invdb12*invdgb2*mmsb22*pow2(s2b)) + 4*lb1u*lb2u
     *(invdgb1*mmsb2 - invdgb2*mmsb2 + invdb12*invdgb1*mmsb22 - invdb12*invdgb2
     *mmsb22 - (5*invdgb1*mmsb2*pow2(s2b))/4. + (5*invdgb2*mmsb2*pow2(s2b))/4.
     - (5*invdb12*invdgb1*mmsb22*pow2(s2b))/4. + (5*invdb12*invdgb2*mmsb22*pow2
     (s2b))/4.) + 4*(16.680555555555557 - 36*invdgb1*mmsb1 - 14*invdgb1*mmsb2 -
     22*invdgb2*mmsb2 - 14*invdb12*invdgb1*mmsb22 + 14*invdb12*invdgb2*mmsb22 +
     31*mmsb12*pow2(invdgb1) + 31*mmsb22*pow2(invdgb2) + (7*invdgb1*mmsb1*pow2(
     s2b))/4. - 7*invdgb1*mmsb2*pow2(s2b) + (35*invdgb2*mmsb2*pow2(s2b))/4. - 7
     *invdb12*invdgb1*mmsb22*pow2(s2b) + 7*invdb12*invdgb2*mmsb22*pow2(s2b)) +
     (8*lb1u*mmsb12*mmsb2*pow3(invdgb1))/3. - (8*lgu*mmsb12*mmsb2*pow3(invdgb1)
     )/3. + (8*lb1u*mmsb12*mmst1*pow3(invdgb1))/3. - (8*lgu*mmsb12*mmst1*pow3(
     invdgb1))/3. - (4*lb1u*lt1u*mmsb12*mmst1*pow3(invdgb1))/3. + (4*lgu*lt1u*
     mmsb12*mmst1*pow3(invdgb1))/3. + (8*lb1u*mmsb12*mmst2*pow3(invdgb1))/3. -
     (8*lgu*mmsb12*mmst2*pow3(invdgb1))/3. - (4*lb1u*lt2u*mmsb12*mmst2*pow3(
     invdgb1))/3. + (4*lgu*lt2u*mmsb12*mmst2*pow3(invdgb1))/3. + (64*lb1u*
     mmsb12*mmsusy*pow3(invdgb1))/3. - (64*lgu*mmsb12*mmsusy*pow3(invdgb1))/3.
     - (16*lb1u*mmsb12*mmt*pow3(invdgb1))/3. + (16*lgu*mmsb12*mmt*pow3(invdgb1)
     )/3. + (8*lb1u*ltu*mmsb12*mmt*pow3(invdgb1))/3. - (8*lgu*ltu*mmsb12*mmt*
     pow3(invdgb1))/3. - (2*mmsb12*mmsb2*pow2(lb1u)*pow3(invdgb1))/3. - (2*
     mmsb12*mmst1*pow2(lb1u)*pow3(invdgb1))/3. - (2*mmsb12*mmst2*pow2(lb1u)*
     pow3(invdgb1))/3. - (16*mmsb12*mmsusy*pow2(lb1u)*pow3(invdgb1))/3. + (4*
     mmsb12*mmt*pow2(lb1u)*pow3(invdgb1))/3. + (2*mmsb12*mmsb2*pow2(lgu)*pow3(
     invdgb1))/3. + (2*mmsb12*mmst1*pow2(lgu)*pow3(invdgb1))/3. + (2*mmsb12*
     mmst2*pow2(lgu)*pow3(invdgb1))/3. + (16*mmsb12*mmsusy*pow2(lgu)*pow3(
     invdgb1))/3. - (4*mmsb12*mmt*pow2(lgu)*pow3(invdgb1))/3. + (4*lb1u*lsu*(6*
     mmsb1*mmsusy*pow2(invdgb1) - 8*mmsb12*mmsusy*pow3(invdgb1)))/3. + Fin20(
     mmsb1,mmsusy,mmu)*(s2b*((32*mmgl*mmsb1*pow2(invdgb1))/(3.*mb*mgl) - (32*
     mmgl*mmsusy*pow2(invdgb1))/(3.*mb*mgl)) + (4*(4*invdgb1 - 14*mmsb1*pow2(
     invdgb1) + 6*mmsusy*pow2(invdgb1) + 8*mmsb12*pow3(invdgb1) - 8*mmsb1*
     mmsusy*pow3(invdgb1)))/3.) + (8*lb2u*mmsb1*mmsb22*pow3(invdgb2))/3. - (8*
     lgu*mmsb1*mmsb22*pow3(invdgb2))/3. + (8*lb2u*mmsb22*mmst1*pow3(invdgb2))/
     3. - (8*lgu*mmsb22*mmst1*pow3(invdgb2))/3. - (4*lb2u*lt1u*mmsb22*mmst1*
     pow3(invdgb2))/3. + (4*lgu*lt1u*mmsb22*mmst1*pow3(invdgb2))/3. + (8*lb2u*
     mmsb22*mmst2*pow3(invdgb2))/3. - (8*lgu*mmsb22*mmst2*pow3(invdgb2))/3. - (
     4*lb2u*lt2u*mmsb22*mmst2*pow3(invdgb2))/3. + (4*lgu*lt2u*mmsb22*mmst2*pow3
     (invdgb2))/3. + (64*lb2u*mmsb22*mmsusy*pow3(invdgb2))/3. - (64*lgu*mmsb22*
     mmsusy*pow3(invdgb2))/3. - (16*lb2u*mmsb22*mmt*pow3(invdgb2))/3. + (16*lgu
     *mmsb22*mmt*pow3(invdgb2))/3. + (8*lb2u*ltu*mmsb22*mmt*pow3(invdgb2))/3. -
     (8*lgu*ltu*mmsb22*mmt*pow3(invdgb2))/3. - (2*mmsb1*mmsb22*pow2(lb2u)*pow3(
     invdgb2))/3. - (2*mmsb22*mmst1*pow2(lb2u)*pow3(invdgb2))/3. - (2*mmsb22*
     mmst2*pow2(lb2u)*pow3(invdgb2))/3. - (16*mmsb22*mmsusy*pow2(lb2u)*pow3(
     invdgb2))/3. + (4*mmsb22*mmt*pow2(lb2u)*pow3(invdgb2))/3. + (2*mmsb1*
     mmsb22*pow2(lgu)*pow3(invdgb2))/3. + (2*mmsb22*mmst1*pow2(lgu)*pow3(
     invdgb2))/3. + (2*mmsb22*mmst2*pow2(lgu)*pow3(invdgb2))/3. + (16*mmsb22*
     mmsusy*pow2(lgu)*pow3(invdgb2))/3. - (4*mmsb22*mmt*pow2(lgu)*pow3(invdgb2)
     )/3. + (4*lb1u*lb2u*((3*mmsb1*mmsb2*pow2(invdgb1))/4. + (3*mmsb1*mmsb2*
     pow2(invdgb2))/4. - mmsb12*mmsb2*pow3(invdgb1) - mmsb1*mmsb22*pow3(invdgb2
     )))/3. + (4*lb2u*lsu*(6*mmsb2*mmsusy*pow2(invdgb2) - 8*mmsb22*mmsusy*pow3(
     invdgb2)))/3. + (4*lgu*lsu*(6*invdgb1*mmsusy + 6*invdgb2*mmsusy - 14*mmsb1
     *mmsusy*pow2(invdgb1) - 14*mmsb2*mmsusy*pow2(invdgb2) + 8*mmsb12*mmsusy*
     pow3(invdgb1) + 8*mmsb22*mmsusy*pow3(invdgb2)))/3. + Fin20(mmsb2,mmgl,mmu)
     *(s2b*((-16*invdgb1*mmgl)/(9.*mb*mgl) + (88*invdgb2*mmgl)/(9.*mb*mgl) - (4
     *mmgl*mmsb1*pow2(invdgb1))/(3.*mb*mgl) + (4*mmgl*mmsb2*pow2(invdgb1))/(3.*
     mb*mgl)) + 4*(-invdgb1 - 2*invdb12*invdgb1*mmsb2 + 2*invdb12*invdgb2*mmsb2
      + 3*mmsb2*pow2(invdgb2) - (3*invdgb1*pow2(s2b))/4. + (3*invdgb2*pow2(s2b)
     )/4. - (invdb12*invdgb1*mmsb2*pow2(s2b))/2. + (invdb12*invdgb2*mmsb2*pow2(
     s2b))/2. - (mmsb1*pow2(invdgb1)*pow2(s2b))/2.) + (4*(invdgb1/4. + invdgb2/
     4. + (3*mmsb1*pow2(invdgb1))/4. - (3*mmsb2*pow2(invdgb1))/4. - invdgb1*
     pow2(s2b) + invdgb2*pow2(s2b) - 2*invdb12*invdgb1*mmsb2*pow2(s2b) + 2*
     invdb12*invdgb2*mmsb2*pow2(s2b) - 2*mmsb1*pow2(invdgb1)*pow2(s2b) - mmsb12
     *pow3(invdgb1) + mmsb1*mmsb2*pow3(invdgb1)))/3. + (16*(8*invdb12 + 2*
     invdgb1 - (5*invdgb2)/2. + 4*invdb12*invdgb1*mmsb2 - 12*invdb12*invdgb2*
     mmsb2 - 8*invdb12*pow2(s2b) + (3*invdgb1*pow2(s2b))/2. - (3*invdgb2*pow2(
     s2b))/2. + invdb12*invdgb1*mmsb2*pow2(s2b) + 7*invdb12*invdgb2*mmsb2*pow2(
     s2b) + mmsb1*pow2(invdgb1)*pow2(s2b) - 3*mmsb22*pow3(invdgb2)))/9.) +
     Fin20(mmsb1,mmgl,mmu)*(s2b*((-88*invdgb1*mmgl)/(9.*mb*mgl) + (16*invdgb2*
     mmgl)/(9.*mb*mgl) - (4*mmgl*mmsb1*pow2(invdgb2))/(3.*mb*mgl) + (4*mmgl*
     mmsb2*pow2(invdgb2))/(3.*mb*mgl)) + 4*(-2*invdgb1 + invdgb2 - 2*invdb12*
     invdgb1*mmsb2 + 2*invdb12*invdgb2*mmsb2 + 3*mmsb1*pow2(invdgb1) + (invdgb1
     *pow2(s2b))/4. - (invdgb2*pow2(s2b))/4. - (invdb12*invdgb1*mmsb2*pow2(s2b)
     )/2. + (invdb12*invdgb2*mmsb2*pow2(s2b))/2. - (mmsb2*pow2(invdgb2)*pow2(
     s2b))/2.) + (16*(-8*invdb12 + (19*invdgb1)/2. - 2*invdgb2 + 12*invdb12*
     invdgb1*mmsb2 - 4*invdb12*invdgb2*mmsb2 + 8*invdb12*pow2(s2b) - (17*
     invdgb1*pow2(s2b))/2. + (invdgb2*pow2(s2b))/2. - 7*invdb12*invdgb1*mmsb2*
     pow2(s2b) - invdb12*invdgb2*mmsb2*pow2(s2b) + mmsb2*pow2(invdgb2)*pow2(s2b
     ) - 3*mmsb12*pow3(invdgb1)))/9. + (4*(invdgb1/4. + invdgb2/4. - (3*mmsb1*
     pow2(invdgb2))/4. + (3*mmsb2*pow2(invdgb2))/4. - invdgb1*pow2(s2b) +
     invdgb2*pow2(s2b) - 2*invdb12*invdgb1*mmsb2*pow2(s2b) + 2*invdb12*invdgb2*
     mmsb2*pow2(s2b) - 2*mmsb2*pow2(invdgb2)*pow2(s2b) + mmsb1*mmsb2*pow3(
     invdgb2) - mmsb22*pow3(invdgb2)))/3.) + Fin20(mmsb1,mmsb2,mmu)*(s2b*((4*
     invdgb1*mmgl)/(9.*mb*mgl) - (4*invdgb2*mmgl)/(9.*mb*mgl) + (4*mmgl*mmsb1*
     pow2(invdgb1))/(3.*mb*mgl) - (4*mmgl*mmsb2*pow2(invdgb1))/(3.*mb*mgl) + (4
     *mmgl*mmsb1*pow2(invdgb2))/(3.*mb*mgl) - (4*mmgl*mmsb2*pow2(invdgb2))/(3.*
     mb*mgl)) + (16*(-2*invdgb1 + 2*invdgb2 - 4*invdb12*invdgb1*mmsb2 + 4*
     invdb12*invdgb2*mmsb2 - (invdgb1*pow2(s2b))/2. - (invdgb2*pow2(s2b))/2. -
     mmsb1*pow2(invdgb1)*pow2(s2b) - mmsb2*pow2(invdgb2)*pow2(s2b)))/9. + 4*(
     invdgb1 - invdgb2 + 2*invdb12*invdgb1*mmsb2 - 2*invdb12*invdgb2*mmsb2 + (
     invdgb1*pow2(s2b))/4. + (invdgb2*pow2(s2b))/4. + (mmsb1*pow2(invdgb1)*pow2
     (s2b))/2. + (mmsb2*pow2(invdgb2)*pow2(s2b))/2.) + (4*(invdgb1/2. + invdgb2
     /2. - (7*mmsb1*pow2(invdgb1))/4. + (3*mmsb2*pow2(invdgb1))/4. + (3*mmsb1*
     pow2(invdgb2))/4. - (7*mmsb2*pow2(invdgb2))/4. - invdgb1*pow2(s2b) -
     invdgb2*pow2(s2b) + 2*mmsb1*pow2(invdgb1)*pow2(s2b) + 2*mmsb2*pow2(invdgb2
     )*pow2(s2b) + mmsb12*pow3(invdgb1) - mmsb1*mmsb2*pow3(invdgb1) - mmsb1*
     mmsb2*pow3(invdgb2) + mmsb22*pow3(invdgb2)))/3.) + Fin20(mmsb2,mmsusy,mmu)
     *(s2b*((-32*mmgl*mmsb2*pow2(invdgb2))/(3.*mb*mgl) + (32*mmgl*mmsusy*pow2(
     invdgb2))/(3.*mb*mgl)) + (4*(4*invdgb2 - 14*mmsb2*pow2(invdgb2) + 6*mmsusy
     *pow2(invdgb2) + 8*mmsb22*pow3(invdgb2) - 8*mmsb2*mmsusy*pow3(invdgb2)))/
     3.) + Fin20(mmgl,mmsusy,mmu)*(s2b*((-32*invdgb1*mmgl)/(3.*mb*mgl) + (32*
     invdgb2*mmgl)/(3.*mb*mgl) - (32*mmgl*mmsb1*pow2(invdgb1))/(3.*mb*mgl) + (
     32*mmgl*mmsusy*pow2(invdgb1))/(3.*mb*mgl) + (32*mmgl*mmsb2*pow2(invdgb2))/
     (3.*mb*mgl) - (32*mmgl*mmsusy*pow2(invdgb2))/(3.*mb*mgl)) + (4*(2*invdgb1
     + 2*invdgb2 + 6*mmsb1*pow2(invdgb1) - 6*mmsusy*pow2(invdgb1) + 6*mmsb2*
     pow2(invdgb2) - 6*mmsusy*pow2(invdgb2) - 8*mmsb12*pow3(invdgb1) + 8*mmsb1*
     mmsusy*pow3(invdgb1) - 8*mmsb22*pow3(invdgb2) + 8*mmsb2*mmsusy*pow3(
     invdgb2)))/3.) - (40*lb1u*pow3(invdgb1)*pow3(mmsb1))/3. + (40*lgu*pow3(
     invdgb1)*pow3(mmsb1))/3. + 6*pow2(lb1u)*pow3(invdgb1)*pow3(mmsb1) - (22*
     pow2(lgu)*pow3(invdgb1)*pow3(mmsb1))/3. + 4*lb1u*((3*invdgb1*mmsb1)/2. + (
     11*mmsb12*pow2(invdgb1))/2. - (3*invdgb1*mmsb1*pow2(s2b))/2. + (3*invdgb1*
     mmsb2*pow2(s2b))/2. - (3*invdgb2*mmsb2*pow2(s2b))/2. + (3*invdb12*invdgb1*
     mmsb22*pow2(s2b))/2. - (3*invdb12*invdgb2*mmsb22*pow2(s2b))/2. - (3*mmsb12
     *pow2(invdgb1)*pow2(s2b))/2. - 10*pow3(invdgb1)*pow3(mmsb1)) + (4*lb1u*lgu
     *((3*invdgb1*mmsb1)/4. + (3*invdgb2*mmsb1)/4. - (7*mmsb12*pow2(invdgb1))/
     4. - (7*mmsb1*mmsb2*pow2(invdgb2))/4. - 2*invdgb1*mmsb2*pow2(s2b) + 2*
     invdgb2*mmsb2*pow2(s2b) - 2*invdb12*invdgb1*mmsb22*pow2(s2b) + 2*invdb12*
     invdgb2*mmsb22*pow2(s2b) - 2*mmsb12*pow2(invdgb1)*pow2(s2b) + mmsb1*mmsb22
     *pow3(invdgb2) + pow3(invdgb1)*pow3(mmsb1)))/3. + 4*pow2(lb1u)*((3*invdgb1
     *mmsb1)/2. - (7*mmsb12*pow2(invdgb1))/4. + (invdgb1*mmsb1*pow2(s2b))/2. -
     (invdgb1*mmsb2*pow2(s2b))/4. + (invdgb2*mmsb2*pow2(s2b))/4. - (invdb12*
     invdgb1*mmsb22*pow2(s2b))/4. + (invdb12*invdgb2*mmsb22*pow2(s2b))/4. + (
     mmsb12*pow2(invdgb1)*pow2(s2b))/8. + 2*pow3(invdgb1)*pow3(mmsb1)) + 4*lb1u
     *lgu*(-2*invdgb1*mmsb1 - invdgb1*mmsb2 + invdgb2*mmsb2 - invdb12*invdgb1*
     mmsb22 + invdb12*invdgb2*mmsb22 - mmsb12*pow2(invdgb1) + (3*invdgb1*mmsb2*
     pow2(s2b))/4. - (3*invdgb2*mmsb2*pow2(s2b))/4. + (3*invdb12*invdgb1*mmsb22
     *pow2(s2b))/4. - (3*invdb12*invdgb2*mmsb22*pow2(s2b))/4. + (3*mmsb12*pow2(
     invdgb1)*pow2(s2b))/4. + 2*pow3(invdgb1)*pow3(mmsb1)) + Fin3(mmt,mmsb1,
     mmst1,mmu)*(-(mmsb1*pow2(invdgb1)) + mmst1*pow2(invdgb1) - mmt*pow2(
     invdgb1) + s2b*((4*mmgl*mmsb1*pow2(invdgb1))/(3.*mb*mgl) - (4*mmgl*mmst1*
     pow2(invdgb1))/(3.*mb*mgl) + (4*mmgl*mmt*pow2(invdgb1))/(3.*mb*mgl)) + (4*
     mmsb12*pow3(invdgb1))/3. - (4*mmsb1*mmst1*pow3(invdgb1))/3. + (4*mmsb1*mmt
     *pow3(invdgb1))/3. + s2t*((4*mmgl*mmt*pow2(invdgb1))/(3.*mgl*mt) + s2b*((4
     *invdgb1*mmt)/(3.*mb*mt) - (8*mmsb1*mmt*pow2(invdgb1))/(3.*mb*mt)) - (8*
     mmgl*mmsb1*mmt*pow3(invdgb1))/(3.*mgl*mt)) + DeltaInv(mmt,mmsb1,mmst1)*((2
     *invdgb1*mmsb12)/3. - (4*invdgb1*mmsb1*mmst1)/3. - (2*invdgb1*mmsb1*mmt)/
     3. - (2*invdgb1*mmst1*mmt)/3. + (8*mmsb12*mmst1*pow2(invdgb1))/3. + (4*
     mmsb12*mmt*pow2(invdgb1))/3. + (4*mmsb1*mmst1*mmt*pow2(invdgb1))/3. + (2*
     invdgb1*pow2(mmst1))/3. - (4*mmsb1*pow2(invdgb1)*pow2(mmst1))/3. + s2t*((4
     *mmgl*mmsb12*mmt*pow2(invdgb1))/(3.*mgl*mt) - (4*mmgl*mmsb1*mmst1*mmt*pow2
     (invdgb1))/(3.*mgl*mt) - (4*mmgl*mmsb1*pow2(invdgb1)*pow2(mmt))/(3.*mgl*mt
     )) - (4*pow2(invdgb1)*pow3(mmsb1))/3.)) + Fin3(mmt,mmsb1,mmst2,mmu)*(-(
     mmsb1*pow2(invdgb1)) + mmst2*pow2(invdgb1) - mmt*pow2(invdgb1) + s2b*((4*
     mmgl*mmsb1*pow2(invdgb1))/(3.*mb*mgl) - (4*mmgl*mmst2*pow2(invdgb1))/(3.*
     mb*mgl) + (4*mmgl*mmt*pow2(invdgb1))/(3.*mb*mgl)) + (4*mmsb12*pow3(invdgb1
     ))/3. - (4*mmsb1*mmst2*pow3(invdgb1))/3. + (4*mmsb1*mmt*pow3(invdgb1))/3.
     + s2t*((-4*mmgl*mmt*pow2(invdgb1))/(3.*mgl*mt) + s2b*((-4*invdgb1*mmt)/(3.
     *mb*mt) + (8*mmsb1*mmt*pow2(invdgb1))/(3.*mb*mt)) + (8*mmgl*mmsb1*mmt*pow3
     (invdgb1))/(3.*mgl*mt)) + DeltaInv(mmt,mmsb1,mmst2)*((2*invdgb1*mmsb12)/3.
      - (4*invdgb1*mmsb1*mmst2)/3. - (2*invdgb1*mmsb1*mmt)/3. - (2*invdgb1*
     mmst2*mmt)/3. + (8*mmsb12*mmst2*pow2(invdgb1))/3. + (4*mmsb12*mmt*pow2(
     invdgb1))/3. + (4*mmsb1*mmst2*mmt*pow2(invdgb1))/3. + (2*invdgb1*pow2(
     mmst2))/3. - (4*mmsb1*pow2(invdgb1)*pow2(mmst2))/3. + s2t*((-4*mmgl*mmsb12
     *mmt*pow2(invdgb1))/(3.*mgl*mt) + (4*mmgl*mmsb1*mmst2*mmt*pow2(invdgb1))/(
     3.*mgl*mt) + (4*mmgl*mmsb1*pow2(invdgb1)*pow2(mmt))/(3.*mgl*mt)) - (4*pow2
     (invdgb1)*pow3(mmsb1))/3.)) - (40*lb2u*pow3(invdgb2)*pow3(mmsb2))/3. + (40
     *lgu*pow3(invdgb2)*pow3(mmsb2))/3. + 6*pow2(lb2u)*pow3(invdgb2)*pow3(mmsb2
     ) - (22*pow2(lgu)*pow3(invdgb2)*pow3(mmsb2))/3. + (16*lb1u*(4 + (17*
     invdgb1*mmsb1)/2. + 8*invdb12*mmsb2 - (33*invdgb1*mmsb2)/4. + (invdgb2*
     mmsb2)/4. - 8*invdb12*invdgb1*mmsb22 - (31*mmsb12*pow2(invdgb1))/2. + (
     mmsb1*mmsb2*pow2(invdgb1))/4. + (mmsb22*pow2(invdgb1))/4. - 12*pow2(s2b) +
     (45*invdgb1*mmsb1*pow2(s2b))/2. + (3*invdgb2*mmsb1*pow2(s2b))/2. - 8*
     invdb12*mmsb2*pow2(s2b) + (15*invdgb1*mmsb2*pow2(s2b))/4. + (13*invdgb2*
     mmsb2*pow2(s2b))/4. + 5*invdb12*invdgb1*mmsb22*pow2(s2b) + 3*invdb12*
     invdgb2*mmsb22*pow2(s2b) - (13*mmsb12*pow2(invdgb1)*pow2(s2b))/2. + (9*
     mmsb1*mmsb2*pow2(invdgb1)*pow2(s2b))/4. + (mmsb22*pow2(invdgb1)*pow2(s2b))
     /4. - mmsb1*mmsb2*pow2(invdgb2)*pow2(s2b) - mmsb12*mmsb2*pow2(s2b)*pow3(
     invdgb1) + (65*pow3(invdgb1)*pow3(mmsb1))/4. + (5*pow2(s2b)*pow3(invdgb1)*
     pow3(mmsb1))/4. + (invdgb1*pow2(invdb12)*pow3(mmsb2))/4. - (invdgb2*pow2(
     invdb12)*pow3(mmsb2))/4. + (invdb12*pow2(invdgb1)*pow3(mmsb2))/4. + (
     invdgb1*pow2(invdb12)*pow2(s2b)*pow3(mmsb2))/4. - (invdgb2*pow2(invdb12)*
     pow2(s2b)*pow3(mmsb2))/4. + (invdb12*pow2(invdgb1)*pow2(s2b)*pow3(mmsb2))/
     4.))/9. + (16*(-85.125 + (159*invdgb1*mmsb1)/2. + (335*invdgb1*mmsb2)/4. -
     (17*invdgb2*mmsb2)/4. + (335*invdb12*invdgb1*mmsb22)/4. - (335*invdb12*
     invdgb2*mmsb22)/4. + (23*mmsb12*pow2(invdgb1))/8. + (23*mmsb22*pow2(
     invdgb2))/8. + 40*pow2(s2b) - 50*invdgb1*mmsb1*pow2(s2b) - (3*invdgb2*
     mmsb1*pow2(s2b))/2. - (175*invdgb1*mmsb2*pow2(s2b))/4. - (31*invdgb2*mmsb2
     *pow2(s2b))/4. - (169*invdb12*invdgb1*mmsb22*pow2(s2b))/4. + (169*invdb12*
     invdgb2*mmsb22*pow2(s2b))/4. - (9*mmsb12*pow2(invdgb1)*pow2(s2b))/8. +
     mmsb1*mmsb2*pow2(invdgb1)*pow2(s2b) + mmsb1*mmsb2*pow2(invdgb2)*pow2(s2b)
     - (9*mmsb22*pow2(invdgb2)*pow2(s2b))/8. - 21*pow3(invdgb1)*pow3(mmsb1) -
     21*pow3(invdgb2)*pow3(mmsb2)))/9. + 4*lb2u*((3*invdgb2*mmsb2)/2. + (11*
     mmsb22*pow2(invdgb2))/2. + (3*invdgb1*mmsb2*pow2(s2b))/2. - 3*invdgb2*
     mmsb2*pow2(s2b) + (3*invdb12*invdgb1*mmsb22*pow2(s2b))/2. - (3*invdb12*
     invdgb2*mmsb22*pow2(s2b))/2. - (3*mmsb22*pow2(invdgb2)*pow2(s2b))/2. - 10*
     pow3(invdgb2)*pow3(mmsb2)) + 4*pow2(lgu)*(1.5 - (13*invdgb1*mmsb1)/2. - 3*
     invdgb1*mmsb2 - (7*invdgb2*mmsb2)/2. - 3*invdb12*invdgb1*mmsb22 + 3*
     invdb12*invdgb2*mmsb22 + (35*mmsb12*pow2(invdgb1))/4. + (35*mmsb22*pow2(
     invdgb2))/4. - (7*invdgb1*mmsb2*pow2(s2b))/4. + (7*invdgb2*mmsb2*pow2(s2b)
     )/4. - (7*invdb12*invdgb1*mmsb22*pow2(s2b))/4. + (7*invdb12*invdgb2*mmsb22
     *pow2(s2b))/4. - (7*mmsb12*pow2(invdgb1)*pow2(s2b))/8. - (7*mmsb22*pow2(
     invdgb2)*pow2(s2b))/8. - 4*pow3(invdgb1)*pow3(mmsb1) - 4*pow3(invdgb2)*
     pow3(mmsb2)) + (16*zt2*(-7.5 + (19*invdgb1*mmsb1)/2. + 12*invdgb1*mmsb2 -
     (5*invdgb2*mmsb2)/2. + 12*invdb12*invdgb1*mmsb22 - 12*invdb12*invdgb2*
     mmsb22 + (3*mmsb12*pow2(invdgb1))/2. + (3*mmsb22*pow2(invdgb2))/2. + 8*
     pow2(s2b) - (17*invdgb1*mmsb1*pow2(s2b))/2. - 6*invdgb1*mmsb2*pow2(s2b) -
     (5*invdgb2*mmsb2*pow2(s2b))/2. - 6*invdb12*invdgb1*mmsb22*pow2(s2b) + 6*
     invdb12*invdgb2*mmsb22*pow2(s2b) - 3*pow3(invdgb1)*pow3(mmsb1) - 3*pow3(
     invdgb2)*pow3(mmsb2)))/9. + (4*lb2u*lgu*((3*invdgb1*mmsb2)/4. + (3*invdgb2
     *mmsb2)/4. - (7*mmsb1*mmsb2*pow2(invdgb1))/4. - (7*mmsb22*pow2(invdgb2))/
     4. - 2*invdgb1*mmsb2*pow2(s2b) + 2*invdgb2*mmsb2*pow2(s2b) - 2*invdb12*
     invdgb1*mmsb22*pow2(s2b) + 2*invdb12*invdgb2*mmsb22*pow2(s2b) - 2*mmsb22*
     pow2(invdgb2)*pow2(s2b) + mmsb12*mmsb2*pow3(invdgb1) + pow3(invdgb2)*pow3(
     mmsb2)))/3. + 4*pow2(lb2u)*((3*invdgb2*mmsb2)/2. - (7*mmsb22*pow2(invdgb2)
     )/4. - (invdgb1*mmsb2*pow2(s2b))/4. + (3*invdgb2*mmsb2*pow2(s2b))/4. - (
     invdb12*invdgb1*mmsb22*pow2(s2b))/4. + (invdb12*invdgb2*mmsb22*pow2(s2b))/
     4. + (mmsb22*pow2(invdgb2)*pow2(s2b))/8. + 2*pow3(invdgb2)*pow3(mmsb2)) +
     4*lb2u*lgu*(-(invdgb1*mmsb2) - invdgb2*mmsb2 - invdb12*invdgb1*mmsb22 +
     invdb12*invdgb2*mmsb22 - mmsb22*pow2(invdgb2) + (3*invdgb1*mmsb2*pow2(s2b)
     )/4. - (3*invdgb2*mmsb2*pow2(s2b))/4. + (3*invdb12*invdgb1*mmsb22*pow2(s2b
     ))/4. - (3*invdb12*invdgb2*mmsb22*pow2(s2b))/4. + (3*mmsb22*pow2(invdgb2)*
     pow2(s2b))/4. + 2*pow3(invdgb2)*pow3(mmsb2)) + 4*lgu*(-1.3333333333333333
     + 27*invdgb1*mmsb1 + 12*invdgb1*mmsb2 + 15*invdgb2*mmsb2 + 12*invdb12*
     invdgb1*mmsb22 - 12*invdb12*invdgb2*mmsb22 - (59*mmsb12*pow2(invdgb1))/2.
     - (59*mmsb22*pow2(invdgb2))/2. + 3*invdgb1*mmsb2*pow2(s2b) - 3*invdgb2*
     mmsb2*pow2(s2b) + 3*invdb12*invdgb1*mmsb22*pow2(s2b) - 3*invdb12*invdgb2*
     mmsb22*pow2(s2b) + (3*mmsb12*pow2(invdgb1)*pow2(s2b))/2. + (3*mmsb22*pow2(
     invdgb2)*pow2(s2b))/2. + 10*pow3(invdgb1)*pow3(mmsb1) + 10*pow3(invdgb2)*
     pow3(mmsb2)) + s2b*((64*lb1u*mmgl)/(9.*mb*mgl) - (64*lb2u*mmgl)/(9.*mb*mgl
     ) - (64*lb1u*lgu*mmgl)/(9.*mb*mgl) + (64*lb2u*lgu*mmgl)/(9.*mb*mgl) - (36*
     invdgb1*mmgl*mmsb1)/(mb*mgl) + (8*invdgb2*mmgl*mmsb1)/(mb*mgl) - (56*
     invdgb1*lb1u*mmgl*mmsb1)/(9.*mb*mgl) - (16*invdgb2*lb1u*mmgl*mmsb1)/(3.*mb
     *mgl) + (48*invdgb1*lgu*mmgl*mmsb1)/(mb*mgl) - (8*invdgb2*lgu*mmgl*mmsb1)/
     (3.*mb*mgl) - (100*invdgb1*lb1u*lgu*mmgl*mmsb1)/(9.*mb*mgl) + (4*invdgb2*
     lb1u*lgu*mmgl*mmsb1)/(3.*mb*mgl) - (16*invdgb1*ltu*mmgl*mmsb1)/(mb*mgl) +
     (16*invdgb1*lgu*ltu*mmgl*mmsb1)/(3.*mb*mgl) - (8*invdgb1*mmgl*mmsb2)/(mb*
     mgl) + (36*invdgb2*mmgl*mmsb2)/(mb*mgl) + (8*invdgb1*lb1u*mmgl*mmsb2)/(9.*
     mb*mgl) - (8*invdgb2*lb1u*mmgl*mmsb2)/(9.*mb*mgl) + (40*invdgb1*lb2u*mmgl*
     mmsb2)/(9.*mb*mgl) + (64*invdgb2*lb2u*mmgl*mmsb2)/(9.*mb*mgl) - (8*invdgb1
     *lb1u*lb2u*mmgl*mmsb2)/(3.*mb*mgl) - (40*invdgb2*lb1u*lb2u*mmgl*mmsb2)/(9.
     *mb*mgl) + (8*invdgb1*lgu*mmgl*mmsb2)/(3.*mb*mgl) - (48*invdgb2*lgu*mmgl*
     mmsb2)/(mb*mgl) - (40*invdgb1*lb1u*lgu*mmgl*mmsb2)/(9.*mb*mgl) + (40*
     invdgb2*lb1u*lgu*mmgl*mmsb2)/(9.*mb*mgl) + (4*invdgb1*lb2u*lgu*mmgl*mmsb2)
     /(3.*mb*mgl) + (76*invdgb2*lb2u*lgu*mmgl*mmsb2)/(9.*mb*mgl) + (16*invdgb2*
     ltu*mmgl*mmsb2)/(mb*mgl) - (16*invdgb2*lgu*ltu*mmgl*mmsb2)/(3.*mb*mgl) + (
     8*invdb12*invdgb1*lb1u*mmgl*mmsb22)/(9.*mb*mgl) - (8*invdb12*invdgb2*lb1u*
     mmgl*mmsb22)/(9.*mb*mgl) - (8*invdb12*invdgb1*lb2u*mmgl*mmsb22)/(9.*mb*mgl
     ) + (8*invdb12*invdgb2*lb2u*mmgl*mmsb22)/(9.*mb*mgl) - (8*invdb12*invdgb1*
     lb1u*lb2u*mmgl*mmsb22)/(9.*mb*mgl) - (56*invdb12*invdgb2*lb1u*lb2u*mmgl*
     mmsb22)/(9.*mb*mgl) - (56*invdb12*invdgb1*lb1u*lgu*mmgl*mmsb22)/(9.*mb*mgl
     ) + (56*invdb12*invdgb2*lb1u*lgu*mmgl*mmsb22)/(9.*mb*mgl) + (8*invdb12*
     invdgb1*lb2u*lgu*mmgl*mmsb22)/(9.*mb*mgl) - (8*invdb12*invdgb2*lb2u*lgu*
     mmgl*mmsb22)/(9.*mb*mgl) + (32*invdgb1*mmgl*mmst1)/(3.*mb*mgl) - (32*
     invdgb2*mmgl*mmst1)/(3.*mb*mgl) + (8*invdgb1*lgu*mmgl*mmst1)/(3.*mb*mgl) -
     (8*invdgb2*lgu*mmgl*mmst1)/(3.*mb*mgl) - (8*invdgb1*lt1u*mmgl*mmst1)/(3.*
     mb*mgl) + (8*invdgb2*lt1u*mmgl*mmst1)/(3.*mb*mgl) - (4*invdgb1*lgu*lt1u*
     mmgl*mmst1)/(3.*mb*mgl) + (4*invdgb2*lgu*lt1u*mmgl*mmst1)/(3.*mb*mgl) - (8
     *invdgb1*ltu*mmgl*mmst1)/(mb*mgl) + (8*invdgb2*ltu*mmgl*mmst1)/(mb*mgl) +
     (8*invdgb1*lt1u*ltu*mmgl*mmst1)/(3.*mb*mgl) - (8*invdgb2*lt1u*ltu*mmgl*
     mmst1)/(3.*mb*mgl) + (32*invdgb1*mmgl*mmst2)/(3.*mb*mgl) - (32*invdgb2*
     mmgl*mmst2)/(3.*mb*mgl) + (8*invdgb1*lgu*mmgl*mmst2)/(3.*mb*mgl) - (8*
     invdgb2*lgu*mmgl*mmst2)/(3.*mb*mgl) - (8*invdgb1*lt2u*mmgl*mmst2)/(3.*mb*
     mgl) + (8*invdgb2*lt2u*mmgl*mmst2)/(3.*mb*mgl) - (4*invdgb1*lgu*lt2u*mmgl*
     mmst2)/(3.*mb*mgl) + (4*invdgb2*lgu*lt2u*mmgl*mmst2)/(3.*mb*mgl) - (8*
     invdgb1*ltu*mmgl*mmst2)/(mb*mgl) + (8*invdgb2*ltu*mmgl*mmst2)/(mb*mgl) + (
     8*invdgb1*lt2u*ltu*mmgl*mmst2)/(3.*mb*mgl) - (8*invdgb2*lt2u*ltu*mmgl*
     mmst2)/(3.*mb*mgl) - (64*invdgb1*mmgl*mmsusy)/(mb*mgl) + (64*invdgb2*mmgl*
     mmsusy)/(mb*mgl) + (64*invdgb1*lgu*mmgl*mmsusy)/(3.*mb*mgl) - (64*invdgb2*
     lgu*mmgl*mmsusy)/(3.*mb*mgl) + (128*invdgb1*lsu*mmgl*mmsusy)/(3.*mb*mgl) -
     (128*invdgb2*lsu*mmgl*mmsusy)/(3.*mb*mgl) - (32*invdgb1*lgu*lsu*mmgl*
     mmsusy)/(3.*mb*mgl) + (32*invdgb2*lgu*lsu*mmgl*mmsusy)/(3.*mb*mgl) + (16*
     invdgb1*mmgl*mmt)/(mb*mgl) - (16*invdgb2*mmgl*mmt)/(mb*mgl) - (16*invdgb1*
     lgu*mmgl*mmt)/(3.*mb*mgl) + (16*invdgb2*lgu*mmgl*mmt)/(3.*mb*mgl) - (32*
     invdgb1*ltu*mmgl*mmt)/(3.*mb*mgl) + (32*invdgb2*ltu*mmgl*mmt)/(3.*mb*mgl)
     + (8*invdgb1*lgu*ltu*mmgl*mmt)/(3.*mb*mgl) - (8*invdgb2*lgu*ltu*mmgl*mmt)/
     (3.*mb*mgl) - (40*invdgb1*mmgl*mmsb1*zt2)/(9.*mb*mgl) + (4*invdgb2*mmgl*
     mmsb1*zt2)/(3.*mb*mgl) - (4*invdgb1*mmgl*mmsb2*zt2)/(3.*mb*mgl) + (40*
     invdgb2*mmgl*mmsb2*zt2)/(9.*mb*mgl) + (4*invdgb1*mmgl*mmst1*zt2)/(3.*mb*
     mgl) - (4*invdgb2*mmgl*mmst1*zt2)/(3.*mb*mgl) + (4*invdgb1*mmgl*mmst2*zt2)
     /(3.*mb*mgl) - (4*invdgb2*mmgl*mmst2*zt2)/(3.*mb*mgl) - (32*invdgb1*mmgl*
     mmsusy*zt2)/(3.*mb*mgl) + (32*invdgb2*mmgl*mmsusy*zt2)/(3.*mb*mgl) + (8*
     invdgb1*mmgl*mmt*zt2)/(3.*mb*mgl) - (8*invdgb2*mmgl*mmt*zt2)/(3.*mb*mgl) -
     (40*lb1u*mmgl*mmsb12*pow2(invdgb1))/(mb*mgl) + (40*lgu*mmgl*mmsb12*pow2(
     invdgb1))/(mb*mgl) + (148*lb1u*lgu*mmgl*mmsb12*pow2(invdgb1))/(9.*mb*mgl)
     + (8*lb1u*mmgl*mmsb1*mmsb2*pow2(invdgb1))/(3.*mb*mgl) - (4*lb1u*lb2u*mmgl*
     mmsb1*mmsb2*pow2(invdgb1))/(9.*mb*mgl) - (8*lgu*mmgl*mmsb1*mmsb2*pow2(
     invdgb1))/(3.*mb*mgl) - (8*lb1u*lgu*mmgl*mmsb1*mmsb2*pow2(invdgb1))/(9.*mb
     *mgl) + (4*lb2u*lgu*mmgl*mmsb1*mmsb2*pow2(invdgb1))/(9.*mb*mgl) + (8*lb1u*
     lb2u*mmgl*mmsb22*pow2(invdgb1))/(9.*mb*mgl) - (8*lb1u*lgu*mmgl*mmsb22*pow2
     (invdgb1))/(9.*mb*mgl) - (8*lb2u*lgu*mmgl*mmsb22*pow2(invdgb1))/(9.*mb*mgl
     ) + (8*lb1u*mmgl*mmsb1*mmst1*pow2(invdgb1))/(3.*mb*mgl) - (8*lgu*mmgl*
     mmsb1*mmst1*pow2(invdgb1))/(3.*mb*mgl) - (4*lb1u*lt1u*mmgl*mmsb1*mmst1*
     pow2(invdgb1))/(3.*mb*mgl) + (4*lgu*lt1u*mmgl*mmsb1*mmst1*pow2(invdgb1))/(
     3.*mb*mgl) + (8*lb1u*mmgl*mmsb1*mmst2*pow2(invdgb1))/(3.*mb*mgl) - (8*lgu*
     mmgl*mmsb1*mmst2*pow2(invdgb1))/(3.*mb*mgl) - (4*lb1u*lt2u*mmgl*mmsb1*
     mmst2*pow2(invdgb1))/(3.*mb*mgl) + (4*lgu*lt2u*mmgl*mmsb1*mmst2*pow2(
     invdgb1))/(3.*mb*mgl) + (64*lb1u*mmgl*mmsb1*mmsusy*pow2(invdgb1))/(3.*mb*
     mgl) - (64*lgu*mmgl*mmsb1*mmsusy*pow2(invdgb1))/(3.*mb*mgl) - (32*lb1u*lsu
     *mmgl*mmsb1*mmsusy*pow2(invdgb1))/(3.*mb*mgl) + (32*lgu*lsu*mmgl*mmsb1*
     mmsusy*pow2(invdgb1))/(3.*mb*mgl) - (16*lb1u*mmgl*mmsb1*mmt*pow2(invdgb1))
     /(3.*mb*mgl) + (16*lgu*mmgl*mmsb1*mmt*pow2(invdgb1))/(3.*mb*mgl) + (8*lb1u
     *ltu*mmgl*mmsb1*mmt*pow2(invdgb1))/(3.*mb*mgl) - (8*lgu*ltu*mmgl*mmsb1*mmt
     *pow2(invdgb1))/(3.*mb*mgl) - (8*lb2u*mmgl*mmsb1*mmsb2*pow2(invdgb2))/(3.*
     mb*mgl) + (4*lb1u*lb2u*mmgl*mmsb1*mmsb2*pow2(invdgb2))/(3.*mb*mgl) + (8*
     lgu*mmgl*mmsb1*mmsb2*pow2(invdgb2))/(3.*mb*mgl) - (4*lb1u*lgu*mmgl*mmsb1*
     mmsb2*pow2(invdgb2))/(3.*mb*mgl) + (40*lb2u*mmgl*mmsb22*pow2(invdgb2))/(mb
     *mgl) + (8*lb1u*lb2u*mmgl*mmsb22*pow2(invdgb2))/(9.*mb*mgl) - (40*lgu*mmgl
     *mmsb22*pow2(invdgb2))/(mb*mgl) - (8*lb1u*lgu*mmgl*mmsb22*pow2(invdgb2))/(
     9.*mb*mgl) - (52*lb2u*lgu*mmgl*mmsb22*pow2(invdgb2))/(3.*mb*mgl) - (8*lb2u
     *mmgl*mmsb2*mmst1*pow2(invdgb2))/(3.*mb*mgl) + (8*lgu*mmgl*mmsb2*mmst1*
     pow2(invdgb2))/(3.*mb*mgl) + (4*lb2u*lt1u*mmgl*mmsb2*mmst1*pow2(invdgb2))/
     (3.*mb*mgl) - (4*lgu*lt1u*mmgl*mmsb2*mmst1*pow2(invdgb2))/(3.*mb*mgl) - (8
     *lb2u*mmgl*mmsb2*mmst2*pow2(invdgb2))/(3.*mb*mgl) + (8*lgu*mmgl*mmsb2*
     mmst2*pow2(invdgb2))/(3.*mb*mgl) + (4*lb2u*lt2u*mmgl*mmsb2*mmst2*pow2(
     invdgb2))/(3.*mb*mgl) - (4*lgu*lt2u*mmgl*mmsb2*mmst2*pow2(invdgb2))/(3.*mb
     *mgl) - (64*lb2u*mmgl*mmsb2*mmsusy*pow2(invdgb2))/(3.*mb*mgl) + (64*lgu*
     mmgl*mmsb2*mmsusy*pow2(invdgb2))/(3.*mb*mgl) + (32*lb2u*lsu*mmgl*mmsb2*
     mmsusy*pow2(invdgb2))/(3.*mb*mgl) - (32*lgu*lsu*mmgl*mmsb2*mmsusy*pow2(
     invdgb2))/(3.*mb*mgl) + (16*lb2u*mmgl*mmsb2*mmt*pow2(invdgb2))/(3.*mb*mgl)
     - (16*lgu*mmgl*mmsb2*mmt*pow2(invdgb2))/(3.*mb*mgl) - (8*lb2u*ltu*mmgl*
     mmsb2*mmt*pow2(invdgb2))/(3.*mb*mgl) + (8*lgu*ltu*mmgl*mmsb2*mmt*pow2(
     invdgb2))/(3.*mb*mgl) + (22*invdgb1*mmgl*mmsb1*pow2(lb1u))/(9.*mb*mgl) + (
     2*invdgb2*mmgl*mmsb1*pow2(lb1u))/(3.*mb*mgl) + (32*invdgb1*mmgl*mmsb2*pow2
     (lb1u))/(9.*mb*mgl) + (32*invdb12*invdgb1*mmgl*mmsb22*pow2(lb1u))/(9.*mb*
     mgl) + (62*mmgl*mmsb12*pow2(invdgb1)*pow2(lb1u))/(9.*mb*mgl) - (2*mmgl*
     mmsb1*mmsb2*pow2(invdgb1)*pow2(lb1u))/(3.*mb*mgl) - (2*mmgl*mmsb1*mmst1*
     pow2(invdgb1)*pow2(lb1u))/(3.*mb*mgl) - (2*mmgl*mmsb1*mmst2*pow2(invdgb1)*
     pow2(lb1u))/(3.*mb*mgl) - (16*mmgl*mmsb1*mmsusy*pow2(invdgb1)*pow2(lb1u))/
     (3.*mb*mgl) + (4*mmgl*mmsb1*mmt*pow2(invdgb1)*pow2(lb1u))/(3.*mb*mgl) - (2
     *invdgb1*mmgl*mmsb2*pow2(lb2u))/(3.*mb*mgl) + (10*invdgb2*mmgl*mmsb2*pow2(
     lb2u))/(9.*mb*mgl) + (32*invdb12*invdgb2*mmgl*mmsb22*pow2(lb2u))/(9.*mb*
     mgl) + (2*mmgl*mmsb1*mmsb2*pow2(invdgb2)*pow2(lb2u))/(3.*mb*mgl) - (62*
     mmgl*mmsb22*pow2(invdgb2)*pow2(lb2u))/(9.*mb*mgl) + (2*mmgl*mmsb2*mmst1*
     pow2(invdgb2)*pow2(lb2u))/(3.*mb*mgl) + (2*mmgl*mmsb2*mmst2*pow2(invdgb2)*
     pow2(lb2u))/(3.*mb*mgl) + (16*mmgl*mmsb2*mmsusy*pow2(invdgb2)*pow2(lb2u))/
     (3.*mb*mgl) - (4*mmgl*mmsb2*mmt*pow2(invdgb2)*pow2(lb2u))/(3.*mb*mgl) - (
     74*invdgb1*mmgl*mmsb1*pow2(lgu))/(9.*mb*mgl) + (2*invdgb2*mmgl*mmsb1*pow2(
     lgu))/(3.*mb*mgl) + (2*invdgb1*mmgl*mmsb2*pow2(lgu))/(9.*mb*mgl) + (22*
     invdgb2*mmgl*mmsb2*pow2(lgu))/(3.*mb*mgl) + (8*invdb12*invdgb1*mmgl*mmsb22
     *pow2(lgu))/(3.*mb*mgl) - (8*invdb12*invdgb2*mmgl*mmsb22*pow2(lgu))/(3.*mb
     *mgl) - (2*invdgb1*mmgl*mmst1*pow2(lgu))/(3.*mb*mgl) + (2*invdgb2*mmgl*
     mmst1*pow2(lgu))/(3.*mb*mgl) - (2*invdgb1*mmgl*mmst2*pow2(lgu))/(3.*mb*mgl
     ) + (2*invdgb2*mmgl*mmst2*pow2(lgu))/(3.*mb*mgl) - (16*invdgb1*mmgl*mmsusy
     *pow2(lgu))/(3.*mb*mgl) + (16*invdgb2*mmgl*mmsusy*pow2(lgu))/(3.*mb*mgl) +
     (4*invdgb1*mmgl*mmt*pow2(lgu))/(3.*mb*mgl) - (4*invdgb2*mmgl*mmt*pow2(lgu)
     )/(3.*mb*mgl) - (70*mmgl*mmsb12*pow2(invdgb1)*pow2(lgu))/(3.*mb*mgl) + (14
     *mmgl*mmsb1*mmsb2*pow2(invdgb1)*pow2(lgu))/(9.*mb*mgl) + (8*mmgl*mmsb22*
     pow2(invdgb1)*pow2(lgu))/(9.*mb*mgl) + (2*mmgl*mmsb1*mmst1*pow2(invdgb1)*
     pow2(lgu))/(3.*mb*mgl) + (2*mmgl*mmsb1*mmst2*pow2(invdgb1)*pow2(lgu))/(3.*
     mb*mgl) + (16*mmgl*mmsb1*mmsusy*pow2(invdgb1)*pow2(lgu))/(3.*mb*mgl) - (4*
     mmgl*mmsb1*mmt*pow2(invdgb1)*pow2(lgu))/(3.*mb*mgl) - (2*mmgl*mmsb1*mmsb2*
     pow2(invdgb2)*pow2(lgu))/(3.*mb*mgl) + (218*mmgl*mmsb22*pow2(invdgb2)*pow2
     (lgu))/(9.*mb*mgl) - (2*mmgl*mmsb2*mmst1*pow2(invdgb2)*pow2(lgu))/(3.*mb*
     mgl) - (2*mmgl*mmsb2*mmst2*pow2(invdgb2)*pow2(lgu))/(3.*mb*mgl) - (16*mmgl
     *mmsb2*mmsusy*pow2(invdgb2)*pow2(lgu))/(3.*mb*mgl) + (4*mmgl*mmsb2*mmt*
     pow2(invdgb2)*pow2(lgu))/(3.*mb*mgl) - (16*invdgb1*mmgl*mmsusy*pow2(lsu))/
     (3.*mb*mgl) + (16*invdgb2*mmgl*mmsusy*pow2(lsu))/(3.*mb*mgl) + (2*invdgb1*
     mmgl*mmst1*pow2(lt1u))/(3.*mb*mgl) - (2*invdgb2*mmgl*mmst1*pow2(lt1u))/(3.
     *mb*mgl) + (2*invdgb1*mmgl*mmst2*pow2(lt2u))/(3.*mb*mgl) - (2*invdgb2*mmgl
     *mmst2*pow2(lt2u))/(3.*mb*mgl) + (8*invdgb1*mmgl*mmsb1*pow2(ltu))/(3.*mb*
     mgl) - (8*invdgb2*mmgl*mmsb2*pow2(ltu))/(3.*mb*mgl) + (4*invdgb1*mmgl*
     mmst1*pow2(ltu))/(3.*mb*mgl) - (4*invdgb2*mmgl*mmst1*pow2(ltu))/(3.*mb*mgl
     ) + (4*invdgb1*mmgl*mmst2*pow2(ltu))/(3.*mb*mgl) - (4*invdgb2*mmgl*mmst2*
     pow2(ltu))/(3.*mb*mgl) + (4*invdgb1*mmgl*mmt*pow2(ltu))/(3.*mb*mgl) - (4*
     invdgb2*mmgl*mmt*pow2(ltu))/(3.*mb*mgl) - (16*lb1u*lgu*mmgl*pow3(invdgb1)*
     pow3(mmsb1))/(9.*mb*mgl) + (8*mmgl*pow2(lb1u)*pow3(invdgb1)*pow3(mmsb1))/(
     9.*mb*mgl) + (8*mmgl*pow2(lgu)*pow3(invdgb1)*pow3(mmsb1))/(9.*mb*mgl) + (
     16*invdgb1*lb1u*lb2u*mmgl*pow2(invdb12)*pow3(mmsb2))/(9.*mb*mgl) - (16*
     invdgb2*lb1u*lb2u*mmgl*pow2(invdb12)*pow3(mmsb2))/(9.*mb*mgl) - (16*
     invdgb1*lb1u*lgu*mmgl*pow2(invdb12)*pow3(mmsb2))/(9.*mb*mgl) + (16*invdgb2
     *lb1u*lgu*mmgl*pow2(invdb12)*pow3(mmsb2))/(9.*mb*mgl) - (16*invdgb1*lb2u*
     lgu*mmgl*pow2(invdb12)*pow3(mmsb2))/(9.*mb*mgl) + (16*invdgb2*lb2u*lgu*
     mmgl*pow2(invdb12)*pow3(mmsb2))/(9.*mb*mgl) + (8*invdb12*lb1u*lb2u*mmgl*
     pow2(invdgb1)*pow3(mmsb2))/(9.*mb*mgl) - (8*invdb12*lb1u*lgu*mmgl*pow2(
     invdgb1)*pow3(mmsb2))/(9.*mb*mgl) - (8*invdb12*lb2u*lgu*mmgl*pow2(invdgb1)
     *pow3(mmsb2))/(9.*mb*mgl) + (8*invdb12*lb1u*lb2u*mmgl*pow2(invdgb2)*pow3(
     mmsb2))/(9.*mb*mgl) - (8*invdb12*lb1u*lgu*mmgl*pow2(invdgb2)*pow3(mmsb2))/
     (9.*mb*mgl) - (8*invdb12*lb2u*lgu*mmgl*pow2(invdgb2)*pow3(mmsb2))/(9.*mb*
     mgl) + (16*invdgb1*mmgl*pow2(invdb12)*pow2(lgu)*pow3(mmsb2))/(9.*mb*mgl) -
     (16*invdgb2*mmgl*pow2(invdb12)*pow2(lgu)*pow3(mmsb2))/(9.*mb*mgl) + (8*
     invdb12*mmgl*pow2(invdgb1)*pow2(lgu)*pow3(mmsb2))/(9.*mb*mgl) + (8*invdb12
     *mmgl*pow2(invdgb2)*pow2(lgu)*pow3(mmsb2))/(9.*mb*mgl) + (16*lb2u*lgu*mmgl
     *pow3(invdgb2)*pow3(mmsb2))/(9.*mb*mgl) - (8*mmgl*pow2(lb2u)*pow3(invdgb2)
     *pow3(mmsb2))/(9.*mb*mgl) - (8*mmgl*pow2(lgu)*pow3(invdgb2)*pow3(mmsb2))/(
     9.*mb*mgl)) + (16*lgu*(62 - 75*invdgb1*mmsb1 - (253*invdgb1*mmsb2)/4. - (
     47*invdgb2*mmsb2)/4. - (253*invdb12*invdgb1*mmsb22)/4. + (253*invdb12*
     invdgb2*mmsb22)/4. + (21*mmsb12*pow2(invdgb1))/2. - (mmsb1*mmsb2*pow2(
     invdgb1))/4. - (mmsb22*pow2(invdgb1))/4. + (43*mmsb22*pow2(invdgb2))/4. -
     24*pow2(s2b) + 23*invdgb1*mmsb1*pow2(s2b) + invdgb2*mmsb1*pow2(s2b) + (143
     *invdgb1*mmsb2*pow2(s2b))/4. - (47*invdgb2*mmsb2*pow2(s2b))/4. + (139*
     invdb12*invdgb1*mmsb22*pow2(s2b))/4. - (139*invdb12*invdgb2*mmsb22*pow2(
     s2b))/4. + (15*mmsb12*pow2(invdgb1)*pow2(s2b))/2. - (9*mmsb1*mmsb2*pow2(
     invdgb1)*pow2(s2b))/4. - (mmsb22*pow2(invdgb1)*pow2(s2b))/4. - 2*mmsb1*
     mmsb2*pow2(invdgb2)*pow2(s2b) + (31*mmsb22*pow2(invdgb2)*pow2(s2b))/4. +
     mmsb12*mmsb2*pow2(s2b)*pow3(invdgb1) + mmsb1*mmsb22*pow2(s2b)*pow3(invdgb2
     ) + (7*pow3(invdgb1)*pow3(mmsb1))/4. - (5*pow2(s2b)*pow3(invdgb1)*pow3(
     mmsb1))/4. - (invdb12*pow2(invdgb1)*pow3(mmsb2))/4. + (invdb12*pow2(
     invdgb2)*pow3(mmsb2))/4. - (invdb12*pow2(invdgb1)*pow2(s2b)*pow3(mmsb2))/
     4. + (invdb12*pow2(invdgb2)*pow2(s2b)*pow3(mmsb2))/4. + (7*pow3(invdgb2)*
     pow3(mmsb2))/4. - (5*pow2(s2b)*pow3(invdgb2)*pow3(mmsb2))/4.))/9. + (16*
     lb2u*(-4 - 8*invdb12*mmsb2 - (invdgb1*mmsb2)/2. + 17*invdgb2*mmsb2 - (3*
     invdb12*invdgb1*mmsb22)/4. + (35*invdb12*invdgb2*mmsb22)/4. - (63*mmsb22*
     pow2(invdgb2))/4. - 4*pow2(s2b) - invdgb2*mmsb1*pow2(s2b) + 8*invdb12*
     mmsb2*pow2(s2b) - 2*invdgb1*mmsb2*pow2(s2b) + 18*invdgb2*mmsb2*pow2(s2b) -
     (15*invdb12*invdgb1*mmsb22*pow2(s2b))/4. - (17*invdb12*invdgb2*mmsb22*pow2
     (s2b))/4. - mmsb1*mmsb2*pow2(invdgb1)*pow2(s2b) + 2*mmsb1*mmsb2*pow2(
     invdgb2)*pow2(s2b) - (27*mmsb22*pow2(invdgb2)*pow2(s2b))/4. - mmsb1*mmsb22
     *pow2(s2b)*pow3(invdgb2) - (invdgb1*pow2(invdb12)*pow3(mmsb2))/4. + (
     invdgb2*pow2(invdb12)*pow3(mmsb2))/4. - (invdb12*pow2(invdgb2)*pow3(mmsb2)
     )/4. - (invdgb1*pow2(invdb12)*pow2(s2b)*pow3(mmsb2))/4. + (invdgb2*pow2(
     invdb12)*pow2(s2b)*pow3(mmsb2))/4. - (invdb12*pow2(invdgb2)*pow2(s2b)*pow3
     (mmsb2))/4. + (65*pow3(invdgb2)*pow3(mmsb2))/4. + (5*pow2(s2b)*pow3(
     invdgb2)*pow3(mmsb2))/4.))/9. + Fin3(mmt,mmsb2,mmst1,mmu)*(-(mmsb2*pow2(
     invdgb2)) + mmst1*pow2(invdgb2) - mmt*pow2(invdgb2) + s2b*((-4*mmgl*mmsb2*
     pow2(invdgb2))/(3.*mb*mgl) + (4*mmgl*mmst1*pow2(invdgb2))/(3.*mb*mgl) - (4
     *mmgl*mmt*pow2(invdgb2))/(3.*mb*mgl)) + (4*mmsb22*pow3(invdgb2))/3. - (4*
     mmsb2*mmst1*pow3(invdgb2))/3. + (4*mmsb2*mmt*pow3(invdgb2))/3. + s2t*((4*
     mmgl*mmt*pow2(invdgb2))/(3.*mgl*mt) + s2b*((-4*invdgb2*mmt)/(3.*mb*mt) + (
     8*mmsb2*mmt*pow2(invdgb2))/(3.*mb*mt)) - (8*mmgl*mmsb2*mmt*pow3(invdgb2))/
     (3.*mgl*mt)) + DeltaInv(mmt,mmsb2,mmst1)*((2*invdgb2*mmsb22)/3. - (4*
     invdgb2*mmsb2*mmst1)/3. - (2*invdgb2*mmsb2*mmt)/3. - (2*invdgb2*mmst1*mmt)
     /3. + (8*mmsb22*mmst1*pow2(invdgb2))/3. + (4*mmsb22*mmt*pow2(invdgb2))/3.
     + (4*mmsb2*mmst1*mmt*pow2(invdgb2))/3. + (2*invdgb2*pow2(mmst1))/3. - (4*
     mmsb2*pow2(invdgb2)*pow2(mmst1))/3. + s2t*((4*mmgl*mmsb22*mmt*pow2(invdgb2
     ))/(3.*mgl*mt) - (4*mmgl*mmsb2*mmst1*mmt*pow2(invdgb2))/(3.*mgl*mt) - (4*
     mmgl*mmsb2*pow2(invdgb2)*pow2(mmt))/(3.*mgl*mt)) - (4*pow2(invdgb2)*pow3(
     mmsb2))/3.)) + Fin3(mmt,mmsb2,mmst2,mmu)*(-(mmsb2*pow2(invdgb2)) + mmst2*
     pow2(invdgb2) - mmt*pow2(invdgb2) + s2b*((-4*mmgl*mmsb2*pow2(invdgb2))/(3.
     *mb*mgl) + (4*mmgl*mmst2*pow2(invdgb2))/(3.*mb*mgl) - (4*mmgl*mmt*pow2(
     invdgb2))/(3.*mb*mgl)) + (4*mmsb22*pow3(invdgb2))/3. - (4*mmsb2*mmst2*pow3
     (invdgb2))/3. + (4*mmsb2*mmt*pow3(invdgb2))/3. + s2t*((-4*mmgl*mmt*pow2(
     invdgb2))/(3.*mgl*mt) + s2b*((4*invdgb2*mmt)/(3.*mb*mt) - (8*mmsb2*mmt*
     pow2(invdgb2))/(3.*mb*mt)) + (8*mmgl*mmsb2*mmt*pow3(invdgb2))/(3.*mgl*mt))
     + DeltaInv(mmt,mmsb2,mmst2)*((2*invdgb2*mmsb22)/3. - (4*invdgb2*mmsb2*
     mmst2)/3. - (2*invdgb2*mmsb2*mmt)/3. - (2*invdgb2*mmst2*mmt)/3. + (8*
     mmsb22*mmst2*pow2(invdgb2))/3. + (4*mmsb22*mmt*pow2(invdgb2))/3. + (4*
     mmsb2*mmst2*mmt*pow2(invdgb2))/3. + (2*invdgb2*pow2(mmst2))/3. - (4*mmsb2*
     pow2(invdgb2)*pow2(mmst2))/3. + s2t*((-4*mmgl*mmsb22*mmt*pow2(invdgb2))/(
     3.*mgl*mt) + (4*mmgl*mmsb2*mmst2*mmt*pow2(invdgb2))/(3.*mgl*mt) + (4*mmgl*
     mmsb2*pow2(invdgb2)*pow2(mmt))/(3.*mgl*mt)) - (4*pow2(invdgb2)*pow3(mmsb2)
     )/3.)) + Fin3(mmt,mmst2,mmgl,mmu)*(-invdgb1 - invdgb2 + (7*mmsb1*pow2(
     invdgb1))/3. - mmst2*pow2(invdgb1) + mmt*pow2(invdgb1) + (7*mmsb2*pow2(
     invdgb2))/3. - mmst2*pow2(invdgb2) + mmt*pow2(invdgb2) + s2b*((4*invdgb1*
     mmgl)/(3.*mb*mgl) - (4*invdgb2*mmgl)/(3.*mb*mgl) - (4*mmgl*mmsb1*pow2(
     invdgb1))/(3.*mb*mgl) + (4*mmgl*mmst2*pow2(invdgb1))/(3.*mb*mgl) - (4*mmgl
     *mmt*pow2(invdgb1))/(3.*mb*mgl) + (4*mmgl*mmsb2*pow2(invdgb2))/(3.*mb*mgl)
     - (4*mmgl*mmst2*pow2(invdgb2))/(3.*mb*mgl) + (4*mmgl*mmt*pow2(invdgb2))/(
     3.*mb*mgl)) - (4*mmsb12*pow3(invdgb1))/3. + (4*mmsb1*mmst2*pow3(invdgb1))/
     3. - (4*mmsb1*mmt*pow3(invdgb1))/3. - (4*mmsb22*pow3(invdgb2))/3. + (4*
     mmsb2*mmst2*pow3(invdgb2))/3. - (4*mmsb2*mmt*pow3(invdgb2))/3. + s2t*((4*
     mmgl*mmt*pow2(invdgb1))/(3.*mgl*mt) + (4*mmgl*mmt*pow2(invdgb2))/(3.*mgl*
     mt) + s2b*((4*invdgb1*mmt)/(3.*mb*mt) - (4*invdgb2*mmt)/(3.*mb*mt) - (8*
     mmsb1*mmt*pow2(invdgb1))/(3.*mb*mt) + (8*mmsb2*mmt*pow2(invdgb2))/(3.*mb*
     mt)) - (8*mmgl*mmsb1*mmt*pow3(invdgb1))/(3.*mgl*mt) - (8*mmgl*mmsb2*mmt*
     pow3(invdgb2))/(3.*mgl*mt)) + DeltaInv(mmt,mmst2,mmgl)*((-8*mmgl)/3. - (8*
     mmsb1)/3. + 4*invdgb1*mmsb12 - (8*mmsb2)/3. + 4*invdgb2*mmsb22 + (16*mmst2
     )/3. - (16*invdgb1*mmsb1*mmst2)/3. - (16*invdgb2*mmsb2*mmst2)/3. + (8*mmt)
     /3. - (8*invdgb1*mmsb1*mmt)/3. - (8*invdgb2*mmsb2*mmt)/3. - (4*invdgb1*
     mmst2*mmt)/3. - (4*invdgb2*mmst2*mmt)/3. + (8*mmsb12*mmst2*pow2(invdgb1))/
     3. + (4*mmsb12*mmt*pow2(invdgb1))/3. + (4*mmsb1*mmst2*mmt*pow2(invdgb1))/
     3. + (8*mmsb22*mmst2*pow2(invdgb2))/3. + (4*mmsb22*mmt*pow2(invdgb2))/3. +
     (4*mmsb2*mmst2*mmt*pow2(invdgb2))/3. + (4*invdgb1*pow2(mmst2))/3. + (4*
     invdgb2*pow2(mmst2))/3. - (4*mmsb1*pow2(invdgb1)*pow2(mmst2))/3. - (4*
     mmsb2*pow2(invdgb2)*pow2(mmst2))/3. + s2b*((8*mmgl*mmsb1)/(3.*mb*mgl) - (8
     *invdgb1*mmgl*mmsb12)/(3.*mb*mgl) - (8*mmgl*mmsb2)/(3.*mb*mgl) + (8*
     invdgb2*mmgl*mmsb22)/(3.*mb*mgl) + (16*invdgb1*mmgl*mmsb1*mmst2)/(3.*mb*
     mgl) - (16*invdgb2*mmgl*mmsb2*mmst2)/(3.*mb*mgl) + (8*invdgb1*mmgl*mmsb1*
     mmt)/(3.*mb*mgl) - (8*invdgb2*mmgl*mmsb2*mmt)/(3.*mb*mgl) + (8*invdgb1*
     mmgl*mmst2*mmt)/(3.*mb*mgl) - (8*invdgb2*mmgl*mmst2*mmt)/(3.*mb*mgl) - (8*
     invdgb1*mmgl*pow2(mmst2))/(3.*mb*mgl) + (8*invdgb2*mmgl*pow2(mmst2))/(3.*
     mb*mgl)) + s2t*((-8*mmgl*mmt)/(3.*mgl*mt) + (8*invdgb1*mmgl*mmsb1*mmt)/(3.
     *mgl*mt) + (8*invdgb2*mmgl*mmsb2*mmt)/(3.*mgl*mt) - (4*invdgb1*mmgl*mmst2*
     mmt)/(3.*mgl*mt) - (4*invdgb2*mmgl*mmst2*mmt)/(3.*mgl*mt) - (4*mmgl*mmsb12
     *mmt*pow2(invdgb1))/(3.*mgl*mt) + (4*mmgl*mmsb1*mmst2*mmt*pow2(invdgb1))/(
     3.*mgl*mt) - (4*mmgl*mmsb22*mmt*pow2(invdgb2))/(3.*mgl*mt) + (4*mmgl*mmsb2
     *mmst2*mmt*pow2(invdgb2))/(3.*mgl*mt) - (4*invdgb1*mmgl*pow2(mmt))/(3.*mgl
     *mt) - (4*invdgb2*mmgl*pow2(mmt))/(3.*mgl*mt) + (4*mmgl*mmsb1*pow2(invdgb1
     )*pow2(mmt))/(3.*mgl*mt) + (4*mmgl*mmsb2*pow2(invdgb2)*pow2(mmt))/(3.*mgl*
     mt) + s2b*((8*mmsb1*mmt)/(3.*mb*mt) - (8*invdgb1*mmsb12*mmt)/(3.*mb*mt) -
     (8*mmsb2*mmt)/(3.*mb*mt) + (8*invdgb2*mmsb22*mmt)/(3.*mb*mt) + (8*invdgb1*
     mmsb1*mmst2*mmt)/(3.*mb*mt) - (8*invdgb2*mmsb2*mmst2*mmt)/(3.*mb*mt) + (8*
     invdgb1*mmsb1*pow2(mmt))/(3.*mb*mt) - (8*invdgb2*mmsb2*pow2(mmt))/(3.*mb*
     mt))) - (4*pow2(invdgb1)*pow3(mmsb1))/3. - (4*pow2(invdgb2)*pow3(mmsb2))/
     3.)) + Fin3(mmt,mmst1,mmgl,mmu)*(-invdgb1 - invdgb2 + (7*mmsb1*pow2(
     invdgb1))/3. - mmst1*pow2(invdgb1) + mmt*pow2(invdgb1) + (7*mmsb2*pow2(
     invdgb2))/3. - mmst1*pow2(invdgb2) + mmt*pow2(invdgb2) + s2b*((4*invdgb1*
     mmgl)/(3.*mb*mgl) - (4*invdgb2*mmgl)/(3.*mb*mgl) - (4*mmgl*mmsb1*pow2(
     invdgb1))/(3.*mb*mgl) + (4*mmgl*mmst1*pow2(invdgb1))/(3.*mb*mgl) - (4*mmgl
     *mmt*pow2(invdgb1))/(3.*mb*mgl) + (4*mmgl*mmsb2*pow2(invdgb2))/(3.*mb*mgl)
     - (4*mmgl*mmst1*pow2(invdgb2))/(3.*mb*mgl) + (4*mmgl*mmt*pow2(invdgb2))/(
     3.*mb*mgl)) - (4*mmsb12*pow3(invdgb1))/3. + (4*mmsb1*mmst1*pow3(invdgb1))/
     3. - (4*mmsb1*mmt*pow3(invdgb1))/3. - (4*mmsb22*pow3(invdgb2))/3. + (4*
     mmsb2*mmst1*pow3(invdgb2))/3. - (4*mmsb2*mmt*pow3(invdgb2))/3. + s2t*((-4*
     mmgl*mmt*pow2(invdgb1))/(3.*mgl*mt) - (4*mmgl*mmt*pow2(invdgb2))/(3.*mgl*
     mt) + s2b*((-4*invdgb1*mmt)/(3.*mb*mt) + (4*invdgb2*mmt)/(3.*mb*mt) + (8*
     mmsb1*mmt*pow2(invdgb1))/(3.*mb*mt) - (8*mmsb2*mmt*pow2(invdgb2))/(3.*mb*
     mt)) + (8*mmgl*mmsb1*mmt*pow3(invdgb1))/(3.*mgl*mt) + (8*mmgl*mmsb2*mmt*
     pow3(invdgb2))/(3.*mgl*mt)) + DeltaInv(mmt,mmst1,mmgl)*((-8*mmgl)/3. - (8*
     mmsb1)/3. + 4*invdgb1*mmsb12 - (8*mmsb2)/3. + 4*invdgb2*mmsb22 + (16*mmst1
     )/3. - (16*invdgb1*mmsb1*mmst1)/3. - (16*invdgb2*mmsb2*mmst1)/3. + (8*mmt)
     /3. - (8*invdgb1*mmsb1*mmt)/3. - (8*invdgb2*mmsb2*mmt)/3. - (4*invdgb1*
     mmst1*mmt)/3. - (4*invdgb2*mmst1*mmt)/3. + (8*mmsb12*mmst1*pow2(invdgb1))/
     3. + (4*mmsb12*mmt*pow2(invdgb1))/3. + (4*mmsb1*mmst1*mmt*pow2(invdgb1))/
     3. + (8*mmsb22*mmst1*pow2(invdgb2))/3. + (4*mmsb22*mmt*pow2(invdgb2))/3. +
     (4*mmsb2*mmst1*mmt*pow2(invdgb2))/3. + (4*invdgb1*pow2(mmst1))/3. + (4*
     invdgb2*pow2(mmst1))/3. - (4*mmsb1*pow2(invdgb1)*pow2(mmst1))/3. - (4*
     mmsb2*pow2(invdgb2)*pow2(mmst1))/3. + s2b*((8*mmgl*mmsb1)/(3.*mb*mgl) - (8
     *invdgb1*mmgl*mmsb12)/(3.*mb*mgl) - (8*mmgl*mmsb2)/(3.*mb*mgl) + (8*
     invdgb2*mmgl*mmsb22)/(3.*mb*mgl) + (16*invdgb1*mmgl*mmsb1*mmst1)/(3.*mb*
     mgl) - (16*invdgb2*mmgl*mmsb2*mmst1)/(3.*mb*mgl) + (8*invdgb1*mmgl*mmsb1*
     mmt)/(3.*mb*mgl) - (8*invdgb2*mmgl*mmsb2*mmt)/(3.*mb*mgl) + (8*invdgb1*
     mmgl*mmst1*mmt)/(3.*mb*mgl) - (8*invdgb2*mmgl*mmst1*mmt)/(3.*mb*mgl) - (8*
     invdgb1*mmgl*pow2(mmst1))/(3.*mb*mgl) + (8*invdgb2*mmgl*pow2(mmst1))/(3.*
     mb*mgl)) + s2t*((8*mmgl*mmt)/(3.*mgl*mt) - (8*invdgb1*mmgl*mmsb1*mmt)/(3.*
     mgl*mt) - (8*invdgb2*mmgl*mmsb2*mmt)/(3.*mgl*mt) + (4*invdgb1*mmgl*mmst1*
     mmt)/(3.*mgl*mt) + (4*invdgb2*mmgl*mmst1*mmt)/(3.*mgl*mt) + (4*mmgl*mmsb12
     *mmt*pow2(invdgb1))/(3.*mgl*mt) - (4*mmgl*mmsb1*mmst1*mmt*pow2(invdgb1))/(
     3.*mgl*mt) + (4*mmgl*mmsb22*mmt*pow2(invdgb2))/(3.*mgl*mt) - (4*mmgl*mmsb2
     *mmst1*mmt*pow2(invdgb2))/(3.*mgl*mt) + (4*invdgb1*mmgl*pow2(mmt))/(3.*mgl
     *mt) + (4*invdgb2*mmgl*pow2(mmt))/(3.*mgl*mt) - (4*mmgl*mmsb1*pow2(invdgb1
     )*pow2(mmt))/(3.*mgl*mt) - (4*mmgl*mmsb2*pow2(invdgb2)*pow2(mmt))/(3.*mgl*
     mt) + s2b*((-8*mmsb1*mmt)/(3.*mb*mt) + (8*invdgb1*mmsb12*mmt)/(3.*mb*mt) +
     (8*mmsb2*mmt)/(3.*mb*mt) - (8*invdgb2*mmsb22*mmt)/(3.*mb*mt) - (8*invdgb1*
     mmsb1*mmst1*mmt)/(3.*mb*mt) + (8*invdgb2*mmsb2*mmst1*mmt)/(3.*mb*mt) - (8*
     invdgb1*mmsb1*pow2(mmt))/(3.*mb*mt) + (8*invdgb2*mmsb2*pow2(mmt))/(3.*mb*
     mt))) - (4*pow2(invdgb1)*pow3(mmsb1))/3. - (4*pow2(invdgb2)*pow3(mmsb2))/
     3.)) - (16*invdgb1*mmgl*mmsb1*pow3(s2b))/(9.*mb*mgl) - (16*invdgb2*mmgl*
     mmsb1*pow3(s2b))/(9.*mb*mgl) + (32*invdgb1*lb1u*mmgl*mmsb1*pow3(s2b))/(9.*
     mb*mgl) + (16*invdgb2*lb1u*mmgl*mmsb1*pow3(s2b))/(9.*mb*mgl) - (16*invdgb2
     *lb2u*mmgl*mmsb1*pow3(s2b))/(9.*mb*mgl) + (16*invdgb2*lb1u*lb2u*mmgl*mmsb1
     *pow3(s2b))/(9.*mb*mgl) - (16*invdgb1*lgu*mmgl*mmsb1*pow3(s2b))/(9.*mb*mgl
     ) + (16*invdgb2*lgu*mmgl*mmsb1*pow3(s2b))/(9.*mb*mgl) + (16*invdgb1*lb1u*
     lgu*mmgl*mmsb1*pow3(s2b))/(9.*mb*mgl) - (16*invdgb2*lb1u*lgu*mmgl*mmsb1*
     pow3(s2b))/(9.*mb*mgl) + (16*invdgb1*mmgl*mmsb2*pow3(s2b))/(9.*mb*mgl) + (
     16*invdgb2*mmgl*mmsb2*pow3(s2b))/(9.*mb*mgl) + (16*invdgb1*lb1u*mmgl*mmsb2
     *pow3(s2b))/(9.*mb*mgl) - (16*invdgb1*lb2u*mmgl*mmsb2*pow3(s2b))/(9.*mb*
     mgl) - (32*invdgb2*lb2u*mmgl*mmsb2*pow3(s2b))/(9.*mb*mgl) + (16*invdgb1*
     lb1u*lb2u*mmgl*mmsb2*pow3(s2b))/(9.*mb*mgl) + (32*invdgb2*lb1u*lb2u*mmgl*
     mmsb2*pow3(s2b))/(9.*mb*mgl) - (16*invdgb1*lgu*mmgl*mmsb2*pow3(s2b))/(9.*
     mb*mgl) + (16*invdgb2*lgu*mmgl*mmsb2*pow3(s2b))/(9.*mb*mgl) + (32*invdgb1*
     lb1u*lgu*mmgl*mmsb2*pow3(s2b))/(9.*mb*mgl) - (32*invdgb2*lb1u*lgu*mmgl*
     mmsb2*pow3(s2b))/(9.*mb*mgl) - (16*invdgb1*lb2u*lgu*mmgl*mmsb2*pow3(s2b))/
     (9.*mb*mgl) + (16*invdgb2*lb2u*lgu*mmgl*mmsb2*pow3(s2b))/(9.*mb*mgl) + (32
     *invdb12*invdgb1*lb1u*lb2u*mmgl*mmsb22*pow3(s2b))/(9.*mb*mgl) + (32*
     invdb12*invdgb2*lb1u*lb2u*mmgl*mmsb22*pow3(s2b))/(9.*mb*mgl) + (32*invdb12
     *invdgb1*lb1u*lgu*mmgl*mmsb22*pow3(s2b))/(9.*mb*mgl) - (32*invdb12*invdgb2
     *lb1u*lgu*mmgl*mmsb22*pow3(s2b))/(9.*mb*mgl) - (32*invdb12*invdgb1*lb2u*
     lgu*mmgl*mmsb22*pow3(s2b))/(9.*mb*mgl) + (32*invdb12*invdgb2*lb2u*lgu*mmgl
     *mmsb22*pow3(s2b))/(9.*mb*mgl) + (16*lb1u*mmgl*mmsb12*pow2(invdgb1)*pow3(
     s2b))/(9.*mb*mgl) - (16*lgu*mmgl*mmsb12*pow2(invdgb1)*pow3(s2b))/(9.*mb*
     mgl) + (16*lb1u*lgu*mmgl*mmsb12*pow2(invdgb1)*pow3(s2b))/(9.*mb*mgl) - (16
     *lb1u*mmgl*mmsb1*mmsb2*pow2(invdgb1)*pow3(s2b))/(9.*mb*mgl) + (16*lb1u*
     lb2u*mmgl*mmsb1*mmsb2*pow2(invdgb1)*pow3(s2b))/(9.*mb*mgl) + (16*lgu*mmgl*
     mmsb1*mmsb2*pow2(invdgb1)*pow3(s2b))/(9.*mb*mgl) - (16*lb2u*lgu*mmgl*mmsb1
     *mmsb2*pow2(invdgb1)*pow3(s2b))/(9.*mb*mgl) + (16*lb2u*mmgl*mmsb1*mmsb2*
     pow2(invdgb2)*pow3(s2b))/(9.*mb*mgl) - (16*lb1u*lb2u*mmgl*mmsb1*mmsb2*pow2
     (invdgb2)*pow3(s2b))/(9.*mb*mgl) - (16*lgu*mmgl*mmsb1*mmsb2*pow2(invdgb2)*
     pow3(s2b))/(9.*mb*mgl) + (16*lb1u*lgu*mmgl*mmsb1*mmsb2*pow2(invdgb2)*pow3(
     s2b))/(9.*mb*mgl) - (16*lb2u*mmgl*mmsb22*pow2(invdgb2)*pow3(s2b))/(9.*mb*
     mgl) + (16*lgu*mmgl*mmsb22*pow2(invdgb2)*pow3(s2b))/(9.*mb*mgl) - (16*lb2u
     *lgu*mmgl*mmsb22*pow2(invdgb2)*pow3(s2b))/(9.*mb*mgl) - (16*invdgb1*mmgl*
     mmsb1*pow2(lb1u)*pow3(s2b))/(9.*mb*mgl) - (32*invdgb1*mmgl*mmsb2*pow2(lb1u
     )*pow3(s2b))/(9.*mb*mgl) - (32*invdb12*invdgb1*mmgl*mmsb22*pow2(lb1u)*pow3
     (s2b))/(9.*mb*mgl) - (16*mmgl*mmsb12*pow2(invdgb1)*pow2(lb1u)*pow3(s2b))/(
     9.*mb*mgl) - (16*invdgb2*mmgl*mmsb2*pow2(lb2u)*pow3(s2b))/(9.*mb*mgl) - (
     32*invdb12*invdgb2*mmgl*mmsb22*pow2(lb2u)*pow3(s2b))/(9.*mb*mgl) + (16*
     mmgl*mmsb22*pow2(invdgb2)*pow2(lb2u)*pow3(s2b))/(9.*mb*mgl) + DeltaInv(mmt
     ,mmsb1,mmst1)*((-14*invdgb1*mmsb12*mmst1)/3. + 4*invdgb1*lb1u*mmsb12*mmst1
      - 2*invdgb1*lt1u*mmsb12*mmst1 + 2*invdgb1*ltu*mmsb12*mmst1 - (4*invdgb1*
     lb1u*ltu*mmsb12*mmst1)/3. + (2*invdgb1*lt1u*ltu*mmsb12*mmst1)/3. - (14*
     invdgb1*mmsb12*mmt)/3. + 2*invdgb1*lb1u*mmsb12*mmt + 2*invdgb1*ltu*mmsb12*
     mmt - (2*invdgb1*lb1u*ltu*mmsb12*mmt)/3. - (56*invdgb1*mmsb1*mmst1*mmt)/3.
      + 2*invdgb1*lb1u*mmsb1*mmst1*mmt + 2*invdgb1*lt1u*mmsb1*mmst1*mmt + (4*
     invdgb1*lb1u*lt1u*mmsb1*mmst1*mmt)/3. + 12*invdgb1*ltu*mmsb1*mmst1*mmt - 2
     *invdgb1*lb1u*ltu*mmsb1*mmst1*mmt - 2*invdgb1*lt1u*ltu*mmsb1*mmst1*mmt - (
     2*invdgb1*mmsb12*mmst1*zt2)/3. - (2*invdgb1*mmsb12*mmt*zt2)/3. - (8*
     invdgb1*mmsb1*mmst1*mmt*zt2)/3. + (112*mmsb12*mmst1*mmt*pow2(invdgb1))/3.
     - 4*lb1u*mmsb12*mmst1*mmt*pow2(invdgb1) - 4*lt1u*mmsb12*mmst1*mmt*pow2(
     invdgb1) - (8*lb1u*lt1u*mmsb12*mmst1*mmt*pow2(invdgb1))/3. - 24*ltu*mmsb12
     *mmst1*mmt*pow2(invdgb1) + 4*lb1u*ltu*mmsb12*mmst1*mmt*pow2(invdgb1) + 4*
     lt1u*ltu*mmsb12*mmst1*mmt*pow2(invdgb1) + (16*mmsb12*mmst1*mmt*zt2*pow2(
     invdgb1))/3. - (2*invdgb1*mmsb12*mmst1*pow2(lb1u))/3. - (invdgb1*mmsb12*
     mmt*pow2(lb1u))/3. - (invdgb1*mmsb1*mmst1*mmt*pow2(lb1u))/3. + (2*mmsb12*
     mmst1*mmt*pow2(invdgb1)*pow2(lb1u))/3. + (invdgb1*mmsb12*mmst1*pow2(lt1u))
     /3. - (invdgb1*mmsb1*mmst1*mmt*pow2(lt1u))/3. + (2*mmsb12*mmst1*mmt*pow2(
     invdgb1)*pow2(lt1u))/3. - (invdgb1*mmsb12*mmst1*pow2(ltu))/3. - (invdgb1*
     mmsb12*mmt*pow2(ltu))/3. - 2*invdgb1*mmsb1*mmst1*mmt*pow2(ltu) + 4*mmsb12*
     mmst1*mmt*pow2(invdgb1)*pow2(ltu) - (14*invdgb1*mmsb1*pow2(mmst1))/3. - 2*
     invdgb1*lb1u*mmsb1*pow2(mmst1) + 4*invdgb1*lt1u*mmsb1*pow2(mmst1) + 2*
     invdgb1*ltu*mmsb1*pow2(mmst1) + (2*invdgb1*lb1u*ltu*mmsb1*pow2(mmst1))/3.
     - (4*invdgb1*lt1u*ltu*mmsb1*pow2(mmst1))/3. - (14*invdgb1*mmt*pow2(mmst1))
     /3. + 2*invdgb1*lt1u*mmt*pow2(mmst1) + 2*invdgb1*ltu*mmt*pow2(mmst1) - (2*
     invdgb1*lt1u*ltu*mmt*pow2(mmst1))/3. - (2*invdgb1*mmsb1*zt2*pow2(mmst1))/
     3. - (2*invdgb1*mmt*zt2*pow2(mmst1))/3. + (28*mmsb12*pow2(invdgb1)*pow2(
     mmst1))/3. + 4*lb1u*mmsb12*pow2(invdgb1)*pow2(mmst1) - 8*lt1u*mmsb12*pow2(
     invdgb1)*pow2(mmst1) - 4*ltu*mmsb12*pow2(invdgb1)*pow2(mmst1) - (4*lb1u*
     ltu*mmsb12*pow2(invdgb1)*pow2(mmst1))/3. + (8*lt1u*ltu*mmsb12*pow2(invdgb1
     )*pow2(mmst1))/3. + (28*mmsb1*mmt*pow2(invdgb1)*pow2(mmst1))/3. - 4*lt1u*
     mmsb1*mmt*pow2(invdgb1)*pow2(mmst1) - 4*ltu*mmsb1*mmt*pow2(invdgb1)*pow2(
     mmst1) + (4*lt1u*ltu*mmsb1*mmt*pow2(invdgb1)*pow2(mmst1))/3. + (4*mmsb12*
     zt2*pow2(invdgb1)*pow2(mmst1))/3. + (4*mmsb1*mmt*zt2*pow2(invdgb1)*pow2(
     mmst1))/3. + (invdgb1*mmsb1*pow2(lb1u)*pow2(mmst1))/3. - (2*mmsb12*pow2(
     invdgb1)*pow2(lb1u)*pow2(mmst1))/3. - (2*invdgb1*mmsb1*pow2(lt1u)*pow2(
     mmst1))/3. - (invdgb1*mmt*pow2(lt1u)*pow2(mmst1))/3. + (4*mmsb12*pow2(
     invdgb1)*pow2(lt1u)*pow2(mmst1))/3. + (2*mmsb1*mmt*pow2(invdgb1)*pow2(lt1u
     )*pow2(mmst1))/3. - (invdgb1*mmsb1*pow2(ltu)*pow2(mmst1))/3. - (invdgb1*
     mmt*pow2(ltu)*pow2(mmst1))/3. + (2*mmsb12*pow2(invdgb1)*pow2(ltu)*pow2(
     mmst1))/3. + (2*mmsb1*mmt*pow2(invdgb1)*pow2(ltu)*pow2(mmst1))/3. + (14*
     invdgb1*pow3(mmsb1))/3. - 2*invdgb1*lb1u*pow3(mmsb1) - 2*invdgb1*ltu*pow3(
     mmsb1) + (2*invdgb1*lb1u*ltu*pow3(mmsb1))/3. + (2*invdgb1*zt2*pow3(mmsb1))
     /3. + (28*mmst1*pow2(invdgb1)*pow3(mmsb1))/3. - 8*lb1u*mmst1*pow2(invdgb1)
     *pow3(mmsb1) + 4*lt1u*mmst1*pow2(invdgb1)*pow3(mmsb1) - 4*ltu*mmst1*pow2(
     invdgb1)*pow3(mmsb1) + (8*lb1u*ltu*mmst1*pow2(invdgb1)*pow3(mmsb1))/3. - (
     4*lt1u*ltu*mmst1*pow2(invdgb1)*pow3(mmsb1))/3. + (28*mmt*pow2(invdgb1)*
     pow3(mmsb1))/3. - 4*lb1u*mmt*pow2(invdgb1)*pow3(mmsb1) - 4*ltu*mmt*pow2(
     invdgb1)*pow3(mmsb1) + (4*lb1u*ltu*mmt*pow2(invdgb1)*pow3(mmsb1))/3. + (4*
     mmst1*zt2*pow2(invdgb1)*pow3(mmsb1))/3. + (4*mmt*zt2*pow2(invdgb1)*pow3(
     mmsb1))/3. + (invdgb1*pow2(lb1u)*pow3(mmsb1))/3. + (4*mmst1*pow2(invdgb1)*
     pow2(lb1u)*pow3(mmsb1))/3. + (2*mmt*pow2(invdgb1)*pow2(lb1u)*pow3(mmsb1))/
     3. - (2*mmst1*pow2(invdgb1)*pow2(lt1u)*pow3(mmsb1))/3. + (invdgb1*pow2(ltu
     )*pow3(mmsb1))/3. + (2*mmst1*pow2(invdgb1)*pow2(ltu)*pow3(mmsb1))/3. + (2*
     mmt*pow2(invdgb1)*pow2(ltu)*pow3(mmsb1))/3. + s2t*((-28*mmgl*mmsb12*mmst1*
     mmt*pow2(invdgb1))/(3.*mgl*mt) + (4*lb1u*mmgl*mmsb12*mmst1*mmt*pow2(
     invdgb1))/(mgl*mt) - (4*lt1u*mmgl*mmsb12*mmst1*mmt*pow2(invdgb1))/(mgl*mt)
     + (4*lb1u*lt1u*mmgl*mmsb12*mmst1*mmt*pow2(invdgb1))/(3.*mgl*mt) + (8*ltu*
     mmgl*mmsb12*mmst1*mmt*pow2(invdgb1))/(mgl*mt) - (8*lb1u*ltu*mmgl*mmsb12*
     mmst1*mmt*pow2(invdgb1))/(3.*mgl*mt) - (4*mmgl*mmsb12*mmst1*mmt*zt2*pow2(
     invdgb1))/(3.*mgl*mt) - (2*mmgl*mmsb12*mmst1*mmt*pow2(invdgb1)*pow2(lb1u))
     /(3.*mgl*mt) + (2*mmgl*mmsb12*mmst1*mmt*pow2(invdgb1)*pow2(lt1u))/(3.*mgl*
     mt) - (4*mmgl*mmsb12*mmst1*mmt*pow2(invdgb1)*pow2(ltu))/(3.*mgl*mt) + (4*
     lt1u*mmgl*mmsb1*mmt*pow2(invdgb1)*pow2(mmst1))/(mgl*mt) - (4*lb1u*lt1u*
     mmgl*mmsb1*mmt*pow2(invdgb1)*pow2(mmst1))/(3.*mgl*mt) - (4*ltu*mmgl*mmsb1*
     mmt*pow2(invdgb1)*pow2(mmst1))/(mgl*mt) + (4*lb1u*ltu*mmgl*mmsb1*mmt*pow2(
     invdgb1)*pow2(mmst1))/(3.*mgl*mt) - (2*mmgl*mmsb1*mmt*pow2(invdgb1)*pow2(
     lt1u)*pow2(mmst1))/(3.*mgl*mt) + (2*mmgl*mmsb1*mmt*pow2(invdgb1)*pow2(ltu)
     *pow2(mmst1))/(3.*mgl*mt) - (28*mmgl*mmsb12*pow2(invdgb1)*pow2(mmt))/(3.*
     mgl*mt) + (4*lb1u*mmgl*mmsb12*pow2(invdgb1)*pow2(mmt))/(mgl*mt) + (4*ltu*
     mmgl*mmsb12*pow2(invdgb1)*pow2(mmt))/(mgl*mt) - (4*lb1u*ltu*mmgl*mmsb12*
     pow2(invdgb1)*pow2(mmt))/(3.*mgl*mt) - (56*mmgl*mmsb1*mmst1*pow2(invdgb1)*
     pow2(mmt))/(3.*mgl*mt) + (4*lt1u*mmgl*mmsb1*mmst1*pow2(invdgb1)*pow2(mmt))
     /(mgl*mt) + (4*lb1u*lt1u*mmgl*mmsb1*mmst1*pow2(invdgb1)*pow2(mmt))/(3.*mgl
     *mt) + (12*ltu*mmgl*mmsb1*mmst1*pow2(invdgb1)*pow2(mmt))/(mgl*mt) - (4*
     lb1u*ltu*mmgl*mmsb1*mmst1*pow2(invdgb1)*pow2(mmt))/(3.*mgl*mt) - (8*lt1u*
     ltu*mmgl*mmsb1*mmst1*pow2(invdgb1)*pow2(mmt))/(3.*mgl*mt) - (4*mmgl*mmsb12
     *zt2*pow2(invdgb1)*pow2(mmt))/(3.*mgl*mt) - (8*mmgl*mmsb1*mmst1*zt2*pow2(
     invdgb1)*pow2(mmt))/(3.*mgl*mt) - (2*mmgl*mmsb12*pow2(invdgb1)*pow2(lb1u)*
     pow2(mmt))/(3.*mgl*mt) - (2*mmgl*mmsb1*mmst1*pow2(invdgb1)*pow2(lt1u)*pow2
     (mmt))/(3.*mgl*mt) - (2*mmgl*mmsb12*pow2(invdgb1)*pow2(ltu)*pow2(mmt))/(3.
     *mgl*mt) - (2*mmgl*mmsb1*mmst1*pow2(invdgb1)*pow2(ltu)*pow2(mmt))/(mgl*mt)
     + (28*mmgl*mmt*pow2(invdgb1)*pow3(mmsb1))/(3.*mgl*mt) - (4*lb1u*mmgl*mmt*
     pow2(invdgb1)*pow3(mmsb1))/(mgl*mt) - (4*ltu*mmgl*mmt*pow2(invdgb1)*pow3(
     mmsb1))/(mgl*mt) + (4*lb1u*ltu*mmgl*mmt*pow2(invdgb1)*pow3(mmsb1))/(3.*mgl
     *mt) + (4*mmgl*mmt*zt2*pow2(invdgb1)*pow3(mmsb1))/(3.*mgl*mt) + (2*mmgl*
     mmt*pow2(invdgb1)*pow2(lb1u)*pow3(mmsb1))/(3.*mgl*mt) + (2*mmgl*mmt*pow2(
     invdgb1)*pow2(ltu)*pow3(mmsb1))/(3.*mgl*mt)) + (14*invdgb1*pow3(mmst1))/3.
      - 2*invdgb1*lt1u*pow3(mmst1) - 2*invdgb1*ltu*pow3(mmst1) + (2*invdgb1*
     lt1u*ltu*pow3(mmst1))/3. + (2*invdgb1*zt2*pow3(mmst1))/3. - (28*mmsb1*pow2
     (invdgb1)*pow3(mmst1))/3. + 4*lt1u*mmsb1*pow2(invdgb1)*pow3(mmst1) + 4*ltu
     *mmsb1*pow2(invdgb1)*pow3(mmst1) - (4*lt1u*ltu*mmsb1*pow2(invdgb1)*pow3(
     mmst1))/3. - (4*mmsb1*zt2*pow2(invdgb1)*pow3(mmst1))/3. + (invdgb1*pow2(
     lt1u)*pow3(mmst1))/3. - (2*mmsb1*pow2(invdgb1)*pow2(lt1u)*pow3(mmst1))/3.
     + (invdgb1*pow2(ltu)*pow3(mmst1))/3. - (2*mmsb1*pow2(invdgb1)*pow2(ltu)*
     pow3(mmst1))/3. - (28*pow2(invdgb1)*pow4(mmsb1))/3. + 4*lb1u*pow2(invdgb1)
     *pow4(mmsb1) + 4*ltu*pow2(invdgb1)*pow4(mmsb1) - (4*lb1u*ltu*pow2(invdgb1)
     *pow4(mmsb1))/3. - (4*zt2*pow2(invdgb1)*pow4(mmsb1))/3. - (2*pow2(invdgb1)
     *pow2(lb1u)*pow4(mmsb1))/3. - (2*pow2(invdgb1)*pow2(ltu)*pow4(mmsb1))/3.)
     + DeltaInv(mmt,mmsb1,mmst2)*((-14*invdgb1*mmsb12*mmst2)/3. + 4*invdgb1*
     lb1u*mmsb12*mmst2 - 2*invdgb1*lt2u*mmsb12*mmst2 + 2*invdgb1*ltu*mmsb12*
     mmst2 - (4*invdgb1*lb1u*ltu*mmsb12*mmst2)/3. + (2*invdgb1*lt2u*ltu*mmsb12*
     mmst2)/3. - (14*invdgb1*mmsb12*mmt)/3. + 2*invdgb1*lb1u*mmsb12*mmt + 2*
     invdgb1*ltu*mmsb12*mmt - (2*invdgb1*lb1u*ltu*mmsb12*mmt)/3. - (56*invdgb1*
     mmsb1*mmst2*mmt)/3. + 2*invdgb1*lb1u*mmsb1*mmst2*mmt + 2*invdgb1*lt2u*
     mmsb1*mmst2*mmt + (4*invdgb1*lb1u*lt2u*mmsb1*mmst2*mmt)/3. + 12*invdgb1*
     ltu*mmsb1*mmst2*mmt - 2*invdgb1*lb1u*ltu*mmsb1*mmst2*mmt - 2*invdgb1*lt2u*
     ltu*mmsb1*mmst2*mmt - (2*invdgb1*mmsb12*mmst2*zt2)/3. - (2*invdgb1*mmsb12*
     mmt*zt2)/3. - (8*invdgb1*mmsb1*mmst2*mmt*zt2)/3. + (112*mmsb12*mmst2*mmt*
     pow2(invdgb1))/3. - 4*lb1u*mmsb12*mmst2*mmt*pow2(invdgb1) - 4*lt2u*mmsb12*
     mmst2*mmt*pow2(invdgb1) - (8*lb1u*lt2u*mmsb12*mmst2*mmt*pow2(invdgb1))/3.
     - 24*ltu*mmsb12*mmst2*mmt*pow2(invdgb1) + 4*lb1u*ltu*mmsb12*mmst2*mmt*pow2
     (invdgb1) + 4*lt2u*ltu*mmsb12*mmst2*mmt*pow2(invdgb1) + (16*mmsb12*mmst2*
     mmt*zt2*pow2(invdgb1))/3. - (2*invdgb1*mmsb12*mmst2*pow2(lb1u))/3. - (
     invdgb1*mmsb12*mmt*pow2(lb1u))/3. - (invdgb1*mmsb1*mmst2*mmt*pow2(lb1u))/
     3. + (2*mmsb12*mmst2*mmt*pow2(invdgb1)*pow2(lb1u))/3. + (invdgb1*mmsb12*
     mmst2*pow2(lt2u))/3. - (invdgb1*mmsb1*mmst2*mmt*pow2(lt2u))/3. + (2*mmsb12
     *mmst2*mmt*pow2(invdgb1)*pow2(lt2u))/3. - (invdgb1*mmsb12*mmst2*pow2(ltu))
     /3. - (invdgb1*mmsb12*mmt*pow2(ltu))/3. - 2*invdgb1*mmsb1*mmst2*mmt*pow2(
     ltu) + 4*mmsb12*mmst2*mmt*pow2(invdgb1)*pow2(ltu) - (14*invdgb1*mmsb1*pow2
     (mmst2))/3. - 2*invdgb1*lb1u*mmsb1*pow2(mmst2) + 4*invdgb1*lt2u*mmsb1*pow2
     (mmst2) + 2*invdgb1*ltu*mmsb1*pow2(mmst2) + (2*invdgb1*lb1u*ltu*mmsb1*pow2
     (mmst2))/3. - (4*invdgb1*lt2u*ltu*mmsb1*pow2(mmst2))/3. - (14*invdgb1*mmt*
     pow2(mmst2))/3. + 2*invdgb1*lt2u*mmt*pow2(mmst2) + 2*invdgb1*ltu*mmt*pow2(
     mmst2) - (2*invdgb1*lt2u*ltu*mmt*pow2(mmst2))/3. - (2*invdgb1*mmsb1*zt2*
     pow2(mmst2))/3. - (2*invdgb1*mmt*zt2*pow2(mmst2))/3. + (28*mmsb12*pow2(
     invdgb1)*pow2(mmst2))/3. + 4*lb1u*mmsb12*pow2(invdgb1)*pow2(mmst2) - 8*
     lt2u*mmsb12*pow2(invdgb1)*pow2(mmst2) - 4*ltu*mmsb12*pow2(invdgb1)*pow2(
     mmst2) - (4*lb1u*ltu*mmsb12*pow2(invdgb1)*pow2(mmst2))/3. + (8*lt2u*ltu*
     mmsb12*pow2(invdgb1)*pow2(mmst2))/3. + (28*mmsb1*mmt*pow2(invdgb1)*pow2(
     mmst2))/3. - 4*lt2u*mmsb1*mmt*pow2(invdgb1)*pow2(mmst2) - 4*ltu*mmsb1*mmt*
     pow2(invdgb1)*pow2(mmst2) + (4*lt2u*ltu*mmsb1*mmt*pow2(invdgb1)*pow2(mmst2
     ))/3. + (4*mmsb12*zt2*pow2(invdgb1)*pow2(mmst2))/3. + (4*mmsb1*mmt*zt2*
     pow2(invdgb1)*pow2(mmst2))/3. + (invdgb1*mmsb1*pow2(lb1u)*pow2(mmst2))/3.
     - (2*mmsb12*pow2(invdgb1)*pow2(lb1u)*pow2(mmst2))/3. - (2*invdgb1*mmsb1*
     pow2(lt2u)*pow2(mmst2))/3. - (invdgb1*mmt*pow2(lt2u)*pow2(mmst2))/3. + (4*
     mmsb12*pow2(invdgb1)*pow2(lt2u)*pow2(mmst2))/3. + (2*mmsb1*mmt*pow2(
     invdgb1)*pow2(lt2u)*pow2(mmst2))/3. - (invdgb1*mmsb1*pow2(ltu)*pow2(mmst2)
     )/3. - (invdgb1*mmt*pow2(ltu)*pow2(mmst2))/3. + (2*mmsb12*pow2(invdgb1)*
     pow2(ltu)*pow2(mmst2))/3. + (2*mmsb1*mmt*pow2(invdgb1)*pow2(ltu)*pow2(
     mmst2))/3. + (14*invdgb1*pow3(mmsb1))/3. - 2*invdgb1*lb1u*pow3(mmsb1) - 2*
     invdgb1*ltu*pow3(mmsb1) + (2*invdgb1*lb1u*ltu*pow3(mmsb1))/3. + (2*invdgb1
     *zt2*pow3(mmsb1))/3. + (28*mmst2*pow2(invdgb1)*pow3(mmsb1))/3. - 8*lb1u*
     mmst2*pow2(invdgb1)*pow3(mmsb1) + 4*lt2u*mmst2*pow2(invdgb1)*pow3(mmsb1) -
     4*ltu*mmst2*pow2(invdgb1)*pow3(mmsb1) + (8*lb1u*ltu*mmst2*pow2(invdgb1)*
     pow3(mmsb1))/3. - (4*lt2u*ltu*mmst2*pow2(invdgb1)*pow3(mmsb1))/3. + (28*
     mmt*pow2(invdgb1)*pow3(mmsb1))/3. - 4*lb1u*mmt*pow2(invdgb1)*pow3(mmsb1) -
     4*ltu*mmt*pow2(invdgb1)*pow3(mmsb1) + (4*lb1u*ltu*mmt*pow2(invdgb1)*pow3(
     mmsb1))/3. + (4*mmst2*zt2*pow2(invdgb1)*pow3(mmsb1))/3. + (4*mmt*zt2*pow2(
     invdgb1)*pow3(mmsb1))/3. + (invdgb1*pow2(lb1u)*pow3(mmsb1))/3. + (4*mmst2*
     pow2(invdgb1)*pow2(lb1u)*pow3(mmsb1))/3. + (2*mmt*pow2(invdgb1)*pow2(lb1u)
     *pow3(mmsb1))/3. - (2*mmst2*pow2(invdgb1)*pow2(lt2u)*pow3(mmsb1))/3. + (
     invdgb1*pow2(ltu)*pow3(mmsb1))/3. + (2*mmst2*pow2(invdgb1)*pow2(ltu)*pow3(
     mmsb1))/3. + (2*mmt*pow2(invdgb1)*pow2(ltu)*pow3(mmsb1))/3. + s2t*((28*
     mmgl*mmsb12*mmst2*mmt*pow2(invdgb1))/(3.*mgl*mt) - (4*lb1u*mmgl*mmsb12*
     mmst2*mmt*pow2(invdgb1))/(mgl*mt) + (4*lt2u*mmgl*mmsb12*mmst2*mmt*pow2(
     invdgb1))/(mgl*mt) - (4*lb1u*lt2u*mmgl*mmsb12*mmst2*mmt*pow2(invdgb1))/(3.
     *mgl*mt) - (8*ltu*mmgl*mmsb12*mmst2*mmt*pow2(invdgb1))/(mgl*mt) + (8*lb1u*
     ltu*mmgl*mmsb12*mmst2*mmt*pow2(invdgb1))/(3.*mgl*mt) + (4*mmgl*mmsb12*
     mmst2*mmt*zt2*pow2(invdgb1))/(3.*mgl*mt) + (2*mmgl*mmsb12*mmst2*mmt*pow2(
     invdgb1)*pow2(lb1u))/(3.*mgl*mt) - (2*mmgl*mmsb12*mmst2*mmt*pow2(invdgb1)*
     pow2(lt2u))/(3.*mgl*mt) + (4*mmgl*mmsb12*mmst2*mmt*pow2(invdgb1)*pow2(ltu)
     )/(3.*mgl*mt) - (4*lt2u*mmgl*mmsb1*mmt*pow2(invdgb1)*pow2(mmst2))/(mgl*mt)
     + (4*lb1u*lt2u*mmgl*mmsb1*mmt*pow2(invdgb1)*pow2(mmst2))/(3.*mgl*mt) + (4*
     ltu*mmgl*mmsb1*mmt*pow2(invdgb1)*pow2(mmst2))/(mgl*mt) - (4*lb1u*ltu*mmgl*
     mmsb1*mmt*pow2(invdgb1)*pow2(mmst2))/(3.*mgl*mt) + (2*mmgl*mmsb1*mmt*pow2(
     invdgb1)*pow2(lt2u)*pow2(mmst2))/(3.*mgl*mt) - (2*mmgl*mmsb1*mmt*pow2(
     invdgb1)*pow2(ltu)*pow2(mmst2))/(3.*mgl*mt) + (28*mmgl*mmsb12*pow2(invdgb1
     )*pow2(mmt))/(3.*mgl*mt) - (4*lb1u*mmgl*mmsb12*pow2(invdgb1)*pow2(mmt))/(
     mgl*mt) - (4*ltu*mmgl*mmsb12*pow2(invdgb1)*pow2(mmt))/(mgl*mt) + (4*lb1u*
     ltu*mmgl*mmsb12*pow2(invdgb1)*pow2(mmt))/(3.*mgl*mt) + (56*mmgl*mmsb1*
     mmst2*pow2(invdgb1)*pow2(mmt))/(3.*mgl*mt) - (4*lt2u*mmgl*mmsb1*mmst2*pow2
     (invdgb1)*pow2(mmt))/(mgl*mt) - (4*lb1u*lt2u*mmgl*mmsb1*mmst2*pow2(invdgb1
     )*pow2(mmt))/(3.*mgl*mt) - (12*ltu*mmgl*mmsb1*mmst2*pow2(invdgb1)*pow2(mmt
     ))/(mgl*mt) + (4*lb1u*ltu*mmgl*mmsb1*mmst2*pow2(invdgb1)*pow2(mmt))/(3.*
     mgl*mt) + (8*lt2u*ltu*mmgl*mmsb1*mmst2*pow2(invdgb1)*pow2(mmt))/(3.*mgl*mt
     ) + (4*mmgl*mmsb12*zt2*pow2(invdgb1)*pow2(mmt))/(3.*mgl*mt) + (8*mmgl*
     mmsb1*mmst2*zt2*pow2(invdgb1)*pow2(mmt))/(3.*mgl*mt) + (2*mmgl*mmsb12*pow2
     (invdgb1)*pow2(lb1u)*pow2(mmt))/(3.*mgl*mt) + (2*mmgl*mmsb1*mmst2*pow2(
     invdgb1)*pow2(lt2u)*pow2(mmt))/(3.*mgl*mt) + (2*mmgl*mmsb12*pow2(invdgb1)*
     pow2(ltu)*pow2(mmt))/(3.*mgl*mt) + (2*mmgl*mmsb1*mmst2*pow2(invdgb1)*pow2(
     ltu)*pow2(mmt))/(mgl*mt) - (28*mmgl*mmt*pow2(invdgb1)*pow3(mmsb1))/(3.*mgl
     *mt) + (4*lb1u*mmgl*mmt*pow2(invdgb1)*pow3(mmsb1))/(mgl*mt) + (4*ltu*mmgl*
     mmt*pow2(invdgb1)*pow3(mmsb1))/(mgl*mt) - (4*lb1u*ltu*mmgl*mmt*pow2(
     invdgb1)*pow3(mmsb1))/(3.*mgl*mt) - (4*mmgl*mmt*zt2*pow2(invdgb1)*pow3(
     mmsb1))/(3.*mgl*mt) - (2*mmgl*mmt*pow2(invdgb1)*pow2(lb1u)*pow3(mmsb1))/(
     3.*mgl*mt) - (2*mmgl*mmt*pow2(invdgb1)*pow2(ltu)*pow3(mmsb1))/(3.*mgl*mt))
     + (14*invdgb1*pow3(mmst2))/3. - 2*invdgb1*lt2u*pow3(mmst2) - 2*invdgb1*ltu
     *pow3(mmst2) + (2*invdgb1*lt2u*ltu*pow3(mmst2))/3. + (2*invdgb1*zt2*pow3(
     mmst2))/3. - (28*mmsb1*pow2(invdgb1)*pow3(mmst2))/3. + 4*lt2u*mmsb1*pow2(
     invdgb1)*pow3(mmst2) + 4*ltu*mmsb1*pow2(invdgb1)*pow3(mmst2) - (4*lt2u*ltu
     *mmsb1*pow2(invdgb1)*pow3(mmst2))/3. - (4*mmsb1*zt2*pow2(invdgb1)*pow3(
     mmst2))/3. + (invdgb1*pow2(lt2u)*pow3(mmst2))/3. - (2*mmsb1*pow2(invdgb1)*
     pow2(lt2u)*pow3(mmst2))/3. + (invdgb1*pow2(ltu)*pow3(mmst2))/3. - (2*mmsb1
     *pow2(invdgb1)*pow2(ltu)*pow3(mmst2))/3. - (28*pow2(invdgb1)*pow4(mmsb1))/
     3. + 4*lb1u*pow2(invdgb1)*pow4(mmsb1) + 4*ltu*pow2(invdgb1)*pow4(mmsb1) -
     (4*lb1u*ltu*pow2(invdgb1)*pow4(mmsb1))/3. - (4*zt2*pow2(invdgb1)*pow4(
     mmsb1))/3. - (2*pow2(invdgb1)*pow2(lb1u)*pow4(mmsb1))/3. - (2*pow2(invdgb1
     )*pow2(ltu)*pow4(mmsb1))/3.) + (16*pow2(lb1u)*(-4 - (9*invdgb1*mmsb1)/4. -
     4*invdb12*mmsb2 + 4*invdgb1*mmsb2 + 4*invdb12*invdgb1*mmsb22 + (15*mmsb12*
     pow2(invdgb1))/2. + 4*pow2(s2b) - 6*invdgb1*mmsb1*pow2(s2b) + 4*invdb12*
     mmsb2*pow2(s2b) - (7*invdgb1*mmsb2*pow2(s2b))/2. - (invdgb2*mmsb2*pow2(s2b
     ))/2. - (7*invdb12*invdgb1*mmsb22*pow2(s2b))/2. - (invdb12*invdgb2*mmsb22*
     pow2(s2b))/2. + (5*mmsb12*pow2(invdgb1)*pow2(s2b))/4. - 6*pow3(invdgb1)*
     pow3(mmsb1) - (pow2(s2b)*pow3(invdgb1)*pow3(mmsb1))/2. + (3*pow4(invdgb1)*
     pow4(mmsb1))/8. - (pow2(s2b)*pow4(invdgb1)*pow4(mmsb1))/8.))/9. + DeltaInv
     (mmt,mmsb2,mmst1)*((-14*invdgb2*mmsb22*mmst1)/3. + 4*invdgb2*lb2u*mmsb22*
     mmst1 - 2*invdgb2*lt1u*mmsb22*mmst1 + 2*invdgb2*ltu*mmsb22*mmst1 - (4*
     invdgb2*lb2u*ltu*mmsb22*mmst1)/3. + (2*invdgb2*lt1u*ltu*mmsb22*mmst1)/3. -
     (14*invdgb2*mmsb22*mmt)/3. + 2*invdgb2*lb2u*mmsb22*mmt + 2*invdgb2*ltu*
     mmsb22*mmt - (2*invdgb2*lb2u*ltu*mmsb22*mmt)/3. - (56*invdgb2*mmsb2*mmst1*
     mmt)/3. + 2*invdgb2*lb2u*mmsb2*mmst1*mmt + 2*invdgb2*lt1u*mmsb2*mmst1*mmt
     + (4*invdgb2*lb2u*lt1u*mmsb2*mmst1*mmt)/3. + 12*invdgb2*ltu*mmsb2*mmst1*
     mmt - 2*invdgb2*lb2u*ltu*mmsb2*mmst1*mmt - 2*invdgb2*lt1u*ltu*mmsb2*mmst1*
     mmt - (2*invdgb2*mmsb22*mmst1*zt2)/3. - (2*invdgb2*mmsb22*mmt*zt2)/3. - (8
     *invdgb2*mmsb2*mmst1*mmt*zt2)/3. + (112*mmsb22*mmst1*mmt*pow2(invdgb2))/3.
      - 4*lb2u*mmsb22*mmst1*mmt*pow2(invdgb2) - 4*lt1u*mmsb22*mmst1*mmt*pow2(
     invdgb2) - (8*lb2u*lt1u*mmsb22*mmst1*mmt*pow2(invdgb2))/3. - 24*ltu*mmsb22
     *mmst1*mmt*pow2(invdgb2) + 4*lb2u*ltu*mmsb22*mmst1*mmt*pow2(invdgb2) + 4*
     lt1u*ltu*mmsb22*mmst1*mmt*pow2(invdgb2) + (16*mmsb22*mmst1*mmt*zt2*pow2(
     invdgb2))/3. - (2*invdgb2*mmsb22*mmst1*pow2(lb2u))/3. - (invdgb2*mmsb22*
     mmt*pow2(lb2u))/3. - (invdgb2*mmsb2*mmst1*mmt*pow2(lb2u))/3. + (2*mmsb22*
     mmst1*mmt*pow2(invdgb2)*pow2(lb2u))/3. + (invdgb2*mmsb22*mmst1*pow2(lt1u))
     /3. - (invdgb2*mmsb2*mmst1*mmt*pow2(lt1u))/3. + (2*mmsb22*mmst1*mmt*pow2(
     invdgb2)*pow2(lt1u))/3. - (invdgb2*mmsb22*mmst1*pow2(ltu))/3. - (invdgb2*
     mmsb22*mmt*pow2(ltu))/3. - 2*invdgb2*mmsb2*mmst1*mmt*pow2(ltu) + 4*mmsb22*
     mmst1*mmt*pow2(invdgb2)*pow2(ltu) - (14*invdgb2*mmsb2*pow2(mmst1))/3. - 2*
     invdgb2*lb2u*mmsb2*pow2(mmst1) + 4*invdgb2*lt1u*mmsb2*pow2(mmst1) + 2*
     invdgb2*ltu*mmsb2*pow2(mmst1) + (2*invdgb2*lb2u*ltu*mmsb2*pow2(mmst1))/3.
     - (4*invdgb2*lt1u*ltu*mmsb2*pow2(mmst1))/3. - (14*invdgb2*mmt*pow2(mmst1))
     /3. + 2*invdgb2*lt1u*mmt*pow2(mmst1) + 2*invdgb2*ltu*mmt*pow2(mmst1) - (2*
     invdgb2*lt1u*ltu*mmt*pow2(mmst1))/3. - (2*invdgb2*mmsb2*zt2*pow2(mmst1))/
     3. - (2*invdgb2*mmt*zt2*pow2(mmst1))/3. + (28*mmsb22*pow2(invdgb2)*pow2(
     mmst1))/3. + 4*lb2u*mmsb22*pow2(invdgb2)*pow2(mmst1) - 8*lt1u*mmsb22*pow2(
     invdgb2)*pow2(mmst1) - 4*ltu*mmsb22*pow2(invdgb2)*pow2(mmst1) - (4*lb2u*
     ltu*mmsb22*pow2(invdgb2)*pow2(mmst1))/3. + (8*lt1u*ltu*mmsb22*pow2(invdgb2
     )*pow2(mmst1))/3. + (28*mmsb2*mmt*pow2(invdgb2)*pow2(mmst1))/3. - 4*lt1u*
     mmsb2*mmt*pow2(invdgb2)*pow2(mmst1) - 4*ltu*mmsb2*mmt*pow2(invdgb2)*pow2(
     mmst1) + (4*lt1u*ltu*mmsb2*mmt*pow2(invdgb2)*pow2(mmst1))/3. + (4*mmsb22*
     zt2*pow2(invdgb2)*pow2(mmst1))/3. + (4*mmsb2*mmt*zt2*pow2(invdgb2)*pow2(
     mmst1))/3. + (invdgb2*mmsb2*pow2(lb2u)*pow2(mmst1))/3. - (2*mmsb22*pow2(
     invdgb2)*pow2(lb2u)*pow2(mmst1))/3. - (2*invdgb2*mmsb2*pow2(lt1u)*pow2(
     mmst1))/3. - (invdgb2*mmt*pow2(lt1u)*pow2(mmst1))/3. + (4*mmsb22*pow2(
     invdgb2)*pow2(lt1u)*pow2(mmst1))/3. + (2*mmsb2*mmt*pow2(invdgb2)*pow2(lt1u
     )*pow2(mmst1))/3. - (invdgb2*mmsb2*pow2(ltu)*pow2(mmst1))/3. - (invdgb2*
     mmt*pow2(ltu)*pow2(mmst1))/3. + (2*mmsb22*pow2(invdgb2)*pow2(ltu)*pow2(
     mmst1))/3. + (2*mmsb2*mmt*pow2(invdgb2)*pow2(ltu)*pow2(mmst1))/3. + (14*
     invdgb2*pow3(mmsb2))/3. - 2*invdgb2*lb2u*pow3(mmsb2) - 2*invdgb2*ltu*pow3(
     mmsb2) + (2*invdgb2*lb2u*ltu*pow3(mmsb2))/3. + (2*invdgb2*zt2*pow3(mmsb2))
     /3. + (28*mmst1*pow2(invdgb2)*pow3(mmsb2))/3. - 8*lb2u*mmst1*pow2(invdgb2)
     *pow3(mmsb2) + 4*lt1u*mmst1*pow2(invdgb2)*pow3(mmsb2) - 4*ltu*mmst1*pow2(
     invdgb2)*pow3(mmsb2) + (8*lb2u*ltu*mmst1*pow2(invdgb2)*pow3(mmsb2))/3. - (
     4*lt1u*ltu*mmst1*pow2(invdgb2)*pow3(mmsb2))/3. + (28*mmt*pow2(invdgb2)*
     pow3(mmsb2))/3. - 4*lb2u*mmt*pow2(invdgb2)*pow3(mmsb2) - 4*ltu*mmt*pow2(
     invdgb2)*pow3(mmsb2) + (4*lb2u*ltu*mmt*pow2(invdgb2)*pow3(mmsb2))/3. + (4*
     mmst1*zt2*pow2(invdgb2)*pow3(mmsb2))/3. + (4*mmt*zt2*pow2(invdgb2)*pow3(
     mmsb2))/3. + (invdgb2*pow2(lb2u)*pow3(mmsb2))/3. + (4*mmst1*pow2(invdgb2)*
     pow2(lb2u)*pow3(mmsb2))/3. + (2*mmt*pow2(invdgb2)*pow2(lb2u)*pow3(mmsb2))/
     3. - (2*mmst1*pow2(invdgb2)*pow2(lt1u)*pow3(mmsb2))/3. + (invdgb2*pow2(ltu
     )*pow3(mmsb2))/3. + (2*mmst1*pow2(invdgb2)*pow2(ltu)*pow3(mmsb2))/3. + (2*
     mmt*pow2(invdgb2)*pow2(ltu)*pow3(mmsb2))/3. + s2t*((-28*mmgl*mmsb22*mmst1*
     mmt*pow2(invdgb2))/(3.*mgl*mt) + (4*lb2u*mmgl*mmsb22*mmst1*mmt*pow2(
     invdgb2))/(mgl*mt) - (4*lt1u*mmgl*mmsb22*mmst1*mmt*pow2(invdgb2))/(mgl*mt)
     + (4*lb2u*lt1u*mmgl*mmsb22*mmst1*mmt*pow2(invdgb2))/(3.*mgl*mt) + (8*ltu*
     mmgl*mmsb22*mmst1*mmt*pow2(invdgb2))/(mgl*mt) - (8*lb2u*ltu*mmgl*mmsb22*
     mmst1*mmt*pow2(invdgb2))/(3.*mgl*mt) - (4*mmgl*mmsb22*mmst1*mmt*zt2*pow2(
     invdgb2))/(3.*mgl*mt) - (2*mmgl*mmsb22*mmst1*mmt*pow2(invdgb2)*pow2(lb2u))
     /(3.*mgl*mt) + (2*mmgl*mmsb22*mmst1*mmt*pow2(invdgb2)*pow2(lt1u))/(3.*mgl*
     mt) - (4*mmgl*mmsb22*mmst1*mmt*pow2(invdgb2)*pow2(ltu))/(3.*mgl*mt) + (4*
     lt1u*mmgl*mmsb2*mmt*pow2(invdgb2)*pow2(mmst1))/(mgl*mt) - (4*lb2u*lt1u*
     mmgl*mmsb2*mmt*pow2(invdgb2)*pow2(mmst1))/(3.*mgl*mt) - (4*ltu*mmgl*mmsb2*
     mmt*pow2(invdgb2)*pow2(mmst1))/(mgl*mt) + (4*lb2u*ltu*mmgl*mmsb2*mmt*pow2(
     invdgb2)*pow2(mmst1))/(3.*mgl*mt) - (2*mmgl*mmsb2*mmt*pow2(invdgb2)*pow2(
     lt1u)*pow2(mmst1))/(3.*mgl*mt) + (2*mmgl*mmsb2*mmt*pow2(invdgb2)*pow2(ltu)
     *pow2(mmst1))/(3.*mgl*mt) - (28*mmgl*mmsb22*pow2(invdgb2)*pow2(mmt))/(3.*
     mgl*mt) + (4*lb2u*mmgl*mmsb22*pow2(invdgb2)*pow2(mmt))/(mgl*mt) + (4*ltu*
     mmgl*mmsb22*pow2(invdgb2)*pow2(mmt))/(mgl*mt) - (4*lb2u*ltu*mmgl*mmsb22*
     pow2(invdgb2)*pow2(mmt))/(3.*mgl*mt) - (56*mmgl*mmsb2*mmst1*pow2(invdgb2)*
     pow2(mmt))/(3.*mgl*mt) + (4*lt1u*mmgl*mmsb2*mmst1*pow2(invdgb2)*pow2(mmt))
     /(mgl*mt) + (4*lb2u*lt1u*mmgl*mmsb2*mmst1*pow2(invdgb2)*pow2(mmt))/(3.*mgl
     *mt) + (12*ltu*mmgl*mmsb2*mmst1*pow2(invdgb2)*pow2(mmt))/(mgl*mt) - (4*
     lb2u*ltu*mmgl*mmsb2*mmst1*pow2(invdgb2)*pow2(mmt))/(3.*mgl*mt) - (8*lt1u*
     ltu*mmgl*mmsb2*mmst1*pow2(invdgb2)*pow2(mmt))/(3.*mgl*mt) - (4*mmgl*mmsb22
     *zt2*pow2(invdgb2)*pow2(mmt))/(3.*mgl*mt) - (8*mmgl*mmsb2*mmst1*zt2*pow2(
     invdgb2)*pow2(mmt))/(3.*mgl*mt) - (2*mmgl*mmsb22*pow2(invdgb2)*pow2(lb2u)*
     pow2(mmt))/(3.*mgl*mt) - (2*mmgl*mmsb2*mmst1*pow2(invdgb2)*pow2(lt1u)*pow2
     (mmt))/(3.*mgl*mt) - (2*mmgl*mmsb22*pow2(invdgb2)*pow2(ltu)*pow2(mmt))/(3.
     *mgl*mt) - (2*mmgl*mmsb2*mmst1*pow2(invdgb2)*pow2(ltu)*pow2(mmt))/(mgl*mt)
     + (28*mmgl*mmt*pow2(invdgb2)*pow3(mmsb2))/(3.*mgl*mt) - (4*lb2u*mmgl*mmt*
     pow2(invdgb2)*pow3(mmsb2))/(mgl*mt) - (4*ltu*mmgl*mmt*pow2(invdgb2)*pow3(
     mmsb2))/(mgl*mt) + (4*lb2u*ltu*mmgl*mmt*pow2(invdgb2)*pow3(mmsb2))/(3.*mgl
     *mt) + (4*mmgl*mmt*zt2*pow2(invdgb2)*pow3(mmsb2))/(3.*mgl*mt) + (2*mmgl*
     mmt*pow2(invdgb2)*pow2(lb2u)*pow3(mmsb2))/(3.*mgl*mt) + (2*mmgl*mmt*pow2(
     invdgb2)*pow2(ltu)*pow3(mmsb2))/(3.*mgl*mt)) + (14*invdgb2*pow3(mmst1))/3.
      - 2*invdgb2*lt1u*pow3(mmst1) - 2*invdgb2*ltu*pow3(mmst1) + (2*invdgb2*
     lt1u*ltu*pow3(mmst1))/3. + (2*invdgb2*zt2*pow3(mmst1))/3. - (28*mmsb2*pow2
     (invdgb2)*pow3(mmst1))/3. + 4*lt1u*mmsb2*pow2(invdgb2)*pow3(mmst1) + 4*ltu
     *mmsb2*pow2(invdgb2)*pow3(mmst1) - (4*lt1u*ltu*mmsb2*pow2(invdgb2)*pow3(
     mmst1))/3. - (4*mmsb2*zt2*pow2(invdgb2)*pow3(mmst1))/3. + (invdgb2*pow2(
     lt1u)*pow3(mmst1))/3. - (2*mmsb2*pow2(invdgb2)*pow2(lt1u)*pow3(mmst1))/3.
     + (invdgb2*pow2(ltu)*pow3(mmst1))/3. - (2*mmsb2*pow2(invdgb2)*pow2(ltu)*
     pow3(mmst1))/3. - (28*pow2(invdgb2)*pow4(mmsb2))/3. + 4*lb2u*pow2(invdgb2)
     *pow4(mmsb2) + 4*ltu*pow2(invdgb2)*pow4(mmsb2) - (4*lb2u*ltu*pow2(invdgb2)
     *pow4(mmsb2))/3. - (4*zt2*pow2(invdgb2)*pow4(mmsb2))/3. - (2*pow2(invdgb2)
     *pow2(lb2u)*pow4(mmsb2))/3. - (2*pow2(invdgb2)*pow2(ltu)*pow4(mmsb2))/3.)
     + DeltaInv(mmt,mmsb2,mmst2)*((-14*invdgb2*mmsb22*mmst2)/3. + 4*invdgb2*
     lb2u*mmsb22*mmst2 - 2*invdgb2*lt2u*mmsb22*mmst2 + 2*invdgb2*ltu*mmsb22*
     mmst2 - (4*invdgb2*lb2u*ltu*mmsb22*mmst2)/3. + (2*invdgb2*lt2u*ltu*mmsb22*
     mmst2)/3. - (14*invdgb2*mmsb22*mmt)/3. + 2*invdgb2*lb2u*mmsb22*mmt + 2*
     invdgb2*ltu*mmsb22*mmt - (2*invdgb2*lb2u*ltu*mmsb22*mmt)/3. - (56*invdgb2*
     mmsb2*mmst2*mmt)/3. + 2*invdgb2*lb2u*mmsb2*mmst2*mmt + 2*invdgb2*lt2u*
     mmsb2*mmst2*mmt + (4*invdgb2*lb2u*lt2u*mmsb2*mmst2*mmt)/3. + 12*invdgb2*
     ltu*mmsb2*mmst2*mmt - 2*invdgb2*lb2u*ltu*mmsb2*mmst2*mmt - 2*invdgb2*lt2u*
     ltu*mmsb2*mmst2*mmt - (2*invdgb2*mmsb22*mmst2*zt2)/3. - (2*invdgb2*mmsb22*
     mmt*zt2)/3. - (8*invdgb2*mmsb2*mmst2*mmt*zt2)/3. + (112*mmsb22*mmst2*mmt*
     pow2(invdgb2))/3. - 4*lb2u*mmsb22*mmst2*mmt*pow2(invdgb2) - 4*lt2u*mmsb22*
     mmst2*mmt*pow2(invdgb2) - (8*lb2u*lt2u*mmsb22*mmst2*mmt*pow2(invdgb2))/3.
     - 24*ltu*mmsb22*mmst2*mmt*pow2(invdgb2) + 4*lb2u*ltu*mmsb22*mmst2*mmt*pow2
     (invdgb2) + 4*lt2u*ltu*mmsb22*mmst2*mmt*pow2(invdgb2) + (16*mmsb22*mmst2*
     mmt*zt2*pow2(invdgb2))/3. - (2*invdgb2*mmsb22*mmst2*pow2(lb2u))/3. - (
     invdgb2*mmsb22*mmt*pow2(lb2u))/3. - (invdgb2*mmsb2*mmst2*mmt*pow2(lb2u))/
     3. + (2*mmsb22*mmst2*mmt*pow2(invdgb2)*pow2(lb2u))/3. + (invdgb2*mmsb22*
     mmst2*pow2(lt2u))/3. - (invdgb2*mmsb2*mmst2*mmt*pow2(lt2u))/3. + (2*mmsb22
     *mmst2*mmt*pow2(invdgb2)*pow2(lt2u))/3. - (invdgb2*mmsb22*mmst2*pow2(ltu))
     /3. - (invdgb2*mmsb22*mmt*pow2(ltu))/3. - 2*invdgb2*mmsb2*mmst2*mmt*pow2(
     ltu) + 4*mmsb22*mmst2*mmt*pow2(invdgb2)*pow2(ltu) - (14*invdgb2*mmsb2*pow2
     (mmst2))/3. - 2*invdgb2*lb2u*mmsb2*pow2(mmst2) + 4*invdgb2*lt2u*mmsb2*pow2
     (mmst2) + 2*invdgb2*ltu*mmsb2*pow2(mmst2) + (2*invdgb2*lb2u*ltu*mmsb2*pow2
     (mmst2))/3. - (4*invdgb2*lt2u*ltu*mmsb2*pow2(mmst2))/3. - (14*invdgb2*mmt*
     pow2(mmst2))/3. + 2*invdgb2*lt2u*mmt*pow2(mmst2) + 2*invdgb2*ltu*mmt*pow2(
     mmst2) - (2*invdgb2*lt2u*ltu*mmt*pow2(mmst2))/3. - (2*invdgb2*mmsb2*zt2*
     pow2(mmst2))/3. - (2*invdgb2*mmt*zt2*pow2(mmst2))/3. + (28*mmsb22*pow2(
     invdgb2)*pow2(mmst2))/3. + 4*lb2u*mmsb22*pow2(invdgb2)*pow2(mmst2) - 8*
     lt2u*mmsb22*pow2(invdgb2)*pow2(mmst2) - 4*ltu*mmsb22*pow2(invdgb2)*pow2(
     mmst2) - (4*lb2u*ltu*mmsb22*pow2(invdgb2)*pow2(mmst2))/3. + (8*lt2u*ltu*
     mmsb22*pow2(invdgb2)*pow2(mmst2))/3. + (28*mmsb2*mmt*pow2(invdgb2)*pow2(
     mmst2))/3. - 4*lt2u*mmsb2*mmt*pow2(invdgb2)*pow2(mmst2) - 4*ltu*mmsb2*mmt*
     pow2(invdgb2)*pow2(mmst2) + (4*lt2u*ltu*mmsb2*mmt*pow2(invdgb2)*pow2(mmst2
     ))/3. + (4*mmsb22*zt2*pow2(invdgb2)*pow2(mmst2))/3. + (4*mmsb2*mmt*zt2*
     pow2(invdgb2)*pow2(mmst2))/3. + (invdgb2*mmsb2*pow2(lb2u)*pow2(mmst2))/3.
     - (2*mmsb22*pow2(invdgb2)*pow2(lb2u)*pow2(mmst2))/3. - (2*invdgb2*mmsb2*
     pow2(lt2u)*pow2(mmst2))/3. - (invdgb2*mmt*pow2(lt2u)*pow2(mmst2))/3. + (4*
     mmsb22*pow2(invdgb2)*pow2(lt2u)*pow2(mmst2))/3. + (2*mmsb2*mmt*pow2(
     invdgb2)*pow2(lt2u)*pow2(mmst2))/3. - (invdgb2*mmsb2*pow2(ltu)*pow2(mmst2)
     )/3. - (invdgb2*mmt*pow2(ltu)*pow2(mmst2))/3. + (2*mmsb22*pow2(invdgb2)*
     pow2(ltu)*pow2(mmst2))/3. + (2*mmsb2*mmt*pow2(invdgb2)*pow2(ltu)*pow2(
     mmst2))/3. + (14*invdgb2*pow3(mmsb2))/3. - 2*invdgb2*lb2u*pow3(mmsb2) - 2*
     invdgb2*ltu*pow3(mmsb2) + (2*invdgb2*lb2u*ltu*pow3(mmsb2))/3. + (2*invdgb2
     *zt2*pow3(mmsb2))/3. + (28*mmst2*pow2(invdgb2)*pow3(mmsb2))/3. - 8*lb2u*
     mmst2*pow2(invdgb2)*pow3(mmsb2) + 4*lt2u*mmst2*pow2(invdgb2)*pow3(mmsb2) -
     4*ltu*mmst2*pow2(invdgb2)*pow3(mmsb2) + (8*lb2u*ltu*mmst2*pow2(invdgb2)*
     pow3(mmsb2))/3. - (4*lt2u*ltu*mmst2*pow2(invdgb2)*pow3(mmsb2))/3. + (28*
     mmt*pow2(invdgb2)*pow3(mmsb2))/3. - 4*lb2u*mmt*pow2(invdgb2)*pow3(mmsb2) -
     4*ltu*mmt*pow2(invdgb2)*pow3(mmsb2) + (4*lb2u*ltu*mmt*pow2(invdgb2)*pow3(
     mmsb2))/3. + (4*mmst2*zt2*pow2(invdgb2)*pow3(mmsb2))/3. + (4*mmt*zt2*pow2(
     invdgb2)*pow3(mmsb2))/3. + (invdgb2*pow2(lb2u)*pow3(mmsb2))/3. + (4*mmst2*
     pow2(invdgb2)*pow2(lb2u)*pow3(mmsb2))/3. + (2*mmt*pow2(invdgb2)*pow2(lb2u)
     *pow3(mmsb2))/3. - (2*mmst2*pow2(invdgb2)*pow2(lt2u)*pow3(mmsb2))/3. + (
     invdgb2*pow2(ltu)*pow3(mmsb2))/3. + (2*mmst2*pow2(invdgb2)*pow2(ltu)*pow3(
     mmsb2))/3. + (2*mmt*pow2(invdgb2)*pow2(ltu)*pow3(mmsb2))/3. + s2t*((28*
     mmgl*mmsb22*mmst2*mmt*pow2(invdgb2))/(3.*mgl*mt) - (4*lb2u*mmgl*mmsb22*
     mmst2*mmt*pow2(invdgb2))/(mgl*mt) + (4*lt2u*mmgl*mmsb22*mmst2*mmt*pow2(
     invdgb2))/(mgl*mt) - (4*lb2u*lt2u*mmgl*mmsb22*mmst2*mmt*pow2(invdgb2))/(3.
     *mgl*mt) - (8*ltu*mmgl*mmsb22*mmst2*mmt*pow2(invdgb2))/(mgl*mt) + (8*lb2u*
     ltu*mmgl*mmsb22*mmst2*mmt*pow2(invdgb2))/(3.*mgl*mt) + (4*mmgl*mmsb22*
     mmst2*mmt*zt2*pow2(invdgb2))/(3.*mgl*mt) + (2*mmgl*mmsb22*mmst2*mmt*pow2(
     invdgb2)*pow2(lb2u))/(3.*mgl*mt) - (2*mmgl*mmsb22*mmst2*mmt*pow2(invdgb2)*
     pow2(lt2u))/(3.*mgl*mt) + (4*mmgl*mmsb22*mmst2*mmt*pow2(invdgb2)*pow2(ltu)
     )/(3.*mgl*mt) - (4*lt2u*mmgl*mmsb2*mmt*pow2(invdgb2)*pow2(mmst2))/(mgl*mt)
     + (4*lb2u*lt2u*mmgl*mmsb2*mmt*pow2(invdgb2)*pow2(mmst2))/(3.*mgl*mt) + (4*
     ltu*mmgl*mmsb2*mmt*pow2(invdgb2)*pow2(mmst2))/(mgl*mt) - (4*lb2u*ltu*mmgl*
     mmsb2*mmt*pow2(invdgb2)*pow2(mmst2))/(3.*mgl*mt) + (2*mmgl*mmsb2*mmt*pow2(
     invdgb2)*pow2(lt2u)*pow2(mmst2))/(3.*mgl*mt) - (2*mmgl*mmsb2*mmt*pow2(
     invdgb2)*pow2(ltu)*pow2(mmst2))/(3.*mgl*mt) + (28*mmgl*mmsb22*pow2(invdgb2
     )*pow2(mmt))/(3.*mgl*mt) - (4*lb2u*mmgl*mmsb22*pow2(invdgb2)*pow2(mmt))/(
     mgl*mt) - (4*ltu*mmgl*mmsb22*pow2(invdgb2)*pow2(mmt))/(mgl*mt) + (4*lb2u*
     ltu*mmgl*mmsb22*pow2(invdgb2)*pow2(mmt))/(3.*mgl*mt) + (56*mmgl*mmsb2*
     mmst2*pow2(invdgb2)*pow2(mmt))/(3.*mgl*mt) - (4*lt2u*mmgl*mmsb2*mmst2*pow2
     (invdgb2)*pow2(mmt))/(mgl*mt) - (4*lb2u*lt2u*mmgl*mmsb2*mmst2*pow2(invdgb2
     )*pow2(mmt))/(3.*mgl*mt) - (12*ltu*mmgl*mmsb2*mmst2*pow2(invdgb2)*pow2(mmt
     ))/(mgl*mt) + (4*lb2u*ltu*mmgl*mmsb2*mmst2*pow2(invdgb2)*pow2(mmt))/(3.*
     mgl*mt) + (8*lt2u*ltu*mmgl*mmsb2*mmst2*pow2(invdgb2)*pow2(mmt))/(3.*mgl*mt
     ) + (4*mmgl*mmsb22*zt2*pow2(invdgb2)*pow2(mmt))/(3.*mgl*mt) + (8*mmgl*
     mmsb2*mmst2*zt2*pow2(invdgb2)*pow2(mmt))/(3.*mgl*mt) + (2*mmgl*mmsb22*pow2
     (invdgb2)*pow2(lb2u)*pow2(mmt))/(3.*mgl*mt) + (2*mmgl*mmsb2*mmst2*pow2(
     invdgb2)*pow2(lt2u)*pow2(mmt))/(3.*mgl*mt) + (2*mmgl*mmsb22*pow2(invdgb2)*
     pow2(ltu)*pow2(mmt))/(3.*mgl*mt) + (2*mmgl*mmsb2*mmst2*pow2(invdgb2)*pow2(
     ltu)*pow2(mmt))/(mgl*mt) - (28*mmgl*mmt*pow2(invdgb2)*pow3(mmsb2))/(3.*mgl
     *mt) + (4*lb2u*mmgl*mmt*pow2(invdgb2)*pow3(mmsb2))/(mgl*mt) + (4*ltu*mmgl*
     mmt*pow2(invdgb2)*pow3(mmsb2))/(mgl*mt) - (4*lb2u*ltu*mmgl*mmt*pow2(
     invdgb2)*pow3(mmsb2))/(3.*mgl*mt) - (4*mmgl*mmt*zt2*pow2(invdgb2)*pow3(
     mmsb2))/(3.*mgl*mt) - (2*mmgl*mmt*pow2(invdgb2)*pow2(lb2u)*pow3(mmsb2))/(
     3.*mgl*mt) - (2*mmgl*mmt*pow2(invdgb2)*pow2(ltu)*pow3(mmsb2))/(3.*mgl*mt))
     + (14*invdgb2*pow3(mmst2))/3. - 2*invdgb2*lt2u*pow3(mmst2) - 2*invdgb2*ltu
     *pow3(mmst2) + (2*invdgb2*lt2u*ltu*pow3(mmst2))/3. + (2*invdgb2*zt2*pow3(
     mmst2))/3. - (28*mmsb2*pow2(invdgb2)*pow3(mmst2))/3. + 4*lt2u*mmsb2*pow2(
     invdgb2)*pow3(mmst2) + 4*ltu*mmsb2*pow2(invdgb2)*pow3(mmst2) - (4*lt2u*ltu
     *mmsb2*pow2(invdgb2)*pow3(mmst2))/3. - (4*mmsb2*zt2*pow2(invdgb2)*pow3(
     mmst2))/3. + (invdgb2*pow2(lt2u)*pow3(mmst2))/3. - (2*mmsb2*pow2(invdgb2)*
     pow2(lt2u)*pow3(mmst2))/3. + (invdgb2*pow2(ltu)*pow3(mmst2))/3. - (2*mmsb2
     *pow2(invdgb2)*pow2(ltu)*pow3(mmst2))/3. - (28*pow2(invdgb2)*pow4(mmsb2))/
     3. + 4*lb2u*pow2(invdgb2)*pow4(mmsb2) + 4*ltu*pow2(invdgb2)*pow4(mmsb2) -
     (4*lb2u*ltu*pow2(invdgb2)*pow4(mmsb2))/3. - (4*zt2*pow2(invdgb2)*pow4(
     mmsb2))/3. - (2*pow2(invdgb2)*pow2(lb2u)*pow4(mmsb2))/3. - (2*pow2(invdgb2
     )*pow2(ltu)*pow4(mmsb2))/3.) + DeltaInv(mmt,mmst1,mmgl)*((-56*mmgl2)/3. +
     8*lgu*mmgl2 + 8*ltu*mmgl2 - (8*lgu*ltu*mmgl2)/3. - (56*mmgl*mmsb1)/3. + 8*
     lgu*mmgl*mmsb1 + 8*ltu*mmgl*mmsb1 - (8*lgu*ltu*mmgl*mmsb1)/3. - 28*mmsb12
     + 12*lgu*mmsb12 + 12*ltu*mmsb12 - 4*lgu*ltu*mmsb12 - (56*mmgl*mmsb2)/3. +
     8*lgu*mmgl*mmsb2 + 8*ltu*mmgl*mmsb2 - (8*lgu*ltu*mmgl*mmsb2)/3. - 28*
     mmsb22 + 12*lgu*mmsb22 + 12*ltu*mmsb22 - 4*lgu*ltu*mmsb22 + (56*mmgl*mmst1
     )/3. - 16*lgu*mmgl*mmst1 + 8*lt1u*mmgl*mmst1 - 8*ltu*mmgl*mmst1 + (16*lgu*
     ltu*mmgl*mmst1)/3. - (8*lt1u*ltu*mmgl*mmst1)/3. + (56*mmsb1*mmst1)/3. - 16
     *lgu*mmsb1*mmst1 + 8*lt1u*mmsb1*mmst1 - 8*ltu*mmsb1*mmst1 + (16*lgu*ltu*
     mmsb1*mmst1)/3. - (8*lt1u*ltu*mmsb1*mmst1)/3. - 28*invdgb1*mmsb12*mmst1 +
     24*invdgb1*lgu*mmsb12*mmst1 - 12*invdgb1*lt1u*mmsb12*mmst1 + 12*invdgb1*
     ltu*mmsb12*mmst1 - 8*invdgb1*lgu*ltu*mmsb12*mmst1 + 4*invdgb1*lt1u*ltu*
     mmsb12*mmst1 + (56*mmsb2*mmst1)/3. - 16*lgu*mmsb2*mmst1 + 8*lt1u*mmsb2*
     mmst1 - 8*ltu*mmsb2*mmst1 + (16*lgu*ltu*mmsb2*mmst1)/3. - (8*lt1u*ltu*
     mmsb2*mmst1)/3. - 28*invdgb2*mmsb22*mmst1 + 24*invdgb2*lgu*mmsb22*mmst1 -
     12*invdgb2*lt1u*mmsb22*mmst1 + 12*invdgb2*ltu*mmsb22*mmst1 - 8*invdgb2*lgu
     *ltu*mmsb22*mmst1 + 4*invdgb2*lt1u*ltu*mmsb22*mmst1 + (56*mmgl*mmt)/3. - 8
     *lgu*mmgl*mmt - 8*ltu*mmgl*mmt + (8*lgu*ltu*mmgl*mmt)/3. + (56*mmsb1*mmt)/
     3. - 8*lgu*mmsb1*mmt - 8*ltu*mmsb1*mmt + (8*lgu*ltu*mmsb1*mmt)/3. - 28*
     invdgb1*mmsb12*mmt + 12*invdgb1*lgu*mmsb12*mmt + 12*invdgb1*ltu*mmsb12*mmt
      - 4*invdgb1*lgu*ltu*mmsb12*mmt + (56*mmsb2*mmt)/3. - 8*lgu*mmsb2*mmt - 8*
     ltu*mmsb2*mmt + (8*lgu*ltu*mmsb2*mmt)/3. - 28*invdgb2*mmsb22*mmt + 12*
     invdgb2*lgu*mmsb22*mmt + 12*invdgb2*ltu*mmsb22*mmt - 4*invdgb2*lgu*ltu*
     mmsb22*mmt + (224*mmst1*mmt)/3. - 8*lgu*mmst1*mmt - 8*lt1u*mmst1*mmt - (16
     *lgu*lt1u*mmst1*mmt)/3. - 48*ltu*mmst1*mmt + 8*lgu*ltu*mmst1*mmt + 8*lt1u*
     ltu*mmst1*mmt - (224*invdgb1*mmsb1*mmst1*mmt)/3. + 8*invdgb1*lgu*mmsb1*
     mmst1*mmt + 8*invdgb1*lt1u*mmsb1*mmst1*mmt + (16*invdgb1*lgu*lt1u*mmsb1*
     mmst1*mmt)/3. + 48*invdgb1*ltu*mmsb1*mmst1*mmt - 8*invdgb1*lgu*ltu*mmsb1*
     mmst1*mmt - 8*invdgb1*lt1u*ltu*mmsb1*mmst1*mmt - (224*invdgb2*mmsb2*mmst1*
     mmt)/3. + 8*invdgb2*lgu*mmsb2*mmst1*mmt + 8*invdgb2*lt1u*mmsb2*mmst1*mmt +
     (16*invdgb2*lgu*lt1u*mmsb2*mmst1*mmt)/3. + 48*invdgb2*ltu*mmsb2*mmst1*mmt
     - 8*invdgb2*lgu*ltu*mmsb2*mmst1*mmt - 8*invdgb2*lt1u*ltu*mmsb2*mmst1*mmt -
     (8*mmgl2*zt2)/3. - (8*mmgl*mmsb1*zt2)/3. - 4*mmsb12*zt2 - (8*mmgl*mmsb2*
     zt2)/3. - 4*mmsb22*zt2 + (8*mmgl*mmst1*zt2)/3. + (8*mmsb1*mmst1*zt2)/3. -
     4*invdgb1*mmsb12*mmst1*zt2 + (8*mmsb2*mmst1*zt2)/3. - 4*invdgb2*mmsb22*
     mmst1*zt2 + (8*mmgl*mmt*zt2)/3. + (8*mmsb1*mmt*zt2)/3. - 4*invdgb1*mmsb12*
     mmt*zt2 + (8*mmsb2*mmt*zt2)/3. - 4*invdgb2*mmsb22*mmt*zt2 + (32*mmst1*mmt*
     zt2)/3. - (32*invdgb1*mmsb1*mmst1*mmt*zt2)/3. - (32*invdgb2*mmsb2*mmst1*
     mmt*zt2)/3. + (112*mmsb12*mmst1*mmt*pow2(invdgb1))/3. - 4*lgu*mmsb12*mmst1
     *mmt*pow2(invdgb1) - 4*lt1u*mmsb12*mmst1*mmt*pow2(invdgb1) - (8*lgu*lt1u*
     mmsb12*mmst1*mmt*pow2(invdgb1))/3. - 24*ltu*mmsb12*mmst1*mmt*pow2(invdgb1)
     + 4*lgu*ltu*mmsb12*mmst1*mmt*pow2(invdgb1) + 4*lt1u*ltu*mmsb12*mmst1*mmt*
     pow2(invdgb1) + (16*mmsb12*mmst1*mmt*zt2*pow2(invdgb1))/3. + (112*mmsb22*
     mmst1*mmt*pow2(invdgb2))/3. - 4*lgu*mmsb22*mmst1*mmt*pow2(invdgb2) - 4*
     lt1u*mmsb22*mmst1*mmt*pow2(invdgb2) - (8*lgu*lt1u*mmsb22*mmst1*mmt*pow2(
     invdgb2))/3. - 24*ltu*mmsb22*mmst1*mmt*pow2(invdgb2) + 4*lgu*ltu*mmsb22*
     mmst1*mmt*pow2(invdgb2) + 4*lt1u*ltu*mmsb22*mmst1*mmt*pow2(invdgb2) + (16*
     mmsb22*mmst1*mmt*zt2*pow2(invdgb2))/3. - (4*mmgl2*pow2(lgu))/3. - (4*mmgl*
     mmsb1*pow2(lgu))/3. - 2*mmsb12*pow2(lgu) - (4*mmgl*mmsb2*pow2(lgu))/3. - 2
     *mmsb22*pow2(lgu) + (8*mmgl*mmst1*pow2(lgu))/3. + (8*mmsb1*mmst1*pow2(lgu)
     )/3. - 4*invdgb1*mmsb12*mmst1*pow2(lgu) + (8*mmsb2*mmst1*pow2(lgu))/3. - 4
     *invdgb2*mmsb22*mmst1*pow2(lgu) + (4*mmgl*mmt*pow2(lgu))/3. + (4*mmsb1*mmt
     *pow2(lgu))/3. - 2*invdgb1*mmsb12*mmt*pow2(lgu) + (4*mmsb2*mmt*pow2(lgu))/
     3. - 2*invdgb2*mmsb22*mmt*pow2(lgu) + (4*mmst1*mmt*pow2(lgu))/3. - (4*
     invdgb1*mmsb1*mmst1*mmt*pow2(lgu))/3. - (4*invdgb2*mmsb2*mmst1*mmt*pow2(
     lgu))/3. + (2*mmsb12*mmst1*mmt*pow2(invdgb1)*pow2(lgu))/3. + (2*mmsb22*
     mmst1*mmt*pow2(invdgb2)*pow2(lgu))/3. - (4*mmgl*mmst1*pow2(lt1u))/3. - (4*
     mmsb1*mmst1*pow2(lt1u))/3. + 2*invdgb1*mmsb12*mmst1*pow2(lt1u) - (4*mmsb2*
     mmst1*pow2(lt1u))/3. + 2*invdgb2*mmsb22*mmst1*pow2(lt1u) + (4*mmst1*mmt*
     pow2(lt1u))/3. - (4*invdgb1*mmsb1*mmst1*mmt*pow2(lt1u))/3. - (4*invdgb2*
     mmsb2*mmst1*mmt*pow2(lt1u))/3. + (2*mmsb12*mmst1*mmt*pow2(invdgb1)*pow2(
     lt1u))/3. + (2*mmsb22*mmst1*mmt*pow2(invdgb2)*pow2(lt1u))/3. - (4*mmgl2*
     pow2(ltu))/3. - (4*mmgl*mmsb1*pow2(ltu))/3. - 2*mmsb12*pow2(ltu) - (4*mmgl
     *mmsb2*pow2(ltu))/3. - 2*mmsb22*pow2(ltu) + (4*mmgl*mmst1*pow2(ltu))/3. +
     (4*mmsb1*mmst1*pow2(ltu))/3. - 2*invdgb1*mmsb12*mmst1*pow2(ltu) + (4*mmsb2
     *mmst1*pow2(ltu))/3. - 2*invdgb2*mmsb22*mmst1*pow2(ltu) + (4*mmgl*mmt*pow2
     (ltu))/3. + (4*mmsb1*mmt*pow2(ltu))/3. - 2*invdgb1*mmsb12*mmt*pow2(ltu) +
     (4*mmsb2*mmt*pow2(ltu))/3. - 2*invdgb2*mmsb22*mmt*pow2(ltu) + 8*mmst1*mmt*
     pow2(ltu) - 8*invdgb1*mmsb1*mmst1*mmt*pow2(ltu) - 8*invdgb2*mmsb2*mmst1*
     mmt*pow2(ltu) + 4*mmsb12*mmst1*mmt*pow2(invdgb1)*pow2(ltu) + 4*mmsb22*
     mmst1*mmt*pow2(invdgb2)*pow2(ltu) + (56*pow2(mmst1))/3. + 8*lgu*pow2(mmst1
     ) - 16*lt1u*pow2(mmst1) - 8*ltu*pow2(mmst1) - (8*lgu*ltu*pow2(mmst1))/3. +
     (16*lt1u*ltu*pow2(mmst1))/3. - (56*invdgb1*mmsb1*pow2(mmst1))/3. - 8*
     invdgb1*lgu*mmsb1*pow2(mmst1) + 16*invdgb1*lt1u*mmsb1*pow2(mmst1) + 8*
     invdgb1*ltu*mmsb1*pow2(mmst1) + (8*invdgb1*lgu*ltu*mmsb1*pow2(mmst1))/3. -
     (16*invdgb1*lt1u*ltu*mmsb1*pow2(mmst1))/3. - (56*invdgb2*mmsb2*pow2(mmst1)
     )/3. - 8*invdgb2*lgu*mmsb2*pow2(mmst1) + 16*invdgb2*lt1u*mmsb2*pow2(mmst1)
     + 8*invdgb2*ltu*mmsb2*pow2(mmst1) + (8*invdgb2*lgu*ltu*mmsb2*pow2(mmst1))/
     3. - (16*invdgb2*lt1u*ltu*mmsb2*pow2(mmst1))/3. - (28*invdgb1*mmt*pow2(
     mmst1))/3. - (28*invdgb2*mmt*pow2(mmst1))/3. + 4*invdgb1*lt1u*mmt*pow2(
     mmst1) + 4*invdgb2*lt1u*mmt*pow2(mmst1) + 4*invdgb1*ltu*mmt*pow2(mmst1) +
     4*invdgb2*ltu*mmt*pow2(mmst1) - (4*invdgb1*lt1u*ltu*mmt*pow2(mmst1))/3. -
     (4*invdgb2*lt1u*ltu*mmt*pow2(mmst1))/3. + (8*zt2*pow2(mmst1))/3. - (8*
     invdgb1*mmsb1*zt2*pow2(mmst1))/3. - (8*invdgb2*mmsb2*zt2*pow2(mmst1))/3. -
     (4*invdgb1*mmt*zt2*pow2(mmst1))/3. - (4*invdgb2*mmt*zt2*pow2(mmst1))/3. +
     (28*mmsb12*pow2(invdgb1)*pow2(mmst1))/3. + 4*lgu*mmsb12*pow2(invdgb1)*pow2
     (mmst1) - 8*lt1u*mmsb12*pow2(invdgb1)*pow2(mmst1) - 4*ltu*mmsb12*pow2(
     invdgb1)*pow2(mmst1) - (4*lgu*ltu*mmsb12*pow2(invdgb1)*pow2(mmst1))/3. + (
     8*lt1u*ltu*mmsb12*pow2(invdgb1)*pow2(mmst1))/3. + (28*mmsb1*mmt*pow2(
     invdgb1)*pow2(mmst1))/3. - 4*lt1u*mmsb1*mmt*pow2(invdgb1)*pow2(mmst1) - 4*
     ltu*mmsb1*mmt*pow2(invdgb1)*pow2(mmst1) + (4*lt1u*ltu*mmsb1*mmt*pow2(
     invdgb1)*pow2(mmst1))/3. + (4*mmsb12*zt2*pow2(invdgb1)*pow2(mmst1))/3. + (
     4*mmsb1*mmt*zt2*pow2(invdgb1)*pow2(mmst1))/3. + (28*mmsb22*pow2(invdgb2)*
     pow2(mmst1))/3. + 4*lgu*mmsb22*pow2(invdgb2)*pow2(mmst1) - 8*lt1u*mmsb22*
     pow2(invdgb2)*pow2(mmst1) - 4*ltu*mmsb22*pow2(invdgb2)*pow2(mmst1) - (4*
     lgu*ltu*mmsb22*pow2(invdgb2)*pow2(mmst1))/3. + (8*lt1u*ltu*mmsb22*pow2(
     invdgb2)*pow2(mmst1))/3. + (28*mmsb2*mmt*pow2(invdgb2)*pow2(mmst1))/3. - 4
     *lt1u*mmsb2*mmt*pow2(invdgb2)*pow2(mmst1) - 4*ltu*mmsb2*mmt*pow2(invdgb2)*
     pow2(mmst1) + (4*lt1u*ltu*mmsb2*mmt*pow2(invdgb2)*pow2(mmst1))/3. + (4*
     mmsb22*zt2*pow2(invdgb2)*pow2(mmst1))/3. + (4*mmsb2*mmt*zt2*pow2(invdgb2)*
     pow2(mmst1))/3. - (4*pow2(lgu)*pow2(mmst1))/3. + (4*invdgb1*mmsb1*pow2(lgu
     )*pow2(mmst1))/3. + (4*invdgb2*mmsb2*pow2(lgu)*pow2(mmst1))/3. - (2*mmsb12
     *pow2(invdgb1)*pow2(lgu)*pow2(mmst1))/3. - (2*mmsb22*pow2(invdgb2)*pow2(
     lgu)*pow2(mmst1))/3. + (8*pow2(lt1u)*pow2(mmst1))/3. - (8*invdgb1*mmsb1*
     pow2(lt1u)*pow2(mmst1))/3. - (8*invdgb2*mmsb2*pow2(lt1u)*pow2(mmst1))/3. -
     (2*invdgb1*mmt*pow2(lt1u)*pow2(mmst1))/3. - (2*invdgb2*mmt*pow2(lt1u)*pow2
     (mmst1))/3. + (4*mmsb12*pow2(invdgb1)*pow2(lt1u)*pow2(mmst1))/3. + (2*
     mmsb1*mmt*pow2(invdgb1)*pow2(lt1u)*pow2(mmst1))/3. + (4*mmsb22*pow2(
     invdgb2)*pow2(lt1u)*pow2(mmst1))/3. + (2*mmsb2*mmt*pow2(invdgb2)*pow2(lt1u
     )*pow2(mmst1))/3. + (4*pow2(ltu)*pow2(mmst1))/3. - (4*invdgb1*mmsb1*pow2(
     ltu)*pow2(mmst1))/3. - (4*invdgb2*mmsb2*pow2(ltu)*pow2(mmst1))/3. - (2*
     invdgb1*mmt*pow2(ltu)*pow2(mmst1))/3. - (2*invdgb2*mmt*pow2(ltu)*pow2(
     mmst1))/3. + (2*mmsb12*pow2(invdgb1)*pow2(ltu)*pow2(mmst1))/3. + (2*mmsb1*
     mmt*pow2(invdgb1)*pow2(ltu)*pow2(mmst1))/3. + (2*mmsb22*pow2(invdgb2)*pow2
     (ltu)*pow2(mmst1))/3. + (2*mmsb2*mmt*pow2(invdgb2)*pow2(ltu)*pow2(mmst1))/
     3. + (112*invdgb1*pow3(mmsb1))/3. - 16*invdgb1*lgu*pow3(mmsb1) - 16*
     invdgb1*ltu*pow3(mmsb1) + (16*invdgb1*lgu*ltu*pow3(mmsb1))/3. + (16*
     invdgb1*zt2*pow3(mmsb1))/3. + (28*mmst1*pow2(invdgb1)*pow3(mmsb1))/3. - 8*
     lgu*mmst1*pow2(invdgb1)*pow3(mmsb1) + 4*lt1u*mmst1*pow2(invdgb1)*pow3(
     mmsb1) - 4*ltu*mmst1*pow2(invdgb1)*pow3(mmsb1) + (8*lgu*ltu*mmst1*pow2(
     invdgb1)*pow3(mmsb1))/3. - (4*lt1u*ltu*mmst1*pow2(invdgb1)*pow3(mmsb1))/3.
      + (28*mmt*pow2(invdgb1)*pow3(mmsb1))/3. - 4*lgu*mmt*pow2(invdgb1)*pow3(
     mmsb1) - 4*ltu*mmt*pow2(invdgb1)*pow3(mmsb1) + (4*lgu*ltu*mmt*pow2(invdgb1
     )*pow3(mmsb1))/3. + (4*mmst1*zt2*pow2(invdgb1)*pow3(mmsb1))/3. + (4*mmt*
     zt2*pow2(invdgb1)*pow3(mmsb1))/3. + (8*invdgb1*pow2(lgu)*pow3(mmsb1))/3. +
     (4*mmst1*pow2(invdgb1)*pow2(lgu)*pow3(mmsb1))/3. + (2*mmt*pow2(invdgb1)*
     pow2(lgu)*pow3(mmsb1))/3. - (2*mmst1*pow2(invdgb1)*pow2(lt1u)*pow3(mmsb1))
     /3. + (8*invdgb1*pow2(ltu)*pow3(mmsb1))/3. + (2*mmst1*pow2(invdgb1)*pow2(
     ltu)*pow3(mmsb1))/3. + (2*mmt*pow2(invdgb1)*pow2(ltu)*pow3(mmsb1))/3. + (
     112*invdgb2*pow3(mmsb2))/3. - 16*invdgb2*lgu*pow3(mmsb2) - 16*invdgb2*ltu*
     pow3(mmsb2) + (16*invdgb2*lgu*ltu*pow3(mmsb2))/3. + (16*invdgb2*zt2*pow3(
     mmsb2))/3. + (28*mmst1*pow2(invdgb2)*pow3(mmsb2))/3. - 8*lgu*mmst1*pow2(
     invdgb2)*pow3(mmsb2) + 4*lt1u*mmst1*pow2(invdgb2)*pow3(mmsb2) - 4*ltu*
     mmst1*pow2(invdgb2)*pow3(mmsb2) + (8*lgu*ltu*mmst1*pow2(invdgb2)*pow3(
     mmsb2))/3. - (4*lt1u*ltu*mmst1*pow2(invdgb2)*pow3(mmsb2))/3. + (28*mmt*
     pow2(invdgb2)*pow3(mmsb2))/3. - 4*lgu*mmt*pow2(invdgb2)*pow3(mmsb2) - 4*
     ltu*mmt*pow2(invdgb2)*pow3(mmsb2) + (4*lgu*ltu*mmt*pow2(invdgb2)*pow3(
     mmsb2))/3. + (4*mmst1*zt2*pow2(invdgb2)*pow3(mmsb2))/3. + (4*mmt*zt2*pow2(
     invdgb2)*pow3(mmsb2))/3. + (8*invdgb2*pow2(lgu)*pow3(mmsb2))/3. + (4*mmst1
     *pow2(invdgb2)*pow2(lgu)*pow3(mmsb2))/3. + (2*mmt*pow2(invdgb2)*pow2(lgu)*
     pow3(mmsb2))/3. - (2*mmst1*pow2(invdgb2)*pow2(lt1u)*pow3(mmsb2))/3. + (8*
     invdgb2*pow2(ltu)*pow3(mmsb2))/3. + (2*mmst1*pow2(invdgb2)*pow2(ltu)*pow3(
     mmsb2))/3. + (2*mmt*pow2(invdgb2)*pow2(ltu)*pow3(mmsb2))/3. + s2t*((56*
     mmgl2*mmt)/(3.*mgl*mt) - (8*lgu*mmgl2*mmt)/(mgl*mt) - (8*ltu*mmgl2*mmt)/(
     mgl*mt) + (8*lgu*ltu*mmgl2*mmt)/(3.*mgl*mt) + (56*mmgl*mmsb1*mmt)/(3.*mgl*
     mt) - (8*lgu*mmgl*mmsb1*mmt)/(mgl*mt) - (8*ltu*mmgl*mmsb1*mmt)/(mgl*mt) +
     (8*lgu*ltu*mmgl*mmsb1*mmt)/(3.*mgl*mt) - (28*invdgb1*mmgl*mmsb12*mmt)/(mgl
     *mt) + (12*invdgb1*lgu*mmgl*mmsb12*mmt)/(mgl*mt) + (12*invdgb1*ltu*mmgl*
     mmsb12*mmt)/(mgl*mt) - (4*invdgb1*lgu*ltu*mmgl*mmsb12*mmt)/(mgl*mt) + (56*
     mmgl*mmsb2*mmt)/(3.*mgl*mt) - (8*lgu*mmgl*mmsb2*mmt)/(mgl*mt) - (8*ltu*
     mmgl*mmsb2*mmt)/(mgl*mt) + (8*lgu*ltu*mmgl*mmsb2*mmt)/(3.*mgl*mt) - (28*
     invdgb2*mmgl*mmsb22*mmt)/(mgl*mt) + (12*invdgb2*lgu*mmgl*mmsb22*mmt)/(mgl*
     mt) + (12*invdgb2*ltu*mmgl*mmsb22*mmt)/(mgl*mt) - (4*invdgb2*lgu*ltu*mmgl*
     mmsb22*mmt)/(mgl*mt) - (56*mmgl*mmst1*mmt)/(3.*mgl*mt) + (8*lgu*mmgl*mmst1
     *mmt)/(mgl*mt) - (8*lt1u*mmgl*mmst1*mmt)/(mgl*mt) + (8*lgu*lt1u*mmgl*mmst1
     *mmt)/(3.*mgl*mt) + (16*ltu*mmgl*mmst1*mmt)/(mgl*mt) - (16*lgu*ltu*mmgl*
     mmst1*mmt)/(3.*mgl*mt) + (56*invdgb1*mmgl*mmsb1*mmst1*mmt)/(3.*mgl*mt) - (
     8*invdgb1*lgu*mmgl*mmsb1*mmst1*mmt)/(mgl*mt) + (8*invdgb1*lt1u*mmgl*mmsb1*
     mmst1*mmt)/(mgl*mt) - (8*invdgb1*lgu*lt1u*mmgl*mmsb1*mmst1*mmt)/(3.*mgl*mt
     ) - (16*invdgb1*ltu*mmgl*mmsb1*mmst1*mmt)/(mgl*mt) + (16*invdgb1*lgu*ltu*
     mmgl*mmsb1*mmst1*mmt)/(3.*mgl*mt) + (56*invdgb2*mmgl*mmsb2*mmst1*mmt)/(3.*
     mgl*mt) - (8*invdgb2*lgu*mmgl*mmsb2*mmst1*mmt)/(mgl*mt) + (8*invdgb2*lt1u*
     mmgl*mmsb2*mmst1*mmt)/(mgl*mt) - (8*invdgb2*lgu*lt1u*mmgl*mmsb2*mmst1*mmt)
     /(3.*mgl*mt) - (16*invdgb2*ltu*mmgl*mmsb2*mmst1*mmt)/(mgl*mt) + (16*
     invdgb2*lgu*ltu*mmgl*mmsb2*mmst1*mmt)/(3.*mgl*mt) + (8*mmgl2*mmt*zt2)/(3.*
     mgl*mt) + (8*mmgl*mmsb1*mmt*zt2)/(3.*mgl*mt) - (4*invdgb1*mmgl*mmsb12*mmt*
     zt2)/(mgl*mt) + (8*mmgl*mmsb2*mmt*zt2)/(3.*mgl*mt) - (4*invdgb2*mmgl*
     mmsb22*mmt*zt2)/(mgl*mt) - (8*mmgl*mmst1*mmt*zt2)/(3.*mgl*mt) + (8*invdgb1
     *mmgl*mmsb1*mmst1*mmt*zt2)/(3.*mgl*mt) + (8*invdgb2*mmgl*mmsb2*mmst1*mmt*
     zt2)/(3.*mgl*mt) - (28*mmgl*mmsb12*mmst1*mmt*pow2(invdgb1))/(3.*mgl*mt) +
     (4*lgu*mmgl*mmsb12*mmst1*mmt*pow2(invdgb1))/(mgl*mt) - (4*lt1u*mmgl*mmsb12
     *mmst1*mmt*pow2(invdgb1))/(mgl*mt) + (4*lgu*lt1u*mmgl*mmsb12*mmst1*mmt*
     pow2(invdgb1))/(3.*mgl*mt) + (8*ltu*mmgl*mmsb12*mmst1*mmt*pow2(invdgb1))/(
     mgl*mt) - (8*lgu*ltu*mmgl*mmsb12*mmst1*mmt*pow2(invdgb1))/(3.*mgl*mt) - (4
     *mmgl*mmsb12*mmst1*mmt*zt2*pow2(invdgb1))/(3.*mgl*mt) - (28*mmgl*mmsb22*
     mmst1*mmt*pow2(invdgb2))/(3.*mgl*mt) + (4*lgu*mmgl*mmsb22*mmst1*mmt*pow2(
     invdgb2))/(mgl*mt) - (4*lt1u*mmgl*mmsb22*mmst1*mmt*pow2(invdgb2))/(mgl*mt)
     + (4*lgu*lt1u*mmgl*mmsb22*mmst1*mmt*pow2(invdgb2))/(3.*mgl*mt) + (8*ltu*
     mmgl*mmsb22*mmst1*mmt*pow2(invdgb2))/(mgl*mt) - (8*lgu*ltu*mmgl*mmsb22*
     mmst1*mmt*pow2(invdgb2))/(3.*mgl*mt) - (4*mmgl*mmsb22*mmst1*mmt*zt2*pow2(
     invdgb2))/(3.*mgl*mt) + (4*mmgl2*mmt*pow2(lgu))/(3.*mgl*mt) + (4*mmgl*
     mmsb1*mmt*pow2(lgu))/(3.*mgl*mt) - (2*invdgb1*mmgl*mmsb12*mmt*pow2(lgu))/(
     mgl*mt) + (4*mmgl*mmsb2*mmt*pow2(lgu))/(3.*mgl*mt) - (2*invdgb2*mmgl*
     mmsb22*mmt*pow2(lgu))/(mgl*mt) - (4*mmgl*mmst1*mmt*pow2(lgu))/(3.*mgl*mt)
     + (4*invdgb1*mmgl*mmsb1*mmst1*mmt*pow2(lgu))/(3.*mgl*mt) + (4*invdgb2*mmgl
     *mmsb2*mmst1*mmt*pow2(lgu))/(3.*mgl*mt) - (2*mmgl*mmsb12*mmst1*mmt*pow2(
     invdgb1)*pow2(lgu))/(3.*mgl*mt) - (2*mmgl*mmsb22*mmst1*mmt*pow2(invdgb2)*
     pow2(lgu))/(3.*mgl*mt) + (4*mmgl*mmst1*mmt*pow2(lt1u))/(3.*mgl*mt) - (4*
     invdgb1*mmgl*mmsb1*mmst1*mmt*pow2(lt1u))/(3.*mgl*mt) - (4*invdgb2*mmgl*
     mmsb2*mmst1*mmt*pow2(lt1u))/(3.*mgl*mt) + (2*mmgl*mmsb12*mmst1*mmt*pow2(
     invdgb1)*pow2(lt1u))/(3.*mgl*mt) + (2*mmgl*mmsb22*mmst1*mmt*pow2(invdgb2)*
     pow2(lt1u))/(3.*mgl*mt) + (4*mmgl2*mmt*pow2(ltu))/(3.*mgl*mt) + (4*mmgl*
     mmsb1*mmt*pow2(ltu))/(3.*mgl*mt) - (2*invdgb1*mmgl*mmsb12*mmt*pow2(ltu))/(
     mgl*mt) + (4*mmgl*mmsb2*mmt*pow2(ltu))/(3.*mgl*mt) - (2*invdgb2*mmgl*
     mmsb22*mmt*pow2(ltu))/(mgl*mt) - (8*mmgl*mmst1*mmt*pow2(ltu))/(3.*mgl*mt)
     + (8*invdgb1*mmgl*mmsb1*mmst1*mmt*pow2(ltu))/(3.*mgl*mt) + (8*invdgb2*mmgl
     *mmsb2*mmst1*mmt*pow2(ltu))/(3.*mgl*mt) - (4*mmgl*mmsb12*mmst1*mmt*pow2(
     invdgb1)*pow2(ltu))/(3.*mgl*mt) - (4*mmgl*mmsb22*mmst1*mmt*pow2(invdgb2)*
     pow2(ltu))/(3.*mgl*mt) - (4*invdgb1*lt1u*mmgl*mmt*pow2(mmst1))/(mgl*mt) -
     (4*invdgb2*lt1u*mmgl*mmt*pow2(mmst1))/(mgl*mt) + (4*invdgb1*lgu*lt1u*mmgl*
     mmt*pow2(mmst1))/(3.*mgl*mt) + (4*invdgb2*lgu*lt1u*mmgl*mmt*pow2(mmst1))/(
     3.*mgl*mt) + (4*invdgb1*ltu*mmgl*mmt*pow2(mmst1))/(mgl*mt) + (4*invdgb2*
     ltu*mmgl*mmt*pow2(mmst1))/(mgl*mt) - (4*invdgb1*lgu*ltu*mmgl*mmt*pow2(
     mmst1))/(3.*mgl*mt) - (4*invdgb2*lgu*ltu*mmgl*mmt*pow2(mmst1))/(3.*mgl*mt)
     + (4*lt1u*mmgl*mmsb1*mmt*pow2(invdgb1)*pow2(mmst1))/(mgl*mt) - (4*lgu*lt1u
     *mmgl*mmsb1*mmt*pow2(invdgb1)*pow2(mmst1))/(3.*mgl*mt) - (4*ltu*mmgl*mmsb1
     *mmt*pow2(invdgb1)*pow2(mmst1))/(mgl*mt) + (4*lgu*ltu*mmgl*mmsb1*mmt*pow2(
     invdgb1)*pow2(mmst1))/(3.*mgl*mt) + (4*lt1u*mmgl*mmsb2*mmt*pow2(invdgb2)*
     pow2(mmst1))/(mgl*mt) - (4*lgu*lt1u*mmgl*mmsb2*mmt*pow2(invdgb2)*pow2(
     mmst1))/(3.*mgl*mt) - (4*ltu*mmgl*mmsb2*mmt*pow2(invdgb2)*pow2(mmst1))/(
     mgl*mt) + (4*lgu*ltu*mmgl*mmsb2*mmt*pow2(invdgb2)*pow2(mmst1))/(3.*mgl*mt)
     + (2*invdgb1*mmgl*mmt*pow2(lt1u)*pow2(mmst1))/(3.*mgl*mt) + (2*invdgb2*
     mmgl*mmt*pow2(lt1u)*pow2(mmst1))/(3.*mgl*mt) - (2*mmgl*mmsb1*mmt*pow2(
     invdgb1)*pow2(lt1u)*pow2(mmst1))/(3.*mgl*mt) - (2*mmgl*mmsb2*mmt*pow2(
     invdgb2)*pow2(lt1u)*pow2(mmst1))/(3.*mgl*mt) - (2*invdgb1*mmgl*mmt*pow2(
     ltu)*pow2(mmst1))/(3.*mgl*mt) - (2*invdgb2*mmgl*mmt*pow2(ltu)*pow2(mmst1))
     /(3.*mgl*mt) + (2*mmgl*mmsb1*mmt*pow2(invdgb1)*pow2(ltu)*pow2(mmst1))/(3.*
     mgl*mt) + (2*mmgl*mmsb2*mmt*pow2(invdgb2)*pow2(ltu)*pow2(mmst1))/(3.*mgl*
     mt) - (56*mmgl*pow2(mmt))/(3.*mgl*mt) + (8*lgu*mmgl*pow2(mmt))/(mgl*mt) +
     (8*ltu*mmgl*pow2(mmt))/(mgl*mt) - (8*lgu*ltu*mmgl*pow2(mmt))/(3.*mgl*mt) +
     (56*invdgb1*mmgl*mmsb1*pow2(mmt))/(3.*mgl*mt) - (8*invdgb1*lgu*mmgl*mmsb1*
     pow2(mmt))/(mgl*mt) - (8*invdgb1*ltu*mmgl*mmsb1*pow2(mmt))/(mgl*mt) + (8*
     invdgb1*lgu*ltu*mmgl*mmsb1*pow2(mmt))/(3.*mgl*mt) + (56*invdgb2*mmgl*mmsb2
     *pow2(mmt))/(3.*mgl*mt) - (8*invdgb2*lgu*mmgl*mmsb2*pow2(mmt))/(mgl*mt) -
     (8*invdgb2*ltu*mmgl*mmsb2*pow2(mmt))/(mgl*mt) + (8*invdgb2*lgu*ltu*mmgl*
     mmsb2*pow2(mmt))/(3.*mgl*mt) + (56*invdgb1*mmgl*mmst1*pow2(mmt))/(3.*mgl*
     mt) + (56*invdgb2*mmgl*mmst1*pow2(mmt))/(3.*mgl*mt) - (4*invdgb1*lt1u*mmgl
     *mmst1*pow2(mmt))/(mgl*mt) - (4*invdgb2*lt1u*mmgl*mmst1*pow2(mmt))/(mgl*mt
     ) - (4*invdgb1*lgu*lt1u*mmgl*mmst1*pow2(mmt))/(3.*mgl*mt) - (4*invdgb2*lgu
     *lt1u*mmgl*mmst1*pow2(mmt))/(3.*mgl*mt) - (12*invdgb1*ltu*mmgl*mmst1*pow2(
     mmt))/(mgl*mt) - (12*invdgb2*ltu*mmgl*mmst1*pow2(mmt))/(mgl*mt) + (4*
     invdgb1*lgu*ltu*mmgl*mmst1*pow2(mmt))/(3.*mgl*mt) + (4*invdgb2*lgu*ltu*
     mmgl*mmst1*pow2(mmt))/(3.*mgl*mt) + (8*invdgb1*lt1u*ltu*mmgl*mmst1*pow2(
     mmt))/(3.*mgl*mt) + (8*invdgb2*lt1u*ltu*mmgl*mmst1*pow2(mmt))/(3.*mgl*mt)
     - (8*mmgl*zt2*pow2(mmt))/(3.*mgl*mt) + (8*invdgb1*mmgl*mmsb1*zt2*pow2(mmt)
     )/(3.*mgl*mt) + (8*invdgb2*mmgl*mmsb2*zt2*pow2(mmt))/(3.*mgl*mt) + (8*
     invdgb1*mmgl*mmst1*zt2*pow2(mmt))/(3.*mgl*mt) + (8*invdgb2*mmgl*mmst1*zt2*
     pow2(mmt))/(3.*mgl*mt) - (28*mmgl*mmsb12*pow2(invdgb1)*pow2(mmt))/(3.*mgl*
     mt) + (4*lgu*mmgl*mmsb12*pow2(invdgb1)*pow2(mmt))/(mgl*mt) + (4*ltu*mmgl*
     mmsb12*pow2(invdgb1)*pow2(mmt))/(mgl*mt) - (4*lgu*ltu*mmgl*mmsb12*pow2(
     invdgb1)*pow2(mmt))/(3.*mgl*mt) - (56*mmgl*mmsb1*mmst1*pow2(invdgb1)*pow2(
     mmt))/(3.*mgl*mt) + (4*lt1u*mmgl*mmsb1*mmst1*pow2(invdgb1)*pow2(mmt))/(mgl
     *mt) + (4*lgu*lt1u*mmgl*mmsb1*mmst1*pow2(invdgb1)*pow2(mmt))/(3.*mgl*mt) +
     (12*ltu*mmgl*mmsb1*mmst1*pow2(invdgb1)*pow2(mmt))/(mgl*mt) - (4*lgu*ltu*
     mmgl*mmsb1*mmst1*pow2(invdgb1)*pow2(mmt))/(3.*mgl*mt) - (8*lt1u*ltu*mmgl*
     mmsb1*mmst1*pow2(invdgb1)*pow2(mmt))/(3.*mgl*mt) - (4*mmgl*mmsb12*zt2*pow2
     (invdgb1)*pow2(mmt))/(3.*mgl*mt) - (8*mmgl*mmsb1*mmst1*zt2*pow2(invdgb1)*
     pow2(mmt))/(3.*mgl*mt) - (28*mmgl*mmsb22*pow2(invdgb2)*pow2(mmt))/(3.*mgl*
     mt) + (4*lgu*mmgl*mmsb22*pow2(invdgb2)*pow2(mmt))/(mgl*mt) + (4*ltu*mmgl*
     mmsb22*pow2(invdgb2)*pow2(mmt))/(mgl*mt) - (4*lgu*ltu*mmgl*mmsb22*pow2(
     invdgb2)*pow2(mmt))/(3.*mgl*mt) - (56*mmgl*mmsb2*mmst1*pow2(invdgb2)*pow2(
     mmt))/(3.*mgl*mt) + (4*lt1u*mmgl*mmsb2*mmst1*pow2(invdgb2)*pow2(mmt))/(mgl
     *mt) + (4*lgu*lt1u*mmgl*mmsb2*mmst1*pow2(invdgb2)*pow2(mmt))/(3.*mgl*mt) +
     (12*ltu*mmgl*mmsb2*mmst1*pow2(invdgb2)*pow2(mmt))/(mgl*mt) - (4*lgu*ltu*
     mmgl*mmsb2*mmst1*pow2(invdgb2)*pow2(mmt))/(3.*mgl*mt) - (8*lt1u*ltu*mmgl*
     mmsb2*mmst1*pow2(invdgb2)*pow2(mmt))/(3.*mgl*mt) - (4*mmgl*mmsb22*zt2*pow2
     (invdgb2)*pow2(mmt))/(3.*mgl*mt) - (8*mmgl*mmsb2*mmst1*zt2*pow2(invdgb2)*
     pow2(mmt))/(3.*mgl*mt) - (4*mmgl*pow2(lgu)*pow2(mmt))/(3.*mgl*mt) + (4*
     invdgb1*mmgl*mmsb1*pow2(lgu)*pow2(mmt))/(3.*mgl*mt) + (4*invdgb2*mmgl*
     mmsb2*pow2(lgu)*pow2(mmt))/(3.*mgl*mt) - (2*mmgl*mmsb12*pow2(invdgb1)*pow2
     (lgu)*pow2(mmt))/(3.*mgl*mt) - (2*mmgl*mmsb22*pow2(invdgb2)*pow2(lgu)*pow2
     (mmt))/(3.*mgl*mt) + (2*invdgb1*mmgl*mmst1*pow2(lt1u)*pow2(mmt))/(3.*mgl*
     mt) + (2*invdgb2*mmgl*mmst1*pow2(lt1u)*pow2(mmt))/(3.*mgl*mt) - (2*mmgl*
     mmsb1*mmst1*pow2(invdgb1)*pow2(lt1u)*pow2(mmt))/(3.*mgl*mt) - (2*mmgl*
     mmsb2*mmst1*pow2(invdgb2)*pow2(lt1u)*pow2(mmt))/(3.*mgl*mt) - (4*mmgl*pow2
     (ltu)*pow2(mmt))/(3.*mgl*mt) + (4*invdgb1*mmgl*mmsb1*pow2(ltu)*pow2(mmt))/
     (3.*mgl*mt) + (4*invdgb2*mmgl*mmsb2*pow2(ltu)*pow2(mmt))/(3.*mgl*mt) + (2*
     invdgb1*mmgl*mmst1*pow2(ltu)*pow2(mmt))/(mgl*mt) + (2*invdgb2*mmgl*mmst1*
     pow2(ltu)*pow2(mmt))/(mgl*mt) - (2*mmgl*mmsb12*pow2(invdgb1)*pow2(ltu)*
     pow2(mmt))/(3.*mgl*mt) - (2*mmgl*mmsb1*mmst1*pow2(invdgb1)*pow2(ltu)*pow2(
     mmt))/(mgl*mt) - (2*mmgl*mmsb22*pow2(invdgb2)*pow2(ltu)*pow2(mmt))/(3.*mgl
     *mt) - (2*mmgl*mmsb2*mmst1*pow2(invdgb2)*pow2(ltu)*pow2(mmt))/(mgl*mt) + (
     28*mmgl*mmt*pow2(invdgb1)*pow3(mmsb1))/(3.*mgl*mt) - (4*lgu*mmgl*mmt*pow2(
     invdgb1)*pow3(mmsb1))/(mgl*mt) - (4*ltu*mmgl*mmt*pow2(invdgb1)*pow3(mmsb1)
     )/(mgl*mt) + (4*lgu*ltu*mmgl*mmt*pow2(invdgb1)*pow3(mmsb1))/(3.*mgl*mt) +
     (4*mmgl*mmt*zt2*pow2(invdgb1)*pow3(mmsb1))/(3.*mgl*mt) + (2*mmgl*mmt*pow2(
     invdgb1)*pow2(lgu)*pow3(mmsb1))/(3.*mgl*mt) + (2*mmgl*mmt*pow2(invdgb1)*
     pow2(ltu)*pow3(mmsb1))/(3.*mgl*mt) + (28*mmgl*mmt*pow2(invdgb2)*pow3(mmsb2
     ))/(3.*mgl*mt) - (4*lgu*mmgl*mmt*pow2(invdgb2)*pow3(mmsb2))/(mgl*mt) - (4*
     ltu*mmgl*mmt*pow2(invdgb2)*pow3(mmsb2))/(mgl*mt) + (4*lgu*ltu*mmgl*mmt*
     pow2(invdgb2)*pow3(mmsb2))/(3.*mgl*mt) + (4*mmgl*mmt*zt2*pow2(invdgb2)*
     pow3(mmsb2))/(3.*mgl*mt) + (2*mmgl*mmt*pow2(invdgb2)*pow2(lgu)*pow3(mmsb2)
     )/(3.*mgl*mt) + (2*mmgl*mmt*pow2(invdgb2)*pow2(ltu)*pow3(mmsb2))/(3.*mgl*
     mt) + s2b*((-56*mmgl*mmsb1*mmt)/(3.*mb*mt) + (8*lgu*mmgl*mmsb1*mmt)/(mb*mt
     ) + (8*ltu*mmgl*mmsb1*mmt)/(mb*mt) - (8*lgu*ltu*mmgl*mmsb1*mmt)/(3.*mb*mt)
     - (56*mmsb12*mmt)/(3.*mb*mt) + (8*lgu*mmsb12*mmt)/(mb*mt) + (8*ltu*mmsb12*
     mmt)/(mb*mt) - (8*lgu*ltu*mmsb12*mmt)/(3.*mb*mt) + (56*mmgl*mmsb2*mmt)/(3.
     *mb*mt) - (8*lgu*mmgl*mmsb2*mmt)/(mb*mt) - (8*ltu*mmgl*mmsb2*mmt)/(mb*mt)
     + (8*lgu*ltu*mmgl*mmsb2*mmt)/(3.*mb*mt) + (56*mmsb22*mmt)/(3.*mb*mt) - (8*
     lgu*mmsb22*mmt)/(mb*mt) - (8*ltu*mmsb22*mmt)/(mb*mt) + (8*lgu*ltu*mmsb22*
     mmt)/(3.*mb*mt) + (56*mmsb1*mmst1*mmt)/(3.*mb*mt) - (8*lgu*mmsb1*mmst1*mmt
     )/(mb*mt) + (8*lt1u*mmsb1*mmst1*mmt)/(mb*mt) - (8*lgu*lt1u*mmsb1*mmst1*mmt
     )/(3.*mb*mt) - (16*ltu*mmsb1*mmst1*mmt)/(mb*mt) + (16*lgu*ltu*mmsb1*mmst1*
     mmt)/(3.*mb*mt) - (56*invdgb1*mmsb12*mmst1*mmt)/(3.*mb*mt) + (8*invdgb1*
     lgu*mmsb12*mmst1*mmt)/(mb*mt) - (8*invdgb1*lt1u*mmsb12*mmst1*mmt)/(mb*mt)
     + (8*invdgb1*lgu*lt1u*mmsb12*mmst1*mmt)/(3.*mb*mt) + (16*invdgb1*ltu*
     mmsb12*mmst1*mmt)/(mb*mt) - (16*invdgb1*lgu*ltu*mmsb12*mmst1*mmt)/(3.*mb*
     mt) - (56*mmsb2*mmst1*mmt)/(3.*mb*mt) + (8*lgu*mmsb2*mmst1*mmt)/(mb*mt) -
     (8*lt1u*mmsb2*mmst1*mmt)/(mb*mt) + (8*lgu*lt1u*mmsb2*mmst1*mmt)/(3.*mb*mt)
     + (16*ltu*mmsb2*mmst1*mmt)/(mb*mt) - (16*lgu*ltu*mmsb2*mmst1*mmt)/(3.*mb*
     mt) + (56*invdgb2*mmsb22*mmst1*mmt)/(3.*mb*mt) - (8*invdgb2*lgu*mmsb22*
     mmst1*mmt)/(mb*mt) + (8*invdgb2*lt1u*mmsb22*mmst1*mmt)/(mb*mt) - (8*
     invdgb2*lgu*lt1u*mmsb22*mmst1*mmt)/(3.*mb*mt) - (16*invdgb2*ltu*mmsb22*
     mmst1*mmt)/(mb*mt) + (16*invdgb2*lgu*ltu*mmsb22*mmst1*mmt)/(3.*mb*mt) - (8
     *mmgl*mmsb1*mmt*zt2)/(3.*mb*mt) - (8*mmsb12*mmt*zt2)/(3.*mb*mt) + (8*mmgl*
     mmsb2*mmt*zt2)/(3.*mb*mt) + (8*mmsb22*mmt*zt2)/(3.*mb*mt) + (8*mmsb1*mmst1
     *mmt*zt2)/(3.*mb*mt) - (8*invdgb1*mmsb12*mmst1*mmt*zt2)/(3.*mb*mt) - (8*
     mmsb2*mmst1*mmt*zt2)/(3.*mb*mt) + (8*invdgb2*mmsb22*mmst1*mmt*zt2)/(3.*mb*
     mt) - (4*mmgl*mmsb1*mmt*pow2(lgu))/(3.*mb*mt) - (4*mmsb12*mmt*pow2(lgu))/(
     3.*mb*mt) + (4*mmgl*mmsb2*mmt*pow2(lgu))/(3.*mb*mt) + (4*mmsb22*mmt*pow2(
     lgu))/(3.*mb*mt) + (4*mmsb1*mmst1*mmt*pow2(lgu))/(3.*mb*mt) - (4*invdgb1*
     mmsb12*mmst1*mmt*pow2(lgu))/(3.*mb*mt) - (4*mmsb2*mmst1*mmt*pow2(lgu))/(3.
     *mb*mt) + (4*invdgb2*mmsb22*mmst1*mmt*pow2(lgu))/(3.*mb*mt) - (4*mmsb1*
     mmst1*mmt*pow2(lt1u))/(3.*mb*mt) + (4*invdgb1*mmsb12*mmst1*mmt*pow2(lt1u))
     /(3.*mb*mt) + (4*mmsb2*mmst1*mmt*pow2(lt1u))/(3.*mb*mt) - (4*invdgb2*
     mmsb22*mmst1*mmt*pow2(lt1u))/(3.*mb*mt) - (4*mmgl*mmsb1*mmt*pow2(ltu))/(3.
     *mb*mt) - (4*mmsb12*mmt*pow2(ltu))/(3.*mb*mt) + (4*mmgl*mmsb2*mmt*pow2(ltu
     ))/(3.*mb*mt) + (4*mmsb22*mmt*pow2(ltu))/(3.*mb*mt) + (8*mmsb1*mmst1*mmt*
     pow2(ltu))/(3.*mb*mt) - (8*invdgb1*mmsb12*mmst1*mmt*pow2(ltu))/(3.*mb*mt)
     - (8*mmsb2*mmst1*mmt*pow2(ltu))/(3.*mb*mt) + (8*invdgb2*mmsb22*mmst1*mmt*
     pow2(ltu))/(3.*mb*mt) + (8*invdgb1*lt1u*mmsb1*mmt*pow2(mmst1))/(mb*mt) - (
     8*invdgb1*lgu*lt1u*mmsb1*mmt*pow2(mmst1))/(3.*mb*mt) - (8*invdgb1*ltu*
     mmsb1*mmt*pow2(mmst1))/(mb*mt) + (8*invdgb1*lgu*ltu*mmsb1*mmt*pow2(mmst1))
     /(3.*mb*mt) - (8*invdgb2*lt1u*mmsb2*mmt*pow2(mmst1))/(mb*mt) + (8*invdgb2*
     lgu*lt1u*mmsb2*mmt*pow2(mmst1))/(3.*mb*mt) + (8*invdgb2*ltu*mmsb2*mmt*pow2
     (mmst1))/(mb*mt) - (8*invdgb2*lgu*ltu*mmsb2*mmt*pow2(mmst1))/(3.*mb*mt) -
     (4*invdgb1*mmsb1*mmt*pow2(lt1u)*pow2(mmst1))/(3.*mb*mt) + (4*invdgb2*mmsb2
     *mmt*pow2(lt1u)*pow2(mmst1))/(3.*mb*mt) + (4*invdgb1*mmsb1*mmt*pow2(ltu)*
     pow2(mmst1))/(3.*mb*mt) - (4*invdgb2*mmsb2*mmt*pow2(ltu)*pow2(mmst1))/(3.*
     mb*mt) + (56*mmsb1*pow2(mmt))/(3.*mb*mt) - (8*lgu*mmsb1*pow2(mmt))/(mb*mt)
     - (8*ltu*mmsb1*pow2(mmt))/(mb*mt) + (8*lgu*ltu*mmsb1*pow2(mmt))/(3.*mb*mt)
     - (56*invdgb1*mmsb12*pow2(mmt))/(3.*mb*mt) + (8*invdgb1*lgu*mmsb12*pow2(
     mmt))/(mb*mt) + (8*invdgb1*ltu*mmsb12*pow2(mmt))/(mb*mt) - (8*invdgb1*lgu*
     ltu*mmsb12*pow2(mmt))/(3.*mb*mt) - (56*mmsb2*pow2(mmt))/(3.*mb*mt) + (8*
     lgu*mmsb2*pow2(mmt))/(mb*mt) + (8*ltu*mmsb2*pow2(mmt))/(mb*mt) - (8*lgu*
     ltu*mmsb2*pow2(mmt))/(3.*mb*mt) + (56*invdgb2*mmsb22*pow2(mmt))/(3.*mb*mt)
     - (8*invdgb2*lgu*mmsb22*pow2(mmt))/(mb*mt) - (8*invdgb2*ltu*mmsb22*pow2(
     mmt))/(mb*mt) + (8*invdgb2*lgu*ltu*mmsb22*pow2(mmt))/(3.*mb*mt) - (112*
     invdgb1*mmsb1*mmst1*pow2(mmt))/(3.*mb*mt) + (8*invdgb1*lt1u*mmsb1*mmst1*
     pow2(mmt))/(mb*mt) + (8*invdgb1*lgu*lt1u*mmsb1*mmst1*pow2(mmt))/(3.*mb*mt)
     + (24*invdgb1*ltu*mmsb1*mmst1*pow2(mmt))/(mb*mt) - (8*invdgb1*lgu*ltu*
     mmsb1*mmst1*pow2(mmt))/(3.*mb*mt) - (16*invdgb1*lt1u*ltu*mmsb1*mmst1*pow2(
     mmt))/(3.*mb*mt) + (112*invdgb2*mmsb2*mmst1*pow2(mmt))/(3.*mb*mt) - (8*
     invdgb2*lt1u*mmsb2*mmst1*pow2(mmt))/(mb*mt) - (8*invdgb2*lgu*lt1u*mmsb2*
     mmst1*pow2(mmt))/(3.*mb*mt) - (24*invdgb2*ltu*mmsb2*mmst1*pow2(mmt))/(mb*
     mt) + (8*invdgb2*lgu*ltu*mmsb2*mmst1*pow2(mmt))/(3.*mb*mt) + (16*invdgb2*
     lt1u*ltu*mmsb2*mmst1*pow2(mmt))/(3.*mb*mt) + (8*mmsb1*zt2*pow2(mmt))/(3.*
     mb*mt) - (8*invdgb1*mmsb12*zt2*pow2(mmt))/(3.*mb*mt) - (8*mmsb2*zt2*pow2(
     mmt))/(3.*mb*mt) + (8*invdgb2*mmsb22*zt2*pow2(mmt))/(3.*mb*mt) - (16*
     invdgb1*mmsb1*mmst1*zt2*pow2(mmt))/(3.*mb*mt) + (16*invdgb2*mmsb2*mmst1*
     zt2*pow2(mmt))/(3.*mb*mt) + (4*mmsb1*pow2(lgu)*pow2(mmt))/(3.*mb*mt) - (4*
     invdgb1*mmsb12*pow2(lgu)*pow2(mmt))/(3.*mb*mt) - (4*mmsb2*pow2(lgu)*pow2(
     mmt))/(3.*mb*mt) + (4*invdgb2*mmsb22*pow2(lgu)*pow2(mmt))/(3.*mb*mt) - (4*
     invdgb1*mmsb1*mmst1*pow2(lt1u)*pow2(mmt))/(3.*mb*mt) + (4*invdgb2*mmsb2*
     mmst1*pow2(lt1u)*pow2(mmt))/(3.*mb*mt) + (4*mmsb1*pow2(ltu)*pow2(mmt))/(3.
     *mb*mt) - (4*invdgb1*mmsb12*pow2(ltu)*pow2(mmt))/(3.*mb*mt) - (4*mmsb2*
     pow2(ltu)*pow2(mmt))/(3.*mb*mt) + (4*invdgb2*mmsb22*pow2(ltu)*pow2(mmt))/(
     3.*mb*mt) - (4*invdgb1*mmsb1*mmst1*pow2(ltu)*pow2(mmt))/(mb*mt) + (4*
     invdgb2*mmsb2*mmst1*pow2(ltu)*pow2(mmt))/(mb*mt) + (56*invdgb1*mmt*pow3(
     mmsb1))/(3.*mb*mt) - (8*invdgb1*lgu*mmt*pow3(mmsb1))/(mb*mt) - (8*invdgb1*
     ltu*mmt*pow3(mmsb1))/(mb*mt) + (8*invdgb1*lgu*ltu*mmt*pow3(mmsb1))/(3.*mb*
     mt) + (8*invdgb1*mmt*zt2*pow3(mmsb1))/(3.*mb*mt) + (4*invdgb1*mmt*pow2(lgu
     )*pow3(mmsb1))/(3.*mb*mt) + (4*invdgb1*mmt*pow2(ltu)*pow3(mmsb1))/(3.*mb*
     mt) - (56*invdgb2*mmt*pow3(mmsb2))/(3.*mb*mt) + (8*invdgb2*lgu*mmt*pow3(
     mmsb2))/(mb*mt) + (8*invdgb2*ltu*mmt*pow3(mmsb2))/(mb*mt) - (8*invdgb2*lgu
     *ltu*mmt*pow3(mmsb2))/(3.*mb*mt) - (8*invdgb2*mmt*zt2*pow3(mmsb2))/(3.*mb*
     mt) - (4*invdgb2*mmt*pow2(lgu)*pow3(mmsb2))/(3.*mb*mt) - (4*invdgb2*mmt*
     pow2(ltu)*pow3(mmsb2))/(3.*mb*mt))) + (28*invdgb1*pow3(mmst1))/3. + (28*
     invdgb2*pow3(mmst1))/3. - 4*invdgb1*lt1u*pow3(mmst1) - 4*invdgb2*lt1u*pow3
     (mmst1) - 4*invdgb1*ltu*pow3(mmst1) - 4*invdgb2*ltu*pow3(mmst1) + (4*
     invdgb1*lt1u*ltu*pow3(mmst1))/3. + (4*invdgb2*lt1u*ltu*pow3(mmst1))/3. + (
     4*invdgb1*zt2*pow3(mmst1))/3. + (4*invdgb2*zt2*pow3(mmst1))/3. - (28*mmsb1
     *pow2(invdgb1)*pow3(mmst1))/3. + 4*lt1u*mmsb1*pow2(invdgb1)*pow3(mmst1) +
     4*ltu*mmsb1*pow2(invdgb1)*pow3(mmst1) - (4*lt1u*ltu*mmsb1*pow2(invdgb1)*
     pow3(mmst1))/3. - (4*mmsb1*zt2*pow2(invdgb1)*pow3(mmst1))/3. - (28*mmsb2*
     pow2(invdgb2)*pow3(mmst1))/3. + 4*lt1u*mmsb2*pow2(invdgb2)*pow3(mmst1) + 4
     *ltu*mmsb2*pow2(invdgb2)*pow3(mmst1) - (4*lt1u*ltu*mmsb2*pow2(invdgb2)*
     pow3(mmst1))/3. - (4*mmsb2*zt2*pow2(invdgb2)*pow3(mmst1))/3. + (2*invdgb1*
     pow2(lt1u)*pow3(mmst1))/3. + (2*invdgb2*pow2(lt1u)*pow3(mmst1))/3. - (2*
     mmsb1*pow2(invdgb1)*pow2(lt1u)*pow3(mmst1))/3. - (2*mmsb2*pow2(invdgb2)*
     pow2(lt1u)*pow3(mmst1))/3. + (2*invdgb1*pow2(ltu)*pow3(mmst1))/3. + (2*
     invdgb2*pow2(ltu)*pow3(mmst1))/3. - (2*mmsb1*pow2(invdgb1)*pow2(ltu)*pow3(
     mmst1))/3. - (2*mmsb2*pow2(invdgb2)*pow2(ltu)*pow3(mmst1))/3. + s2b*((56*
     mmgl2*mmsb1)/(3.*mb*mgl) - (8*lgu*mmgl2*mmsb1)/(mb*mgl) - (8*ltu*mmgl2*
     mmsb1)/(mb*mgl) + (8*lgu*ltu*mmgl2*mmsb1)/(3.*mb*mgl) + (56*mmgl*mmsb12)/(
     3.*mb*mgl) - (8*lgu*mmgl*mmsb12)/(mb*mgl) - (8*ltu*mmgl*mmsb12)/(mb*mgl) +
     (8*lgu*ltu*mmgl*mmsb12)/(3.*mb*mgl) - (56*mmgl2*mmsb2)/(3.*mb*mgl) + (8*
     lgu*mmgl2*mmsb2)/(mb*mgl) + (8*ltu*mmgl2*mmsb2)/(mb*mgl) - (8*lgu*ltu*
     mmgl2*mmsb2)/(3.*mb*mgl) - (56*mmgl*mmsb22)/(3.*mb*mgl) + (8*lgu*mmgl*
     mmsb22)/(mb*mgl) + (8*ltu*mmgl*mmsb22)/(mb*mgl) - (8*lgu*ltu*mmgl*mmsb22)/
     (3.*mb*mgl) - (56*mmgl*mmsb1*mmst1)/(3.*mb*mgl) + (16*lgu*mmgl*mmsb1*mmst1
     )/(mb*mgl) - (8*lt1u*mmgl*mmsb1*mmst1)/(mb*mgl) + (8*ltu*mmgl*mmsb1*mmst1)
     /(mb*mgl) - (16*lgu*ltu*mmgl*mmsb1*mmst1)/(3.*mb*mgl) + (8*lt1u*ltu*mmgl*
     mmsb1*mmst1)/(3.*mb*mgl) + (56*invdgb1*mmgl*mmsb12*mmst1)/(3.*mb*mgl) - (
     16*invdgb1*lgu*mmgl*mmsb12*mmst1)/(mb*mgl) + (8*invdgb1*lt1u*mmgl*mmsb12*
     mmst1)/(mb*mgl) - (8*invdgb1*ltu*mmgl*mmsb12*mmst1)/(mb*mgl) + (16*invdgb1
     *lgu*ltu*mmgl*mmsb12*mmst1)/(3.*mb*mgl) - (8*invdgb1*lt1u*ltu*mmgl*mmsb12*
     mmst1)/(3.*mb*mgl) + (56*mmgl*mmsb2*mmst1)/(3.*mb*mgl) - (16*lgu*mmgl*
     mmsb2*mmst1)/(mb*mgl) + (8*lt1u*mmgl*mmsb2*mmst1)/(mb*mgl) - (8*ltu*mmgl*
     mmsb2*mmst1)/(mb*mgl) + (16*lgu*ltu*mmgl*mmsb2*mmst1)/(3.*mb*mgl) - (8*
     lt1u*ltu*mmgl*mmsb2*mmst1)/(3.*mb*mgl) - (56*invdgb2*mmgl*mmsb22*mmst1)/(
     3.*mb*mgl) + (16*invdgb2*lgu*mmgl*mmsb22*mmst1)/(mb*mgl) - (8*invdgb2*lt1u
     *mmgl*mmsb22*mmst1)/(mb*mgl) + (8*invdgb2*ltu*mmgl*mmsb22*mmst1)/(mb*mgl)
     - (16*invdgb2*lgu*ltu*mmgl*mmsb22*mmst1)/(3.*mb*mgl) + (8*invdgb2*lt1u*ltu
     *mmgl*mmsb22*mmst1)/(3.*mb*mgl) - (56*mmgl*mmsb1*mmt)/(3.*mb*mgl) + (8*lgu
     *mmgl*mmsb1*mmt)/(mb*mgl) + (8*ltu*mmgl*mmsb1*mmt)/(mb*mgl) - (8*lgu*ltu*
     mmgl*mmsb1*mmt)/(3.*mb*mgl) + (56*invdgb1*mmgl*mmsb12*mmt)/(3.*mb*mgl) - (
     8*invdgb1*lgu*mmgl*mmsb12*mmt)/(mb*mgl) - (8*invdgb1*ltu*mmgl*mmsb12*mmt)/
     (mb*mgl) + (8*invdgb1*lgu*ltu*mmgl*mmsb12*mmt)/(3.*mb*mgl) + (56*mmgl*
     mmsb2*mmt)/(3.*mb*mgl) - (8*lgu*mmgl*mmsb2*mmt)/(mb*mgl) - (8*ltu*mmgl*
     mmsb2*mmt)/(mb*mgl) + (8*lgu*ltu*mmgl*mmsb2*mmt)/(3.*mb*mgl) - (56*invdgb2
     *mmgl*mmsb22*mmt)/(3.*mb*mgl) + (8*invdgb2*lgu*mmgl*mmsb22*mmt)/(mb*mgl) +
     (8*invdgb2*ltu*mmgl*mmsb22*mmt)/(mb*mgl) - (8*invdgb2*lgu*ltu*mmgl*mmsb22*
     mmt)/(3.*mb*mgl) + (224*invdgb1*mmgl*mmsb1*mmst1*mmt)/(3.*mb*mgl) - (8*
     invdgb1*lgu*mmgl*mmsb1*mmst1*mmt)/(mb*mgl) - (8*invdgb1*lt1u*mmgl*mmsb1*
     mmst1*mmt)/(mb*mgl) - (16*invdgb1*lgu*lt1u*mmgl*mmsb1*mmst1*mmt)/(3.*mb*
     mgl) - (48*invdgb1*ltu*mmgl*mmsb1*mmst1*mmt)/(mb*mgl) + (8*invdgb1*lgu*ltu
     *mmgl*mmsb1*mmst1*mmt)/(mb*mgl) + (8*invdgb1*lt1u*ltu*mmgl*mmsb1*mmst1*mmt
     )/(mb*mgl) - (224*invdgb2*mmgl*mmsb2*mmst1*mmt)/(3.*mb*mgl) + (8*invdgb2*
     lgu*mmgl*mmsb2*mmst1*mmt)/(mb*mgl) + (8*invdgb2*lt1u*mmgl*mmsb2*mmst1*mmt)
     /(mb*mgl) + (16*invdgb2*lgu*lt1u*mmgl*mmsb2*mmst1*mmt)/(3.*mb*mgl) + (48*
     invdgb2*ltu*mmgl*mmsb2*mmst1*mmt)/(mb*mgl) - (8*invdgb2*lgu*ltu*mmgl*mmsb2
     *mmst1*mmt)/(mb*mgl) - (8*invdgb2*lt1u*ltu*mmgl*mmsb2*mmst1*mmt)/(mb*mgl)
     + (8*mmgl2*mmsb1*zt2)/(3.*mb*mgl) + (8*mmgl*mmsb12*zt2)/(3.*mb*mgl) - (8*
     mmgl2*mmsb2*zt2)/(3.*mb*mgl) - (8*mmgl*mmsb22*zt2)/(3.*mb*mgl) - (8*mmgl*
     mmsb1*mmst1*zt2)/(3.*mb*mgl) + (8*invdgb1*mmgl*mmsb12*mmst1*zt2)/(3.*mb*
     mgl) + (8*mmgl*mmsb2*mmst1*zt2)/(3.*mb*mgl) - (8*invdgb2*mmgl*mmsb22*mmst1
     *zt2)/(3.*mb*mgl) - (8*mmgl*mmsb1*mmt*zt2)/(3.*mb*mgl) + (8*invdgb1*mmgl*
     mmsb12*mmt*zt2)/(3.*mb*mgl) + (8*mmgl*mmsb2*mmt*zt2)/(3.*mb*mgl) - (8*
     invdgb2*mmgl*mmsb22*mmt*zt2)/(3.*mb*mgl) + (32*invdgb1*mmgl*mmsb1*mmst1*
     mmt*zt2)/(3.*mb*mgl) - (32*invdgb2*mmgl*mmsb2*mmst1*mmt*zt2)/(3.*mb*mgl) +
     (4*mmgl2*mmsb1*pow2(lgu))/(3.*mb*mgl) + (4*mmgl*mmsb12*pow2(lgu))/(3.*mb*
     mgl) - (4*mmgl2*mmsb2*pow2(lgu))/(3.*mb*mgl) - (4*mmgl*mmsb22*pow2(lgu))/(
     3.*mb*mgl) - (8*mmgl*mmsb1*mmst1*pow2(lgu))/(3.*mb*mgl) + (8*invdgb1*mmgl*
     mmsb12*mmst1*pow2(lgu))/(3.*mb*mgl) + (8*mmgl*mmsb2*mmst1*pow2(lgu))/(3.*
     mb*mgl) - (8*invdgb2*mmgl*mmsb22*mmst1*pow2(lgu))/(3.*mb*mgl) - (4*mmgl*
     mmsb1*mmt*pow2(lgu))/(3.*mb*mgl) + (4*invdgb1*mmgl*mmsb12*mmt*pow2(lgu))/(
     3.*mb*mgl) + (4*mmgl*mmsb2*mmt*pow2(lgu))/(3.*mb*mgl) - (4*invdgb2*mmgl*
     mmsb22*mmt*pow2(lgu))/(3.*mb*mgl) + (4*invdgb1*mmgl*mmsb1*mmst1*mmt*pow2(
     lgu))/(3.*mb*mgl) - (4*invdgb2*mmgl*mmsb2*mmst1*mmt*pow2(lgu))/(3.*mb*mgl)
     + (4*mmgl*mmsb1*mmst1*pow2(lt1u))/(3.*mb*mgl) - (4*invdgb1*mmgl*mmsb12*
     mmst1*pow2(lt1u))/(3.*mb*mgl) - (4*mmgl*mmsb2*mmst1*pow2(lt1u))/(3.*mb*mgl
     ) + (4*invdgb2*mmgl*mmsb22*mmst1*pow2(lt1u))/(3.*mb*mgl) + (4*invdgb1*mmgl
     *mmsb1*mmst1*mmt*pow2(lt1u))/(3.*mb*mgl) - (4*invdgb2*mmgl*mmsb2*mmst1*mmt
     *pow2(lt1u))/(3.*mb*mgl) + (4*mmgl2*mmsb1*pow2(ltu))/(3.*mb*mgl) + (4*mmgl
     *mmsb12*pow2(ltu))/(3.*mb*mgl) - (4*mmgl2*mmsb2*pow2(ltu))/(3.*mb*mgl) - (
     4*mmgl*mmsb22*pow2(ltu))/(3.*mb*mgl) - (4*mmgl*mmsb1*mmst1*pow2(ltu))/(3.*
     mb*mgl) + (4*invdgb1*mmgl*mmsb12*mmst1*pow2(ltu))/(3.*mb*mgl) + (4*mmgl*
     mmsb2*mmst1*pow2(ltu))/(3.*mb*mgl) - (4*invdgb2*mmgl*mmsb22*mmst1*pow2(ltu
     ))/(3.*mb*mgl) - (4*mmgl*mmsb1*mmt*pow2(ltu))/(3.*mb*mgl) + (4*invdgb1*
     mmgl*mmsb12*mmt*pow2(ltu))/(3.*mb*mgl) + (4*mmgl*mmsb2*mmt*pow2(ltu))/(3.*
     mb*mgl) - (4*invdgb2*mmgl*mmsb22*mmt*pow2(ltu))/(3.*mb*mgl) + (8*invdgb1*
     mmgl*mmsb1*mmst1*mmt*pow2(ltu))/(mb*mgl) - (8*invdgb2*mmgl*mmsb2*mmst1*mmt
     *pow2(ltu))/(mb*mgl) + (56*invdgb1*mmgl*mmsb1*pow2(mmst1))/(3.*mb*mgl) + (
     8*invdgb1*lgu*mmgl*mmsb1*pow2(mmst1))/(mb*mgl) - (16*invdgb1*lt1u*mmgl*
     mmsb1*pow2(mmst1))/(mb*mgl) - (8*invdgb1*ltu*mmgl*mmsb1*pow2(mmst1))/(mb*
     mgl) - (8*invdgb1*lgu*ltu*mmgl*mmsb1*pow2(mmst1))/(3.*mb*mgl) + (16*
     invdgb1*lt1u*ltu*mmgl*mmsb1*pow2(mmst1))/(3.*mb*mgl) - (56*invdgb2*mmgl*
     mmsb2*pow2(mmst1))/(3.*mb*mgl) - (8*invdgb2*lgu*mmgl*mmsb2*pow2(mmst1))/(
     mb*mgl) + (16*invdgb2*lt1u*mmgl*mmsb2*pow2(mmst1))/(mb*mgl) + (8*invdgb2*
     ltu*mmgl*mmsb2*pow2(mmst1))/(mb*mgl) + (8*invdgb2*lgu*ltu*mmgl*mmsb2*pow2(
     mmst1))/(3.*mb*mgl) - (16*invdgb2*lt1u*ltu*mmgl*mmsb2*pow2(mmst1))/(3.*mb*
     mgl) + (56*invdgb1*mmgl*mmt*pow2(mmst1))/(3.*mb*mgl) - (56*invdgb2*mmgl*
     mmt*pow2(mmst1))/(3.*mb*mgl) - (8*invdgb1*lt1u*mmgl*mmt*pow2(mmst1))/(mb*
     mgl) + (8*invdgb2*lt1u*mmgl*mmt*pow2(mmst1))/(mb*mgl) - (8*invdgb1*ltu*
     mmgl*mmt*pow2(mmst1))/(mb*mgl) + (8*invdgb2*ltu*mmgl*mmt*pow2(mmst1))/(mb*
     mgl) + (8*invdgb1*lt1u*ltu*mmgl*mmt*pow2(mmst1))/(3.*mb*mgl) - (8*invdgb2*
     lt1u*ltu*mmgl*mmt*pow2(mmst1))/(3.*mb*mgl) + (8*invdgb1*mmgl*mmsb1*zt2*
     pow2(mmst1))/(3.*mb*mgl) - (8*invdgb2*mmgl*mmsb2*zt2*pow2(mmst1))/(3.*mb*
     mgl) + (8*invdgb1*mmgl*mmt*zt2*pow2(mmst1))/(3.*mb*mgl) - (8*invdgb2*mmgl*
     mmt*zt2*pow2(mmst1))/(3.*mb*mgl) - (4*invdgb1*mmgl*mmsb1*pow2(lgu)*pow2(
     mmst1))/(3.*mb*mgl) + (4*invdgb2*mmgl*mmsb2*pow2(lgu)*pow2(mmst1))/(3.*mb*
     mgl) + (8*invdgb1*mmgl*mmsb1*pow2(lt1u)*pow2(mmst1))/(3.*mb*mgl) - (8*
     invdgb2*mmgl*mmsb2*pow2(lt1u)*pow2(mmst1))/(3.*mb*mgl) + (4*invdgb1*mmgl*
     mmt*pow2(lt1u)*pow2(mmst1))/(3.*mb*mgl) - (4*invdgb2*mmgl*mmt*pow2(lt1u)*
     pow2(mmst1))/(3.*mb*mgl) + (4*invdgb1*mmgl*mmsb1*pow2(ltu)*pow2(mmst1))/(
     3.*mb*mgl) - (4*invdgb2*mmgl*mmsb2*pow2(ltu)*pow2(mmst1))/(3.*mb*mgl) + (4
     *invdgb1*mmgl*mmt*pow2(ltu)*pow2(mmst1))/(3.*mb*mgl) - (4*invdgb2*mmgl*mmt
     *pow2(ltu)*pow2(mmst1))/(3.*mb*mgl) - (56*invdgb1*mmgl*pow3(mmsb1))/(3.*mb
     *mgl) + (8*invdgb1*lgu*mmgl*pow3(mmsb1))/(mb*mgl) + (8*invdgb1*ltu*mmgl*
     pow3(mmsb1))/(mb*mgl) - (8*invdgb1*lgu*ltu*mmgl*pow3(mmsb1))/(3.*mb*mgl) -
     (8*invdgb1*mmgl*zt2*pow3(mmsb1))/(3.*mb*mgl) - (4*invdgb1*mmgl*pow2(lgu)*
     pow3(mmsb1))/(3.*mb*mgl) - (4*invdgb1*mmgl*pow2(ltu)*pow3(mmsb1))/(3.*mb*
     mgl) + (56*invdgb2*mmgl*pow3(mmsb2))/(3.*mb*mgl) - (8*invdgb2*lgu*mmgl*
     pow3(mmsb2))/(mb*mgl) - (8*invdgb2*ltu*mmgl*pow3(mmsb2))/(mb*mgl) + (8*
     invdgb2*lgu*ltu*mmgl*pow3(mmsb2))/(3.*mb*mgl) + (8*invdgb2*mmgl*zt2*pow3(
     mmsb2))/(3.*mb*mgl) + (4*invdgb2*mmgl*pow2(lgu)*pow3(mmsb2))/(3.*mb*mgl) +
     (4*invdgb2*mmgl*pow2(ltu)*pow3(mmsb2))/(3.*mb*mgl) - (56*invdgb1*mmgl*pow3
     (mmst1))/(3.*mb*mgl) + (56*invdgb2*mmgl*pow3(mmst1))/(3.*mb*mgl) + (8*
     invdgb1*lt1u*mmgl*pow3(mmst1))/(mb*mgl) - (8*invdgb2*lt1u*mmgl*pow3(mmst1)
     )/(mb*mgl) + (8*invdgb1*ltu*mmgl*pow3(mmst1))/(mb*mgl) - (8*invdgb2*ltu*
     mmgl*pow3(mmst1))/(mb*mgl) - (8*invdgb1*lt1u*ltu*mmgl*pow3(mmst1))/(3.*mb*
     mgl) + (8*invdgb2*lt1u*ltu*mmgl*pow3(mmst1))/(3.*mb*mgl) - (8*invdgb1*mmgl
     *zt2*pow3(mmst1))/(3.*mb*mgl) + (8*invdgb2*mmgl*zt2*pow3(mmst1))/(3.*mb*
     mgl) - (4*invdgb1*mmgl*pow2(lt1u)*pow3(mmst1))/(3.*mb*mgl) + (4*invdgb2*
     mmgl*pow2(lt1u)*pow3(mmst1))/(3.*mb*mgl) - (4*invdgb1*mmgl*pow2(ltu)*pow3(
     mmst1))/(3.*mb*mgl) + (4*invdgb2*mmgl*pow2(ltu)*pow3(mmst1))/(3.*mb*mgl))
     - (28*pow2(invdgb1)*pow4(mmsb1))/3. + 4*lgu*pow2(invdgb1)*pow4(mmsb1) + 4*
     ltu*pow2(invdgb1)*pow4(mmsb1) - (4*lgu*ltu*pow2(invdgb1)*pow4(mmsb1))/3. -
     (4*zt2*pow2(invdgb1)*pow4(mmsb1))/3. - (2*pow2(invdgb1)*pow2(lgu)*pow4(
     mmsb1))/3. - (2*pow2(invdgb1)*pow2(ltu)*pow4(mmsb1))/3. - (28*pow2(invdgb2
     )*pow4(mmsb2))/3. + 4*lgu*pow2(invdgb2)*pow4(mmsb2) + 4*ltu*pow2(invdgb2)*
     pow4(mmsb2) - (4*lgu*ltu*pow2(invdgb2)*pow4(mmsb2))/3. - (4*zt2*pow2(
     invdgb2)*pow4(mmsb2))/3. - (2*pow2(invdgb2)*pow2(lgu)*pow4(mmsb2))/3. - (2
     *pow2(invdgb2)*pow2(ltu)*pow4(mmsb2))/3.) + DeltaInv(mmt,mmst2,mmgl)*((-56
     *mmgl2)/3. + 8*lgu*mmgl2 + 8*ltu*mmgl2 - (8*lgu*ltu*mmgl2)/3. - (56*mmgl*
     mmsb1)/3. + 8*lgu*mmgl*mmsb1 + 8*ltu*mmgl*mmsb1 - (8*lgu*ltu*mmgl*mmsb1)/
     3. - 28*mmsb12 + 12*lgu*mmsb12 + 12*ltu*mmsb12 - 4*lgu*ltu*mmsb12 - (56*
     mmgl*mmsb2)/3. + 8*lgu*mmgl*mmsb2 + 8*ltu*mmgl*mmsb2 - (8*lgu*ltu*mmgl*
     mmsb2)/3. - 28*mmsb22 + 12*lgu*mmsb22 + 12*ltu*mmsb22 - 4*lgu*ltu*mmsb22 +
     (56*mmgl*mmst2)/3. - 16*lgu*mmgl*mmst2 + 8*lt2u*mmgl*mmst2 - 8*ltu*mmgl*
     mmst2 + (16*lgu*ltu*mmgl*mmst2)/3. - (8*lt2u*ltu*mmgl*mmst2)/3. + (56*
     mmsb1*mmst2)/3. - 16*lgu*mmsb1*mmst2 + 8*lt2u*mmsb1*mmst2 - 8*ltu*mmsb1*
     mmst2 + (16*lgu*ltu*mmsb1*mmst2)/3. - (8*lt2u*ltu*mmsb1*mmst2)/3. - 28*
     invdgb1*mmsb12*mmst2 + 24*invdgb1*lgu*mmsb12*mmst2 - 12*invdgb1*lt2u*
     mmsb12*mmst2 + 12*invdgb1*ltu*mmsb12*mmst2 - 8*invdgb1*lgu*ltu*mmsb12*
     mmst2 + 4*invdgb1*lt2u*ltu*mmsb12*mmst2 + (56*mmsb2*mmst2)/3. - 16*lgu*
     mmsb2*mmst2 + 8*lt2u*mmsb2*mmst2 - 8*ltu*mmsb2*mmst2 + (16*lgu*ltu*mmsb2*
     mmst2)/3. - (8*lt2u*ltu*mmsb2*mmst2)/3. - 28*invdgb2*mmsb22*mmst2 + 24*
     invdgb2*lgu*mmsb22*mmst2 - 12*invdgb2*lt2u*mmsb22*mmst2 + 12*invdgb2*ltu*
     mmsb22*mmst2 - 8*invdgb2*lgu*ltu*mmsb22*mmst2 + 4*invdgb2*lt2u*ltu*mmsb22*
     mmst2 + (56*mmgl*mmt)/3. - 8*lgu*mmgl*mmt - 8*ltu*mmgl*mmt + (8*lgu*ltu*
     mmgl*mmt)/3. + (56*mmsb1*mmt)/3. - 8*lgu*mmsb1*mmt - 8*ltu*mmsb1*mmt + (8*
     lgu*ltu*mmsb1*mmt)/3. - 28*invdgb1*mmsb12*mmt + 12*invdgb1*lgu*mmsb12*mmt
     + 12*invdgb1*ltu*mmsb12*mmt - 4*invdgb1*lgu*ltu*mmsb12*mmt + (56*mmsb2*mmt
     )/3. - 8*lgu*mmsb2*mmt - 8*ltu*mmsb2*mmt + (8*lgu*ltu*mmsb2*mmt)/3. - 28*
     invdgb2*mmsb22*mmt + 12*invdgb2*lgu*mmsb22*mmt + 12*invdgb2*ltu*mmsb22*mmt
      - 4*invdgb2*lgu*ltu*mmsb22*mmt + (224*mmst2*mmt)/3. - 8*lgu*mmst2*mmt - 8
     *lt2u*mmst2*mmt - (16*lgu*lt2u*mmst2*mmt)/3. - 48*ltu*mmst2*mmt + 8*lgu*
     ltu*mmst2*mmt + 8*lt2u*ltu*mmst2*mmt - (224*invdgb1*mmsb1*mmst2*mmt)/3. +
     8*invdgb1*lgu*mmsb1*mmst2*mmt + 8*invdgb1*lt2u*mmsb1*mmst2*mmt + (16*
     invdgb1*lgu*lt2u*mmsb1*mmst2*mmt)/3. + 48*invdgb1*ltu*mmsb1*mmst2*mmt - 8*
     invdgb1*lgu*ltu*mmsb1*mmst2*mmt - 8*invdgb1*lt2u*ltu*mmsb1*mmst2*mmt - (
     224*invdgb2*mmsb2*mmst2*mmt)/3. + 8*invdgb2*lgu*mmsb2*mmst2*mmt + 8*
     invdgb2*lt2u*mmsb2*mmst2*mmt + (16*invdgb2*lgu*lt2u*mmsb2*mmst2*mmt)/3. +
     48*invdgb2*ltu*mmsb2*mmst2*mmt - 8*invdgb2*lgu*ltu*mmsb2*mmst2*mmt - 8*
     invdgb2*lt2u*ltu*mmsb2*mmst2*mmt - (8*mmgl2*zt2)/3. - (8*mmgl*mmsb1*zt2)/
     3. - 4*mmsb12*zt2 - (8*mmgl*mmsb2*zt2)/3. - 4*mmsb22*zt2 + (8*mmgl*mmst2*
     zt2)/3. + (8*mmsb1*mmst2*zt2)/3. - 4*invdgb1*mmsb12*mmst2*zt2 + (8*mmsb2*
     mmst2*zt2)/3. - 4*invdgb2*mmsb22*mmst2*zt2 + (8*mmgl*mmt*zt2)/3. + (8*
     mmsb1*mmt*zt2)/3. - 4*invdgb1*mmsb12*mmt*zt2 + (8*mmsb2*mmt*zt2)/3. - 4*
     invdgb2*mmsb22*mmt*zt2 + (32*mmst2*mmt*zt2)/3. - (32*invdgb1*mmsb1*mmst2*
     mmt*zt2)/3. - (32*invdgb2*mmsb2*mmst2*mmt*zt2)/3. + (112*mmsb12*mmst2*mmt*
     pow2(invdgb1))/3. - 4*lgu*mmsb12*mmst2*mmt*pow2(invdgb1) - 4*lt2u*mmsb12*
     mmst2*mmt*pow2(invdgb1) - (8*lgu*lt2u*mmsb12*mmst2*mmt*pow2(invdgb1))/3. -
     24*ltu*mmsb12*mmst2*mmt*pow2(invdgb1) + 4*lgu*ltu*mmsb12*mmst2*mmt*pow2(
     invdgb1) + 4*lt2u*ltu*mmsb12*mmst2*mmt*pow2(invdgb1) + (16*mmsb12*mmst2*
     mmt*zt2*pow2(invdgb1))/3. + (112*mmsb22*mmst2*mmt*pow2(invdgb2))/3. - 4*
     lgu*mmsb22*mmst2*mmt*pow2(invdgb2) - 4*lt2u*mmsb22*mmst2*mmt*pow2(invdgb2)
     - (8*lgu*lt2u*mmsb22*mmst2*mmt*pow2(invdgb2))/3. - 24*ltu*mmsb22*mmst2*mmt
     *pow2(invdgb2) + 4*lgu*ltu*mmsb22*mmst2*mmt*pow2(invdgb2) + 4*lt2u*ltu*
     mmsb22*mmst2*mmt*pow2(invdgb2) + (16*mmsb22*mmst2*mmt*zt2*pow2(invdgb2))/
     3. - (4*mmgl2*pow2(lgu))/3. - (4*mmgl*mmsb1*pow2(lgu))/3. - 2*mmsb12*pow2(
     lgu) - (4*mmgl*mmsb2*pow2(lgu))/3. - 2*mmsb22*pow2(lgu) + (8*mmgl*mmst2*
     pow2(lgu))/3. + (8*mmsb1*mmst2*pow2(lgu))/3. - 4*invdgb1*mmsb12*mmst2*pow2
     (lgu) + (8*mmsb2*mmst2*pow2(lgu))/3. - 4*invdgb2*mmsb22*mmst2*pow2(lgu) +
     (4*mmgl*mmt*pow2(lgu))/3. + (4*mmsb1*mmt*pow2(lgu))/3. - 2*invdgb1*mmsb12*
     mmt*pow2(lgu) + (4*mmsb2*mmt*pow2(lgu))/3. - 2*invdgb2*mmsb22*mmt*pow2(lgu
     ) + (4*mmst2*mmt*pow2(lgu))/3. - (4*invdgb1*mmsb1*mmst2*mmt*pow2(lgu))/3.
     - (4*invdgb2*mmsb2*mmst2*mmt*pow2(lgu))/3. + (2*mmsb12*mmst2*mmt*pow2(
     invdgb1)*pow2(lgu))/3. + (2*mmsb22*mmst2*mmt*pow2(invdgb2)*pow2(lgu))/3. -
     (4*mmgl*mmst2*pow2(lt2u))/3. - (4*mmsb1*mmst2*pow2(lt2u))/3. + 2*invdgb1*
     mmsb12*mmst2*pow2(lt2u) - (4*mmsb2*mmst2*pow2(lt2u))/3. + 2*invdgb2*mmsb22
     *mmst2*pow2(lt2u) + (4*mmst2*mmt*pow2(lt2u))/3. - (4*invdgb1*mmsb1*mmst2*
     mmt*pow2(lt2u))/3. - (4*invdgb2*mmsb2*mmst2*mmt*pow2(lt2u))/3. + (2*mmsb12
     *mmst2*mmt*pow2(invdgb1)*pow2(lt2u))/3. + (2*mmsb22*mmst2*mmt*pow2(invdgb2
     )*pow2(lt2u))/3. - (4*mmgl2*pow2(ltu))/3. - (4*mmgl*mmsb1*pow2(ltu))/3. -
     2*mmsb12*pow2(ltu) - (4*mmgl*mmsb2*pow2(ltu))/3. - 2*mmsb22*pow2(ltu) + (4
     *mmgl*mmst2*pow2(ltu))/3. + (4*mmsb1*mmst2*pow2(ltu))/3. - 2*invdgb1*
     mmsb12*mmst2*pow2(ltu) + (4*mmsb2*mmst2*pow2(ltu))/3. - 2*invdgb2*mmsb22*
     mmst2*pow2(ltu) + (4*mmgl*mmt*pow2(ltu))/3. + (4*mmsb1*mmt*pow2(ltu))/3. -
     2*invdgb1*mmsb12*mmt*pow2(ltu) + (4*mmsb2*mmt*pow2(ltu))/3. - 2*invdgb2*
     mmsb22*mmt*pow2(ltu) + 8*mmst2*mmt*pow2(ltu) - 8*invdgb1*mmsb1*mmst2*mmt*
     pow2(ltu) - 8*invdgb2*mmsb2*mmst2*mmt*pow2(ltu) + 4*mmsb12*mmst2*mmt*pow2(
     invdgb1)*pow2(ltu) + 4*mmsb22*mmst2*mmt*pow2(invdgb2)*pow2(ltu) + (56*pow2
     (mmst2))/3. + 8*lgu*pow2(mmst2) - 16*lt2u*pow2(mmst2) - 8*ltu*pow2(mmst2)
     - (8*lgu*ltu*pow2(mmst2))/3. + (16*lt2u*ltu*pow2(mmst2))/3. - (56*invdgb1*
     mmsb1*pow2(mmst2))/3. - 8*invdgb1*lgu*mmsb1*pow2(mmst2) + 16*invdgb1*lt2u*
     mmsb1*pow2(mmst2) + 8*invdgb1*ltu*mmsb1*pow2(mmst2) + (8*invdgb1*lgu*ltu*
     mmsb1*pow2(mmst2))/3. - (16*invdgb1*lt2u*ltu*mmsb1*pow2(mmst2))/3. - (56*
     invdgb2*mmsb2*pow2(mmst2))/3. - 8*invdgb2*lgu*mmsb2*pow2(mmst2) + 16*
     invdgb2*lt2u*mmsb2*pow2(mmst2) + 8*invdgb2*ltu*mmsb2*pow2(mmst2) + (8*
     invdgb2*lgu*ltu*mmsb2*pow2(mmst2))/3. - (16*invdgb2*lt2u*ltu*mmsb2*pow2(
     mmst2))/3. - (28*invdgb1*mmt*pow2(mmst2))/3. - (28*invdgb2*mmt*pow2(mmst2)
     )/3. + 4*invdgb1*lt2u*mmt*pow2(mmst2) + 4*invdgb2*lt2u*mmt*pow2(mmst2) + 4
     *invdgb1*ltu*mmt*pow2(mmst2) + 4*invdgb2*ltu*mmt*pow2(mmst2) - (4*invdgb1*
     lt2u*ltu*mmt*pow2(mmst2))/3. - (4*invdgb2*lt2u*ltu*mmt*pow2(mmst2))/3. + (
     8*zt2*pow2(mmst2))/3. - (8*invdgb1*mmsb1*zt2*pow2(mmst2))/3. - (8*invdgb2*
     mmsb2*zt2*pow2(mmst2))/3. - (4*invdgb1*mmt*zt2*pow2(mmst2))/3. - (4*
     invdgb2*mmt*zt2*pow2(mmst2))/3. + (28*mmsb12*pow2(invdgb1)*pow2(mmst2))/3.
      + 4*lgu*mmsb12*pow2(invdgb1)*pow2(mmst2) - 8*lt2u*mmsb12*pow2(invdgb1)*
     pow2(mmst2) - 4*ltu*mmsb12*pow2(invdgb1)*pow2(mmst2) - (4*lgu*ltu*mmsb12*
     pow2(invdgb1)*pow2(mmst2))/3. + (8*lt2u*ltu*mmsb12*pow2(invdgb1)*pow2(
     mmst2))/3. + (28*mmsb1*mmt*pow2(invdgb1)*pow2(mmst2))/3. - 4*lt2u*mmsb1*
     mmt*pow2(invdgb1)*pow2(mmst2) - 4*ltu*mmsb1*mmt*pow2(invdgb1)*pow2(mmst2)
     + (4*lt2u*ltu*mmsb1*mmt*pow2(invdgb1)*pow2(mmst2))/3. + (4*mmsb12*zt2*pow2
     (invdgb1)*pow2(mmst2))/3. + (4*mmsb1*mmt*zt2*pow2(invdgb1)*pow2(mmst2))/3.
      + (28*mmsb22*pow2(invdgb2)*pow2(mmst2))/3. + 4*lgu*mmsb22*pow2(invdgb2)*
     pow2(mmst2) - 8*lt2u*mmsb22*pow2(invdgb2)*pow2(mmst2) - 4*ltu*mmsb22*pow2(
     invdgb2)*pow2(mmst2) - (4*lgu*ltu*mmsb22*pow2(invdgb2)*pow2(mmst2))/3. + (
     8*lt2u*ltu*mmsb22*pow2(invdgb2)*pow2(mmst2))/3. + (28*mmsb2*mmt*pow2(
     invdgb2)*pow2(mmst2))/3. - 4*lt2u*mmsb2*mmt*pow2(invdgb2)*pow2(mmst2) - 4*
     ltu*mmsb2*mmt*pow2(invdgb2)*pow2(mmst2) + (4*lt2u*ltu*mmsb2*mmt*pow2(
     invdgb2)*pow2(mmst2))/3. + (4*mmsb22*zt2*pow2(invdgb2)*pow2(mmst2))/3. + (
     4*mmsb2*mmt*zt2*pow2(invdgb2)*pow2(mmst2))/3. - (4*pow2(lgu)*pow2(mmst2))/
     3. + (4*invdgb1*mmsb1*pow2(lgu)*pow2(mmst2))/3. + (4*invdgb2*mmsb2*pow2(
     lgu)*pow2(mmst2))/3. - (2*mmsb12*pow2(invdgb1)*pow2(lgu)*pow2(mmst2))/3. -
     (2*mmsb22*pow2(invdgb2)*pow2(lgu)*pow2(mmst2))/3. + (8*pow2(lt2u)*pow2(
     mmst2))/3. - (8*invdgb1*mmsb1*pow2(lt2u)*pow2(mmst2))/3. - (8*invdgb2*
     mmsb2*pow2(lt2u)*pow2(mmst2))/3. - (2*invdgb1*mmt*pow2(lt2u)*pow2(mmst2))/
     3. - (2*invdgb2*mmt*pow2(lt2u)*pow2(mmst2))/3. + (4*mmsb12*pow2(invdgb1)*
     pow2(lt2u)*pow2(mmst2))/3. + (2*mmsb1*mmt*pow2(invdgb1)*pow2(lt2u)*pow2(
     mmst2))/3. + (4*mmsb22*pow2(invdgb2)*pow2(lt2u)*pow2(mmst2))/3. + (2*mmsb2
     *mmt*pow2(invdgb2)*pow2(lt2u)*pow2(mmst2))/3. + (4*pow2(ltu)*pow2(mmst2))/
     3. - (4*invdgb1*mmsb1*pow2(ltu)*pow2(mmst2))/3. - (4*invdgb2*mmsb2*pow2(
     ltu)*pow2(mmst2))/3. - (2*invdgb1*mmt*pow2(ltu)*pow2(mmst2))/3. - (2*
     invdgb2*mmt*pow2(ltu)*pow2(mmst2))/3. + (2*mmsb12*pow2(invdgb1)*pow2(ltu)*
     pow2(mmst2))/3. + (2*mmsb1*mmt*pow2(invdgb1)*pow2(ltu)*pow2(mmst2))/3. + (
     2*mmsb22*pow2(invdgb2)*pow2(ltu)*pow2(mmst2))/3. + (2*mmsb2*mmt*pow2(
     invdgb2)*pow2(ltu)*pow2(mmst2))/3. + (112*invdgb1*pow3(mmsb1))/3. - 16*
     invdgb1*lgu*pow3(mmsb1) - 16*invdgb1*ltu*pow3(mmsb1) + (16*invdgb1*lgu*ltu
     *pow3(mmsb1))/3. + (16*invdgb1*zt2*pow3(mmsb1))/3. + (28*mmst2*pow2(
     invdgb1)*pow3(mmsb1))/3. - 8*lgu*mmst2*pow2(invdgb1)*pow3(mmsb1) + 4*lt2u*
     mmst2*pow2(invdgb1)*pow3(mmsb1) - 4*ltu*mmst2*pow2(invdgb1)*pow3(mmsb1) +
     (8*lgu*ltu*mmst2*pow2(invdgb1)*pow3(mmsb1))/3. - (4*lt2u*ltu*mmst2*pow2(
     invdgb1)*pow3(mmsb1))/3. + (28*mmt*pow2(invdgb1)*pow3(mmsb1))/3. - 4*lgu*
     mmt*pow2(invdgb1)*pow3(mmsb1) - 4*ltu*mmt*pow2(invdgb1)*pow3(mmsb1) + (4*
     lgu*ltu*mmt*pow2(invdgb1)*pow3(mmsb1))/3. + (4*mmst2*zt2*pow2(invdgb1)*
     pow3(mmsb1))/3. + (4*mmt*zt2*pow2(invdgb1)*pow3(mmsb1))/3. + (8*invdgb1*
     pow2(lgu)*pow3(mmsb1))/3. + (4*mmst2*pow2(invdgb1)*pow2(lgu)*pow3(mmsb1))/
     3. + (2*mmt*pow2(invdgb1)*pow2(lgu)*pow3(mmsb1))/3. - (2*mmst2*pow2(
     invdgb1)*pow2(lt2u)*pow3(mmsb1))/3. + (8*invdgb1*pow2(ltu)*pow3(mmsb1))/3.
      + (2*mmst2*pow2(invdgb1)*pow2(ltu)*pow3(mmsb1))/3. + (2*mmt*pow2(invdgb1)
     *pow2(ltu)*pow3(mmsb1))/3. + (112*invdgb2*pow3(mmsb2))/3. - 16*invdgb2*lgu
     *pow3(mmsb2) - 16*invdgb2*ltu*pow3(mmsb2) + (16*invdgb2*lgu*ltu*pow3(mmsb2
     ))/3. + (16*invdgb2*zt2*pow3(mmsb2))/3. + (28*mmst2*pow2(invdgb2)*pow3(
     mmsb2))/3. - 8*lgu*mmst2*pow2(invdgb2)*pow3(mmsb2) + 4*lt2u*mmst2*pow2(
     invdgb2)*pow3(mmsb2) - 4*ltu*mmst2*pow2(invdgb2)*pow3(mmsb2) + (8*lgu*ltu*
     mmst2*pow2(invdgb2)*pow3(mmsb2))/3. - (4*lt2u*ltu*mmst2*pow2(invdgb2)*pow3
     (mmsb2))/3. + (28*mmt*pow2(invdgb2)*pow3(mmsb2))/3. - 4*lgu*mmt*pow2(
     invdgb2)*pow3(mmsb2) - 4*ltu*mmt*pow2(invdgb2)*pow3(mmsb2) + (4*lgu*ltu*
     mmt*pow2(invdgb2)*pow3(mmsb2))/3. + (4*mmst2*zt2*pow2(invdgb2)*pow3(mmsb2)
     )/3. + (4*mmt*zt2*pow2(invdgb2)*pow3(mmsb2))/3. + (8*invdgb2*pow2(lgu)*
     pow3(mmsb2))/3. + (4*mmst2*pow2(invdgb2)*pow2(lgu)*pow3(mmsb2))/3. + (2*
     mmt*pow2(invdgb2)*pow2(lgu)*pow3(mmsb2))/3. - (2*mmst2*pow2(invdgb2)*pow2(
     lt2u)*pow3(mmsb2))/3. + (8*invdgb2*pow2(ltu)*pow3(mmsb2))/3. + (2*mmst2*
     pow2(invdgb2)*pow2(ltu)*pow3(mmsb2))/3. + (2*mmt*pow2(invdgb2)*pow2(ltu)*
     pow3(mmsb2))/3. + s2t*((-56*mmgl2*mmt)/(3.*mgl*mt) + (8*lgu*mmgl2*mmt)/(
     mgl*mt) + (8*ltu*mmgl2*mmt)/(mgl*mt) - (8*lgu*ltu*mmgl2*mmt)/(3.*mgl*mt) -
     (56*mmgl*mmsb1*mmt)/(3.*mgl*mt) + (8*lgu*mmgl*mmsb1*mmt)/(mgl*mt) + (8*ltu
     *mmgl*mmsb1*mmt)/(mgl*mt) - (8*lgu*ltu*mmgl*mmsb1*mmt)/(3.*mgl*mt) + (28*
     invdgb1*mmgl*mmsb12*mmt)/(mgl*mt) - (12*invdgb1*lgu*mmgl*mmsb12*mmt)/(mgl*
     mt) - (12*invdgb1*ltu*mmgl*mmsb12*mmt)/(mgl*mt) + (4*invdgb1*lgu*ltu*mmgl*
     mmsb12*mmt)/(mgl*mt) - (56*mmgl*mmsb2*mmt)/(3.*mgl*mt) + (8*lgu*mmgl*mmsb2
     *mmt)/(mgl*mt) + (8*ltu*mmgl*mmsb2*mmt)/(mgl*mt) - (8*lgu*ltu*mmgl*mmsb2*
     mmt)/(3.*mgl*mt) + (28*invdgb2*mmgl*mmsb22*mmt)/(mgl*mt) - (12*invdgb2*lgu
     *mmgl*mmsb22*mmt)/(mgl*mt) - (12*invdgb2*ltu*mmgl*mmsb22*mmt)/(mgl*mt) + (
     4*invdgb2*lgu*ltu*mmgl*mmsb22*mmt)/(mgl*mt) + (56*mmgl*mmst2*mmt)/(3.*mgl*
     mt) - (8*lgu*mmgl*mmst2*mmt)/(mgl*mt) + (8*lt2u*mmgl*mmst2*mmt)/(mgl*mt) -
     (8*lgu*lt2u*mmgl*mmst2*mmt)/(3.*mgl*mt) - (16*ltu*mmgl*mmst2*mmt)/(mgl*mt)
     + (16*lgu*ltu*mmgl*mmst2*mmt)/(3.*mgl*mt) - (56*invdgb1*mmgl*mmsb1*mmst2*
     mmt)/(3.*mgl*mt) + (8*invdgb1*lgu*mmgl*mmsb1*mmst2*mmt)/(mgl*mt) - (8*
     invdgb1*lt2u*mmgl*mmsb1*mmst2*mmt)/(mgl*mt) + (8*invdgb1*lgu*lt2u*mmgl*
     mmsb1*mmst2*mmt)/(3.*mgl*mt) + (16*invdgb1*ltu*mmgl*mmsb1*mmst2*mmt)/(mgl*
     mt) - (16*invdgb1*lgu*ltu*mmgl*mmsb1*mmst2*mmt)/(3.*mgl*mt) - (56*invdgb2*
     mmgl*mmsb2*mmst2*mmt)/(3.*mgl*mt) + (8*invdgb2*lgu*mmgl*mmsb2*mmst2*mmt)/(
     mgl*mt) - (8*invdgb2*lt2u*mmgl*mmsb2*mmst2*mmt)/(mgl*mt) + (8*invdgb2*lgu*
     lt2u*mmgl*mmsb2*mmst2*mmt)/(3.*mgl*mt) + (16*invdgb2*ltu*mmgl*mmsb2*mmst2*
     mmt)/(mgl*mt) - (16*invdgb2*lgu*ltu*mmgl*mmsb2*mmst2*mmt)/(3.*mgl*mt) - (8
     *mmgl2*mmt*zt2)/(3.*mgl*mt) - (8*mmgl*mmsb1*mmt*zt2)/(3.*mgl*mt) + (4*
     invdgb1*mmgl*mmsb12*mmt*zt2)/(mgl*mt) - (8*mmgl*mmsb2*mmt*zt2)/(3.*mgl*mt)
     + (4*invdgb2*mmgl*mmsb22*mmt*zt2)/(mgl*mt) + (8*mmgl*mmst2*mmt*zt2)/(3.*
     mgl*mt) - (8*invdgb1*mmgl*mmsb1*mmst2*mmt*zt2)/(3.*mgl*mt) - (8*invdgb2*
     mmgl*mmsb2*mmst2*mmt*zt2)/(3.*mgl*mt) + (28*mmgl*mmsb12*mmst2*mmt*pow2(
     invdgb1))/(3.*mgl*mt) - (4*lgu*mmgl*mmsb12*mmst2*mmt*pow2(invdgb1))/(mgl*
     mt) + (4*lt2u*mmgl*mmsb12*mmst2*mmt*pow2(invdgb1))/(mgl*mt) - (4*lgu*lt2u*
     mmgl*mmsb12*mmst2*mmt*pow2(invdgb1))/(3.*mgl*mt) - (8*ltu*mmgl*mmsb12*
     mmst2*mmt*pow2(invdgb1))/(mgl*mt) + (8*lgu*ltu*mmgl*mmsb12*mmst2*mmt*pow2(
     invdgb1))/(3.*mgl*mt) + (4*mmgl*mmsb12*mmst2*mmt*zt2*pow2(invdgb1))/(3.*
     mgl*mt) + (28*mmgl*mmsb22*mmst2*mmt*pow2(invdgb2))/(3.*mgl*mt) - (4*lgu*
     mmgl*mmsb22*mmst2*mmt*pow2(invdgb2))/(mgl*mt) + (4*lt2u*mmgl*mmsb22*mmst2*
     mmt*pow2(invdgb2))/(mgl*mt) - (4*lgu*lt2u*mmgl*mmsb22*mmst2*mmt*pow2(
     invdgb2))/(3.*mgl*mt) - (8*ltu*mmgl*mmsb22*mmst2*mmt*pow2(invdgb2))/(mgl*
     mt) + (8*lgu*ltu*mmgl*mmsb22*mmst2*mmt*pow2(invdgb2))/(3.*mgl*mt) + (4*
     mmgl*mmsb22*mmst2*mmt*zt2*pow2(invdgb2))/(3.*mgl*mt) - (4*mmgl2*mmt*pow2(
     lgu))/(3.*mgl*mt) - (4*mmgl*mmsb1*mmt*pow2(lgu))/(3.*mgl*mt) + (2*invdgb1*
     mmgl*mmsb12*mmt*pow2(lgu))/(mgl*mt) - (4*mmgl*mmsb2*mmt*pow2(lgu))/(3.*mgl
     *mt) + (2*invdgb2*mmgl*mmsb22*mmt*pow2(lgu))/(mgl*mt) + (4*mmgl*mmst2*mmt*
     pow2(lgu))/(3.*mgl*mt) - (4*invdgb1*mmgl*mmsb1*mmst2*mmt*pow2(lgu))/(3.*
     mgl*mt) - (4*invdgb2*mmgl*mmsb2*mmst2*mmt*pow2(lgu))/(3.*mgl*mt) + (2*mmgl
     *mmsb12*mmst2*mmt*pow2(invdgb1)*pow2(lgu))/(3.*mgl*mt) + (2*mmgl*mmsb22*
     mmst2*mmt*pow2(invdgb2)*pow2(lgu))/(3.*mgl*mt) - (4*mmgl*mmst2*mmt*pow2(
     lt2u))/(3.*mgl*mt) + (4*invdgb1*mmgl*mmsb1*mmst2*mmt*pow2(lt2u))/(3.*mgl*
     mt) + (4*invdgb2*mmgl*mmsb2*mmst2*mmt*pow2(lt2u))/(3.*mgl*mt) - (2*mmgl*
     mmsb12*mmst2*mmt*pow2(invdgb1)*pow2(lt2u))/(3.*mgl*mt) - (2*mmgl*mmsb22*
     mmst2*mmt*pow2(invdgb2)*pow2(lt2u))/(3.*mgl*mt) - (4*mmgl2*mmt*pow2(ltu))/
     (3.*mgl*mt) - (4*mmgl*mmsb1*mmt*pow2(ltu))/(3.*mgl*mt) + (2*invdgb1*mmgl*
     mmsb12*mmt*pow2(ltu))/(mgl*mt) - (4*mmgl*mmsb2*mmt*pow2(ltu))/(3.*mgl*mt)
     + (2*invdgb2*mmgl*mmsb22*mmt*pow2(ltu))/(mgl*mt) + (8*mmgl*mmst2*mmt*pow2(
     ltu))/(3.*mgl*mt) - (8*invdgb1*mmgl*mmsb1*mmst2*mmt*pow2(ltu))/(3.*mgl*mt)
     - (8*invdgb2*mmgl*mmsb2*mmst2*mmt*pow2(ltu))/(3.*mgl*mt) + (4*mmgl*mmsb12*
     mmst2*mmt*pow2(invdgb1)*pow2(ltu))/(3.*mgl*mt) + (4*mmgl*mmsb22*mmst2*mmt*
     pow2(invdgb2)*pow2(ltu))/(3.*mgl*mt) + (4*invdgb1*lt2u*mmgl*mmt*pow2(mmst2
     ))/(mgl*mt) + (4*invdgb2*lt2u*mmgl*mmt*pow2(mmst2))/(mgl*mt) - (4*invdgb1*
     lgu*lt2u*mmgl*mmt*pow2(mmst2))/(3.*mgl*mt) - (4*invdgb2*lgu*lt2u*mmgl*mmt*
     pow2(mmst2))/(3.*mgl*mt) - (4*invdgb1*ltu*mmgl*mmt*pow2(mmst2))/(mgl*mt) -
     (4*invdgb2*ltu*mmgl*mmt*pow2(mmst2))/(mgl*mt) + (4*invdgb1*lgu*ltu*mmgl*
     mmt*pow2(mmst2))/(3.*mgl*mt) + (4*invdgb2*lgu*ltu*mmgl*mmt*pow2(mmst2))/(
     3.*mgl*mt) - (4*lt2u*mmgl*mmsb1*mmt*pow2(invdgb1)*pow2(mmst2))/(mgl*mt) +
     (4*lgu*lt2u*mmgl*mmsb1*mmt*pow2(invdgb1)*pow2(mmst2))/(3.*mgl*mt) + (4*ltu
     *mmgl*mmsb1*mmt*pow2(invdgb1)*pow2(mmst2))/(mgl*mt) - (4*lgu*ltu*mmgl*
     mmsb1*mmt*pow2(invdgb1)*pow2(mmst2))/(3.*mgl*mt) - (4*lt2u*mmgl*mmsb2*mmt*
     pow2(invdgb2)*pow2(mmst2))/(mgl*mt) + (4*lgu*lt2u*mmgl*mmsb2*mmt*pow2(
     invdgb2)*pow2(mmst2))/(3.*mgl*mt) + (4*ltu*mmgl*mmsb2*mmt*pow2(invdgb2)*
     pow2(mmst2))/(mgl*mt) - (4*lgu*ltu*mmgl*mmsb2*mmt*pow2(invdgb2)*pow2(mmst2
     ))/(3.*mgl*mt) - (2*invdgb1*mmgl*mmt*pow2(lt2u)*pow2(mmst2))/(3.*mgl*mt) -
     (2*invdgb2*mmgl*mmt*pow2(lt2u)*pow2(mmst2))/(3.*mgl*mt) + (2*mmgl*mmsb1*
     mmt*pow2(invdgb1)*pow2(lt2u)*pow2(mmst2))/(3.*mgl*mt) + (2*mmgl*mmsb2*mmt*
     pow2(invdgb2)*pow2(lt2u)*pow2(mmst2))/(3.*mgl*mt) + (2*invdgb1*mmgl*mmt*
     pow2(ltu)*pow2(mmst2))/(3.*mgl*mt) + (2*invdgb2*mmgl*mmt*pow2(ltu)*pow2(
     mmst2))/(3.*mgl*mt) - (2*mmgl*mmsb1*mmt*pow2(invdgb1)*pow2(ltu)*pow2(mmst2
     ))/(3.*mgl*mt) - (2*mmgl*mmsb2*mmt*pow2(invdgb2)*pow2(ltu)*pow2(mmst2))/(
     3.*mgl*mt) + (56*mmgl*pow2(mmt))/(3.*mgl*mt) - (8*lgu*mmgl*pow2(mmt))/(mgl
     *mt) - (8*ltu*mmgl*pow2(mmt))/(mgl*mt) + (8*lgu*ltu*mmgl*pow2(mmt))/(3.*
     mgl*mt) - (56*invdgb1*mmgl*mmsb1*pow2(mmt))/(3.*mgl*mt) + (8*invdgb1*lgu*
     mmgl*mmsb1*pow2(mmt))/(mgl*mt) + (8*invdgb1*ltu*mmgl*mmsb1*pow2(mmt))/(mgl
     *mt) - (8*invdgb1*lgu*ltu*mmgl*mmsb1*pow2(mmt))/(3.*mgl*mt) - (56*invdgb2*
     mmgl*mmsb2*pow2(mmt))/(3.*mgl*mt) + (8*invdgb2*lgu*mmgl*mmsb2*pow2(mmt))/(
     mgl*mt) + (8*invdgb2*ltu*mmgl*mmsb2*pow2(mmt))/(mgl*mt) - (8*invdgb2*lgu*
     ltu*mmgl*mmsb2*pow2(mmt))/(3.*mgl*mt) - (56*invdgb1*mmgl*mmst2*pow2(mmt))/
     (3.*mgl*mt) - (56*invdgb2*mmgl*mmst2*pow2(mmt))/(3.*mgl*mt) + (4*invdgb1*
     lt2u*mmgl*mmst2*pow2(mmt))/(mgl*mt) + (4*invdgb2*lt2u*mmgl*mmst2*pow2(mmt)
     )/(mgl*mt) + (4*invdgb1*lgu*lt2u*mmgl*mmst2*pow2(mmt))/(3.*mgl*mt) + (4*
     invdgb2*lgu*lt2u*mmgl*mmst2*pow2(mmt))/(3.*mgl*mt) + (12*invdgb1*ltu*mmgl*
     mmst2*pow2(mmt))/(mgl*mt) + (12*invdgb2*ltu*mmgl*mmst2*pow2(mmt))/(mgl*mt)
     - (4*invdgb1*lgu*ltu*mmgl*mmst2*pow2(mmt))/(3.*mgl*mt) - (4*invdgb2*lgu*
     ltu*mmgl*mmst2*pow2(mmt))/(3.*mgl*mt) - (8*invdgb1*lt2u*ltu*mmgl*mmst2*
     pow2(mmt))/(3.*mgl*mt) - (8*invdgb2*lt2u*ltu*mmgl*mmst2*pow2(mmt))/(3.*mgl
     *mt) + (8*mmgl*zt2*pow2(mmt))/(3.*mgl*mt) - (8*invdgb1*mmgl*mmsb1*zt2*pow2
     (mmt))/(3.*mgl*mt) - (8*invdgb2*mmgl*mmsb2*zt2*pow2(mmt))/(3.*mgl*mt) - (8
     *invdgb1*mmgl*mmst2*zt2*pow2(mmt))/(3.*mgl*mt) - (8*invdgb2*mmgl*mmst2*zt2
     *pow2(mmt))/(3.*mgl*mt) + (28*mmgl*mmsb12*pow2(invdgb1)*pow2(mmt))/(3.*mgl
     *mt) - (4*lgu*mmgl*mmsb12*pow2(invdgb1)*pow2(mmt))/(mgl*mt) - (4*ltu*mmgl*
     mmsb12*pow2(invdgb1)*pow2(mmt))/(mgl*mt) + (4*lgu*ltu*mmgl*mmsb12*pow2(
     invdgb1)*pow2(mmt))/(3.*mgl*mt) + (56*mmgl*mmsb1*mmst2*pow2(invdgb1)*pow2(
     mmt))/(3.*mgl*mt) - (4*lt2u*mmgl*mmsb1*mmst2*pow2(invdgb1)*pow2(mmt))/(mgl
     *mt) - (4*lgu*lt2u*mmgl*mmsb1*mmst2*pow2(invdgb1)*pow2(mmt))/(3.*mgl*mt) -
     (12*ltu*mmgl*mmsb1*mmst2*pow2(invdgb1)*pow2(mmt))/(mgl*mt) + (4*lgu*ltu*
     mmgl*mmsb1*mmst2*pow2(invdgb1)*pow2(mmt))/(3.*mgl*mt) + (8*lt2u*ltu*mmgl*
     mmsb1*mmst2*pow2(invdgb1)*pow2(mmt))/(3.*mgl*mt) + (4*mmgl*mmsb12*zt2*pow2
     (invdgb1)*pow2(mmt))/(3.*mgl*mt) + (8*mmgl*mmsb1*mmst2*zt2*pow2(invdgb1)*
     pow2(mmt))/(3.*mgl*mt) + (28*mmgl*mmsb22*pow2(invdgb2)*pow2(mmt))/(3.*mgl*
     mt) - (4*lgu*mmgl*mmsb22*pow2(invdgb2)*pow2(mmt))/(mgl*mt) - (4*ltu*mmgl*
     mmsb22*pow2(invdgb2)*pow2(mmt))/(mgl*mt) + (4*lgu*ltu*mmgl*mmsb22*pow2(
     invdgb2)*pow2(mmt))/(3.*mgl*mt) + (56*mmgl*mmsb2*mmst2*pow2(invdgb2)*pow2(
     mmt))/(3.*mgl*mt) - (4*lt2u*mmgl*mmsb2*mmst2*pow2(invdgb2)*pow2(mmt))/(mgl
     *mt) - (4*lgu*lt2u*mmgl*mmsb2*mmst2*pow2(invdgb2)*pow2(mmt))/(3.*mgl*mt) -
     (12*ltu*mmgl*mmsb2*mmst2*pow2(invdgb2)*pow2(mmt))/(mgl*mt) + (4*lgu*ltu*
     mmgl*mmsb2*mmst2*pow2(invdgb2)*pow2(mmt))/(3.*mgl*mt) + (8*lt2u*ltu*mmgl*
     mmsb2*mmst2*pow2(invdgb2)*pow2(mmt))/(3.*mgl*mt) + (4*mmgl*mmsb22*zt2*pow2
     (invdgb2)*pow2(mmt))/(3.*mgl*mt) + (8*mmgl*mmsb2*mmst2*zt2*pow2(invdgb2)*
     pow2(mmt))/(3.*mgl*mt) + (4*mmgl*pow2(lgu)*pow2(mmt))/(3.*mgl*mt) - (4*
     invdgb1*mmgl*mmsb1*pow2(lgu)*pow2(mmt))/(3.*mgl*mt) - (4*invdgb2*mmgl*
     mmsb2*pow2(lgu)*pow2(mmt))/(3.*mgl*mt) + (2*mmgl*mmsb12*pow2(invdgb1)*pow2
     (lgu)*pow2(mmt))/(3.*mgl*mt) + (2*mmgl*mmsb22*pow2(invdgb2)*pow2(lgu)*pow2
     (mmt))/(3.*mgl*mt) - (2*invdgb1*mmgl*mmst2*pow2(lt2u)*pow2(mmt))/(3.*mgl*
     mt) - (2*invdgb2*mmgl*mmst2*pow2(lt2u)*pow2(mmt))/(3.*mgl*mt) + (2*mmgl*
     mmsb1*mmst2*pow2(invdgb1)*pow2(lt2u)*pow2(mmt))/(3.*mgl*mt) + (2*mmgl*
     mmsb2*mmst2*pow2(invdgb2)*pow2(lt2u)*pow2(mmt))/(3.*mgl*mt) + (4*mmgl*pow2
     (ltu)*pow2(mmt))/(3.*mgl*mt) - (4*invdgb1*mmgl*mmsb1*pow2(ltu)*pow2(mmt))/
     (3.*mgl*mt) - (4*invdgb2*mmgl*mmsb2*pow2(ltu)*pow2(mmt))/(3.*mgl*mt) - (2*
     invdgb1*mmgl*mmst2*pow2(ltu)*pow2(mmt))/(mgl*mt) - (2*invdgb2*mmgl*mmst2*
     pow2(ltu)*pow2(mmt))/(mgl*mt) + (2*mmgl*mmsb12*pow2(invdgb1)*pow2(ltu)*
     pow2(mmt))/(3.*mgl*mt) + (2*mmgl*mmsb1*mmst2*pow2(invdgb1)*pow2(ltu)*pow2(
     mmt))/(mgl*mt) + (2*mmgl*mmsb22*pow2(invdgb2)*pow2(ltu)*pow2(mmt))/(3.*mgl
     *mt) + (2*mmgl*mmsb2*mmst2*pow2(invdgb2)*pow2(ltu)*pow2(mmt))/(mgl*mt) - (
     28*mmgl*mmt*pow2(invdgb1)*pow3(mmsb1))/(3.*mgl*mt) + (4*lgu*mmgl*mmt*pow2(
     invdgb1)*pow3(mmsb1))/(mgl*mt) + (4*ltu*mmgl*mmt*pow2(invdgb1)*pow3(mmsb1)
     )/(mgl*mt) - (4*lgu*ltu*mmgl*mmt*pow2(invdgb1)*pow3(mmsb1))/(3.*mgl*mt) -
     (4*mmgl*mmt*zt2*pow2(invdgb1)*pow3(mmsb1))/(3.*mgl*mt) - (2*mmgl*mmt*pow2(
     invdgb1)*pow2(lgu)*pow3(mmsb1))/(3.*mgl*mt) - (2*mmgl*mmt*pow2(invdgb1)*
     pow2(ltu)*pow3(mmsb1))/(3.*mgl*mt) - (28*mmgl*mmt*pow2(invdgb2)*pow3(mmsb2
     ))/(3.*mgl*mt) + (4*lgu*mmgl*mmt*pow2(invdgb2)*pow3(mmsb2))/(mgl*mt) + (4*
     ltu*mmgl*mmt*pow2(invdgb2)*pow3(mmsb2))/(mgl*mt) - (4*lgu*ltu*mmgl*mmt*
     pow2(invdgb2)*pow3(mmsb2))/(3.*mgl*mt) - (4*mmgl*mmt*zt2*pow2(invdgb2)*
     pow3(mmsb2))/(3.*mgl*mt) - (2*mmgl*mmt*pow2(invdgb2)*pow2(lgu)*pow3(mmsb2)
     )/(3.*mgl*mt) - (2*mmgl*mmt*pow2(invdgb2)*pow2(ltu)*pow3(mmsb2))/(3.*mgl*
     mt) + s2b*((56*mmgl*mmsb1*mmt)/(3.*mb*mt) - (8*lgu*mmgl*mmsb1*mmt)/(mb*mt)
     - (8*ltu*mmgl*mmsb1*mmt)/(mb*mt) + (8*lgu*ltu*mmgl*mmsb1*mmt)/(3.*mb*mt) +
     (56*mmsb12*mmt)/(3.*mb*mt) - (8*lgu*mmsb12*mmt)/(mb*mt) - (8*ltu*mmsb12*
     mmt)/(mb*mt) + (8*lgu*ltu*mmsb12*mmt)/(3.*mb*mt) - (56*mmgl*mmsb2*mmt)/(3.
     *mb*mt) + (8*lgu*mmgl*mmsb2*mmt)/(mb*mt) + (8*ltu*mmgl*mmsb2*mmt)/(mb*mt)
     - (8*lgu*ltu*mmgl*mmsb2*mmt)/(3.*mb*mt) - (56*mmsb22*mmt)/(3.*mb*mt) + (8*
     lgu*mmsb22*mmt)/(mb*mt) + (8*ltu*mmsb22*mmt)/(mb*mt) - (8*lgu*ltu*mmsb22*
     mmt)/(3.*mb*mt) - (56*mmsb1*mmst2*mmt)/(3.*mb*mt) + (8*lgu*mmsb1*mmst2*mmt
     )/(mb*mt) - (8*lt2u*mmsb1*mmst2*mmt)/(mb*mt) + (8*lgu*lt2u*mmsb1*mmst2*mmt
     )/(3.*mb*mt) + (16*ltu*mmsb1*mmst2*mmt)/(mb*mt) - (16*lgu*ltu*mmsb1*mmst2*
     mmt)/(3.*mb*mt) + (56*invdgb1*mmsb12*mmst2*mmt)/(3.*mb*mt) - (8*invdgb1*
     lgu*mmsb12*mmst2*mmt)/(mb*mt) + (8*invdgb1*lt2u*mmsb12*mmst2*mmt)/(mb*mt)
     - (8*invdgb1*lgu*lt2u*mmsb12*mmst2*mmt)/(3.*mb*mt) - (16*invdgb1*ltu*
     mmsb12*mmst2*mmt)/(mb*mt) + (16*invdgb1*lgu*ltu*mmsb12*mmst2*mmt)/(3.*mb*
     mt) + (56*mmsb2*mmst2*mmt)/(3.*mb*mt) - (8*lgu*mmsb2*mmst2*mmt)/(mb*mt) +
     (8*lt2u*mmsb2*mmst2*mmt)/(mb*mt) - (8*lgu*lt2u*mmsb2*mmst2*mmt)/(3.*mb*mt)
     - (16*ltu*mmsb2*mmst2*mmt)/(mb*mt) + (16*lgu*ltu*mmsb2*mmst2*mmt)/(3.*mb*
     mt) - (56*invdgb2*mmsb22*mmst2*mmt)/(3.*mb*mt) + (8*invdgb2*lgu*mmsb22*
     mmst2*mmt)/(mb*mt) - (8*invdgb2*lt2u*mmsb22*mmst2*mmt)/(mb*mt) + (8*
     invdgb2*lgu*lt2u*mmsb22*mmst2*mmt)/(3.*mb*mt) + (16*invdgb2*ltu*mmsb22*
     mmst2*mmt)/(mb*mt) - (16*invdgb2*lgu*ltu*mmsb22*mmst2*mmt)/(3.*mb*mt) + (8
     *mmgl*mmsb1*mmt*zt2)/(3.*mb*mt) + (8*mmsb12*mmt*zt2)/(3.*mb*mt) - (8*mmgl*
     mmsb2*mmt*zt2)/(3.*mb*mt) - (8*mmsb22*mmt*zt2)/(3.*mb*mt) - (8*mmsb1*mmst2
     *mmt*zt2)/(3.*mb*mt) + (8*invdgb1*mmsb12*mmst2*mmt*zt2)/(3.*mb*mt) + (8*
     mmsb2*mmst2*mmt*zt2)/(3.*mb*mt) - (8*invdgb2*mmsb22*mmst2*mmt*zt2)/(3.*mb*
     mt) + (4*mmgl*mmsb1*mmt*pow2(lgu))/(3.*mb*mt) + (4*mmsb12*mmt*pow2(lgu))/(
     3.*mb*mt) - (4*mmgl*mmsb2*mmt*pow2(lgu))/(3.*mb*mt) - (4*mmsb22*mmt*pow2(
     lgu))/(3.*mb*mt) - (4*mmsb1*mmst2*mmt*pow2(lgu))/(3.*mb*mt) + (4*invdgb1*
     mmsb12*mmst2*mmt*pow2(lgu))/(3.*mb*mt) + (4*mmsb2*mmst2*mmt*pow2(lgu))/(3.
     *mb*mt) - (4*invdgb2*mmsb22*mmst2*mmt*pow2(lgu))/(3.*mb*mt) + (4*mmsb1*
     mmst2*mmt*pow2(lt2u))/(3.*mb*mt) - (4*invdgb1*mmsb12*mmst2*mmt*pow2(lt2u))
     /(3.*mb*mt) - (4*mmsb2*mmst2*mmt*pow2(lt2u))/(3.*mb*mt) + (4*invdgb2*
     mmsb22*mmst2*mmt*pow2(lt2u))/(3.*mb*mt) + (4*mmgl*mmsb1*mmt*pow2(ltu))/(3.
     *mb*mt) + (4*mmsb12*mmt*pow2(ltu))/(3.*mb*mt) - (4*mmgl*mmsb2*mmt*pow2(ltu
     ))/(3.*mb*mt) - (4*mmsb22*mmt*pow2(ltu))/(3.*mb*mt) - (8*mmsb1*mmst2*mmt*
     pow2(ltu))/(3.*mb*mt) + (8*invdgb1*mmsb12*mmst2*mmt*pow2(ltu))/(3.*mb*mt)
     + (8*mmsb2*mmst2*mmt*pow2(ltu))/(3.*mb*mt) - (8*invdgb2*mmsb22*mmst2*mmt*
     pow2(ltu))/(3.*mb*mt) - (8*invdgb1*lt2u*mmsb1*mmt*pow2(mmst2))/(mb*mt) + (
     8*invdgb1*lgu*lt2u*mmsb1*mmt*pow2(mmst2))/(3.*mb*mt) + (8*invdgb1*ltu*
     mmsb1*mmt*pow2(mmst2))/(mb*mt) - (8*invdgb1*lgu*ltu*mmsb1*mmt*pow2(mmst2))
     /(3.*mb*mt) + (8*invdgb2*lt2u*mmsb2*mmt*pow2(mmst2))/(mb*mt) - (8*invdgb2*
     lgu*lt2u*mmsb2*mmt*pow2(mmst2))/(3.*mb*mt) - (8*invdgb2*ltu*mmsb2*mmt*pow2
     (mmst2))/(mb*mt) + (8*invdgb2*lgu*ltu*mmsb2*mmt*pow2(mmst2))/(3.*mb*mt) +
     (4*invdgb1*mmsb1*mmt*pow2(lt2u)*pow2(mmst2))/(3.*mb*mt) - (4*invdgb2*mmsb2
     *mmt*pow2(lt2u)*pow2(mmst2))/(3.*mb*mt) - (4*invdgb1*mmsb1*mmt*pow2(ltu)*
     pow2(mmst2))/(3.*mb*mt) + (4*invdgb2*mmsb2*mmt*pow2(ltu)*pow2(mmst2))/(3.*
     mb*mt) - (56*mmsb1*pow2(mmt))/(3.*mb*mt) + (8*lgu*mmsb1*pow2(mmt))/(mb*mt)
     + (8*ltu*mmsb1*pow2(mmt))/(mb*mt) - (8*lgu*ltu*mmsb1*pow2(mmt))/(3.*mb*mt)
     + (56*invdgb1*mmsb12*pow2(mmt))/(3.*mb*mt) - (8*invdgb1*lgu*mmsb12*pow2(
     mmt))/(mb*mt) - (8*invdgb1*ltu*mmsb12*pow2(mmt))/(mb*mt) + (8*invdgb1*lgu*
     ltu*mmsb12*pow2(mmt))/(3.*mb*mt) + (56*mmsb2*pow2(mmt))/(3.*mb*mt) - (8*
     lgu*mmsb2*pow2(mmt))/(mb*mt) - (8*ltu*mmsb2*pow2(mmt))/(mb*mt) + (8*lgu*
     ltu*mmsb2*pow2(mmt))/(3.*mb*mt) - (56*invdgb2*mmsb22*pow2(mmt))/(3.*mb*mt)
     + (8*invdgb2*lgu*mmsb22*pow2(mmt))/(mb*mt) + (8*invdgb2*ltu*mmsb22*pow2(
     mmt))/(mb*mt) - (8*invdgb2*lgu*ltu*mmsb22*pow2(mmt))/(3.*mb*mt) + (112*
     invdgb1*mmsb1*mmst2*pow2(mmt))/(3.*mb*mt) - (8*invdgb1*lt2u*mmsb1*mmst2*
     pow2(mmt))/(mb*mt) - (8*invdgb1*lgu*lt2u*mmsb1*mmst2*pow2(mmt))/(3.*mb*mt)
     - (24*invdgb1*ltu*mmsb1*mmst2*pow2(mmt))/(mb*mt) + (8*invdgb1*lgu*ltu*
     mmsb1*mmst2*pow2(mmt))/(3.*mb*mt) + (16*invdgb1*lt2u*ltu*mmsb1*mmst2*pow2(
     mmt))/(3.*mb*mt) - (112*invdgb2*mmsb2*mmst2*pow2(mmt))/(3.*mb*mt) + (8*
     invdgb2*lt2u*mmsb2*mmst2*pow2(mmt))/(mb*mt) + (8*invdgb2*lgu*lt2u*mmsb2*
     mmst2*pow2(mmt))/(3.*mb*mt) + (24*invdgb2*ltu*mmsb2*mmst2*pow2(mmt))/(mb*
     mt) - (8*invdgb2*lgu*ltu*mmsb2*mmst2*pow2(mmt))/(3.*mb*mt) - (16*invdgb2*
     lt2u*ltu*mmsb2*mmst2*pow2(mmt))/(3.*mb*mt) - (8*mmsb1*zt2*pow2(mmt))/(3.*
     mb*mt) + (8*invdgb1*mmsb12*zt2*pow2(mmt))/(3.*mb*mt) + (8*mmsb2*zt2*pow2(
     mmt))/(3.*mb*mt) - (8*invdgb2*mmsb22*zt2*pow2(mmt))/(3.*mb*mt) + (16*
     invdgb1*mmsb1*mmst2*zt2*pow2(mmt))/(3.*mb*mt) - (16*invdgb2*mmsb2*mmst2*
     zt2*pow2(mmt))/(3.*mb*mt) - (4*mmsb1*pow2(lgu)*pow2(mmt))/(3.*mb*mt) + (4*
     invdgb1*mmsb12*pow2(lgu)*pow2(mmt))/(3.*mb*mt) + (4*mmsb2*pow2(lgu)*pow2(
     mmt))/(3.*mb*mt) - (4*invdgb2*mmsb22*pow2(lgu)*pow2(mmt))/(3.*mb*mt) + (4*
     invdgb1*mmsb1*mmst2*pow2(lt2u)*pow2(mmt))/(3.*mb*mt) - (4*invdgb2*mmsb2*
     mmst2*pow2(lt2u)*pow2(mmt))/(3.*mb*mt) - (4*mmsb1*pow2(ltu)*pow2(mmt))/(3.
     *mb*mt) + (4*invdgb1*mmsb12*pow2(ltu)*pow2(mmt))/(3.*mb*mt) + (4*mmsb2*
     pow2(ltu)*pow2(mmt))/(3.*mb*mt) - (4*invdgb2*mmsb22*pow2(ltu)*pow2(mmt))/(
     3.*mb*mt) + (4*invdgb1*mmsb1*mmst2*pow2(ltu)*pow2(mmt))/(mb*mt) - (4*
     invdgb2*mmsb2*mmst2*pow2(ltu)*pow2(mmt))/(mb*mt) - (56*invdgb1*mmt*pow3(
     mmsb1))/(3.*mb*mt) + (8*invdgb1*lgu*mmt*pow3(mmsb1))/(mb*mt) + (8*invdgb1*
     ltu*mmt*pow3(mmsb1))/(mb*mt) - (8*invdgb1*lgu*ltu*mmt*pow3(mmsb1))/(3.*mb*
     mt) - (8*invdgb1*mmt*zt2*pow3(mmsb1))/(3.*mb*mt) - (4*invdgb1*mmt*pow2(lgu
     )*pow3(mmsb1))/(3.*mb*mt) - (4*invdgb1*mmt*pow2(ltu)*pow3(mmsb1))/(3.*mb*
     mt) + (56*invdgb2*mmt*pow3(mmsb2))/(3.*mb*mt) - (8*invdgb2*lgu*mmt*pow3(
     mmsb2))/(mb*mt) - (8*invdgb2*ltu*mmt*pow3(mmsb2))/(mb*mt) + (8*invdgb2*lgu
     *ltu*mmt*pow3(mmsb2))/(3.*mb*mt) + (8*invdgb2*mmt*zt2*pow3(mmsb2))/(3.*mb*
     mt) + (4*invdgb2*mmt*pow2(lgu)*pow3(mmsb2))/(3.*mb*mt) + (4*invdgb2*mmt*
     pow2(ltu)*pow3(mmsb2))/(3.*mb*mt))) + (28*invdgb1*pow3(mmst2))/3. + (28*
     invdgb2*pow3(mmst2))/3. - 4*invdgb1*lt2u*pow3(mmst2) - 4*invdgb2*lt2u*pow3
     (mmst2) - 4*invdgb1*ltu*pow3(mmst2) - 4*invdgb2*ltu*pow3(mmst2) + (4*
     invdgb1*lt2u*ltu*pow3(mmst2))/3. + (4*invdgb2*lt2u*ltu*pow3(mmst2))/3. + (
     4*invdgb1*zt2*pow3(mmst2))/3. + (4*invdgb2*zt2*pow3(mmst2))/3. - (28*mmsb1
     *pow2(invdgb1)*pow3(mmst2))/3. + 4*lt2u*mmsb1*pow2(invdgb1)*pow3(mmst2) +
     4*ltu*mmsb1*pow2(invdgb1)*pow3(mmst2) - (4*lt2u*ltu*mmsb1*pow2(invdgb1)*
     pow3(mmst2))/3. - (4*mmsb1*zt2*pow2(invdgb1)*pow3(mmst2))/3. - (28*mmsb2*
     pow2(invdgb2)*pow3(mmst2))/3. + 4*lt2u*mmsb2*pow2(invdgb2)*pow3(mmst2) + 4
     *ltu*mmsb2*pow2(invdgb2)*pow3(mmst2) - (4*lt2u*ltu*mmsb2*pow2(invdgb2)*
     pow3(mmst2))/3. - (4*mmsb2*zt2*pow2(invdgb2)*pow3(mmst2))/3. + (2*invdgb1*
     pow2(lt2u)*pow3(mmst2))/3. + (2*invdgb2*pow2(lt2u)*pow3(mmst2))/3. - (2*
     mmsb1*pow2(invdgb1)*pow2(lt2u)*pow3(mmst2))/3. - (2*mmsb2*pow2(invdgb2)*
     pow2(lt2u)*pow3(mmst2))/3. + (2*invdgb1*pow2(ltu)*pow3(mmst2))/3. + (2*
     invdgb2*pow2(ltu)*pow3(mmst2))/3. - (2*mmsb1*pow2(invdgb1)*pow2(ltu)*pow3(
     mmst2))/3. - (2*mmsb2*pow2(invdgb2)*pow2(ltu)*pow3(mmst2))/3. + s2b*((56*
     mmgl2*mmsb1)/(3.*mb*mgl) - (8*lgu*mmgl2*mmsb1)/(mb*mgl) - (8*ltu*mmgl2*
     mmsb1)/(mb*mgl) + (8*lgu*ltu*mmgl2*mmsb1)/(3.*mb*mgl) + (56*mmgl*mmsb12)/(
     3.*mb*mgl) - (8*lgu*mmgl*mmsb12)/(mb*mgl) - (8*ltu*mmgl*mmsb12)/(mb*mgl) +
     (8*lgu*ltu*mmgl*mmsb12)/(3.*mb*mgl) - (56*mmgl2*mmsb2)/(3.*mb*mgl) + (8*
     lgu*mmgl2*mmsb2)/(mb*mgl) + (8*ltu*mmgl2*mmsb2)/(mb*mgl) - (8*lgu*ltu*
     mmgl2*mmsb2)/(3.*mb*mgl) - (56*mmgl*mmsb22)/(3.*mb*mgl) + (8*lgu*mmgl*
     mmsb22)/(mb*mgl) + (8*ltu*mmgl*mmsb22)/(mb*mgl) - (8*lgu*ltu*mmgl*mmsb22)/
     (3.*mb*mgl) - (56*mmgl*mmsb1*mmst2)/(3.*mb*mgl) + (16*lgu*mmgl*mmsb1*mmst2
     )/(mb*mgl) - (8*lt2u*mmgl*mmsb1*mmst2)/(mb*mgl) + (8*ltu*mmgl*mmsb1*mmst2)
     /(mb*mgl) - (16*lgu*ltu*mmgl*mmsb1*mmst2)/(3.*mb*mgl) + (8*lt2u*ltu*mmgl*
     mmsb1*mmst2)/(3.*mb*mgl) + (56*invdgb1*mmgl*mmsb12*mmst2)/(3.*mb*mgl) - (
     16*invdgb1*lgu*mmgl*mmsb12*mmst2)/(mb*mgl) + (8*invdgb1*lt2u*mmgl*mmsb12*
     mmst2)/(mb*mgl) - (8*invdgb1*ltu*mmgl*mmsb12*mmst2)/(mb*mgl) + (16*invdgb1
     *lgu*ltu*mmgl*mmsb12*mmst2)/(3.*mb*mgl) - (8*invdgb1*lt2u*ltu*mmgl*mmsb12*
     mmst2)/(3.*mb*mgl) + (56*mmgl*mmsb2*mmst2)/(3.*mb*mgl) - (16*lgu*mmgl*
     mmsb2*mmst2)/(mb*mgl) + (8*lt2u*mmgl*mmsb2*mmst2)/(mb*mgl) - (8*ltu*mmgl*
     mmsb2*mmst2)/(mb*mgl) + (16*lgu*ltu*mmgl*mmsb2*mmst2)/(3.*mb*mgl) - (8*
     lt2u*ltu*mmgl*mmsb2*mmst2)/(3.*mb*mgl) - (56*invdgb2*mmgl*mmsb22*mmst2)/(
     3.*mb*mgl) + (16*invdgb2*lgu*mmgl*mmsb22*mmst2)/(mb*mgl) - (8*invdgb2*lt2u
     *mmgl*mmsb22*mmst2)/(mb*mgl) + (8*invdgb2*ltu*mmgl*mmsb22*mmst2)/(mb*mgl)
     - (16*invdgb2*lgu*ltu*mmgl*mmsb22*mmst2)/(3.*mb*mgl) + (8*invdgb2*lt2u*ltu
     *mmgl*mmsb22*mmst2)/(3.*mb*mgl) - (56*mmgl*mmsb1*mmt)/(3.*mb*mgl) + (8*lgu
     *mmgl*mmsb1*mmt)/(mb*mgl) + (8*ltu*mmgl*mmsb1*mmt)/(mb*mgl) - (8*lgu*ltu*
     mmgl*mmsb1*mmt)/(3.*mb*mgl) + (56*invdgb1*mmgl*mmsb12*mmt)/(3.*mb*mgl) - (
     8*invdgb1*lgu*mmgl*mmsb12*mmt)/(mb*mgl) - (8*invdgb1*ltu*mmgl*mmsb12*mmt)/
     (mb*mgl) + (8*invdgb1*lgu*ltu*mmgl*mmsb12*mmt)/(3.*mb*mgl) + (56*mmgl*
     mmsb2*mmt)/(3.*mb*mgl) - (8*lgu*mmgl*mmsb2*mmt)/(mb*mgl) - (8*ltu*mmgl*
     mmsb2*mmt)/(mb*mgl) + (8*lgu*ltu*mmgl*mmsb2*mmt)/(3.*mb*mgl) - (56*invdgb2
     *mmgl*mmsb22*mmt)/(3.*mb*mgl) + (8*invdgb2*lgu*mmgl*mmsb22*mmt)/(mb*mgl) +
     (8*invdgb2*ltu*mmgl*mmsb22*mmt)/(mb*mgl) - (8*invdgb2*lgu*ltu*mmgl*mmsb22*
     mmt)/(3.*mb*mgl) + (224*invdgb1*mmgl*mmsb1*mmst2*mmt)/(3.*mb*mgl) - (8*
     invdgb1*lgu*mmgl*mmsb1*mmst2*mmt)/(mb*mgl) - (8*invdgb1*lt2u*mmgl*mmsb1*
     mmst2*mmt)/(mb*mgl) - (16*invdgb1*lgu*lt2u*mmgl*mmsb1*mmst2*mmt)/(3.*mb*
     mgl) - (48*invdgb1*ltu*mmgl*mmsb1*mmst2*mmt)/(mb*mgl) + (8*invdgb1*lgu*ltu
     *mmgl*mmsb1*mmst2*mmt)/(mb*mgl) + (8*invdgb1*lt2u*ltu*mmgl*mmsb1*mmst2*mmt
     )/(mb*mgl) - (224*invdgb2*mmgl*mmsb2*mmst2*mmt)/(3.*mb*mgl) + (8*invdgb2*
     lgu*mmgl*mmsb2*mmst2*mmt)/(mb*mgl) + (8*invdgb2*lt2u*mmgl*mmsb2*mmst2*mmt)
     /(mb*mgl) + (16*invdgb2*lgu*lt2u*mmgl*mmsb2*mmst2*mmt)/(3.*mb*mgl) + (48*
     invdgb2*ltu*mmgl*mmsb2*mmst2*mmt)/(mb*mgl) - (8*invdgb2*lgu*ltu*mmgl*mmsb2
     *mmst2*mmt)/(mb*mgl) - (8*invdgb2*lt2u*ltu*mmgl*mmsb2*mmst2*mmt)/(mb*mgl)
     + (8*mmgl2*mmsb1*zt2)/(3.*mb*mgl) + (8*mmgl*mmsb12*zt2)/(3.*mb*mgl) - (8*
     mmgl2*mmsb2*zt2)/(3.*mb*mgl) - (8*mmgl*mmsb22*zt2)/(3.*mb*mgl) - (8*mmgl*
     mmsb1*mmst2*zt2)/(3.*mb*mgl) + (8*invdgb1*mmgl*mmsb12*mmst2*zt2)/(3.*mb*
     mgl) + (8*mmgl*mmsb2*mmst2*zt2)/(3.*mb*mgl) - (8*invdgb2*mmgl*mmsb22*mmst2
     *zt2)/(3.*mb*mgl) - (8*mmgl*mmsb1*mmt*zt2)/(3.*mb*mgl) + (8*invdgb1*mmgl*
     mmsb12*mmt*zt2)/(3.*mb*mgl) + (8*mmgl*mmsb2*mmt*zt2)/(3.*mb*mgl) - (8*
     invdgb2*mmgl*mmsb22*mmt*zt2)/(3.*mb*mgl) + (32*invdgb1*mmgl*mmsb1*mmst2*
     mmt*zt2)/(3.*mb*mgl) - (32*invdgb2*mmgl*mmsb2*mmst2*mmt*zt2)/(3.*mb*mgl) +
     (4*mmgl2*mmsb1*pow2(lgu))/(3.*mb*mgl) + (4*mmgl*mmsb12*pow2(lgu))/(3.*mb*
     mgl) - (4*mmgl2*mmsb2*pow2(lgu))/(3.*mb*mgl) - (4*mmgl*mmsb22*pow2(lgu))/(
     3.*mb*mgl) - (8*mmgl*mmsb1*mmst2*pow2(lgu))/(3.*mb*mgl) + (8*invdgb1*mmgl*
     mmsb12*mmst2*pow2(lgu))/(3.*mb*mgl) + (8*mmgl*mmsb2*mmst2*pow2(lgu))/(3.*
     mb*mgl) - (8*invdgb2*mmgl*mmsb22*mmst2*pow2(lgu))/(3.*mb*mgl) - (4*mmgl*
     mmsb1*mmt*pow2(lgu))/(3.*mb*mgl) + (4*invdgb1*mmgl*mmsb12*mmt*pow2(lgu))/(
     3.*mb*mgl) + (4*mmgl*mmsb2*mmt*pow2(lgu))/(3.*mb*mgl) - (4*invdgb2*mmgl*
     mmsb22*mmt*pow2(lgu))/(3.*mb*mgl) + (4*invdgb1*mmgl*mmsb1*mmst2*mmt*pow2(
     lgu))/(3.*mb*mgl) - (4*invdgb2*mmgl*mmsb2*mmst2*mmt*pow2(lgu))/(3.*mb*mgl)
     + (4*mmgl*mmsb1*mmst2*pow2(lt2u))/(3.*mb*mgl) - (4*invdgb1*mmgl*mmsb12*
     mmst2*pow2(lt2u))/(3.*mb*mgl) - (4*mmgl*mmsb2*mmst2*pow2(lt2u))/(3.*mb*mgl
     ) + (4*invdgb2*mmgl*mmsb22*mmst2*pow2(lt2u))/(3.*mb*mgl) + (4*invdgb1*mmgl
     *mmsb1*mmst2*mmt*pow2(lt2u))/(3.*mb*mgl) - (4*invdgb2*mmgl*mmsb2*mmst2*mmt
     *pow2(lt2u))/(3.*mb*mgl) + (4*mmgl2*mmsb1*pow2(ltu))/(3.*mb*mgl) + (4*mmgl
     *mmsb12*pow2(ltu))/(3.*mb*mgl) - (4*mmgl2*mmsb2*pow2(ltu))/(3.*mb*mgl) - (
     4*mmgl*mmsb22*pow2(ltu))/(3.*mb*mgl) - (4*mmgl*mmsb1*mmst2*pow2(ltu))/(3.*
     mb*mgl) + (4*invdgb1*mmgl*mmsb12*mmst2*pow2(ltu))/(3.*mb*mgl) + (4*mmgl*
     mmsb2*mmst2*pow2(ltu))/(3.*mb*mgl) - (4*invdgb2*mmgl*mmsb22*mmst2*pow2(ltu
     ))/(3.*mb*mgl) - (4*mmgl*mmsb1*mmt*pow2(ltu))/(3.*mb*mgl) + (4*invdgb1*
     mmgl*mmsb12*mmt*pow2(ltu))/(3.*mb*mgl) + (4*mmgl*mmsb2*mmt*pow2(ltu))/(3.*
     mb*mgl) - (4*invdgb2*mmgl*mmsb22*mmt*pow2(ltu))/(3.*mb*mgl) + (8*invdgb1*
     mmgl*mmsb1*mmst2*mmt*pow2(ltu))/(mb*mgl) - (8*invdgb2*mmgl*mmsb2*mmst2*mmt
     *pow2(ltu))/(mb*mgl) + (56*invdgb1*mmgl*mmsb1*pow2(mmst2))/(3.*mb*mgl) + (
     8*invdgb1*lgu*mmgl*mmsb1*pow2(mmst2))/(mb*mgl) - (16*invdgb1*lt2u*mmgl*
     mmsb1*pow2(mmst2))/(mb*mgl) - (8*invdgb1*ltu*mmgl*mmsb1*pow2(mmst2))/(mb*
     mgl) - (8*invdgb1*lgu*ltu*mmgl*mmsb1*pow2(mmst2))/(3.*mb*mgl) + (16*
     invdgb1*lt2u*ltu*mmgl*mmsb1*pow2(mmst2))/(3.*mb*mgl) - (56*invdgb2*mmgl*
     mmsb2*pow2(mmst2))/(3.*mb*mgl) - (8*invdgb2*lgu*mmgl*mmsb2*pow2(mmst2))/(
     mb*mgl) + (16*invdgb2*lt2u*mmgl*mmsb2*pow2(mmst2))/(mb*mgl) + (8*invdgb2*
     ltu*mmgl*mmsb2*pow2(mmst2))/(mb*mgl) + (8*invdgb2*lgu*ltu*mmgl*mmsb2*pow2(
     mmst2))/(3.*mb*mgl) - (16*invdgb2*lt2u*ltu*mmgl*mmsb2*pow2(mmst2))/(3.*mb*
     mgl) + (56*invdgb1*mmgl*mmt*pow2(mmst2))/(3.*mb*mgl) - (56*invdgb2*mmgl*
     mmt*pow2(mmst2))/(3.*mb*mgl) - (8*invdgb1*lt2u*mmgl*mmt*pow2(mmst2))/(mb*
     mgl) + (8*invdgb2*lt2u*mmgl*mmt*pow2(mmst2))/(mb*mgl) - (8*invdgb1*ltu*
     mmgl*mmt*pow2(mmst2))/(mb*mgl) + (8*invdgb2*ltu*mmgl*mmt*pow2(mmst2))/(mb*
     mgl) + (8*invdgb1*lt2u*ltu*mmgl*mmt*pow2(mmst2))/(3.*mb*mgl) - (8*invdgb2*
     lt2u*ltu*mmgl*mmt*pow2(mmst2))/(3.*mb*mgl) + (8*invdgb1*mmgl*mmsb1*zt2*
     pow2(mmst2))/(3.*mb*mgl) - (8*invdgb2*mmgl*mmsb2*zt2*pow2(mmst2))/(3.*mb*
     mgl) + (8*invdgb1*mmgl*mmt*zt2*pow2(mmst2))/(3.*mb*mgl) - (8*invdgb2*mmgl*
     mmt*zt2*pow2(mmst2))/(3.*mb*mgl) - (4*invdgb1*mmgl*mmsb1*pow2(lgu)*pow2(
     mmst2))/(3.*mb*mgl) + (4*invdgb2*mmgl*mmsb2*pow2(lgu)*pow2(mmst2))/(3.*mb*
     mgl) + (8*invdgb1*mmgl*mmsb1*pow2(lt2u)*pow2(mmst2))/(3.*mb*mgl) - (8*
     invdgb2*mmgl*mmsb2*pow2(lt2u)*pow2(mmst2))/(3.*mb*mgl) + (4*invdgb1*mmgl*
     mmt*pow2(lt2u)*pow2(mmst2))/(3.*mb*mgl) - (4*invdgb2*mmgl*mmt*pow2(lt2u)*
     pow2(mmst2))/(3.*mb*mgl) + (4*invdgb1*mmgl*mmsb1*pow2(ltu)*pow2(mmst2))/(
     3.*mb*mgl) - (4*invdgb2*mmgl*mmsb2*pow2(ltu)*pow2(mmst2))/(3.*mb*mgl) + (4
     *invdgb1*mmgl*mmt*pow2(ltu)*pow2(mmst2))/(3.*mb*mgl) - (4*invdgb2*mmgl*mmt
     *pow2(ltu)*pow2(mmst2))/(3.*mb*mgl) - (56*invdgb1*mmgl*pow3(mmsb1))/(3.*mb
     *mgl) + (8*invdgb1*lgu*mmgl*pow3(mmsb1))/(mb*mgl) + (8*invdgb1*ltu*mmgl*
     pow3(mmsb1))/(mb*mgl) - (8*invdgb1*lgu*ltu*mmgl*pow3(mmsb1))/(3.*mb*mgl) -
     (8*invdgb1*mmgl*zt2*pow3(mmsb1))/(3.*mb*mgl) - (4*invdgb1*mmgl*pow2(lgu)*
     pow3(mmsb1))/(3.*mb*mgl) - (4*invdgb1*mmgl*pow2(ltu)*pow3(mmsb1))/(3.*mb*
     mgl) + (56*invdgb2*mmgl*pow3(mmsb2))/(3.*mb*mgl) - (8*invdgb2*lgu*mmgl*
     pow3(mmsb2))/(mb*mgl) - (8*invdgb2*ltu*mmgl*pow3(mmsb2))/(mb*mgl) + (8*
     invdgb2*lgu*ltu*mmgl*pow3(mmsb2))/(3.*mb*mgl) + (8*invdgb2*mmgl*zt2*pow3(
     mmsb2))/(3.*mb*mgl) + (4*invdgb2*mmgl*pow2(lgu)*pow3(mmsb2))/(3.*mb*mgl) +
     (4*invdgb2*mmgl*pow2(ltu)*pow3(mmsb2))/(3.*mb*mgl) - (56*invdgb1*mmgl*pow3
     (mmst2))/(3.*mb*mgl) + (56*invdgb2*mmgl*pow3(mmst2))/(3.*mb*mgl) + (8*
     invdgb1*lt2u*mmgl*pow3(mmst2))/(mb*mgl) - (8*invdgb2*lt2u*mmgl*pow3(mmst2)
     )/(mb*mgl) + (8*invdgb1*ltu*mmgl*pow3(mmst2))/(mb*mgl) - (8*invdgb2*ltu*
     mmgl*pow3(mmst2))/(mb*mgl) - (8*invdgb1*lt2u*ltu*mmgl*pow3(mmst2))/(3.*mb*
     mgl) + (8*invdgb2*lt2u*ltu*mmgl*pow3(mmst2))/(3.*mb*mgl) - (8*invdgb1*mmgl
     *zt2*pow3(mmst2))/(3.*mb*mgl) + (8*invdgb2*mmgl*zt2*pow3(mmst2))/(3.*mb*
     mgl) - (4*invdgb1*mmgl*pow2(lt2u)*pow3(mmst2))/(3.*mb*mgl) + (4*invdgb2*
     mmgl*pow2(lt2u)*pow3(mmst2))/(3.*mb*mgl) - (4*invdgb1*mmgl*pow2(ltu)*pow3(
     mmst2))/(3.*mb*mgl) + (4*invdgb2*mmgl*pow2(ltu)*pow3(mmst2))/(3.*mb*mgl))
     - (28*pow2(invdgb1)*pow4(mmsb1))/3. + 4*lgu*pow2(invdgb1)*pow4(mmsb1) + 4*
     ltu*pow2(invdgb1)*pow4(mmsb1) - (4*lgu*ltu*pow2(invdgb1)*pow4(mmsb1))/3. -
     (4*zt2*pow2(invdgb1)*pow4(mmsb1))/3. - (2*pow2(invdgb1)*pow2(lgu)*pow4(
     mmsb1))/3. - (2*pow2(invdgb1)*pow2(ltu)*pow4(mmsb1))/3. - (28*pow2(invdgb2
     )*pow4(mmsb2))/3. + 4*lgu*pow2(invdgb2)*pow4(mmsb2) + 4*ltu*pow2(invdgb2)*
     pow4(mmsb2) - (4*lgu*ltu*pow2(invdgb2)*pow4(mmsb2))/3. - (4*zt2*pow2(
     invdgb2)*pow4(mmsb2))/3. - (2*pow2(invdgb2)*pow2(lgu)*pow4(mmsb2))/3. - (2
     *pow2(invdgb2)*pow2(ltu)*pow4(mmsb2))/3.) + (16*lb1u*lb2u*((-5*invdgb1*
     mmsb2)/2. + (5*invdgb2*mmsb2)/2. - 2*invdb12*invdgb1*mmsb22 + 2*invdb12*
     invdgb2*mmsb22 + (mmsb1*mmsb2*pow2(invdgb1))/2. + (3*mmsb22*pow2(invdgb1))
     /4. - (mmsb22*pow2(invdgb2))/4. + invdgb2*mmsb1*pow2(s2b) + 3*invdgb1*
     mmsb2*pow2(s2b) - 2*invdgb2*mmsb2*pow2(s2b) + (5*invdb12*invdgb1*mmsb22*
     pow2(s2b))/2. - (5*invdb12*invdgb2*mmsb22*pow2(s2b))/2. - (3*mmsb1*mmsb2*
     pow2(invdgb1)*pow2(s2b))/2. + (3*mmsb22*pow2(invdgb1)*pow2(s2b))/4. - 2*
     mmsb1*mmsb2*pow2(invdgb2)*pow2(s2b) - (mmsb22*pow2(invdgb2)*pow2(s2b))/4.
     + mmsb12*mmsb2*pow2(s2b)*pow3(invdgb1) + mmsb1*mmsb22*pow2(s2b)*pow3(
     invdgb2) + invdgb1*pow2(invdb12)*pow3(mmsb2) - invdgb2*pow2(invdb12)*pow3(
     mmsb2) + invdb12*pow2(invdgb1)*pow3(mmsb2) + invdgb1*pow2(invdb12)*pow2(
     s2b)*pow3(mmsb2) - invdgb2*pow2(invdb12)*pow2(s2b)*pow3(mmsb2) + invdb12*
     pow2(invdgb1)*pow2(s2b)*pow3(mmsb2) + (pow2(invdb12)*pow2(invdgb1)*pow4(
     mmsb2))/4. + (pow2(invdb12)*pow2(invdgb2)*pow4(mmsb2))/4. + (pow2(invdb12)
     *pow2(invdgb1)*pow2(s2b)*pow4(mmsb2))/4. + (pow2(invdb12)*pow2(invdgb2)*
     pow2(s2b)*pow4(mmsb2))/4. + (invdgb1*pow3(invdb12)*pow4(mmsb2))/2. - (
     invdgb2*pow3(invdb12)*pow4(mmsb2))/2. + (invdgb1*pow2(s2b)*pow3(invdb12)*
     pow4(mmsb2))/2. - (invdgb2*pow2(s2b)*pow3(invdb12)*pow4(mmsb2))/2.))/9. +
     (16*lb1u*lgu*(4 + (invdgb1*mmsb1)/2. + (5*invdgb1*mmsb2)/2. - (5*invdgb2*
     mmsb2)/2. + 2*invdb12*invdgb1*mmsb22 - 2*invdb12*invdgb2*mmsb22 - (7*
     mmsb12*pow2(invdgb1))/2. - (mmsb1*mmsb2*pow2(invdgb1))/2. - (3*mmsb22*pow2
     (invdgb1))/4. + (mmsb22*pow2(invdgb2))/4. + 4*pow2(s2b) - 7*invdgb1*mmsb1*
     pow2(s2b) - invdgb2*mmsb1*pow2(s2b) - invdgb1*mmsb2*pow2(s2b) + invdgb2*
     mmsb2*pow2(s2b) - (3*invdb12*invdgb1*mmsb22*pow2(s2b))/2. + (3*invdb12*
     invdgb2*mmsb22*pow2(s2b))/2. + (3*mmsb12*pow2(invdgb1)*pow2(s2b))/2. - (
     mmsb1*mmsb2*pow2(invdgb1)*pow2(s2b))/2. - (3*mmsb22*pow2(invdgb1)*pow2(s2b
     ))/4. + 2*mmsb1*mmsb2*pow2(invdgb2)*pow2(s2b) + (mmsb22*pow2(invdgb2)*pow2
     (s2b))/4. - mmsb1*mmsb22*pow2(s2b)*pow3(invdgb2) + 2*pow3(invdgb1)*pow3(
     mmsb1) - invdgb1*pow2(invdb12)*pow3(mmsb2) + invdgb2*pow2(invdb12)*pow3(
     mmsb2) - invdb12*pow2(invdgb1)*pow3(mmsb2) - invdgb1*pow2(invdb12)*pow2(
     s2b)*pow3(mmsb2) + invdgb2*pow2(invdb12)*pow2(s2b)*pow3(mmsb2) - invdb12*
     pow2(invdgb1)*pow2(s2b)*pow3(mmsb2) - (3*pow4(invdgb1)*pow4(mmsb1))/4. + (
     pow2(s2b)*pow4(invdgb1)*pow4(mmsb1))/4. - (pow2(invdb12)*pow2(invdgb1)*
     pow4(mmsb2))/4. - (pow2(invdb12)*pow2(invdgb2)*pow4(mmsb2))/4. - (pow2(
     invdb12)*pow2(invdgb1)*pow2(s2b)*pow4(mmsb2))/4. - (pow2(invdb12)*pow2(
     invdgb2)*pow2(s2b)*pow4(mmsb2))/4. - (invdgb1*pow3(invdb12)*pow4(mmsb2))/
     2. + (invdgb2*pow3(invdb12)*pow4(mmsb2))/2. - (invdgb1*pow2(s2b)*pow3(
     invdb12)*pow4(mmsb2))/2. + (invdgb2*pow2(s2b)*pow3(invdb12)*pow4(mmsb2))/
     2.))/9. + (16*pow2(lb2u)*(4*invdb12*mmsb2 - (25*invdgb2*mmsb2)/4. - 4*
     invdb12*invdgb2*mmsb22 + (15*mmsb22*pow2(invdgb2))/2. - 4*invdb12*mmsb2*
     pow2(s2b) + (invdgb1*mmsb2*pow2(s2b))/2. - (5*invdgb2*mmsb2*pow2(s2b))/2.
     + (invdb12*invdgb1*mmsb22*pow2(s2b))/2. + (7*invdb12*invdgb2*mmsb22*pow2(
     s2b))/2. + (5*mmsb22*pow2(invdgb2)*pow2(s2b))/4. - 6*pow3(invdgb2)*pow3(
     mmsb2) - (pow2(s2b)*pow3(invdgb2)*pow3(mmsb2))/2. + (3*pow4(invdgb2)*pow4(
     mmsb2))/8. - (pow2(s2b)*pow4(invdgb2)*pow4(mmsb2))/8.))/9. + (16*pow2(lgu)
     *(-18.5 + (83*invdgb1*mmsb1)/4. + (35*invdgb1*mmsb2)/2. + (13*invdgb2*
     mmsb2)/4. + 18*invdb12*invdgb1*mmsb22 - 18*invdb12*invdgb2*mmsb22 - mmsb12
     *pow2(invdgb1) + (mmsb1*mmsb2*pow2(invdgb1))/2. + (3*mmsb22*pow2(invdgb1))
     /4. - (5*mmsb22*pow2(invdgb2))/4. + 4*pow2(s2b) - 4*invdgb1*mmsb1*pow2(s2b
     ) - 9*invdgb1*mmsb2*pow2(s2b) + 5*invdgb2*mmsb2*pow2(s2b) - (17*invdb12*
     invdgb1*mmsb22*pow2(s2b))/2. + (17*invdb12*invdgb2*mmsb22*pow2(s2b))/2. -
     (11*mmsb12*pow2(invdgb1)*pow2(s2b))/4. + (mmsb1*mmsb2*pow2(invdgb1)*pow2(
     s2b))/2. + (3*mmsb22*pow2(invdgb1)*pow2(s2b))/4. - 3*mmsb22*pow2(invdgb2)*
     pow2(s2b) - 2*pow3(invdgb1)*pow3(mmsb1) + (pow2(s2b)*pow3(invdgb1)*pow3(
     mmsb1))/2. + invdgb1*pow2(invdb12)*pow3(mmsb2) - invdgb2*pow2(invdb12)*
     pow3(mmsb2) + invdb12*pow2(invdgb1)*pow3(mmsb2) + invdgb1*pow2(invdb12)*
     pow2(s2b)*pow3(mmsb2) - invdgb2*pow2(invdb12)*pow2(s2b)*pow3(mmsb2) +
     invdb12*pow2(invdgb1)*pow2(s2b)*pow3(mmsb2) - 2*pow3(invdgb2)*pow3(mmsb2)
     + (pow2(s2b)*pow3(invdgb2)*pow3(mmsb2))/2. + (3*pow4(invdgb1)*pow4(mmsb1))
     /8. - (pow2(s2b)*pow4(invdgb1)*pow4(mmsb1))/8. + (pow2(invdb12)*pow2(
     invdgb1)*pow4(mmsb2))/4. + (pow2(invdb12)*pow2(invdgb2)*pow4(mmsb2))/4. +
     (pow2(invdb12)*pow2(invdgb1)*pow2(s2b)*pow4(mmsb2))/4. + (pow2(invdb12)*
     pow2(invdgb2)*pow2(s2b)*pow4(mmsb2))/4. + (invdgb1*pow3(invdb12)*pow4(
     mmsb2))/2. - (invdgb2*pow3(invdb12)*pow4(mmsb2))/2. + (invdgb1*pow2(s2b)*
     pow3(invdb12)*pow4(mmsb2))/2. - (invdgb2*pow2(s2b)*pow3(invdb12)*pow4(
     mmsb2))/2. + (3*pow4(invdgb2)*pow4(mmsb2))/8. - (pow2(s2b)*pow4(invdgb2)*
     pow4(mmsb2))/8.))/9. + (16*lb2u*lgu*(4 + (5*invdgb1*mmsb2)/2. - 2*invdgb2*
     mmsb2 + 2*invdb12*invdgb1*mmsb22 - 2*invdb12*invdgb2*mmsb22 - (mmsb1*mmsb2
     *pow2(invdgb1))/2. - (3*mmsb22*pow2(invdgb1))/4. - (13*mmsb22*pow2(invdgb2
     ))/4. + 4*pow2(s2b) - 2*invdgb1*mmsb2*pow2(s2b) - 6*invdgb2*mmsb2*pow2(s2b
     ) - (3*invdb12*invdgb1*mmsb22*pow2(s2b))/2. + (3*invdb12*invdgb2*mmsb22*
     pow2(s2b))/2. + (3*mmsb1*mmsb2*pow2(invdgb1)*pow2(s2b))/2. - (3*mmsb22*
     pow2(invdgb1)*pow2(s2b))/4. + (7*mmsb22*pow2(invdgb2)*pow2(s2b))/4. -
     mmsb12*mmsb2*pow2(s2b)*pow3(invdgb1) - invdgb1*pow2(invdb12)*pow3(mmsb2) +
     invdgb2*pow2(invdb12)*pow3(mmsb2) - invdb12*pow2(invdgb1)*pow3(mmsb2) -
     invdgb1*pow2(invdb12)*pow2(s2b)*pow3(mmsb2) + invdgb2*pow2(invdb12)*pow2(
     s2b)*pow3(mmsb2) - invdb12*pow2(invdgb1)*pow2(s2b)*pow3(mmsb2) + 2*pow3(
     invdgb2)*pow3(mmsb2) - (pow2(invdb12)*pow2(invdgb1)*pow4(mmsb2))/4. - (
     pow2(invdb12)*pow2(invdgb2)*pow4(mmsb2))/4. - (pow2(invdb12)*pow2(invdgb1)
     *pow2(s2b)*pow4(mmsb2))/4. - (pow2(invdb12)*pow2(invdgb2)*pow2(s2b)*pow4(
     mmsb2))/4. - (invdgb1*pow3(invdb12)*pow4(mmsb2))/2. + (invdgb2*pow3(
     invdb12)*pow4(mmsb2))/2. - (invdgb1*pow2(s2b)*pow3(invdb12)*pow4(mmsb2))/
     2. + (invdgb2*pow2(s2b)*pow3(invdb12)*pow4(mmsb2))/2. - (3*pow4(invdgb2)*
     pow4(mmsb2))/4. + (pow2(s2b)*pow4(invdgb2)*pow4(mmsb2))/4.))/9.;

   return pow4(g3) * result * twoLoop;
}

} // namespace mssm_twoloop_mb
} // namespace flexiblesusy
