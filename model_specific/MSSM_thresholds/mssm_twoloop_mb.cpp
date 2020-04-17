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

// This file has been generated at Sat 18 Apr 2020 00:18:47
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
   const double mmt    = pow2(pars.mt);
   const double mmb    = pow2(pars.mb);
   const double mmgl   = pow2(pars.mg);
   const double mmu    = pow2(pars.Q);
   const double mmst1  = pow2(pars.mst1);
   const double mmst2  = pow2(pars.mst2);
   const double mmsb1  = pow2(pars.msb1);
   const double mmsb2  = pow2(pars.msb2);
   const double mmsusy = pow2(pars.msusy);
   const double lgu    = std::log(mmgl/mmu);
   const double lt1u   = std::log(mmst1/mmu);
   const double lt2u   = std::log(mmst2/mmu);
   const double lb1u   = std::log(mmsb1/mmu);
   const double lb2u   = std::log(mmsb2/mmu);
   const double lsu    = std::log(mmsusy/mmu);
   const double ltu    = std::log(mmt/mmu);

   const double result =
   (Xb*((128*lb1u*mmgl*pow4(g3))/(9.*mgl*(mmsb1 - mmsb2)) - (128*lb2u*mmgl*pow4
     (g3))/(9.*mgl*(mmsb1 - mmsb2)) - (128*lb1u*lgu*mmgl*pow4(g3))/(9.*mgl*(
     mmsb1 - mmsb2)) + (128*lb2u*lgu*mmgl*pow4(g3))/(9.*mgl*(mmsb1 - mmsb2)) -
     (72*mmgl*mmsb1*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (112*lb1u
     *mmgl*mmsb1*pow4(g3))/(9.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (96*lgu*
     mmgl*mmsb1*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (200*lb1u*lgu
     *mmgl*mmsb1*pow4(g3))/(9.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (32*ltu*
     mmgl*mmsb1*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (32*lgu*ltu*
     mmgl*mmsb1*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*
     mmsb2*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*lb1u*mmgl*
     mmsb2*pow4(g3))/(9.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (80*lb2u*mmgl*
     mmsb2*pow4(g3))/(9.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*lb1u*lb2u*
     mmgl*mmsb2*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*lgu*
     mmgl*mmsb2*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (80*lb1u*
     lgu*mmgl*mmsb2*pow4(g3))/(9.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (8*
     lb2u*lgu*mmgl*mmsb2*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (
     16*mmgl*mmsb1*pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (32*lb1u*
     mmgl*mmsb1*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*lgu*
     mmgl*mmsb1*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (8*lb1u*
     lgu*mmgl*mmsb1*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (72*
     mmgl*mmsb2*pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*lb1u*mmgl
     *mmsb2*pow4(g3))/(9.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (128*lb2u*mmgl
     *mmsb2*pow4(g3))/(9.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (80*lb1u*lb2u*
     mmgl*mmsb2*pow4(g3))/(9.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (96*lgu*
     mmgl*mmsb2*pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (80*lb1u*lgu*
     mmgl*mmsb2*pow4(g3))/(9.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (152*lb2u*
     lgu*mmgl*mmsb2*pow4(g3))/(9.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (32*
     ltu*mmgl*mmsb2*pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (32*lgu*
     ltu*mmgl*mmsb2*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (64*
     mmgl*mmst1*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*lgu*
     mmgl*mmst1*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*lt1u*
     mmgl*mmst1*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (8*lgu*
     lt1u*mmgl*mmst1*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*
     ltu*mmgl*mmst1*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*lt1u*
     ltu*mmgl*mmst1*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (64*
     mmgl*mmst1*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*lgu*
     mmgl*mmst1*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*lt1u*
     mmgl*mmst1*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (8*lgu*
     lt1u*mmgl*mmst1*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*
     ltu*mmgl*mmst1*pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*lt1u*
     ltu*mmgl*mmst1*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (64*
     mmgl*mmst2*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*lgu*
     mmgl*mmst2*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*lt2u*
     mmgl*mmst2*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (8*lgu*
     lt2u*mmgl*mmst2*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*
     ltu*mmgl*mmst2*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*lt2u*
     ltu*mmgl*mmst2*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (64*
     mmgl*mmst2*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*lgu*
     mmgl*mmst2*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*lt2u*
     mmgl*mmst2*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (8*lgu*
     lt2u*mmgl*mmst2*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*
     ltu*mmgl*mmst2*pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*lt2u*
     ltu*mmgl*mmst2*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (128*
     mmgl*mmsusy*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (128*lgu*
     mmgl*mmsusy*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (256*lsu*
     mmgl*mmsusy*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (64*lgu*
     lsu*mmgl*mmsusy*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (128*
     mmgl*mmsusy*pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (128*lgu*
     mmgl*mmsusy*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (256*lsu*
     mmgl*mmsusy*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (64*lgu*
     lsu*mmgl*mmsusy*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (32*
     mmgl*mmt*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (32*lgu*mmgl*
     mmt*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (64*ltu*mmgl*mmt*
     pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*lgu*ltu*mmgl*mmt*
     pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (32*mmgl*mmt*pow4(g3)
     )/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (32*lgu*mmgl*mmt*pow4(g3))/(3.*
     mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (64*ltu*mmgl*mmt*pow4(g3))/(3.*mgl*
     (mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*lgu*ltu*mmgl*mmt*pow4(g3))/(3.*mgl*
     (mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (80*mmgl*mmsb1*zt2*pow4(g3))/(9.*mgl*(-
     mmgl + mmsb1)*(mmsb1 - mmsb2)) - (8*mmgl*mmsb2*zt2*pow4(g3))/(3.*mgl*(-
     mmgl + mmsb1)*(mmsb1 - mmsb2)) + (8*mmgl*mmsb1*zt2*pow4(g3))/(3.*mgl*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (80*mmgl*mmsb2*zt2*pow4(g3))/(9.*mgl*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (8*mmgl*mmst1*zt2*pow4(g3))/(3.*mgl*(-
     mmgl + mmsb1)*(mmsb1 - mmsb2)) - (8*mmgl*mmst1*zt2*pow4(g3))/(3.*mgl*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (8*mmgl*mmst2*zt2*pow4(g3))/(3.*mgl*(-
     mmgl + mmsb1)*(mmsb1 - mmsb2)) - (8*mmgl*mmst2*zt2*pow4(g3))/(3.*mgl*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (64*mmgl*mmsusy*zt2*pow4(g3))/(3.*mgl*(-
     mmgl + mmsb1)*(mmsb1 - mmsb2)) + (64*mmgl*mmsusy*zt2*pow4(g3))/(3.*mgl*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*mmt*zt2*pow4(g3))/(3.*mgl*(-
     mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*mmt*zt2*pow4(g3))/(3.*mgl*(mmsb1
      - mmsb2)*(-mmgl + mmsb2)) - (112*mmgl*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl
     )*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) + (32*lgu*mmgl*mmsb1*mmst1*DeltaInv(
     mmt,mmst1,mmgl)*pow4(g3))/(mgl*(mmsb1 - mmsb2)) - (16*lt1u*mmgl*mmsb1*
     mmst1*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(mgl*(mmsb1 - mmsb2)) + (16*ltu*
     mmgl*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(mgl*(mmsb1 - mmsb2))
     - (32*lgu*ltu*mmgl*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(3.*mgl*
     (mmsb1 - mmsb2)) + (16*lt1u*ltu*mmgl*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*
     pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) + (112*mmgl*mmsb2*mmst1*DeltaInv(mmt,
     mmst1,mmgl)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) - (32*lgu*mmgl*mmsb2*mmst1*
     DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(mgl*(mmsb1 - mmsb2)) + (16*lt1u*mmgl*
     mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(mgl*(mmsb1 - mmsb2)) - (16
     *ltu*mmgl*mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(mgl*(mmsb1 -
     mmsb2)) + (32*lgu*ltu*mmgl*mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/
     (3.*mgl*(mmsb1 - mmsb2)) - (16*lt1u*ltu*mmgl*mmsb2*mmst1*DeltaInv(mmt,
     mmst1,mmgl)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) - (112*mmgl*mmsb1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) + (16*lgu*mmgl
     *mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(mgl*(mmsb1 - mmsb2)) + (16*
     ltu*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(mgl*(mmsb1 - mmsb2)
     ) - (16*lgu*ltu*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(3.*mgl*
     (mmsb1 - mmsb2)) + (112*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/
     (3.*mgl*(mmsb1 - mmsb2)) - (16*lgu*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)
     *pow4(g3))/(mgl*(mmsb1 - mmsb2)) - (16*ltu*mmgl*mmsb2*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow4(g3))/(mgl*(mmsb1 - mmsb2)) + (16*lgu*ltu*mmgl*mmsb2*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) + (448*mmgl*
     mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)
     *(mmsb1 - mmsb2)) - (16*lgu*mmgl*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*lt1u*mmgl*mmsb1*
     mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) - (32*lgu*lt1u*mmgl*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(
     g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (96*ltu*mmgl*mmsb1*mmst1*
     mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2
     )) + (16*lgu*ltu*mmgl*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(
     mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*lt1u*ltu*mmgl*mmsb1*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) -
     (448*mmgl*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(3.*mgl*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*lgu*mmgl*mmsb2*mmst1*mmt*DeltaInv(
     mmt,mmst1,mmgl)*pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*lt1u
     *mmgl*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) + (32*lgu*lt1u*mmgl*mmsb2*mmst1*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (96*ltu*
     mmgl*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) - (16*lgu*ltu*mmgl*mmsb2*mmst1*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*lt1u*ltu
     *mmgl*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*mmsb1*mmst1*zt2*DeltaInv(mmt,mmst1,mmgl
     )*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) + (16*mmgl*mmsb2*mmst1*zt2*DeltaInv(
     mmt,mmst1,mmgl)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) - (16*mmgl*mmsb1*mmt*
     zt2*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) + (16*mmgl
     *mmsb2*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2))
     + (64*mmgl*mmsb1*mmst1*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(3.*mgl*
     (-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (64*mmgl*mmsb2*mmst1*mmt*zt2*DeltaInv(
     mmt,mmst1,mmgl)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (112*
     mmgl*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2
     )) + (32*lgu*mmgl*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(mgl*(
     mmsb1 - mmsb2)) - (16*lt2u*mmgl*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow4(
     g3))/(mgl*(mmsb1 - mmsb2)) + (16*ltu*mmgl*mmsb1*mmst2*DeltaInv(mmt,mmst2,
     mmgl)*pow4(g3))/(mgl*(mmsb1 - mmsb2)) - (32*lgu*ltu*mmgl*mmsb1*mmst2*
     DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) + (16*lt2u*ltu
     *mmgl*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*mgl*(mmsb1 -
     mmsb2)) + (112*mmgl*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*mgl
     *(mmsb1 - mmsb2)) - (32*lgu*mmgl*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow4
     (g3))/(mgl*(mmsb1 - mmsb2)) + (16*lt2u*mmgl*mmsb2*mmst2*DeltaInv(mmt,mmst2
     ,mmgl)*pow4(g3))/(mgl*(mmsb1 - mmsb2)) - (16*ltu*mmgl*mmsb2*mmst2*DeltaInv
     (mmt,mmst2,mmgl)*pow4(g3))/(mgl*(mmsb1 - mmsb2)) + (32*lgu*ltu*mmgl*mmsb2*
     mmst2*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) - (16*
     lt2u*ltu*mmgl*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*mgl*(
     mmsb1 - mmsb2)) - (112*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(
     3.*mgl*(mmsb1 - mmsb2)) + (16*lgu*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*
     pow4(g3))/(mgl*(mmsb1 - mmsb2)) + (16*ltu*mmgl*mmsb1*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow4(g3))/(mgl*(mmsb1 - mmsb2)) - (16*lgu*ltu*mmgl*mmsb1*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) + (112*mmgl*
     mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) - (
     16*lgu*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(mgl*(mmsb1 -
     mmsb2)) - (16*ltu*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(mgl*(
     mmsb1 - mmsb2)) + (16*lgu*ltu*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4
     (g3))/(3.*mgl*(mmsb1 - mmsb2)) + (448*mmgl*mmsb1*mmst2*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*lgu*
     mmgl*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) - (16*lt2u*mmgl*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2
     ,mmgl)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (32*lgu*lt2u*mmgl
     *mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1
     )*(mmsb1 - mmsb2)) - (96*ltu*mmgl*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)
     *pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*lgu*ltu*mmgl*mmsb1*
     mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) + (16*lt2u*ltu*mmgl*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(
     g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (448*mmgl*mmsb2*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)
     ) + (16*lgu*mmgl*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(mgl*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*lt2u*mmgl*mmsb2*mmst2*mmt*DeltaInv(
     mmt,mmst2,mmgl)*pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (32*lgu*
     lt2u*mmgl*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*mgl*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (96*ltu*mmgl*mmsb2*mmst2*mmt*DeltaInv(
     mmt,mmst2,mmgl)*pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*lgu*
     ltu*mmgl*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) - (16*lt2u*ltu*mmgl*mmsb2*mmst2*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*
     mmsb1*mmst2*zt2*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)
     ) + (16*mmgl*mmsb2*mmst2*zt2*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*mgl*(
     mmsb1 - mmsb2)) - (16*mmgl*mmsb1*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow4(g3)
     )/(3.*mgl*(mmsb1 - mmsb2)) + (16*mmgl*mmsb2*mmt*zt2*DeltaInv(mmt,mmst2,
     mmgl)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) + (64*mmgl*mmsb1*mmst2*mmt*zt2*
     DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)
     ) - (64*mmgl*mmsb2*mmst2*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*
     mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (64*mmgl*Fin20(mmgl,mmsusy,mmu)*
     pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (64*mmgl*Fin20(mmgl,
     mmsusy,mmu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (176*mmgl
     *Fin20(mmsb1,mmgl,mmu)*pow4(g3))/(9.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2))
     + (32*mmgl*Fin20(mmsb1,mmgl,mmu)*pow4(g3))/(9.*mgl*(mmsb1 - mmsb2)*(-mmgl
     + mmsb2)) + (8*mmgl*Fin20(mmsb1,mmsb2,mmu)*pow4(g3))/(9.*mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) - (8*mmgl*Fin20(mmsb1,mmsb2,mmu)*pow4(g3))/(9.*mgl
     *(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (32*mmgl*Fin20(mmsb2,mmgl,mmu)*pow4(g3
     ))/(9.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (176*mmgl*Fin20(mmsb2,mmgl,
     mmu)*pow4(g3))/(9.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (8*mmgl*Fin3(mmt
     ,mmst1,mmgl,mmu)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (8*
     mmgl*Fin3(mmt,mmst1,mmgl,mmu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) + (16*mmgl*mmsb1*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu)
     *pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) - (16*mmgl*mmsb2*DeltaInv(mmt,mmst1,
     mmgl)*Fin3(mmt,mmst1,mmgl,mmu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) + (32*
     mmgl*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu)*pow4(g3
     ))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (32*mmgl*mmsb2*mmst1*
     DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu)*pow4(g3))/(3.*mgl*(mmsb1
      - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*
     Fin3(mmt,mmst1,mmgl,mmu)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)
     ) - (16*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu)*
     pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu)*pow4(g3))/(3.*mgl*(-mmgl
      + mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*
     Fin3(mmt,mmst1,mmgl,mmu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)
     ) + (8*mmgl*Fin3(mmt,mmst2,mmgl,mmu)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) - (8*mmgl*Fin3(mmt,mmst2,mmgl,mmu)*pow4(g3))/(3.*mgl*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*mmsb1*DeltaInv(mmt,mmst2,mmgl)*
     Fin3(mmt,mmst2,mmgl,mmu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) - (16*mmgl*
     mmsb2*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*pow4(g3))/(3.*mgl*
     (mmsb1 - mmsb2)) + (32*mmgl*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,
     mmst2,mmgl,mmu)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (32*
     mmgl*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*pow4(g3
     ))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*mmsb1*mmt*DeltaInv(
     mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)
     *(mmsb1 - mmsb2)) - (16*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,
     mmst2,mmgl,mmu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*
     mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*pow4(g3))
     /(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*mmst2*mmt*DeltaInv(
     mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)
     *(-mmgl + mmsb2)) + (44*mmgl*mmsb1*pow2(lb1u)*pow4(g3))/(9.*mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) + (64*mmgl*mmsb2*pow2(lb1u)*pow4(g3))/(9.*mgl*(-
     mmgl + mmsb1)*(mmsb1 - mmsb2)) + (4*mmgl*mmsb1*pow2(lb1u)*pow4(g3))/(3.*
     mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (4*mmgl*mmsb2*pow2(lb2u)*pow4(g3))/
     (3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (20*mmgl*mmsb2*pow2(lb2u)*pow4(
     g3))/(9.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (148*mmgl*mmsb1*pow2(lgu)*
     pow4(g3))/(9.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (4*mmgl*mmsb2*pow2(
     lgu)*pow4(g3))/(9.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (4*mmgl*mmsb1*
     pow2(lgu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (44*mmgl*
     mmsb2*pow2(lgu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (4*
     mmgl*mmst1*pow2(lgu)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) +
     (4*mmgl*mmst1*pow2(lgu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2))
     - (4*mmgl*mmst2*pow2(lgu)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2
     )) + (4*mmgl*mmst2*pow2(lgu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) - (32*mmgl*mmsusy*pow2(lgu)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) + (32*mmgl*mmsusy*pow2(lgu)*pow4(g3))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) + (8*mmgl*mmt*pow2(lgu)*pow4(g3))/(3.*mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) - (8*mmgl*mmt*pow2(lgu)*pow4(g3))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*
     pow2(lgu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) + (16*mmgl*mmsb2*mmst1*
     DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) - (8
     *mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow4(g3))/(3.*mgl*(
     mmsb1 - mmsb2)) + (8*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*
     pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) + (8*mmgl*mmsb1*mmst1*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow2(lgu)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) -
     (8*mmgl*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow4(g3))/(3.*
     mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*mmsb1*mmst2*DeltaInv(mmt,
     mmst2,mmgl)*pow2(lgu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) + (16*mmgl*mmsb2*
     mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)
     ) - (8*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow4(g3))/(3.*mgl
     *(mmsb1 - mmsb2)) + (8*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*
     pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) + (8*mmgl*mmsb1*mmst2*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow2(lgu)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) -
     (8*mmgl*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow4(g3))/(3.*
     mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (32*mmgl*mmsusy*pow2(lsu)*pow4(g3))
     /(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (32*mmgl*mmsusy*pow2(lsu)*pow4
     (g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (4*mmgl*mmst1*pow2(lt1u)*
     pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (4*mmgl*mmst1*pow2(
     lt1u)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (8*mmgl*mmsb1*
     mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(lt1u)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2
     )) - (8*mmgl*mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(lt1u)*pow4(g3))/(3.
     *mgl*(mmsb1 - mmsb2)) + (8*mmgl*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow2(lt1u)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (8*mmgl*
     mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lt1u)*pow4(g3))/(3.*mgl*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (4*mmgl*mmst2*pow2(lt2u)*pow4(g3))/(3.*
     mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (4*mmgl*mmst2*pow2(lt2u)*pow4(g3))/
     (3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (8*mmgl*mmsb1*mmst2*DeltaInv(
     mmt,mmst2,mmgl)*pow2(lt2u)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) - (8*mmgl*
     mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(lt2u)*pow4(g3))/(3.*mgl*(mmsb1 -
     mmsb2)) + (8*mmgl*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lt2u)*pow4
     (g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (8*mmgl*mmsb2*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(lt2u)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-
     mmgl + mmsb2)) + (16*mmgl*mmsb1*pow2(ltu)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1
     )*(mmsb1 - mmsb2)) - (16*mmgl*mmsb2*pow2(ltu)*pow4(g3))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) + (8*mmgl*mmst1*pow2(ltu)*pow4(g3))/(3.*mgl*(-mmgl
      + mmsb1)*(mmsb1 - mmsb2)) - (8*mmgl*mmst1*pow2(ltu)*pow4(g3))/(3.*mgl*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (8*mmgl*mmst2*pow2(ltu)*pow4(g3))/(3.*
     mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (8*mmgl*mmst2*pow2(ltu)*pow4(g3))/(
     3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (8*mmgl*mmt*pow2(ltu)*pow4(g3))/
     (3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (8*mmgl*mmt*pow2(ltu)*pow4(g3))
     /(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (8*mmgl*mmsb1*mmst1*DeltaInv(
     mmt,mmst1,mmgl)*pow2(ltu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) + (8*mmgl*
     mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow4(g3))/(3.*mgl*(mmsb1 -
     mmsb2)) - (8*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow4(g3))/(
     3.*mgl*(mmsb1 - mmsb2)) + (8*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     ltu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) + (16*mmgl*mmsb1*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) - (16*mmgl*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow4
     (g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (8*mmgl*mmsb1*mmst2*DeltaInv
     (mmt,mmst2,mmgl)*pow2(ltu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) + (8*mmgl*
     mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow4(g3))/(3.*mgl*(mmsb1 -
     mmsb2)) - (8*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow4(g3))/(
     3.*mgl*(mmsb1 - mmsb2)) + (8*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(
     ltu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) + (16*mmgl*mmsb1*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) - (16*mmgl*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow4
     (g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (112*mmsb1*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmgl)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) - (16*lgu*mmsb1*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmgl)*pow4(g3))/(mgl*(mmsb1 - mmsb2)) - (16*
     ltu*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmgl)*pow4(g3))/(mgl*(mmsb1 -
     mmsb2)) + (16*lgu*ltu*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmgl)*pow4(g3))/
     (3.*mgl*(mmsb1 - mmsb2)) - (112*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmgl)*
     pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) + (16*lgu*mmsb2*DeltaInv(mmt,mmst1,mmgl
     )*pow2(mmgl)*pow4(g3))/(mgl*(mmsb1 - mmsb2)) + (16*ltu*mmsb2*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmgl)*pow4(g3))/(mgl*(mmsb1 - mmsb2)) - (16*lgu*ltu*mmsb2
     *DeltaInv(mmt,mmst1,mmgl)*pow2(mmgl)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) +
     (16*mmsb1*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmgl)*pow4(g3))/(3.*mgl*(mmsb1
      - mmsb2)) - (16*mmsb2*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmgl)*pow4(g3))/(
     3.*mgl*(mmsb1 - mmsb2)) + (112*mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(mmgl)*
     pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) - (16*lgu*mmsb1*DeltaInv(mmt,mmst2,mmgl
     )*pow2(mmgl)*pow4(g3))/(mgl*(mmsb1 - mmsb2)) - (16*ltu*mmsb1*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmgl)*pow4(g3))/(mgl*(mmsb1 - mmsb2)) + (16*lgu*ltu*mmsb1
     *DeltaInv(mmt,mmst2,mmgl)*pow2(mmgl)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) -
     (112*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmgl)*pow4(g3))/(3.*mgl*(mmsb1 -
     mmsb2)) + (16*lgu*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmgl)*pow4(g3))/(mgl
     *(mmsb1 - mmsb2)) + (16*ltu*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmgl)*pow4
     (g3))/(mgl*(mmsb1 - mmsb2)) - (16*lgu*ltu*mmsb2*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmgl)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) + (16*mmsb1*zt2*DeltaInv(mmt
     ,mmst2,mmgl)*pow2(mmgl)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) - (16*mmsb2*zt2
     *DeltaInv(mmt,mmst2,mmgl)*pow2(mmgl)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) +
     (8*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow2(mmgl)*pow4(g3))/(3.*mgl*(
     mmsb1 - mmsb2)) - (8*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow2(mmgl)*
     pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) + (8*mmsb1*DeltaInv(mmt,mmst2,mmgl)*
     pow2(lgu)*pow2(mmgl)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) - (8*mmsb2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(mmgl)*pow4(g3))/(3.*mgl*(mmsb1 -
     mmsb2)) + (8*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow2(mmgl)*pow4(g3))
     /(3.*mgl*(mmsb1 - mmsb2)) - (8*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*
     pow2(mmgl)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) + (8*mmsb1*DeltaInv(mmt,
     mmst2,mmgl)*pow2(ltu)*pow2(mmgl)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) - (8*
     mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow2(mmgl)*pow4(g3))/(3.*mgl*(
     mmsb1 - mmsb2)) + (112*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow4(g3))
     /(3.*mgl*(mmsb1 - mmsb2)) - (16*lgu*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmsb1)*pow4(g3))/(mgl*(mmsb1 - mmsb2)) - (16*ltu*mmgl*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmsb1)*pow4(g3))/(mgl*(mmsb1 - mmsb2)) + (16*lgu*ltu*mmgl*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) +
     (112*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(-
     mmgl + mmsb1)*(mmsb1 - mmsb2)) - (32*lgu*mmgl*mmst1*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmsb1)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*
     lt1u*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow4(g3))/(mgl*(-mmgl
      + mmsb1)*(mmsb1 - mmsb2)) - (16*ltu*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmsb1)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (32*lgu*ltu*
     mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) - (16*lt1u*ltu*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)
     *pow2(mmsb1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (112*
     mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) - (16*lgu*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmsb1)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*ltu*mmgl*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1
      - mmsb2)) + (16*lgu*ltu*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*
     pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*mmgl*zt2*DeltaInv
     (mmt,mmst1,mmgl)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) + (16*mmgl
     *mmst1*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) + (16*mmgl*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmsb1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (112*mmgl*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) -
     (16*lgu*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3))/(mgl*(mmsb1 -
     mmsb2)) - (16*ltu*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3))/(mgl
     *(mmsb1 - mmsb2)) + (16*lgu*ltu*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*
     pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) + (112*mmgl*mmst2*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (32
     *lgu*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3))/(mgl*(-mmgl
      + mmsb1)*(mmsb1 - mmsb2)) + (16*lt2u*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmsb1)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*ltu*mmgl
     *mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3))/(mgl*(-mmgl + mmsb1)
     *(mmsb1 - mmsb2)) + (32*lgu*ltu*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmsb1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*lt2u*ltu*
     mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) + (112*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmsb1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*lgu*mmgl*
     mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) - (16*ltu*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*
     pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*lgu*ltu*mmgl*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) + (16*mmgl*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(
     g3))/(3.*mgl*(mmsb1 - mmsb2)) + (16*mmgl*mmst2*zt2*DeltaInv(mmt,mmst2,mmgl
     )*pow2(mmsb1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*
     mmgl*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(-mmgl
      + mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,
     mmst1,mmgl,mmu)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) - (16*mmgl*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*pow2(
     mmsb1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (8*mmgl*
     DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(mmsb1 -
     mmsb2)) + (16*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow2(mmsb1)*
     pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (8*mmgl*mmt*DeltaInv(
     mmt,mmst1,mmgl)*pow2(lgu)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) + (8*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(mmsb1)*
     pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) + (16*mmgl*mmst2*DeltaInv(mmt,mmst2,
     mmgl)*pow2(lgu)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) + (8*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(mmsb1)*pow4(
     g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (8*mmgl*mmst1*DeltaInv(mmt
     ,mmst1,mmgl)*pow2(lt1u)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) - (8*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(lt2u)*pow2(
     mmsb1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (8*mmgl*
     DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(mmsb1 -
     mmsb2)) + (8*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow2(mmsb1)*
     pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (8*mmgl*mmt*DeltaInv(
     mmt,mmst1,mmgl)*pow2(ltu)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) + (8*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow2(mmsb1)*
     pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) + (8*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl
     )*pow2(ltu)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2))
     + (8*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow2(mmsb1)*pow4(g3))/(3.
     *mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*lb1u*mmgl*mmsb1*mmsb2*pow4(g3)
     )/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (8*lb1u*lb2u*mmgl*mmsb1*
     mmsb2*pow4(g3))/(9.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (16*lgu*
     mmgl*mmsb1*mmsb2*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) -
     (16*lb1u*lgu*mmgl*mmsb1*mmsb2*pow4(g3))/(9.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl
      + mmsb1)) + (8*lb2u*lgu*mmgl*mmsb1*mmsb2*pow4(g3))/(9.*mgl*(mmsb1 - mmsb2
     )*pow2(-mmgl + mmsb1)) + (16*lb1u*mmgl*mmsb1*mmst1*pow4(g3))/(3.*mgl*(
     mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (16*lgu*mmgl*mmsb1*mmst1*pow4(g3))/(
     3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (8*lb1u*lt1u*mmgl*mmsb1*
     mmst1*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) + (8*lgu*lt1u
     *mmgl*mmsb1*mmst1*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) +
     (16*lb1u*mmgl*mmsb1*mmst2*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl +
     mmsb1)) - (16*lgu*mmgl*mmsb1*mmst2*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(
     -mmgl + mmsb1)) - (8*lb1u*lt2u*mmgl*mmsb1*mmst2*pow4(g3))/(3.*mgl*(mmsb1 -
     mmsb2)*pow2(-mmgl + mmsb1)) + (8*lgu*lt2u*mmgl*mmsb1*mmst2*pow4(g3))/(3.*
     mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) + (128*lb1u*mmgl*mmsb1*mmsusy*
     pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (128*lgu*mmgl*
     mmsb1*mmsusy*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (64*
     lb1u*lsu*mmgl*mmsb1*mmsusy*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl +
     mmsb1)) + (64*lgu*lsu*mmgl*mmsb1*mmsusy*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*
     pow2(-mmgl + mmsb1)) - (32*lb1u*mmgl*mmsb1*mmt*pow4(g3))/(3.*mgl*(mmsb1 -
     mmsb2)*pow2(-mmgl + mmsb1)) + (32*lgu*mmgl*mmsb1*mmt*pow4(g3))/(3.*mgl*(
     mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) + (16*lb1u*ltu*mmgl*mmsb1*mmt*pow4(g3)
     )/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (16*lgu*ltu*mmgl*mmsb1*
     mmt*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (64*mmgl*
     mmsb1*Fin20(mmgl,mmsusy,mmu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl
     + mmsb1)) + (64*mmgl*mmsusy*Fin20(mmgl,mmsusy,mmu)*pow4(g3))/(3.*mgl*(
     mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmsb1*Fin20(mmsb1,mmsb2,mmu)
     *pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmsb2*
     Fin20(mmsb1,mmsb2,mmu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl +
     mmsb1)) + (64*mmgl*mmsb1*Fin20(mmsb1,mmsusy,mmu)*pow4(g3))/(3.*mgl*(mmsb1
     - mmsb2)*pow2(-mmgl + mmsb1)) - (64*mmgl*mmsusy*Fin20(mmsb1,mmsusy,mmu)*
     pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmsb1*
     Fin20(mmsb2,mmgl,mmu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1
     )) + (8*mmgl*mmsb2*Fin20(mmsb2,mmgl,mmu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)
     *pow2(-mmgl + mmsb1)) + (8*mmgl*mmsb1*Fin3(mmt,mmsb1,mmst1,mmu)*pow4(g3))/
     (3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmst1*Fin3(mmt,
     mmsb1,mmst1,mmu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) +
     (8*mmgl*mmt*Fin3(mmt,mmsb1,mmst1,mmu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*
     pow2(-mmgl + mmsb1)) + (8*mmgl*mmsb1*Fin3(mmt,mmsb1,mmst2,mmu)*pow4(g3))/(
     3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmst2*Fin3(mmt,mmsb1
     ,mmst2,mmu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) + (8*
     mmgl*mmt*Fin3(mmt,mmsb1,mmst2,mmu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(
     -mmgl + mmsb1)) - (8*mmgl*mmsb1*Fin3(mmt,mmst1,mmgl,mmu)*pow4(g3))/(3.*mgl
     *(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmst1*Fin3(mmt,mmst1,mmgl,
     mmu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmt*
     Fin3(mmt,mmst1,mmgl,mmu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl +
     mmsb1)) - (8*mmgl*mmsb1*Fin3(mmt,mmst2,mmgl,mmu)*pow4(g3))/(3.*mgl*(mmsb1
     - mmsb2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmst2*Fin3(mmt,mmst2,mmgl,mmu)*
     pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmt*Fin3(
     mmt,mmst2,mmgl,mmu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1))
     - (4*mmgl*mmsb1*mmsb2*pow2(lb1u)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-
     mmgl + mmsb1)) - (4*mmgl*mmsb1*mmst1*pow2(lb1u)*pow4(g3))/(3.*mgl*(mmsb1 -
     mmsb2)*pow2(-mmgl + mmsb1)) - (4*mmgl*mmsb1*mmst2*pow2(lb1u)*pow4(g3))/(3.
     *mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (32*mmgl*mmsb1*mmsusy*pow2(
     lb1u)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) + (8*mmgl*
     mmsb1*mmt*pow2(lb1u)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)
     ) + (28*mmgl*mmsb1*mmsb2*pow2(lgu)*pow4(g3))/(9.*mgl*(mmsb1 - mmsb2)*pow2(
     -mmgl + mmsb1)) + (4*mmgl*mmsb1*mmst1*pow2(lgu)*pow4(g3))/(3.*mgl*(mmsb1 -
     mmsb2)*pow2(-mmgl + mmsb1)) + (4*mmgl*mmsb1*mmst2*pow2(lgu)*pow4(g3))/(3.*
     mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) + (32*mmgl*mmsb1*mmsusy*pow2(lgu)
     *pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmsb1*
     mmt*pow2(lgu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (80
     *lb1u*mmgl*pow2(mmsb1)*pow4(g3))/(mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1))
     + (80*lgu*mmgl*pow2(mmsb1)*pow4(g3))/(mgl*(mmsb1 - mmsb2)*pow2(-mmgl +
     mmsb1)) + (296*lb1u*lgu*mmgl*pow2(mmsb1)*pow4(g3))/(9.*mgl*(mmsb1 - mmsb2)
     *pow2(-mmgl + mmsb1)) + (124*mmgl*pow2(lb1u)*pow2(mmsb1)*pow4(g3))/(9.*mgl
     *(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (140*mmgl*pow2(lgu)*pow2(mmsb1)*
     pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (112*mmgl*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) +
     (16*lgu*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(g3))/(mgl*(mmsb1 -
     mmsb2)) + (16*ltu*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(g3))/(mgl
     *(mmsb1 - mmsb2)) - (16*lgu*ltu*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*
     pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) - (112*mmgl*mmst1*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (32
     *lgu*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(g3))/(mgl*(mmsb1
      - mmsb2)*(-mmgl + mmsb2)) - (16*lt1u*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmsb2)*pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*ltu*mmgl
     *mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(g3))/(mgl*(mmsb1 - mmsb2)
     *(-mmgl + mmsb2)) - (32*lgu*ltu*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*lt1u*ltu*
     mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) - (112*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*lgu*mmgl*
     mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-
     mmgl + mmsb2)) + (16*ltu*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*
     pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*lgu*ltu*mmgl*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-
     mmgl + mmsb2)) - (16*mmgl*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(g3
     ))/(3.*mgl*(mmsb1 - mmsb2)) - (16*mmgl*mmst1*zt2*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*
     mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) - (112*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*
     pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) + (16*lgu*mmgl*DeltaInv(mmt,mmst2,mmgl)
     *pow2(mmsb2)*pow4(g3))/(mgl*(mmsb1 - mmsb2)) + (16*ltu*mmgl*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmsb2)*pow4(g3))/(mgl*(mmsb1 - mmsb2)) - (16*lgu*ltu*mmgl
     *DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) -
     (112*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (32*lgu*mmgl*mmst2*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmsb2)*pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*
     lt2u*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow4(g3))/(mgl*(mmsb1
      - mmsb2)*(-mmgl + mmsb2)) + (16*ltu*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmsb2)*pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (32*lgu*ltu*
     mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) + (16*lt2u*ltu*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)
     *pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (112*
     mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) + (16*lgu*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmsb2)*pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*ltu*mmgl*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl
      + mmsb2)) - (16*lgu*ltu*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*
     pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*zt2*DeltaInv
     (mmt,mmst2,mmgl)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) - (16*mmgl
     *mmst2*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*
     DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu)*pow2(mmsb2)*pow4(g3))/(
     3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*DeltaInv(mmt,mmst2,mmgl
     )*Fin3(mmt,mmst2,mmgl,mmu)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(
     -mmgl + mmsb2)) - (8*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow2(mmsb2)*
     pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) - (16*mmgl*mmst1*DeltaInv(mmt,mmst1,
     mmgl)*pow2(lgu)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) - (8*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow2(mmsb2)*pow4(
     g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (8*mmgl*DeltaInv(mmt,mmst2
     ,mmgl)*pow2(lgu)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) - (16*mmgl
     *mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (8*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*
     pow2(lgu)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) +
     (8*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(lt1u)*pow2(mmsb2)*pow4(g3))/(
     3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (8*mmgl*mmst2*DeltaInv(mmt,mmst2
     ,mmgl)*pow2(lt2u)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) - (8*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow2(mmsb2)*pow4(g3))
     /(3.*mgl*(mmsb1 - mmsb2)) - (8*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(
     ltu)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (8*
     mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow2(mmsb2)*pow4(g3))/(3.*mgl*
     (mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (8*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(
     ltu)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)) - (8*mmgl*mmst2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) - (8*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*
     pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*lb1u*
     lb2u*mmgl*pow2(mmsb2)*pow4(g3))/(9.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1
     )) - (16*lb1u*lgu*mmgl*pow2(mmsb2)*pow4(g3))/(9.*mgl*(mmsb1 - mmsb2)*pow2(
     -mmgl + mmsb1)) - (16*lb2u*lgu*mmgl*pow2(mmsb2)*pow4(g3))/(9.*mgl*(mmsb1 -
     mmsb2)*pow2(-mmgl + mmsb1)) + (16*mmgl*pow2(lgu)*pow2(mmsb2)*pow4(g3))/(9.
     *mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) + (16*lb1u*mmgl*pow2(mmsb2)*pow4
     (g3))/(9.*mgl*(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (16*lb2u*mmgl*pow2(
     mmsb2)*pow4(g3))/(9.*mgl*(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (16*lb1u*
     lb2u*mmgl*pow2(mmsb2)*pow4(g3))/(9.*mgl*(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2
     )) - (112*lb1u*lgu*mmgl*pow2(mmsb2)*pow4(g3))/(9.*mgl*(-mmgl + mmsb1)*pow2
     (mmsb1 - mmsb2)) + (16*lb2u*lgu*mmgl*pow2(mmsb2)*pow4(g3))/(9.*mgl*(-mmgl
     + mmsb1)*pow2(mmsb1 - mmsb2)) - (16*lb1u*mmgl*pow2(mmsb2)*pow4(g3))/(9.*
     mgl*(-mmgl + mmsb2)*pow2(mmsb1 - mmsb2)) + (16*lb2u*mmgl*pow2(mmsb2)*pow4(
     g3))/(9.*mgl*(-mmgl + mmsb2)*pow2(mmsb1 - mmsb2)) - (112*lb1u*lb2u*mmgl*
     pow2(mmsb2)*pow4(g3))/(9.*mgl*(-mmgl + mmsb2)*pow2(mmsb1 - mmsb2)) + (112*
     lb1u*lgu*mmgl*pow2(mmsb2)*pow4(g3))/(9.*mgl*(-mmgl + mmsb2)*pow2(mmsb1 -
     mmsb2)) - (16*lb2u*lgu*mmgl*pow2(mmsb2)*pow4(g3))/(9.*mgl*(-mmgl + mmsb2)*
     pow2(mmsb1 - mmsb2)) + (64*mmgl*pow2(lb1u)*pow2(mmsb2)*pow4(g3))/(9.*mgl*(
     -mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) + (64*mmgl*pow2(lb2u)*pow2(mmsb2)*pow4
     (g3))/(9.*mgl*(-mmgl + mmsb2)*pow2(mmsb1 - mmsb2)) + (16*mmgl*pow2(lgu)*
     pow2(mmsb2)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (16*
     mmgl*pow2(lgu)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*pow2(mmsb1 -
     mmsb2)) - (16*lb2u*mmgl*mmsb1*mmsb2*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2
     (-mmgl + mmsb2)) + (8*lb1u*lb2u*mmgl*mmsb1*mmsb2*pow4(g3))/(3.*mgl*(mmsb1
     - mmsb2)*pow2(-mmgl + mmsb2)) + (16*lgu*mmgl*mmsb1*mmsb2*pow4(g3))/(3.*mgl
     *(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (8*lb1u*lgu*mmgl*mmsb1*mmsb2*pow4(
     g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (16*lb2u*mmgl*mmsb2*
     mmst1*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (16*lgu*
     mmgl*mmsb2*mmst1*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) +
     (8*lb2u*lt1u*mmgl*mmsb2*mmst1*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl
      + mmsb2)) - (8*lgu*lt1u*mmgl*mmsb2*mmst1*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2
     )*pow2(-mmgl + mmsb2)) - (16*lb2u*mmgl*mmsb2*mmst2*pow4(g3))/(3.*mgl*(
     mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (16*lgu*mmgl*mmsb2*mmst2*pow4(g3))/(
     3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (8*lb2u*lt2u*mmgl*mmsb2*
     mmst2*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (8*lgu*lt2u
     *mmgl*mmsb2*mmst2*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) -
     (128*lb2u*mmgl*mmsb2*mmsusy*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl +
     mmsb2)) + (128*lgu*mmgl*mmsb2*mmsusy*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*
     pow2(-mmgl + mmsb2)) + (64*lb2u*lsu*mmgl*mmsb2*mmsusy*pow4(g3))/(3.*mgl*(
     mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (64*lgu*lsu*mmgl*mmsb2*mmsusy*pow4(
     g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (32*lb2u*mmgl*mmsb2*
     mmt*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (32*lgu*mmgl*
     mmsb2*mmt*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (16*
     lb2u*ltu*mmgl*mmsb2*mmt*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl +
     mmsb2)) + (16*lgu*ltu*mmgl*mmsb2*mmt*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*
     pow2(-mmgl + mmsb2)) + (64*mmgl*mmsb2*Fin20(mmgl,mmsusy,mmu)*pow4(g3))/(3.
     *mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (64*mmgl*mmsusy*Fin20(mmgl,
     mmsusy,mmu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (8*
     mmgl*mmsb1*Fin20(mmsb1,mmgl,mmu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-
     mmgl + mmsb2)) + (8*mmgl*mmsb2*Fin20(mmsb1,mmgl,mmu)*pow4(g3))/(3.*mgl*(
     mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmsb1*Fin20(mmsb1,mmsb2,mmu)
     *pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmsb2*
     Fin20(mmsb1,mmsb2,mmu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl +
     mmsb2)) - (64*mmgl*mmsb2*Fin20(mmsb2,mmsusy,mmu)*pow4(g3))/(3.*mgl*(mmsb1
     - mmsb2)*pow2(-mmgl + mmsb2)) + (64*mmgl*mmsusy*Fin20(mmsb2,mmsusy,mmu)*
     pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmsb2*
     Fin3(mmt,mmsb2,mmst1,mmu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl +
     mmsb2)) + (8*mmgl*mmst1*Fin3(mmt,mmsb2,mmst1,mmu)*pow4(g3))/(3.*mgl*(mmsb1
      - mmsb2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmt*Fin3(mmt,mmsb2,mmst1,mmu)*
     pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmsb2*
     Fin3(mmt,mmsb2,mmst2,mmu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl +
     mmsb2)) + (8*mmgl*mmst2*Fin3(mmt,mmsb2,mmst2,mmu)*pow4(g3))/(3.*mgl*(mmsb1
      - mmsb2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmt*Fin3(mmt,mmsb2,mmst2,mmu)*
     pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmsb2*
     Fin3(mmt,mmst1,mmgl,mmu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl +
     mmsb2)) - (8*mmgl*mmst1*Fin3(mmt,mmst1,mmgl,mmu)*pow4(g3))/(3.*mgl*(mmsb1
     - mmsb2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmt*Fin3(mmt,mmst1,mmgl,mmu)*pow4(
     g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmsb2*Fin3(mmt
     ,mmst2,mmgl,mmu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) -
     (8*mmgl*mmst2*Fin3(mmt,mmst2,mmgl,mmu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*
     pow2(-mmgl + mmsb2)) + (8*mmgl*mmt*Fin3(mmt,mmst2,mmgl,mmu)*pow4(g3))/(3.*
     mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (4*mmgl*mmsb1*mmsb2*pow2(lb2u)*
     pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (4*mmgl*mmsb2*
     mmst1*pow2(lb2u)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) +
     (4*mmgl*mmsb2*mmst2*pow2(lb2u)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-
     mmgl + mmsb2)) + (32*mmgl*mmsb2*mmsusy*pow2(lb2u)*pow4(g3))/(3.*mgl*(mmsb1
      - mmsb2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmsb2*mmt*pow2(lb2u)*pow4(g3))/(
     3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (4*mmgl*mmsb1*mmsb2*pow2(lgu
     )*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (4*mmgl*mmsb2*
     mmst1*pow2(lgu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (
     4*mmgl*mmsb2*mmst2*pow2(lgu)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl
     + mmsb2)) - (32*mmgl*mmsb2*mmsusy*pow2(lgu)*pow4(g3))/(3.*mgl*(mmsb1 -
     mmsb2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmsb2*mmt*pow2(lgu)*pow4(g3))/(3.*
     mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (80*lb2u*mmgl*pow2(mmsb2)*pow4(
     g3))/(mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (16*lb1u*lb2u*mmgl*pow2(
     mmsb2)*pow4(g3))/(9.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (80*lgu*
     mmgl*pow2(mmsb2)*pow4(g3))/(mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (16
     *lb1u*lgu*mmgl*pow2(mmsb2)*pow4(g3))/(9.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl +
     mmsb2)) - (104*lb2u*lgu*mmgl*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)
     *pow2(-mmgl + mmsb2)) - (124*mmgl*pow2(lb2u)*pow2(mmsb2)*pow4(g3))/(9.*mgl
     *(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (436*mmgl*pow2(lgu)*pow2(mmsb2)*
     pow4(g3))/(9.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (112*mmgl*mmsb1*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) + (16*lgu*mmgl*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*
     pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (32*lt1u*mmgl*mmsb1*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1
      - mmsb2)) - (16*ltu*mmgl*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow4(
     g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*lgu*ltu*mmgl*mmsb1*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) + (32*lt1u*ltu*mmgl*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmst1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (112*mmgl*
     mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow4(g3))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) - (16*lgu*mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2
     (mmst1)*pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (32*lt1u*mmgl*
     mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow4(g3))/(mgl*(mmsb1 - mmsb2)*
     (-mmgl + mmsb2)) + (16*ltu*mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)
     *pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*lgu*ltu*mmgl*mmsb2*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-
     mmgl + mmsb2)) - (32*lt1u*ltu*mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmst1)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (112*mmgl*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) - (16*lt1u*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*
     pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*ltu*mmgl*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1
      - mmsb2)) + (16*lt1u*ltu*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*
     pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (112*mmgl*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-
     mmgl + mmsb2)) + (16*lt1u*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*
     pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*ltu*mmgl*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl
      + mmsb2)) - (16*lt1u*ltu*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*
     pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*mmsb1*zt2*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) - (16*mmgl*mmsb2*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*
     pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*mmt*zt2*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) - (16*mmgl*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*
     pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*DeltaInv(mmt
     ,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu)*pow2(mmst1)*pow4(g3))/(3.*mgl*(-mmgl
      + mmsb1)*(mmsb1 - mmsb2)) + (16*mmgl*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,
     mmst1,mmgl,mmu)*pow2(mmst1)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) - (8*mmgl*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow2(mmst1)*
     pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (8*mmgl*mmsb2*
     DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow2(mmst1)*pow4(g3))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(
     lt1u)*pow2(mmst1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16
     *mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(lt1u)*pow2(mmst1)*pow4(g3))/(3.*
     mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (8*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl
     )*pow2(lt1u)*pow2(mmst1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)
     ) - (8*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lt1u)*pow2(mmst1)*pow4(g3))/
     (3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (8*mmgl*mmsb1*DeltaInv(mmt,
     mmst1,mmgl)*pow2(ltu)*pow2(mmst1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1
      - mmsb2)) - (8*mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow2(mmst1)*
     pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (8*mmgl*mmt*DeltaInv(
     mmt,mmst1,mmgl)*pow2(ltu)*pow2(mmst1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) - (8*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow2(
     mmst1)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (112*mmgl*
     mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(3.*mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) + (16*lgu*mmgl*mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2
     (mmst2)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (32*lt2u*mmgl*
     mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(mgl*(-mmgl + mmsb1)*
     (mmsb1 - mmsb2)) - (16*ltu*mmgl*mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)
     *pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*lgu*ltu*mmgl*mmsb1*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) + (32*lt2u*ltu*mmgl*mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmst2)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (112*mmgl*
     mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) - (16*lgu*mmgl*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2
     (mmst2)*pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (32*lt2u*mmgl*
     mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(mgl*(mmsb1 - mmsb2)*
     (-mmgl + mmsb2)) + (16*ltu*mmgl*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)
     *pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*lgu*ltu*mmgl*mmsb2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-
     mmgl + mmsb2)) - (32*lt2u*ltu*mmgl*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmst2)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (112*mmgl*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) - (16*lt2u*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*
     pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*ltu*mmgl*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1
      - mmsb2)) + (16*lt2u*ltu*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*
     pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (112*mmgl*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-
     mmgl + mmsb2)) + (16*lt2u*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*
     pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*ltu*mmgl*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl
      + mmsb2)) - (16*lt2u*ltu*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*
     pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*mmsb1*zt2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) - (16*mmgl*mmsb2*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*
     pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*mmt*zt2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) - (16*mmgl*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*
     pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*DeltaInv(mmt
     ,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*pow2(mmst2)*pow4(g3))/(3.*mgl*(-mmgl
      + mmsb1)*(mmsb1 - mmsb2)) + (16*mmgl*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,
     mmst2,mmgl,mmu)*pow2(mmst2)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) - (8*mmgl*mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(mmst2)*
     pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (8*mmgl*mmsb2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(mmst2)*pow4(g3))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(
     lt2u)*pow2(mmst2)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16
     *mmgl*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(lt2u)*pow2(mmst2)*pow4(g3))/(3.*
     mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (8*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl
     )*pow2(lt2u)*pow2(mmst2)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)
     ) - (8*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lt2u)*pow2(mmst2)*pow4(g3))/
     (3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (8*mmgl*mmsb1*DeltaInv(mmt,
     mmst2,mmgl)*pow2(ltu)*pow2(mmst2)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1
      - mmsb2)) - (8*mmgl*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow2(mmst2)*
     pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (8*mmgl*mmt*DeltaInv(
     mmt,mmst2,mmgl)*pow2(ltu)*pow2(mmst2)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) - (8*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow2(
     mmst2)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (112*mmgl*
     DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) + (16*lgu*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1)*pow4(
     g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*ltu*mmgl*DeltaInv(mmt,
     mmst1,mmgl)*pow3(mmsb1)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) -
     (16*lgu*ltu*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1)*pow4(g3))/(3.*mgl*(-
     mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*zt2*DeltaInv(mmt,mmst1,mmgl)*
     pow3(mmsb1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (112*mmgl
     *DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) + (16*lgu*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1)*pow4(
     g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*ltu*mmgl*DeltaInv(mmt,
     mmst2,mmgl)*pow3(mmsb1)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) -
     (16*lgu*ltu*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1)*pow4(g3))/(3.*mgl*(-
     mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*zt2*DeltaInv(mmt,mmst2,mmgl)*
     pow3(mmsb1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (8*mmgl*
     DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow3(mmsb1)*pow4(g3))/(3.*mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) - (8*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow3(
     mmsb1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (8*mmgl*
     DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow3(mmsb1)*pow4(g3))/(3.*mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) - (8*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow3(
     mmsb1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (32*lb1u*lgu*
     mmgl*pow3(mmsb1)*pow4(g3))/(9.*mgl*(mmsb1 - mmsb2)*pow3(-mmgl + mmsb1)) +
     (16*mmgl*pow2(lb1u)*pow3(mmsb1)*pow4(g3))/(9.*mgl*(mmsb1 - mmsb2)*pow3(-
     mmgl + mmsb1)) + (16*mmgl*pow2(lgu)*pow3(mmsb1)*pow4(g3))/(9.*mgl*(mmsb1 -
     mmsb2)*pow3(-mmgl + mmsb1)) + (112*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow3(
     mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*lgu*mmgl*
     DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2)*pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl
      + mmsb2)) - (16*ltu*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2)*pow4(g3))/(
     mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*lgu*ltu*mmgl*DeltaInv(mmt,mmst1
     ,mmgl)*pow3(mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (
     16*mmgl*zt2*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2)*pow4(g3))/(3.*mgl*(mmsb1
     - mmsb2)*(-mmgl + mmsb2)) + (112*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2)
     *pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*lgu*mmgl*
     DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2)*pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl
      + mmsb2)) - (16*ltu*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2)*pow4(g3))/(
     mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*lgu*ltu*mmgl*DeltaInv(mmt,mmst2
     ,mmgl)*pow3(mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (
     16*mmgl*zt2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2)*pow4(g3))/(3.*mgl*(mmsb1
     - mmsb2)*(-mmgl + mmsb2)) + (8*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*
     pow3(mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (8*mmgl*
     DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow3(mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) + (8*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow3(
     mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (8*mmgl*
     DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow3(mmsb2)*pow4(g3))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) + (16*lb1u*lb2u*mmgl*pow3(mmsb2)*pow4(g3))/(9.*mgl
     *pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (16*lb1u*lgu*mmgl*pow3(mmsb2)*
     pow4(g3))/(9.*mgl*pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (16*lb2u*lgu*
     mmgl*pow3(mmsb2)*pow4(g3))/(9.*mgl*pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)
     ) + (16*mmgl*pow2(lgu)*pow3(mmsb2)*pow4(g3))/(9.*mgl*pow2(-mmgl + mmsb1)*
     pow2(mmsb1 - mmsb2)) + (16*lb1u*lb2u*mmgl*pow3(mmsb2)*pow4(g3))/(9.*mgl*
     pow2(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (16*lb1u*lgu*mmgl*pow3(mmsb2)*
     pow4(g3))/(9.*mgl*pow2(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (16*lb2u*lgu*
     mmgl*pow3(mmsb2)*pow4(g3))/(9.*mgl*pow2(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)
     ) + (16*mmgl*pow2(lgu)*pow3(mmsb2)*pow4(g3))/(9.*mgl*pow2(mmsb1 - mmsb2)*
     pow2(-mmgl + mmsb2)) + (32*lb1u*lb2u*mmgl*pow3(mmsb2)*pow4(g3))/(9.*mgl*(-
     mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) - (32*lb1u*lgu*mmgl*pow3(mmsb2)*pow4(g3
     ))/(9.*mgl*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) - (32*lb2u*lgu*mmgl*pow3(
     mmsb2)*pow4(g3))/(9.*mgl*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) - (32*lb1u*
     lb2u*mmgl*pow3(mmsb2)*pow4(g3))/(9.*mgl*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2
     )) + (32*lb1u*lgu*mmgl*pow3(mmsb2)*pow4(g3))/(9.*mgl*(-mmgl + mmsb2)*pow3(
     mmsb1 - mmsb2)) + (32*lb2u*lgu*mmgl*pow3(mmsb2)*pow4(g3))/(9.*mgl*(-mmgl +
     mmsb2)*pow3(mmsb1 - mmsb2)) + (32*mmgl*pow2(lgu)*pow3(mmsb2)*pow4(g3))/(9.
     *mgl*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) - (32*mmgl*pow2(lgu)*pow3(mmsb2)
     *pow4(g3))/(9.*mgl*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) + (32*lb2u*lgu*
     mmgl*pow3(mmsb2)*pow4(g3))/(9.*mgl*(mmsb1 - mmsb2)*pow3(-mmgl + mmsb2)) -
     (16*mmgl*pow2(lb2u)*pow3(mmsb2)*pow4(g3))/(9.*mgl*(mmsb1 - mmsb2)*pow3(-
     mmgl + mmsb2)) - (16*mmgl*pow2(lgu)*pow3(mmsb2)*pow4(g3))/(9.*mgl*(mmsb1 -
     mmsb2)*pow3(-mmgl + mmsb2)) - (112*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow3(
     mmst1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*lt1u*mmgl*
     DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1
      - mmsb2)) + (16*ltu*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1)*pow4(g3))/(
     mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*lt1u*ltu*mmgl*DeltaInv(mmt,
     mmst1,mmgl)*pow3(mmst1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2))
     + (112*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1)*pow4(g3))/(3.*mgl*(mmsb1
     - mmsb2)*(-mmgl + mmsb2)) - (16*lt1u*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow3(
     mmst1)*pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*ltu*mmgl*
     DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1)*pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl
      + mmsb2)) + (16*lt1u*ltu*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1)*pow4(
     g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*zt2*DeltaInv(mmt,
     mmst1,mmgl)*pow3(mmst1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2))
     + (16*mmgl*zt2*DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1)*pow4(g3))/(3.*mgl*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (8*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(
     lt1u)*pow3(mmst1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (8*
     mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(lt1u)*pow3(mmst1)*pow4(g3))/(3.*mgl*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (8*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(
     ltu)*pow3(mmst1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (8*
     mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow3(mmst1)*pow4(g3))/(3.*mgl*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (112*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow3(
     mmst2)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*lt2u*mmgl*
     DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmsb1
      - mmsb2)) + (16*ltu*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2)*pow4(g3))/(
     mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*lt2u*ltu*mmgl*DeltaInv(mmt,
     mmst2,mmgl)*pow3(mmst2)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2))
     + (112*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2)*pow4(g3))/(3.*mgl*(mmsb1
     - mmsb2)*(-mmgl + mmsb2)) - (16*lt2u*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow3(
     mmst2)*pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*ltu*mmgl*
     DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2)*pow4(g3))/(mgl*(mmsb1 - mmsb2)*(-mmgl
      + mmsb2)) + (16*lt2u*ltu*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2)*pow4(
     g3))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*zt2*DeltaInv(mmt,
     mmst2,mmgl)*pow3(mmst2)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2))
     + (16*mmgl*zt2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2)*pow4(g3))/(3.*mgl*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (8*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(
     lt2u)*pow3(mmst2)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (8*
     mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(lt2u)*pow3(mmst2)*pow4(g3))/(3.*mgl*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (8*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(
     ltu)*pow3(mmst2)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (8*
     mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow3(mmst2)*pow4(g3))/(3.*mgl*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2))))/pow4(g3) + Xt*((Xb*((-224*mmgl*mmsb1*mmt
     *DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) +
     (32*lgu*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/((mmsb1 - mmsb2)
     *(mmst1 - mmst2)) + (32*ltu*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(
     g3))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*lgu*ltu*mmgl*mmsb1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) +
     (224*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(3.*(mmsb1 - mmsb2)
     *(mmst1 - mmst2)) - (32*lgu*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(
     g3))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*ltu*mmgl*mmsb2*mmt*DeltaInv(
     mmt,mmst1,mmgl)*pow4(g3))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*lgu*ltu*
     mmgl*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) + (224*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/
     (3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*lgu*mmsb1*mmst1*mmt*DeltaInv(
     mmt,mmst1,mmgl)*pow4(g3))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*lt1u*
     mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/((mmsb1 - mmsb2)*(mmst1
      - mmst2)) - (32*lgu*lt1u*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3
     ))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (64*ltu*mmsb1*mmst1*mmt*DeltaInv
     (mmt,mmst1,mmgl)*pow4(g3))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) + (64*lgu*ltu
     *mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) - (224*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/
     (3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*lgu*mmsb2*mmst1*mmt*DeltaInv(
     mmt,mmst1,mmgl)*pow4(g3))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*lt1u*
     mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/((mmsb1 - mmsb2)*(mmst1
      - mmst2)) + (32*lgu*lt1u*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3
     ))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (64*ltu*mmsb2*mmst1*mmt*DeltaInv
     (mmt,mmst1,mmgl)*pow4(g3))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) - (64*lgu*ltu
     *mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) - (32*mmgl*mmsb1*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow4(g3)
     )/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*mmgl*mmsb2*mmt*zt2*DeltaInv(
     mmt,mmst1,mmgl)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*mmsb1
     *mmst1*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) - (32*mmsb2*mmst1*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow4(g3
     ))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (224*mmgl*mmsb1*mmt*DeltaInv(mmt
     ,mmst2,mmgl)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*lgu*mmgl
     *mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/((mmsb1 - mmsb2)*(mmst1 -
     mmst2)) - (32*ltu*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/((
     mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*lgu*ltu*mmgl*mmsb1*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (224*mmgl*
     mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) + (32*lgu*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/((
     mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*ltu*mmgl*mmsb2*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow4(g3))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*lgu*ltu*mmgl
     *mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) - (224*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*(
     mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*lgu*mmsb1*mmst2*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow4(g3))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*lt2u*mmsb1*
     mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/((mmsb1 - mmsb2)*(mmst1 -
     mmst2)) + (32*lgu*lt2u*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/
     (3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (64*ltu*mmsb1*mmst2*mmt*DeltaInv(
     mmt,mmst2,mmgl)*pow4(g3))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) - (64*lgu*ltu*
     mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) + (224*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/
     (3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*lgu*mmsb2*mmst2*mmt*DeltaInv(
     mmt,mmst2,mmgl)*pow4(g3))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*lt2u*
     mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/((mmsb1 - mmsb2)*(mmst1
      - mmst2)) - (32*lgu*lt2u*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3
     ))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (64*ltu*mmsb2*mmst2*mmt*DeltaInv
     (mmt,mmst2,mmgl)*pow4(g3))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) + (64*lgu*ltu
     *mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) + (32*mmgl*mmsb1*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow4(g3)
     )/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmgl*mmsb2*mmt*zt2*DeltaInv(
     mmt,mmst2,mmgl)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmsb1
     *mmst2*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) + (32*mmsb2*mmst2*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow4(g3
     ))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*mmt*Fin3(mmt,mmsb1,mmst1,mmu
     )*pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (16*mmt
     *Fin3(mmt,mmsb1,mmst2,mmu)*pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) - (16*mmt*Fin3(mmt,mmsb2,mmst1,mmu)*pow4(g3))/(3.*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (16*mmt*Fin3(mmt,mmsb2,mmst2,mmu
     )*pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (16*mmt
     *Fin3(mmt,mmst1,mmgl,mmu)*pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) + (16*mmt*Fin3(mmt,mmst1,mmgl,mmu)*pow4(g3))/(3.*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (32*mmsb1*mmt*DeltaInv(mmt,mmst1
     ,mmgl)*Fin3(mmt,mmst1,mmgl,mmu)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) + (32*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu)*
     pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmsb1*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu)*pow4(g3))/(3.*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*mmsb2*mmst1*mmt*DeltaInv(mmt
     ,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl
      + mmsb2)*(mmst1 - mmst2)) + (16*mmt*Fin3(mmt,mmst2,mmgl,mmu)*pow4(g3))/(
     3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (16*mmt*Fin3(mmt,
     mmst2,mmgl,mmu)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 -
     mmst2)) + (32*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*
     pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmsb2*mmt*DeltaInv(
     mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) + (32*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,
     mmst2,mmgl,mmu)*pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) - (32*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl
     ,mmu)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (16
     *mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow4(g3))/(3.*(mmsb1 -
     mmsb2)*(mmst1 - mmst2)) + (16*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2
     (lgu)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*mmsb1*mmst1*mmt
     *DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) - (16*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow4(g3))
     /(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*mmgl*mmsb1*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow2(lgu)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (16
     *mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow4(g3))/(3.*(mmsb1 -
     mmsb2)*(mmst1 - mmst2)) - (16*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*
     pow2(lgu)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*mmsb2*mmst2
     *mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) - (16*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lt1u)*
     pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*mmsb2*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(lt1u)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) + (16*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lt2u)*pow4(g3)
     )/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (16*mmsb2*mmst2*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow2(lt2u)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (
     16*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow4(g3))/(3.*(mmsb1
     - mmsb2)*(mmst1 - mmst2)) + (16*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow2(ltu)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*mmsb1*mmst1
     *mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) - (32*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*
     pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*mmgl*mmsb1*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) - (16*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow4(g3))/
     (3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmsb1*mmst2*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow2(ltu)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32
     *mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow4(g3))/(3.*(mmsb1 -
     mmsb2)*(mmst1 - mmst2)) - (224*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*
     pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*lgu*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmsb1)*pow4(g3))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*
     ltu*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow4(g3))/((mmsb1 - mmsb2)*(
     mmst1 - mmst2)) - (32*lgu*ltu*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*
     pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (224*mmst1*mmt*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmsb1)*pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*
     (mmst1 - mmst2)) + (32*lgu*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*
     pow4(g3))/((-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*lt1u*
     mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow4(g3))/((-mmgl + mmsb1)*
     (mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*lgu*lt1u*mmst1*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmsb1)*pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) + (64*ltu*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*
     pow4(g3))/((-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (64*lgu*ltu*
     mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow4(g3))/(3.*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmt*zt2*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmsb1)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*
     mmst1*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow4(g3))/(3.*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (224*mmt*DeltaInv(mmt,mmst2,mmgl
     )*pow2(mmsb1)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*lgu*mmt
     *DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3))/((mmsb1 - mmsb2)*(mmst1 -
     mmst2)) - (32*ltu*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3))/((
     mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*lgu*ltu*mmt*DeltaInv(mmt,mmst2,mmgl)
     *pow2(mmsb1)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (224*mmst2*
     mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3))/(3.*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*lgu*mmst2*mmt*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmsb1)*pow4(g3))/((-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) + (32*lt2u*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3)
     )/((-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*lgu*lt2u*mmst2*
     mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3))/(3.*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)*(mmst1 - mmst2)) - (64*ltu*mmst2*mmt*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmsb1)*pow4(g3))/((-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) + (64*lgu*ltu*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(
     g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*mmt*zt2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1
     - mmst2)) + (32*mmst2*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3
     ))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*mmt*DeltaInv
     (mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu)*pow2(mmsb1)*pow4(g3))/(3.*(-mmgl
      + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmt*DeltaInv(mmt,mmst2,
     mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*pow2(mmsb1)*pow4(g3))/(3.*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)*(mmst1 - mmst2)) - (16*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     lgu)*pow2(mmsb1)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (16*
     mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow2(mmsb1)*pow4(g3))/(3.*(-
     mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow2(lgu)*pow2(mmsb1)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) + (16*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(mmsb1)*
     pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*mmst1
     *mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lt1u)*pow2(mmsb1)*pow4(g3))/(3.*(-mmgl
     + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (16*mmst2*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow2(lt2u)*pow2(mmsb1)*pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)*(mmst1 - mmst2)) - (16*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow2(
     mmsb1)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow2(mmsb1)*pow4(g3))/(3.*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*mmt*DeltaInv(mmt,mmst2,mmgl)
     *pow2(ltu)*pow2(mmsb1)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (
     32*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow2(mmsb1)*pow4(g3))/(3.*
     (-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmsb1*mmt*Fin3(mmt,
     mmsb1,mmst1,mmu)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)*pow2(-mmgl
     + mmsb1)) + (32*mmsb1*mmt*Fin3(mmt,mmsb1,mmst2,mmu)*pow4(g3))/(3.*(mmsb1 -
     mmsb2)*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (32*mmsb1*mmt*Fin3(mmt,mmst1
     ,mmgl,mmu)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) - (32*mmsb1*mmt*Fin3(mmt,mmst2,mmgl,mmu)*pow4(g3))/(3.*(mmsb1 -
     mmsb2)*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (224*mmt*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmsb2)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*lgu
     *mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(g3))/((mmsb1 - mmsb2)*(
     mmst1 - mmst2)) - (32*ltu*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(g3
     ))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*lgu*ltu*mmt*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmsb2)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (224*
     mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(g3))/(3.*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (32*lgu*mmst1*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmsb2)*pow4(g3))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1
     - mmst2)) + (32*lt1u*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(
     g3))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (32*lgu*lt1u*
     mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(g3))/(3.*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (64*ltu*mmst1*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmsb2)*pow4(g3))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1
     - mmst2)) + (64*lgu*ltu*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*
     pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (32*mmt*
     zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) + (32*mmst1*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*
     pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (224*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1
     - mmst2)) + (32*lgu*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow4(g3))/((
     mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*ltu*mmt*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmsb2)*pow4(g3))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*lgu*ltu*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1
     - mmst2)) - (224*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow4(g3))/
     (3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (32*lgu*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow4(g3))/((mmsb1 - mmsb2)*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) - (32*lt2u*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2
     (mmsb2)*pow4(g3))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (32*
     lgu*lt2u*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow4(g3))/(3.*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (64*ltu*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow4(g3))/((mmsb1 - mmsb2)*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) - (64*lgu*ltu*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmsb2)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2))
     - (32*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow4(g3))/(3.*(mmsb1 -
     mmsb2)*(mmst1 - mmst2)) - (32*mmst2*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmsb2)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (
     32*mmt*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu)*pow2(mmsb2)*pow4(
     g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (32*mmt*
     DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*pow2(mmsb2)*pow4(g3))/(
     3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (16*mmt*DeltaInv(mmt
     ,mmst1,mmgl)*pow2(lgu)*pow2(mmsb2)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) + (16*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow2(mmsb2)*
     pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (16*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(mmsb2)*pow4(g3))/(3.*(mmsb1 -
     mmsb2)*(mmst1 - mmst2)) - (16*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)
     *pow2(mmsb2)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)
     ) - (16*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lt1u)*pow2(mmsb2)*pow4(g3)
     )/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (16*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(lt2u)*pow2(mmsb2)*pow4(g3))/(3.*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (16*mmt*DeltaInv(mmt,mmst1,mmgl)
     *pow2(ltu)*pow2(mmsb2)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (
     32*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow2(mmsb2)*pow4(g3))/(3.*
     (mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (16*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow2(ltu)*pow2(mmsb2)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) - (32*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow2(mmsb2)*
     pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (32*mmsb2
     *mmt*Fin3(mmt,mmsb2,mmst1,mmu)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb2)) - (32*mmsb2*mmt*Fin3(mmt,mmsb2,mmst2,mmu)*pow4
     (g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (32*mmsb2
     *mmt*Fin3(mmt,mmst1,mmgl,mmu)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2
     )*pow2(-mmgl + mmsb2)) + (32*mmsb2*mmt*Fin3(mmt,mmst2,mmgl,mmu)*pow4(g3))/
     (3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (32*lt1u*mmsb1*
     mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow4(g3))/((-mmgl + mmsb1)*(mmsb1
      - mmsb2)*(mmst1 - mmst2)) - (32*lgu*lt1u*mmsb1*mmt*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmst1)*pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) - (32*ltu*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow4(g3))
     /((-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*lgu*ltu*mmsb1*mmt
     *DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1
      - mmsb2)*(mmst1 - mmst2)) - (32*lt1u*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmst1)*pow4(g3))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) +
     (32*lgu*lt1u*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow4(g3))/(3.*
     (mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (32*ltu*mmsb2*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow4(g3))/((mmsb1 - mmsb2)*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) - (32*lgu*ltu*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmst1)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2))
     - (16*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lt1u)*pow2(mmst1)*pow4(g3))/
     (3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*mmsb2*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(lt1u)*pow2(mmst1)*pow4(g3))/(3.*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (16*mmsb1*mmt*DeltaInv(mmt,mmst1
     ,mmgl)*pow2(ltu)*pow2(mmst1)*pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)
     *(mmst1 - mmst2)) - (16*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow2(
     mmst1)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (
     32*lt2u*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/((-mmgl +
     mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*lgu*lt2u*mmsb1*mmt*DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*
     (mmst1 - mmst2)) + (32*ltu*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*
     pow4(g3))/((-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*lgu*ltu*
     mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(3.*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*lt2u*mmsb2*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmst2)*pow4(g3))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1
     - mmst2)) - (32*lgu*lt2u*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*
     pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (32*ltu*
     mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/((mmsb1 - mmsb2)*
     (-mmgl + mmsb2)*(mmst1 - mmst2)) + (32*lgu*ltu*mmsb2*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(
     mmst1 - mmst2)) + (16*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lt2u)*pow2(
     mmst2)*pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (
     16*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lt2u)*pow2(mmst2)*pow4(g3))/(3.
     *(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (16*mmsb1*mmt*DeltaInv
     (mmt,mmst2,mmgl)*pow2(ltu)*pow2(mmst2)*pow4(g3))/(3.*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*
     pow2(ltu)*pow2(mmst2)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1
      - mmst2)) + (224*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(3.*(
     mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*lgu*mmsb1*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmt)*pow4(g3))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*ltu*mmsb1*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/((mmsb1 - mmsb2)*(mmst1 -
     mmst2)) + (32*lgu*ltu*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(
     3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (224*mmsb2*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmt)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*lgu*mmsb2*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/((mmsb1 - mmsb2)*(mmst1 -
     mmst2)) + (32*ltu*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/((
     mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*lgu*ltu*mmsb2*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmt)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (448*
     mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(3.*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*lt1u*mmsb1*mmst1*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/((-mmgl + mmsb1)*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) + (32*lgu*lt1u*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmt)*pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (96*
     ltu*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/((-mmgl +
     mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*lgu*ltu*mmsb1*mmst1*DeltaInv
     (mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) - (64*lt1u*ltu*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmt)*pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (448
     *mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(3.*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (32*lt1u*mmsb2*mmst1*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(
     mmst1 - mmst2)) - (32*lgu*lt1u*mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmt)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (96*
     ltu*mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/((mmsb1 -
     mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (32*lgu*ltu*mmsb2*mmst1*DeltaInv
     (mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(
     mmst1 - mmst2)) + (64*lt1u*ltu*mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmt)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (32*
     mmsb1*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(3.*(mmsb1 - mmsb2)
     *(mmst1 - mmst2)) - (32*mmsb2*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(
     g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (64*mmsb1*mmst1*zt2*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) + (64*mmsb2*mmst1*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*
     pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (224*
     mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) + (32*lgu*mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3
     ))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*ltu*mmsb1*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmt)*pow4(g3))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*lgu*ltu*
     mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) + (224*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/
     (3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*lgu*mmsb2*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmt)*pow4(g3))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*ltu*
     mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/((mmsb1 - mmsb2)*(mmst1
      - mmst2)) + (32*lgu*ltu*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3)
     )/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (448*mmsb1*mmst2*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmt)*pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1
      - mmst2)) - (32*lt2u*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(
     g3))/((-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*lgu*lt2u*
     mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/(3.*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (96*ltu*mmsb1*mmst2*DeltaInv(mmt
     ,mmst2,mmgl)*pow2(mmt)*pow4(g3))/((-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) + (32*lgu*ltu*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(
     g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (64*lt2u*ltu*
     mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/(3.*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (448*mmsb2*mmst2*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmt)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1
      - mmst2)) + (32*lt2u*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(
     g3))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (32*lgu*lt2u*
     mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/(3.*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (96*ltu*mmsb2*mmst2*DeltaInv(mmt
     ,mmst2,mmgl)*pow2(mmt)*pow4(g3))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 -
     mmst2)) - (32*lgu*ltu*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(
     g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (64*lt2u*ltu*
     mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/(3.*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (32*mmsb1*zt2*DeltaInv(mmt,mmst2
     ,mmgl)*pow2(mmt)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*
     mmsb2*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/(3.*(mmsb1 - mmsb2)
     *(mmst1 - mmst2)) + (64*mmsb1*mmst2*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)
     *pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (64*
     mmsb2*mmst2*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/(3.*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (32*mmsb1*DeltaInv(mmt,mmst1,
     mmgl)*Fin3(mmt,mmst1,mmgl,mmu)*pow2(mmt)*pow4(g3))/(3.*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*mmsb2*DeltaInv(mmt,mmst1,mmgl)*Fin3(
     mmt,mmst1,mmgl,mmu)*pow2(mmt)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2
     )*(mmst1 - mmst2)) + (32*mmsb1*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,
     mmgl,mmu)*pow2(mmt)*pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) - (32*mmsb2*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*pow2
     (mmt)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (16
     *mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow2(mmt)*pow4(g3))/(3.*(mmsb1 -
     mmsb2)*(mmst1 - mmst2)) - (16*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*
     pow2(mmt)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (16*mmsb1*
     DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(mmt)*pow4(g3))/(3.*(mmsb1 - mmsb2)
     *(mmst1 - mmst2)) + (16*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(mmt)
     *pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (16*mmsb1*mmst1*DeltaInv
     (mmt,mmst1,mmgl)*pow2(lt1u)*pow2(mmt)*pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1
      - mmsb2)*(mmst1 - mmst2)) + (16*mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2
     (lt1u)*pow2(mmt)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 -
     mmst2)) + (16*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(lt2u)*pow2(mmt)*
     pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (16*mmsb2
     *mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(lt2u)*pow2(mmt)*pow4(g3))/(3.*(mmsb1
     - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (16*mmsb1*DeltaInv(mmt,mmst1,
     mmgl)*pow2(ltu)*pow2(mmt)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) -
     (16*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow2(mmt)*pow4(g3))/(3.*(
     mmsb1 - mmsb2)*(mmst1 - mmst2)) - (16*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)
     *pow2(ltu)*pow2(mmt)*pow4(g3))/((-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) + (16*mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow2(mmt)*
     pow4(g3))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (16*mmsb1*
     DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow2(mmt)*pow4(g3))/(3.*(mmsb1 - mmsb2)
     *(mmst1 - mmst2)) + (16*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow2(mmt)
     *pow4(g3))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*mmsb1*mmst2*DeltaInv
     (mmt,mmst2,mmgl)*pow2(ltu)*pow2(mmt)*pow4(g3))/((-mmgl + mmsb1)*(mmsb1 -
     mmsb2)*(mmst1 - mmst2)) - (16*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(
     ltu)*pow2(mmt)*pow4(g3))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2))
     - (224*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(mmt)*pow4(g3))/(3.*(-mmgl
      + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*lgu*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmsb1)*pow2(mmt)*pow4(g3))/((-mmgl + mmsb1)*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) + (32*ltu*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(mmt)*
     pow4(g3))/((-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*lgu*ltu*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(mmt)*pow4(g3))/(3.*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*zt2*DeltaInv(mmt,mmst1,mmgl)
     *pow2(mmsb1)*pow2(mmt)*pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) + (224*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(mmt)*pow4
     (g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*lgu*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(mmt)*pow4(g3))/((-mmgl + mmsb1)*
     (mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*ltu*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmsb1)*pow2(mmt)*pow4(g3))/((-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2
     )) + (32*lgu*ltu*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(mmt)*pow4(g3))/
     (3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*zt2*DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(mmt)*pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1
      - mmsb2)*(mmst1 - mmst2)) - (16*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow2(
     mmsb1)*pow2(mmt)*pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) + (16*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(mmsb1)*pow2(mmt)*
     pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (16*
     DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow2(mmsb1)*pow2(mmt)*pow4(g3))/(3.*(-
     mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*DeltaInv(mmt,mmst2,
     mmgl)*pow2(ltu)*pow2(mmsb1)*pow2(mmt)*pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1
      - mmsb2)*(mmst1 - mmst2)) + (224*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*
     pow2(mmt)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) -
     (32*lgu*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(mmt)*pow4(g3))/((mmsb1 -
     mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (32*ltu*DeltaInv(mmt,mmst1,mmgl)
     *pow2(mmsb2)*pow2(mmt)*pow4(g3))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 -
     mmst2)) + (32*lgu*ltu*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(mmt)*pow4(
     g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (32*zt2*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(mmt)*pow4(g3))/(3.*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (224*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmsb2)*pow2(mmt)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1
      - mmst2)) + (32*lgu*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(mmt)*pow4(
     g3))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (32*ltu*DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(mmt)*pow4(g3))/((mmsb1 - mmsb2)*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) - (32*lgu*ltu*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)
     *pow2(mmt)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2))
     - (32*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(mmt)*pow4(g3))/(3.*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (16*DeltaInv(mmt,mmst1,
     mmgl)*pow2(lgu)*pow2(mmsb2)*pow2(mmt)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl
      + mmsb2)*(mmst1 - mmst2)) - (16*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(
     mmsb2)*pow2(mmt)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 -
     mmst2)) + (16*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow2(mmsb2)*pow2(mmt)*
     pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (16*
     DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow2(mmsb2)*pow2(mmt)*pow4(g3))/(3.*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (224*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow3(mmsb1)*pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) - (32*lgu*mmt*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1)*pow4(g3
     ))/((-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*ltu*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1)*pow4(g3))/((-mmgl + mmsb1)*(mmsb1 -
     mmsb2)*(mmst1 - mmst2)) + (32*lgu*ltu*mmt*DeltaInv(mmt,mmst1,mmgl)*pow3(
     mmsb1)*pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (
     32*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1)*pow4(g3))/(3.*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (224*mmt*DeltaInv(mmt,mmst2,mmgl
     )*pow3(mmsb1)*pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2
     )) + (32*lgu*mmt*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1)*pow4(g3))/((-mmgl +
     mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*ltu*mmt*DeltaInv(mmt,mmst2,
     mmgl)*pow3(mmsb1)*pow4(g3))/((-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) - (32*lgu*ltu*mmt*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1)*pow4(g3))/(
     3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmt*zt2*DeltaInv
     (mmt,mmst2,mmgl)*pow3(mmsb1)*pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)
     *(mmst1 - mmst2)) + (16*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow3(mmsb1)
     *pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (16*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow3(mmsb1)*pow4(g3))/(3.*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*mmt*DeltaInv(mmt,mmst1,mmgl)
     *pow2(ltu)*pow3(mmsb1)*pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) - (16*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow3(mmsb1)*
     pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (224*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl
     + mmsb2)*(mmst1 - mmst2)) + (32*lgu*mmt*DeltaInv(mmt,mmst1,mmgl)*pow3(
     mmsb2)*pow4(g3))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (32*
     ltu*mmt*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2)*pow4(g3))/((mmsb1 - mmsb2)*(-
     mmgl + mmsb2)*(mmst1 - mmst2)) - (32*lgu*ltu*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow3(mmsb2)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2))
     - (32*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2)*pow4(g3))/(3.*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (224*mmt*DeltaInv(mmt,mmst2,mmgl
     )*pow3(mmsb2)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2
     )) - (32*lgu*mmt*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2)*pow4(g3))/((mmsb1 -
     mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (32*ltu*mmt*DeltaInv(mmt,mmst2,
     mmgl)*pow3(mmsb2)*pow4(g3))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 -
     mmst2)) + (32*lgu*ltu*mmt*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2)*pow4(g3))/(
     3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (32*mmt*zt2*DeltaInv
     (mmt,mmst2,mmgl)*pow3(mmsb2)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)
     *(mmst1 - mmst2)) - (16*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow3(mmsb2)
     *pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (16*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow3(mmsb2)*pow4(g3))/(3.*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (16*mmt*DeltaInv(mmt,mmst1,mmgl)
     *pow2(ltu)*pow3(mmsb2)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(
     mmst1 - mmst2)) + (16*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow3(mmsb2)*
     pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2))))/pow4(g3)
     + ((112*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(3.*mgl*(mmst1 -
     mmst2)) - (16*lgu*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(mgl*(
     mmst1 - mmst2)) - (16*ltu*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3)
     )/(mgl*(mmst1 - mmst2)) + (16*lgu*ltu*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,
     mmgl)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)) + (112*mmgl*mmsb2*mmt*DeltaInv(
     mmt,mmst1,mmgl)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)) - (16*lgu*mmgl*mmsb2*
     mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(mgl*(mmst1 - mmst2)) - (16*ltu*
     mmgl*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(mgl*(mmst1 - mmst2)) +
     (16*lgu*ltu*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(3.*mgl*(
     mmst1 - mmst2)) - (112*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(
     3.*mgl*(mmst1 - mmst2)) + (16*lgu*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow4(g3))/(mgl*(mmst1 - mmst2)) - (16*lt1u*mmgl*mmst1*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow4(g3))/(mgl*(mmst1 - mmst2)) + (16*lgu*lt1u*mmgl*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)) + (32*ltu*mmgl
     *mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(mgl*(mmst1 - mmst2)) - (32*
     lgu*ltu*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(3.*mgl*(mmst1 -
     mmst2)) + (112*mmgl*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(3.
     *mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (16*lgu*mmgl*mmsb1*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) +
     (16*lt1u*mmgl*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(mgl*(-
     mmgl + mmsb1)*(mmst1 - mmst2)) - (16*lgu*lt1u*mmgl*mmsb1*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)
     ) - (32*ltu*mmgl*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(mgl*(
     -mmgl + mmsb1)*(mmst1 - mmst2)) + (32*lgu*ltu*mmgl*mmsb1*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)
     ) + (112*mmgl*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(3.*mgl*(
     -mmgl + mmsb2)*(mmst1 - mmst2)) - (16*lgu*mmgl*mmsb2*mmst1*mmt*DeltaInv(
     mmt,mmst1,mmgl)*pow4(g3))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (16*lt1u
     *mmgl*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(mgl*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) - (16*lgu*lt1u*mmgl*mmsb2*mmst1*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (32*ltu*
     mmgl*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(mgl*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) + (32*lgu*ltu*mmgl*mmsb2*mmst1*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (16*mmgl*
     mmsb1*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(3.*mgl*(mmst1 - mmst2))
     + (16*mmgl*mmsb2*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(3.*mgl*(mmst1
      - mmst2)) - (16*mmgl*mmst1*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(3.
     *mgl*(mmst1 - mmst2)) + (16*mmgl*mmsb1*mmst1*mmt*zt2*DeltaInv(mmt,mmst1,
     mmgl)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (16*mmgl*mmsb2*
     mmst1*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(
     mmst1 - mmst2)) - (112*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(
     3.*mgl*(mmst1 - mmst2)) + (16*lgu*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*
     pow4(g3))/(mgl*(mmst1 - mmst2)) + (16*ltu*mmgl*mmsb1*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow4(g3))/(mgl*(mmst1 - mmst2)) - (16*lgu*ltu*mmgl*mmsb1*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)) - (112*mmgl*
     mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)) + (
     16*lgu*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(mgl*(mmst1 -
     mmst2)) + (16*ltu*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(mgl*(
     mmst1 - mmst2)) - (16*lgu*ltu*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4
     (g3))/(3.*mgl*(mmst1 - mmst2)) + (112*mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,
     mmgl)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)) - (16*lgu*mmgl*mmst2*mmt*DeltaInv
     (mmt,mmst2,mmgl)*pow4(g3))/(mgl*(mmst1 - mmst2)) + (16*lt2u*mmgl*mmst2*mmt
     *DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(mgl*(mmst1 - mmst2)) - (16*lgu*lt2u*
     mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*mgl*(mmst1 - mmst2))
     - (32*ltu*mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(mgl*(mmst1 -
     mmst2)) + (32*lgu*ltu*mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(
     3.*mgl*(mmst1 - mmst2)) - (112*mmgl*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,
     mmgl)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (16*lgu*mmgl*
     mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(
     mmst1 - mmst2)) - (16*lt2u*mmgl*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*
     pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (16*lgu*lt2u*mmgl*mmsb1*
     mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(
     mmst1 - mmst2)) + (32*ltu*mmgl*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*
     pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (32*lgu*ltu*mmgl*mmsb1*
     mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(
     mmst1 - mmst2)) - (112*mmgl*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(
     g3))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (16*lgu*mmgl*mmsb2*mmst2*
     mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2
     )) - (16*lt2u*mmgl*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(mgl
     *(-mmgl + mmsb2)*(mmst1 - mmst2)) + (16*lgu*lt2u*mmgl*mmsb2*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)
     ) + (32*ltu*mmgl*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(mgl*(
     -mmgl + mmsb2)*(mmst1 - mmst2)) - (32*lgu*ltu*mmgl*mmsb2*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)
     ) - (16*mmgl*mmsb1*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*mgl*(
     mmst1 - mmst2)) - (16*mmgl*mmsb2*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow4(g3)
     )/(3.*mgl*(mmst1 - mmst2)) + (16*mmgl*mmst2*mmt*zt2*DeltaInv(mmt,mmst2,
     mmgl)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)) - (16*mmgl*mmsb1*mmst2*mmt*zt2*
     DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)
     ) - (16*mmgl*mmsb2*mmst2*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*
     mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (16*mmgl*mmt*DeltaInv(mmt,mmst1,
     mmgl)*Fin3(mmt,mmst1,mmgl,mmu)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)) - (16*
     mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu)*pow4(g3))
     /(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (16*mmgl*mmsb2*mmt*DeltaInv(
     mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)
     *(mmst1 - mmst2)) + (8*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,
     mmst1,mmgl,mmu)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (8*
     mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu)*pow4(g3))
     /(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (16*mmgl*mmt*DeltaInv(mmt,
     mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)) +
     (16*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*pow4(
     g3))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (16*mmgl*mmsb2*mmt*
     DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*pow4(g3))/(3.*mgl*(-mmgl
      + mmsb2)*(mmst1 - mmst2)) - (8*mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*
     Fin3(mmt,mmst2,mmgl,mmu)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)
     ) - (8*mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*
     pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (8*mmgl*mmsb1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)) + (8
     *mmgl*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow4(g3))/(3.*mgl*(
     mmst1 - mmst2)) - (8*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*
     pow4(g3))/(3.*mgl*(mmst1 - mmst2)) + (8*mmgl*mmsb1*mmst1*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow2(lgu)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) +
     (8*mmgl*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow4(g3))/(3.*
     mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (8*mmgl*mmsb1*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow2(lgu)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)) - (8*mmgl*mmsb2*
     mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow4(g3))/(3.*mgl*(mmst1 - mmst2))
     + (8*mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow4(g3))/(3.*mgl*(
     mmst1 - mmst2)) - (8*mmgl*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(
     lgu)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (8*mmgl*mmsb2*
     mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow4(g3))/(3.*mgl*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) + (8*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     lt1u)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)) - (8*mmgl*mmsb1*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(lt1u)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(
     mmst1 - mmst2)) - (8*mmgl*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     lt1u)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (8*mmgl*mmst2*
     mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lt2u)*pow4(g3))/(3.*mgl*(mmst1 - mmst2))
     + (8*mmgl*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lt2u)*pow4(g3))/(
     3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (8*mmgl*mmsb2*mmst2*mmt*DeltaInv
     (mmt,mmst2,mmgl)*pow2(lt2u)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 -
     mmst2)) + (8*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow4(g3))/(
     3.*mgl*(mmst1 - mmst2)) + (8*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     ltu)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)) - (16*mmgl*mmst1*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow2(ltu)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)) + (16*mmgl*mmsb1*
     mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow4(g3))/(3.*mgl*(-mmgl +
     mmsb1)*(mmst1 - mmst2)) + (16*mmgl*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl
     )*pow2(ltu)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (8*mmgl*
     mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow4(g3))/(3.*mgl*(mmst1 -
     mmst2)) - (8*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow4(g3))/(
     3.*mgl*(mmst1 - mmst2)) + (16*mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2
     (ltu)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)) - (16*mmgl*mmsb1*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(
     mmst1 - mmst2)) - (16*mmgl*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(
     ltu)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (112*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmgl)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)) - (
     16*lgu*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmgl)*pow4(g3))/(mgl*(mmst1 -
     mmst2)) - (16*ltu*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmgl)*pow4(g3))/(mgl*(
     mmst1 - mmst2)) + (16*lgu*ltu*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmgl)*pow4
     (g3))/(3.*mgl*(mmst1 - mmst2)) + (16*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2
     (mmgl)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)) - (112*mmt*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmgl)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)) + (16*lgu*mmt*DeltaInv
     (mmt,mmst2,mmgl)*pow2(mmgl)*pow4(g3))/(mgl*(mmst1 - mmst2)) + (16*ltu*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmgl)*pow4(g3))/(mgl*(mmst1 - mmst2)) - (16*
     lgu*ltu*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmgl)*pow4(g3))/(3.*mgl*(mmst1 -
     mmst2)) - (16*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmgl)*pow4(g3))/(3.*
     mgl*(mmst1 - mmst2)) + (8*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow2(mmgl
     )*pow4(g3))/(3.*mgl*(mmst1 - mmst2)) - (8*mmt*DeltaInv(mmt,mmst2,mmgl)*
     pow2(lgu)*pow2(mmgl)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)) + (8*mmt*DeltaInv(
     mmt,mmst1,mmgl)*pow2(ltu)*pow2(mmgl)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)) -
     (8*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow2(mmgl)*pow4(g3))/(3.*mgl*(
     mmst1 - mmst2)) - (56*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow4(
     g3))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (24*lgu*mmgl*mmt*DeltaInv(mmt
     ,mmst1,mmgl)*pow2(mmsb1)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) +
     (24*ltu*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow4(g3))/(mgl*(-
     mmgl + mmsb1)*(mmst1 - mmst2)) - (8*lgu*ltu*mmgl*mmt*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmsb1)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (8*
     mmgl*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow4(g3))/(mgl*(-mmgl +
     mmsb1)*(mmst1 - mmst2)) + (56*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1
     )*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (24*lgu*mmgl*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmst1
      - mmst2)) - (24*ltu*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3
     ))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (8*lgu*ltu*mmgl*mmt*DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)
     ) + (8*mmgl*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3))/(mgl*(-
     mmgl + mmsb1)*(mmst1 - mmst2)) - (4*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2
     (lgu)*pow2(mmsb1)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (4*
     mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(mmsb1)*pow4(g3))/(mgl*(-
     mmgl + mmsb1)*(mmst1 - mmst2)) - (4*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2
     (ltu)*pow2(mmsb1)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (4*
     mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow2(mmsb1)*pow4(g3))/(mgl*(-
     mmgl + mmsb1)*(mmst1 - mmst2)) + (8*mmgl*mmt*Fin3(mmt,mmsb1,mmst1,mmu)*
     pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmsb1*
     mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1)*Fin3(mmt,mmsb1,mmst1,mmu)*pow4(g3))/(
     3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmt*Fin3(mmt,mmsb1,
     mmst2,mmu)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*
     mmgl*mmsb1*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2)*Fin3(mmt,mmsb1,mmst2,mmu)*
     pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmt*Fin3(
     mmt,mmst1,mmgl,mmu)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1))
     - (8*mmgl*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu
     )*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmt*
     Fin3(mmt,mmst2,mmgl,mmu)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) + (8*mmgl*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,
     mmgl,mmu)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (56*
     mmgl*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*lb1u*mmgl*mmst1*mmt*DeltaInv(mmt,
     mmsb1,mmst1)*pow2(mmsb1)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1
     )) - (8*lt1u*mmgl*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow4(g3)
     )/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*lb1u*lt1u*mmgl*mmst1*mmt*
     DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb1)) + (16*ltu*mmgl*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1)*
     pow2(mmsb1)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (16*lb1u
     *ltu*mmgl*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow4(g3))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmst1*mmt*zt2*DeltaInv(
     mmt,mmsb1,mmst1)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl
     + mmsb1)) + (56*mmgl*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*pow4(
     g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*lb1u*mmgl*mmst2*mmt
     *DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2
     (-mmgl + mmsb1)) + (8*lt2u*mmgl*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(
     mmsb1)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*lb1u*lt2u*
     mmgl*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (16*ltu*mmgl*mmst2*mmt*DeltaInv(mmt,
     mmsb1,mmst2)*pow2(mmsb1)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1
     )) + (16*lb1u*ltu*mmgl*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*
     pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmst2*mmt
     *zt2*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb1)) - (56*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmsb1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*
     lgu*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow4(g3))/(mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*lt1u*mmgl*mmst1*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmsb1)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)
     ) + (8*lgu*lt1u*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow4(
     g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (16*ltu*mmgl*mmst1*mmt
     *DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(
     -mmgl + mmsb1)) - (16*lgu*ltu*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2
     (mmsb1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*
     mmst1*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (56*mmgl*mmst2*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) - (8*lgu*mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(
     g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*lt2u*mmgl*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-
     mmgl + mmsb1)) - (8*lgu*lt2u*mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmsb1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (16*ltu*
     mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3))/(mgl*(mmst1
     - mmst2)*pow2(-mmgl + mmsb1)) + (16*lgu*ltu*mmgl*mmst2*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) + (8*mmgl*mmst2*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(
     g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmt*DeltaInv(
     mmt,mmsb1,mmst1)*Fin3(mmt,mmsb1,mmst1,mmu)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmt*DeltaInv(mmt,mmsb1,mmst2
     )*Fin3(mmt,mmsb1,mmst2,mmu)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb1)) + (8*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1
     ,mmgl,mmu)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) - (8*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*
     pow2(mmsb1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (4*
     mmgl*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(lb1u)*pow2(mmsb1)*pow4(g3))/
     (3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (4*mmgl*mmst2*mmt*DeltaInv(
     mmt,mmsb1,mmst2)*pow2(lb1u)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb1)) - (4*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu
     )*pow2(mmsb1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (4*
     mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(mmsb1)*pow4(g3))/(
     3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (4*mmgl*mmst1*mmt*DeltaInv(
     mmt,mmsb1,mmst1)*pow2(lt1u)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb1)) + (4*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     lt1u)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) -
     (4*mmgl*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(lt2u)*pow2(mmsb1)*pow4(g3
     ))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (4*mmgl*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(lt2u)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1)*
     pow2(ltu)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1
     )) + (8*mmgl*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(ltu)*pow2(mmsb1)*
     pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmst1*mmt
     *DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*
     pow2(ltu)*pow2(mmsb1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1
     )) - (56*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(g3))/(mgl*(-
     mmgl + mmsb2)*(mmst1 - mmst2)) + (24*lgu*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)
     *pow2(mmsb2)*pow4(g3))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (24*ltu*
     mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(g3))/(mgl*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) - (8*lgu*ltu*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmsb2)*pow4(g3))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (8*mmgl*mmt*
     zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(g3))/(mgl*(-mmgl + mmsb2)*(
     mmst1 - mmst2)) + (56*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow4(
     g3))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (24*lgu*mmgl*mmt*DeltaInv(mmt
     ,mmst2,mmgl)*pow2(mmsb2)*pow4(g3))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) -
     (24*ltu*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow4(g3))/(mgl*(-
     mmgl + mmsb2)*(mmst1 - mmst2)) + (8*lgu*ltu*mmgl*mmt*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmsb2)*pow4(g3))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (8*
     mmgl*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow4(g3))/(mgl*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) - (4*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*
     pow2(mmsb2)*pow4(g3))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (4*mmgl*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(mmsb2)*pow4(g3))/(mgl*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) - (4*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*
     pow2(mmsb2)*pow4(g3))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (4*mmgl*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow2(mmsb2)*pow4(g3))/(mgl*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) + (8*mmgl*mmt*Fin3(mmt,mmsb2,mmst1,mmu)*pow4(g3))/
     (3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmsb2*mmst1*mmt*
     DeltaInv(mmt,mmsb2,mmst1)*Fin3(mmt,mmsb2,mmst1,mmu)*pow4(g3))/(3.*mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmt*Fin3(mmt,mmsb2,mmst2,mmu
     )*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmsb2*
     mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2)*Fin3(mmt,mmsb2,mmst2,mmu)*pow4(g3))/(
     3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmt*Fin3(mmt,mmst1,
     mmgl,mmu)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl
     *mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu)*pow4(g3
     ))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmt*Fin3(mmt,
     mmst2,mmgl,mmu)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (
     8*mmgl*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*
     pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (56*mmgl*mmst1*
     mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2
     )*pow2(-mmgl + mmsb2)) + (8*lb2u*mmgl*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1)*
     pow2(mmsb2)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*lt1u*
     mmgl*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2)*pow4(g3))/(mgl*(mmst1
      - mmst2)*pow2(-mmgl + mmsb2)) + (8*lb2u*lt1u*mmgl*mmst1*mmt*DeltaInv(mmt,
     mmsb2,mmst1)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb2)) + (16*ltu*mmgl*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2)*
     pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (16*lb2u*ltu*mmgl*
     mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmst1*mmt*zt2*DeltaInv(mmt,mmsb2,
     mmst1)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2))
     + (56*mmgl*mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow4(g3))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*lb2u*mmgl*mmst2*mmt*DeltaInv
     (mmt,mmsb2,mmst2)*pow2(mmsb2)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb2)) + (8*lt2u*mmgl*mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*
     pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*lb2u*lt2u*mmgl*
     mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb2)) - (16*ltu*mmgl*mmst2*mmt*DeltaInv(mmt,mmsb2,
     mmst2)*pow2(mmsb2)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (
     16*lb2u*ltu*mmgl*mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow4(g3))
     /(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmst2*mmt*zt2*
     DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) - (56*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmsb2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*lgu*
     mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(g3))/(mgl*(mmst1
     - mmst2)*pow2(-mmgl + mmsb2)) - (8*lt1u*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmsb2)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8
     *lgu*lt1u*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(g3))/(
     3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (16*ltu*mmgl*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-
     mmgl + mmsb2)) - (16*lgu*ltu*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmsb2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*
     mmst1*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (56*mmgl*mmst2*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb2)) - (8*lgu*mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow4(
     g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*lt2u*mmgl*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-
     mmgl + mmsb2)) - (8*lgu*lt2u*mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmsb2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (16*ltu*
     mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow4(g3))/(mgl*(mmst1
     - mmst2)*pow2(-mmgl + mmsb2)) + (16*lgu*ltu*mmgl*mmst2*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb2)) + (8*mmgl*mmst2*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow4(
     g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmt*DeltaInv(
     mmt,mmsb2,mmst1)*Fin3(mmt,mmsb2,mmst1,mmu)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmt*DeltaInv(mmt,mmsb2,mmst2
     )*Fin3(mmt,mmsb2,mmst2,mmu)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) + (8*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1
     ,mmgl,mmu)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb2)) - (8*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*
     pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (4*
     mmgl*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(lb2u)*pow2(mmsb2)*pow4(g3))/
     (3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (4*mmgl*mmst2*mmt*DeltaInv(
     mmt,mmsb2,mmst2)*pow2(lb2u)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) - (4*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu
     )*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (4*
     mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(mmsb2)*pow4(g3))/(
     3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (4*mmgl*mmst1*mmt*DeltaInv(
     mmt,mmsb2,mmst1)*pow2(lt1u)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) + (4*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     lt1u)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) -
     (4*mmgl*mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(lt2u)*pow2(mmsb2)*pow4(g3
     ))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (4*mmgl*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(lt2u)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1)*
     pow2(ltu)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2
     )) + (8*mmgl*mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(ltu)*pow2(mmsb2)*
     pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmst1*mmt
     *DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*
     pow2(ltu)*pow2(mmsb2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2
     )) - (8*lt1u*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow4(g3))/(mgl*
     (-mmgl + mmsb1)*(mmst1 - mmst2)) + (8*lgu*lt1u*mmgl*mmt*DeltaInv(mmt,mmst1
     ,mmgl)*pow2(mmst1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (8
     *ltu*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow4(g3))/(mgl*(-mmgl +
     mmsb1)*(mmst1 - mmst2)) - (8*lgu*ltu*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmst1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (8*lt1u*
     mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow4(g3))/(mgl*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) + (8*lgu*lt1u*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmst1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (8*ltu*
     mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow4(g3))/(mgl*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) - (8*lgu*ltu*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmst1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (4*mmgl*
     mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lt1u)*pow2(mmst1)*pow4(g3))/(3.*mgl*(-
     mmgl + mmsb1)*(mmst1 - mmst2)) + (4*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2
     (lt1u)*pow2(mmst1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (4
     *mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow2(mmst1)*pow4(g3))/(3.*mgl
     *(-mmgl + mmsb1)*(mmst1 - mmst2)) - (4*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow2(ltu)*pow2(mmst1)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) +
     (8*lt1u*mmgl*mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1)*pow4(g3))/(
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*lb1u*lt1u*mmgl*mmsb1*mmt*
     DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb1)) - (8*ltu*mmgl*mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst1)*
     pow2(mmst1)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*lb1u*
     ltu*mmgl*mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1)*pow4(g3))/(3.*mgl
     *(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*lt1u*mmgl*mmsb1*mmt*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmst1)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) - (8*lgu*lt1u*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*
     pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*ltu*mmgl*mmsb1
     *mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow4(g3))/(mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb1)) + (8*lgu*ltu*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmst1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (4*
     mmgl*mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(lt1u)*pow2(mmst1)*pow4(g3))/
     (3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (4*mmgl*mmsb1*mmt*DeltaInv(
     mmt,mmst1,mmgl)*pow2(lt1u)*pow2(mmst1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb1)) + (4*mmgl*mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(
     ltu)*pow2(mmst1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) +
     (4*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow2(mmst1)*pow4(g3))
     /(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*lt1u*mmgl*mmsb2*mmt*
     DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(
     -mmgl + mmsb2)) - (8*lb2u*lt1u*mmgl*mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst1)*
     pow2(mmst1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*
     ltu*mmgl*mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1)*pow4(g3))/(mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*lb2u*ltu*mmgl*mmsb2*mmt*DeltaInv(
     mmt,mmsb2,mmst1)*pow2(mmst1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl
     + mmsb2)) + (8*lt1u*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*
     pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*lgu*lt1u*mmgl*
     mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow4(g3))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb2)) - (8*ltu*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmst1)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8
     *lgu*ltu*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow4(g3))/(3.
     *mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (4*mmgl*mmsb2*mmt*DeltaInv(mmt
     ,mmsb2,mmst1)*pow2(lt1u)*pow2(mmst1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) - (4*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     lt1u)*pow2(mmst1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) +
     (4*mmgl*mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(ltu)*pow2(mmst1)*pow4(g3)
     )/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (4*mmgl*mmsb2*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow2(mmst1)*pow4(g3))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb2)) + (8*lt2u*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmst2)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (8*lgu*lt2u*
     mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(3.*mgl*(-mmgl +
     mmsb1)*(mmst1 - mmst2)) - (8*ltu*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmst2)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (8*lgu*ltu*mmgl*
     mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)
     *(mmst1 - mmst2)) + (8*lt2u*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*
     pow4(g3))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (8*lgu*lt2u*mmgl*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(
     mmst1 - mmst2)) - (8*ltu*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*
     pow4(g3))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (8*lgu*ltu*mmgl*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(
     mmst1 - mmst2)) - (4*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lt2u)*pow2(
     mmst2)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (4*mmgl*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(lt2u)*pow2(mmst2)*pow4(g3))/(3.*mgl*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) + (4*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*
     pow2(mmst2)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (4*mmgl*
     mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow2(mmst2)*pow4(g3))/(3.*mgl*(-
     mmgl + mmsb2)*(mmst1 - mmst2)) - (8*lt2u*mmgl*mmsb1*mmt*DeltaInv(mmt,mmsb1
     ,mmst2)*pow2(mmst2)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) +
     (8*lb1u*lt2u*mmgl*mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2)*pow4(g3)
     )/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*ltu*mmgl*mmsb1*mmt*
     DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(
     -mmgl + mmsb1)) - (8*lb1u*ltu*mmgl*mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst2)*
     pow2(mmst2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*
     lt2u*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*lgu*lt2u*mmgl*mmsb1*mmt*DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) + (8*ltu*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(
     g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*lgu*ltu*mmgl*mmsb1*mmt
     *DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb1)) + (4*mmgl*mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(
     lt2u)*pow2(mmst2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) +
     (4*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lt2u)*pow2(mmst2)*pow4(g3)
     )/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (4*mmgl*mmsb1*mmt*
     DeltaInv(mmt,mmsb1,mmst2)*pow2(ltu)*pow2(mmst2)*pow4(g3))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb1)) - (4*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*
     pow2(ltu)*pow2(mmst2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1
     )) - (8*lt2u*mmgl*mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2)*pow4(g3)
     )/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*lb2u*lt2u*mmgl*mmsb2*mmt*
     DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) + (8*ltu*mmgl*mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst2)*
     pow2(mmst2)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*lb2u*
     ltu*mmgl*mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2)*pow4(g3))/(3.*mgl
     *(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*lt2u*mmgl*mmsb2*mmt*DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb2)) + (8*lgu*lt2u*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*
     pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*ltu*mmgl*mmsb2
     *mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) - (8*lgu*ltu*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmst2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (4*
     mmgl*mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(lt2u)*pow2(mmst2)*pow4(g3))/
     (3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (4*mmgl*mmsb2*mmt*DeltaInv(
     mmt,mmst2,mmgl)*pow2(lt2u)*pow2(mmst2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) - (4*mmgl*mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(
     ltu)*pow2(mmst2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) -
     (4*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow2(mmst2)*pow4(g3))
     /(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (112*mmgl*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)) + (16*lgu*mmgl*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(mgl*(mmst1 - mmst2)) + (16*
     ltu*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(mgl*(mmst1 - mmst2)
     ) - (16*lgu*ltu*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*
     (mmst1 - mmst2)) + (112*mmgl*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4
     (g3))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (16*lgu*mmgl*mmsb1*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmst1 -
     mmst2)) - (16*ltu*mmgl*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/
     (mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (16*lgu*ltu*mmgl*mmsb1*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2
     )) + (112*mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*
     (-mmgl + mmsb2)*(mmst1 - mmst2)) - (16*lgu*mmgl*mmsb2*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmt)*pow4(g3))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (16*ltu*
     mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(mgl*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) + (16*lgu*ltu*mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (112*mmgl*
     mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)
     *(mmst1 - mmst2)) - (8*lt1u*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*
     pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (8*lgu*lt1u*mmgl*mmst1*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(
     mmst1 - mmst2)) - (24*ltu*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*
     pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (8*lgu*ltu*mmgl*mmst1*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(
     mmst1 - mmst2)) + (16*lt1u*ltu*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmt)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (112*mmgl*mmst1*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(
     mmst1 - mmst2)) - (8*lt1u*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*
     pow4(g3))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (8*lgu*lt1u*mmgl*mmst1*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(
     mmst1 - mmst2)) - (24*ltu*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*
     pow4(g3))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (8*lgu*ltu*mmgl*mmst1*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(
     mmst1 - mmst2)) + (16*lt1u*ltu*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmt)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (16*mmgl*zt2*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)) + (
     16*mmgl*mmsb1*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(-
     mmgl + mmsb1)*(mmst1 - mmst2)) + (16*mmgl*mmsb2*zt2*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (16*
     mmgl*mmst1*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl
      + mmsb1)*(mmst1 - mmst2)) + (16*mmgl*mmst1*zt2*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (112*mmgl*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)) - (
     16*lgu*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/(mgl*(mmst1 -
     mmst2)) - (16*ltu*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/(mgl*(
     mmst1 - mmst2)) + (16*lgu*ltu*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4
     (g3))/(3.*mgl*(mmst1 - mmst2)) - (112*mmgl*mmsb1*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (16*lgu*
     mmgl*mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/(mgl*(-mmgl +
     mmsb1)*(mmst1 - mmst2)) + (16*ltu*mmgl*mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2
     (mmt)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (16*lgu*ltu*mmgl*
     mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)
     *(mmst1 - mmst2)) - (112*mmgl*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*
     pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (16*lgu*mmgl*mmsb2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/(mgl*(-mmgl + mmsb2)*(mmst1 -
     mmst2)) + (16*ltu*mmgl*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/
     (mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (16*lgu*ltu*mmgl*mmsb2*DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2
     )) - (112*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*
     (-mmgl + mmsb1)*(mmst1 - mmst2)) + (8*lt2u*mmgl*mmst2*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmt)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (8*lgu*
     lt2u*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(-
     mmgl + mmsb1)*(mmst1 - mmst2)) + (24*ltu*mmgl*mmst2*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmt)*pow4(g3))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (8*lgu*
     ltu*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl
      + mmsb1)*(mmst1 - mmst2)) - (16*lt2u*ltu*mmgl*mmst2*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (112*
     mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) + (8*lt2u*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2
     (mmt)*pow4(g3))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (8*lgu*lt2u*mmgl*
     mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)
     *(mmst1 - mmst2)) + (24*ltu*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*
     pow4(g3))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (8*lgu*ltu*mmgl*mmst2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(
     mmst1 - mmst2)) - (16*lt2u*ltu*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmt)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (16*mmgl*zt2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)) - (
     16*mmgl*mmsb1*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(-
     mmgl + mmsb1)*(mmst1 - mmst2)) - (16*mmgl*mmsb2*zt2*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (16*
     mmgl*mmst2*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl
      + mmsb1)*(mmst1 - mmst2)) - (16*mmgl*mmst2*zt2*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (8*mmgl*
     DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu)*pow2(mmt)*pow4(g3))/(3.*
     mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (8*mmgl*DeltaInv(mmt,mmst1,mmgl)*
     Fin3(mmt,mmst1,mmgl,mmu)*pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(
     mmst1 - mmst2)) - (8*mmgl*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu
     )*pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (8*mmgl*
     DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*pow2(mmt)*pow4(g3))/(3.*
     mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (8*mmgl*DeltaInv(mmt,mmst1,mmgl)*
     pow2(lgu)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)) + (8*mmgl*mmsb1*
     DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl +
     mmsb1)*(mmst1 - mmst2)) + (8*mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)
     *pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (8*mmgl*
     DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 -
     mmst2)) - (8*mmgl*mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(mmt)*pow4(
     g3))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (8*mmgl*mmsb2*DeltaInv(mmt
     ,mmst2,mmgl)*pow2(lgu)*pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(mmst1
     - mmst2)) + (4*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(lt1u)*pow2(mmt)*
     pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (4*mmgl*mmst1*
     DeltaInv(mmt,mmst1,mmgl)*pow2(lt1u)*pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) - (4*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(lt2u
     )*pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (4*mmgl*
     mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(lt2u)*pow2(mmt)*pow4(g3))/(3.*mgl*(-
     mmgl + mmsb2)*(mmst1 - mmst2)) - (8*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu
     )*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)) + (8*mmgl*mmsb1*DeltaInv(
     mmt,mmst1,mmgl)*pow2(ltu)*pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl + mmsb1)*(
     mmst1 - mmst2)) + (8*mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow2(
     mmt)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (4*mmgl*mmst1*
     DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow2(mmt)*pow4(g3))/(mgl*(-mmgl + mmsb1
     )*(mmst1 - mmst2)) + (4*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow2
     (mmt)*pow4(g3))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (8*mmgl*DeltaInv(
     mmt,mmst2,mmgl)*pow2(ltu)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)) - (
     8*mmgl*mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow2(mmt)*pow4(g3))/(3.*
     mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (8*mmgl*mmsb2*DeltaInv(mmt,mmst2,
     mmgl)*pow2(ltu)*pow2(mmt)*pow4(g3))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2
     )) - (4*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow2(mmt)*pow4(g3))/
     (mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (4*mmgl*mmst2*DeltaInv(mmt,mmst2,
     mmgl)*pow2(ltu)*pow2(mmt)*pow4(g3))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2))
     - (112*mmgl*mmsb1*mmst1*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmt)*pow4(g3))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*lt1u*mmgl*mmsb1*mmst1*
     DeltaInv(mmt,mmsb1,mmst1)*pow2(mmt)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-
     mmgl + mmsb1)) + (8*lb1u*lt1u*mmgl*mmsb1*mmst1*DeltaInv(mmt,mmsb1,mmst1)*
     pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (24*ltu
     *mmgl*mmsb1*mmst1*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmt)*pow4(g3))/(mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*lb1u*ltu*mmgl*mmsb1*mmst1*
     DeltaInv(mmt,mmsb1,mmst1)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2
     (-mmgl + mmsb1)) - (16*lt1u*ltu*mmgl*mmsb1*mmst1*DeltaInv(mmt,mmsb1,mmst1)
     *pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (16*
     mmgl*mmsb1*mmst1*zt2*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmt)*pow4(g3))/(3.*mgl
     *(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (112*mmgl*mmsb1*mmst2*DeltaInv(mmt
     ,mmsb1,mmst2)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) - (8*lt2u*mmgl*mmsb1*mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmt)*
     pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*lb1u*lt2u*mmgl*
     mmsb1*mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb1)) - (24*ltu*mmgl*mmsb1*mmst2*DeltaInv(mmt,mmsb1,
     mmst2)*pow2(mmt)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*
     lb1u*ltu*mmgl*mmsb1*mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmt)*pow4(g3))/(
     3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (16*lt2u*ltu*mmgl*mmsb1*
     mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2
     )*pow2(-mmgl + mmsb1)) + (16*mmgl*mmsb1*mmst2*zt2*DeltaInv(mmt,mmsb1,mmst2
     )*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (112*
     mmgl*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*lt1u*mmgl*mmsb1*mmst1*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) + (8*lgu*lt1u*mmgl*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*
     pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (24*ltu*mmgl*
     mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb1)) - (8*lgu*ltu*mmgl*mmsb1*mmst1*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1
     )) - (16*lt1u*ltu*mmgl*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4
     (g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (16*mmgl*mmsb1*mmst1*
     zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb1)) + (112*mmgl*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2
     (mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*lt2u*
     mmgl*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/(mgl*(mmst1
     - mmst2)*pow2(-mmgl + mmsb1)) - (8*lgu*lt2u*mmgl*mmsb1*mmst2*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1
     )) - (24*ltu*mmgl*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))
     /(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*lgu*ltu*mmgl*mmsb1*mmst2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(
     -mmgl + mmsb1)) + (16*lt2u*ltu*mmgl*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (16*
     mmgl*mmsb1*mmst2*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*
     (mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmsb1*DeltaInv(mmt,mmsb1,
     mmst1)*Fin3(mmt,mmsb1,mmst1,mmu)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmsb1*DeltaInv(mmt,mmsb1,mmst2)*Fin3
     (mmt,mmsb1,mmst2,mmu)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-
     mmgl + mmsb1)) - (8*mmgl*mmsb1*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,
     mmgl,mmu)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1))
     + (8*mmgl*mmsb1*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*pow2(mmt
     )*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (4*mmgl*mmsb1*
     mmst1*DeltaInv(mmt,mmsb1,mmst1)*pow2(lt1u)*pow2(mmt)*pow4(g3))/(3.*mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (4*mmgl*mmsb1*mmst1*DeltaInv(mmt,
     mmst1,mmgl)*pow2(lt1u)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-
     mmgl + mmsb1)) + (4*mmgl*mmsb1*mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow2(lt2u)*
     pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (4*mmgl
     *mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(lt2u)*pow2(mmt)*pow4(g3))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (4*mmgl*mmsb1*mmst1*DeltaInv(
     mmt,mmsb1,mmst1)*pow2(ltu)*pow2(mmt)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-
     mmgl + mmsb1)) + (4*mmgl*mmsb1*mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow2(ltu)*
     pow2(mmt)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (4*mmgl*
     mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow2(mmt)*pow4(g3))/(mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (4*mmgl*mmsb1*mmst2*DeltaInv(mmt,
     mmst2,mmgl)*pow2(ltu)*pow2(mmt)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl
     + mmsb1)) - (56*mmgl*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow2(mmt)*pow4(
     g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*lb1u*mmgl*DeltaInv(
     mmt,mmsb1,mmst1)*pow2(mmsb1)*pow2(mmt)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2
     (-mmgl + mmsb1)) + (8*ltu*mmgl*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow2(
     mmt)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*lb1u*ltu*
     mmgl*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow2(mmt)*pow4(g3))/(3.*mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*zt2*DeltaInv(mmt,mmsb1,mmst1
     )*pow2(mmsb1)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) + (56*mmgl*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*pow2(mmt)*pow4(g3
     ))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*lb1u*mmgl*DeltaInv(
     mmt,mmsb1,mmst2)*pow2(mmsb1)*pow2(mmt)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2
     (-mmgl + mmsb1)) - (8*ltu*mmgl*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*pow2(
     mmt)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*lb1u*ltu*
     mmgl*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*pow2(mmt)*pow4(g3))/(3.*mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*zt2*DeltaInv(mmt,mmsb1,mmst2
     )*pow2(mmsb1)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) - (56*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(mmt)*pow4(g3)
     )/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*lgu*mmgl*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmsb1)*pow2(mmt)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-
     mmgl + mmsb1)) + (8*ltu*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(mmt
     )*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*lgu*ltu*mmgl*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmsb1)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) +
     (56*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(mmt)*pow4(g3))/(3.*mgl*
     (mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*lgu*mmgl*DeltaInv(mmt,mmst2,mmgl
     )*pow2(mmsb1)*pow2(mmt)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)
     ) - (8*ltu*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(mmt)*pow4(g3))/(
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*lgu*ltu*mmgl*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmsb1)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-
     mmgl + mmsb1)) + (8*mmgl*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(mmt
     )*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (4*mmgl*
     DeltaInv(mmt,mmsb1,mmst1)*pow2(lb1u)*pow2(mmsb1)*pow2(mmt)*pow4(g3))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (4*mmgl*DeltaInv(mmt,mmsb1,
     mmst2)*pow2(lb1u)*pow2(mmsb1)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb1)) - (4*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow2(
     mmsb1)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) +
     (4*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(mmsb1)*pow2(mmt)*pow4(g3))
     /(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (4*mmgl*DeltaInv(mmt,mmsb1
     ,mmst1)*pow2(ltu)*pow2(mmsb1)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb1)) + (4*mmgl*DeltaInv(mmt,mmsb1,mmst2)*pow2(ltu)*pow2(
     mmsb1)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) -
     (4*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow2(mmsb1)*pow2(mmt)*pow4(g3))
     /(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (4*mmgl*DeltaInv(mmt,mmst2
     ,mmgl)*pow2(ltu)*pow2(mmsb1)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb1)) - (112*mmgl*mmsb2*mmst1*DeltaInv(mmt,mmsb2,mmst1)*
     pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*lt1u
     *mmgl*mmsb2*mmst1*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmt)*pow4(g3))/(mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*lb2u*lt1u*mmgl*mmsb2*mmst1*
     DeltaInv(mmt,mmsb2,mmst1)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2
     (-mmgl + mmsb2)) + (24*ltu*mmgl*mmsb2*mmst1*DeltaInv(mmt,mmsb2,mmst1)*pow2
     (mmt)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*lb2u*ltu*
     mmgl*mmsb2*mmst1*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmt)*pow4(g3))/(3.*mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (16*lt1u*ltu*mmgl*mmsb2*mmst1*
     DeltaInv(mmt,mmsb2,mmst1)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2
     (-mmgl + mmsb2)) - (16*mmgl*mmsb2*mmst1*zt2*DeltaInv(mmt,mmsb2,mmst1)*pow2
     (mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (112*mmgl*
     mmsb2*mmst2*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb2)) - (8*lt2u*mmgl*mmsb2*mmst2*DeltaInv(mmt,mmsb2,
     mmst2)*pow2(mmt)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*
     lb2u*lt2u*mmgl*mmsb2*mmst2*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmt)*pow4(g3))/(
     3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (24*ltu*mmgl*mmsb2*mmst2*
     DeltaInv(mmt,mmsb2,mmst2)*pow2(mmt)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-
     mmgl + mmsb2)) + (8*lb2u*ltu*mmgl*mmsb2*mmst2*DeltaInv(mmt,mmsb2,mmst2)*
     pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (16*
     lt2u*ltu*mmgl*mmsb2*mmst2*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmt)*pow4(g3))/(
     3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (16*mmgl*mmsb2*mmst2*zt2*
     DeltaInv(mmt,mmsb2,mmst2)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2
     (-mmgl + mmsb2)) - (112*mmgl*mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt
     )*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*lt1u*mmgl*
     mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb2)) + (8*lgu*lt1u*mmgl*mmsb2*mmst1*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2
     )) + (24*ltu*mmgl*mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))
     /(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*lgu*ltu*mmgl*mmsb2*mmst1*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(
     -mmgl + mmsb2)) - (16*lt1u*ltu*mmgl*mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (16*
     mmgl*mmsb2*mmst1*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*
     (mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (112*mmgl*mmsb2*mmst2*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2
     )) - (8*lt2u*mmgl*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))
     /(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*lgu*lt2u*mmgl*mmsb2*mmst2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(
     -mmgl + mmsb2)) - (24*ltu*mmgl*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmt)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*lgu*ltu*mmgl
     *mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb2)) + (16*lt2u*ltu*mmgl*mmsb2*mmst2*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2
     )) + (16*mmgl*mmsb2*mmst2*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow4(g3))
     /(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmsb2*DeltaInv(mmt
     ,mmsb2,mmst1)*Fin3(mmt,mmsb2,mmst1,mmu)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1
      - mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmsb2*DeltaInv(mmt,mmsb2,mmst2)*
     Fin3(mmt,mmsb2,mmst2,mmu)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2
     (-mmgl + mmsb2)) - (8*mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,
     mmgl,mmu)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2))
     + (8*mmgl*mmsb2*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*pow2(mmt
     )*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (4*mmgl*mmsb2*
     mmst1*DeltaInv(mmt,mmsb2,mmst1)*pow2(lt1u)*pow2(mmt)*pow4(g3))/(3.*mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (4*mmgl*mmsb2*mmst1*DeltaInv(mmt,
     mmst1,mmgl)*pow2(lt1u)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-
     mmgl + mmsb2)) + (4*mmgl*mmsb2*mmst2*DeltaInv(mmt,mmsb2,mmst2)*pow2(lt2u)*
     pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (4*mmgl
     *mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(lt2u)*pow2(mmt)*pow4(g3))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (4*mmgl*mmsb2*mmst1*DeltaInv(
     mmt,mmsb2,mmst1)*pow2(ltu)*pow2(mmt)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-
     mmgl + mmsb2)) + (4*mmgl*mmsb2*mmst2*DeltaInv(mmt,mmsb2,mmst2)*pow2(ltu)*
     pow2(mmt)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (4*mmgl*
     mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow2(mmt)*pow4(g3))/(mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (4*mmgl*mmsb2*mmst2*DeltaInv(mmt,
     mmst2,mmgl)*pow2(ltu)*pow2(mmt)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl
     + mmsb2)) - (56*mmgl*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2)*pow2(mmt)*pow4(
     g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*lb2u*mmgl*DeltaInv(
     mmt,mmsb2,mmst1)*pow2(mmsb2)*pow2(mmt)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2
     (-mmgl + mmsb2)) + (8*ltu*mmgl*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2)*pow2(
     mmt)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*lb2u*ltu*
     mmgl*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2)*pow2(mmt)*pow4(g3))/(3.*mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*zt2*DeltaInv(mmt,mmsb2,mmst1
     )*pow2(mmsb2)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb2)) + (56*mmgl*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow2(mmt)*pow4(g3
     ))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*lb2u*mmgl*DeltaInv(
     mmt,mmsb2,mmst2)*pow2(mmsb2)*pow2(mmt)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2
     (-mmgl + mmsb2)) - (8*ltu*mmgl*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow2(
     mmt)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*lb2u*ltu*
     mmgl*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow2(mmt)*pow4(g3))/(3.*mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*zt2*DeltaInv(mmt,mmsb2,mmst2
     )*pow2(mmsb2)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb2)) - (56*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(mmt)*pow4(g3)
     )/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*lgu*mmgl*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmsb2)*pow2(mmt)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-
     mmgl + mmsb2)) + (8*ltu*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(mmt
     )*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*lgu*ltu*mmgl*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmsb2)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) +
     (56*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(mmt)*pow4(g3))/(3.*mgl*
     (mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*lgu*mmgl*DeltaInv(mmt,mmst2,mmgl
     )*pow2(mmsb2)*pow2(mmt)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)
     ) - (8*ltu*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(mmt)*pow4(g3))/(
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*lgu*ltu*mmgl*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmsb2)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-
     mmgl + mmsb2)) + (8*mmgl*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(mmt
     )*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (4*mmgl*
     DeltaInv(mmt,mmsb2,mmst1)*pow2(lb2u)*pow2(mmsb2)*pow2(mmt)*pow4(g3))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (4*mmgl*DeltaInv(mmt,mmsb2,
     mmst2)*pow2(lb2u)*pow2(mmsb2)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) - (4*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow2(
     mmsb2)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) +
     (4*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(mmsb2)*pow2(mmt)*pow4(g3))
     /(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (4*mmgl*DeltaInv(mmt,mmsb2
     ,mmst1)*pow2(ltu)*pow2(mmsb2)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) + (4*mmgl*DeltaInv(mmt,mmsb2,mmst2)*pow2(ltu)*pow2(
     mmsb2)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) -
     (4*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow2(mmsb2)*pow2(mmt)*pow4(g3))
     /(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (4*mmgl*DeltaInv(mmt,mmst2
     ,mmgl)*pow2(ltu)*pow2(mmsb2)*pow2(mmt)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) + (56*mmgl*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1)*
     pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*lb1u*mmgl*mmt*
     DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(
     -mmgl + mmsb1)) - (8*ltu*mmgl*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1)*
     pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*lb1u*ltu*mmgl*mmt
     *DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb1)) + (8*mmgl*mmt*zt2*DeltaInv(mmt,mmsb1,mmst1)*pow3(
     mmsb1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (56*mmgl*
     mmt*DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2
     )*pow2(-mmgl + mmsb1)) + (8*lb1u*mmgl*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow3(
     mmsb1)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*ltu*mmgl*
     mmt*DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1)*pow4(g3))/(mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb1)) - (8*lb1u*ltu*mmgl*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow3
     (mmsb1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*
     mmt*zt2*DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1)*pow4(g3))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb1)) + (56*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow3(
     mmsb1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*lgu*
     mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1)*pow4(g3))/(mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb1)) - (8*ltu*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow3(mmsb1)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*lgu*
     ltu*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1)*pow4(g3))/(3.*mgl*(mmst1
      - mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*
     pow3(mmsb1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (56*
     mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1)*pow4(g3))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb1)) + (8*lgu*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*
     pow3(mmsb1)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*ltu*
     mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1)*pow4(g3))/(mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb1)) - (8*lgu*ltu*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)
     *pow3(mmsb1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*
     mmgl*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1)*pow4(g3))/(3.*mgl*(mmst1
      - mmst2)*pow2(-mmgl + mmsb1)) + (4*mmgl*mmt*DeltaInv(mmt,mmsb1,mmst1)*
     pow2(lb1u)*pow3(mmsb1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) - (4*mmgl*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(lb1u)*pow3(mmsb1)*
     pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (4*mmgl*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow3(mmsb1)*pow4(g3))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb1)) - (4*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(
     lgu)*pow3(mmsb1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) +
     (4*mmgl*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(ltu)*pow3(mmsb1)*pow4(g3))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (4*mmgl*mmt*DeltaInv(mmt,mmsb1,
     mmst2)*pow2(ltu)*pow3(mmsb1)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl
     + mmsb1)) + (4*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow3(mmsb1)*
     pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (4*mmgl*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(ltu)*pow3(mmsb1)*pow4(g3))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb1)) - (16*mmgl*mmsb1*mmt*Fin3(mmt,mmsb1,mmst1,mmu)
     *pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow3(-mmgl + mmsb1)) + (16*mmgl*mmsb1*
     mmt*Fin3(mmt,mmsb1,mmst2,mmu)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow3(-mmgl
      + mmsb1)) + (16*mmgl*mmsb1*mmt*Fin3(mmt,mmst1,mmgl,mmu)*pow4(g3))/(3.*mgl
     *(mmst1 - mmst2)*pow3(-mmgl + mmsb1)) - (16*mmgl*mmsb1*mmt*Fin3(mmt,mmst2,
     mmgl,mmu)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow3(-mmgl + mmsb1)) + (56*
     mmgl*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmsb2)*pow4(g3))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb2)) - (8*lb2u*mmgl*mmt*DeltaInv(mmt,mmsb2,mmst1)*
     pow3(mmsb2)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*ltu*
     mmgl*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmsb2)*pow4(g3))/(mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb2)) + (8*lb2u*ltu*mmgl*mmt*DeltaInv(mmt,mmsb2,
     mmst1)*pow3(mmsb2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2))
     + (8*mmgl*mmt*zt2*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmsb2)*pow4(g3))/(3.*mgl*
     (mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (56*mmgl*mmt*DeltaInv(mmt,mmsb2,
     mmst2)*pow3(mmsb2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2))
     + (8*lb2u*mmgl*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmsb2)*pow4(g3))/(mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*ltu*mmgl*mmt*DeltaInv(mmt,mmsb2,
     mmst2)*pow3(mmsb2)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (
     8*lb2u*ltu*mmgl*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmsb2)*pow4(g3))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmt*zt2*DeltaInv(mmt,
     mmsb2,mmst2)*pow3(mmsb2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb2)) + (56*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2)*pow4(g3))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*lgu*mmgl*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow3(mmsb2)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)
     ) - (8*ltu*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2)*pow4(g3))/(mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*lgu*ltu*mmgl*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow3(mmsb2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb2)) + (8*mmgl*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2)*pow4(g3))/(
     3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (56*mmgl*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow3(mmsb2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb2)) + (8*lgu*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2)*pow4(g3))/(
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*ltu*mmgl*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow3(mmsb2)*pow4(g3))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)
     ) - (8*lgu*ltu*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2)*pow4(g3))/(3.
     *mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmt*zt2*DeltaInv(mmt,
     mmst2,mmgl)*pow3(mmsb2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb2)) + (4*mmgl*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(lb2u)*pow3(mmsb2)*
     pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (4*mmgl*mmt*
     DeltaInv(mmt,mmsb2,mmst2)*pow2(lb2u)*pow3(mmsb2)*pow4(g3))/(3.*mgl*(mmst1
     - mmst2)*pow2(-mmgl + mmsb2)) + (4*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     lgu)*pow3(mmsb2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) -
     (4*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow3(mmsb2)*pow4(g3))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (4*mmgl*mmt*DeltaInv(mmt,mmsb2,
     mmst1)*pow2(ltu)*pow3(mmsb2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl
     + mmsb2)) - (4*mmgl*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(ltu)*pow3(mmsb2)*
     pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (4*mmgl*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(ltu)*pow3(mmsb2)*pow4(g3))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb2)) - (4*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(
     ltu)*pow3(mmsb2)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) -
     (16*mmgl*mmsb2*mmt*Fin3(mmt,mmsb2,mmst1,mmu)*pow4(g3))/(3.*mgl*(mmst1 -
     mmst2)*pow3(-mmgl + mmsb2)) + (16*mmgl*mmsb2*mmt*Fin3(mmt,mmsb2,mmst2,mmu)
     *pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow3(-mmgl + mmsb2)) + (16*mmgl*mmsb2*
     mmt*Fin3(mmt,mmst1,mmgl,mmu)*pow4(g3))/(3.*mgl*(mmst1 - mmst2)*pow3(-mmgl
     + mmsb2)) - (16*mmgl*mmsb2*mmt*Fin3(mmt,mmst2,mmgl,mmu)*pow4(g3))/(3.*mgl*
     (mmst1 - mmst2)*pow3(-mmgl + mmsb2)))/pow4(g3)) + (pow3(Xb)*((-128*mmb*
     mmgl*mmsb1*pow4(g3))/(9.*mgl*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) + (256*
     lb1u*mmb*mmgl*mmsb1*pow4(g3))/(9.*mgl*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2))
     - (128*lgu*mmb*mmgl*mmsb1*pow4(g3))/(9.*mgl*(-mmgl + mmsb1)*pow3(mmsb1 -
     mmsb2)) + (128*lb1u*lgu*mmb*mmgl*mmsb1*pow4(g3))/(9.*mgl*(-mmgl + mmsb1)*
     pow3(mmsb1 - mmsb2)) + (128*mmb*mmgl*mmsb2*pow4(g3))/(9.*mgl*(-mmgl +
     mmsb1)*pow3(mmsb1 - mmsb2)) + (128*lb1u*mmb*mmgl*mmsb2*pow4(g3))/(9.*mgl*(
     -mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) - (128*lb2u*mmb*mmgl*mmsb2*pow4(g3))/(
     9.*mgl*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) + (128*lb1u*lb2u*mmb*mmgl*
     mmsb2*pow4(g3))/(9.*mgl*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) - (128*lgu*
     mmb*mmgl*mmsb2*pow4(g3))/(9.*mgl*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) + (
     256*lb1u*lgu*mmb*mmgl*mmsb2*pow4(g3))/(9.*mgl*(-mmgl + mmsb1)*pow3(mmsb1 -
     mmsb2)) - (128*lb2u*lgu*mmb*mmgl*mmsb2*pow4(g3))/(9.*mgl*(-mmgl + mmsb1)*
     pow3(mmsb1 - mmsb2)) - (128*mmb*mmgl*mmsb1*pow4(g3))/(9.*mgl*(-mmgl +
     mmsb2)*pow3(mmsb1 - mmsb2)) + (128*lb1u*mmb*mmgl*mmsb1*pow4(g3))/(9.*mgl*(
     -mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) - (128*lb2u*mmb*mmgl*mmsb1*pow4(g3))/(
     9.*mgl*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) + (128*lb1u*lb2u*mmb*mmgl*
     mmsb1*pow4(g3))/(9.*mgl*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) + (128*lgu*
     mmb*mmgl*mmsb1*pow4(g3))/(9.*mgl*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) - (
     128*lb1u*lgu*mmb*mmgl*mmsb1*pow4(g3))/(9.*mgl*(-mmgl + mmsb2)*pow3(mmsb1 -
     mmsb2)) + (128*mmb*mmgl*mmsb2*pow4(g3))/(9.*mgl*(-mmgl + mmsb2)*pow3(mmsb1
      - mmsb2)) - (256*lb2u*mmb*mmgl*mmsb2*pow4(g3))/(9.*mgl*(-mmgl + mmsb2)*
     pow3(mmsb1 - mmsb2)) + (256*lb1u*lb2u*mmb*mmgl*mmsb2*pow4(g3))/(9.*mgl*(-
     mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) + (128*lgu*mmb*mmgl*mmsb2*pow4(g3))/(9.
     *mgl*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) - (256*lb1u*lgu*mmb*mmgl*mmsb2*
     pow4(g3))/(9.*mgl*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) + (128*lb2u*lgu*mmb
     *mmgl*mmsb2*pow4(g3))/(9.*mgl*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) - (128*
     mmb*mmgl*mmsb1*pow2(lb1u)*pow4(g3))/(9.*mgl*(-mmgl + mmsb1)*pow3(mmsb1 -
     mmsb2)) - (256*mmb*mmgl*mmsb2*pow2(lb1u)*pow4(g3))/(9.*mgl*(-mmgl + mmsb1)
     *pow3(mmsb1 - mmsb2)) - (128*mmb*mmgl*mmsb2*pow2(lb2u)*pow4(g3))/(9.*mgl*(
     -mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) - (128*lb1u*mmb*mmgl*mmsb1*mmsb2*pow4(
     g3))/(9.*mgl*pow2(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) + (128*lb1u*lb2u*mmb
     *mmgl*mmsb1*mmsb2*pow4(g3))/(9.*mgl*pow2(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2
     )) + (128*lgu*mmb*mmgl*mmsb1*mmsb2*pow4(g3))/(9.*mgl*pow2(-mmgl + mmsb1)*
     pow3(mmsb1 - mmsb2)) - (128*lb2u*lgu*mmb*mmgl*mmsb1*mmsb2*pow4(g3))/(9.*
     mgl*pow2(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) + (128*lb1u*mmb*mmgl*pow2(
     mmsb1)*pow4(g3))/(9.*mgl*pow2(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) - (128*
     lgu*mmb*mmgl*pow2(mmsb1)*pow4(g3))/(9.*mgl*pow2(-mmgl + mmsb1)*pow3(mmsb1
     - mmsb2)) + (128*lb1u*lgu*mmb*mmgl*pow2(mmsb1)*pow4(g3))/(9.*mgl*pow2(-
     mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) - (128*mmb*mmgl*pow2(lb1u)*pow2(mmsb1)*
     pow4(g3))/(9.*mgl*pow2(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) + (128*lb2u*mmb
     *mmgl*mmsb1*mmsb2*pow4(g3))/(9.*mgl*pow2(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2
     )) - (128*lb1u*lb2u*mmb*mmgl*mmsb1*mmsb2*pow4(g3))/(9.*mgl*pow2(-mmgl +
     mmsb2)*pow3(mmsb1 - mmsb2)) - (128*lgu*mmb*mmgl*mmsb1*mmsb2*pow4(g3))/(9.*
     mgl*pow2(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) + (128*lb1u*lgu*mmb*mmgl*
     mmsb1*mmsb2*pow4(g3))/(9.*mgl*pow2(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) - (
     128*lb2u*mmb*mmgl*pow2(mmsb2)*pow4(g3))/(9.*mgl*pow2(-mmgl + mmsb2)*pow3(
     mmsb1 - mmsb2)) + (128*lgu*mmb*mmgl*pow2(mmsb2)*pow4(g3))/(9.*mgl*pow2(-
     mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) - (128*lb2u*lgu*mmb*mmgl*pow2(mmsb2)*
     pow4(g3))/(9.*mgl*pow2(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) + (128*mmb*mmgl
     *pow2(lb2u)*pow2(mmsb2)*pow4(g3))/(9.*mgl*pow2(-mmgl + mmsb2)*pow3(mmsb1 -
     mmsb2)) + (256*lb1u*lb2u*mmb*mmgl*pow2(mmsb2)*pow4(g3))/(9.*mgl*(-mmgl +
     mmsb1)*pow4(mmsb1 - mmsb2)) + (256*lb1u*lgu*mmb*mmgl*pow2(mmsb2)*pow4(g3))
     /(9.*mgl*(-mmgl + mmsb1)*pow4(mmsb1 - mmsb2)) - (256*lb2u*lgu*mmb*mmgl*
     pow2(mmsb2)*pow4(g3))/(9.*mgl*(-mmgl + mmsb1)*pow4(mmsb1 - mmsb2)) + (256*
     lb1u*lb2u*mmb*mmgl*pow2(mmsb2)*pow4(g3))/(9.*mgl*(-mmgl + mmsb2)*pow4(
     mmsb1 - mmsb2)) - (256*lb1u*lgu*mmb*mmgl*pow2(mmsb2)*pow4(g3))/(9.*mgl*(-
     mmgl + mmsb2)*pow4(mmsb1 - mmsb2)) + (256*lb2u*lgu*mmb*mmgl*pow2(mmsb2)*
     pow4(g3))/(9.*mgl*(-mmgl + mmsb2)*pow4(mmsb1 - mmsb2)) - (256*mmb*mmgl*
     pow2(lb1u)*pow2(mmsb2)*pow4(g3))/(9.*mgl*(-mmgl + mmsb1)*pow4(mmsb1 -
     mmsb2)) - (256*mmb*mmgl*pow2(lb2u)*pow2(mmsb2)*pow4(g3))/(9.*mgl*(-mmgl +
     mmsb2)*pow4(mmsb1 - mmsb2))))/pow4(g3) + ((-3659*pow4(g3))/54. + (25*lb1u*
     pow4(g3))/3. - (53*lb2u*pow4(g3))/9. + (1016*lgu*pow4(g3))/9. + (64*lb1u*
     lgu*pow4(g3))/9. + (64*lb2u*lgu*pow4(g3))/9. - (8*mmsb1*pow4(g3))/(-mmgl +
     mmsb1) + (118*lb1u*mmsb1*pow4(g3))/(9.*(-mmgl + mmsb1)) - (58*lgu*mmsb1*
     pow4(g3))/(3.*(-mmgl + mmsb1)) - (55*lb1u*lgu*mmsb1*pow4(g3))/(9.*(-mmgl +
     mmsb1)) + (890*mmsb2*pow4(g3))/(9.*(-mmgl + mmsb1)) - (44*lb1u*mmsb2*pow4(
     g3))/(3.*(-mmgl + mmsb1)) - (44*lb2u*mmsb2*pow4(g3))/(9.*(-mmgl + mmsb1))
     - (4*lb1u*lb2u*mmsb2*pow4(g3))/(9.*(-mmgl + mmsb1)) - (598*lgu*mmsb2*pow4(
     g3))/(9.*(-mmgl + mmsb1)) + (4*lb1u*lgu*mmsb2*pow4(g3))/(9.*(-mmgl + mmsb1
     )) + (13*lb2u*lgu*mmsb2*pow4(g3))/(9.*(-mmgl + mmsb1)) + (128*lb1u*mmsb2*
     pow4(g3))/(9.*(mmsb1 - mmsb2)) - (128*lb2u*mmsb2*pow4(g3))/(9.*(mmsb1 -
     mmsb2)) + (6*mmsb1*pow4(g3))/(-mmgl + mmsb2) - (4*lb1u*mmsb1*pow4(g3))/(-
     mmgl + mmsb2) - (2*lgu*mmsb1*pow4(g3))/(-mmgl + mmsb2) + (lb1u*lgu*mmsb1*
     pow4(g3))/(-mmgl + mmsb2) - (908*mmsb2*pow4(g3))/(9.*(-mmgl + mmsb2)) + (4
     *lb1u*mmsb2*pow4(g3))/(9.*(-mmgl + mmsb2)) + (254*lb2u*mmsb2*pow4(g3))/(9.
     *(-mmgl + mmsb2)) + (4*lb1u*lb2u*mmsb2*pow4(g3))/(9.*(-mmgl + mmsb2)) + (
     406*lgu*mmsb2*pow4(g3))/(9.*(-mmgl + mmsb2)) - (4*lb1u*lgu*mmsb2*pow4(g3))
     /(9.*(-mmgl + mmsb2)) - (59*lb2u*lgu*mmsb2*pow4(g3))/(9.*(-mmgl + mmsb2))
     - (8*mmst1*pow4(g3))/(-mmgl + mmsb1) - (2*lgu*mmst1*pow4(g3))/(-mmgl +
     mmsb1) - (8*mmst1*pow4(g3))/(-mmgl + mmsb2) - (2*lgu*mmst1*pow4(g3))/(-
     mmgl + mmsb2) - (8*mmst2*pow4(g3))/(-mmgl + mmsb1) - (2*lgu*mmst2*pow4(g3)
     )/(-mmgl + mmsb1) - (8*mmst2*pow4(g3))/(-mmgl + mmsb2) - (2*lgu*mmst2*pow4
     (g3))/(-mmgl + mmsb2) + (48*mmsusy*pow4(g3))/(-mmgl + mmsb1) - (16*lgu*
     mmsusy*pow4(g3))/(-mmgl + mmsb1) + (48*mmsusy*pow4(g3))/(-mmgl + mmsb2) -
     (16*lgu*mmsusy*pow4(g3))/(-mmgl + mmsb2) - (12*mmt*pow4(g3))/(-mmgl +
     mmsb1) + (4*lgu*mmt*pow4(g3))/(-mmgl + mmsb1) - (12*mmt*pow4(g3))/(-mmgl +
     mmsb2) + (4*lgu*mmt*pow4(g3))/(-mmgl + mmsb2) - 8*zt2*pow4(g3) - (31*mmsb1
     *zt2*pow4(g3))/(9.*(-mmgl + mmsb1)) + (43*mmsb2*zt2*pow4(g3))/(3.*(-mmgl +
     mmsb1)) + (mmsb1*zt2*pow4(g3))/(-mmgl + mmsb2) - (151*mmsb2*zt2*pow4(g3))/
     (9.*(-mmgl + mmsb2)) - (mmst1*zt2*pow4(g3))/(-mmgl + mmsb1) - (mmst1*zt2*
     pow4(g3))/(-mmgl + mmsb2) - (mmst2*zt2*pow4(g3))/(-mmgl + mmsb1) - (mmst2*
     zt2*pow4(g3))/(-mmgl + mmsb2) + (8*mmsusy*zt2*pow4(g3))/(-mmgl + mmsb1) +
     (8*mmsusy*zt2*pow4(g3))/(-mmgl + mmsb2) - (2*mmt*zt2*pow4(g3))/(-mmgl +
     mmsb1) - (2*mmt*zt2*pow4(g3))/(-mmgl + mmsb2) - (56*mmsb1*mmst1*mmt*
     DeltaInv(mmt,mmsb1,mmst1)*pow4(g3))/(3.*(-mmgl + mmsb1)) + (2*lb1u*mmsb1*
     mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow4(g3))/(-mmgl + mmsb1) - (8*mmsb1*
     mmst1*mmt*zt2*DeltaInv(mmt,mmsb1,mmst1)*pow4(g3))/(3.*(-mmgl + mmsb1)) - (
     56*mmsb1*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow4(g3))/(3.*(-mmgl + mmsb1)
     ) + (2*lb1u*mmsb1*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow4(g3))/(-mmgl +
     mmsb1) - (8*mmsb1*mmst2*mmt*zt2*DeltaInv(mmt,mmsb1,mmst2)*pow4(g3))/(3.*(-
     mmgl + mmsb1)) - (56*mmsb2*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow4(g3))/(
     3.*(-mmgl + mmsb2)) + (2*lb2u*mmsb2*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1)*
     pow4(g3))/(-mmgl + mmsb2) - (8*mmsb2*mmst1*mmt*zt2*DeltaInv(mmt,mmsb2,
     mmst1)*pow4(g3))/(3.*(-mmgl + mmsb2)) - (56*mmsb2*mmst2*mmt*DeltaInv(mmt,
     mmsb2,mmst2)*pow4(g3))/(3.*(-mmgl + mmsb2)) + (2*lb2u*mmsb2*mmst2*mmt*
     DeltaInv(mmt,mmsb2,mmst2)*pow4(g3))/(-mmgl + mmsb2) - (8*mmsb2*mmst2*mmt*
     zt2*DeltaInv(mmt,mmsb2,mmst2)*pow4(g3))/(3.*(-mmgl + mmsb2)) - (56*mmgl*
     mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/3. + 8*lgu*mmgl*mmsb1*DeltaInv(
     mmt,mmst1,mmgl)*pow4(g3) - (56*mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow4(g3
     ))/3. + 8*lgu*mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow4(g3) + (56*mmgl*
     mmst1*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/3. - 16*lgu*mmgl*mmst1*DeltaInv(
     mmt,mmst1,mmgl)*pow4(g3) + (56*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow4(
     g3))/3. - 16*lgu*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow4(g3) + (56*mmsb2
     *mmst1*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/3. - 16*lgu*mmsb2*mmst1*DeltaInv
     (mmt,mmst1,mmgl)*pow4(g3) + (56*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3)
     )/3. - 8*lgu*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3) + (56*mmsb1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/3. - 8*lgu*mmsb1*mmt*DeltaInv(mmt,mmst1
     ,mmgl)*pow4(g3) + (56*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/3. - 8*
     lgu*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3) + (224*mmst1*mmt*DeltaInv(
     mmt,mmst1,mmgl)*pow4(g3))/3. - 8*lgu*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow4(g3) - (224*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(3.*(-
     mmgl + mmsb1)) + (8*lgu*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))
     /(-mmgl + mmsb1) - (224*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))
     /(3.*(-mmgl + mmsb2)) + (8*lgu*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow4(g3))/(-mmgl + mmsb2) - (8*mmgl*mmsb1*zt2*DeltaInv(mmt,mmst1,mmgl)*
     pow4(g3))/3. - (8*mmgl*mmsb2*zt2*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/3. + (
     8*mmgl*mmst1*zt2*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/3. + (8*mmsb1*mmst1*
     zt2*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/3. + (8*mmsb2*mmst1*zt2*DeltaInv(
     mmt,mmst1,mmgl)*pow4(g3))/3. + (8*mmgl*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*
     pow4(g3))/3. + (8*mmsb1*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/3. + (8
     *mmsb2*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/3. + (32*mmst1*mmt*zt2*
     DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/3. - (32*mmsb1*mmst1*mmt*zt2*DeltaInv(
     mmt,mmst1,mmgl)*pow4(g3))/(3.*(-mmgl + mmsb1)) - (32*mmsb2*mmst1*mmt*zt2*
     DeltaInv(mmt,mmst1,mmgl)*pow4(g3))/(3.*(-mmgl + mmsb2)) - (56*mmgl*mmsb1*
     DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/3. + 8*lgu*mmgl*mmsb1*DeltaInv(mmt,
     mmst2,mmgl)*pow4(g3) - (56*mmgl*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/
     3. + 8*lgu*mmgl*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow4(g3) + (56*mmgl*mmst2*
     DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/3. - 16*lgu*mmgl*mmst2*DeltaInv(mmt,
     mmst2,mmgl)*pow4(g3) + (56*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/
     3. - 16*lgu*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow4(g3) + (56*mmsb2*
     mmst2*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/3. - 16*lgu*mmsb2*mmst2*DeltaInv(
     mmt,mmst2,mmgl)*pow4(g3) + (56*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))
     /3. - 8*lgu*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3) + (56*mmsb1*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/3. - 8*lgu*mmsb1*mmt*DeltaInv(mmt,mmst2
     ,mmgl)*pow4(g3) + (56*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/3. - 8*
     lgu*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3) + (224*mmst2*mmt*DeltaInv(
     mmt,mmst2,mmgl)*pow4(g3))/3. - 8*lgu*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*
     pow4(g3) - (224*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*(-
     mmgl + mmsb1)) + (8*lgu*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))
     /(-mmgl + mmsb1) - (224*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))
     /(3.*(-mmgl + mmsb2)) + (8*lgu*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*
     pow4(g3))/(-mmgl + mmsb2) - (8*mmgl*mmsb1*zt2*DeltaInv(mmt,mmst2,mmgl)*
     pow4(g3))/3. - (8*mmgl*mmsb2*zt2*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/3. + (
     8*mmgl*mmst2*zt2*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/3. + (8*mmsb1*mmst2*
     zt2*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/3. + (8*mmsb2*mmst2*zt2*DeltaInv(
     mmt,mmst2,mmgl)*pow4(g3))/3. + (8*mmgl*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*
     pow4(g3))/3. + (8*mmsb1*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/3. + (8
     *mmsb2*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/3. + (32*mmst2*mmt*zt2*
     DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/3. - (32*mmsb1*mmst2*mmt*zt2*DeltaInv(
     mmt,mmst2,mmgl)*pow4(g3))/(3.*(-mmgl + mmsb1)) - (32*mmsb2*mmst2*mmt*zt2*
     DeltaInv(mmt,mmst2,mmgl)*pow4(g3))/(3.*(-mmgl + mmsb2)) + (83*Fin20(mmsb1,
     mmgl,mmu)*pow4(g3))/(9.*(-mmgl + mmsb1)) - (128*Fin20(mmsb1,mmgl,mmu)*pow4
     (g3))/(9.*(mmsb1 - mmsb2)) + (40*mmsb2*Fin20(mmsb1,mmgl,mmu)*pow4(g3))/(3.
     *(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (7*Fin20(mmsb1,mmgl,mmu)*pow4(g3))/(9.
     *(-mmgl + mmsb2)) + (8*mmsb2*Fin20(mmsb1,mmgl,mmu)*pow4(g3))/(9.*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) + (10*Fin20(mmsb1,mmsb2,mmu)*pow4(g3))/(9.*(-mmgl
     + mmsb1)) + (8*mmsb2*Fin20(mmsb1,mmsb2,mmu)*pow4(g3))/(9.*(-mmgl + mmsb1)*
     (mmsb1 - mmsb2)) + (2*Fin20(mmsb1,mmsb2,mmu)*pow4(g3))/(9.*(-mmgl + mmsb2)
     ) - (8*mmsb2*Fin20(mmsb1,mmsb2,mmu)*pow4(g3))/(9.*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) - (Fin20(mmsb2,mmgl,mmu)*pow4(g3))/(9.*(-mmgl + mmsb1)) + (128*
     Fin20(mmsb2,mmgl,mmu)*pow4(g3))/(9.*(mmsb1 - mmsb2)) - (8*mmsb2*Fin20(
     mmsb2,mmgl,mmu)*pow4(g3))/(9.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (37*Fin20
     (mmsb2,mmgl,mmu)*pow4(g3))/(9.*(-mmgl + mmsb2)) - (40*mmsb2*Fin20(mmsb2,
     mmgl,mmu)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (61*pow2(lb1u)*
     pow4(g3))/9. + (9*mmsb1*pow2(lb1u)*pow4(g3))/(2.*(-mmgl + mmsb1)) + (64*
     mmsb2*pow2(lb1u)*pow4(g3))/(9.*(-mmgl + mmsb1)) - (64*mmsb2*pow2(lb1u)*
     pow4(g3))/(9.*(mmsb1 - mmsb2)) + (mmsb1*pow2(lb1u)*pow4(g3))/(2.*(-mmgl +
     mmsb2)) - (mmsb1*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(lb1u)*pow4(g3))/
     (3.*(-mmgl + mmsb1)) - (mmsb1*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(
     lb1u)*pow4(g3))/(3.*(-mmgl + mmsb1)) + (pow2(lb2u)*pow4(g3))/3. + (mmsb2*
     pow2(lb2u)*pow4(g3))/(2.*(-mmgl + mmsb1)) + (64*mmsb2*pow2(lb2u)*pow4(g3))
     /(9.*(mmsb1 - mmsb2)) - (47*mmsb2*pow2(lb2u)*pow4(g3))/(18.*(-mmgl + mmsb2
     )) - (mmsb2*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(lb2u)*pow4(g3))/(3.*(
     -mmgl + mmsb2)) - (mmsb2*mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(lb2u)*
     pow4(g3))/(3.*(-mmgl + mmsb2)) - (254*pow2(lgu)*pow4(g3))/9. + (85*mmsb1*
     pow2(lgu)*pow4(g3))/(18.*(-mmgl + mmsb1)) + (353*mmsb2*pow2(lgu)*pow4(g3))
     /(18.*(-mmgl + mmsb1)) + (mmsb1*pow2(lgu)*pow4(g3))/(2.*(-mmgl + mmsb2)) -
     (259*mmsb2*pow2(lgu)*pow4(g3))/(18.*(-mmgl + mmsb2)) + (mmst1*pow2(lgu)*
     pow4(g3))/(2.*(-mmgl + mmsb1)) + (mmst1*pow2(lgu)*pow4(g3))/(2.*(-mmgl +
     mmsb2)) + (mmst2*pow2(lgu)*pow4(g3))/(2.*(-mmgl + mmsb1)) + (mmst2*pow2(
     lgu)*pow4(g3))/(2.*(-mmgl + mmsb2)) + (4*mmsusy*pow2(lgu)*pow4(g3))/(-mmgl
      + mmsb1) + (4*mmsusy*pow2(lgu)*pow4(g3))/(-mmgl + mmsb2) - (mmt*pow2(lgu)
     *pow4(g3))/(-mmgl + mmsb1) - (mmt*pow2(lgu)*pow4(g3))/(-mmgl + mmsb2) - (4
     *mmgl*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow4(g3))/3. - (4*mmgl*
     mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow4(g3))/3. + (8*mmgl*mmst1*
     DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow4(g3))/3. + (8*mmsb1*mmst1*DeltaInv(
     mmt,mmst1,mmgl)*pow2(lgu)*pow4(g3))/3. + (8*mmsb2*mmst1*DeltaInv(mmt,mmst1
     ,mmgl)*pow2(lgu)*pow4(g3))/3. + (4*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     lgu)*pow4(g3))/3. + (4*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow4(
     g3))/3. + (4*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow4(g3))/3. + (
     4*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow4(g3))/3. - (4*mmsb1*
     mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow4(g3))/(3.*(-mmgl + mmsb1)
     ) - (4*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow4(g3))/(3.*(-
     mmgl + mmsb2)) - (4*mmgl*mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow4(g3)
     )/3. - (4*mmgl*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow4(g3))/3. + (8*
     mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow4(g3))/3. + (8*mmsb1*
     mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow4(g3))/3. + (8*mmsb2*mmst2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow4(g3))/3. + (4*mmgl*mmt*DeltaInv(mmt
     ,mmst2,mmgl)*pow2(lgu)*pow4(g3))/3. + (4*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl
     )*pow2(lgu)*pow4(g3))/3. + (4*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)
     *pow4(g3))/3. + (4*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow4(g3))/
     3. - (4*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow4(g3))/(3.*(
     -mmgl + mmsb1)) - (4*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*
     pow4(g3))/(3.*(-mmgl + mmsb2)) - (56*DeltaInv(mmt,mmst1,mmgl)*pow2(mmgl)*
     pow4(g3))/3. + 8*lgu*DeltaInv(mmt,mmst1,mmgl)*pow2(mmgl)*pow4(g3) - (8*zt2
     *DeltaInv(mmt,mmst1,mmgl)*pow2(mmgl)*pow4(g3))/3. - (56*DeltaInv(mmt,mmst2
     ,mmgl)*pow2(mmgl)*pow4(g3))/3. + 8*lgu*DeltaInv(mmt,mmst2,mmgl)*pow2(mmgl)
     *pow4(g3) - (8*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmgl)*pow4(g3))/3. - (4*
     DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow2(mmgl)*pow4(g3))/3. - (4*DeltaInv(
     mmt,mmst2,mmgl)*pow2(lgu)*pow2(mmgl)*pow4(g3))/3. - (14*mmst1*DeltaInv(mmt
     ,mmsb1,mmst1)*pow2(mmsb1)*pow4(g3))/(3.*(-mmgl + mmsb1)) + (4*lb1u*mmst1*
     DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow4(g3))/(-mmgl + mmsb1) - (14*mmt*
     DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow4(g3))/(3.*(-mmgl + mmsb1)) + (2*
     lb1u*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow4(g3))/(-mmgl + mmsb1) -
     (2*mmst1*zt2*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow4(g3))/(3.*(-mmgl +
     mmsb1)) - (2*mmt*zt2*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow4(g3))/(3.*(
     -mmgl + mmsb1)) - (14*mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*pow4(g3)
     )/(3.*(-mmgl + mmsb1)) + (4*lb1u*mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow2(
     mmsb1)*pow4(g3))/(-mmgl + mmsb1) - (14*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(
     mmsb1)*pow4(g3))/(3.*(-mmgl + mmsb1)) + (2*lb1u*mmt*DeltaInv(mmt,mmsb1,
     mmst2)*pow2(mmsb1)*pow4(g3))/(-mmgl + mmsb1) - (2*mmst2*zt2*DeltaInv(mmt,
     mmsb1,mmst2)*pow2(mmsb1)*pow4(g3))/(3.*(-mmgl + mmsb1)) - (2*mmt*zt2*
     DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*pow4(g3))/(3.*(-mmgl + mmsb1)) - 28*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow4(g3) + 12*lgu*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmsb1)*pow4(g3) - (28*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1
     )*pow4(g3))/(-mmgl + mmsb1) + (24*lgu*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmsb1)*pow4(g3))/(-mmgl + mmsb1) - (28*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmsb1)*pow4(g3))/(-mmgl + mmsb1) + (12*lgu*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmsb1)*pow4(g3))/(-mmgl + mmsb1) - 4*zt2*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmsb1)*pow4(g3) - (4*mmst1*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*
     pow4(g3))/(-mmgl + mmsb1) - (4*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1
     )*pow4(g3))/(-mmgl + mmsb1) - 28*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4
     (g3) + 12*lgu*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3) - (28*mmst2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3))/(-mmgl + mmsb1) + (24*lgu*
     mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3))/(-mmgl + mmsb1) - (28
     *mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3))/(-mmgl + mmsb1) + (12*
     lgu*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3))/(-mmgl + mmsb1) - 4
     *zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3) - (4*mmst2*zt2*DeltaInv
     (mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3))/(-mmgl + mmsb1) - (4*mmt*zt2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3))/(-mmgl + mmsb1) - (2*mmst1*
     DeltaInv(mmt,mmsb1,mmst1)*pow2(lb1u)*pow2(mmsb1)*pow4(g3))/(3.*(-mmgl +
     mmsb1)) - (mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(lb1u)*pow2(mmsb1)*pow4(g3))/
     (3.*(-mmgl + mmsb1)) - (2*mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow2(lb1u)*pow2(
     mmsb1)*pow4(g3))/(3.*(-mmgl + mmsb1)) - (mmt*DeltaInv(mmt,mmsb1,mmst2)*
     pow2(lb1u)*pow2(mmsb1)*pow4(g3))/(3.*(-mmgl + mmsb1)) - 2*DeltaInv(mmt,
     mmst1,mmgl)*pow2(lgu)*pow2(mmsb1)*pow4(g3) - (4*mmst1*DeltaInv(mmt,mmst1,
     mmgl)*pow2(lgu)*pow2(mmsb1)*pow4(g3))/(-mmgl + mmsb1) - (2*mmt*DeltaInv(
     mmt,mmst1,mmgl)*pow2(lgu)*pow2(mmsb1)*pow4(g3))/(-mmgl + mmsb1) - 2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(mmsb1)*pow4(g3) - (4*mmst2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(mmsb1)*pow4(g3))/(-mmgl + mmsb1) -
     (2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(mmsb1)*pow4(g3))/(-mmgl +
     mmsb1) - (8*mmsb1*mmsb2*pow4(g3))/pow2(-mmgl + mmsb1) - (14*lb1u*mmsb1*
     mmsb2*pow4(g3))/(9.*pow2(-mmgl + mmsb1)) + (16*lb2u*mmsb1*mmsb2*pow4(g3))/
     (3.*pow2(-mmgl + mmsb1)) + (8*lb1u*lb2u*mmsb1*mmsb2*pow4(g3))/(9.*pow2(-
     mmgl + mmsb1)) + (38*lgu*mmsb1*mmsb2*pow4(g3))/(9.*pow2(-mmgl + mmsb1)) -
     (8*lb1u*lgu*mmsb1*mmsb2*pow4(g3))/(9.*pow2(-mmgl + mmsb1)) - (29*lb2u*lgu*
     mmsb1*mmsb2*pow4(g3))/(9.*pow2(-mmgl + mmsb1)) + (32*mmsb1*mmst1*pow4(g3))
     /(3.*pow2(-mmgl + mmsb1)) - (2*lb1u*mmsb1*mmst1*pow4(g3))/pow2(-mmgl +
     mmsb1) + (14*lgu*mmsb1*mmst1*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) + (32*
     mmsb1*mmst2*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) - (2*lb1u*mmsb1*mmst2*pow4(
     g3))/pow2(-mmgl + mmsb1) + (14*lgu*mmsb1*mmst2*pow4(g3))/(3.*pow2(-mmgl +
     mmsb1)) - (64*mmsb1*mmsusy*pow4(g3))/pow2(-mmgl + mmsb1) - (16*lb1u*mmsb1*
     mmsusy*pow4(g3))/pow2(-mmgl + mmsb1) + (112*lgu*mmsb1*mmsusy*pow4(g3))/(3.
     *pow2(-mmgl + mmsb1)) + (16*mmsb1*mmt*pow4(g3))/pow2(-mmgl + mmsb1) + (4*
     lb1u*mmsb1*mmt*pow4(g3))/pow2(-mmgl + mmsb1) - (28*lgu*mmsb1*mmt*pow4(g3))
     /(3.*pow2(-mmgl + mmsb1)) - (4*mmsb1*mmsb2*zt2*pow4(g3))/(3.*pow2(-mmgl +
     mmsb1)) + (4*mmsb1*mmst1*zt2*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) + (4*mmsb1
     *mmst2*zt2*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) - (32*mmsb1*mmsusy*zt2*pow4(
     g3))/(3.*pow2(-mmgl + mmsb1)) + (8*mmsb1*mmt*zt2*pow4(g3))/(3.*pow2(-mmgl
     + mmsb1)) + (12*mmsb1*Fin20(mmsb1,mmgl,mmu)*pow4(g3))/pow2(-mmgl + mmsb1)
     - (7*mmsb1*Fin20(mmsb1,mmsb2,mmu)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) + (
     mmsb2*Fin20(mmsb1,mmsb2,mmu)*pow4(g3))/pow2(-mmgl + mmsb1) + (mmsb1*Fin20(
     mmsb2,mmgl,mmu)*pow4(g3))/pow2(-mmgl + mmsb1) - (mmsb2*Fin20(mmsb2,mmgl,
     mmu)*pow4(g3))/pow2(-mmgl + mmsb1) + (mmsb1*mmsb2*pow2(lb1u)*pow4(g3))/(2.
     *pow2(-mmgl + mmsb1)) + (mmsb1*mmst1*pow2(lb1u)*pow4(g3))/(2.*pow2(-mmgl +
     mmsb1)) + (mmsb1*mmst2*pow2(lb1u)*pow4(g3))/(2.*pow2(-mmgl + mmsb1)) + (4*
     mmsb1*mmsusy*pow2(lb1u)*pow4(g3))/pow2(-mmgl + mmsb1) - (mmsb1*mmt*pow2(
     lb1u)*pow4(g3))/pow2(-mmgl + mmsb1) - (2*mmsb1*mmsb2*pow2(lb2u)*pow4(g3))/
     (3.*pow2(-mmgl + mmsb1)) - (5*mmsb1*mmsb2*pow2(lgu)*pow4(g3))/(18.*pow2(-
     mmgl + mmsb1)) - (7*mmsb1*mmst1*pow2(lgu)*pow4(g3))/(6.*pow2(-mmgl + mmsb1
     )) - (7*mmsb1*mmst2*pow2(lgu)*pow4(g3))/(6.*pow2(-mmgl + mmsb1)) - (28*
     mmsb1*mmsusy*pow2(lgu)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) + (7*mmsb1*mmt*
     pow2(lgu)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) + (1222*pow2(mmsb1)*pow4(g3))
     /(9.*pow2(-mmgl + mmsb1)) + (160*lb1u*pow2(mmsb1)*pow4(g3))/(9.*pow2(-mmgl
      + mmsb1)) - (368*lgu*pow2(mmsb1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) - (
     113*lb1u*lgu*pow2(mmsb1)*pow4(g3))/(9.*pow2(-mmgl + mmsb1)) + (20*zt2*pow2
     (mmsb1)*pow4(g3))/pow2(-mmgl + mmsb1) + (112*mmst1*mmt*DeltaInv(mmt,mmsb1,
     mmst1)*pow2(mmsb1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) - (4*lb1u*mmst1*mmt*
     DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow4(g3))/pow2(-mmgl + mmsb1) + (16*
     mmst1*mmt*zt2*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow4(g3))/(3.*pow2(-
     mmgl + mmsb1)) + (112*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*pow4
     (g3))/(3.*pow2(-mmgl + mmsb1)) - (4*lb1u*mmst2*mmt*DeltaInv(mmt,mmsb1,
     mmst2)*pow2(mmsb1)*pow4(g3))/pow2(-mmgl + mmsb1) + (16*mmst2*mmt*zt2*
     DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) +
     (112*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow4(g3))/(3.*pow2(-
     mmgl + mmsb1)) - (4*lgu*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*
     pow4(g3))/pow2(-mmgl + mmsb1) + (16*mmst1*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)
     *pow2(mmsb1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) + (112*mmst2*mmt*DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) - (4*lgu*
     mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3))/pow2(-mmgl +
     mmsb1) + (16*mmst2*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow4(g3))/
     (3.*pow2(-mmgl + mmsb1)) - (17*pow2(lb1u)*pow2(mmsb1)*pow4(g3))/(6.*pow2(-
     mmgl + mmsb1)) + (2*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(lb1u)*pow2(
     mmsb1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) + (2*mmst2*mmt*DeltaInv(mmt,
     mmsb1,mmst2)*pow2(lb1u)*pow2(mmsb1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) + (
     853*pow2(lgu)*pow2(mmsb1)*pow4(g3))/(18.*pow2(-mmgl + mmsb1)) + (2*mmst1*
     mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow2(mmsb1)*pow4(g3))/(3.*pow2(-
     mmgl + mmsb1)) + (2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(
     mmsb1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) + (836*pow2(mmsb2)*pow4(g3))/(9.
     *(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (128*lb1u*pow2(mmsb2)*pow4(g3))/(9.*(-
     mmgl + mmsb1)*(mmsb1 - mmsb2)) - (4*lb2u*pow2(mmsb2)*pow4(g3))/(3.*(-mmgl
     + mmsb1)*(mmsb1 - mmsb2)) + (4*lb1u*lb2u*pow2(mmsb2)*pow4(g3))/(9.*(-mmgl
     + mmsb1)*(mmsb1 - mmsb2)) - (580*lgu*pow2(mmsb2)*pow4(g3))/(9.*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) - (4*lb1u*lgu*pow2(mmsb2)*pow4(g3))/(9.*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) - (4*lb2u*lgu*pow2(mmsb2)*pow4(g3))/(9.*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) - (836*pow2(mmsb2)*pow4(g3))/(9.*(mmsb1 - mmsb2)*(
     -mmgl + mmsb2)) + (140*lb2u*pow2(mmsb2)*pow4(g3))/(9.*(mmsb1 - mmsb2)*(-
     mmgl + mmsb2)) - (4*lb1u*lb2u*pow2(mmsb2)*pow4(g3))/(9.*(mmsb1 - mmsb2)*(-
     mmgl + mmsb2)) + (580*lgu*pow2(mmsb2)*pow4(g3))/(9.*(mmsb1 - mmsb2)*(-mmgl
      + mmsb2)) + (4*lb1u*lgu*pow2(mmsb2)*pow4(g3))/(9.*(mmsb1 - mmsb2)*(-mmgl
     + mmsb2)) + (4*lb2u*lgu*pow2(mmsb2)*pow4(g3))/(9.*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) + (40*zt2*pow2(mmsb2)*pow4(g3))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2
     )) - (40*zt2*pow2(mmsb2)*pow4(g3))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) -
     (14*mmst1*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2)*pow4(g3))/(3.*(-mmgl +
     mmsb2)) + (4*lb2u*mmst1*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2)*pow4(g3))/(-
     mmgl + mmsb2) - (14*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2)*pow4(g3))/(
     3.*(-mmgl + mmsb2)) + (2*lb2u*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2)*
     pow4(g3))/(-mmgl + mmsb2) - (2*mmst1*zt2*DeltaInv(mmt,mmsb2,mmst1)*pow2(
     mmsb2)*pow4(g3))/(3.*(-mmgl + mmsb2)) - (2*mmt*zt2*DeltaInv(mmt,mmsb2,
     mmst1)*pow2(mmsb2)*pow4(g3))/(3.*(-mmgl + mmsb2)) - (14*mmst2*DeltaInv(mmt
     ,mmsb2,mmst2)*pow2(mmsb2)*pow4(g3))/(3.*(-mmgl + mmsb2)) + (4*lb2u*mmst2*
     DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow4(g3))/(-mmgl + mmsb2) - (14*mmt*
     DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow4(g3))/(3.*(-mmgl + mmsb2)) + (2*
     lb2u*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow4(g3))/(-mmgl + mmsb2) -
     (2*mmst2*zt2*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow4(g3))/(3.*(-mmgl +
     mmsb2)) - (2*mmt*zt2*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow4(g3))/(3.*(
     -mmgl + mmsb2)) - 28*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(g3) + 12*
     lgu*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(g3) - (28*mmst1*DeltaInv(mmt
     ,mmst1,mmgl)*pow2(mmsb2)*pow4(g3))/(-mmgl + mmsb2) + (24*lgu*mmst1*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(g3))/(-mmgl + mmsb2) - (28*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(g3))/(-mmgl + mmsb2) + (12*lgu*
     mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(g3))/(-mmgl + mmsb2) - 4*zt2
     *DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(g3) - (4*mmst1*zt2*DeltaInv(mmt
     ,mmst1,mmgl)*pow2(mmsb2)*pow4(g3))/(-mmgl + mmsb2) - (4*mmt*zt2*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(g3))/(-mmgl + mmsb2) - 28*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmsb2)*pow4(g3) + 12*lgu*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmsb2)*pow4(g3) - (28*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow4(g3))
     /(-mmgl + mmsb2) + (24*lgu*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow4
     (g3))/(-mmgl + mmsb2) - (28*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow4(
     g3))/(-mmgl + mmsb2) + (12*lgu*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*
     pow4(g3))/(-mmgl + mmsb2) - 4*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*
     pow4(g3) - (4*mmst2*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow4(g3))/(-
     mmgl + mmsb2) - (4*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow4(g3))/
     (-mmgl + mmsb2) + (64*pow2(lb1u)*pow2(mmsb2)*pow4(g3))/(9.*(-mmgl + mmsb1)
     *(mmsb1 - mmsb2)) - (64*pow2(lb2u)*pow2(mmsb2)*pow4(g3))/(9.*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) - (2*mmst1*DeltaInv(mmt,mmsb2,mmst1)*pow2(lb2u)*
     pow2(mmsb2)*pow4(g3))/(3.*(-mmgl + mmsb2)) - (mmt*DeltaInv(mmt,mmsb2,mmst1
     )*pow2(lb2u)*pow2(mmsb2)*pow4(g3))/(3.*(-mmgl + mmsb2)) - (2*mmst2*
     DeltaInv(mmt,mmsb2,mmst2)*pow2(lb2u)*pow2(mmsb2)*pow4(g3))/(3.*(-mmgl +
     mmsb2)) - (mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(lb2u)*pow2(mmsb2)*pow4(g3))/
     (3.*(-mmgl + mmsb2)) + (20*pow2(lgu)*pow2(mmsb2)*pow4(g3))/((-mmgl + mmsb1
     )*(mmsb1 - mmsb2)) - (20*pow2(lgu)*pow2(mmsb2)*pow4(g3))/((mmsb1 - mmsb2)*
     (-mmgl + mmsb2)) - 2*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow2(mmsb2)*pow4(
     g3) - (4*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow2(mmsb2)*pow4(g3))/(-
     mmgl + mmsb2) - (2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow2(mmsb2)*pow4
     (g3))/(-mmgl + mmsb2) - 2*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(mmsb2)*
     pow4(g3) - (4*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(mmsb2)*pow4(g3
     ))/(-mmgl + mmsb2) - (2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(mmsb2)
     *pow4(g3))/(-mmgl + mmsb2) + (4*lb1u*pow2(mmsb2)*pow4(g3))/(9.*pow2(-mmgl
     + mmsb1)) + (4*lb1u*lb2u*pow2(mmsb2)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) -
     (4*lgu*pow2(mmsb2)*pow4(g3))/(9.*pow2(-mmgl + mmsb1)) - (4*lb1u*lgu*pow2(
     mmsb2)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) - (4*lb2u*lgu*pow2(mmsb2)*pow4(
     g3))/(3.*pow2(-mmgl + mmsb1)) + (4*pow2(lgu)*pow2(mmsb2)*pow4(g3))/(3.*
     pow2(-mmgl + mmsb1)) + (4*pow2(lsu)*(2 + (3*mmsusy)/(-mmgl + mmsb1) + (3*
     mmsusy)/(-mmgl + mmsb2) - (4*mmsb1*mmsusy)/pow2(-mmgl + mmsb1) - (4*mmsb2*
     mmsusy)/pow2(-mmgl + mmsb2))*pow4(g3))/3. + (4*lsu*(7.333333333333333 - (
     24*mmsusy)/(-mmgl + mmsb1) - (24*mmsusy)/(-mmgl + mmsb2) + (32*mmsb1*
     mmsusy)/pow2(-mmgl + mmsb1) + (32*mmsb2*mmsusy)/pow2(-mmgl + mmsb2))*pow4(
     g3))/3. - (8*mmsb1*mmsb2*pow4(g3))/pow2(-mmgl + mmsb2) + (16*lb1u*mmsb1*
     mmsb2*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) - (2*lb2u*mmsb1*mmsb2*pow4(g3))/
     pow2(-mmgl + mmsb2) + (14*lgu*mmsb1*mmsb2*pow4(g3))/(3.*pow2(-mmgl + mmsb2
     )) - (7*lb1u*lgu*mmsb1*mmsb2*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) + (32*
     mmsb2*mmst1*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) - (2*lb2u*mmsb2*mmst1*pow4(
     g3))/pow2(-mmgl + mmsb2) + (14*lgu*mmsb2*mmst1*pow4(g3))/(3.*pow2(-mmgl +
     mmsb2)) + (32*mmsb2*mmst2*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) - (2*lb2u*
     mmsb2*mmst2*pow4(g3))/pow2(-mmgl + mmsb2) + (14*lgu*mmsb2*mmst2*pow4(g3))/
     (3.*pow2(-mmgl + mmsb2)) - (64*mmsb2*mmsusy*pow4(g3))/pow2(-mmgl + mmsb2)
     - (16*lb2u*mmsb2*mmsusy*pow4(g3))/pow2(-mmgl + mmsb2) + (112*lgu*mmsb2*
     mmsusy*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) + (16*mmsb2*mmt*pow4(g3))/pow2(-
     mmgl + mmsb2) + (4*lb2u*mmsb2*mmt*pow4(g3))/pow2(-mmgl + mmsb2) - (28*lgu*
     mmsb2*mmt*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) - (4*mmsb1*mmsb2*zt2*pow4(g3)
     )/(3.*pow2(-mmgl + mmsb2)) + (4*mmsb2*mmst1*zt2*pow4(g3))/(3.*pow2(-mmgl +
     mmsb2)) + (4*mmsb2*mmst2*zt2*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) - (32*
     mmsb2*mmsusy*zt2*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) + (8*mmsb2*mmt*zt2*
     pow4(g3))/(3.*pow2(-mmgl + mmsb2)) - (mmsb1*Fin20(mmsb1,mmgl,mmu)*pow4(g3)
     )/pow2(-mmgl + mmsb2) + (mmsb2*Fin20(mmsb1,mmgl,mmu)*pow4(g3))/pow2(-mmgl
     + mmsb2) + (mmsb1*Fin20(mmsb1,mmsb2,mmu)*pow4(g3))/pow2(-mmgl + mmsb2) - (
     7*mmsb2*Fin20(mmsb1,mmsb2,mmu)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) + (12*
     mmsb2*Fin20(mmsb2,mmgl,mmu)*pow4(g3))/pow2(-mmgl + mmsb2) - (2*mmsb1*mmsb2
     *pow2(lb1u)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) + (mmsb1*mmsb2*pow2(lb2u)*
     pow4(g3))/(2.*pow2(-mmgl + mmsb2)) + (mmsb2*mmst1*pow2(lb2u)*pow4(g3))/(2.
     *pow2(-mmgl + mmsb2)) + (mmsb2*mmst2*pow2(lb2u)*pow4(g3))/(2.*pow2(-mmgl +
     mmsb2)) + (4*mmsb2*mmsusy*pow2(lb2u)*pow4(g3))/pow2(-mmgl + mmsb2) - (
     mmsb2*mmt*pow2(lb2u)*pow4(g3))/pow2(-mmgl + mmsb2) - (7*mmsb1*mmsb2*pow2(
     lgu)*pow4(g3))/(6.*pow2(-mmgl + mmsb2)) - (7*mmsb2*mmst1*pow2(lgu)*pow4(g3
     ))/(6.*pow2(-mmgl + mmsb2)) - (7*mmsb2*mmst2*pow2(lgu)*pow4(g3))/(6.*pow2(
     -mmgl + mmsb2)) - (28*mmsb2*mmsusy*pow2(lgu)*pow4(g3))/(3.*pow2(-mmgl +
     mmsb2)) + (7*mmsb2*mmt*pow2(lgu)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) + (
     1222*pow2(mmsb2)*pow4(g3))/(9.*pow2(-mmgl + mmsb2)) + (52*lb2u*pow2(mmsb2)
     *pow4(g3))/(3.*pow2(-mmgl + mmsb2)) - (4*lb1u*lb2u*pow2(mmsb2)*pow4(g3))/(
     9.*pow2(-mmgl + mmsb2)) - (1100*lgu*pow2(mmsb2)*pow4(g3))/(9.*pow2(-mmgl +
     mmsb2)) + (4*lb1u*lgu*pow2(mmsb2)*pow4(g3))/(9.*pow2(-mmgl + mmsb2)) - (
     109*lb2u*lgu*pow2(mmsb2)*pow4(g3))/(9.*pow2(-mmgl + mmsb2)) + (20*zt2*pow2
     (mmsb2)*pow4(g3))/pow2(-mmgl + mmsb2) + (112*mmst1*mmt*DeltaInv(mmt,mmsb2,
     mmst1)*pow2(mmsb2)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) - (4*lb2u*mmst1*mmt*
     DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2)*pow4(g3))/pow2(-mmgl + mmsb2) + (16*
     mmst1*mmt*zt2*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2)*pow4(g3))/(3.*pow2(-
     mmgl + mmsb2)) + (112*mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow4
     (g3))/(3.*pow2(-mmgl + mmsb2)) - (4*lb2u*mmst2*mmt*DeltaInv(mmt,mmsb2,
     mmst2)*pow2(mmsb2)*pow4(g3))/pow2(-mmgl + mmsb2) + (16*mmst2*mmt*zt2*
     DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) +
     (112*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow4(g3))/(3.*pow2(-
     mmgl + mmsb2)) - (4*lgu*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*
     pow4(g3))/pow2(-mmgl + mmsb2) + (16*mmst1*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)
     *pow2(mmsb2)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) + (112*mmst2*mmt*DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmsb2)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) - (4*lgu*
     mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow4(g3))/pow2(-mmgl +
     mmsb2) + (16*mmst2*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow4(g3))/
     (3.*pow2(-mmgl + mmsb2)) - (17*pow2(lb2u)*pow2(mmsb2)*pow4(g3))/(6.*pow2(-
     mmgl + mmsb2)) + (2*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(lb2u)*pow2(
     mmsb2)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) + (2*mmst2*mmt*DeltaInv(mmt,
     mmsb2,mmst2)*pow2(lb2u)*pow2(mmsb2)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) + (
     845*pow2(lgu)*pow2(mmsb2)*pow4(g3))/(18.*pow2(-mmgl + mmsb2)) + (2*mmst1*
     mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*pow2(mmsb2)*pow4(g3))/(3.*pow2(-
     mmgl + mmsb2)) + (2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(
     mmsb2)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) - (14*mmsb1*DeltaInv(mmt,mmsb1,
     mmst1)*pow2(mmst1)*pow4(g3))/(3.*(-mmgl + mmsb1)) - (2*lb1u*mmsb1*DeltaInv
     (mmt,mmsb1,mmst1)*pow2(mmst1)*pow4(g3))/(-mmgl + mmsb1) - (14*mmt*DeltaInv
     (mmt,mmsb1,mmst1)*pow2(mmst1)*pow4(g3))/(3.*(-mmgl + mmsb1)) - (2*mmsb1*
     zt2*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1)*pow4(g3))/(3.*(-mmgl + mmsb1)) -
     (2*mmt*zt2*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1)*pow4(g3))/(3.*(-mmgl +
     mmsb1)) - (14*mmsb2*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1)*pow4(g3))/(3.*(-
     mmgl + mmsb2)) - (2*lb2u*mmsb2*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1)*pow4(
     g3))/(-mmgl + mmsb2) - (14*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1)*pow4(
     g3))/(3.*(-mmgl + mmsb2)) - (2*mmsb2*zt2*DeltaInv(mmt,mmsb2,mmst1)*pow2(
     mmst1)*pow4(g3))/(3.*(-mmgl + mmsb2)) - (2*mmt*zt2*DeltaInv(mmt,mmsb2,
     mmst1)*pow2(mmst1)*pow4(g3))/(3.*(-mmgl + mmsb2)) + (56*DeltaInv(mmt,mmst1
     ,mmgl)*pow2(mmst1)*pow4(g3))/3. + 8*lgu*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmst1)*pow4(g3) - (56*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow4(g3))
     /(3.*(-mmgl + mmsb1)) - (8*lgu*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*
     pow4(g3))/(-mmgl + mmsb1) - (56*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)
     *pow4(g3))/(3.*(-mmgl + mmsb2)) - (8*lgu*mmsb2*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmst1)*pow4(g3))/(-mmgl + mmsb2) - (28*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmst1)*pow4(g3))/(3.*(-mmgl + mmsb1)) - (28*mmt*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmst1)*pow4(g3))/(3.*(-mmgl + mmsb2)) + (8*zt2*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmst1)*pow4(g3))/3. - (8*mmsb1*zt2*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmst1)*pow4(g3))/(3.*(-mmgl + mmsb1)) - (8*mmsb2*zt2*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmst1)*pow4(g3))/(3.*(-mmgl + mmsb2)) - (4*mmt*zt2*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow4(g3))/(3.*(-mmgl + mmsb1)) - (4*
     mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow4(g3))/(3.*(-mmgl + mmsb2)
     ) + (mmsb1*DeltaInv(mmt,mmsb1,mmst1)*pow2(lb1u)*pow2(mmst1)*pow4(g3))/(3.*
     (-mmgl + mmsb1)) + (mmsb2*DeltaInv(mmt,mmsb2,mmst1)*pow2(lb2u)*pow2(mmst1)
     *pow4(g3))/(3.*(-mmgl + mmsb2)) - (4*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*
     pow2(mmst1)*pow4(g3))/3. + (4*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*
     pow2(mmst1)*pow4(g3))/(3.*(-mmgl + mmsb1)) + (4*mmsb2*DeltaInv(mmt,mmst1,
     mmgl)*pow2(lgu)*pow2(mmst1)*pow4(g3))/(3.*(-mmgl + mmsb2)) + (28*mmsb1*mmt
     *DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1))
     + (4*mmsb1*mmt*zt2*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1)*pow4(g3))/(3.*
     pow2(-mmgl + mmsb1)) + (28*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*
     pow4(g3))/(3.*pow2(-mmgl + mmsb1)) + (4*mmsb1*mmt*zt2*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmst1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) + (28*DeltaInv(mmt,
     mmsb1,mmst1)*pow2(mmsb1)*pow2(mmst1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) +
     (4*lb1u*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow2(mmst1)*pow4(g3))/pow2(-
     mmgl + mmsb1) + (4*zt2*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow2(mmst1)*
     pow4(g3))/(3.*pow2(-mmgl + mmsb1)) + (28*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmsb1)*pow2(mmst1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) + (4*lgu*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(mmst1)*pow4(g3))/pow2(-mmgl + mmsb1) + (4
     *zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(mmst1)*pow4(g3))/(3.*pow2(-
     mmgl + mmsb1)) - (2*DeltaInv(mmt,mmsb1,mmst1)*pow2(lb1u)*pow2(mmsb1)*pow2(
     mmst1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) - (2*DeltaInv(mmt,mmst1,mmgl)*
     pow2(lgu)*pow2(mmsb1)*pow2(mmst1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) + (28
     *mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1)*pow4(g3))/(3.*pow2(-mmgl
     + mmsb2)) + (4*mmsb2*mmt*zt2*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1)*pow4(g3
     ))/(3.*pow2(-mmgl + mmsb2)) + (28*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmst1)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) + (4*mmsb2*mmt*zt2*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmst1)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) + (28*DeltaInv(
     mmt,mmsb2,mmst1)*pow2(mmsb2)*pow2(mmst1)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)
     ) + (4*lb2u*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2)*pow2(mmst1)*pow4(g3))/
     pow2(-mmgl + mmsb2) + (4*zt2*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2)*pow2(
     mmst1)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) + (28*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmsb2)*pow2(mmst1)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) + (4*lgu*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(mmst1)*pow4(g3))/pow2(-mmgl +
     mmsb2) + (4*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(mmst1)*pow4(g3))
     /(3.*pow2(-mmgl + mmsb2)) - (2*DeltaInv(mmt,mmsb2,mmst1)*pow2(lb2u)*pow2(
     mmsb2)*pow2(mmst1)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) - (2*DeltaInv(mmt,
     mmst1,mmgl)*pow2(lgu)*pow2(mmsb2)*pow2(mmst1)*pow4(g3))/(3.*pow2(-mmgl +
     mmsb2)) - (14*mmsb1*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2)*pow4(g3))/(3.*(-
     mmgl + mmsb1)) - (2*lb1u*mmsb1*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2)*pow4(
     g3))/(-mmgl + mmsb1) - (14*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2)*pow4(
     g3))/(3.*(-mmgl + mmsb1)) - (2*mmsb1*zt2*DeltaInv(mmt,mmsb1,mmst2)*pow2(
     mmst2)*pow4(g3))/(3.*(-mmgl + mmsb1)) - (2*mmt*zt2*DeltaInv(mmt,mmsb1,
     mmst2)*pow2(mmst2)*pow4(g3))/(3.*(-mmgl + mmsb1)) - (14*mmsb2*DeltaInv(mmt
     ,mmsb2,mmst2)*pow2(mmst2)*pow4(g3))/(3.*(-mmgl + mmsb2)) - (2*lb2u*mmsb2*
     DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2)*pow4(g3))/(-mmgl + mmsb2) - (14*mmt*
     DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2)*pow4(g3))/(3.*(-mmgl + mmsb2)) - (2*
     mmsb2*zt2*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2)*pow4(g3))/(3.*(-mmgl +
     mmsb2)) - (2*mmt*zt2*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2)*pow4(g3))/(3.*(
     -mmgl + mmsb2)) + (56*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/3. +
     8*lgu*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3) - (56*mmsb1*DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(3.*(-mmgl + mmsb1)) - (8*lgu*mmsb1*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(-mmgl + mmsb1) - (56*mmsb2
     *DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(3.*(-mmgl + mmsb2)) - (8*
     lgu*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(-mmgl + mmsb2) -
     (28*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(3.*(-mmgl + mmsb1)
     ) - (28*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(3.*(-mmgl +
     mmsb2)) + (8*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/3. - (8*
     mmsb1*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(3.*(-mmgl +
     mmsb1)) - (8*mmsb2*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(3.*
     (-mmgl + mmsb2)) - (4*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3
     ))/(3.*(-mmgl + mmsb1)) - (4*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*
     pow4(g3))/(3.*(-mmgl + mmsb2)) + (mmsb1*DeltaInv(mmt,mmsb1,mmst2)*pow2(
     lb1u)*pow2(mmst2)*pow4(g3))/(3.*(-mmgl + mmsb1)) + (mmsb2*DeltaInv(mmt,
     mmsb2,mmst2)*pow2(lb2u)*pow2(mmst2)*pow4(g3))/(3.*(-mmgl + mmsb2)) - (4*
     DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(mmst2)*pow4(g3))/3. + (4*mmsb1*
     DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(mmst2)*pow4(g3))/(3.*(-mmgl +
     mmsb1)) + (4*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(mmst2)*pow4(g3)
     )/(3.*(-mmgl + mmsb2)) + (28*mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(
     mmst2)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) + (4*mmsb1*mmt*zt2*DeltaInv(mmt,
     mmsb1,mmst2)*pow2(mmst2)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) + (28*mmsb1*
     mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)
     ) + (4*mmsb1*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(3.*
     pow2(-mmgl + mmsb1)) + (28*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*pow2(
     mmst2)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) + (4*lb1u*DeltaInv(mmt,mmsb1,
     mmst2)*pow2(mmsb1)*pow2(mmst2)*pow4(g3))/pow2(-mmgl + mmsb1) + (4*zt2*
     DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*pow2(mmst2)*pow4(g3))/(3.*pow2(-mmgl
      + mmsb1)) + (28*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(mmst2)*pow4(g3)
     )/(3.*pow2(-mmgl + mmsb1)) + (4*lgu*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*
     pow2(mmst2)*pow4(g3))/pow2(-mmgl + mmsb1) + (4*zt2*DeltaInv(mmt,mmst2,mmgl
     )*pow2(mmsb1)*pow2(mmst2)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) - (2*DeltaInv
     (mmt,mmsb1,mmst2)*pow2(lb1u)*pow2(mmsb1)*pow2(mmst2)*pow4(g3))/(3.*pow2(-
     mmgl + mmsb1)) - (2*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow2(mmsb1)*pow2(
     mmst2)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) + (28*mmsb2*mmt*DeltaInv(mmt,
     mmsb2,mmst2)*pow2(mmst2)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) + (4*mmsb2*mmt
     *zt2*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2)*pow4(g3))/(3.*pow2(-mmgl +
     mmsb2)) + (28*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow4(g3))/(3.
     *pow2(-mmgl + mmsb2)) + (4*mmsb2*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmst2)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) + (28*DeltaInv(mmt,mmsb2,mmst2)*
     pow2(mmsb2)*pow2(mmst2)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) + (4*lb2u*
     DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow2(mmst2)*pow4(g3))/pow2(-mmgl +
     mmsb2) + (4*zt2*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow2(mmst2)*pow4(g3)
     )/(3.*pow2(-mmgl + mmsb2)) + (28*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2
     (mmst2)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) + (4*lgu*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmsb2)*pow2(mmst2)*pow4(g3))/pow2(-mmgl + mmsb2) + (4*zt2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(mmst2)*pow4(g3))/(3.*pow2(-mmgl
     + mmsb2)) - (2*DeltaInv(mmt,mmsb2,mmst2)*pow2(lb2u)*pow2(mmsb2)*pow2(mmst2
     )*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) - (2*DeltaInv(mmt,mmst2,mmgl)*pow2(
     lgu)*pow2(mmsb2)*pow2(mmst2)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) + (14*
     DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1)*pow4(g3))/(3.*(-mmgl + mmsb1)) - (2*
     lb1u*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1)*pow4(g3))/(-mmgl + mmsb1) + (2*
     zt2*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1)*pow4(g3))/(3.*(-mmgl + mmsb1)) +
     (14*DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1)*pow4(g3))/(3.*(-mmgl + mmsb1)) -
     (2*lb1u*DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1)*pow4(g3))/(-mmgl + mmsb1) +
     (2*zt2*DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1)*pow4(g3))/(3.*(-mmgl + mmsb1)
     ) + (112*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1)*pow4(g3))/(3.*(-mmgl + mmsb1
     )) - (16*lgu*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1)*pow4(g3))/(-mmgl + mmsb1
     ) + (16*zt2*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1)*pow4(g3))/(3.*(-mmgl +
     mmsb1)) + (112*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1)*pow4(g3))/(3.*(-mmgl +
     mmsb1)) - (16*lgu*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1)*pow4(g3))/(-mmgl +
     mmsb1) + (16*zt2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1)*pow4(g3))/(3.*(-mmgl
      + mmsb1)) + (DeltaInv(mmt,mmsb1,mmst1)*pow2(lb1u)*pow3(mmsb1)*pow4(g3))/(
     3.*(-mmgl + mmsb1)) + (DeltaInv(mmt,mmsb1,mmst2)*pow2(lb1u)*pow3(mmsb1)*
     pow4(g3))/(3.*(-mmgl + mmsb1)) + (8*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*
     pow3(mmsb1)*pow4(g3))/(3.*(-mmgl + mmsb1)) + (8*DeltaInv(mmt,mmst2,mmgl)*
     pow2(lgu)*pow3(mmsb1)*pow4(g3))/(3.*(-mmgl + mmsb1)) + (28*mmst1*DeltaInv(
     mmt,mmsb1,mmst1)*pow3(mmsb1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) - (8*lb1u*
     mmst1*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1)*pow4(g3))/pow2(-mmgl + mmsb1)
     + (28*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1)*pow4(g3))/(3.*pow2(-mmgl +
     mmsb1)) - (4*lb1u*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1)*pow4(g3))/pow2
     (-mmgl + mmsb1) + (4*mmst1*zt2*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1)*pow4(
     g3))/(3.*pow2(-mmgl + mmsb1)) + (4*mmt*zt2*DeltaInv(mmt,mmsb1,mmst1)*pow3(
     mmsb1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) + (28*mmst2*DeltaInv(mmt,mmsb1,
     mmst2)*pow3(mmsb1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) - (8*lb1u*mmst2*
     DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1)*pow4(g3))/pow2(-mmgl + mmsb1) + (28*
     mmt*DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1
     )) - (4*lb1u*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1)*pow4(g3))/pow2(-
     mmgl + mmsb1) + (4*mmst2*zt2*DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1)*pow4(g3
     ))/(3.*pow2(-mmgl + mmsb1)) + (4*mmt*zt2*DeltaInv(mmt,mmsb1,mmst2)*pow3(
     mmsb1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) + (28*mmst1*DeltaInv(mmt,mmst1,
     mmgl)*pow3(mmsb1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) - (8*lgu*mmst1*
     DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1)*pow4(g3))/pow2(-mmgl + mmsb1) + (28*
     mmt*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)
     ) - (4*lgu*mmt*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1)*pow4(g3))/pow2(-mmgl +
     mmsb1) + (4*mmst1*zt2*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1)*pow4(g3))/(3.*
     pow2(-mmgl + mmsb1)) + (4*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1)*
     pow4(g3))/(3.*pow2(-mmgl + mmsb1)) + (28*mmst2*DeltaInv(mmt,mmst2,mmgl)*
     pow3(mmsb1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) - (8*lgu*mmst2*DeltaInv(mmt
     ,mmst2,mmgl)*pow3(mmsb1)*pow4(g3))/pow2(-mmgl + mmsb1) + (28*mmt*DeltaInv(
     mmt,mmst2,mmgl)*pow3(mmsb1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) - (4*lgu*
     mmt*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1)*pow4(g3))/pow2(-mmgl + mmsb1) + (
     4*mmst2*zt2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1)*pow4(g3))/(3.*pow2(-mmgl
     + mmsb1)) + (4*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1)*pow4(g3))/(3.*
     pow2(-mmgl + mmsb1)) + (4*mmst1*DeltaInv(mmt,mmsb1,mmst1)*pow2(lb1u)*pow3(
     mmsb1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) + (2*mmt*DeltaInv(mmt,mmsb1,
     mmst1)*pow2(lb1u)*pow3(mmsb1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) + (4*
     mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow2(lb1u)*pow3(mmsb1)*pow4(g3))/(3.*pow2(
     -mmgl + mmsb1)) + (2*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(lb1u)*pow3(mmsb1)*
     pow4(g3))/(3.*pow2(-mmgl + mmsb1)) + (4*mmst1*DeltaInv(mmt,mmst1,mmgl)*
     pow2(lgu)*pow3(mmsb1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) + (2*mmt*DeltaInv
     (mmt,mmst1,mmgl)*pow2(lgu)*pow3(mmsb1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1))
     + (4*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow3(mmsb1)*pow4(g3))/(3.*
     pow2(-mmgl + mmsb1)) + (2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow3(
     mmsb1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) + (4*Fin3(mmt,mmsb1,mmst1,mmu)*(
     -((mmsb1*mmst1*DeltaInv(mmt,mmsb1,mmst1))/(-mmgl + mmsb1)) - (mmsb1*mmt*
     DeltaInv(mmt,mmsb1,mmst1))/(2.*(-mmgl + mmsb1)) - (mmst1*mmt*DeltaInv(mmt,
     mmsb1,mmst1))/(2.*(-mmgl + mmsb1)) + (DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1
     ))/(2.*(-mmgl + mmsb1)) - (3*mmsb1)/(4.*pow2(-mmgl + mmsb1)) + (3*mmst1)/(
     4.*pow2(-mmgl + mmsb1)) - (3*mmt)/(4.*pow2(-mmgl + mmsb1)) + (mmsb1*mmst1*
     mmt*DeltaInv(mmt,mmsb1,mmst1))/pow2(-mmgl + mmsb1) + (2*mmst1*DeltaInv(mmt
     ,mmsb1,mmst1)*pow2(mmsb1))/pow2(-mmgl + mmsb1) + (mmt*DeltaInv(mmt,mmsb1,
     mmst1)*pow2(mmsb1))/pow2(-mmgl + mmsb1) + (DeltaInv(mmt,mmsb1,mmst1)*pow2(
     mmst1))/(2.*(-mmgl + mmsb1)) - (mmsb1*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1
     ))/pow2(-mmgl + mmsb1) - (DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1))/pow2(-
     mmgl + mmsb1) - (mmsb1*mmst1)/pow3(-mmgl + mmsb1) + (mmsb1*mmt)/pow3(-mmgl
      + mmsb1) + pow2(mmsb1)/pow3(-mmgl + mmsb1))*pow4(g3))/3. + (4*Fin3(mmt,
     mmsb1,mmst2,mmu)*(-((mmsb1*mmst2*DeltaInv(mmt,mmsb1,mmst2))/(-mmgl + mmsb1
     )) - (mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst2))/(2.*(-mmgl + mmsb1)) - (mmst2*
     mmt*DeltaInv(mmt,mmsb1,mmst2))/(2.*(-mmgl + mmsb1)) + (DeltaInv(mmt,mmsb1,
     mmst2)*pow2(mmsb1))/(2.*(-mmgl + mmsb1)) - (3*mmsb1)/(4.*pow2(-mmgl +
     mmsb1)) + (3*mmst2)/(4.*pow2(-mmgl + mmsb1)) - (3*mmt)/(4.*pow2(-mmgl +
     mmsb1)) + (mmsb1*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2))/pow2(-mmgl + mmsb1)
     + (2*mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))/pow2(-mmgl + mmsb1) + (
     mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))/pow2(-mmgl + mmsb1) + (DeltaInv
     (mmt,mmsb1,mmst2)*pow2(mmst2))/(2.*(-mmgl + mmsb1)) - (mmsb1*DeltaInv(mmt,
     mmsb1,mmst2)*pow2(mmst2))/pow2(-mmgl + mmsb1) - (DeltaInv(mmt,mmsb1,mmst2)
     *pow3(mmsb1))/pow2(-mmgl + mmsb1) - (mmsb1*mmst2)/pow3(-mmgl + mmsb1) + (
     mmsb1*mmt)/pow3(-mmgl + mmsb1) + pow2(mmsb1)/pow3(-mmgl + mmsb1))*pow4(g3)
     )/3. + (4*Fin20(mmsb1,mmsusy,mmu)*(4/(-mmgl + mmsb1) - (14*mmsb1)/pow2(-
     mmgl + mmsb1) + (6*mmsusy)/pow2(-mmgl + mmsb1) - (8*mmsb1*mmsusy)/pow3(-
     mmgl + mmsb1) + (8*pow2(mmsb1))/pow3(-mmgl + mmsb1))*pow4(g3))/3. + (4*
     lb1u*lt1u*((mmsb1*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1))/(-mmgl + mmsb1) + (
     3*mmsb1*mmst1)/(4.*pow2(-mmgl + mmsb1)) - (2*mmst1*mmt*DeltaInv(mmt,mmsb1,
     mmst1)*pow2(mmsb1))/pow2(-mmgl + mmsb1) - (mmst1*pow2(mmsb1))/pow3(-mmgl +
     mmsb1))*pow4(g3))/3. + (4*lb1u*lt2u*((mmsb1*mmst2*mmt*DeltaInv(mmt,mmsb1,
     mmst2))/(-mmgl + mmsb1) + (3*mmsb1*mmst2)/(4.*pow2(-mmgl + mmsb1)) - (2*
     mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))/pow2(-mmgl + mmsb1) - (
     mmst2*pow2(mmsb1))/pow3(-mmgl + mmsb1))*pow4(g3))/3. + (4*lb1u*lsu*((6*
     mmsb1*mmsusy)/pow2(-mmgl + mmsb1) - (8*mmsusy*pow2(mmsb1))/pow3(-mmgl +
     mmsb1))*pow4(g3))/3. - (4*mmsb1*mmsb2*Fin20(mmsb1,mmsb2,mmu)*pow4(g3))/(3.
     *pow3(-mmgl + mmsb1)) + (4*mmsb1*mmsb2*Fin20(mmsb2,mmgl,mmu)*pow4(g3))/(3.
     *pow3(-mmgl + mmsb1)) + (8*lb1u*mmsb2*pow2(mmsb1)*pow4(g3))/(3.*pow3(-mmgl
      + mmsb1)) - (8*lgu*mmsb2*pow2(mmsb1)*pow4(g3))/(3.*pow3(-mmgl + mmsb1)) +
     (4*lb2u*lgu*mmsb2*pow2(mmsb1)*pow4(g3))/(3.*pow3(-mmgl + mmsb1)) + (8*lb1u
     *mmst1*pow2(mmsb1)*pow4(g3))/(3.*pow3(-mmgl + mmsb1)) - (8*lgu*mmst1*pow2(
     mmsb1)*pow4(g3))/(3.*pow3(-mmgl + mmsb1)) + (8*lb1u*mmst2*pow2(mmsb1)*pow4
     (g3))/(3.*pow3(-mmgl + mmsb1)) - (8*lgu*mmst2*pow2(mmsb1)*pow4(g3))/(3.*
     pow3(-mmgl + mmsb1)) + (64*lb1u*mmsusy*pow2(mmsb1)*pow4(g3))/(3.*pow3(-
     mmgl + mmsb1)) - (64*lgu*mmsusy*pow2(mmsb1)*pow4(g3))/(3.*pow3(-mmgl +
     mmsb1)) - (16*lb1u*mmt*pow2(mmsb1)*pow4(g3))/(3.*pow3(-mmgl + mmsb1)) + (
     16*lgu*mmt*pow2(mmsb1)*pow4(g3))/(3.*pow3(-mmgl + mmsb1)) - (16*Fin20(
     mmsb1,mmgl,mmu)*pow2(mmsb1)*pow4(g3))/(3.*pow3(-mmgl + mmsb1)) + (4*Fin20(
     mmsb1,mmsb2,mmu)*pow2(mmsb1)*pow4(g3))/(3.*pow3(-mmgl + mmsb1)) - (4*Fin20
     (mmsb2,mmgl,mmu)*pow2(mmsb1)*pow4(g3))/(3.*pow3(-mmgl + mmsb1)) - (2*mmsb2
     *pow2(lb1u)*pow2(mmsb1)*pow4(g3))/(3.*pow3(-mmgl + mmsb1)) - (2*mmst1*pow2
     (lb1u)*pow2(mmsb1)*pow4(g3))/(3.*pow3(-mmgl + mmsb1)) - (2*mmst2*pow2(lb1u
     )*pow2(mmsb1)*pow4(g3))/(3.*pow3(-mmgl + mmsb1)) - (16*mmsusy*pow2(lb1u)*
     pow2(mmsb1)*pow4(g3))/(3.*pow3(-mmgl + mmsb1)) + (4*mmt*pow2(lb1u)*pow2(
     mmsb1)*pow4(g3))/(3.*pow3(-mmgl + mmsb1)) + (2*mmsb2*pow2(lgu)*pow2(mmsb1)
     *pow4(g3))/(3.*pow3(-mmgl + mmsb1)) + (2*mmst1*pow2(lgu)*pow2(mmsb1)*pow4(
     g3))/(3.*pow3(-mmgl + mmsb1)) + (2*mmst2*pow2(lgu)*pow2(mmsb1)*pow4(g3))/(
     3.*pow3(-mmgl + mmsb1)) + (16*mmsusy*pow2(lgu)*pow2(mmsb1)*pow4(g3))/(3.*
     pow3(-mmgl + mmsb1)) - (4*mmt*pow2(lgu)*pow2(mmsb1)*pow4(g3))/(3.*pow3(-
     mmgl + mmsb1)) - (112*pow3(mmsb1)*pow4(g3))/(3.*pow3(-mmgl + mmsb1)) - (
     220*lb1u*pow3(mmsb1)*pow4(g3))/(9.*pow3(-mmgl + mmsb1)) + (508*lgu*pow3(
     mmsb1)*pow4(g3))/(9.*pow3(-mmgl + mmsb1)) + (116*lb1u*lgu*pow3(mmsb1)*pow4
     (g3))/(9.*pow3(-mmgl + mmsb1)) - (16*zt2*pow3(mmsb1)*pow4(g3))/(3.*pow3(-
     mmgl + mmsb1)) + (10*pow2(lb1u)*pow3(mmsb1)*pow4(g3))/(3.*pow3(-mmgl +
     mmsb1)) - (242*pow2(lgu)*pow3(mmsb1)*pow4(g3))/(9.*pow3(-mmgl + mmsb1)) +
     (14*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmsb2)*pow4(g3))/(3.*(-mmgl + mmsb2)) -
     (2*lb2u*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmsb2)*pow4(g3))/(-mmgl + mmsb2) +
     (2*zt2*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmsb2)*pow4(g3))/(3.*(-mmgl + mmsb2)
     ) + (14*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmsb2)*pow4(g3))/(3.*(-mmgl + mmsb2
     )) - (2*lb2u*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmsb2)*pow4(g3))/(-mmgl +
     mmsb2) + (2*zt2*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmsb2)*pow4(g3))/(3.*(-mmgl
      + mmsb2)) + (112*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2)*pow4(g3))/(3.*(-
     mmgl + mmsb2)) - (16*lgu*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2)*pow4(g3))/(-
     mmgl + mmsb2) + (16*zt2*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2)*pow4(g3))/(3.
     *(-mmgl + mmsb2)) + (112*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2)*pow4(g3))/(
     3.*(-mmgl + mmsb2)) - (16*lgu*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2)*pow4(g3
     ))/(-mmgl + mmsb2) + (16*zt2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2)*pow4(g3)
     )/(3.*(-mmgl + mmsb2)) + (DeltaInv(mmt,mmsb2,mmst1)*pow2(lb2u)*pow3(mmsb2)
     *pow4(g3))/(3.*(-mmgl + mmsb2)) + (DeltaInv(mmt,mmsb2,mmst2)*pow2(lb2u)*
     pow3(mmsb2)*pow4(g3))/(3.*(-mmgl + mmsb2)) + (8*DeltaInv(mmt,mmst1,mmgl)*
     pow2(lgu)*pow3(mmsb2)*pow4(g3))/(3.*(-mmgl + mmsb2)) + (8*DeltaInv(mmt,
     mmst2,mmgl)*pow2(lgu)*pow3(mmsb2)*pow4(g3))/(3.*(-mmgl + mmsb2)) + (4*lb1u
     *pow3(mmsb2)*pow4(g3))/(9.*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) + (16*lb1u
     *lb2u*pow3(mmsb2)*pow4(g3))/(9.*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (4*
     lgu*pow3(mmsb2)*pow4(g3))/(9.*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (16*
     lb1u*lgu*pow3(mmsb2)*pow4(g3))/(9.*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) -
     (16*lb2u*lgu*pow3(mmsb2)*pow4(g3))/(9.*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)
     ) + (16*pow2(lgu)*pow3(mmsb2)*pow4(g3))/(9.*(mmsb1 - mmsb2)*pow2(-mmgl +
     mmsb1)) + (4*lb1u*pow3(mmsb2)*pow4(g3))/(9.*(-mmgl + mmsb1)*pow2(mmsb1 -
     mmsb2)) - (4*lb2u*pow3(mmsb2)*pow4(g3))/(9.*(-mmgl + mmsb1)*pow2(mmsb1 -
     mmsb2)) + (16*lb1u*lb2u*pow3(mmsb2)*pow4(g3))/(9.*(-mmgl + mmsb1)*pow2(
     mmsb1 - mmsb2)) - (16*lb1u*lgu*pow3(mmsb2)*pow4(g3))/(9.*(-mmgl + mmsb1)*
     pow2(mmsb1 - mmsb2)) - (16*lb2u*lgu*pow3(mmsb2)*pow4(g3))/(9.*(-mmgl +
     mmsb1)*pow2(mmsb1 - mmsb2)) - (4*lb1u*pow3(mmsb2)*pow4(g3))/(9.*(-mmgl +
     mmsb2)*pow2(mmsb1 - mmsb2)) + (4*lb2u*pow3(mmsb2)*pow4(g3))/(9.*(-mmgl +
     mmsb2)*pow2(mmsb1 - mmsb2)) - (16*lb1u*lb2u*pow3(mmsb2)*pow4(g3))/(9.*(-
     mmgl + mmsb2)*pow2(mmsb1 - mmsb2)) + (16*lb1u*lgu*pow3(mmsb2)*pow4(g3))/(
     9.*(-mmgl + mmsb2)*pow2(mmsb1 - mmsb2)) + (16*lb2u*lgu*pow3(mmsb2)*pow4(g3
     ))/(9.*(-mmgl + mmsb2)*pow2(mmsb1 - mmsb2)) + (16*pow2(lgu)*pow3(mmsb2)*
     pow4(g3))/(9.*(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (16*pow2(lgu)*pow3(
     mmsb2)*pow4(g3))/(9.*(-mmgl + mmsb2)*pow2(mmsb1 - mmsb2)) - (4*lb2u*pow3(
     mmsb2)*pow4(g3))/(9.*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (4*lgu*pow3(
     mmsb2)*pow4(g3))/(9.*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (28*mmst1*
     DeltaInv(mmt,mmsb2,mmst1)*pow3(mmsb2)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) -
     (8*lb2u*mmst1*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmsb2)*pow4(g3))/pow2(-mmgl +
     mmsb2) + (28*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmsb2)*pow4(g3))/(3.*pow2(
     -mmgl + mmsb2)) - (4*lb2u*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmsb2)*pow4(
     g3))/pow2(-mmgl + mmsb2) + (4*mmst1*zt2*DeltaInv(mmt,mmsb2,mmst1)*pow3(
     mmsb2)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) + (4*mmt*zt2*DeltaInv(mmt,mmsb2,
     mmst1)*pow3(mmsb2)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) + (28*mmst2*DeltaInv
     (mmt,mmsb2,mmst2)*pow3(mmsb2)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) - (8*lb2u
     *mmst2*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmsb2)*pow4(g3))/pow2(-mmgl + mmsb2)
     + (28*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmsb2)*pow4(g3))/(3.*pow2(-mmgl +
     mmsb2)) - (4*lb2u*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmsb2)*pow4(g3))/pow2
     (-mmgl + mmsb2) + (4*mmst2*zt2*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmsb2)*pow4(
     g3))/(3.*pow2(-mmgl + mmsb2)) + (4*mmt*zt2*DeltaInv(mmt,mmsb2,mmst2)*pow3(
     mmsb2)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) + (28*mmst1*DeltaInv(mmt,mmst1,
     mmgl)*pow3(mmsb2)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) - (8*lgu*mmst1*
     DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2)*pow4(g3))/pow2(-mmgl + mmsb2) + (28*
     mmt*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)
     ) - (4*lgu*mmt*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2)*pow4(g3))/pow2(-mmgl +
     mmsb2) + (4*mmst1*zt2*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2)*pow4(g3))/(3.*
     pow2(-mmgl + mmsb2)) + (4*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2)*
     pow4(g3))/(3.*pow2(-mmgl + mmsb2)) + (28*mmst2*DeltaInv(mmt,mmst2,mmgl)*
     pow3(mmsb2)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) - (8*lgu*mmst2*DeltaInv(mmt
     ,mmst2,mmgl)*pow3(mmsb2)*pow4(g3))/pow2(-mmgl + mmsb2) + (28*mmt*DeltaInv(
     mmt,mmst2,mmgl)*pow3(mmsb2)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) - (4*lgu*
     mmt*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2)*pow4(g3))/pow2(-mmgl + mmsb2) + (
     4*mmst2*zt2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2)*pow4(g3))/(3.*pow2(-mmgl
     + mmsb2)) + (4*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2)*pow4(g3))/(3.*
     pow2(-mmgl + mmsb2)) + (4*mmst1*DeltaInv(mmt,mmsb2,mmst1)*pow2(lb2u)*pow3(
     mmsb2)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) + (2*mmt*DeltaInv(mmt,mmsb2,
     mmst1)*pow2(lb2u)*pow3(mmsb2)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) + (4*
     mmst2*DeltaInv(mmt,mmsb2,mmst2)*pow2(lb2u)*pow3(mmsb2)*pow4(g3))/(3.*pow2(
     -mmgl + mmsb2)) + (2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(lb2u)*pow3(mmsb2)*
     pow4(g3))/(3.*pow2(-mmgl + mmsb2)) + (4*mmst1*DeltaInv(mmt,mmst1,mmgl)*
     pow2(lgu)*pow3(mmsb2)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) + (2*mmt*DeltaInv
     (mmt,mmst1,mmgl)*pow2(lgu)*pow3(mmsb2)*pow4(g3))/(3.*pow2(-mmgl + mmsb2))
     + (4*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow3(mmsb2)*pow4(g3))/(3.*
     pow2(-mmgl + mmsb2)) + (2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow3(
     mmsb2)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) + (4*Fin20(mmgl,mmsusy,mmu)*(2/(
     -mmgl + mmsb1) + 2/(-mmgl + mmsb2) + (6*mmsb1)/pow2(-mmgl + mmsb1) - (6*
     mmsusy)/pow2(-mmgl + mmsb1) + (6*mmsb2)/pow2(-mmgl + mmsb2) - (6*mmsusy)/
     pow2(-mmgl + mmsb2) + (8*mmsb1*mmsusy)/pow3(-mmgl + mmsb1) - (8*pow2(mmsb1
     ))/pow3(-mmgl + mmsb1) + (8*mmsb2*mmsusy)/pow3(-mmgl + mmsb2) - (8*pow2(
     mmsb2))/pow3(-mmgl + mmsb2))*pow4(g3))/3. + (4*Fin3(mmt,mmst1,mmgl,mmu)*(-
     3/(4.*(-mmgl + mmsb1)) - 3/(4.*(-mmgl + mmsb2)) - 2*mmgl*DeltaInv(mmt,
     mmst1,mmgl) - 2*mmsb1*DeltaInv(mmt,mmst1,mmgl) - 2*mmsb2*DeltaInv(mmt,
     mmst1,mmgl) + 4*mmst1*DeltaInv(mmt,mmst1,mmgl) - (4*mmsb1*mmst1*DeltaInv(
     mmt,mmst1,mmgl))/(-mmgl + mmsb1) - (4*mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)
     )/(-mmgl + mmsb2) + 2*mmt*DeltaInv(mmt,mmst1,mmgl) - (2*mmsb1*mmt*DeltaInv
     (mmt,mmst1,mmgl))/(-mmgl + mmsb1) - (2*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl))
     /(-mmgl + mmsb2) - (mmst1*mmt*DeltaInv(mmt,mmst1,mmgl))/(-mmgl + mmsb1) -
     (mmst1*mmt*DeltaInv(mmt,mmst1,mmgl))/(-mmgl + mmsb2) + (3*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmsb1))/(-mmgl + mmsb1) + (7*mmsb1)/(4.*pow2(-mmgl +
     mmsb1)) - (3*mmst1)/(4.*pow2(-mmgl + mmsb1)) + (3*mmt)/(4.*pow2(-mmgl +
     mmsb1)) + (mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl))/pow2(-mmgl + mmsb1) +
     (2*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1))/pow2(-mmgl + mmsb1) + (mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1))/pow2(-mmgl + mmsb1) + (3*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmsb2))/(-mmgl + mmsb2) + (7*mmsb2)/(4.*pow2(-mmgl +
     mmsb2)) - (3*mmst1)/(4.*pow2(-mmgl + mmsb2)) + (3*mmt)/(4.*pow2(-mmgl +
     mmsb2)) + (mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl))/pow2(-mmgl + mmsb2) +
     (2*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2))/pow2(-mmgl + mmsb2) + (mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2))/pow2(-mmgl + mmsb2) + (DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmst1))/(-mmgl + mmsb1) + (DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmst1))/(-mmgl + mmsb2) - (mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/
     pow2(-mmgl + mmsb1) - (mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/pow2(-
     mmgl + mmsb2) - (DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1))/pow2(-mmgl + mmsb1)
     + (mmsb1*mmst1)/pow3(-mmgl + mmsb1) - (mmsb1*mmt)/pow3(-mmgl + mmsb1) -
     pow2(mmsb1)/pow3(-mmgl + mmsb1) - (DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2))/
     pow2(-mmgl + mmsb2) + (mmsb2*mmst1)/pow3(-mmgl + mmsb2) - (mmsb2*mmt)/pow3
     (-mmgl + mmsb2) - pow2(mmsb2)/pow3(-mmgl + mmsb2))*pow4(g3))/3. + (4*Fin3(
     mmt,mmst2,mmgl,mmu)*(-3/(4.*(-mmgl + mmsb1)) - 3/(4.*(-mmgl + mmsb2)) - 2*
     mmgl*DeltaInv(mmt,mmst2,mmgl) - 2*mmsb1*DeltaInv(mmt,mmst2,mmgl) - 2*mmsb2
     *DeltaInv(mmt,mmst2,mmgl) + 4*mmst2*DeltaInv(mmt,mmst2,mmgl) - (4*mmsb1*
     mmst2*DeltaInv(mmt,mmst2,mmgl))/(-mmgl + mmsb1) - (4*mmsb2*mmst2*DeltaInv(
     mmt,mmst2,mmgl))/(-mmgl + mmsb2) + 2*mmt*DeltaInv(mmt,mmst2,mmgl) - (2*
     mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl))/(-mmgl + mmsb1) - (2*mmsb2*mmt*
     DeltaInv(mmt,mmst2,mmgl))/(-mmgl + mmsb2) - (mmst2*mmt*DeltaInv(mmt,mmst2,
     mmgl))/(-mmgl + mmsb1) - (mmst2*mmt*DeltaInv(mmt,mmst2,mmgl))/(-mmgl +
     mmsb2) + (3*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/(-mmgl + mmsb1) + (7*
     mmsb1)/(4.*pow2(-mmgl + mmsb1)) - (3*mmst2)/(4.*pow2(-mmgl + mmsb1)) + (3*
     mmt)/(4.*pow2(-mmgl + mmsb1)) + (mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl))
     /pow2(-mmgl + mmsb1) + (2*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/pow2
     (-mmgl + mmsb1) + (mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/pow2(-mmgl +
     mmsb1) + (3*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/(-mmgl + mmsb2) + (7*
     mmsb2)/(4.*pow2(-mmgl + mmsb2)) - (3*mmst2)/(4.*pow2(-mmgl + mmsb2)) + (3*
     mmt)/(4.*pow2(-mmgl + mmsb2)) + (mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl))
     /pow2(-mmgl + mmsb2) + (2*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/pow2
     (-mmgl + mmsb2) + (mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/pow2(-mmgl +
     mmsb2) + (DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(-mmgl + mmsb1) + (
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(-mmgl + mmsb2) - (mmsb1*DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmst2))/pow2(-mmgl + mmsb1) - (mmsb2*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmst2))/pow2(-mmgl + mmsb2) - (DeltaInv(mmt,mmst2,mmgl)*
     pow3(mmsb1))/pow2(-mmgl + mmsb1) + (mmsb1*mmst2)/pow3(-mmgl + mmsb1) - (
     mmsb1*mmt)/pow3(-mmgl + mmsb1) - pow2(mmsb1)/pow3(-mmgl + mmsb1) - (
     DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2))/pow2(-mmgl + mmsb2) + (mmsb2*mmst2)/
     pow3(-mmgl + mmsb2) - (mmsb2*mmt)/pow3(-mmgl + mmsb2) - pow2(mmsb2)/pow3(-
     mmgl + mmsb2))*pow4(g3))/3. + (4*Fin3(mmt,mmsb2,mmst1,mmu)*(-((mmsb2*mmst1
     *DeltaInv(mmt,mmsb2,mmst1))/(-mmgl + mmsb2)) - (mmsb2*mmt*DeltaInv(mmt,
     mmsb2,mmst1))/(2.*(-mmgl + mmsb2)) - (mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1))
     /(2.*(-mmgl + mmsb2)) + (DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2))/(2.*(-mmgl
      + mmsb2)) - (3*mmsb2)/(4.*pow2(-mmgl + mmsb2)) + (3*mmst1)/(4.*pow2(-mmgl
      + mmsb2)) - (3*mmt)/(4.*pow2(-mmgl + mmsb2)) + (mmsb2*mmst1*mmt*DeltaInv(
     mmt,mmsb2,mmst1))/pow2(-mmgl + mmsb2) + (2*mmst1*DeltaInv(mmt,mmsb2,mmst1)
     *pow2(mmsb2))/pow2(-mmgl + mmsb2) + (mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(
     mmsb2))/pow2(-mmgl + mmsb2) + (DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/(2.*
     (-mmgl + mmsb2)) - (mmsb2*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/pow2(-
     mmgl + mmsb2) - (DeltaInv(mmt,mmsb2,mmst1)*pow3(mmsb2))/pow2(-mmgl + mmsb2
     ) - (mmsb2*mmst1)/pow3(-mmgl + mmsb2) + (mmsb2*mmt)/pow3(-mmgl + mmsb2) +
     pow2(mmsb2)/pow3(-mmgl + mmsb2))*pow4(g3))/3. + (4*Fin3(mmt,mmsb2,mmst2,
     mmu)*(-((mmsb2*mmst2*DeltaInv(mmt,mmsb2,mmst2))/(-mmgl + mmsb2)) - (mmsb2*
     mmt*DeltaInv(mmt,mmsb2,mmst2))/(2.*(-mmgl + mmsb2)) - (mmst2*mmt*DeltaInv(
     mmt,mmsb2,mmst2))/(2.*(-mmgl + mmsb2)) + (DeltaInv(mmt,mmsb2,mmst2)*pow2(
     mmsb2))/(2.*(-mmgl + mmsb2)) - (3*mmsb2)/(4.*pow2(-mmgl + mmsb2)) + (3*
     mmst2)/(4.*pow2(-mmgl + mmsb2)) - (3*mmt)/(4.*pow2(-mmgl + mmsb2)) + (
     mmsb2*mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2))/pow2(-mmgl + mmsb2) + (2*mmst2*
     DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2))/pow2(-mmgl + mmsb2) + (mmt*DeltaInv
     (mmt,mmsb2,mmst2)*pow2(mmsb2))/pow2(-mmgl + mmsb2) + (DeltaInv(mmt,mmsb2,
     mmst2)*pow2(mmst2))/(2.*(-mmgl + mmsb2)) - (mmsb2*DeltaInv(mmt,mmsb2,mmst2
     )*pow2(mmst2))/pow2(-mmgl + mmsb2) - (DeltaInv(mmt,mmsb2,mmst2)*pow3(mmsb2
     ))/pow2(-mmgl + mmsb2) - (mmsb2*mmst2)/pow3(-mmgl + mmsb2) + (mmsb2*mmt)/
     pow3(-mmgl + mmsb2) + pow2(mmsb2)/pow3(-mmgl + mmsb2))*pow4(g3))/3. + (4*
     Fin20(mmsb2,mmsusy,mmu)*(4/(-mmgl + mmsb2) - (14*mmsb2)/pow2(-mmgl + mmsb2
     ) + (6*mmsusy)/pow2(-mmgl + mmsb2) - (8*mmsb2*mmsusy)/pow3(-mmgl + mmsb2)
     + (8*pow2(mmsb2))/pow3(-mmgl + mmsb2))*pow4(g3))/3. + (4*lb1u*lb2u*((3*
     mmsb1*mmsb2)/(4.*pow2(-mmgl + mmsb1)) + (3*mmsb1*mmsb2)/(4.*pow2(-mmgl +
     mmsb2)) - (mmsb2*pow2(mmsb1))/pow3(-mmgl + mmsb1) - (mmsb1*pow2(mmsb2))/
     pow3(-mmgl + mmsb2))*pow4(g3))/3. + (4*lb2u*lt1u*((mmsb2*mmst1*mmt*
     DeltaInv(mmt,mmsb2,mmst1))/(-mmgl + mmsb2) + (3*mmsb2*mmst1)/(4.*pow2(-
     mmgl + mmsb2)) - (2*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2))/pow2(
     -mmgl + mmsb2) - (mmst1*pow2(mmsb2))/pow3(-mmgl + mmsb2))*pow4(g3))/3. + (
     4*lgu*lt1u*((3*mmst1)/(4.*(-mmgl + mmsb1)) + (3*mmst1)/(4.*(-mmgl + mmsb2)
     ) - 4*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl) + (4*mmsb1*mmst1*mmt*DeltaInv(mmt
     ,mmst1,mmgl))/(-mmgl + mmsb1) + (4*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl
     ))/(-mmgl + mmsb2) - (7*mmsb1*mmst1)/(4.*pow2(-mmgl + mmsb1)) - (2*mmst1*
     mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1))/pow2(-mmgl + mmsb1) - (7*mmsb2*
     mmst1)/(4.*pow2(-mmgl + mmsb2)) - (2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmsb2))/pow2(-mmgl + mmsb2) + (mmst1*pow2(mmsb1))/pow3(-mmgl + mmsb1)
     + (mmst1*pow2(mmsb2))/pow3(-mmgl + mmsb2))*pow4(g3))/3. + (4*lb2u*lt2u*((
     mmsb2*mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2))/(-mmgl + mmsb2) + (3*mmsb2*
     mmst2)/(4.*pow2(-mmgl + mmsb2)) - (2*mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2)*
     pow2(mmsb2))/pow2(-mmgl + mmsb2) - (mmst2*pow2(mmsb2))/pow3(-mmgl + mmsb2)
     )*pow4(g3))/3. + (4*lgu*lt2u*((3*mmst2)/(4.*(-mmgl + mmsb1)) + (3*mmst2)/(
     4.*(-mmgl + mmsb2)) - 4*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl) + (4*mmsb1*
     mmst2*mmt*DeltaInv(mmt,mmst2,mmgl))/(-mmgl + mmsb1) + (4*mmsb2*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl))/(-mmgl + mmsb2) - (7*mmsb1*mmst2)/(4.*pow2(-mmgl
      + mmsb1)) - (2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/pow2(-mmgl
      + mmsb1) - (7*mmsb2*mmst2)/(4.*pow2(-mmgl + mmsb2)) - (2*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/pow2(-mmgl + mmsb2) + (mmst2*pow2(
     mmsb1))/pow3(-mmgl + mmsb1) + (mmst2*pow2(mmsb2))/pow3(-mmgl + mmsb2))*
     pow4(g3))/3. + (4*lb2u*lsu*((6*mmsb2*mmsusy)/pow2(-mmgl + mmsb2) - (8*
     mmsusy*pow2(mmsb2))/pow3(-mmgl + mmsb2))*pow4(g3))/3. + (4*lgu*lsu*((6*
     mmsusy)/(-mmgl + mmsb1) + (6*mmsusy)/(-mmgl + mmsb2) - (14*mmsb1*mmsusy)/
     pow2(-mmgl + mmsb1) - (14*mmsb2*mmsusy)/pow2(-mmgl + mmsb2) + (8*mmsusy*
     pow2(mmsb1))/pow3(-mmgl + mmsb1) + (8*mmsusy*pow2(mmsb2))/pow3(-mmgl +
     mmsb2))*pow4(g3))/3. + (4*mmsb1*mmsb2*Fin20(mmsb1,mmgl,mmu)*pow4(g3))/(3.*
     pow3(-mmgl + mmsb2)) - (4*mmsb1*mmsb2*Fin20(mmsb1,mmsb2,mmu)*pow4(g3))/(3.
     *pow3(-mmgl + mmsb2)) + (8*lb2u*mmsb1*pow2(mmsb2)*pow4(g3))/(3.*pow3(-mmgl
      + mmsb2)) - (8*lgu*mmsb1*pow2(mmsb2)*pow4(g3))/(3.*pow3(-mmgl + mmsb2)) +
     (4*lb1u*lgu*mmsb1*pow2(mmsb2)*pow4(g3))/(3.*pow3(-mmgl + mmsb2)) + (8*lb2u
     *mmst1*pow2(mmsb2)*pow4(g3))/(3.*pow3(-mmgl + mmsb2)) - (8*lgu*mmst1*pow2(
     mmsb2)*pow4(g3))/(3.*pow3(-mmgl + mmsb2)) + (8*lb2u*mmst2*pow2(mmsb2)*pow4
     (g3))/(3.*pow3(-mmgl + mmsb2)) - (8*lgu*mmst2*pow2(mmsb2)*pow4(g3))/(3.*
     pow3(-mmgl + mmsb2)) + (64*lb2u*mmsusy*pow2(mmsb2)*pow4(g3))/(3.*pow3(-
     mmgl + mmsb2)) - (64*lgu*mmsusy*pow2(mmsb2)*pow4(g3))/(3.*pow3(-mmgl +
     mmsb2)) - (16*lb2u*mmt*pow2(mmsb2)*pow4(g3))/(3.*pow3(-mmgl + mmsb2)) + (
     16*lgu*mmt*pow2(mmsb2)*pow4(g3))/(3.*pow3(-mmgl + mmsb2)) - (4*Fin20(mmsb1
     ,mmgl,mmu)*pow2(mmsb2)*pow4(g3))/(3.*pow3(-mmgl + mmsb2)) + (4*Fin20(mmsb1
     ,mmsb2,mmu)*pow2(mmsb2)*pow4(g3))/(3.*pow3(-mmgl + mmsb2)) - (16*Fin20(
     mmsb2,mmgl,mmu)*pow2(mmsb2)*pow4(g3))/(3.*pow3(-mmgl + mmsb2)) - (2*mmsb1*
     pow2(lb2u)*pow2(mmsb2)*pow4(g3))/(3.*pow3(-mmgl + mmsb2)) - (2*mmst1*pow2(
     lb2u)*pow2(mmsb2)*pow4(g3))/(3.*pow3(-mmgl + mmsb2)) - (2*mmst2*pow2(lb2u)
     *pow2(mmsb2)*pow4(g3))/(3.*pow3(-mmgl + mmsb2)) - (16*mmsusy*pow2(lb2u)*
     pow2(mmsb2)*pow4(g3))/(3.*pow3(-mmgl + mmsb2)) + (4*mmt*pow2(lb2u)*pow2(
     mmsb2)*pow4(g3))/(3.*pow3(-mmgl + mmsb2)) + (2*mmsb1*pow2(lgu)*pow2(mmsb2)
     *pow4(g3))/(3.*pow3(-mmgl + mmsb2)) + (2*mmst1*pow2(lgu)*pow2(mmsb2)*pow4(
     g3))/(3.*pow3(-mmgl + mmsb2)) + (2*mmst2*pow2(lgu)*pow2(mmsb2)*pow4(g3))/(
     3.*pow3(-mmgl + mmsb2)) + (16*mmsusy*pow2(lgu)*pow2(mmsb2)*pow4(g3))/(3.*
     pow3(-mmgl + mmsb2)) - (4*mmt*pow2(lgu)*pow2(mmsb2)*pow4(g3))/(3.*pow3(-
     mmgl + mmsb2)) - (112*pow3(mmsb2)*pow4(g3))/(3.*pow3(-mmgl + mmsb2)) - (
     220*lb2u*pow3(mmsb2)*pow4(g3))/(9.*pow3(-mmgl + mmsb2)) + (508*lgu*pow3(
     mmsb2)*pow4(g3))/(9.*pow3(-mmgl + mmsb2)) + (116*lb2u*lgu*pow3(mmsb2)*pow4
     (g3))/(9.*pow3(-mmgl + mmsb2)) - (16*zt2*pow3(mmsb2)*pow4(g3))/(3.*pow3(-
     mmgl + mmsb2)) + (10*pow2(lb2u)*pow3(mmsb2)*pow4(g3))/(3.*pow3(-mmgl +
     mmsb2)) - (242*pow2(lgu)*pow3(mmsb2)*pow4(g3))/(9.*pow3(-mmgl + mmsb2)) +
     (14*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmst1)*pow4(g3))/(3.*(-mmgl + mmsb1)) +
     (2*zt2*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmst1)*pow4(g3))/(3.*(-mmgl + mmsb1)
     ) + (14*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmst1)*pow4(g3))/(3.*(-mmgl + mmsb2
     )) + (2*zt2*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmst1)*pow4(g3))/(3.*(-mmgl +
     mmsb2)) + (28*DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1)*pow4(g3))/(3.*(-mmgl +
     mmsb1)) + (28*DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1)*pow4(g3))/(3.*(-mmgl +
     mmsb2)) + (4*zt2*DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1)*pow4(g3))/(3.*(-mmgl
      + mmsb1)) + (4*zt2*DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1)*pow4(g3))/(3.*(-
     mmgl + mmsb2)) - (28*mmsb1*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmst1)*pow4(g3))
     /(3.*pow2(-mmgl + mmsb1)) - (4*mmsb1*zt2*DeltaInv(mmt,mmsb1,mmst1)*pow3(
     mmst1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) - (28*mmsb1*DeltaInv(mmt,mmst1,
     mmgl)*pow3(mmst1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) - (4*mmsb1*zt2*
     DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) -
     (28*mmsb2*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmst1)*pow4(g3))/(3.*pow2(-mmgl +
     mmsb2)) - (4*mmsb2*zt2*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmst1)*pow4(g3))/(3.
     *pow2(-mmgl + mmsb2)) - (28*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1)*
     pow4(g3))/(3.*pow2(-mmgl + mmsb2)) - (4*mmsb2*zt2*DeltaInv(mmt,mmst1,mmgl)
     *pow3(mmst1)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) + (4*lt1u*ltu*((-3*mmst1)/
     (2.*(-mmgl + mmsb1)) - (3*mmst1)/(2.*(-mmgl + mmsb2)) - (3*mmsb1*mmst1*mmt
     *DeltaInv(mmt,mmsb1,mmst1))/(2.*(-mmgl + mmsb1)) - (3*mmsb2*mmst1*mmt*
     DeltaInv(mmt,mmsb2,mmst1))/(2.*(-mmgl + mmsb2)) - 2*mmgl*mmst1*DeltaInv(
     mmt,mmst1,mmgl) - 2*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl) - 2*mmsb2*mmst1*
     DeltaInv(mmt,mmst1,mmgl) + 6*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl) - (6*mmsb1
     *mmst1*mmt*DeltaInv(mmt,mmst1,mmgl))/(-mmgl + mmsb1) - (6*mmsb2*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl))/(-mmgl + mmsb2) + (mmst1*DeltaInv(mmt,mmsb1,
     mmst1)*pow2(mmsb1))/(2.*(-mmgl + mmsb1)) + (3*mmst1*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmsb1))/(-mmgl + mmsb1) + (2*mmsb1*mmst1)/pow2(-mmgl + mmsb1) +
     (3*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1))/pow2(-mmgl + mmsb1) +
     (3*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1))/pow2(-mmgl + mmsb1) + (
     mmst1*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2))/(2.*(-mmgl + mmsb2)) + (3*
     mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2))/(-mmgl + mmsb2) + (2*mmsb2*
     mmst1)/pow2(-mmgl + mmsb2) + (3*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(
     mmsb2))/pow2(-mmgl + mmsb2) + (3*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmsb2))/pow2(-mmgl + mmsb2) - (mmsb1*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1)
     )/(-mmgl + mmsb1) - (mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1))/(2.*(-mmgl
      + mmsb1)) - (mmsb2*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/(-mmgl + mmsb2)
     - (mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/(2.*(-mmgl + mmsb2)) + 4*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1) - (4*mmsb1*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmst1))/(-mmgl + mmsb1) - (4*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmst1))/(-mmgl + mmsb2) - (mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(-
     mmgl + mmsb1) - (mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(-mmgl + mmsb2)
     + (mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1))/pow2(-mmgl + mmsb1) +
     (mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/pow2(-mmgl + mmsb1) + (2*
     DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow2(mmst1))/pow2(-mmgl + mmsb1) + (
     2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(mmst1))/pow2(-mmgl + mmsb1) +
     (mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/pow2(-mmgl + mmsb2) + (
     mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/pow2(-mmgl + mmsb2) + (2*
     DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2)*pow2(mmst1))/pow2(-mmgl + mmsb2) + (
     2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(mmst1))/pow2(-mmgl + mmsb2) -
     (mmst1*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1))/pow2(-mmgl + mmsb1) - (mmst1
     *DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1))/pow2(-mmgl + mmsb1) - (mmst1*
     DeltaInv(mmt,mmsb2,mmst1)*pow3(mmsb2))/pow2(-mmgl + mmsb2) - (mmst1*
     DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2))/pow2(-mmgl + mmsb2) + (DeltaInv(mmt,
     mmsb1,mmst1)*pow3(mmst1))/(2.*(-mmgl + mmsb1)) + (DeltaInv(mmt,mmsb2,mmst1
     )*pow3(mmst1))/(2.*(-mmgl + mmsb2)) + (DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1
     ))/(-mmgl + mmsb1) + (DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1))/(-mmgl + mmsb2
     ) - (mmsb1*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmst1))/pow2(-mmgl + mmsb1) - (
     mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1))/pow2(-mmgl + mmsb1) - (mmsb2*
     DeltaInv(mmt,mmsb2,mmst1)*pow3(mmst1))/pow2(-mmgl + mmsb2) - (mmsb2*
     DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1))/pow2(-mmgl + mmsb2))*pow4(g3))/3. +
     (4*pow2(lt1u)*(0.25 - (3*mmst1)/(8.*(-mmgl + mmsb1)) - (3*mmst1)/(8.*(-
     mmgl + mmsb2)) - (mmsb1*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1))/(4.*(-mmgl +
     mmsb1)) - (mmsb2*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1))/(4.*(-mmgl + mmsb2))
     - mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl) - mmsb1*mmst1*DeltaInv(mmt,mmst1,
     mmgl) - mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl) + mmst1*mmt*DeltaInv(mmt,
     mmst1,mmgl) - (mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl))/(-mmgl + mmsb1) -
     (mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl))/(-mmgl + mmsb2) + (mmst1*
     DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1))/(4.*(-mmgl + mmsb1)) + (3*mmst1*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1))/(2.*(-mmgl + mmsb1)) + (mmsb1*mmst1)
     /(2.*pow2(-mmgl + mmsb1)) + (mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(
     mmsb1))/(2.*pow2(-mmgl + mmsb1)) + (mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmsb1))/(2.*pow2(-mmgl + mmsb1)) + (mmst1*DeltaInv(mmt,mmsb2,mmst1)*
     pow2(mmsb2))/(4.*(-mmgl + mmsb2)) + (3*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2
     (mmsb2))/(2.*(-mmgl + mmsb2)) + (mmsb2*mmst1)/(2.*pow2(-mmgl + mmsb2)) + (
     mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2))/(2.*pow2(-mmgl + mmsb2))
     + (mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2))/(2.*pow2(-mmgl + mmsb2)
     ) - (mmsb1*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1))/(2.*(-mmgl + mmsb1)) - (
     mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1))/(4.*(-mmgl + mmsb1)) - (mmsb2*
     DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/(2.*(-mmgl + mmsb2)) - (mmt*
     DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/(4.*(-mmgl + mmsb2)) + 2*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmst1) - (2*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1
     ))/(-mmgl + mmsb1) - (2*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(-mmgl
      + mmsb2) - (mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(2.*(-mmgl + mmsb1)
     ) - (mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(2.*(-mmgl + mmsb2)) + (
     mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1))/(2.*pow2(-mmgl + mmsb1))
     + (mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(2.*pow2(-mmgl + mmsb1)
     ) + (DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow2(mmst1))/pow2(-mmgl + mmsb1
     ) + (DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(mmst1))/pow2(-mmgl + mmsb1)
     + (mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/(2.*pow2(-mmgl + mmsb2
     )) + (mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(2.*pow2(-mmgl +
     mmsb2)) + (DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2)*pow2(mmst1))/pow2(-mmgl +
     mmsb2) + (DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(mmst1))/pow2(-mmgl +
     mmsb2) - (mmst1*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1))/(2.*pow2(-mmgl +
     mmsb1)) - (mmst1*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1))/(2.*pow2(-mmgl +
     mmsb1)) - (mmst1*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmsb2))/(2.*pow2(-mmgl +
     mmsb2)) - (mmst1*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2))/(2.*pow2(-mmgl +
     mmsb2)) + (DeltaInv(mmt,mmsb1,mmst1)*pow3(mmst1))/(4.*(-mmgl + mmsb1)) + (
     DeltaInv(mmt,mmsb2,mmst1)*pow3(mmst1))/(4.*(-mmgl + mmsb2)) + (DeltaInv(
     mmt,mmst1,mmgl)*pow3(mmst1))/(2.*(-mmgl + mmsb1)) + (DeltaInv(mmt,mmst1,
     mmgl)*pow3(mmst1))/(2.*(-mmgl + mmsb2)) - (mmsb1*DeltaInv(mmt,mmsb1,mmst1)
     *pow3(mmst1))/(2.*pow2(-mmgl + mmsb1)) - (mmsb1*DeltaInv(mmt,mmst1,mmgl)*
     pow3(mmst1))/(2.*pow2(-mmgl + mmsb1)) - (mmsb2*DeltaInv(mmt,mmsb2,mmst1)*
     pow3(mmst1))/(2.*pow2(-mmgl + mmsb2)) - (mmsb2*DeltaInv(mmt,mmst1,mmgl)*
     pow3(mmst1))/(2.*pow2(-mmgl + mmsb2)))*pow4(g3))/3. + (4*lt1u*(
     0.9166666666666666 + (3*mmst1)/(2.*(-mmgl + mmsb1)) + (3*mmst1)/(2.*(-mmgl
      + mmsb2)) + (3*mmsb1*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1))/(2.*(-mmgl +
     mmsb1)) + (3*mmsb2*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1))/(2.*(-mmgl + mmsb2
     )) + 6*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl) + 6*mmsb1*mmst1*DeltaInv(mmt,
     mmst1,mmgl) + 6*mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl) - 6*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl) + (6*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl))/(-
     mmgl + mmsb1) + (6*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl))/(-mmgl +
     mmsb2) - (3*mmst1*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1))/(2.*(-mmgl +
     mmsb1)) - (9*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1))/(-mmgl + mmsb1) -
     (2*mmsb1*mmst1)/pow2(-mmgl + mmsb1) - (3*mmst1*mmt*DeltaInv(mmt,mmsb1,
     mmst1)*pow2(mmsb1))/pow2(-mmgl + mmsb1) - (3*mmst1*mmt*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmsb1))/pow2(-mmgl + mmsb1) - (3*mmst1*DeltaInv(mmt,mmsb2,mmst1
     )*pow2(mmsb2))/(2.*(-mmgl + mmsb2)) - (9*mmst1*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmsb2))/(-mmgl + mmsb2) - (2*mmsb2*mmst1)/pow2(-mmgl + mmsb2) - (3*
     mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2))/pow2(-mmgl + mmsb2) - (3*
     mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2))/pow2(-mmgl + mmsb2) + (3*
     mmsb1*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1))/(-mmgl + mmsb1) + (3*mmt*
     DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1))/(2.*(-mmgl + mmsb1)) + (3*mmsb2*
     DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/(-mmgl + mmsb2) + (3*mmt*DeltaInv(
     mmt,mmsb2,mmst1)*pow2(mmst1))/(2.*(-mmgl + mmsb2)) - 12*DeltaInv(mmt,mmst1
     ,mmgl)*pow2(mmst1) + (12*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(-
     mmgl + mmsb1) + (12*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(-mmgl +
     mmsb2) + (3*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(-mmgl + mmsb1) + (3
     *mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(-mmgl + mmsb2) - (3*mmsb1*mmt*
     DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1))/pow2(-mmgl + mmsb1) - (3*mmsb1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/pow2(-mmgl + mmsb1) - (6*DeltaInv(
     mmt,mmsb1,mmst1)*pow2(mmsb1)*pow2(mmst1))/pow2(-mmgl + mmsb1) - (6*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(mmst1))/pow2(-mmgl + mmsb1) - (3
     *mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/pow2(-mmgl + mmsb2) - (3
     *mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/pow2(-mmgl + mmsb2) - (6*
     DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2)*pow2(mmst1))/pow2(-mmgl + mmsb2) - (
     6*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(mmst1))/pow2(-mmgl + mmsb2) +
     (3*mmst1*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1))/pow2(-mmgl + mmsb1) + (3*
     mmst1*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1))/pow2(-mmgl + mmsb1) + (3*mmst1
     *DeltaInv(mmt,mmsb2,mmst1)*pow3(mmsb2))/pow2(-mmgl + mmsb2) + (3*mmst1*
     DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2))/pow2(-mmgl + mmsb2) - (3*DeltaInv(
     mmt,mmsb1,mmst1)*pow3(mmst1))/(2.*(-mmgl + mmsb1)) - (3*DeltaInv(mmt,mmsb2
     ,mmst1)*pow3(mmst1))/(2.*(-mmgl + mmsb2)) - (3*DeltaInv(mmt,mmst1,mmgl)*
     pow3(mmst1))/(-mmgl + mmsb1) - (3*DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1))/(-
     mmgl + mmsb2) + (3*mmsb1*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmst1))/pow2(-mmgl
      + mmsb1) + (3*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1))/pow2(-mmgl +
     mmsb1) + (3*mmsb2*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmst1))/pow2(-mmgl +
     mmsb2) + (3*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1))/pow2(-mmgl + mmsb2
     ))*pow4(g3))/3. + (14*DeltaInv(mmt,mmsb1,mmst2)*pow3(mmst2)*pow4(g3))/(3.*
     (-mmgl + mmsb1)) + (2*zt2*DeltaInv(mmt,mmsb1,mmst2)*pow3(mmst2)*pow4(g3))/
     (3.*(-mmgl + mmsb1)) + (14*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmst2)*pow4(g3))
     /(3.*(-mmgl + mmsb2)) + (2*zt2*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmst2)*pow4(
     g3))/(3.*(-mmgl + mmsb2)) + (28*DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2)*pow4(
     g3))/(3.*(-mmgl + mmsb1)) + (28*DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2)*pow4(
     g3))/(3.*(-mmgl + mmsb2)) + (4*zt2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2)*
     pow4(g3))/(3.*(-mmgl + mmsb1)) + (4*zt2*DeltaInv(mmt,mmst2,mmgl)*pow3(
     mmst2)*pow4(g3))/(3.*(-mmgl + mmsb2)) - (28*mmsb1*DeltaInv(mmt,mmsb1,mmst2
     )*pow3(mmst2)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) - (4*mmsb1*zt2*DeltaInv(
     mmt,mmsb1,mmst2)*pow3(mmst2)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)) - (28*
     mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2)*pow4(g3))/(3.*pow2(-mmgl +
     mmsb1)) - (4*mmsb1*zt2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2)*pow4(g3))/(3.*
     pow2(-mmgl + mmsb1)) - (28*mmsb2*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmst2)*
     pow4(g3))/(3.*pow2(-mmgl + mmsb2)) - (4*mmsb2*zt2*DeltaInv(mmt,mmsb2,mmst2
     )*pow3(mmst2)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) - (28*mmsb2*DeltaInv(mmt,
     mmst2,mmgl)*pow3(mmst2)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) - (4*mmsb2*zt2*
     DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2)*pow4(g3))/(3.*pow2(-mmgl + mmsb2)) +
     (4*lt2u*ltu*((-3*mmst2)/(2.*(-mmgl + mmsb1)) - (3*mmst2)/(2.*(-mmgl +
     mmsb2)) - (3*mmsb1*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2))/(2.*(-mmgl + mmsb1
     )) - (3*mmsb2*mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2))/(2.*(-mmgl + mmsb2)) -
     2*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl) - 2*mmsb1*mmst2*DeltaInv(mmt,mmst2,
     mmgl) - 2*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl) + 6*mmst2*mmt*DeltaInv(mmt,
     mmst2,mmgl) - (6*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl))/(-mmgl + mmsb1)
     - (6*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl))/(-mmgl + mmsb2) + (mmst2*
     DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))/(2.*(-mmgl + mmsb1)) + (3*mmst2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/(-mmgl + mmsb1) + (2*mmsb1*mmst2)/
     pow2(-mmgl + mmsb1) + (3*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))/
     pow2(-mmgl + mmsb1) + (3*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/
     pow2(-mmgl + mmsb1) + (mmst2*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2))/(2.*(-
     mmgl + mmsb2)) + (3*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/(-mmgl +
     mmsb2) + (2*mmsb2*mmst2)/pow2(-mmgl + mmsb2) + (3*mmst2*mmt*DeltaInv(mmt,
     mmsb2,mmst2)*pow2(mmsb2))/pow2(-mmgl + mmsb2) + (3*mmst2*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmsb2))/pow2(-mmgl + mmsb2) - (mmsb1*DeltaInv(mmt,mmsb1,
     mmst2)*pow2(mmst2))/(-mmgl + mmsb1) - (mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(
     mmst2))/(2.*(-mmgl + mmsb1)) - (mmsb2*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2
     ))/(-mmgl + mmsb2) - (mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2))/(2.*(-
     mmgl + mmsb2)) + 4*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2) - (4*mmsb1*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(-mmgl + mmsb1) - (4*mmsb2*DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmst2))/(-mmgl + mmsb2) - (mmt*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmst2))/(-mmgl + mmsb1) - (mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmst2))/(-mmgl + mmsb2) + (mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2)
     )/pow2(-mmgl + mmsb1) + (mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/
     pow2(-mmgl + mmsb1) + (2*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*pow2(mmst2)
     )/pow2(-mmgl + mmsb1) + (2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(mmst2
     ))/pow2(-mmgl + mmsb1) + (mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2))
     /pow2(-mmgl + mmsb2) + (mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/
     pow2(-mmgl + mmsb2) + (2*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow2(mmst2)
     )/pow2(-mmgl + mmsb2) + (2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(mmst2
     ))/pow2(-mmgl + mmsb2) - (mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1))/
     pow2(-mmgl + mmsb1) - (mmst2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1))/pow2(-
     mmgl + mmsb1) - (mmst2*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmsb2))/pow2(-mmgl +
     mmsb2) - (mmst2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2))/pow2(-mmgl + mmsb2)
     + (DeltaInv(mmt,mmsb1,mmst2)*pow3(mmst2))/(2.*(-mmgl + mmsb1)) + (DeltaInv
     (mmt,mmsb2,mmst2)*pow3(mmst2))/(2.*(-mmgl + mmsb2)) + (DeltaInv(mmt,mmst2,
     mmgl)*pow3(mmst2))/(-mmgl + mmsb1) + (DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2)
     )/(-mmgl + mmsb2) - (mmsb1*DeltaInv(mmt,mmsb1,mmst2)*pow3(mmst2))/pow2(-
     mmgl + mmsb1) - (mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2))/pow2(-mmgl +
     mmsb1) - (mmsb2*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmst2))/pow2(-mmgl + mmsb2)
     - (mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2))/pow2(-mmgl + mmsb2))*pow4(
     g3))/3. + (4*pow2(lt2u)*(0.25 - (3*mmst2)/(8.*(-mmgl + mmsb1)) - (3*mmst2)
     /(8.*(-mmgl + mmsb2)) - (mmsb1*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2))/(4.*(-
     mmgl + mmsb1)) - (mmsb2*mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2))/(4.*(-mmgl +
     mmsb2)) - mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl) - mmsb1*mmst2*DeltaInv(mmt,
     mmst2,mmgl) - mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl) + mmst2*mmt*DeltaInv(
     mmt,mmst2,mmgl) - (mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl))/(-mmgl +
     mmsb1) - (mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl))/(-mmgl + mmsb2) + (
     mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))/(4.*(-mmgl + mmsb1)) + (3*
     mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/(2.*(-mmgl + mmsb1)) + (mmsb1*
     mmst2)/(2.*pow2(-mmgl + mmsb1)) + (mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2)*
     pow2(mmsb1))/(2.*pow2(-mmgl + mmsb1)) + (mmst2*mmt*DeltaInv(mmt,mmst2,mmgl
     )*pow2(mmsb1))/(2.*pow2(-mmgl + mmsb1)) + (mmst2*DeltaInv(mmt,mmsb2,mmst2)
     *pow2(mmsb2))/(4.*(-mmgl + mmsb2)) + (3*mmst2*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmsb2))/(2.*(-mmgl + mmsb2)) + (mmsb2*mmst2)/(2.*pow2(-mmgl + mmsb2))
     + (mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2))/(2.*pow2(-mmgl + mmsb2
     )) + (mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/(2.*pow2(-mmgl +
     mmsb2)) - (mmsb1*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2))/(2.*(-mmgl + mmsb1
     )) - (mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2))/(4.*(-mmgl + mmsb1)) - (
     mmsb2*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2))/(2.*(-mmgl + mmsb2)) - (mmt*
     DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2))/(4.*(-mmgl + mmsb2)) + 2*DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmst2) - (2*mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2
     ))/(-mmgl + mmsb1) - (2*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(-mmgl
      + mmsb2) - (mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(2.*(-mmgl + mmsb1)
     ) - (mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(2.*(-mmgl + mmsb2)) + (
     mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2))/(2.*pow2(-mmgl + mmsb1))
     + (mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(2.*pow2(-mmgl + mmsb1)
     ) + (DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*pow2(mmst2))/pow2(-mmgl + mmsb1
     ) + (DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(mmst2))/pow2(-mmgl + mmsb1)
     + (mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2))/(2.*pow2(-mmgl + mmsb2
     )) + (mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(2.*pow2(-mmgl +
     mmsb2)) + (DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow2(mmst2))/pow2(-mmgl +
     mmsb2) + (DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(mmst2))/pow2(-mmgl +
     mmsb2) - (mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1))/(2.*pow2(-mmgl +
     mmsb1)) - (mmst2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1))/(2.*pow2(-mmgl +
     mmsb1)) - (mmst2*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmsb2))/(2.*pow2(-mmgl +
     mmsb2)) - (mmst2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2))/(2.*pow2(-mmgl +
     mmsb2)) + (DeltaInv(mmt,mmsb1,mmst2)*pow3(mmst2))/(4.*(-mmgl + mmsb1)) + (
     DeltaInv(mmt,mmsb2,mmst2)*pow3(mmst2))/(4.*(-mmgl + mmsb2)) + (DeltaInv(
     mmt,mmst2,mmgl)*pow3(mmst2))/(2.*(-mmgl + mmsb1)) + (DeltaInv(mmt,mmst2,
     mmgl)*pow3(mmst2))/(2.*(-mmgl + mmsb2)) - (mmsb1*DeltaInv(mmt,mmsb1,mmst2)
     *pow3(mmst2))/(2.*pow2(-mmgl + mmsb1)) - (mmsb1*DeltaInv(mmt,mmst2,mmgl)*
     pow3(mmst2))/(2.*pow2(-mmgl + mmsb1)) - (mmsb2*DeltaInv(mmt,mmsb2,mmst2)*
     pow3(mmst2))/(2.*pow2(-mmgl + mmsb2)) - (mmsb2*DeltaInv(mmt,mmst2,mmgl)*
     pow3(mmst2))/(2.*pow2(-mmgl + mmsb2)))*pow4(g3))/3. + (4*lt2u*(
     0.9166666666666666 + (3*mmst2)/(2.*(-mmgl + mmsb1)) + (3*mmst2)/(2.*(-mmgl
      + mmsb2)) + (3*mmsb1*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2))/(2.*(-mmgl +
     mmsb1)) + (3*mmsb2*mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2))/(2.*(-mmgl + mmsb2
     )) + 6*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl) + 6*mmsb1*mmst2*DeltaInv(mmt,
     mmst2,mmgl) + 6*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl) - 6*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl) + (6*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl))/(-
     mmgl + mmsb1) + (6*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl))/(-mmgl +
     mmsb2) - (3*mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))/(2.*(-mmgl +
     mmsb1)) - (9*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/(-mmgl + mmsb1) -
     (2*mmsb1*mmst2)/pow2(-mmgl + mmsb1) - (3*mmst2*mmt*DeltaInv(mmt,mmsb1,
     mmst2)*pow2(mmsb1))/pow2(-mmgl + mmsb1) - (3*mmst2*mmt*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmsb1))/pow2(-mmgl + mmsb1) - (3*mmst2*DeltaInv(mmt,mmsb2,mmst2
     )*pow2(mmsb2))/(2.*(-mmgl + mmsb2)) - (9*mmst2*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmsb2))/(-mmgl + mmsb2) - (2*mmsb2*mmst2)/pow2(-mmgl + mmsb2) - (3*
     mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2))/pow2(-mmgl + mmsb2) - (3*
     mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/pow2(-mmgl + mmsb2) + (3*
     mmsb1*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2))/(-mmgl + mmsb1) + (3*mmt*
     DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2))/(2.*(-mmgl + mmsb1)) + (3*mmsb2*
     DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2))/(-mmgl + mmsb2) + (3*mmt*DeltaInv(
     mmt,mmsb2,mmst2)*pow2(mmst2))/(2.*(-mmgl + mmsb2)) - 12*DeltaInv(mmt,mmst2
     ,mmgl)*pow2(mmst2) + (12*mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(-
     mmgl + mmsb1) + (12*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(-mmgl +
     mmsb2) + (3*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(-mmgl + mmsb1) + (3
     *mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(-mmgl + mmsb2) - (3*mmsb1*mmt*
     DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2))/pow2(-mmgl + mmsb1) - (3*mmsb1*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/pow2(-mmgl + mmsb1) - (6*DeltaInv(
     mmt,mmsb1,mmst2)*pow2(mmsb1)*pow2(mmst2))/pow2(-mmgl + mmsb1) - (6*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(mmst2))/pow2(-mmgl + mmsb1) - (3
     *mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2))/pow2(-mmgl + mmsb2) - (3
     *mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/pow2(-mmgl + mmsb2) - (6*
     DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow2(mmst2))/pow2(-mmgl + mmsb2) - (
     6*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(mmst2))/pow2(-mmgl + mmsb2) +
     (3*mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1))/pow2(-mmgl + mmsb1) + (3*
     mmst2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1))/pow2(-mmgl + mmsb1) + (3*mmst2
     *DeltaInv(mmt,mmsb2,mmst2)*pow3(mmsb2))/pow2(-mmgl + mmsb2) + (3*mmst2*
     DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2))/pow2(-mmgl + mmsb2) - (3*DeltaInv(
     mmt,mmsb1,mmst2)*pow3(mmst2))/(2.*(-mmgl + mmsb1)) - (3*DeltaInv(mmt,mmsb2
     ,mmst2)*pow3(mmst2))/(2.*(-mmgl + mmsb2)) - (3*DeltaInv(mmt,mmst2,mmgl)*
     pow3(mmst2))/(-mmgl + mmsb1) - (3*DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2))/(-
     mmgl + mmsb2) + (3*mmsb1*DeltaInv(mmt,mmsb1,mmst2)*pow3(mmst2))/pow2(-mmgl
      + mmsb1) + (3*mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2))/pow2(-mmgl +
     mmsb1) + (3*mmsb2*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmst2))/pow2(-mmgl +
     mmsb2) + (3*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2))/pow2(-mmgl + mmsb2
     ))*pow4(g3))/3. - (28*DeltaInv(mmt,mmsb1,mmst1)*pow4(g3)*pow4(mmsb1))/(3.*
     pow2(-mmgl + mmsb1)) + (4*lb1u*DeltaInv(mmt,mmsb1,mmst1)*pow4(g3)*pow4(
     mmsb1))/pow2(-mmgl + mmsb1) - (4*zt2*DeltaInv(mmt,mmsb1,mmst1)*pow4(g3)*
     pow4(mmsb1))/(3.*pow2(-mmgl + mmsb1)) - (28*DeltaInv(mmt,mmsb1,mmst2)*pow4
     (g3)*pow4(mmsb1))/(3.*pow2(-mmgl + mmsb1)) + (4*lb1u*DeltaInv(mmt,mmsb1,
     mmst2)*pow4(g3)*pow4(mmsb1))/pow2(-mmgl + mmsb1) - (4*zt2*DeltaInv(mmt,
     mmsb1,mmst2)*pow4(g3)*pow4(mmsb1))/(3.*pow2(-mmgl + mmsb1)) - (28*DeltaInv
     (mmt,mmst1,mmgl)*pow4(g3)*pow4(mmsb1))/(3.*pow2(-mmgl + mmsb1)) + (4*lgu*
     DeltaInv(mmt,mmst1,mmgl)*pow4(g3)*pow4(mmsb1))/pow2(-mmgl + mmsb1) - (4*
     zt2*DeltaInv(mmt,mmst1,mmgl)*pow4(g3)*pow4(mmsb1))/(3.*pow2(-mmgl + mmsb1)
     ) - (28*DeltaInv(mmt,mmst2,mmgl)*pow4(g3)*pow4(mmsb1))/(3.*pow2(-mmgl +
     mmsb1)) + (4*lgu*DeltaInv(mmt,mmst2,mmgl)*pow4(g3)*pow4(mmsb1))/pow2(-mmgl
      + mmsb1) - (4*zt2*DeltaInv(mmt,mmst2,mmgl)*pow4(g3)*pow4(mmsb1))/(3.*pow2
     (-mmgl + mmsb1)) - (2*DeltaInv(mmt,mmsb1,mmst1)*pow2(lb1u)*pow4(g3)*pow4(
     mmsb1))/(3.*pow2(-mmgl + mmsb1)) - (2*DeltaInv(mmt,mmsb1,mmst2)*pow2(lb1u)
     *pow4(g3)*pow4(mmsb1))/(3.*pow2(-mmgl + mmsb1)) - (2*DeltaInv(mmt,mmst1,
     mmgl)*pow2(lgu)*pow4(g3)*pow4(mmsb1))/(3.*pow2(-mmgl + mmsb1)) - (2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(lgu)*pow4(g3)*pow4(mmsb1))/(3.*pow2(-mmgl +
     mmsb1)) + (4*lb1u*ltu*pow4(g3)*(-(mmsb1/(-mmgl + mmsb1)) - (3*mmsb1*mmst1*
     mmt*DeltaInv(mmt,mmsb1,mmst1))/(2.*(-mmgl + mmsb1)) - (3*mmsb1*mmst2*mmt*
     DeltaInv(mmt,mmsb1,mmst2))/(2.*(-mmgl + mmsb1)) - (mmst1*DeltaInv(mmt,
     mmsb1,mmst1)*pow2(mmsb1))/(-mmgl + mmsb1) - (mmt*DeltaInv(mmt,mmsb1,mmst1)
     *pow2(mmsb1))/(2.*(-mmgl + mmsb1)) - (mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow2
     (mmsb1))/(-mmgl + mmsb1) - (mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))/(2.
     *(-mmgl + mmsb1)) - (3*mmsb1*mmt)/(2.*pow2(-mmgl + mmsb1)) + (2*pow2(mmsb1
     ))/pow2(-mmgl + mmsb1) + (3*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1
     ))/pow2(-mmgl + mmsb1) + (3*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1
     ))/pow2(-mmgl + mmsb1) + (mmsb1*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1))/(2.
     *(-mmgl + mmsb1)) - (DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow2(mmst1))/
     pow2(-mmgl + mmsb1) + (mmsb1*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2))/(2.*(-
     mmgl + mmsb1)) - (DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*pow2(mmst2))/pow2(
     -mmgl + mmsb1) + (DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1))/(2.*(-mmgl +
     mmsb1)) + (DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1))/(2.*(-mmgl + mmsb1)) + (
     2*mmst1*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1))/pow2(-mmgl + mmsb1) + (mmt*
     DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1))/pow2(-mmgl + mmsb1) + (2*mmst2*
     DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1))/pow2(-mmgl + mmsb1) + (mmt*DeltaInv
     (mmt,mmsb1,mmst2)*pow3(mmsb1))/pow2(-mmgl + mmsb1) + (2*mmt*pow2(mmsb1))/
     pow3(-mmgl + mmsb1) - (DeltaInv(mmt,mmsb1,mmst1)*pow4(mmsb1))/pow2(-mmgl +
     mmsb1) - (DeltaInv(mmt,mmsb1,mmst2)*pow4(mmsb1))/pow2(-mmgl + mmsb1)))/3.
     - (4*lb1u*lgu*pow4(g3)*pow4(mmsb1))/(3.*pow4(-mmgl + mmsb1)) + (2*pow2(
     lb1u)*pow4(g3)*pow4(mmsb1))/(3.*pow4(-mmgl + mmsb1)) + (2*pow2(lgu)*pow4(
     g3)*pow4(mmsb1))/(3.*pow4(-mmgl + mmsb1)) + (4*lb1u*lb2u*pow4(g3)*pow4(
     mmsb2))/(9.*pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (4*lb1u*lgu*pow4(g3
     )*pow4(mmsb2))/(9.*pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (4*lb2u*lgu*
     pow4(g3)*pow4(mmsb2))/(9.*pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) + (4*
     pow2(lgu)*pow4(g3)*pow4(mmsb2))/(9.*pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2
     )) - (28*DeltaInv(mmt,mmsb2,mmst1)*pow4(g3)*pow4(mmsb2))/(3.*pow2(-mmgl +
     mmsb2)) + (4*lb2u*DeltaInv(mmt,mmsb2,mmst1)*pow4(g3)*pow4(mmsb2))/pow2(-
     mmgl + mmsb2) - (4*zt2*DeltaInv(mmt,mmsb2,mmst1)*pow4(g3)*pow4(mmsb2))/(3.
     *pow2(-mmgl + mmsb2)) - (28*DeltaInv(mmt,mmsb2,mmst2)*pow4(g3)*pow4(mmsb2)
     )/(3.*pow2(-mmgl + mmsb2)) + (4*lb2u*DeltaInv(mmt,mmsb2,mmst2)*pow4(g3)*
     pow4(mmsb2))/pow2(-mmgl + mmsb2) - (4*zt2*DeltaInv(mmt,mmsb2,mmst2)*pow4(
     g3)*pow4(mmsb2))/(3.*pow2(-mmgl + mmsb2)) - (28*DeltaInv(mmt,mmst1,mmgl)*
     pow4(g3)*pow4(mmsb2))/(3.*pow2(-mmgl + mmsb2)) + (4*lgu*DeltaInv(mmt,mmst1
     ,mmgl)*pow4(g3)*pow4(mmsb2))/pow2(-mmgl + mmsb2) - (4*zt2*DeltaInv(mmt,
     mmst1,mmgl)*pow4(g3)*pow4(mmsb2))/(3.*pow2(-mmgl + mmsb2)) - (28*DeltaInv(
     mmt,mmst2,mmgl)*pow4(g3)*pow4(mmsb2))/(3.*pow2(-mmgl + mmsb2)) + (4*lgu*
     DeltaInv(mmt,mmst2,mmgl)*pow4(g3)*pow4(mmsb2))/pow2(-mmgl + mmsb2) - (4*
     zt2*DeltaInv(mmt,mmst2,mmgl)*pow4(g3)*pow4(mmsb2))/(3.*pow2(-mmgl + mmsb2)
     ) - (2*DeltaInv(mmt,mmsb2,mmst1)*pow2(lb2u)*pow4(g3)*pow4(mmsb2))/(3.*pow2
     (-mmgl + mmsb2)) - (2*DeltaInv(mmt,mmsb2,mmst2)*pow2(lb2u)*pow4(g3)*pow4(
     mmsb2))/(3.*pow2(-mmgl + mmsb2)) - (2*DeltaInv(mmt,mmst1,mmgl)*pow2(lgu)*
     pow4(g3)*pow4(mmsb2))/(3.*pow2(-mmgl + mmsb2)) - (2*DeltaInv(mmt,mmst2,
     mmgl)*pow2(lgu)*pow4(g3)*pow4(mmsb2))/(3.*pow2(-mmgl + mmsb2)) + (4*lb1u*
     lb2u*pow4(g3)*pow4(mmsb2))/(9.*pow2(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) -
     (4*lb1u*lgu*pow4(g3)*pow4(mmsb2))/(9.*pow2(mmsb1 - mmsb2)*pow2(-mmgl +
     mmsb2)) - (4*lb2u*lgu*pow4(g3)*pow4(mmsb2))/(9.*pow2(mmsb1 - mmsb2)*pow2(-
     mmgl + mmsb2)) + (4*pow2(lgu)*pow4(g3)*pow4(mmsb2))/(9.*pow2(mmsb1 - mmsb2
     )*pow2(-mmgl + mmsb2)) + (8*lb1u*lb2u*pow4(g3)*pow4(mmsb2))/(9.*(-mmgl +
     mmsb1)*pow3(mmsb1 - mmsb2)) - (8*lb1u*lgu*pow4(g3)*pow4(mmsb2))/(9.*(-mmgl
      + mmsb1)*pow3(mmsb1 - mmsb2)) - (8*lb2u*lgu*pow4(g3)*pow4(mmsb2))/(9.*(-
     mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) - (8*lb1u*lb2u*pow4(g3)*pow4(mmsb2))/(
     9.*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) + (8*lb1u*lgu*pow4(g3)*pow4(mmsb2)
     )/(9.*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) + (8*lb2u*lgu*pow4(g3)*pow4(
     mmsb2))/(9.*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) + (8*pow2(lgu)*pow4(g3)*
     pow4(mmsb2))/(9.*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) - (8*pow2(lgu)*pow4(
     g3)*pow4(mmsb2))/(9.*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) + (4*lb2u*ltu*
     pow4(g3)*(-(mmsb2/(-mmgl + mmsb2)) - (3*mmsb2*mmst1*mmt*DeltaInv(mmt,mmsb2
     ,mmst1))/(2.*(-mmgl + mmsb2)) - (3*mmsb2*mmst2*mmt*DeltaInv(mmt,mmsb2,
     mmst2))/(2.*(-mmgl + mmsb2)) - (mmst1*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2
     ))/(-mmgl + mmsb2) - (mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2))/(2.*(-
     mmgl + mmsb2)) - (mmst2*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2))/(-mmgl +
     mmsb2) - (mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2))/(2.*(-mmgl + mmsb2))
     - (3*mmsb2*mmt)/(2.*pow2(-mmgl + mmsb2)) + (2*pow2(mmsb2))/pow2(-mmgl +
     mmsb2) + (3*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2))/pow2(-mmgl +
     mmsb2) + (3*mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2))/pow2(-mmgl +
     mmsb2) + (mmsb2*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/(2.*(-mmgl + mmsb2)
     ) - (DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2)*pow2(mmst1))/pow2(-mmgl + mmsb2
     ) + (mmsb2*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2))/(2.*(-mmgl + mmsb2)) - (
     DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow2(mmst2))/pow2(-mmgl + mmsb2) + (
     DeltaInv(mmt,mmsb2,mmst1)*pow3(mmsb2))/(2.*(-mmgl + mmsb2)) + (DeltaInv(
     mmt,mmsb2,mmst2)*pow3(mmsb2))/(2.*(-mmgl + mmsb2)) + (2*mmst1*DeltaInv(mmt
     ,mmsb2,mmst1)*pow3(mmsb2))/pow2(-mmgl + mmsb2) + (mmt*DeltaInv(mmt,mmsb2,
     mmst1)*pow3(mmsb2))/pow2(-mmgl + mmsb2) + (2*mmst2*DeltaInv(mmt,mmsb2,
     mmst2)*pow3(mmsb2))/pow2(-mmgl + mmsb2) + (mmt*DeltaInv(mmt,mmsb2,mmst2)*
     pow3(mmsb2))/pow2(-mmgl + mmsb2) + (2*mmt*pow2(mmsb2))/pow3(-mmgl + mmsb2)
     - (DeltaInv(mmt,mmsb2,mmst1)*pow4(mmsb2))/pow2(-mmgl + mmsb2) - (DeltaInv(
     mmt,mmsb2,mmst2)*pow4(mmsb2))/pow2(-mmgl + mmsb2)))/3. + (4*lgu*ltu*pow4(
     g3)*(4 - (4*mmsb1)/(-mmgl + mmsb1) - (4*mmsb2)/(-mmgl + mmsb2) - (3*mmt)/(
     2.*(-mmgl + mmsb1)) - (3*mmt)/(2.*(-mmgl + mmsb2)) - 2*mmgl*mmsb1*DeltaInv
     (mmt,mmst1,mmgl) - 2*mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl) + 4*mmgl*mmst1*
     DeltaInv(mmt,mmst1,mmgl) + 4*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl) + 4*
     mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl) + 2*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)
     + 2*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl) + 2*mmsb2*mmt*DeltaInv(mmt,mmst1,
     mmgl) + 6*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl) - (6*mmsb1*mmst1*mmt*DeltaInv
     (mmt,mmst1,mmgl))/(-mmgl + mmsb1) - (6*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,
     mmgl))/(-mmgl + mmsb2) - 2*mmgl*mmsb1*DeltaInv(mmt,mmst2,mmgl) - 2*mmgl*
     mmsb2*DeltaInv(mmt,mmst2,mmgl) + 4*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl) + 4
     *mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl) + 4*mmsb2*mmst2*DeltaInv(mmt,mmst2,
     mmgl) + 2*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl) + 2*mmsb1*mmt*DeltaInv(mmt,
     mmst2,mmgl) + 2*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl) + 6*mmst2*mmt*DeltaInv(
     mmt,mmst2,mmgl) - (6*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl))/(-mmgl +
     mmsb1) - (6*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl))/(-mmgl + mmsb2) - 2*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmgl) - 2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmgl
     ) - 3*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1) - (6*mmst1*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmsb1))/(-mmgl + mmsb1) - (3*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmsb1))/(-mmgl + mmsb1) - 3*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1) - (6*
     mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/(-mmgl + mmsb1) - (3*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/(-mmgl + mmsb1) + (7*mmsb1*mmt)/(2.*
     pow2(-mmgl + mmsb1)) + (2*pow2(mmsb1))/pow2(-mmgl + mmsb1) + (3*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1))/pow2(-mmgl + mmsb1) + (3*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/pow2(-mmgl + mmsb1) - 3*DeltaInv(mmt
     ,mmst1,mmgl)*pow2(mmsb2) - (6*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2))/
     (-mmgl + mmsb2) - (3*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2))/(-mmgl +
     mmsb2) - 3*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2) - (6*mmst2*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmsb2))/(-mmgl + mmsb2) - (3*mmt*DeltaInv(mmt,mmst2,mmgl)
     *pow2(mmsb2))/(-mmgl + mmsb2) + (7*mmsb2*mmt)/(2.*pow2(-mmgl + mmsb2)) + (
     2*pow2(mmsb2))/pow2(-mmgl + mmsb2) + (3*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)
     *pow2(mmsb2))/pow2(-mmgl + mmsb2) + (3*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmsb2))/pow2(-mmgl + mmsb2) - 2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)
     + (2*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(-mmgl + mmsb1) + (2*
     mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(-mmgl + mmsb2) - (DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(mmst1))/pow2(-mmgl + mmsb1) - (DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(mmst1))/pow2(-mmgl + mmsb2) - 2*DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmst2) + (2*mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2
     ))/(-mmgl + mmsb1) + (2*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(-mmgl
      + mmsb2) - (DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(mmst2))/pow2(-mmgl
     + mmsb1) - (DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(mmst2))/pow2(-mmgl +
     mmsb2) + (4*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1))/(-mmgl + mmsb1) + (4*
     DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1))/(-mmgl + mmsb1) + (2*mmst1*DeltaInv(
     mmt,mmst1,mmgl)*pow3(mmsb1))/pow2(-mmgl + mmsb1) + (mmt*DeltaInv(mmt,mmst1
     ,mmgl)*pow3(mmsb1))/pow2(-mmgl + mmsb1) + (2*mmst2*DeltaInv(mmt,mmst2,mmgl
     )*pow3(mmsb1))/pow2(-mmgl + mmsb1) + (mmt*DeltaInv(mmt,mmst2,mmgl)*pow3(
     mmsb1))/pow2(-mmgl + mmsb1) - (2*mmt*pow2(mmsb1))/pow3(-mmgl + mmsb1) + (4
     *DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2))/(-mmgl + mmsb2) + (4*DeltaInv(mmt,
     mmst2,mmgl)*pow3(mmsb2))/(-mmgl + mmsb2) + (2*mmst1*DeltaInv(mmt,mmst1,
     mmgl)*pow3(mmsb2))/pow2(-mmgl + mmsb2) + (mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow3(mmsb2))/pow2(-mmgl + mmsb2) + (2*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow3(
     mmsb2))/pow2(-mmgl + mmsb2) + (mmt*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2))/
     pow2(-mmgl + mmsb2) - (2*mmt*pow2(mmsb2))/pow3(-mmgl + mmsb2) - (DeltaInv(
     mmt,mmst1,mmgl)*pow4(mmsb1))/pow2(-mmgl + mmsb1) - (DeltaInv(mmt,mmst2,
     mmgl)*pow4(mmsb1))/pow2(-mmgl + mmsb1) - (DeltaInv(mmt,mmst1,mmgl)*pow4(
     mmsb2))/pow2(-mmgl + mmsb2) - (DeltaInv(mmt,mmst2,mmgl)*pow4(mmsb2))/pow2(
     -mmgl + mmsb2)))/3. + (4*pow2(ltu)*pow4(g3)*(3 - (5*mmsb1)/(2.*(-mmgl +
     mmsb1)) - (5*mmsb2)/(2.*(-mmgl + mmsb2)) - (3*mmst1)/(4.*(-mmgl + mmsb1))
     - (3*mmst1)/(4.*(-mmgl + mmsb2)) - (3*mmst2)/(4.*(-mmgl + mmsb1)) - (3*
     mmst2)/(4.*(-mmgl + mmsb2)) - (3*mmt)/(4.*(-mmgl + mmsb1)) - (3*mmt)/(4.*(
     -mmgl + mmsb2)) - (3*mmsb1*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1))/(2.*(-mmgl
      + mmsb1)) - (3*mmsb1*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2))/(2.*(-mmgl +
     mmsb1)) - (3*mmsb2*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1))/(2.*(-mmgl + mmsb2
     )) - (3*mmsb2*mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2))/(2.*(-mmgl + mmsb2)) -
     mmgl*mmsb1*DeltaInv(mmt,mmst1,mmgl) - mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl)
     + mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl) + mmsb1*mmst1*DeltaInv(mmt,mmst1,
     mmgl) + mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl) + mmgl*mmt*DeltaInv(mmt,mmst1
     ,mmgl) + mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl) + mmsb2*mmt*DeltaInv(mmt,mmst1
     ,mmgl) + 6*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl) - (6*mmsb1*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl))/(-mmgl + mmsb1) - (6*mmsb2*mmst1*mmt*DeltaInv(
     mmt,mmst1,mmgl))/(-mmgl + mmsb2) - mmgl*mmsb1*DeltaInv(mmt,mmst2,mmgl) -
     mmgl*mmsb2*DeltaInv(mmt,mmst2,mmgl) + mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)
     + mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl) + mmsb2*mmst2*DeltaInv(mmt,mmst2,
     mmgl) + mmgl*mmt*DeltaInv(mmt,mmst2,mmgl) + mmsb1*mmt*DeltaInv(mmt,mmst2,
     mmgl) + mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl) + 6*mmst2*mmt*DeltaInv(mmt,
     mmst2,mmgl) - (6*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl))/(-mmgl + mmsb1)
     - (6*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl))/(-mmgl + mmsb2) - DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmgl) - DeltaInv(mmt,mmst2,mmgl)*pow2(mmgl) - (mmst1*
     DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1))/(4.*(-mmgl + mmsb1)) - (mmt*
     DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1))/(4.*(-mmgl + mmsb1)) - (mmst2*
     DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))/(4.*(-mmgl + mmsb1)) - (mmt*
     DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))/(4.*(-mmgl + mmsb1)) - (3*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmsb1))/2. - (3*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmsb1))/(2.*(-mmgl + mmsb1)) - (3*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)
     )/(2.*(-mmgl + mmsb1)) - (3*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/2. - (3*
     mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/(2.*(-mmgl + mmsb1)) - (3*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/(2.*(-mmgl + mmsb1)) + (mmsb1*mmst1)
     /pow2(-mmgl + mmsb1) + (mmsb1*mmst2)/pow2(-mmgl + mmsb1) + (mmsb1*mmt)/
     pow2(-mmgl + mmsb1) + (2*pow2(mmsb1))/pow2(-mmgl + mmsb1) + (3*mmst1*mmt*
     DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1))/pow2(-mmgl + mmsb1) + (3*mmst2*mmt*
     DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))/pow2(-mmgl + mmsb1) + (3*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1))/pow2(-mmgl + mmsb1) + (3*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/pow2(-mmgl + mmsb1) - (mmst1*
     DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2))/(4.*(-mmgl + mmsb2)) - (mmt*
     DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2))/(4.*(-mmgl + mmsb2)) - (mmst2*
     DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2))/(4.*(-mmgl + mmsb2)) - (mmt*
     DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2))/(4.*(-mmgl + mmsb2)) - (3*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmsb2))/2. - (3*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmsb2))/(2.*(-mmgl + mmsb2)) - (3*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)
     )/(2.*(-mmgl + mmsb2)) - (3*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/2. - (3*
     mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/(2.*(-mmgl + mmsb2)) - (3*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/(2.*(-mmgl + mmsb2)) + (mmsb2*mmst1)
     /pow2(-mmgl + mmsb2) + (mmsb2*mmst2)/pow2(-mmgl + mmsb2) + (mmsb2*mmt)/
     pow2(-mmgl + mmsb2) + (2*pow2(mmsb2))/pow2(-mmgl + mmsb2) + (3*mmst1*mmt*
     DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2))/pow2(-mmgl + mmsb2) + (3*mmst2*mmt*
     DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2))/pow2(-mmgl + mmsb2) + (3*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2))/pow2(-mmgl + mmsb2) + (3*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/pow2(-mmgl + mmsb2) - (mmsb1*
     DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1))/(4.*(-mmgl + mmsb1)) - (mmt*
     DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1))/(4.*(-mmgl + mmsb1)) - (mmsb2*
     DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/(4.*(-mmgl + mmsb2)) - (mmt*
     DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/(4.*(-mmgl + mmsb2)) + DeltaInv(mmt
     ,mmst1,mmgl)*pow2(mmst1) - (mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(-
     mmgl + mmsb1) - (mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(-mmgl +
     mmsb2) - (mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(2.*(-mmgl + mmsb1)) -
     (mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(2.*(-mmgl + mmsb2)) + (mmsb1*
     mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1))/(2.*pow2(-mmgl + mmsb1)) + (
     mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(2.*pow2(-mmgl + mmsb1)) +
     (DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow2(mmst1))/(2.*pow2(-mmgl + mmsb1
     )) + (DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(mmst1))/(2.*pow2(-mmgl +
     mmsb1)) + (mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/(2.*pow2(-mmgl
      + mmsb2)) + (mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(2.*pow2(-
     mmgl + mmsb2)) + (DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2)*pow2(mmst1))/(2.*
     pow2(-mmgl + mmsb2)) + (DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(mmst1))/
     (2.*pow2(-mmgl + mmsb2)) - (mmsb1*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2))/(
     4.*(-mmgl + mmsb1)) - (mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2))/(4.*(-
     mmgl + mmsb1)) - (mmsb2*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2))/(4.*(-mmgl
     + mmsb2)) - (mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2))/(4.*(-mmgl + mmsb2
     )) + DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2) - (mmsb1*DeltaInv(mmt,mmst2,mmgl
     )*pow2(mmst2))/(-mmgl + mmsb1) - (mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmst2))/(-mmgl + mmsb2) - (mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(2.*(
     -mmgl + mmsb1)) - (mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(2.*(-mmgl +
     mmsb2)) + (mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2))/(2.*pow2(-mmgl
      + mmsb1)) + (mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(2.*pow2(-
     mmgl + mmsb1)) + (DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*pow2(mmst2))/(2.*
     pow2(-mmgl + mmsb1)) + (DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(mmst2))/
     (2.*pow2(-mmgl + mmsb1)) + (mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2
     ))/(2.*pow2(-mmgl + mmsb2)) + (mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmst2))/(2.*pow2(-mmgl + mmsb2)) + (DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*
     pow2(mmst2))/(2.*pow2(-mmgl + mmsb2)) + (DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmsb2)*pow2(mmst2))/(2.*pow2(-mmgl + mmsb2)) + (DeltaInv(mmt,mmsb1,mmst1)*
     pow3(mmsb1))/(4.*(-mmgl + mmsb1)) + (DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1)
     )/(4.*(-mmgl + mmsb1)) + (2*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1))/(-mmgl +
     mmsb1) + (2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1))/(-mmgl + mmsb1) + (mmst1
     *DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1))/(2.*pow2(-mmgl + mmsb1)) + (mmt*
     DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1))/(2.*pow2(-mmgl + mmsb1)) + (mmst2*
     DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1))/(2.*pow2(-mmgl + mmsb1)) + (mmt*
     DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1))/(2.*pow2(-mmgl + mmsb1)) + (mmst1*
     DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1))/(2.*pow2(-mmgl + mmsb1)) + (mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1))/(2.*pow2(-mmgl + mmsb1)) + (mmst2*
     DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1))/(2.*pow2(-mmgl + mmsb1)) + (mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1))/(2.*pow2(-mmgl + mmsb1)) + (DeltaInv
     (mmt,mmsb2,mmst1)*pow3(mmsb2))/(4.*(-mmgl + mmsb2)) + (DeltaInv(mmt,mmsb2,
     mmst2)*pow3(mmsb2))/(4.*(-mmgl + mmsb2)) + (2*DeltaInv(mmt,mmst1,mmgl)*
     pow3(mmsb2))/(-mmgl + mmsb2) + (2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2))/(-
     mmgl + mmsb2) + (mmst1*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmsb2))/(2.*pow2(-
     mmgl + mmsb2)) + (mmt*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmsb2))/(2.*pow2(-
     mmgl + mmsb2)) + (mmst2*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmsb2))/(2.*pow2(-
     mmgl + mmsb2)) + (mmt*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmsb2))/(2.*pow2(-
     mmgl + mmsb2)) + (mmst1*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2))/(2.*pow2(-
     mmgl + mmsb2)) + (mmt*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2))/(2.*pow2(-mmgl
      + mmsb2)) + (mmst2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2))/(2.*pow2(-mmgl +
     mmsb2)) + (mmt*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2))/(2.*pow2(-mmgl +
     mmsb2)) + (DeltaInv(mmt,mmsb1,mmst1)*pow3(mmst1))/(4.*(-mmgl + mmsb1)) + (
     DeltaInv(mmt,mmsb2,mmst1)*pow3(mmst1))/(4.*(-mmgl + mmsb2)) + (DeltaInv(
     mmt,mmst1,mmgl)*pow3(mmst1))/(2.*(-mmgl + mmsb1)) + (DeltaInv(mmt,mmst1,
     mmgl)*pow3(mmst1))/(2.*(-mmgl + mmsb2)) - (mmsb1*DeltaInv(mmt,mmsb1,mmst1)
     *pow3(mmst1))/(2.*pow2(-mmgl + mmsb1)) - (mmsb1*DeltaInv(mmt,mmst1,mmgl)*
     pow3(mmst1))/(2.*pow2(-mmgl + mmsb1)) - (mmsb2*DeltaInv(mmt,mmsb2,mmst1)*
     pow3(mmst1))/(2.*pow2(-mmgl + mmsb2)) - (mmsb2*DeltaInv(mmt,mmst1,mmgl)*
     pow3(mmst1))/(2.*pow2(-mmgl + mmsb2)) + (DeltaInv(mmt,mmsb1,mmst2)*pow3(
     mmst2))/(4.*(-mmgl + mmsb1)) + (DeltaInv(mmt,mmsb2,mmst2)*pow3(mmst2))/(4.
     *(-mmgl + mmsb2)) + (DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2))/(2.*(-mmgl +
     mmsb1)) + (DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2))/(2.*(-mmgl + mmsb2)) - (
     mmsb1*DeltaInv(mmt,mmsb1,mmst2)*pow3(mmst2))/(2.*pow2(-mmgl + mmsb1)) - (
     mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2))/(2.*pow2(-mmgl + mmsb1)) - (
     mmsb2*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmst2))/(2.*pow2(-mmgl + mmsb2)) - (
     mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2))/(2.*pow2(-mmgl + mmsb2)) - (
     DeltaInv(mmt,mmsb1,mmst1)*pow4(mmsb1))/(2.*pow2(-mmgl + mmsb1)) - (
     DeltaInv(mmt,mmsb1,mmst2)*pow4(mmsb1))/(2.*pow2(-mmgl + mmsb1)) - (
     DeltaInv(mmt,mmst1,mmgl)*pow4(mmsb1))/(2.*pow2(-mmgl + mmsb1)) - (DeltaInv
     (mmt,mmst2,mmgl)*pow4(mmsb1))/(2.*pow2(-mmgl + mmsb1)) - (DeltaInv(mmt,
     mmsb2,mmst1)*pow4(mmsb2))/(2.*pow2(-mmgl + mmsb2)) - (DeltaInv(mmt,mmsb2,
     mmst2)*pow4(mmsb2))/(2.*pow2(-mmgl + mmsb2)) - (DeltaInv(mmt,mmst1,mmgl)*
     pow4(mmsb2))/(2.*pow2(-mmgl + mmsb2)) - (DeltaInv(mmt,mmst2,mmgl)*pow4(
     mmsb2))/(2.*pow2(-mmgl + mmsb2))))/3. + (4*ltu*pow4(g3)*(-
     10.333333333333334 + (15*mmsb1)/(-mmgl + mmsb1) + (15*mmsb2)/(-mmgl +
     mmsb2) + (9*mmst1)/(2.*(-mmgl + mmsb1)) + (9*mmst1)/(2.*(-mmgl + mmsb2)) +
     (9*mmst2)/(2.*(-mmgl + mmsb1)) + (9*mmst2)/(2.*(-mmgl + mmsb2)) + (6*mmt)/
     (-mmgl + mmsb1) + (6*mmt)/(-mmgl + mmsb2) + (9*mmsb1*mmst1*mmt*DeltaInv(
     mmt,mmsb1,mmst1))/(-mmgl + mmsb1) + (9*mmsb1*mmst2*mmt*DeltaInv(mmt,mmsb1,
     mmst2))/(-mmgl + mmsb1) + (9*mmsb2*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1))/(-
     mmgl + mmsb2) + (9*mmsb2*mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2))/(-mmgl +
     mmsb2) + 6*mmgl*mmsb1*DeltaInv(mmt,mmst1,mmgl) + 6*mmgl*mmsb2*DeltaInv(mmt
     ,mmst1,mmgl) - 6*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl) - 6*mmsb1*mmst1*
     DeltaInv(mmt,mmst1,mmgl) - 6*mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl) - 6*mmgl
     *mmt*DeltaInv(mmt,mmst1,mmgl) - 6*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl) - 6*
     mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl) - 36*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)
     + (36*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl))/(-mmgl + mmsb1) + (36*
     mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl))/(-mmgl + mmsb2) + 6*mmgl*mmsb1*
     DeltaInv(mmt,mmst2,mmgl) + 6*mmgl*mmsb2*DeltaInv(mmt,mmst2,mmgl) - 6*mmgl*
     mmst2*DeltaInv(mmt,mmst2,mmgl) - 6*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl) -
     6*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl) - 6*mmgl*mmt*DeltaInv(mmt,mmst2,
     mmgl) - 6*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl) - 6*mmsb2*mmt*DeltaInv(mmt,
     mmst2,mmgl) - 36*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl) + (36*mmsb1*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl))/(-mmgl + mmsb1) + (36*mmsb2*mmst2*mmt*DeltaInv(
     mmt,mmst2,mmgl))/(-mmgl + mmsb2) + 6*DeltaInv(mmt,mmst1,mmgl)*pow2(mmgl) +
     6*DeltaInv(mmt,mmst2,mmgl)*pow2(mmgl) + (3*mmst1*DeltaInv(mmt,mmsb1,mmst1)
     *pow2(mmsb1))/(2.*(-mmgl + mmsb1)) + (3*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2
     (mmsb1))/(2.*(-mmgl + mmsb1)) + (3*mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow2(
     mmsb1))/(2.*(-mmgl + mmsb1)) + (3*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1
     ))/(2.*(-mmgl + mmsb1)) + 9*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1) + (9*
     mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1))/(-mmgl + mmsb1) + (9*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1))/(-mmgl + mmsb1) + 9*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmsb1) + (9*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/(
     -mmgl + mmsb1) + (9*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/(-mmgl +
     mmsb1) - (6*mmsb1*mmst1)/pow2(-mmgl + mmsb1) - (6*mmsb1*mmst2)/pow2(-mmgl
     + mmsb1) - (8*mmsb1*mmt)/pow2(-mmgl + mmsb1) - (12*pow2(mmsb1))/pow2(-mmgl
      + mmsb1) - (18*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1))/pow2(-
     mmgl + mmsb1) - (18*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))/pow2(
     -mmgl + mmsb1) - (18*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1))/pow2(
     -mmgl + mmsb1) - (18*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/pow2(
     -mmgl + mmsb1) + (3*mmst1*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2))/(2.*(-
     mmgl + mmsb2)) + (3*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2))/(2.*(-mmgl
     + mmsb2)) + (3*mmst2*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2))/(2.*(-mmgl +
     mmsb2)) + (3*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2))/(2.*(-mmgl + mmsb2
     )) + 9*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2) + (9*mmst1*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmsb2))/(-mmgl + mmsb2) + (9*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmsb2))/(-mmgl + mmsb2) + 9*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2) + (9*
     mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/(-mmgl + mmsb2) + (9*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/(-mmgl + mmsb2) - (6*mmsb2*mmst1)/
     pow2(-mmgl + mmsb2) - (6*mmsb2*mmst2)/pow2(-mmgl + mmsb2) - (8*mmsb2*mmt)/
     pow2(-mmgl + mmsb2) - (12*pow2(mmsb2))/pow2(-mmgl + mmsb2) - (18*mmst1*mmt
     *DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2))/pow2(-mmgl + mmsb2) - (18*mmst2*
     mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2))/pow2(-mmgl + mmsb2) - (18*mmst1
     *mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2))/pow2(-mmgl + mmsb2) - (18*mmst2
     *mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/pow2(-mmgl + mmsb2) + (3*mmsb1*
     DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1))/(2.*(-mmgl + mmsb1)) + (3*mmt*
     DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1))/(2.*(-mmgl + mmsb1)) + (3*mmsb2*
     DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/(2.*(-mmgl + mmsb2)) + (3*mmt*
     DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/(2.*(-mmgl + mmsb2)) - 6*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmst1) + (6*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1
     ))/(-mmgl + mmsb1) + (6*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(-mmgl
      + mmsb2) + (3*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(-mmgl + mmsb1) +
     (3*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(-mmgl + mmsb2) - (3*mmsb1*
     mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1))/pow2(-mmgl + mmsb1) - (3*mmsb1*
     mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/pow2(-mmgl + mmsb1) - (3*
     DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow2(mmst1))/pow2(-mmgl + mmsb1) - (
     3*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(mmst1))/pow2(-mmgl + mmsb1) -
     (3*mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/pow2(-mmgl + mmsb2) -
     (3*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/pow2(-mmgl + mmsb2) - (
     3*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2)*pow2(mmst1))/pow2(-mmgl + mmsb2) -
     (3*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(mmst1))/pow2(-mmgl + mmsb2) +
     (3*mmsb1*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2))/(2.*(-mmgl + mmsb1)) + (3*
     mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2))/(2.*(-mmgl + mmsb1)) + (3*mmsb2
     *DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2))/(2.*(-mmgl + mmsb2)) + (3*mmt*
     DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2))/(2.*(-mmgl + mmsb2)) - 6*DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmst2) + (6*mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2
     ))/(-mmgl + mmsb1) + (6*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(-mmgl
      + mmsb2) + (3*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(-mmgl + mmsb1) +
     (3*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(-mmgl + mmsb2) - (3*mmsb1*
     mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2))/pow2(-mmgl + mmsb1) - (3*mmsb1*
     mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/pow2(-mmgl + mmsb1) - (3*
     DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*pow2(mmst2))/pow2(-mmgl + mmsb1) - (
     3*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(mmst2))/pow2(-mmgl + mmsb1) -
     (3*mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2))/pow2(-mmgl + mmsb2) -
     (3*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/pow2(-mmgl + mmsb2) - (
     3*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow2(mmst2))/pow2(-mmgl + mmsb2) -
     (3*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(mmst2))/pow2(-mmgl + mmsb2) -
     (3*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1))/(2.*(-mmgl + mmsb1)) - (3*
     DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1))/(2.*(-mmgl + mmsb1)) - (12*DeltaInv
     (mmt,mmst1,mmgl)*pow3(mmsb1))/(-mmgl + mmsb1) - (12*DeltaInv(mmt,mmst2,
     mmgl)*pow3(mmsb1))/(-mmgl + mmsb1) - (3*mmst1*DeltaInv(mmt,mmsb1,mmst1)*
     pow3(mmsb1))/pow2(-mmgl + mmsb1) - (3*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow3(
     mmsb1))/pow2(-mmgl + mmsb1) - (3*mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow3(
     mmsb1))/pow2(-mmgl + mmsb1) - (3*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1)
     )/pow2(-mmgl + mmsb1) - (3*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1))/
     pow2(-mmgl + mmsb1) - (3*mmt*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1))/pow2(-
     mmgl + mmsb1) - (3*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1))/pow2(-mmgl
     + mmsb1) - (3*mmt*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1))/pow2(-mmgl + mmsb1
     ) - (3*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmsb2))/(2.*(-mmgl + mmsb2)) - (3*
     DeltaInv(mmt,mmsb2,mmst2)*pow3(mmsb2))/(2.*(-mmgl + mmsb2)) - (12*DeltaInv
     (mmt,mmst1,mmgl)*pow3(mmsb2))/(-mmgl + mmsb2) - (12*DeltaInv(mmt,mmst2,
     mmgl)*pow3(mmsb2))/(-mmgl + mmsb2) - (3*mmst1*DeltaInv(mmt,mmsb2,mmst1)*
     pow3(mmsb2))/pow2(-mmgl + mmsb2) - (3*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow3(
     mmsb2))/pow2(-mmgl + mmsb2) - (3*mmst2*DeltaInv(mmt,mmsb2,mmst2)*pow3(
     mmsb2))/pow2(-mmgl + mmsb2) - (3*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmsb2)
     )/pow2(-mmgl + mmsb2) - (3*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2))/
     pow2(-mmgl + mmsb2) - (3*mmt*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2))/pow2(-
     mmgl + mmsb2) - (3*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2))/pow2(-mmgl
     + mmsb2) - (3*mmt*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2))/pow2(-mmgl + mmsb2
     ) - (3*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmst1))/(2.*(-mmgl + mmsb1)) - (3*
     DeltaInv(mmt,mmsb2,mmst1)*pow3(mmst1))/(2.*(-mmgl + mmsb2)) - (3*DeltaInv(
     mmt,mmst1,mmgl)*pow3(mmst1))/(-mmgl + mmsb1) - (3*DeltaInv(mmt,mmst1,mmgl)
     *pow3(mmst1))/(-mmgl + mmsb2) + (3*mmsb1*DeltaInv(mmt,mmsb1,mmst1)*pow3(
     mmst1))/pow2(-mmgl + mmsb1) + (3*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1
     ))/pow2(-mmgl + mmsb1) + (3*mmsb2*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmst1))/
     pow2(-mmgl + mmsb2) + (3*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1))/pow2(
     -mmgl + mmsb2) - (3*DeltaInv(mmt,mmsb1,mmst2)*pow3(mmst2))/(2.*(-mmgl +
     mmsb1)) - (3*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmst2))/(2.*(-mmgl + mmsb2)) -
     (3*DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2))/(-mmgl + mmsb1) - (3*DeltaInv(mmt
     ,mmst2,mmgl)*pow3(mmst2))/(-mmgl + mmsb2) + (3*mmsb1*DeltaInv(mmt,mmsb1,
     mmst2)*pow3(mmst2))/pow2(-mmgl + mmsb1) + (3*mmsb1*DeltaInv(mmt,mmst2,mmgl
     )*pow3(mmst2))/pow2(-mmgl + mmsb1) + (3*mmsb2*DeltaInv(mmt,mmsb2,mmst2)*
     pow3(mmst2))/pow2(-mmgl + mmsb2) + (3*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow3(
     mmst2))/pow2(-mmgl + mmsb2) + (3*DeltaInv(mmt,mmsb1,mmst1)*pow4(mmsb1))/
     pow2(-mmgl + mmsb1) + (3*DeltaInv(mmt,mmsb1,mmst2)*pow4(mmsb1))/pow2(-mmgl
      + mmsb1) + (3*DeltaInv(mmt,mmst1,mmgl)*pow4(mmsb1))/pow2(-mmgl + mmsb1) +
     (3*DeltaInv(mmt,mmst2,mmgl)*pow4(mmsb1))/pow2(-mmgl + mmsb1) + (3*DeltaInv
     (mmt,mmsb2,mmst1)*pow4(mmsb2))/pow2(-mmgl + mmsb2) + (3*DeltaInv(mmt,mmsb2
     ,mmst2)*pow4(mmsb2))/pow2(-mmgl + mmsb2) + (3*DeltaInv(mmt,mmst1,mmgl)*
     pow4(mmsb2))/pow2(-mmgl + mmsb2) + (3*DeltaInv(mmt,mmst2,mmgl)*pow4(mmsb2)
     )/pow2(-mmgl + mmsb2)))/3. - (4*lb2u*lgu*pow4(g3)*pow4(mmsb2))/(3.*pow4(-
     mmgl + mmsb2)) + (2*pow2(lb2u)*pow4(g3)*pow4(mmsb2))/(3.*pow4(-mmgl +
     mmsb2)) + (2*pow2(lgu)*pow4(g3)*pow4(mmsb2))/(3.*pow4(-mmgl + mmsb2)))/
     pow4(g3) + (pow2(Xb)*((2560*mmb*pow4(g3))/(9.*pow2(mmsb1 - mmsb2)) - (256*
     lb1u*mmb*pow4(g3))/(3.*pow2(mmsb1 - mmsb2)) - (256*lb2u*mmb*pow4(g3))/(9.*
     pow2(mmsb1 - mmsb2)) - (512*lgu*mmb*pow4(g3))/(3.*pow2(mmsb1 - mmsb2)) + (
     256*lb1u*lgu*mmb*pow4(g3))/(9.*pow2(mmsb1 - mmsb2)) + (256*lb2u*lgu*mmb*
     pow4(g3))/(9.*pow2(mmsb1 - mmsb2)) - (3284*mmb*mmsb1*pow4(g3))/(9.*(-mmgl
     + mmsb1)*pow2(mmsb1 - mmsb2)) + (168*lb1u*mmb*mmsb1*pow4(g3))/((-mmgl +
     mmsb1)*pow2(mmsb1 - mmsb2)) + (1472*lgu*mmb*mmsb1*pow4(g3))/(9.*(-mmgl +
     mmsb1)*pow2(mmsb1 - mmsb2)) - (448*lb1u*lgu*mmb*mmsb1*pow4(g3))/(9.*(-mmgl
      + mmsb1)*pow2(mmsb1 - mmsb2)) - (5152*mmb*mmsb2*pow4(g3))/(9.*(-mmgl +
     mmsb1)*pow2(mmsb1 - mmsb2)) + (248*lb1u*mmb*mmsb2*pow4(g3))/(3.*(-mmgl +
     mmsb1)*pow2(mmsb1 - mmsb2)) + (376*lb2u*mmb*mmsb2*pow4(g3))/(9.*(-mmgl +
     mmsb1)*pow2(mmsb1 - mmsb2)) + (4*lb1u*lb2u*mmb*mmsb2*pow4(g3))/(3.*(-mmgl
     + mmsb1)*pow2(mmsb1 - mmsb2)) + (3296*lgu*mmb*mmsb2*pow4(g3))/(9.*(-mmgl +
     mmsb1)*pow2(mmsb1 - mmsb2)) - (52*lb1u*lgu*mmb*mmsb2*pow4(g3))/(9.*(-mmgl
     + mmsb1)*pow2(mmsb1 - mmsb2)) - (116*lb2u*lgu*mmb*mmsb2*pow4(g3))/(9.*(-
     mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (32*mmb*mmsb1*pow4(g3))/(3.*(-mmgl +
     mmsb2)*pow2(mmsb1 - mmsb2)) + (32*lb1u*mmb*mmsb1*pow4(g3))/(3.*(-mmgl +
     mmsb2)*pow2(mmsb1 - mmsb2)) - (64*lb2u*mmb*mmsb1*pow4(g3))/(9.*(-mmgl +
     mmsb2)*pow2(mmsb1 - mmsb2)) + (64*lb1u*lb2u*mmb*mmsb1*pow4(g3))/(9.*(-mmgl
      + mmsb2)*pow2(mmsb1 - mmsb2)) + (64*lgu*mmb*mmsb1*pow4(g3))/(9.*(-mmgl +
     mmsb2)*pow2(mmsb1 - mmsb2)) - (64*lb1u*lgu*mmb*mmsb1*pow4(g3))/(9.*(-mmgl
     + mmsb2)*pow2(mmsb1 - mmsb2)) + (1772*mmb*mmsb2*pow4(g3))/(9.*(-mmgl +
     mmsb2)*pow2(mmsb1 - mmsb2)) - (296*lb1u*mmb*mmsb2*pow4(g3))/(9.*(-mmgl +
     mmsb2)*pow2(mmsb1 - mmsb2)) + (80*lb2u*mmb*mmsb2*pow4(g3))/((-mmgl + mmsb2
     )*pow2(mmsb1 - mmsb2)) + (52*lb1u*lb2u*mmb*mmsb2*pow4(g3))/(9.*(-mmgl +
     mmsb2)*pow2(mmsb1 - mmsb2)) - (1760*lgu*mmb*mmsb2*pow4(g3))/(9.*(-mmgl +
     mmsb2)*pow2(mmsb1 - mmsb2)) + (52*lb1u*lgu*mmb*mmsb2*pow4(g3))/(9.*(-mmgl
     + mmsb2)*pow2(mmsb1 - mmsb2)) - (44*lb2u*lgu*mmb*mmsb2*pow4(g3))/((-mmgl +
     mmsb2)*pow2(mmsb1 - mmsb2)) + (512*mmb*zt2*pow4(g3))/(9.*pow2(mmsb1 -
     mmsb2)) - (556*mmb*mmsb1*zt2*pow4(g3))/(9.*(-mmgl + mmsb1)*pow2(mmsb1 -
     mmsb2)) - (80*mmb*mmsb2*zt2*pow4(g3))/((-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)
     ) + (164*mmb*mmsb2*zt2*pow4(g3))/(9.*(-mmgl + mmsb2)*pow2(mmsb1 - mmsb2))
     - (556*mmb*Fin20(mmsb1,mmgl,mmu)*pow4(g3))/(9.*(-mmgl + mmsb1)*pow2(mmsb1
     - mmsb2)) + (44*mmb*Fin20(mmsb1,mmgl,mmu)*pow4(g3))/(9.*(-mmgl + mmsb2)*
     pow2(mmsb1 - mmsb2)) - (44*mmb*Fin20(mmsb1,mmsb2,mmu)*pow4(g3))/(9.*(-mmgl
      + mmsb1)*pow2(mmsb1 - mmsb2)) - (44*mmb*Fin20(mmsb1,mmsb2,mmu)*pow4(g3))/
     (9.*(-mmgl + mmsb2)*pow2(mmsb1 - mmsb2)) - (20*mmb*Fin20(mmsb2,mmgl,mmu)*
     pow4(g3))/(3.*(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) + (20*mmb*Fin20(mmsb2,
     mmgl,mmu)*pow4(g3))/(3.*(-mmgl + mmsb2)*pow2(mmsb1 - mmsb2)) + (256*mmb*
     pow2(lb1u)*pow4(g3))/(9.*pow2(mmsb1 - mmsb2)) - (136*mmb*mmsb1*pow2(lb1u)*
     pow4(g3))/(3.*(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (308*mmb*mmsb2*pow2(
     lb1u)*pow4(g3))/(9.*(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) + (52*mmb*mmsb2*
     pow2(lb1u)*pow4(g3))/(9.*(-mmgl + mmsb2)*pow2(mmsb1 - mmsb2)) - (52*mmb*
     mmsb2*pow2(lb2u)*pow4(g3))/(9.*(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (100
     *mmb*mmsb2*pow2(lb2u)*pow4(g3))/(9.*(-mmgl + mmsb2)*pow2(mmsb1 - mmsb2)) +
     (256*mmb*pow2(lgu)*pow4(g3))/(9.*pow2(mmsb1 - mmsb2)) - (256*mmb*mmsb1*
     pow2(lgu)*pow4(g3))/(9.*(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (308*mmb*
     mmsb2*pow2(lgu)*pow4(g3))/(3.*(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) + (668*
     mmb*mmsb2*pow2(lgu)*pow4(g3))/(9.*(-mmgl + mmsb2)*pow2(mmsb1 - mmsb2)) + (
     64*mmb*mmsb1*mmsb2*pow4(g3))/(9.*pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2))
     + (16*lb1u*mmb*mmsb1*mmsb2*pow4(g3))/(pow2(-mmgl + mmsb1)*pow2(mmsb1 -
     mmsb2)) - (64*lb2u*mmb*mmsb1*mmsb2*pow4(g3))/(9.*pow2(-mmgl + mmsb1)*pow2(
     mmsb1 - mmsb2)) - (32*lb1u*lb2u*mmb*mmsb1*mmsb2*pow4(g3))/(3.*pow2(-mmgl +
     mmsb1)*pow2(mmsb1 - mmsb2)) - (16*lgu*mmb*mmsb1*mmsb2*pow4(g3))/(pow2(-
     mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (32*lb1u*lgu*mmb*mmsb1*mmsb2*pow4(g3)
     )/(9.*pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) + (32*lb2u*lgu*mmb*mmsb1*
     mmsb2*pow4(g3))/(3.*pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) + (104*mmb*
     mmsb1*Fin20(mmsb1,mmsb2,mmu)*pow4(g3))/(9.*pow2(-mmgl + mmsb1)*pow2(mmsb1
     - mmsb2)) - (104*mmb*mmsb1*Fin20(mmsb2,mmgl,mmu)*pow4(g3))/(9.*pow2(-mmgl
     + mmsb1)*pow2(mmsb1 - mmsb2)) + (32*mmb*mmsb1*mmsb2*pow2(lgu)*pow4(g3))/(
     9.*pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (8*mmb*pow2(mmsb1)*pow4(g3))
     /(pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (920*lb1u*mmb*pow2(mmsb1)*
     pow4(g3))/(9.*pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) + (328*lgu*mmb*pow2
     (mmsb1)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) + (12*lb1u*
     lgu*mmb*pow2(mmsb1)*pow4(g3))/(pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) +
     (242*mmb*pow2(lb1u)*pow2(mmsb1)*pow4(g3))/(9.*pow2(-mmgl + mmsb1)*pow2(
     mmsb1 - mmsb2)) - (350*mmb*pow2(lgu)*pow2(mmsb1)*pow4(g3))/(9.*pow2(-mmgl
     + mmsb1)*pow2(mmsb1 - mmsb2)) + (16*lb1u*mmb*pow2(mmsb2)*pow4(g3))/(9.*
     pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) + (16*lb1u*lb2u*mmb*pow2(mmsb2)*
     pow4(g3))/(3.*pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (16*lgu*mmb*pow2(
     mmsb2)*pow4(g3))/(9.*pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (16*lb1u*
     lgu*mmb*pow2(mmsb2)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2))
     - (16*lb2u*lgu*mmb*pow2(mmsb2)*pow4(g3))/(3.*pow2(-mmgl + mmsb1)*pow2(
     mmsb1 - mmsb2)) + (16*mmb*pow2(lgu)*pow2(mmsb2)*pow4(g3))/(3.*pow2(-mmgl +
     mmsb1)*pow2(mmsb1 - mmsb2)) + (64*mmb*mmsb1*mmsb2*pow4(g3))/(9.*pow2(mmsb1
      - mmsb2)*pow2(-mmgl + mmsb2)) - (64*lb1u*mmb*mmsb1*mmsb2*pow4(g3))/(9.*
     pow2(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (128*lb2u*mmb*mmsb1*mmsb2*pow4(
     g3))/(9.*pow2(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (128*lb1u*lb2u*mmb*
     mmsb1*mmsb2*pow4(g3))/(9.*pow2(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (128*
     lgu*mmb*mmsb1*mmsb2*pow4(g3))/(9.*pow2(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2))
     + (128*lb1u*lgu*mmb*mmsb1*mmsb2*pow4(g3))/(9.*pow2(mmsb1 - mmsb2)*pow2(-
     mmgl + mmsb2)) - (104*mmb*mmsb2*Fin20(mmsb1,mmgl,mmu)*pow4(g3))/(9.*pow2(
     mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (104*mmb*mmsb2*Fin20(mmsb1,mmsb2,mmu
     )*pow4(g3))/(9.*pow2(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (8*mmb*pow2(
     mmsb2)*pow4(g3))/(pow2(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (104*lb2u*mmb
     *pow2(mmsb2)*pow4(g3))/(pow2(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (16*
     lb1u*lb2u*mmb*pow2(mmsb2)*pow4(g3))/(9.*pow2(mmsb1 - mmsb2)*pow2(-mmgl +
     mmsb2)) + (1000*lgu*mmb*pow2(mmsb2)*pow4(g3))/(9.*pow2(mmsb1 - mmsb2)*pow2
     (-mmgl + mmsb2)) + (16*lb1u*lgu*mmb*pow2(mmsb2)*pow4(g3))/(9.*pow2(mmsb1 -
     mmsb2)*pow2(-mmgl + mmsb2)) + (124*lb2u*lgu*mmb*pow2(mmsb2)*pow4(g3))/(9.*
     pow2(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (242*mmb*pow2(lb2u)*pow2(mmsb2)
     *pow4(g3))/(9.*pow2(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (122*mmb*pow2(
     lgu)*pow2(mmsb2)*pow4(g3))/(3.*pow2(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) -
     (64*lb1u*mmb*mmsb2*pow2(mmsb1)*pow4(g3))/(9.*pow2(mmsb1 - mmsb2)*pow3(-
     mmgl + mmsb1)) + (64*lb1u*lb2u*mmb*mmsb2*pow2(mmsb1)*pow4(g3))/(9.*pow2(
     mmsb1 - mmsb2)*pow3(-mmgl + mmsb1)) + (64*lgu*mmb*mmsb2*pow2(mmsb1)*pow4(
     g3))/(9.*pow2(mmsb1 - mmsb2)*pow3(-mmgl + mmsb1)) - (64*lb2u*lgu*mmb*mmsb2
     *pow2(mmsb1)*pow4(g3))/(9.*pow2(mmsb1 - mmsb2)*pow3(-mmgl + mmsb1)) + (80*
     lb1u*mmb*pow3(mmsb1)*pow4(g3))/(9.*pow2(mmsb1 - mmsb2)*pow3(-mmgl + mmsb1)
     ) - (80*lgu*mmb*pow3(mmsb1)*pow4(g3))/(9.*pow2(mmsb1 - mmsb2)*pow3(-mmgl +
     mmsb1)) - (32*mmb*pow2(lb1u)*pow3(mmsb1)*pow4(g3))/(9.*pow2(mmsb1 - mmsb2)
     *pow3(-mmgl + mmsb1)) + (32*mmb*pow2(lgu)*pow3(mmsb1)*pow4(g3))/(9.*pow2(
     mmsb1 - mmsb2)*pow3(-mmgl + mmsb1)) - (512*lb1u*mmb*mmsb2*pow4(g3))/(9.*
     pow3(mmsb1 - mmsb2)) + (512*lb2u*mmb*mmsb2*pow4(g3))/(9.*pow3(mmsb1 -
     mmsb2)) + (512*mmb*Fin20(mmsb1,mmgl,mmu)*pow4(g3))/(9.*pow3(mmsb1 - mmsb2)
     ) - (616*mmb*mmsb2*Fin20(mmsb1,mmgl,mmu)*pow4(g3))/(9.*(-mmgl + mmsb1)*
     pow3(mmsb1 - mmsb2)) + (104*mmb*mmsb2*Fin20(mmsb1,mmgl,mmu)*pow4(g3))/(9.*
     (-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) - (512*mmb*Fin20(mmsb2,mmgl,mmu)*pow4
     (g3))/(9.*pow3(mmsb1 - mmsb2)) - (104*mmb*mmsb2*Fin20(mmsb2,mmgl,mmu)*pow4
     (g3))/(9.*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) + (616*mmb*mmsb2*Fin20(
     mmsb2,mmgl,mmu)*pow4(g3))/(9.*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) + (256*
     mmb*mmsb2*pow2(lb1u)*pow4(g3))/(9.*pow3(mmsb1 - mmsb2)) - (256*mmb*mmsb2*
     pow2(lb2u)*pow4(g3))/(9.*pow3(mmsb1 - mmsb2)) - (5056*mmb*pow2(mmsb2)*pow4
     (g3))/(9.*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) + (824*lb1u*mmb*pow2(mmsb2)
     *pow4(g3))/(9.*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) + (88*lb2u*mmb*pow2(
     mmsb2)*pow4(g3))/(3.*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) - (20*lb1u*lb2u*
     mmb*pow2(mmsb2)*pow4(g3))/(9.*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) + (3232
     *lgu*mmb*pow2(mmsb2)*pow4(g3))/(9.*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) -
     (28*lb1u*lgu*mmb*pow2(mmsb2)*pow4(g3))/(3.*(-mmgl + mmsb1)*pow3(mmsb1 -
     mmsb2)) - (28*lb2u*lgu*mmb*pow2(mmsb2)*pow4(g3))/(3.*(-mmgl + mmsb1)*pow3(
     mmsb1 - mmsb2)) + (5056*mmb*pow2(mmsb2)*pow4(g3))/(9.*(-mmgl + mmsb2)*pow3
     (mmsb1 - mmsb2)) - (104*lb1u*mmb*pow2(mmsb2)*pow4(g3))/(3.*(-mmgl + mmsb2)
     *pow3(mmsb1 - mmsb2)) - (776*lb2u*mmb*pow2(mmsb2)*pow4(g3))/(9.*(-mmgl +
     mmsb2)*pow3(mmsb1 - mmsb2)) + (20*lb1u*lb2u*mmb*pow2(mmsb2)*pow4(g3))/(9.*
     (-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) - (3232*lgu*mmb*pow2(mmsb2)*pow4(g3))
     /(9.*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) + (28*lb1u*lgu*mmb*pow2(mmsb2)*
     pow4(g3))/(3.*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) + (28*lb2u*lgu*mmb*pow2
     (mmsb2)*pow4(g3))/(3.*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) - (80*mmb*zt2*
     pow2(mmsb2)*pow4(g3))/((-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) + (80*mmb*zt2*
     pow2(mmsb2)*pow4(g3))/((-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) - (308*mmb*
     pow2(lb1u)*pow2(mmsb2)*pow4(g3))/(9.*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2))
     + (52*mmb*pow2(lb1u)*pow2(mmsb2)*pow4(g3))/(9.*(-mmgl + mmsb2)*pow3(mmsb1
     - mmsb2)) - (52*mmb*pow2(lb2u)*pow2(mmsb2)*pow4(g3))/(9.*(-mmgl + mmsb1)*
     pow3(mmsb1 - mmsb2)) + (308*mmb*pow2(lb2u)*pow2(mmsb2)*pow4(g3))/(9.*(-
     mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) - (892*mmb*pow2(lgu)*pow2(mmsb2)*pow4(
     g3))/(9.*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) + (892*mmb*pow2(lgu)*pow2(
     mmsb2)*pow4(g3))/(9.*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) + (16*lb1u*mmb*
     pow3(mmsb2)*pow4(g3))/(9.*pow2(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) + (64*
     lb1u*lb2u*mmb*pow3(mmsb2)*pow4(g3))/(9.*pow2(-mmgl + mmsb1)*pow3(mmsb1 -
     mmsb2)) - (16*lgu*mmb*pow3(mmsb2)*pow4(g3))/(9.*pow2(-mmgl + mmsb1)*pow3(
     mmsb1 - mmsb2)) - (64*lb1u*lgu*mmb*pow3(mmsb2)*pow4(g3))/(9.*pow2(-mmgl +
     mmsb1)*pow3(mmsb1 - mmsb2)) - (64*lb2u*lgu*mmb*pow3(mmsb2)*pow4(g3))/(9.*
     pow2(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) + (64*mmb*pow2(lgu)*pow3(mmsb2)*
     pow4(g3))/(9.*pow2(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) - (16*lb2u*mmb*pow3
     (mmsb2)*pow4(g3))/(9.*pow2(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) + (16*lgu*
     mmb*pow3(mmsb2)*pow4(g3))/(9.*pow2(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) - (
     64*lb2u*mmb*mmsb1*pow2(mmsb2)*pow4(g3))/(9.*pow2(mmsb1 - mmsb2)*pow3(-mmgl
      + mmsb2)) + (64*lb1u*lb2u*mmb*mmsb1*pow2(mmsb2)*pow4(g3))/(9.*pow2(mmsb1
     - mmsb2)*pow3(-mmgl + mmsb2)) + (64*lgu*mmb*mmsb1*pow2(mmsb2)*pow4(g3))/(
     9.*pow2(mmsb1 - mmsb2)*pow3(-mmgl + mmsb2)) - (64*lb1u*lgu*mmb*mmsb1*pow2(
     mmsb2)*pow4(g3))/(9.*pow2(mmsb1 - mmsb2)*pow3(-mmgl + mmsb2)) + (80*lb2u*
     mmb*pow3(mmsb2)*pow4(g3))/(9.*pow2(mmsb1 - mmsb2)*pow3(-mmgl + mmsb2)) - (
     80*lgu*mmb*pow3(mmsb2)*pow4(g3))/(9.*pow2(mmsb1 - mmsb2)*pow3(-mmgl +
     mmsb2)) - (32*mmb*pow2(lb2u)*pow3(mmsb2)*pow4(g3))/(9.*pow2(mmsb1 - mmsb2)
     *pow3(-mmgl + mmsb2)) + (32*mmb*pow2(lgu)*pow3(mmsb2)*pow4(g3))/(9.*pow2(
     mmsb1 - mmsb2)*pow3(-mmgl + mmsb2)) + (16*lb1u*lgu*mmb*pow4(g3)*pow4(mmsb1
     ))/(9.*pow2(mmsb1 - mmsb2)*pow4(-mmgl + mmsb1)) - (8*mmb*pow2(lb1u)*pow4(
     g3)*pow4(mmsb1))/(9.*pow2(mmsb1 - mmsb2)*pow4(-mmgl + mmsb1)) - (8*mmb*
     pow2(lgu)*pow4(g3)*pow4(mmsb1))/(9.*pow2(mmsb1 - mmsb2)*pow4(-mmgl + mmsb1
     )) + (16*lb1u*mmb*pow3(mmsb2)*pow4(g3))/(9.*(-mmgl + mmsb1)*pow4(mmsb1 -
     mmsb2)) - (16*lb2u*mmb*pow3(mmsb2)*pow4(g3))/(9.*(-mmgl + mmsb1)*pow4(
     mmsb1 - mmsb2)) + (64*lb1u*lb2u*mmb*pow3(mmsb2)*pow4(g3))/(9.*(-mmgl +
     mmsb1)*pow4(mmsb1 - mmsb2)) - (64*lb1u*lgu*mmb*pow3(mmsb2)*pow4(g3))/(9.*(
     -mmgl + mmsb1)*pow4(mmsb1 - mmsb2)) - (64*lb2u*lgu*mmb*pow3(mmsb2)*pow4(g3
     ))/(9.*(-mmgl + mmsb1)*pow4(mmsb1 - mmsb2)) - (16*lb1u*mmb*pow3(mmsb2)*
     pow4(g3))/(9.*(-mmgl + mmsb2)*pow4(mmsb1 - mmsb2)) + (16*lb2u*mmb*pow3(
     mmsb2)*pow4(g3))/(9.*(-mmgl + mmsb2)*pow4(mmsb1 - mmsb2)) - (64*lb1u*lb2u*
     mmb*pow3(mmsb2)*pow4(g3))/(9.*(-mmgl + mmsb2)*pow4(mmsb1 - mmsb2)) + (64*
     lb1u*lgu*mmb*pow3(mmsb2)*pow4(g3))/(9.*(-mmgl + mmsb2)*pow4(mmsb1 - mmsb2)
     ) + (64*lb2u*lgu*mmb*pow3(mmsb2)*pow4(g3))/(9.*(-mmgl + mmsb2)*pow4(mmsb1
     - mmsb2)) + (64*mmb*pow2(lgu)*pow3(mmsb2)*pow4(g3))/(9.*(-mmgl + mmsb1)*
     pow4(mmsb1 - mmsb2)) - (64*mmb*pow2(lgu)*pow3(mmsb2)*pow4(g3))/(9.*(-mmgl
     + mmsb2)*pow4(mmsb1 - mmsb2)) + (16*lb1u*lb2u*mmb*pow4(g3)*pow4(mmsb2))/(
     9.*pow2(-mmgl + mmsb1)*pow4(mmsb1 - mmsb2)) - (16*lb1u*lgu*mmb*pow4(g3)*
     pow4(mmsb2))/(9.*pow2(-mmgl + mmsb1)*pow4(mmsb1 - mmsb2)) - (16*lb2u*lgu*
     mmb*pow4(g3)*pow4(mmsb2))/(9.*pow2(-mmgl + mmsb1)*pow4(mmsb1 - mmsb2)) + (
     16*mmb*pow2(lgu)*pow4(g3)*pow4(mmsb2))/(9.*pow2(-mmgl + mmsb1)*pow4(mmsb1
     - mmsb2)) + (16*lb1u*lb2u*mmb*pow4(g3)*pow4(mmsb2))/(9.*pow2(-mmgl + mmsb2
     )*pow4(mmsb1 - mmsb2)) - (16*lb1u*lgu*mmb*pow4(g3)*pow4(mmsb2))/(9.*pow2(-
     mmgl + mmsb2)*pow4(mmsb1 - mmsb2)) - (16*lb2u*lgu*mmb*pow4(g3)*pow4(mmsb2)
     )/(9.*pow2(-mmgl + mmsb2)*pow4(mmsb1 - mmsb2)) + (16*mmb*pow2(lgu)*pow4(g3
     )*pow4(mmsb2))/(9.*pow2(-mmgl + mmsb2)*pow4(mmsb1 - mmsb2)) + (16*lb2u*lgu
     *mmb*pow4(g3)*pow4(mmsb2))/(9.*pow2(mmsb1 - mmsb2)*pow4(-mmgl + mmsb2)) -
     (8*mmb*pow2(lb2u)*pow4(g3)*pow4(mmsb2))/(9.*pow2(mmsb1 - mmsb2)*pow4(-mmgl
      + mmsb2)) - (8*mmb*pow2(lgu)*pow4(g3)*pow4(mmsb2))/(9.*pow2(mmsb1 - mmsb2
     )*pow4(-mmgl + mmsb2)) + (32*lb1u*lb2u*mmb*pow4(g3)*pow4(mmsb2))/(9.*(-
     mmgl + mmsb1)*pow5(mmsb1 - mmsb2)) - (32*lb1u*lgu*mmb*pow4(g3)*pow4(mmsb2)
     )/(9.*(-mmgl + mmsb1)*pow5(mmsb1 - mmsb2)) - (32*lb2u*lgu*mmb*pow4(g3)*
     pow4(mmsb2))/(9.*(-mmgl + mmsb1)*pow5(mmsb1 - mmsb2)) - (32*lb1u*lb2u*mmb*
     pow4(g3)*pow4(mmsb2))/(9.*(-mmgl + mmsb2)*pow5(mmsb1 - mmsb2)) + (32*lb1u*
     lgu*mmb*pow4(g3)*pow4(mmsb2))/(9.*(-mmgl + mmsb2)*pow5(mmsb1 - mmsb2)) + (
     32*lb2u*lgu*mmb*pow4(g3)*pow4(mmsb2))/(9.*(-mmgl + mmsb2)*pow5(mmsb1 -
     mmsb2)) + (32*mmb*pow2(lgu)*pow4(g3)*pow4(mmsb2))/(9.*(-mmgl + mmsb1)*pow5
     (mmsb1 - mmsb2)) - (32*mmb*pow2(lgu)*pow4(g3)*pow4(mmsb2))/(9.*(-mmgl +
     mmsb2)*pow5(mmsb1 - mmsb2))))/pow4(g3);

   return pow4(g3) * result * twoLoop;
}

} // namespace mssm_twoloop_mb
} // namespace flexiblesusy
