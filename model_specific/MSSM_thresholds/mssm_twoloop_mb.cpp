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

// This file has been generated at Thu 13 Jul 2017 10:34:54
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

   template <typename T> T pow2(T x)  { return x*x; }
   template <typename T> T pow3(T x)  { return x*x*x; }
   template <typename T> T pow4(T x)  { return x*x*x*x; }
   template <typename T> T pow5(T x)  { return x*x*x*x*x; }

   const double oneLoop = 1./pow2(4*Pi);
   const double twoLoop = pow2(oneLoop);

   template <typename T>
   bool is_zero(T a, T prec = std::numeric_limits<T>::epsilon())
   {
      return std::fabs(a) < prec;
   }

   template <typename T>
   bool is_equal(T a, T b, T prec = std::numeric_limits<T>::epsilon())
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
   double Fin20(double mm1, double mm2, double mmu)
   {
      using std::log;
      const double PI = 3.14159265358979323846264338327950288;

      return (6*(mm1*log(mm1/mmu) + mm2*log(mm2/mmu)) +
         (-mm1 - mm2)*(7 + pow2(PI)/6.) +
         (mm1 - mm2)*(2*dilog(1 - mm1/mm2) +
            pow2(log(mm1/mm2))/2.) +
         ((mm1 + mm2)*pow2(log(mm1/mm2)))/2. -
         2*(mm1*pow2(log(mm1/mmu)) + mm2*pow2(log(mm2/mmu))))/2.;
   }

   double LambdaSquared(double x, double y)
   {
      return pow2(1 - x - y) - 4*x*y;
   }

   /// ClausenCl[2,x]
   double ClausenCl2(double x)
   {
      using std::exp;
      const std::complex<double> img(0.,1.);

      return std::imag(dilog(exp(img*x)));
   }

   /// x < 1 && y < 1, LambdaSquared(x,y) > 0
   double PhiPos(double x, double y)
   {
      const double lambda = std::sqrt(LambdaSquared(x,y));

      return (-(log(x)*log(y))
              + 2*log((1 - lambda + x - y)/2.)*log((1 - lambda - x + y)/2.)
              - 2*dilog((1 - lambda + x - y)/2.)
              - 2*dilog((1 - lambda - x + y)/2.)
              + pow2(Pi)/3.)/lambda;
   }

   /// LambdaSquared(x,y) < 0
   double PhiNeg(double x, double y)
   {
      using std::acos;
      using std::sqrt;
      const double lambda = std::sqrt(-LambdaSquared(x,y));

      return 2*(+ ClausenCl2(2*acos((1 + x - y)/(2.*sqrt(x))))
                + ClausenCl2(2*acos((1 - x + y)/(2.*sqrt(y))))
                + ClausenCl2(2*acos((-1 + x + y)/(2.*sqrt(x*y)))))/lambda;
   }

   double Phi(double x, double y)
   {
      const double lambda = LambdaSquared(x,y);

      if (lambda > 0.)
         return PhiPos(x,y);

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
   double Fin3(double mm1, double mm2, double mm3, double mmu)
   {
      using std::log;

      std::array<double,3> masses = { mm1, mm2, mm3 };
      std::sort(masses.begin(), masses.end());

      const double mm = masses[2];
      const double x = masses[0]/mm;
      const double y = masses[1]/mm;

      const double lambda = LambdaSquared(x,y);

      if (is_zero(lambda, 1e-10)) {
         return -(mm*(2*y*(-3 + 2*log(mm/mmu))*log(y)
                      + log(x)*(2*x*(-3 + 2*log(mm/mmu)) + (-1 + x + y)*log(y))
                      + (1 + x + y)*(7 - 6*log(mm/mmu) + pow2(Pi)/6. + 2*pow2(log(mm/mmu)))
                      + x*pow2(log(x)) + y*pow2(log(y))))/2.;
      }

      return mm*((-7 + 6*log(mm/mmu) + log(x)*log(y)
                  - lambda*Phi(x,y) - pow2(Pi)/6. - 2*pow2(log(mm/mmu)))/2.
                 - (x*(7 - 6*log(mm/mmu) + log(x)*(-6 + 4*log(mm/mmu) + log(y))
                       + pow2(Pi)/6. + 2*pow2(log(mm/mmu)) + pow2(log(x))))/2.
                 - (y*(7 - 6*log(mm/mmu) + (
                     -6 + 4*log(mm/mmu) + log(x))*log(y) + pow2(Pi)/6.
                       + 2*pow2(log(mm/mmu)) + pow2(log(y))))/2.);
   }

   /// Delta[m1,m2,m3,-1]
   double DeltaInv(double m1, double m2, double m3)
   {
      return 1./(pow2(m1) + pow2(m2) + pow2(m3) - 2*(m1*m2 + m1*m3 + m2*m3));
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
   using std::log;
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

   const double result =
   pow4(g3)*(-67.75925925925925 - (8*mmsb1)/(-mmgl + mmsb1) + (890*mmsb2)/(9.*(
     -mmgl + mmsb1)) + (6*mmsb1)/(-mmgl + mmsb2) - (908*mmsb2)/(9.*(-mmgl +
     mmsb2)) - (8*mmst1)/(-mmgl + mmsb1) - (8*mmst1)/(-mmgl + mmsb2) - (8*
     mmst2)/(-mmgl + mmsb1) - (8*mmst2)/(-mmgl + mmsb2) + (48*mmsusy)/(-mmgl +
     mmsb1) + (48*mmsusy)/(-mmgl + mmsb2) - (12*mmt)/(-mmgl + mmsb1) - (12*mmt)
     /(-mmgl + mmsb2) - 8*zt2 - (31*mmsb1*zt2)/(9.*(-mmgl + mmsb1)) + (43*
     mmsb2*zt2)/(3.*(-mmgl + mmsb1)) + (mmsb1*zt2)/(-mmgl + mmsb2) - (151*
     mmsb2*zt2)/(9.*(-mmgl + mmsb2)) - (mmst1*zt2)/(-mmgl + mmsb1) - (mmst1*
     zt2)/(-mmgl + mmsb2) - (mmst2*zt2)/(-mmgl + mmsb1) - (mmst2*zt2)/(-mmgl +
     mmsb2) + (8*mmsusy*zt2)/(-mmgl + mmsb1) + (8*mmsusy*zt2)/(-mmgl + mmsb2) -
     (2*mmt*zt2)/(-mmgl + mmsb1) - (2*mmt*zt2)/(-mmgl + mmsb2) - (56*mmsb1*
     mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1))/(3.*(-mmgl + mmsb1)) - (8*mmsb1*
     mmst1*mmt*zt2*DeltaInv(mmt,mmsb1,mmst1))/(3.*(-mmgl + mmsb1)) - (56*mmsb1*
     mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2))/(3.*(-mmgl + mmsb1)) - (8*mmsb1*
     mmst2*mmt*zt2*DeltaInv(mmt,mmsb1,mmst2))/(3.*(-mmgl + mmsb1)) - (56*mmsb2*
     mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1))/(3.*(-mmgl + mmsb2)) - (8*mmsb2*
     mmst1*mmt*zt2*DeltaInv(mmt,mmsb2,mmst1))/(3.*(-mmgl + mmsb2)) - (56*mmsb2*
     mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2))/(3.*(-mmgl + mmsb2)) - (8*mmsb2*
     mmst2*mmt*zt2*DeltaInv(mmt,mmsb2,mmst2))/(3.*(-mmgl + mmsb2)) - (56*mmgl*
     mmsb1*DeltaInv(mmt,mmst1,mmgl))/3. - (56*mmgl*mmsb2*DeltaInv(mmt,mmst1,
     mmgl))/3. + (56*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl))/3. + (56*mmsb1*mmst1*
     DeltaInv(mmt,mmst1,mmgl))/3. + (56*mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl))/
     3. + (56*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl))/3. + (56*mmsb1*mmt*DeltaInv(
     mmt,mmst1,mmgl))/3. + (56*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl))/3. + (224*
     mmst1*mmt*DeltaInv(mmt,mmst1,mmgl))/3. - (224*mmsb1*mmst1*mmt*DeltaInv(
     mmt,mmst1,mmgl))/(3.*(-mmgl + mmsb1)) - (224*mmsb2*mmst1*mmt*DeltaInv(mmt,
     mmst1,mmgl))/(3.*(-mmgl + mmsb2)) - (8*mmgl*mmsb1*zt2*DeltaInv(mmt,mmst1,
     mmgl))/3. - (8*mmgl*mmsb2*zt2*DeltaInv(mmt,mmst1,mmgl))/3. + (8*mmgl*
     mmst1*zt2*DeltaInv(mmt,mmst1,mmgl))/3. + (8*mmsb1*mmst1*zt2*DeltaInv(mmt,
     mmst1,mmgl))/3. + (8*mmsb2*mmst1*zt2*DeltaInv(mmt,mmst1,mmgl))/3. + (8*
     mmgl*mmt*zt2*DeltaInv(mmt,mmst1,mmgl))/3. + (8*mmsb1*mmt*zt2*DeltaInv(mmt,
     mmst1,mmgl))/3. + (8*mmsb2*mmt*zt2*DeltaInv(mmt,mmst1,mmgl))/3. + (32*
     mmst1*mmt*zt2*DeltaInv(mmt,mmst1,mmgl))/3. - (32*mmsb1*mmst1*mmt*zt2*
     DeltaInv(mmt,mmst1,mmgl))/(3.*(-mmgl + mmsb1)) - (32*mmsb2*mmst1*mmt*zt2*
     DeltaInv(mmt,mmst1,mmgl))/(3.*(-mmgl + mmsb2)) - (56*mmgl*mmsb1*DeltaInv(
     mmt,mmst2,mmgl))/3. - (56*mmgl*mmsb2*DeltaInv(mmt,mmst2,mmgl))/3. + (56*
     mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl))/3. + (56*mmsb1*mmst2*DeltaInv(mmt,
     mmst2,mmgl))/3. + (56*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl))/3. + (56*mmgl*
     mmt*DeltaInv(mmt,mmst2,mmgl))/3. + (56*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl))
     /3. + (56*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl))/3. + (224*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl))/3. - (224*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,
     mmgl))/(3.*(-mmgl + mmsb1)) - (224*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,
     mmgl))/(3.*(-mmgl + mmsb2)) - (8*mmgl*mmsb1*zt2*DeltaInv(mmt,mmst2,mmgl))/
     3. - (8*mmgl*mmsb2*zt2*DeltaInv(mmt,mmst2,mmgl))/3. + (8*mmgl*mmst2*zt2*
     DeltaInv(mmt,mmst2,mmgl))/3. + (8*mmsb1*mmst2*zt2*DeltaInv(mmt,mmst2,mmgl)
     )/3. + (8*mmsb2*mmst2*zt2*DeltaInv(mmt,mmst2,mmgl))/3. + (8*mmgl*mmt*zt2*
     DeltaInv(mmt,mmst2,mmgl))/3. + (8*mmsb1*mmt*zt2*DeltaInv(mmt,mmst2,mmgl))/
     3. + (8*mmsb2*mmt*zt2*DeltaInv(mmt,mmst2,mmgl))/3. + (32*mmst2*mmt*zt2*
     DeltaInv(mmt,mmst2,mmgl))/3. - (32*mmsb1*mmst2*mmt*zt2*DeltaInv(mmt,mmst2,
     mmgl))/(3.*(-mmgl + mmsb1)) - (32*mmsb2*mmst2*mmt*zt2*DeltaInv(mmt,mmst2,
     mmgl))/(3.*(-mmgl + mmsb2)) + (83*Fin20(mmsb1,mmgl,mmu))/(9.*(-mmgl +
     mmsb1)) - (128*Fin20(mmsb1,mmgl,mmu))/(9.*(mmsb1 - mmsb2)) + (40*mmsb2*
     Fin20(mmsb1,mmgl,mmu))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (7*Fin20(
     mmsb1,mmgl,mmu))/(9.*(-mmgl + mmsb2)) + (8*mmsb2*Fin20(mmsb1,mmgl,mmu))/(
     9.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (10*Fin20(mmsb1,mmsb2,mmu))/(9.*(-
     mmgl + mmsb1)) + (8*mmsb2*Fin20(mmsb1,mmsb2,mmu))/(9.*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) + (2*Fin20(mmsb1,mmsb2,mmu))/(9.*(-mmgl + mmsb2)) - (8*
     mmsb2*Fin20(mmsb1,mmsb2,mmu))/(9.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) -
     Fin20(mmsb2,mmgl,mmu)/(9.*(-mmgl + mmsb1)) + (128*Fin20(mmsb2,mmgl,mmu))/(
     9.*(mmsb1 - mmsb2)) - (8*mmsb2*Fin20(mmsb2,mmgl,mmu))/(9.*(-mmgl + mmsb1)*
     (mmsb1 - mmsb2)) - (37*Fin20(mmsb2,mmgl,mmu))/(9.*(-mmgl + mmsb2)) - (40*
     mmsb2*Fin20(mmsb2,mmgl,mmu))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (1016*
     log(mmgl/mmu))/9. - (58*mmsb1*log(mmgl/mmu))/(3.*(-mmgl + mmsb1)) - (598*
     mmsb2*log(mmgl/mmu))/(9.*(-mmgl + mmsb1)) - (2*mmsb1*log(mmgl/mmu))/(-mmgl
     + mmsb2) + (406*mmsb2*log(mmgl/mmu))/(9.*(-mmgl + mmsb2)) - (2*mmst1*log(
     mmgl/mmu))/(-mmgl + mmsb1) - (2*mmst1*log(mmgl/mmu))/(-mmgl + mmsb2) - (2*
     mmst2*log(mmgl/mmu))/(-mmgl + mmsb1) - (2*mmst2*log(mmgl/mmu))/(-mmgl +
     mmsb2) - (16*mmsusy*log(mmgl/mmu))/(-mmgl + mmsb1) - (16*mmsusy*log(mmgl/
     mmu))/(-mmgl + mmsb2) + (4*mmt*log(mmgl/mmu))/(-mmgl + mmsb1) + (4*mmt*
     log(mmgl/mmu))/(-mmgl + mmsb2) + 8*mmgl*mmsb1*DeltaInv(mmt,mmst1,mmgl)*
     log(mmgl/mmu) + 8*mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu) - 16*
     mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu) - 16*mmsb1*mmst1*
     DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu) - 16*mmsb2*mmst1*DeltaInv(mmt,
     mmst1,mmgl)*log(mmgl/mmu) - 8*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/
     mmu) - 8*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu) - 8*mmsb2*mmt*
     DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu) - 8*mmst1*mmt*DeltaInv(mmt,mmst1,
     mmgl)*log(mmgl/mmu) + (8*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(
     mmgl/mmu))/(-mmgl + mmsb1) + (8*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*
     log(mmgl/mmu))/(-mmgl + mmsb2) + 8*mmgl*mmsb1*DeltaInv(mmt,mmst2,mmgl)*
     log(mmgl/mmu) + 8*mmgl*mmsb2*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu) - 16*
     mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu) - 16*mmsb1*mmst2*
     DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu) - 16*mmsb2*mmst2*DeltaInv(mmt,
     mmst2,mmgl)*log(mmgl/mmu) - 8*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/
     mmu) - 8*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu) - 8*mmsb2*mmt*
     DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu) - 8*mmst2*mmt*DeltaInv(mmt,mmst2,
     mmgl)*log(mmgl/mmu) + (8*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(
     mmgl/mmu))/(-mmgl + mmsb1) + (8*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*
     log(mmgl/mmu))/(-mmgl + mmsb2) + (25*log(mmsb1/mmu))/3. + (118*mmsb1*log(
     mmsb1/mmu))/(9.*(-mmgl + mmsb1)) - (44*mmsb2*log(mmsb1/mmu))/(3.*(-mmgl +
     mmsb1)) + (128*mmsb2*log(mmsb1/mmu))/(9.*(mmsb1 - mmsb2)) - (4*mmsb1*log(
     mmsb1/mmu))/(-mmgl + mmsb2) + (4*mmsb2*log(mmsb1/mmu))/(9.*(-mmgl + mmsb2)
     ) + (2*mmsb1*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1)*log(mmsb1/mmu))/(-mmgl +
     mmsb1) + (2*mmsb1*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2)*log(mmsb1/mmu))/(-
     mmgl + mmsb1) + (64*log(mmgl/mmu)*log(mmsb1/mmu))/9. - (55*mmsb1*log(mmgl/
     mmu)*log(mmsb1/mmu))/(9.*(-mmgl + mmsb1)) + (4*mmsb2*log(mmgl/mmu)*log(
     mmsb1/mmu))/(9.*(-mmgl + mmsb1)) + (mmsb1*log(mmgl/mmu)*log(mmsb1/mmu))/(-
     mmgl + mmsb2) - (4*mmsb2*log(mmgl/mmu)*log(mmsb1/mmu))/(9.*(-mmgl + mmsb2)
     ) - (53*log(mmsb2/mmu))/9. - (44*mmsb2*log(mmsb2/mmu))/(9.*(-mmgl + mmsb1)
     ) - (128*mmsb2*log(mmsb2/mmu))/(9.*(mmsb1 - mmsb2)) + (254*mmsb2*log(
     mmsb2/mmu))/(9.*(-mmgl + mmsb2)) + (2*mmsb2*mmst1*mmt*DeltaInv(mmt,mmsb2,
     mmst1)*log(mmsb2/mmu))/(-mmgl + mmsb2) + (2*mmsb2*mmst2*mmt*DeltaInv(mmt,
     mmsb2,mmst2)*log(mmsb2/mmu))/(-mmgl + mmsb2) + (64*log(mmgl/mmu)*log(
     mmsb2/mmu))/9. + (13*mmsb2*log(mmgl/mmu)*log(mmsb2/mmu))/(9.*(-mmgl +
     mmsb1)) - (59*mmsb2*log(mmgl/mmu)*log(mmsb2/mmu))/(9.*(-mmgl + mmsb2)) - (
     4*mmsb2*log(mmsb1/mmu)*log(mmsb2/mmu))/(9.*(-mmgl + mmsb1)) + (4*mmsb2*
     log(mmsb1/mmu)*log(mmsb2/mmu))/(9.*(-mmgl + mmsb2)) - (56*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmgl))/3. - (8*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmgl))/
     3. - (56*DeltaInv(mmt,mmst2,mmgl)*pow2(mmgl))/3. - (8*zt2*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmgl))/3. + 8*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*
     pow2(mmgl) + 8*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*pow2(mmgl) - (14*
     mmst1*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1))/(3.*(-mmgl + mmsb1)) - (14*
     mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1))/(3.*(-mmgl + mmsb1)) - (2*
     mmst1*zt2*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1))/(3.*(-mmgl + mmsb1)) - (
     2*mmt*zt2*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1))/(3.*(-mmgl + mmsb1)) - (
     14*mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))/(3.*(-mmgl + mmsb1)) - (
     14*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))/(3.*(-mmgl + mmsb1)) - (2*
     mmst2*zt2*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))/(3.*(-mmgl + mmsb1)) - (
     2*mmt*zt2*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))/(3.*(-mmgl + mmsb1)) -
     28*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1) - (28*mmst1*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmsb1))/(-mmgl + mmsb1) - (28*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmsb1))/(-mmgl + mmsb1) - 4*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)
     - (4*mmst1*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1))/(-mmgl + mmsb1) - (4*
     mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1))/(-mmgl + mmsb1) - 28*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1) - (28*mmst2*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmsb1))/(-mmgl + mmsb1) - (28*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmsb1))/(-mmgl + mmsb1) - 4*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1) - (4*
     mmst2*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/(-mmgl + mmsb1) - (4*mmt*
     zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/(-mmgl + mmsb1) + 12*DeltaInv(
     mmt,mmst1,mmgl)*log(mmgl/mmu)*pow2(mmsb1) + (24*mmst1*DeltaInv(mmt,mmst1,
     mmgl)*log(mmgl/mmu)*pow2(mmsb1))/(-mmgl + mmsb1) + (12*mmt*DeltaInv(mmt,
     mmst1,mmgl)*log(mmgl/mmu)*pow2(mmsb1))/(-mmgl + mmsb1) + 12*DeltaInv(mmt,
     mmst2,mmgl)*log(mmgl/mmu)*pow2(mmsb1) + (24*mmst2*DeltaInv(mmt,mmst2,mmgl)
     *log(mmgl/mmu)*pow2(mmsb1))/(-mmgl + mmsb1) + (12*mmt*DeltaInv(mmt,mmst2,
     mmgl)*log(mmgl/mmu)*pow2(mmsb1))/(-mmgl + mmsb1) + (4*mmst1*DeltaInv(mmt,
     mmsb1,mmst1)*log(mmsb1/mmu)*pow2(mmsb1))/(-mmgl + mmsb1) + (2*mmt*
     DeltaInv(mmt,mmsb1,mmst1)*log(mmsb1/mmu)*pow2(mmsb1))/(-mmgl + mmsb1) + (
     4*mmst2*DeltaInv(mmt,mmsb1,mmst2)*log(mmsb1/mmu)*pow2(mmsb1))/(-mmgl +
     mmsb1) + (2*mmt*DeltaInv(mmt,mmsb1,mmst2)*log(mmsb1/mmu)*pow2(mmsb1))/(-
     mmgl + mmsb1) - (8*mmsb1*mmsb2)/pow2(-mmgl + mmsb1) + (32*mmsb1*mmst1)/(3.
     *pow2(-mmgl + mmsb1)) + (32*mmsb1*mmst2)/(3.*pow2(-mmgl + mmsb1)) - (64*
     mmsb1*mmsusy)/pow2(-mmgl + mmsb1) + (16*mmsb1*mmt)/pow2(-mmgl + mmsb1) - (
     4*mmsb1*mmsb2*zt2)/(3.*pow2(-mmgl + mmsb1)) + (4*mmsb1*mmst1*zt2)/(3.*
     pow2(-mmgl + mmsb1)) + (4*mmsb1*mmst2*zt2)/(3.*pow2(-mmgl + mmsb1)) - (32*
     mmsb1*mmsusy*zt2)/(3.*pow2(-mmgl + mmsb1)) + (8*mmsb1*mmt*zt2)/(3.*pow2(-
     mmgl + mmsb1)) + (12*mmsb1*Fin20(mmsb1,mmgl,mmu))/pow2(-mmgl + mmsb1) - (
     7*mmsb1*Fin20(mmsb1,mmsb2,mmu))/(3.*pow2(-mmgl + mmsb1)) + (mmsb2*Fin20(
     mmsb1,mmsb2,mmu))/pow2(-mmgl + mmsb1) + (mmsb1*Fin20(mmsb2,mmgl,mmu))/
     pow2(-mmgl + mmsb1) - (mmsb2*Fin20(mmsb2,mmgl,mmu))/pow2(-mmgl + mmsb1) +
     (38*mmsb1*mmsb2*log(mmgl/mmu))/(9.*pow2(-mmgl + mmsb1)) + (14*mmsb1*mmst1*
     log(mmgl/mmu))/(3.*pow2(-mmgl + mmsb1)) + (14*mmsb1*mmst2*log(mmgl/mmu))/(
     3.*pow2(-mmgl + mmsb1)) + (112*mmsb1*mmsusy*log(mmgl/mmu))/(3.*pow2(-mmgl
     + mmsb1)) - (28*mmsb1*mmt*log(mmgl/mmu))/(3.*pow2(-mmgl + mmsb1)) - (14*
     mmsb1*mmsb2*log(mmsb1/mmu))/(9.*pow2(-mmgl + mmsb1)) - (2*mmsb1*mmst1*log(
     mmsb1/mmu))/pow2(-mmgl + mmsb1) - (2*mmsb1*mmst2*log(mmsb1/mmu))/pow2(-
     mmgl + mmsb1) - (16*mmsb1*mmsusy*log(mmsb1/mmu))/pow2(-mmgl + mmsb1) + (4*
     mmsb1*mmt*log(mmsb1/mmu))/pow2(-mmgl + mmsb1) - (8*mmsb1*mmsb2*log(mmgl/
     mmu)*log(mmsb1/mmu))/(9.*pow2(-mmgl + mmsb1)) + (16*mmsb1*mmsb2*log(mmsb2/
     mmu))/(3.*pow2(-mmgl + mmsb1)) - (29*mmsb1*mmsb2*log(mmgl/mmu)*log(mmsb2/
     mmu))/(9.*pow2(-mmgl + mmsb1)) + (8*mmsb1*mmsb2*log(mmsb1/mmu)*log(mmsb2/
     mmu))/(9.*pow2(-mmgl + mmsb1)) + (1222*pow2(mmsb1))/(9.*pow2(-mmgl +
     mmsb1)) + (20*zt2*pow2(mmsb1))/pow2(-mmgl + mmsb1) + (112*mmst1*mmt*
     DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1))/(3.*pow2(-mmgl + mmsb1)) + (16*
     mmst1*mmt*zt2*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1))/(3.*pow2(-mmgl +
     mmsb1)) + (112*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))/(3.*pow2(-
     mmgl + mmsb1)) + (16*mmst2*mmt*zt2*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))/
     (3.*pow2(-mmgl + mmsb1)) + (112*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmsb1))/(3.*pow2(-mmgl + mmsb1)) + (16*mmst1*mmt*zt2*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmsb1))/(3.*pow2(-mmgl + mmsb1)) + (112*mmst2*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmsb1))/(3.*pow2(-mmgl + mmsb1)) + (16*mmst2*mmt*zt2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/(3.*pow2(-mmgl + mmsb1)) - (368*log(
     mmgl/mmu)*pow2(mmsb1))/(3.*pow2(-mmgl + mmsb1)) - (4*mmst1*mmt*DeltaInv(
     mmt,mmst1,mmgl)*log(mmgl/mmu)*pow2(mmsb1))/pow2(-mmgl + mmsb1) - (4*mmst2*
     mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*pow2(mmsb1))/pow2(-mmgl +
     mmsb1) + (160*log(mmsb1/mmu)*pow2(mmsb1))/(9.*pow2(-mmgl + mmsb1)) - (4*
     mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1)*log(mmsb1/mmu)*pow2(mmsb1))/pow2(-mmgl
     + mmsb1) - (4*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2)*log(mmsb1/mmu)*pow2(
     mmsb1))/pow2(-mmgl + mmsb1) - (113*log(mmgl/mmu)*log(mmsb1/mmu)*pow2(
     mmsb1))/(9.*pow2(-mmgl + mmsb1)) + (836*pow2(mmsb2))/(9.*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) - (836*pow2(mmsb2))/(9.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) +
     (40*zt2*pow2(mmsb2))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (40*zt2*pow2(
     mmsb2))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (14*mmst1*DeltaInv(mmt,
     mmsb2,mmst1)*pow2(mmsb2))/(3.*(-mmgl + mmsb2)) - (14*mmt*DeltaInv(mmt,
     mmsb2,mmst1)*pow2(mmsb2))/(3.*(-mmgl + mmsb2)) - (2*mmst1*zt2*DeltaInv(
     mmt,mmsb2,mmst1)*pow2(mmsb2))/(3.*(-mmgl + mmsb2)) - (2*mmt*zt2*DeltaInv(
     mmt,mmsb2,mmst1)*pow2(mmsb2))/(3.*(-mmgl + mmsb2)) - (14*mmst2*DeltaInv(
     mmt,mmsb2,mmst2)*pow2(mmsb2))/(3.*(-mmgl + mmsb2)) - (14*mmt*DeltaInv(mmt,
     mmsb2,mmst2)*pow2(mmsb2))/(3.*(-mmgl + mmsb2)) - (2*mmst2*zt2*DeltaInv(
     mmt,mmsb2,mmst2)*pow2(mmsb2))/(3.*(-mmgl + mmsb2)) - (2*mmt*zt2*DeltaInv(
     mmt,mmsb2,mmst2)*pow2(mmsb2))/(3.*(-mmgl + mmsb2)) - 28*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmsb2) - (28*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2))/
     (-mmgl + mmsb2) - (28*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2))/(-mmgl +
     mmsb2) - 4*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2) - (4*mmst1*zt2*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2))/(-mmgl + mmsb2) - (4*mmt*zt2*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2))/(-mmgl + mmsb2) - 28*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmsb2) - (28*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/
     (-mmgl + mmsb2) - (28*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/(-mmgl +
     mmsb2) - 4*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2) - (4*mmst2*zt2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/(-mmgl + mmsb2) - (4*mmt*zt2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/(-mmgl + mmsb2) - (580*log(mmgl/mmu)
     *pow2(mmsb2))/(9.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (580*log(mmgl/mmu)*
     pow2(mmsb2))/(9.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + 12*DeltaInv(mmt,mmst1,
     mmgl)*log(mmgl/mmu)*pow2(mmsb2) + (24*mmst1*DeltaInv(mmt,mmst1,mmgl)*log(
     mmgl/mmu)*pow2(mmsb2))/(-mmgl + mmsb2) + (12*mmt*DeltaInv(mmt,mmst1,mmgl)*
     log(mmgl/mmu)*pow2(mmsb2))/(-mmgl + mmsb2) + 12*DeltaInv(mmt,mmst2,mmgl)*
     log(mmgl/mmu)*pow2(mmsb2) + (24*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/
     mmu)*pow2(mmsb2))/(-mmgl + mmsb2) + (12*mmt*DeltaInv(mmt,mmst2,mmgl)*log(
     mmgl/mmu)*pow2(mmsb2))/(-mmgl + mmsb2) - (128*log(mmsb1/mmu)*pow2(mmsb2))/
     (9.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (4*log(mmgl/mmu)*log(mmsb1/mmu)*
     pow2(mmsb2))/(9.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (4*log(mmgl/mmu)*log(
     mmsb1/mmu)*pow2(mmsb2))/(9.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (4*log(
     mmsb2/mmu)*pow2(mmsb2))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (140*log(
     mmsb2/mmu)*pow2(mmsb2))/(9.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (4*mmst1*
     DeltaInv(mmt,mmsb2,mmst1)*log(mmsb2/mmu)*pow2(mmsb2))/(-mmgl + mmsb2) + (
     2*mmt*DeltaInv(mmt,mmsb2,mmst1)*log(mmsb2/mmu)*pow2(mmsb2))/(-mmgl +
     mmsb2) + (4*mmst2*DeltaInv(mmt,mmsb2,mmst2)*log(mmsb2/mmu)*pow2(mmsb2))/(-
     mmgl + mmsb2) + (2*mmt*DeltaInv(mmt,mmsb2,mmst2)*log(mmsb2/mmu)*pow2(
     mmsb2))/(-mmgl + mmsb2) - (4*log(mmgl/mmu)*log(mmsb2/mmu)*pow2(mmsb2))/(9.
     *(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (4*log(mmgl/mmu)*log(mmsb2/mmu)*pow2(
     mmsb2))/(9.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (4*log(mmsb1/mmu)*log(
     mmsb2/mmu)*pow2(mmsb2))/(9.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (4*log(
     mmsb1/mmu)*log(mmsb2/mmu)*pow2(mmsb2))/(9.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)
     ) - (4*log(mmgl/mmu)*pow2(mmsb2))/(9.*pow2(-mmgl + mmsb1)) + (4*log(mmsb1/
     mmu)*pow2(mmsb2))/(9.*pow2(-mmgl + mmsb1)) - (4*log(mmgl/mmu)*log(mmsb1/
     mmu)*pow2(mmsb2))/(3.*pow2(-mmgl + mmsb1)) - (4*log(mmgl/mmu)*log(mmsb2/
     mmu)*pow2(mmsb2))/(3.*pow2(-mmgl + mmsb1)) + (4*log(mmsb1/mmu)*log(mmsb2/
     mmu)*pow2(mmsb2))/(3.*pow2(-mmgl + mmsb1)) + (4*log(mmsusy/mmu)*(
     7.333333333333333 - (24*mmsusy)/(-mmgl + mmsb1) - (24*mmsusy)/(-mmgl +
     mmsb2) + (32*mmsb1*mmsusy)/pow2(-mmgl + mmsb1) + (32*mmsb2*mmsusy)/pow2(-
     mmgl + mmsb2)))/3. - (8*mmsb1*mmsb2)/pow2(-mmgl + mmsb2) + (32*mmsb2*
     mmst1)/(3.*pow2(-mmgl + mmsb2)) + (32*mmsb2*mmst2)/(3.*pow2(-mmgl + mmsb2)
     ) - (64*mmsb2*mmsusy)/pow2(-mmgl + mmsb2) + (16*mmsb2*mmt)/pow2(-mmgl +
     mmsb2) - (4*mmsb1*mmsb2*zt2)/(3.*pow2(-mmgl + mmsb2)) + (4*mmsb2*mmst1*
     zt2)/(3.*pow2(-mmgl + mmsb2)) + (4*mmsb2*mmst2*zt2)/(3.*pow2(-mmgl +
     mmsb2)) - (32*mmsb2*mmsusy*zt2)/(3.*pow2(-mmgl + mmsb2)) + (8*mmsb2*mmt*
     zt2)/(3.*pow2(-mmgl + mmsb2)) - (mmsb1*Fin20(mmsb1,mmgl,mmu))/pow2(-mmgl +
     mmsb2) + (mmsb2*Fin20(mmsb1,mmgl,mmu))/pow2(-mmgl + mmsb2) + (mmsb1*Fin20(
     mmsb1,mmsb2,mmu))/pow2(-mmgl + mmsb2) - (7*mmsb2*Fin20(mmsb1,mmsb2,mmu))/(
     3.*pow2(-mmgl + mmsb2)) + (12*mmsb2*Fin20(mmsb2,mmgl,mmu))/pow2(-mmgl +
     mmsb2) + (14*mmsb1*mmsb2*log(mmgl/mmu))/(3.*pow2(-mmgl + mmsb2)) + (14*
     mmsb2*mmst1*log(mmgl/mmu))/(3.*pow2(-mmgl + mmsb2)) + (14*mmsb2*mmst2*log(
     mmgl/mmu))/(3.*pow2(-mmgl + mmsb2)) + (112*mmsb2*mmsusy*log(mmgl/mmu))/(3.
     *pow2(-mmgl + mmsb2)) - (28*mmsb2*mmt*log(mmgl/mmu))/(3.*pow2(-mmgl +
     mmsb2)) + (16*mmsb1*mmsb2*log(mmsb1/mmu))/(3.*pow2(-mmgl + mmsb2)) - (7*
     mmsb1*mmsb2*log(mmgl/mmu)*log(mmsb1/mmu))/(3.*pow2(-mmgl + mmsb2)) - (2*
     mmsb1*mmsb2*log(mmsb2/mmu))/pow2(-mmgl + mmsb2) - (2*mmsb2*mmst1*log(
     mmsb2/mmu))/pow2(-mmgl + mmsb2) - (2*mmsb2*mmst2*log(mmsb2/mmu))/pow2(-
     mmgl + mmsb2) - (16*mmsb2*mmsusy*log(mmsb2/mmu))/pow2(-mmgl + mmsb2) + (4*
     mmsb2*mmt*log(mmsb2/mmu))/pow2(-mmgl + mmsb2) + (1222*pow2(mmsb2))/(9.*
     pow2(-mmgl + mmsb2)) + (20*zt2*pow2(mmsb2))/pow2(-mmgl + mmsb2) + (112*
     mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2))/(3.*pow2(-mmgl + mmsb2))
     + (16*mmst1*mmt*zt2*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2))/(3.*pow2(-mmgl
     + mmsb2)) + (112*mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2))/(3.*
     pow2(-mmgl + mmsb2)) + (16*mmst2*mmt*zt2*DeltaInv(mmt,mmsb2,mmst2)*pow2(
     mmsb2))/(3.*pow2(-mmgl + mmsb2)) + (112*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)
     *pow2(mmsb2))/(3.*pow2(-mmgl + mmsb2)) + (16*mmst1*mmt*zt2*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmsb2))/(3.*pow2(-mmgl + mmsb2)) + (112*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/(3.*pow2(-mmgl + mmsb2)) + (16*
     mmst2*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/(3.*pow2(-mmgl +
     mmsb2)) - (1100*log(mmgl/mmu)*pow2(mmsb2))/(9.*pow2(-mmgl + mmsb2)) - (4*
     mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*pow2(mmsb2))/pow2(-mmgl +
     mmsb2) - (4*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*pow2(mmsb2))/
     pow2(-mmgl + mmsb2) + (4*log(mmgl/mmu)*log(mmsb1/mmu)*pow2(mmsb2))/(9.*
     pow2(-mmgl + mmsb2)) + (52*log(mmsb2/mmu)*pow2(mmsb2))/(3.*pow2(-mmgl +
     mmsb2)) - (4*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1)*log(mmsb2/mmu)*pow2(
     mmsb2))/pow2(-mmgl + mmsb2) - (4*mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2)*log(
     mmsb2/mmu)*pow2(mmsb2))/pow2(-mmgl + mmsb2) - (109*log(mmgl/mmu)*log(
     mmsb2/mmu)*pow2(mmsb2))/(9.*pow2(-mmgl + mmsb2)) - (4*log(mmsb1/mmu)*log(
     mmsb2/mmu)*pow2(mmsb2))/(9.*pow2(-mmgl + mmsb2)) - (14*mmsb1*DeltaInv(mmt,
     mmsb1,mmst1)*pow2(mmst1))/(3.*(-mmgl + mmsb1)) - (14*mmt*DeltaInv(mmt,
     mmsb1,mmst1)*pow2(mmst1))/(3.*(-mmgl + mmsb1)) - (2*mmsb1*zt2*DeltaInv(
     mmt,mmsb1,mmst1)*pow2(mmst1))/(3.*(-mmgl + mmsb1)) - (2*mmt*zt2*DeltaInv(
     mmt,mmsb1,mmst1)*pow2(mmst1))/(3.*(-mmgl + mmsb1)) - (14*mmsb2*DeltaInv(
     mmt,mmsb2,mmst1)*pow2(mmst1))/(3.*(-mmgl + mmsb2)) - (14*mmt*DeltaInv(mmt,
     mmsb2,mmst1)*pow2(mmst1))/(3.*(-mmgl + mmsb2)) - (2*mmsb2*zt2*DeltaInv(
     mmt,mmsb2,mmst1)*pow2(mmst1))/(3.*(-mmgl + mmsb2)) - (2*mmt*zt2*DeltaInv(
     mmt,mmsb2,mmst1)*pow2(mmst1))/(3.*(-mmgl + mmsb2)) + (56*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmst1))/3. - (56*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmst1))/(3.*(-mmgl + mmsb1)) - (56*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmst1))/(3.*(-mmgl + mmsb2)) - (28*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmst1))/(3.*(-mmgl + mmsb1)) - (28*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmst1))/(3.*(-mmgl + mmsb2)) + (8*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)
     )/3. - (8*mmsb1*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(3.*(-mmgl +
     mmsb1)) - (8*mmsb2*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(3.*(-mmgl +
     mmsb2)) - (4*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(3.*(-mmgl +
     mmsb1)) - (4*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(3.*(-mmgl +
     mmsb2)) + 8*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*pow2(mmst1) - (8*mmsb1*
     DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*pow2(mmst1))/(-mmgl + mmsb1) - (8*
     mmsb2*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*pow2(mmst1))/(-mmgl + mmsb2)
     - (2*mmsb1*DeltaInv(mmt,mmsb1,mmst1)*log(mmsb1/mmu)*pow2(mmst1))/(-mmgl +
     mmsb1) - (2*mmsb2*DeltaInv(mmt,mmsb2,mmst1)*log(mmsb2/mmu)*pow2(mmst1))/(-
     mmgl + mmsb2) + (28*mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1))/(3.*
     pow2(-mmgl + mmsb1)) + (4*mmsb1*mmt*zt2*DeltaInv(mmt,mmsb1,mmst1)*pow2(
     mmst1))/(3.*pow2(-mmgl + mmsb1)) + (28*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmst1))/(3.*pow2(-mmgl + mmsb1)) + (4*mmsb1*mmt*zt2*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmst1))/(3.*pow2(-mmgl + mmsb1)) + (28*DeltaInv(mmt,
     mmsb1,mmst1)*pow2(mmsb1)*pow2(mmst1))/(3.*pow2(-mmgl + mmsb1)) + (4*zt2*
     DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow2(mmst1))/(3.*pow2(-mmgl + mmsb1)
     ) + (28*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(mmst1))/(3.*pow2(-mmgl +
     mmsb1)) + (4*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(mmst1))/(3.*
     pow2(-mmgl + mmsb1)) + (4*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*pow2(
     mmsb1)*pow2(mmst1))/pow2(-mmgl + mmsb1) + (4*DeltaInv(mmt,mmsb1,mmst1)*
     log(mmsb1/mmu)*pow2(mmsb1)*pow2(mmst1))/pow2(-mmgl + mmsb1) + (28*mmsb2*
     mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/(3.*pow2(-mmgl + mmsb2)) + (4*
     mmsb2*mmt*zt2*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/(3.*pow2(-mmgl +
     mmsb2)) + (28*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(3.*pow2(-
     mmgl + mmsb2)) + (4*mmsb2*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(
     3.*pow2(-mmgl + mmsb2)) + (28*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2)*pow2(
     mmst1))/(3.*pow2(-mmgl + mmsb2)) + (4*zt2*DeltaInv(mmt,mmsb2,mmst1)*pow2(
     mmsb2)*pow2(mmst1))/(3.*pow2(-mmgl + mmsb2)) + (28*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmsb2)*pow2(mmst1))/(3.*pow2(-mmgl + mmsb2)) + (4*zt2*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(mmst1))/(3.*pow2(-mmgl + mmsb2)) + (4*
     DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*pow2(mmsb2)*pow2(mmst1))/pow2(-mmgl
     + mmsb2) + (4*DeltaInv(mmt,mmsb2,mmst1)*log(mmsb2/mmu)*pow2(mmsb2)*pow2(
     mmst1))/pow2(-mmgl + mmsb2) - (14*mmsb1*DeltaInv(mmt,mmsb1,mmst2)*pow2(
     mmst2))/(3.*(-mmgl + mmsb1)) - (14*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(
     mmst2))/(3.*(-mmgl + mmsb1)) - (2*mmsb1*zt2*DeltaInv(mmt,mmsb1,mmst2)*
     pow2(mmst2))/(3.*(-mmgl + mmsb1)) - (2*mmt*zt2*DeltaInv(mmt,mmsb1,mmst2)*
     pow2(mmst2))/(3.*(-mmgl + mmsb1)) - (14*mmsb2*DeltaInv(mmt,mmsb2,mmst2)*
     pow2(mmst2))/(3.*(-mmgl + mmsb2)) - (14*mmt*DeltaInv(mmt,mmsb2,mmst2)*
     pow2(mmst2))/(3.*(-mmgl + mmsb2)) - (2*mmsb2*zt2*DeltaInv(mmt,mmsb2,mmst2)
     *pow2(mmst2))/(3.*(-mmgl + mmsb2)) - (2*mmt*zt2*DeltaInv(mmt,mmsb2,mmst2)*
     pow2(mmst2))/(3.*(-mmgl + mmsb2)) + (56*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmst2))/3. - (56*mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(3.*(-mmgl +
     mmsb1)) - (56*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(3.*(-mmgl +
     mmsb2)) - (28*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(3.*(-mmgl +
     mmsb1)) - (28*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(3.*(-mmgl +
     mmsb2)) + (8*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/3. - (8*mmsb1*zt2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(3.*(-mmgl + mmsb1)) - (8*mmsb2*zt2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(3.*(-mmgl + mmsb2)) - (4*mmt*zt2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(3.*(-mmgl + mmsb1)) - (4*mmt*zt2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(3.*(-mmgl + mmsb2)) + 8*DeltaInv(
     mmt,mmst2,mmgl)*log(mmgl/mmu)*pow2(mmst2) - (8*mmsb1*DeltaInv(mmt,mmst2,
     mmgl)*log(mmgl/mmu)*pow2(mmst2))/(-mmgl + mmsb1) - (8*mmsb2*DeltaInv(mmt,
     mmst2,mmgl)*log(mmgl/mmu)*pow2(mmst2))/(-mmgl + mmsb2) - (2*mmsb1*
     DeltaInv(mmt,mmsb1,mmst2)*log(mmsb1/mmu)*pow2(mmst2))/(-mmgl + mmsb1) - (
     2*mmsb2*DeltaInv(mmt,mmsb2,mmst2)*log(mmsb2/mmu)*pow2(mmst2))/(-mmgl +
     mmsb2) + (28*mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2))/(3.*pow2(-
     mmgl + mmsb1)) + (4*mmsb1*mmt*zt2*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2))/(
     3.*pow2(-mmgl + mmsb1)) + (28*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmst2))/(3.*pow2(-mmgl + mmsb1)) + (4*mmsb1*mmt*zt2*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmst2))/(3.*pow2(-mmgl + mmsb1)) + (28*DeltaInv(mmt,mmsb1,
     mmst2)*pow2(mmsb1)*pow2(mmst2))/(3.*pow2(-mmgl + mmsb1)) + (4*zt2*
     DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*pow2(mmst2))/(3.*pow2(-mmgl + mmsb1)
     ) + (28*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(mmst2))/(3.*pow2(-mmgl +
     mmsb1)) + (4*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(mmst2))/(3.*
     pow2(-mmgl + mmsb1)) + (4*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*pow2(
     mmsb1)*pow2(mmst2))/pow2(-mmgl + mmsb1) + (4*DeltaInv(mmt,mmsb1,mmst2)*
     log(mmsb1/mmu)*pow2(mmsb1)*pow2(mmst2))/pow2(-mmgl + mmsb1) + (28*mmsb2*
     mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2))/(3.*pow2(-mmgl + mmsb2)) + (4*
     mmsb2*mmt*zt2*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2))/(3.*pow2(-mmgl +
     mmsb2)) + (28*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(3.*pow2(-
     mmgl + mmsb2)) + (4*mmsb2*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(
     3.*pow2(-mmgl + mmsb2)) + (28*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow2(
     mmst2))/(3.*pow2(-mmgl + mmsb2)) + (4*zt2*DeltaInv(mmt,mmsb2,mmst2)*pow2(
     mmsb2)*pow2(mmst2))/(3.*pow2(-mmgl + mmsb2)) + (28*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmsb2)*pow2(mmst2))/(3.*pow2(-mmgl + mmsb2)) + (4*zt2*DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(mmst2))/(3.*pow2(-mmgl + mmsb2)) + (4*
     DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*pow2(mmsb2)*pow2(mmst2))/pow2(-mmgl
     + mmsb2) + (4*DeltaInv(mmt,mmsb2,mmst2)*log(mmsb2/mmu)*pow2(mmsb2)*pow2(
     mmst2))/pow2(-mmgl + mmsb2) - (254*pow2(log(mmgl/mmu)))/9. + (85*mmsb1*
     pow2(log(mmgl/mmu)))/(18.*(-mmgl + mmsb1)) + (353*mmsb2*pow2(log(mmgl/mmu)
     ))/(18.*(-mmgl + mmsb1)) + (mmsb1*pow2(log(mmgl/mmu)))/(2.*(-mmgl + mmsb2)
     ) - (259*mmsb2*pow2(log(mmgl/mmu)))/(18.*(-mmgl + mmsb2)) + (mmst1*pow2(
     log(mmgl/mmu)))/(2.*(-mmgl + mmsb1)) + (mmst1*pow2(log(mmgl/mmu)))/(2.*(-
     mmgl + mmsb2)) + (mmst2*pow2(log(mmgl/mmu)))/(2.*(-mmgl + mmsb1)) + (
     mmst2*pow2(log(mmgl/mmu)))/(2.*(-mmgl + mmsb2)) + (4*mmsusy*pow2(log(mmgl/
     mmu)))/(-mmgl + mmsb1) + (4*mmsusy*pow2(log(mmgl/mmu)))/(-mmgl + mmsb2) -
     (mmt*pow2(log(mmgl/mmu)))/(-mmgl + mmsb1) - (mmt*pow2(log(mmgl/mmu)))/(-
     mmgl + mmsb2) - (4*mmgl*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmgl/mmu))
     )/3. - (4*mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmgl/mmu)))/3. + (
     8*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmgl/mmu)))/3. + (8*mmsb1*
     mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmgl/mmu)))/3. + (8*mmsb2*mmst1*
     DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmgl/mmu)))/3. + (4*mmgl*mmt*DeltaInv(
     mmt,mmst1,mmgl)*pow2(log(mmgl/mmu)))/3. + (4*mmsb1*mmt*DeltaInv(mmt,mmst1,
     mmgl)*pow2(log(mmgl/mmu)))/3. + (4*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow2(log(mmgl/mmu)))/3. + (4*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(log(
     mmgl/mmu)))/3. - (4*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(log(
     mmgl/mmu)))/(3.*(-mmgl + mmsb1)) - (4*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,
     mmgl)*pow2(log(mmgl/mmu)))/(3.*(-mmgl + mmsb2)) - (4*mmgl*mmsb1*DeltaInv(
     mmt,mmst2,mmgl)*pow2(log(mmgl/mmu)))/3. - (4*mmgl*mmsb2*DeltaInv(mmt,
     mmst2,mmgl)*pow2(log(mmgl/mmu)))/3. + (8*mmgl*mmst2*DeltaInv(mmt,mmst2,
     mmgl)*pow2(log(mmgl/mmu)))/3. + (8*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*
     pow2(log(mmgl/mmu)))/3. + (8*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(
     log(mmgl/mmu)))/3. + (4*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmgl/
     mmu)))/3. + (4*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmgl/mmu)))/3.
      + (4*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmgl/mmu)))/
     3. + (4*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmgl/mmu)))/3. - (4*
     mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmgl/mmu)))/(3.*(-mmgl +
     mmsb1)) - (4*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmgl/mmu)))
     /(3.*(-mmgl + mmsb2)) - (4*DeltaInv(mmt,mmst1,mmgl)*pow2(mmgl)*pow2(log(
     mmgl/mmu)))/3. - (4*DeltaInv(mmt,mmst2,mmgl)*pow2(mmgl)*pow2(log(mmgl/mmu)
     ))/3. - 2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(log(mmgl/mmu)) - (4*
     mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(log(mmgl/mmu)))/(-mmgl +
     mmsb1) - (2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(log(mmgl/mmu)))/
     (-mmgl + mmsb1) - 2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(log(mmgl/
     mmu)) - (4*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(log(mmgl/mmu)))
     /(-mmgl + mmsb1) - (2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(log(
     mmgl/mmu)))/(-mmgl + mmsb1) - (5*mmsb1*mmsb2*pow2(log(mmgl/mmu)))/(18.*
     pow2(-mmgl + mmsb1)) - (7*mmsb1*mmst1*pow2(log(mmgl/mmu)))/(6.*pow2(-mmgl
     + mmsb1)) - (7*mmsb1*mmst2*pow2(log(mmgl/mmu)))/(6.*pow2(-mmgl + mmsb1)) -
     (28*mmsb1*mmsusy*pow2(log(mmgl/mmu)))/(3.*pow2(-mmgl + mmsb1)) + (7*mmsb1*
     mmt*pow2(log(mmgl/mmu)))/(3.*pow2(-mmgl + mmsb1)) + (853*pow2(mmsb1)*pow2(
     log(mmgl/mmu)))/(18.*pow2(-mmgl + mmsb1)) + (2*mmst1*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmsb1)*pow2(log(mmgl/mmu)))/(3.*pow2(-mmgl + mmsb1)) + (
     2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(log(mmgl/mmu)))/(3.*
     pow2(-mmgl + mmsb1)) + (20*pow2(mmsb2)*pow2(log(mmgl/mmu)))/((-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) - (20*pow2(mmsb2)*pow2(log(mmgl/mmu)))/((mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) - 2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(log(
     mmgl/mmu)) - (4*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(log(mmgl/
     mmu)))/(-mmgl + mmsb2) - (2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(
     log(mmgl/mmu)))/(-mmgl + mmsb2) - 2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*
     pow2(log(mmgl/mmu)) - (4*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(
     log(mmgl/mmu)))/(-mmgl + mmsb2) - (2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmsb2)*pow2(log(mmgl/mmu)))/(-mmgl + mmsb2) + (4*pow2(mmsb2)*pow2(log(
     mmgl/mmu)))/(3.*pow2(-mmgl + mmsb1)) - (7*mmsb1*mmsb2*pow2(log(mmgl/mmu)))
     /(6.*pow2(-mmgl + mmsb2)) - (7*mmsb2*mmst1*pow2(log(mmgl/mmu)))/(6.*pow2(-
     mmgl + mmsb2)) - (7*mmsb2*mmst2*pow2(log(mmgl/mmu)))/(6.*pow2(-mmgl +
     mmsb2)) - (28*mmsb2*mmsusy*pow2(log(mmgl/mmu)))/(3.*pow2(-mmgl + mmsb2)) +
     (7*mmsb2*mmt*pow2(log(mmgl/mmu)))/(3.*pow2(-mmgl + mmsb2)) + (845*pow2(
     mmsb2)*pow2(log(mmgl/mmu)))/(18.*pow2(-mmgl + mmsb2)) + (2*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(log(mmgl/mmu)))/(3.*pow2(-mmgl +
     mmsb2)) + (2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(log(mmgl/
     mmu)))/(3.*pow2(-mmgl + mmsb2)) - (4*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*
     pow2(log(mmgl/mmu)))/3. + (4*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*
     pow2(log(mmgl/mmu)))/(3.*(-mmgl + mmsb1)) + (4*mmsb2*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmst1)*pow2(log(mmgl/mmu)))/(3.*(-mmgl + mmsb2)) - (2*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(mmst1)*pow2(log(mmgl/mmu)))/(3.*pow2(-
     mmgl + mmsb1)) - (2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(mmst1)*pow2(
     log(mmgl/mmu)))/(3.*pow2(-mmgl + mmsb2)) - (4*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmst2)*pow2(log(mmgl/mmu)))/3. + (4*mmsb1*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmst2)*pow2(log(mmgl/mmu)))/(3.*(-mmgl + mmsb1)) + (4*mmsb2*DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmst2)*pow2(log(mmgl/mmu)))/(3.*(-mmgl + mmsb2)) - (
     2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(mmst2)*pow2(log(mmgl/mmu)))/(
     3.*pow2(-mmgl + mmsb1)) - (2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(
     mmst2)*pow2(log(mmgl/mmu)))/(3.*pow2(-mmgl + mmsb2)) - (61*pow2(log(mmsb1/
     mmu)))/9. + (9*mmsb1*pow2(log(mmsb1/mmu)))/(2.*(-mmgl + mmsb1)) + (64*
     mmsb2*pow2(log(mmsb1/mmu)))/(9.*(-mmgl + mmsb1)) - (64*mmsb2*pow2(log(
     mmsb1/mmu)))/(9.*(mmsb1 - mmsb2)) + (mmsb1*pow2(log(mmsb1/mmu)))/(2.*(-
     mmgl + mmsb2)) - (mmsb1*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(log(
     mmsb1/mmu)))/(3.*(-mmgl + mmsb1)) - (mmsb1*mmst2*mmt*DeltaInv(mmt,mmsb1,
     mmst2)*pow2(log(mmsb1/mmu)))/(3.*(-mmgl + mmsb1)) - (2*mmst1*DeltaInv(mmt,
     mmsb1,mmst1)*pow2(mmsb1)*pow2(log(mmsb1/mmu)))/(3.*(-mmgl + mmsb1)) - (
     mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow2(log(mmsb1/mmu)))/(3.*(-mmgl
     + mmsb1)) - (2*mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*pow2(log(mmsb1/
     mmu)))/(3.*(-mmgl + mmsb1)) - (mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*
     pow2(log(mmsb1/mmu)))/(3.*(-mmgl + mmsb1)) + (mmsb1*mmsb2*pow2(log(mmsb1/
     mmu)))/(2.*pow2(-mmgl + mmsb1)) + (mmsb1*mmst1*pow2(log(mmsb1/mmu)))/(2.*
     pow2(-mmgl + mmsb1)) + (mmsb1*mmst2*pow2(log(mmsb1/mmu)))/(2.*pow2(-mmgl +
     mmsb1)) + (4*mmsb1*mmsusy*pow2(log(mmsb1/mmu)))/pow2(-mmgl + mmsb1) - (
     mmsb1*mmt*pow2(log(mmsb1/mmu)))/pow2(-mmgl + mmsb1) - (17*pow2(mmsb1)*
     pow2(log(mmsb1/mmu)))/(6.*pow2(-mmgl + mmsb1)) + (2*mmst1*mmt*DeltaInv(
     mmt,mmsb1,mmst1)*pow2(mmsb1)*pow2(log(mmsb1/mmu)))/(3.*pow2(-mmgl + mmsb1)
     ) + (2*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*pow2(log(mmsb1/mmu)
     ))/(3.*pow2(-mmgl + mmsb1)) + (64*pow2(mmsb2)*pow2(log(mmsb1/mmu)))/(9.*(-
     mmgl + mmsb1)*(mmsb1 - mmsb2)) - (2*mmsb1*mmsb2*pow2(log(mmsb1/mmu)))/(3.*
     pow2(-mmgl + mmsb2)) + (mmsb1*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1)*pow2(
     log(mmsb1/mmu)))/(3.*(-mmgl + mmsb1)) - (2*DeltaInv(mmt,mmsb1,mmst1)*pow2(
     mmsb1)*pow2(mmst1)*pow2(log(mmsb1/mmu)))/(3.*pow2(-mmgl + mmsb1)) + (
     mmsb1*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2)*pow2(log(mmsb1/mmu)))/(3.*(-
     mmgl + mmsb1)) - (2*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*pow2(mmst2)*
     pow2(log(mmsb1/mmu)))/(3.*pow2(-mmgl + mmsb1)) + pow2(log(mmsb2/mmu))/3. +
     (mmsb2*pow2(log(mmsb2/mmu)))/(2.*(-mmgl + mmsb1)) + (64*mmsb2*pow2(log(
     mmsb2/mmu)))/(9.*(mmsb1 - mmsb2)) - (47*mmsb2*pow2(log(mmsb2/mmu)))/(18.*(
     -mmgl + mmsb2)) - (mmsb2*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(log(
     mmsb2/mmu)))/(3.*(-mmgl + mmsb2)) - (mmsb2*mmst2*mmt*DeltaInv(mmt,mmsb2,
     mmst2)*pow2(log(mmsb2/mmu)))/(3.*(-mmgl + mmsb2)) - (2*mmsb1*mmsb2*pow2(
     log(mmsb2/mmu)))/(3.*pow2(-mmgl + mmsb1)) - (64*pow2(mmsb2)*pow2(log(
     mmsb2/mmu)))/(9.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (2*mmst1*DeltaInv(mmt,
     mmsb2,mmst1)*pow2(mmsb2)*pow2(log(mmsb2/mmu)))/(3.*(-mmgl + mmsb2)) - (
     mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2)*pow2(log(mmsb2/mmu)))/(3.*(-mmgl
     + mmsb2)) - (2*mmst2*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow2(log(mmsb2/
     mmu)))/(3.*(-mmgl + mmsb2)) - (mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*
     pow2(log(mmsb2/mmu)))/(3.*(-mmgl + mmsb2)) + (mmsb1*mmsb2*pow2(log(mmsb2/
     mmu)))/(2.*pow2(-mmgl + mmsb2)) + (mmsb2*mmst1*pow2(log(mmsb2/mmu)))/(2.*
     pow2(-mmgl + mmsb2)) + (mmsb2*mmst2*pow2(log(mmsb2/mmu)))/(2.*pow2(-mmgl +
     mmsb2)) + (4*mmsb2*mmsusy*pow2(log(mmsb2/mmu)))/pow2(-mmgl + mmsb2) - (
     mmsb2*mmt*pow2(log(mmsb2/mmu)))/pow2(-mmgl + mmsb2) - (17*pow2(mmsb2)*
     pow2(log(mmsb2/mmu)))/(6.*pow2(-mmgl + mmsb2)) + (2*mmst1*mmt*DeltaInv(
     mmt,mmsb2,mmst1)*pow2(mmsb2)*pow2(log(mmsb2/mmu)))/(3.*pow2(-mmgl + mmsb2)
     ) + (2*mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow2(log(mmsb2/mmu)
     ))/(3.*pow2(-mmgl + mmsb2)) + (mmsb2*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1)
     *pow2(log(mmsb2/mmu)))/(3.*(-mmgl + mmsb2)) - (2*DeltaInv(mmt,mmsb2,mmst1)
     *pow2(mmsb2)*pow2(mmst1)*pow2(log(mmsb2/mmu)))/(3.*pow2(-mmgl + mmsb2)) +
     (mmsb2*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2)*pow2(log(mmsb2/mmu)))/(3.*(-
     mmgl + mmsb2)) - (2*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow2(mmst2)*
     pow2(log(mmsb2/mmu)))/(3.*pow2(-mmgl + mmsb2)) + (4*(2 + (3*mmsusy)/(-mmgl
     + mmsb1) + (3*mmsusy)/(-mmgl + mmsb2) - (4*mmsb1*mmsusy)/pow2(-mmgl +
     mmsb1) - (4*mmsb2*mmsusy)/pow2(-mmgl + mmsb2))*pow2(log(mmsusy/mmu)))/3. +
     (14*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1))/(3.*(-mmgl + mmsb1)) + (2*zt2*
     DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1))/(3.*(-mmgl + mmsb1)) + (14*
     DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1))/(3.*(-mmgl + mmsb1)) + (2*zt2*
     DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1))/(3.*(-mmgl + mmsb1)) + (112*
     DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1))/(3.*(-mmgl + mmsb1)) + (16*zt2*
     DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1))/(3.*(-mmgl + mmsb1)) + (112*
     DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1))/(3.*(-mmgl + mmsb1)) + (16*zt2*
     DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1))/(3.*(-mmgl + mmsb1)) - (16*DeltaInv(
     mmt,mmst1,mmgl)*log(mmgl/mmu)*pow3(mmsb1))/(-mmgl + mmsb1) - (16*DeltaInv(
     mmt,mmst2,mmgl)*log(mmgl/mmu)*pow3(mmsb1))/(-mmgl + mmsb1) - (2*DeltaInv(
     mmt,mmsb1,mmst1)*log(mmsb1/mmu)*pow3(mmsb1))/(-mmgl + mmsb1) - (2*
     DeltaInv(mmt,mmsb1,mmst2)*log(mmsb1/mmu)*pow3(mmsb1))/(-mmgl + mmsb1) + (
     28*mmst1*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1))/(3.*pow2(-mmgl + mmsb1)) +
     (28*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1))/(3.*pow2(-mmgl + mmsb1)) +
     (4*mmst1*zt2*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1))/(3.*pow2(-mmgl +
     mmsb1)) + (4*mmt*zt2*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1))/(3.*pow2(-mmgl
     + mmsb1)) + (28*mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1))/(3.*pow2(-
     mmgl + mmsb1)) + (28*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1))/(3.*pow2(-
     mmgl + mmsb1)) + (4*mmst2*zt2*DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1))/(3.*
     pow2(-mmgl + mmsb1)) + (4*mmt*zt2*DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1))/(
     3.*pow2(-mmgl + mmsb1)) + (28*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1))/
     (3.*pow2(-mmgl + mmsb1)) + (28*mmt*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1))/(
     3.*pow2(-mmgl + mmsb1)) + (4*mmst1*zt2*DeltaInv(mmt,mmst1,mmgl)*pow3(
     mmsb1))/(3.*pow2(-mmgl + mmsb1)) + (4*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*
     pow3(mmsb1))/(3.*pow2(-mmgl + mmsb1)) + (28*mmst2*DeltaInv(mmt,mmst2,mmgl)
     *pow3(mmsb1))/(3.*pow2(-mmgl + mmsb1)) + (28*mmt*DeltaInv(mmt,mmst2,mmgl)*
     pow3(mmsb1))/(3.*pow2(-mmgl + mmsb1)) + (4*mmst2*zt2*DeltaInv(mmt,mmst2,
     mmgl)*pow3(mmsb1))/(3.*pow2(-mmgl + mmsb1)) + (4*mmt*zt2*DeltaInv(mmt,
     mmst2,mmgl)*pow3(mmsb1))/(3.*pow2(-mmgl + mmsb1)) - (8*mmst1*DeltaInv(mmt,
     mmst1,mmgl)*log(mmgl/mmu)*pow3(mmsb1))/pow2(-mmgl + mmsb1) - (4*mmt*
     DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*pow3(mmsb1))/pow2(-mmgl + mmsb1) -
     (8*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*pow3(mmsb1))/pow2(-mmgl +
     mmsb1) - (4*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*pow3(mmsb1))/pow2(-
     mmgl + mmsb1) - (8*mmst1*DeltaInv(mmt,mmsb1,mmst1)*log(mmsb1/mmu)*pow3(
     mmsb1))/pow2(-mmgl + mmsb1) - (4*mmt*DeltaInv(mmt,mmsb1,mmst1)*log(mmsb1/
     mmu)*pow3(mmsb1))/pow2(-mmgl + mmsb1) - (8*mmst2*DeltaInv(mmt,mmsb1,mmst2)
     *log(mmsb1/mmu)*pow3(mmsb1))/pow2(-mmgl + mmsb1) - (4*mmt*DeltaInv(mmt,
     mmsb1,mmst2)*log(mmsb1/mmu)*pow3(mmsb1))/pow2(-mmgl + mmsb1) + (8*
     DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmgl/mmu))*pow3(mmsb1))/(3.*(-mmgl +
     mmsb1)) + (8*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmgl/mmu))*pow3(mmsb1))/(3.
     *(-mmgl + mmsb1)) + (4*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmgl/mmu))*
     pow3(mmsb1))/(3.*pow2(-mmgl + mmsb1)) + (2*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow2(log(mmgl/mmu))*pow3(mmsb1))/(3.*pow2(-mmgl + mmsb1)) + (4*mmst2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmgl/mmu))*pow3(mmsb1))/(3.*pow2(-mmgl +
     mmsb1)) + (2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmgl/mmu))*pow3(mmsb1))
     /(3.*pow2(-mmgl + mmsb1)) + (DeltaInv(mmt,mmsb1,mmst1)*pow2(log(mmsb1/mmu)
     )*pow3(mmsb1))/(3.*(-mmgl + mmsb1)) + (DeltaInv(mmt,mmsb1,mmst2)*pow2(log(
     mmsb1/mmu))*pow3(mmsb1))/(3.*(-mmgl + mmsb1)) + (4*mmst1*DeltaInv(mmt,
     mmsb1,mmst1)*pow2(log(mmsb1/mmu))*pow3(mmsb1))/(3.*pow2(-mmgl + mmsb1)) +
     (2*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(log(mmsb1/mmu))*pow3(mmsb1))/(3.*
     pow2(-mmgl + mmsb1)) + (4*mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow2(log(mmsb1/
     mmu))*pow3(mmsb1))/(3.*pow2(-mmgl + mmsb1)) + (2*mmt*DeltaInv(mmt,mmsb1,
     mmst2)*pow2(log(mmsb1/mmu))*pow3(mmsb1))/(3.*pow2(-mmgl + mmsb1)) + (4*
     Fin3(mmt,mmsb1,mmst1,mmu)*(-((mmsb1*mmst1*DeltaInv(mmt,mmsb1,mmst1))/(-
     mmgl + mmsb1)) - (mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst1))/(2.*(-mmgl + mmsb1)
     ) - (mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1))/(2.*(-mmgl + mmsb1)) + (
     DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1))/(2.*(-mmgl + mmsb1)) - (3*mmsb1)/(
     4.*pow2(-mmgl + mmsb1)) + (3*mmst1)/(4.*pow2(-mmgl + mmsb1)) - (3*mmt)/(4.
     *pow2(-mmgl + mmsb1)) + (mmsb1*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1))/pow2(-
     mmgl + mmsb1) + (2*mmst1*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1))/pow2(-mmgl
     + mmsb1) + (mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1))/pow2(-mmgl + mmsb1)
     + (DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1))/(2.*(-mmgl + mmsb1)) - (mmsb1*
     DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1))/pow2(-mmgl + mmsb1) - (DeltaInv(
     mmt,mmsb1,mmst1)*pow3(mmsb1))/pow2(-mmgl + mmsb1) - (mmsb1*mmst1)/pow3(-
     mmgl + mmsb1) + (mmsb1*mmt)/pow3(-mmgl + mmsb1) + pow2(mmsb1)/pow3(-mmgl +
     mmsb1)))/3. + (4*Fin3(mmt,mmsb1,mmst2,mmu)*(-((mmsb1*mmst2*DeltaInv(mmt,
     mmsb1,mmst2))/(-mmgl + mmsb1)) - (mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst2))/(2.
     *(-mmgl + mmsb1)) - (mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2))/(2.*(-mmgl +
     mmsb1)) + (DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))/(2.*(-mmgl + mmsb1)) - (
     3*mmsb1)/(4.*pow2(-mmgl + mmsb1)) + (3*mmst2)/(4.*pow2(-mmgl + mmsb1)) - (
     3*mmt)/(4.*pow2(-mmgl + mmsb1)) + (mmsb1*mmst2*mmt*DeltaInv(mmt,mmsb1,
     mmst2))/pow2(-mmgl + mmsb1) + (2*mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow2(
     mmsb1))/pow2(-mmgl + mmsb1) + (mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))/
     pow2(-mmgl + mmsb1) + (DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2))/(2.*(-mmgl +
     mmsb1)) - (mmsb1*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2))/pow2(-mmgl +
     mmsb1) - (DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1))/pow2(-mmgl + mmsb1) - (
     mmsb1*mmst2)/pow3(-mmgl + mmsb1) + (mmsb1*mmt)/pow3(-mmgl + mmsb1) + pow2(
     mmsb1)/pow3(-mmgl + mmsb1)))/3. + (4*Fin20(mmsb1,mmsusy,mmu)*(4/(-mmgl +
     mmsb1) - (14*mmsb1)/pow2(-mmgl + mmsb1) + (6*mmsusy)/pow2(-mmgl + mmsb1) -
     (8*mmsb1*mmsusy)/pow3(-mmgl + mmsb1) + (8*pow2(mmsb1))/pow3(-mmgl + mmsb1)
     ))/3. + (4*log(mmsb1/mmu)*log(mmst1/mmu)*((mmsb1*mmst1*mmt*DeltaInv(mmt,
     mmsb1,mmst1))/(-mmgl + mmsb1) + (3*mmsb1*mmst1)/(4.*pow2(-mmgl + mmsb1)) -
     (2*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1))/pow2(-mmgl + mmsb1) -
     (mmst1*pow2(mmsb1))/pow3(-mmgl + mmsb1)))/3. + (4*log(mmsb1/mmu)*log(
     mmst2/mmu)*((mmsb1*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2))/(-mmgl + mmsb1) +
     (3*mmsb1*mmst2)/(4.*pow2(-mmgl + mmsb1)) - (2*mmst2*mmt*DeltaInv(mmt,
     mmsb1,mmst2)*pow2(mmsb1))/pow2(-mmgl + mmsb1) - (mmst2*pow2(mmsb1))/pow3(-
     mmgl + mmsb1)))/3. + (4*log(mmsb1/mmu)*log(mmsusy/mmu)*((6*mmsb1*mmsusy)/
     pow2(-mmgl + mmsb1) - (8*mmsusy*pow2(mmsb1))/pow3(-mmgl + mmsb1)))/3. - (
     4*mmsb1*mmsb2*Fin20(mmsb1,mmsb2,mmu))/(3.*pow3(-mmgl + mmsb1)) + (4*mmsb1*
     mmsb2*Fin20(mmsb2,mmgl,mmu))/(3.*pow3(-mmgl + mmsb1)) - (16*Fin20(mmsb1,
     mmgl,mmu)*pow2(mmsb1))/(3.*pow3(-mmgl + mmsb1)) + (4*Fin20(mmsb1,mmsb2,
     mmu)*pow2(mmsb1))/(3.*pow3(-mmgl + mmsb1)) - (4*Fin20(mmsb2,mmgl,mmu)*
     pow2(mmsb1))/(3.*pow3(-mmgl + mmsb1)) - (8*mmsb2*log(mmgl/mmu)*pow2(mmsb1)
     )/(3.*pow3(-mmgl + mmsb1)) - (8*mmst1*log(mmgl/mmu)*pow2(mmsb1))/(3.*pow3(
     -mmgl + mmsb1)) - (8*mmst2*log(mmgl/mmu)*pow2(mmsb1))/(3.*pow3(-mmgl +
     mmsb1)) - (64*mmsusy*log(mmgl/mmu)*pow2(mmsb1))/(3.*pow3(-mmgl + mmsb1)) +
     (16*mmt*log(mmgl/mmu)*pow2(mmsb1))/(3.*pow3(-mmgl + mmsb1)) + (8*mmsb2*
     log(mmsb1/mmu)*pow2(mmsb1))/(3.*pow3(-mmgl + mmsb1)) + (8*mmst1*log(mmsb1/
     mmu)*pow2(mmsb1))/(3.*pow3(-mmgl + mmsb1)) + (8*mmst2*log(mmsb1/mmu)*pow2(
     mmsb1))/(3.*pow3(-mmgl + mmsb1)) + (64*mmsusy*log(mmsb1/mmu)*pow2(mmsb1))/
     (3.*pow3(-mmgl + mmsb1)) - (16*mmt*log(mmsb1/mmu)*pow2(mmsb1))/(3.*pow3(-
     mmgl + mmsb1)) + (4*mmsb2*log(mmgl/mmu)*log(mmsb2/mmu)*pow2(mmsb1))/(3.*
     pow3(-mmgl + mmsb1)) + (2*mmsb2*pow2(mmsb1)*pow2(log(mmgl/mmu)))/(3.*pow3(
     -mmgl + mmsb1)) + (2*mmst1*pow2(mmsb1)*pow2(log(mmgl/mmu)))/(3.*pow3(-mmgl
     + mmsb1)) + (2*mmst2*pow2(mmsb1)*pow2(log(mmgl/mmu)))/(3.*pow3(-mmgl +
     mmsb1)) + (16*mmsusy*pow2(mmsb1)*pow2(log(mmgl/mmu)))/(3.*pow3(-mmgl +
     mmsb1)) - (4*mmt*pow2(mmsb1)*pow2(log(mmgl/mmu)))/(3.*pow3(-mmgl + mmsb1))
     - (2*mmsb2*pow2(mmsb1)*pow2(log(mmsb1/mmu)))/(3.*pow3(-mmgl + mmsb1)) - (
     2*mmst1*pow2(mmsb1)*pow2(log(mmsb1/mmu)))/(3.*pow3(-mmgl + mmsb1)) - (2*
     mmst2*pow2(mmsb1)*pow2(log(mmsb1/mmu)))/(3.*pow3(-mmgl + mmsb1)) - (16*
     mmsusy*pow2(mmsb1)*pow2(log(mmsb1/mmu)))/(3.*pow3(-mmgl + mmsb1)) + (4*
     mmt*pow2(mmsb1)*pow2(log(mmsb1/mmu)))/(3.*pow3(-mmgl + mmsb1)) - (112*
     pow3(mmsb1))/(3.*pow3(-mmgl + mmsb1)) - (16*zt2*pow3(mmsb1))/(3.*pow3(-
     mmgl + mmsb1)) + (508*log(mmgl/mmu)*pow3(mmsb1))/(9.*pow3(-mmgl + mmsb1))
     - (220*log(mmsb1/mmu)*pow3(mmsb1))/(9.*pow3(-mmgl + mmsb1)) + (116*log(
     mmgl/mmu)*log(mmsb1/mmu)*pow3(mmsb1))/(9.*pow3(-mmgl + mmsb1)) - (242*
     pow2(log(mmgl/mmu))*pow3(mmsb1))/(9.*pow3(-mmgl + mmsb1)) + (10*pow2(log(
     mmsb1/mmu))*pow3(mmsb1))/(3.*pow3(-mmgl + mmsb1)) + (14*DeltaInv(mmt,
     mmsb2,mmst1)*pow3(mmsb2))/(3.*(-mmgl + mmsb2)) + (2*zt2*DeltaInv(mmt,
     mmsb2,mmst1)*pow3(mmsb2))/(3.*(-mmgl + mmsb2)) + (14*DeltaInv(mmt,mmsb2,
     mmst2)*pow3(mmsb2))/(3.*(-mmgl + mmsb2)) + (2*zt2*DeltaInv(mmt,mmsb2,
     mmst2)*pow3(mmsb2))/(3.*(-mmgl + mmsb2)) + (112*DeltaInv(mmt,mmst1,mmgl)*
     pow3(mmsb2))/(3.*(-mmgl + mmsb2)) + (16*zt2*DeltaInv(mmt,mmst1,mmgl)*pow3(
     mmsb2))/(3.*(-mmgl + mmsb2)) + (112*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2))/
     (3.*(-mmgl + mmsb2)) + (16*zt2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2))/(3.*(
     -mmgl + mmsb2)) - (16*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*pow3(mmsb2))/
     (-mmgl + mmsb2) - (16*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*pow3(mmsb2))/
     (-mmgl + mmsb2) - (2*DeltaInv(mmt,mmsb2,mmst1)*log(mmsb2/mmu)*pow3(mmsb2))
     /(-mmgl + mmsb2) - (2*DeltaInv(mmt,mmsb2,mmst2)*log(mmsb2/mmu)*pow3(mmsb2)
     )/(-mmgl + mmsb2) - (4*log(mmgl/mmu)*pow3(mmsb2))/(9.*(mmsb1 - mmsb2)*
     pow2(-mmgl + mmsb1)) + (4*log(mmsb1/mmu)*pow3(mmsb2))/(9.*(mmsb1 - mmsb2)*
     pow2(-mmgl + mmsb1)) - (16*log(mmgl/mmu)*log(mmsb1/mmu)*pow3(mmsb2))/(9.*(
     mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (16*log(mmgl/mmu)*log(mmsb2/mmu)*
     pow3(mmsb2))/(9.*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) + (16*log(mmsb1/mmu)
     *log(mmsb2/mmu)*pow3(mmsb2))/(9.*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) + (
     4*log(mmsb1/mmu)*pow3(mmsb2))/(9.*(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (
     4*log(mmsb1/mmu)*pow3(mmsb2))/(9.*(-mmgl + mmsb2)*pow2(mmsb1 - mmsb2)) - (
     16*log(mmgl/mmu)*log(mmsb1/mmu)*pow3(mmsb2))/(9.*(-mmgl + mmsb1)*pow2(
     mmsb1 - mmsb2)) + (16*log(mmgl/mmu)*log(mmsb1/mmu)*pow3(mmsb2))/(9.*(-mmgl
     + mmsb2)*pow2(mmsb1 - mmsb2)) - (4*log(mmsb2/mmu)*pow3(mmsb2))/(9.*(-mmgl
     + mmsb1)*pow2(mmsb1 - mmsb2)) + (4*log(mmsb2/mmu)*pow3(mmsb2))/(9.*(-mmgl
     + mmsb2)*pow2(mmsb1 - mmsb2)) - (16*log(mmgl/mmu)*log(mmsb2/mmu)*pow3(
     mmsb2))/(9.*(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) + (16*log(mmgl/mmu)*log(
     mmsb2/mmu)*pow3(mmsb2))/(9.*(-mmgl + mmsb2)*pow2(mmsb1 - mmsb2)) + (16*
     log(mmsb1/mmu)*log(mmsb2/mmu)*pow3(mmsb2))/(9.*(-mmgl + mmsb1)*pow2(mmsb1
     - mmsb2)) - (16*log(mmsb1/mmu)*log(mmsb2/mmu)*pow3(mmsb2))/(9.*(-mmgl +
     mmsb2)*pow2(mmsb1 - mmsb2)) + (28*mmst1*DeltaInv(mmt,mmsb2,mmst1)*pow3(
     mmsb2))/(3.*pow2(-mmgl + mmsb2)) + (28*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow3(
     mmsb2))/(3.*pow2(-mmgl + mmsb2)) + (4*mmst1*zt2*DeltaInv(mmt,mmsb2,mmst1)*
     pow3(mmsb2))/(3.*pow2(-mmgl + mmsb2)) + (4*mmt*zt2*DeltaInv(mmt,mmsb2,
     mmst1)*pow3(mmsb2))/(3.*pow2(-mmgl + mmsb2)) + (28*mmst2*DeltaInv(mmt,
     mmsb2,mmst2)*pow3(mmsb2))/(3.*pow2(-mmgl + mmsb2)) + (28*mmt*DeltaInv(mmt,
     mmsb2,mmst2)*pow3(mmsb2))/(3.*pow2(-mmgl + mmsb2)) + (4*mmst2*zt2*
     DeltaInv(mmt,mmsb2,mmst2)*pow3(mmsb2))/(3.*pow2(-mmgl + mmsb2)) + (4*mmt*
     zt2*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmsb2))/(3.*pow2(-mmgl + mmsb2)) + (28*
     mmst1*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2))/(3.*pow2(-mmgl + mmsb2)) + (
     28*mmt*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2))/(3.*pow2(-mmgl + mmsb2)) + (
     4*mmst1*zt2*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2))/(3.*pow2(-mmgl + mmsb2))
     + (4*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2))/(3.*pow2(-mmgl + mmsb2)
     ) + (28*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2))/(3.*pow2(-mmgl +
     mmsb2)) + (28*mmt*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2))/(3.*pow2(-mmgl +
     mmsb2)) + (4*mmst2*zt2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2))/(3.*pow2(-
     mmgl + mmsb2)) + (4*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2))/(3.*
     pow2(-mmgl + mmsb2)) + (4*log(mmgl/mmu)*pow3(mmsb2))/(9.*(mmsb1 - mmsb2)*
     pow2(-mmgl + mmsb2)) - (8*mmst1*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*
     pow3(mmsb2))/pow2(-mmgl + mmsb2) - (4*mmt*DeltaInv(mmt,mmst1,mmgl)*log(
     mmgl/mmu)*pow3(mmsb2))/pow2(-mmgl + mmsb2) - (8*mmst2*DeltaInv(mmt,mmst2,
     mmgl)*log(mmgl/mmu)*pow3(mmsb2))/pow2(-mmgl + mmsb2) - (4*mmt*DeltaInv(
     mmt,mmst2,mmgl)*log(mmgl/mmu)*pow3(mmsb2))/pow2(-mmgl + mmsb2) - (4*log(
     mmsb2/mmu)*pow3(mmsb2))/(9.*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (8*
     mmst1*DeltaInv(mmt,mmsb2,mmst1)*log(mmsb2/mmu)*pow3(mmsb2))/pow2(-mmgl +
     mmsb2) - (4*mmt*DeltaInv(mmt,mmsb2,mmst1)*log(mmsb2/mmu)*pow3(mmsb2))/
     pow2(-mmgl + mmsb2) - (8*mmst2*DeltaInv(mmt,mmsb2,mmst2)*log(mmsb2/mmu)*
     pow3(mmsb2))/pow2(-mmgl + mmsb2) - (4*mmt*DeltaInv(mmt,mmsb2,mmst2)*log(
     mmsb2/mmu)*pow3(mmsb2))/pow2(-mmgl + mmsb2) + (8*DeltaInv(mmt,mmst1,mmgl)*
     pow2(log(mmgl/mmu))*pow3(mmsb2))/(3.*(-mmgl + mmsb2)) + (8*DeltaInv(mmt,
     mmst2,mmgl)*pow2(log(mmgl/mmu))*pow3(mmsb2))/(3.*(-mmgl + mmsb2)) + (16*
     pow2(log(mmgl/mmu))*pow3(mmsb2))/(9.*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1))
     + (16*pow2(log(mmgl/mmu))*pow3(mmsb2))/(9.*(-mmgl + mmsb1)*pow2(mmsb1 -
     mmsb2)) - (16*pow2(log(mmgl/mmu))*pow3(mmsb2))/(9.*(-mmgl + mmsb2)*pow2(
     mmsb1 - mmsb2)) + (4*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmgl/mmu))*
     pow3(mmsb2))/(3.*pow2(-mmgl + mmsb2)) + (2*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow2(log(mmgl/mmu))*pow3(mmsb2))/(3.*pow2(-mmgl + mmsb2)) + (4*mmst2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmgl/mmu))*pow3(mmsb2))/(3.*pow2(-mmgl +
     mmsb2)) + (2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmgl/mmu))*pow3(mmsb2))
     /(3.*pow2(-mmgl + mmsb2)) + (DeltaInv(mmt,mmsb2,mmst1)*pow2(log(mmsb2/mmu)
     )*pow3(mmsb2))/(3.*(-mmgl + mmsb2)) + (DeltaInv(mmt,mmsb2,mmst2)*pow2(log(
     mmsb2/mmu))*pow3(mmsb2))/(3.*(-mmgl + mmsb2)) + (4*mmst1*DeltaInv(mmt,
     mmsb2,mmst1)*pow2(log(mmsb2/mmu))*pow3(mmsb2))/(3.*pow2(-mmgl + mmsb2)) +
     (2*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(log(mmsb2/mmu))*pow3(mmsb2))/(3.*
     pow2(-mmgl + mmsb2)) + (4*mmst2*DeltaInv(mmt,mmsb2,mmst2)*pow2(log(mmsb2/
     mmu))*pow3(mmsb2))/(3.*pow2(-mmgl + mmsb2)) + (2*mmt*DeltaInv(mmt,mmsb2,
     mmst2)*pow2(log(mmsb2/mmu))*pow3(mmsb2))/(3.*pow2(-mmgl + mmsb2)) + Xt*((
     112*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl))/(3.*mgl*(mmst1 - mmst2)) + (
     112*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl))/(3.*mgl*(mmst1 - mmst2)) - (
     112*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl))/(3.*mgl*(mmst1 - mmst2)) + (
     112*mmgl*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl))/(3.*mgl*(-mmgl + mmsb1)
     *(mmst1 - mmst2)) + (112*mmgl*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl))/(
     3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (16*mmgl*mmsb1*mmt*zt2*DeltaInv(
     mmt,mmst1,mmgl))/(3.*mgl*(mmst1 - mmst2)) + (16*mmgl*mmsb2*mmt*zt2*
     DeltaInv(mmt,mmst1,mmgl))/(3.*mgl*(mmst1 - mmst2)) - (16*mmgl*mmst1*mmt*
     zt2*DeltaInv(mmt,mmst1,mmgl))/(3.*mgl*(mmst1 - mmst2)) + (16*mmgl*mmsb1*
     mmst1*mmt*zt2*DeltaInv(mmt,mmst1,mmgl))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 -
     mmst2)) + (16*mmgl*mmsb2*mmst1*mmt*zt2*DeltaInv(mmt,mmst1,mmgl))/(3.*mgl*(
     -mmgl + mmsb2)*(mmst1 - mmst2)) - (112*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst2,
     mmgl))/(3.*mgl*(mmst1 - mmst2)) - (112*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,
     mmgl))/(3.*mgl*(mmst1 - mmst2)) + (112*mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,
     mmgl))/(3.*mgl*(mmst1 - mmst2)) - (112*mmgl*mmsb1*mmst2*mmt*DeltaInv(mmt,
     mmst2,mmgl))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (112*mmgl*mmsb2*
     mmst2*mmt*DeltaInv(mmt,mmst2,mmgl))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 -
     mmst2)) - (16*mmgl*mmsb1*mmt*zt2*DeltaInv(mmt,mmst2,mmgl))/(3.*mgl*(mmst1
     - mmst2)) - (16*mmgl*mmsb2*mmt*zt2*DeltaInv(mmt,mmst2,mmgl))/(3.*mgl*(
     mmst1 - mmst2)) + (16*mmgl*mmst2*mmt*zt2*DeltaInv(mmt,mmst2,mmgl))/(3.*
     mgl*(mmst1 - mmst2)) - (16*mmgl*mmsb1*mmst2*mmt*zt2*DeltaInv(mmt,mmst2,
     mmgl))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (16*mmgl*mmsb2*mmst2*
     mmt*zt2*DeltaInv(mmt,mmst2,mmgl))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2))
     + (16*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu))/(3.*mgl*
     (mmst1 - mmst2)) - (16*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,
     mmst1,mmgl,mmu))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (16*mmgl*
     mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu))/(3.*mgl*(-
     mmgl + mmsb2)*(mmst1 - mmst2)) + (8*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,
     mmgl)*Fin3(mmt,mmst1,mmgl,mmu))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) +
     (8*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu))/(3.*
     mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (16*mmgl*mmt*DeltaInv(mmt,mmst2,
     mmgl)*Fin3(mmt,mmst2,mmgl,mmu))/(3.*mgl*(mmst1 - mmst2)) + (16*mmgl*mmsb1*
     mmt*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu))/(3.*mgl*(-mmgl +
     mmsb1)*(mmst1 - mmst2)) + (16*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*
     Fin3(mmt,mmst2,mmgl,mmu))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (8*
     mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu))/(3.*mgl*
     (-mmgl + mmsb1)*(mmst1 - mmst2)) - (8*mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,
     mmgl)*Fin3(mmt,mmst2,mmgl,mmu))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) -
     (16*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu))/(mgl*(mmst1 -
     mmst2)) - (16*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu))/(mgl*
     (mmst1 - mmst2)) + (16*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/
     mmu))/(mgl*(mmst1 - mmst2)) - (16*mmgl*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,
     mmgl)*log(mmgl/mmu))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (16*mmgl*
     mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu))/(mgl*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) + (16*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*log(
     mmgl/mmu))/(mgl*(mmst1 - mmst2)) + (16*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,
     mmgl)*log(mmgl/mmu))/(mgl*(mmst1 - mmst2)) - (16*mmgl*mmst2*mmt*DeltaInv(
     mmt,mmst2,mmgl)*log(mmgl/mmu))/(mgl*(mmst1 - mmst2)) + (16*mmgl*mmsb1*
     mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu))/(mgl*(-mmgl + mmsb1)*(
     mmst1 - mmst2)) + (16*mmgl*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(
     mmgl/mmu))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (16*mmgl*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*log(mmst1/mmu))/(mgl*(mmst1 - mmst2)) + (16*mmgl*
     mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmst1/mmu))/(mgl*(-mmgl +
     mmsb1)*(mmst1 - mmst2)) + (16*mmgl*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,
     mmgl)*log(mmst1/mmu))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (16*mmgl*
     mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmst1/mmu))/(3.*mgl*(
     mmst1 - mmst2)) - (16*mmgl*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(
     mmgl/mmu)*log(mmst1/mmu))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (16*
     mmgl*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmst1/mmu)
     )/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (16*mmgl*mmst2*mmt*DeltaInv(
     mmt,mmst2,mmgl)*log(mmst2/mmu))/(mgl*(mmst1 - mmst2)) - (16*mmgl*mmsb1*
     mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmst2/mmu))/(mgl*(-mmgl + mmsb1)*(
     mmst1 - mmst2)) - (16*mmgl*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(
     mmst2/mmu))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (16*mmgl*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmst2/mmu))/(3.*mgl*(mmst1 -
     mmst2)) + (16*mmgl*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*
     log(mmst2/mmu))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (16*mmgl*mmsb2*
     mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmst2/mmu))/(3.*mgl*(
     -mmgl + mmsb2)*(mmst1 - mmst2)) - (16*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,
     mmgl)*log(mmt/mmu))/(mgl*(mmst1 - mmst2)) - (16*mmgl*mmsb2*mmt*DeltaInv(
     mmt,mmst1,mmgl)*log(mmt/mmu))/(mgl*(mmst1 - mmst2)) + (32*mmgl*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu))/(mgl*(mmst1 - mmst2)) - (32*mmgl*
     mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu))/(mgl*(-mmgl +
     mmsb1)*(mmst1 - mmst2)) - (32*mmgl*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,
     mmgl)*log(mmt/mmu))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (16*mmgl*
     mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu))/(mgl*(mmst1 - mmst2)) + (
     16*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu))/(mgl*(mmst1 -
     mmst2)) - (32*mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu))/(mgl*(
     mmst1 - mmst2)) + (32*mmgl*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(
     mmt/mmu))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (32*mmgl*mmsb2*mmst2*
     mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu))/(mgl*(-mmgl + mmsb2)*(mmst1 -
     mmst2)) + (16*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(
     mmt/mmu))/(3.*mgl*(mmst1 - mmst2)) + (16*mmgl*mmsb2*mmt*DeltaInv(mmt,
     mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu))/(3.*mgl*(mmst1 - mmst2)) - (32*
     mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu))/(3.*
     mgl*(mmst1 - mmst2)) + (32*mmgl*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*
     log(mmgl/mmu)*log(mmt/mmu))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (
     32*mmgl*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/
     mmu))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (16*mmgl*mmsb1*mmt*
     DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu))/(3.*mgl*(mmst1 -
     mmst2)) - (16*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(
     mmt/mmu))/(3.*mgl*(mmst1 - mmst2)) + (32*mmgl*mmst2*mmt*DeltaInv(mmt,
     mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu))/(3.*mgl*(mmst1 - mmst2)) - (32*
     mmgl*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu))/
     (3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (32*mmgl*mmsb2*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu))/(3.*mgl*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) + (112*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmgl))/(
     3.*mgl*(mmst1 - mmst2)) + (16*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmgl))
     /(3.*mgl*(mmst1 - mmst2)) - (112*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmgl))/
     (3.*mgl*(mmst1 - mmst2)) - (16*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmgl)
     )/(3.*mgl*(mmst1 - mmst2)) - (16*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/
     mmu)*pow2(mmgl))/(mgl*(mmst1 - mmst2)) + (16*mmt*DeltaInv(mmt,mmst2,mmgl)*
     log(mmgl/mmu)*pow2(mmgl))/(mgl*(mmst1 - mmst2)) - (16*mmt*DeltaInv(mmt,
     mmst1,mmgl)*log(mmt/mmu)*pow2(mmgl))/(mgl*(mmst1 - mmst2)) + (16*mmt*
     DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*pow2(mmgl))/(mgl*(mmst1 - mmst2)) +
     (16*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmgl))/(
     3.*mgl*(mmst1 - mmst2)) - (16*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*
     log(mmt/mmu)*pow2(mmgl))/(3.*mgl*(mmst1 - mmst2)) - (56*mmgl*mmt*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmsb1))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (8*
     mmgl*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1))/(mgl*(-mmgl + mmsb1)*(
     mmst1 - mmst2)) + (56*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/(mgl*
     (-mmgl + mmsb1)*(mmst1 - mmst2)) + (8*mmgl*mmt*zt2*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmsb1))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (24*mmgl*mmt*
     DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*pow2(mmsb1))/(mgl*(-mmgl + mmsb1)*(
     mmst1 - mmst2)) - (24*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*
     pow2(mmsb1))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (24*mmgl*mmt*
     DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*pow2(mmsb1))/(mgl*(-mmgl + mmsb1)*(
     mmst1 - mmst2)) - (24*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*pow2(
     mmsb1))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (8*mmgl*mmt*DeltaInv(mmt,
     mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmsb1))/(mgl*(-mmgl + mmsb1)*(
     mmst1 - mmst2)) + (8*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(
     mmt/mmu)*pow2(mmsb1))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (8*mmgl*mmt*
     Fin3(mmt,mmsb1,mmst1,mmu))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) -
     (8*mmgl*mmsb1*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1)*Fin3(mmt,mmsb1,mmst1,
     mmu))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmt*Fin3(mmt,
     mmsb1,mmst2,mmu))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*
     mmsb1*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2)*Fin3(mmt,mmsb1,mmst2,mmu))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmt*Fin3(mmt,mmst1,
     mmgl,mmu))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmsb1*
     mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu))/(3.*mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmt*Fin3(mmt,mmst2,mmgl,mmu)
     )/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmsb1*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu))/(3.*mgl*(mmst1 - mmst2)
     *pow2(-mmgl + mmsb1)) - (56*mmgl*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(
     mmsb1))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmst1*mmt*
     zt2*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1))/(3.*mgl*(mmst1 - mmst2)*pow2(-
     mmgl + mmsb1)) + (56*mmgl*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))
     /(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmst2*mmt*zt2*
     DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl
     + mmsb1)) - (56*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmst1*mmt*zt2*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmsb1))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1))
     + (56*mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/(3.*mgl*(mmst1
     - mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmst2*mmt*zt2*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmsb1))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*
     mmt*DeltaInv(mmt,mmsb1,mmst1)*Fin3(mmt,mmsb1,mmst1,mmu)*pow2(mmsb1))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmt*DeltaInv(mmt,mmsb1,
     mmst2)*Fin3(mmt,mmsb1,mmst2,mmu)*pow2(mmsb1))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb1)) + (8*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,
     mmst1,mmgl,mmu)*pow2(mmsb1))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1))
     - (8*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*pow2(
     mmsb1))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*pow2(mmsb1))/(mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb1)) - (8*mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(
     mmgl/mmu)*pow2(mmsb1))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*
     mmgl*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1)*log(mmsb1/mmu)*pow2(mmsb1))/(mgl*
     (mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmst2*mmt*DeltaInv(mmt,
     mmsb1,mmst2)*log(mmsb1/mmu)*pow2(mmsb1))/(mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) - (8*mmgl*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1)*log(mmst1/mmu)*pow2(
     mmsb1))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*log(mmst1/mmu)*pow2(mmsb1))/(mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb1)) + (8*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(
     mmgl/mmu)*log(mmst1/mmu)*pow2(mmsb1))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) + (8*mmgl*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1)*log(mmsb1/mmu)*log(
     mmst1/mmu)*pow2(mmsb1))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*
     mmgl*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2)*log(mmst2/mmu)*pow2(mmsb1))/(mgl*
     (mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmst2*mmt*DeltaInv(mmt,
     mmst2,mmgl)*log(mmst2/mmu)*pow2(mmsb1))/(mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) - (8*mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(
     mmst2/mmu)*pow2(mmsb1))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*
     mmgl*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2)*log(mmsb1/mmu)*log(mmst2/mmu)*
     pow2(mmsb1))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (16*mmgl*
     mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1)*log(mmt/mmu)*pow2(mmsb1))/(mgl*(mmst1
     - mmst2)*pow2(-mmgl + mmsb1)) - (16*mmgl*mmst2*mmt*DeltaInv(mmt,mmsb1,
     mmst2)*log(mmt/mmu)*pow2(mmsb1))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1))
     + (16*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*pow2(mmsb1))/(
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (16*mmgl*mmst2*mmt*DeltaInv(
     mmt,mmst2,mmgl)*log(mmt/mmu)*pow2(mmsb1))/(mgl*(mmst1 - mmst2)*pow2(-mmgl
     + mmsb1)) - (16*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(
     mmt/mmu)*pow2(mmsb1))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (16*
     mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(
     mmsb1))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (16*mmgl*mmst1*mmt*
     DeltaInv(mmt,mmsb1,mmst1)*log(mmsb1/mmu)*log(mmt/mmu)*pow2(mmsb1))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (16*mmgl*mmst2*mmt*DeltaInv(
     mmt,mmsb1,mmst2)*log(mmsb1/mmu)*log(mmt/mmu)*pow2(mmsb1))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb1)) - (56*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmsb2))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (8*mmgl*mmt*zt2*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmsb2))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (56*
     mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/(mgl*(-mmgl + mmsb2)*(mmst1
     - mmst2)) + (8*mmgl*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/(mgl*(-
     mmgl + mmsb2)*(mmst1 - mmst2)) + (24*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*
     log(mmgl/mmu)*pow2(mmsb2))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (24*
     mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*pow2(mmsb2))/(mgl*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) + (24*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmt/
     mmu)*pow2(mmsb2))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (24*mmgl*mmt*
     DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*pow2(mmsb2))/(mgl*(-mmgl + mmsb2)*(
     mmst1 - mmst2)) - (8*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(
     mmt/mmu)*pow2(mmsb2))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (8*mmgl*mmt*
     DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmsb2))/(mgl*(-
     mmgl + mmsb2)*(mmst1 - mmst2)) + (8*mmgl*mmt*Fin3(mmt,mmsb2,mmst1,mmu))/(
     3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmsb2*mmst1*mmt*
     DeltaInv(mmt,mmsb2,mmst1)*Fin3(mmt,mmsb2,mmst1,mmu))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmt*Fin3(mmt,mmsb2,mmst2,mmu))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmsb2*mmst2*mmt*
     DeltaInv(mmt,mmsb2,mmst2)*Fin3(mmt,mmsb2,mmst2,mmu))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmt*Fin3(mmt,mmst1,mmgl,mmu))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmsb2*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu))/(3.*mgl*(mmst1 - mmst2)
     *pow2(-mmgl + mmsb2)) + (8*mmgl*mmt*Fin3(mmt,mmst2,mmgl,mmu))/(3.*mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmsb2*mmst2*mmt*DeltaInv(
     mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu))/(3.*mgl*(mmst1 - mmst2)*pow2(-
     mmgl + mmsb2)) - (56*mmgl*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2))
     /(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmst1*mmt*zt2*
     DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl
     + mmsb2)) + (56*mmgl*mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmst2*mmt*zt2*DeltaInv(
     mmt,mmsb2,mmst2)*pow2(mmsb2))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2))
     - (56*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2))/(3.*mgl*(mmst1
     - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmst1*mmt*zt2*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmsb2))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (56*
     mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmst2*mmt*zt2*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmsb2))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*
     mmt*DeltaInv(mmt,mmsb2,mmst1)*Fin3(mmt,mmsb2,mmst1,mmu)*pow2(mmsb2))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmt*DeltaInv(mmt,mmsb2,
     mmst2)*Fin3(mmt,mmsb2,mmst2,mmu)*pow2(mmsb2))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) + (8*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,
     mmst1,mmgl,mmu)*pow2(mmsb2))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2))
     - (8*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*pow2(
     mmsb2))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*pow2(mmsb2))/(mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) - (8*mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(
     mmgl/mmu)*pow2(mmsb2))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*
     mmgl*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1)*log(mmsb2/mmu)*pow2(mmsb2))/(mgl*
     (mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmst2*mmt*DeltaInv(mmt,
     mmsb2,mmst2)*log(mmsb2/mmu)*pow2(mmsb2))/(mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb2)) - (8*mmgl*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1)*log(mmst1/mmu)*pow2(
     mmsb2))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*log(mmst1/mmu)*pow2(mmsb2))/(mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) + (8*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(
     mmgl/mmu)*log(mmst1/mmu)*pow2(mmsb2))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb2)) + (8*mmgl*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1)*log(mmsb2/mmu)*log(
     mmst1/mmu)*pow2(mmsb2))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*
     mmgl*mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2)*log(mmst2/mmu)*pow2(mmsb2))/(mgl*
     (mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmst2*mmt*DeltaInv(mmt,
     mmst2,mmgl)*log(mmst2/mmu)*pow2(mmsb2))/(mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb2)) - (8*mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(
     mmst2/mmu)*pow2(mmsb2))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*
     mmgl*mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2)*log(mmsb2/mmu)*log(mmst2/mmu)*
     pow2(mmsb2))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (16*mmgl*
     mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1)*log(mmt/mmu)*pow2(mmsb2))/(mgl*(mmst1
     - mmst2)*pow2(-mmgl + mmsb2)) - (16*mmgl*mmst2*mmt*DeltaInv(mmt,mmsb2,
     mmst2)*log(mmt/mmu)*pow2(mmsb2))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2))
     + (16*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*pow2(mmsb2))/(
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (16*mmgl*mmst2*mmt*DeltaInv(
     mmt,mmst2,mmgl)*log(mmt/mmu)*pow2(mmsb2))/(mgl*(mmst1 - mmst2)*pow2(-mmgl
     + mmsb2)) - (16*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(
     mmt/mmu)*pow2(mmsb2))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (16*
     mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(
     mmsb2))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (16*mmgl*mmst1*mmt*
     DeltaInv(mmt,mmsb2,mmst1)*log(mmsb2/mmu)*log(mmt/mmu)*pow2(mmsb2))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (16*mmgl*mmst2*mmt*DeltaInv(
     mmt,mmsb2,mmst2)*log(mmsb2/mmu)*log(mmt/mmu)*pow2(mmsb2))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*log(
     mmst1/mmu)*pow2(mmst1))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (8*mmgl*
     mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmst1/mmu)*pow2(mmst1))/(mgl*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) + (8*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/
     mmu)*log(mmst1/mmu)*pow2(mmst1))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2))
     + (8*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmst1/mmu)*pow2(
     mmst1))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (8*mmgl*mmt*DeltaInv(
     mmt,mmst1,mmgl)*log(mmt/mmu)*pow2(mmst1))/(mgl*(-mmgl + mmsb1)*(mmst1 -
     mmst2)) + (8*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*pow2(mmst1))/(
     mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (8*mmgl*mmt*DeltaInv(mmt,mmst1,
     mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmst1))/(3.*mgl*(-mmgl + mmsb1)*(
     mmst1 - mmst2)) - (8*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(
     mmt/mmu)*pow2(mmst1))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (8*mmgl*
     mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst1)*log(mmst1/mmu)*pow2(mmst1))/(mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmsb1*mmt*DeltaInv(mmt,
     mmst1,mmgl)*log(mmst1/mmu)*pow2(mmst1))/(mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) - (8*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(
     mmst1/mmu)*pow2(mmst1))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*
     mmgl*mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst1)*log(mmsb1/mmu)*log(mmst1/mmu)*
     pow2(mmst1))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmsb1*
     mmt*DeltaInv(mmt,mmsb1,mmst1)*log(mmt/mmu)*pow2(mmst1))/(mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*
     log(mmt/mmu)*pow2(mmst1))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*
     mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(
     mmst1))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmsb1*mmt*
     DeltaInv(mmt,mmsb1,mmst1)*log(mmsb1/mmu)*log(mmt/mmu)*pow2(mmst1))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmsb2*mmt*DeltaInv(mmt,
     mmsb2,mmst1)*log(mmst1/mmu)*pow2(mmst1))/(mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb2)) + (8*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmst1/mmu)*pow2(
     mmst1))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmsb2*mmt*
     DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmst1/mmu)*pow2(mmst1))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmsb2*mmt*DeltaInv(mmt,
     mmsb2,mmst1)*log(mmsb2/mmu)*log(mmst1/mmu)*pow2(mmst1))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst1)*
     log(mmt/mmu)*pow2(mmst1))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*
     mmgl*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*pow2(mmst1))/(mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmsb2*mmt*DeltaInv(mmt,
     mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmst1))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst1)*
     log(mmsb2/mmu)*log(mmt/mmu)*pow2(mmst1))/(3.*mgl*(mmst1 - mmst2)*pow2(-
     mmgl + mmsb2)) + (8*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmst2/mmu)*pow2(
     mmst2))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (8*mmgl*mmt*DeltaInv(mmt,
     mmst2,mmgl)*log(mmst2/mmu)*pow2(mmst2))/(mgl*(-mmgl + mmsb2)*(mmst1 -
     mmst2)) - (8*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmst2/
     mmu)*pow2(mmst2))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (8*mmgl*mmt*
     DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmst2/mmu)*pow2(mmst2))/(3.*
     mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (8*mmgl*mmt*DeltaInv(mmt,mmst2,
     mmgl)*log(mmt/mmu)*pow2(mmst2))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (
     8*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*pow2(mmst2))/(mgl*(-mmgl
     + mmsb2)*(mmst1 - mmst2)) + (8*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/
     mmu)*log(mmt/mmu)*pow2(mmst2))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) +
     (8*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(
     mmst2))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (8*mmgl*mmsb1*mmt*
     DeltaInv(mmt,mmsb1,mmst2)*log(mmst2/mmu)*pow2(mmst2))/(mgl*(mmst1 - mmst2)
     *pow2(-mmgl + mmsb1)) - (8*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*log(
     mmst2/mmu)*pow2(mmst2))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*
     mmgl*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmst2/mmu)*pow2(
     mmst2))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmsb1*mmt*
     DeltaInv(mmt,mmsb1,mmst2)*log(mmsb1/mmu)*log(mmst2/mmu)*pow2(mmst2))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmsb1*mmt*DeltaInv(mmt,
     mmsb1,mmst2)*log(mmt/mmu)*pow2(mmst2))/(mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) + (8*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*pow2(
     mmst2))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmsb1*mmt*
     DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmst2))/(3.*mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmsb1*mmt*DeltaInv(mmt,
     mmsb1,mmst2)*log(mmsb1/mmu)*log(mmt/mmu)*pow2(mmst2))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst2)*
     log(mmst2/mmu)*pow2(mmst2))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (
     8*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmst2/mmu)*pow2(mmst2))/(
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmsb2*mmt*DeltaInv(mmt,
     mmst2,mmgl)*log(mmgl/mmu)*log(mmst2/mmu)*pow2(mmst2))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst2)*
     log(mmsb2/mmu)*log(mmst2/mmu)*pow2(mmst2))/(3.*mgl*(mmst1 - mmst2)*pow2(-
     mmgl + mmsb2)) + (8*mmgl*mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst2)*log(mmt/mmu)*
     pow2(mmst2))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmsb2*
     mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*pow2(mmst2))/(mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*
     log(mmgl/mmu)*log(mmt/mmu)*pow2(mmst2))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl
     + mmsb2)) - (8*mmgl*mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst2)*log(mmsb2/mmu)*
     log(mmt/mmu)*pow2(mmst2))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (
     112*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)) + (
     112*mmgl*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt))/(3.*mgl*(-mmgl + mmsb1)
     *(mmst1 - mmst2)) + (112*mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt))/(
     3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (112*mmgl*mmst1*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmt))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (112*
     mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt))/(3.*mgl*(-mmgl + mmsb2)*(
     mmst1 - mmst2)) - (16*mmgl*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt))/(3.*
     mgl*(mmst1 - mmst2)) + (16*mmgl*mmsb1*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmt))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (16*mmgl*mmsb2*zt2*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmt))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 -
     mmst2)) + (16*mmgl*mmst1*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt))/(3.*mgl*(
     -mmgl + mmsb1)*(mmst1 - mmst2)) + (16*mmgl*mmst1*zt2*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmt))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (112*mmgl*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)) - (112*mmgl*
     mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 -
     mmst2)) - (112*mmgl*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt))/(3.*mgl*(-
     mmgl + mmsb2)*(mmst1 - mmst2)) - (112*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmt))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (112*mmgl*mmst2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmt))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 -
     mmst2)) + (16*mmgl*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt))/(3.*mgl*(mmst1
     - mmst2)) - (16*mmgl*mmsb1*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt))/(3.*
     mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (16*mmgl*mmsb2*zt2*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmt))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (16*
     mmgl*mmst2*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt))/(3.*mgl*(-mmgl + mmsb1)
     *(mmst1 - mmst2)) - (16*mmgl*mmst2*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt))
     /(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (8*mmgl*DeltaInv(mmt,mmst1,
     mmgl)*Fin3(mmt,mmst1,mmgl,mmu)*pow2(mmt))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 -
     mmst2)) + (8*mmgl*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu)*pow2(
     mmt))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (8*mmgl*DeltaInv(mmt,
     mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*pow2(mmt))/(3.*mgl*(-mmgl + mmsb1)*(
     mmst1 - mmst2)) - (8*mmgl*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,
     mmu)*pow2(mmt))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (16*mmgl*
     DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*pow2(mmt))/(mgl*(mmst1 - mmst2)) -
     (16*mmgl*mmsb1*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*pow2(mmt))/(mgl*(-
     mmgl + mmsb1)*(mmst1 - mmst2)) - (16*mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl)*
     log(mmgl/mmu)*pow2(mmt))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (16*mmgl*
     DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*pow2(mmt))/(mgl*(mmst1 - mmst2)) +
     (16*mmgl*mmsb1*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*pow2(mmt))/(mgl*(-
     mmgl + mmsb1)*(mmst1 - mmst2)) + (16*mmgl*mmsb2*DeltaInv(mmt,mmst2,mmgl)*
     log(mmgl/mmu)*pow2(mmt))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (8*mmgl*
     mmst1*DeltaInv(mmt,mmst1,mmgl)*log(mmst1/mmu)*pow2(mmt))/(mgl*(-mmgl +
     mmsb1)*(mmst1 - mmst2)) - (8*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*log(
     mmst1/mmu)*pow2(mmt))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (8*mmgl*
     mmst1*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmst1/mmu)*pow2(mmt))/(3.
     *mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (8*mmgl*mmst1*DeltaInv(mmt,mmst1,
     mmgl)*log(mmgl/mmu)*log(mmst1/mmu)*pow2(mmt))/(3.*mgl*(-mmgl + mmsb2)*(
     mmst1 - mmst2)) + (8*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(mmst2/mmu)*
     pow2(mmt))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (8*mmgl*mmst2*DeltaInv(
     mmt,mmst2,mmgl)*log(mmst2/mmu)*pow2(mmt))/(mgl*(-mmgl + mmsb2)*(mmst1 -
     mmst2)) + (8*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmst2/
     mmu)*pow2(mmt))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (8*mmgl*mmst2*
     DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmst2/mmu)*pow2(mmt))/(3.*mgl*(
     -mmgl + mmsb2)*(mmst1 - mmst2)) + (16*mmgl*DeltaInv(mmt,mmst1,mmgl)*log(
     mmt/mmu)*pow2(mmt))/(mgl*(mmst1 - mmst2)) - (16*mmgl*mmsb1*DeltaInv(mmt,
     mmst1,mmgl)*log(mmt/mmu)*pow2(mmt))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2))
     - (16*mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*pow2(mmt))/(mgl*(-
     mmgl + mmsb2)*(mmst1 - mmst2)) - (24*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*
     log(mmt/mmu)*pow2(mmt))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (24*mmgl*
     mmst1*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*pow2(mmt))/(mgl*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) - (16*mmgl*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*
     pow2(mmt))/(mgl*(mmst1 - mmst2)) + (16*mmgl*mmsb1*DeltaInv(mmt,mmst2,mmgl)
     *log(mmt/mmu)*pow2(mmt))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (16*mmgl*
     mmsb2*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*pow2(mmt))/(mgl*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) + (24*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(mmt/
     mmu)*pow2(mmt))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (24*mmgl*mmst2*
     DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*pow2(mmt))/(mgl*(-mmgl + mmsb2)*(
     mmst1 - mmst2)) - (16*mmgl*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/
     mmu)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)) + (16*mmgl*mmsb1*DeltaInv(mmt,
     mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmt))/(3.*mgl*(-mmgl + mmsb1)*
     (mmst1 - mmst2)) + (16*mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*
     log(mmt/mmu)*pow2(mmt))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (8*
     mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmt))/
     (3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (8*mmgl*mmst1*DeltaInv(mmt,
     mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmt))/(3.*mgl*(-mmgl + mmsb2)*
     (mmst1 - mmst2)) + (16*mmgl*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(
     mmt/mmu)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)) - (16*mmgl*mmsb1*DeltaInv(
     mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmt))/(3.*mgl*(-mmgl +
     mmsb1)*(mmst1 - mmst2)) - (16*mmgl*mmsb2*DeltaInv(mmt,mmst2,mmgl)*log(
     mmgl/mmu)*log(mmt/mmu)*pow2(mmt))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2))
     - (8*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(
     mmt))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (8*mmgl*mmst2*DeltaInv(
     mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmt))/(3.*mgl*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) + (16*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*log(
     mmst1/mmu)*log(mmt/mmu)*pow2(mmt))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)
     ) + (16*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*log(mmst1/mmu)*log(mmt/mmu)*
     pow2(mmt))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (16*mmgl*mmst2*
     DeltaInv(mmt,mmst2,mmgl)*log(mmst2/mmu)*log(mmt/mmu)*pow2(mmt))/(3.*mgl*(-
     mmgl + mmsb1)*(mmst1 - mmst2)) - (16*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*
     log(mmst2/mmu)*log(mmt/mmu)*pow2(mmt))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 -
     mmst2)) - (112*mmgl*mmsb1*mmst1*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmt))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (16*mmgl*mmsb1*mmst1*zt2*
     DeltaInv(mmt,mmsb1,mmst1)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) + (112*mmgl*mmsb1*mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmt))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (16*mmgl*mmsb1*mmst2*zt2*
     DeltaInv(mmt,mmsb1,mmst2)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) - (112*mmgl*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (16*mmgl*mmsb1*mmst1*zt2*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) + (112*mmgl*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (16*mmgl*mmsb1*mmst2*zt2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) - (8*mmgl*mmsb1*DeltaInv(mmt,mmsb1,mmst1)*Fin3(mmt,mmsb1,mmst1,
     mmu)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*
     mmsb1*DeltaInv(mmt,mmsb1,mmst2)*Fin3(mmt,mmsb1,mmst2,mmu)*pow2(mmt))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmsb1*DeltaInv(mmt,
     mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb1)) + (8*mmgl*mmsb1*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,
     mmst2,mmgl,mmu)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) +
     (8*mmgl*mmsb1*mmst1*DeltaInv(mmt,mmsb1,mmst1)*log(mmst1/mmu)*pow2(mmt))/(
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmsb1*mmst1*DeltaInv(
     mmt,mmst1,mmgl)*log(mmst1/mmu)*pow2(mmt))/(mgl*(mmst1 - mmst2)*pow2(-mmgl
     + mmsb1)) + (8*mmgl*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*
     log(mmst1/mmu)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (
     8*mmgl*mmsb1*mmst1*DeltaInv(mmt,mmsb1,mmst1)*log(mmsb1/mmu)*log(mmst1/mmu)
     *pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmsb1*
     mmst2*DeltaInv(mmt,mmsb1,mmst2)*log(mmst2/mmu)*pow2(mmt))/(mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)
     *log(mmst2/mmu)*pow2(mmt))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*
     mmgl*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmst2/mmu)*
     pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmsb1*
     mmst2*DeltaInv(mmt,mmsb1,mmst2)*log(mmsb1/mmu)*log(mmst2/mmu)*pow2(mmt))/(
     3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (24*mmgl*mmsb1*mmst1*
     DeltaInv(mmt,mmsb1,mmst1)*log(mmt/mmu)*pow2(mmt))/(mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb1)) - (24*mmgl*mmsb1*mmst2*DeltaInv(mmt,mmsb1,mmst2)*log(
     mmt/mmu)*pow2(mmt))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (24*mmgl*
     mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*pow2(mmt))/(mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb1)) - (24*mmgl*mmsb1*mmst2*DeltaInv(mmt,mmst2,
     mmgl)*log(mmt/mmu)*pow2(mmt))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) -
     (8*mmgl*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*
     pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmsb1*
     mmst2*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmt))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmsb1*mmst1*DeltaInv(
     mmt,mmsb1,mmst1)*log(mmsb1/mmu)*log(mmt/mmu)*pow2(mmt))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmsb1*mmst2*DeltaInv(mmt,mmsb1,
     mmst2)*log(mmsb1/mmu)*log(mmt/mmu)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb1)) - (16*mmgl*mmsb1*mmst1*DeltaInv(mmt,mmsb1,mmst1)*log(
     mmst1/mmu)*log(mmt/mmu)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) - (16*mmgl*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*log(mmst1/mmu)*
     log(mmt/mmu)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (
     16*mmgl*mmsb1*mmst2*DeltaInv(mmt,mmsb1,mmst2)*log(mmst2/mmu)*log(mmt/mmu)*
     pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (16*mmgl*mmsb1*
     mmst2*DeltaInv(mmt,mmst2,mmgl)*log(mmst2/mmu)*log(mmt/mmu)*pow2(mmt))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (56*mmgl*DeltaInv(mmt,mmsb1,
     mmst1)*pow2(mmsb1)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1))
     - (8*mmgl*zt2*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow2(mmt))/(3.*mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (56*mmgl*DeltaInv(mmt,mmsb1,mmst2)*
     pow2(mmsb1)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*
     mmgl*zt2*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*pow2(mmt))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb1)) - (56*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmsb1)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*
     zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(mmt))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb1)) + (56*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmsb1)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*
     zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(mmt))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/
     mmu)*pow2(mmsb1)*pow2(mmt))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (
     8*mmgl*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*pow2(mmsb1)*pow2(mmt))/(mgl*
     (mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*DeltaInv(mmt,mmsb1,mmst1)*
     log(mmsb1/mmu)*pow2(mmsb1)*pow2(mmt))/(mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) - (8*mmgl*DeltaInv(mmt,mmsb1,mmst2)*log(mmsb1/mmu)*pow2(mmsb1)*
     pow2(mmt))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*DeltaInv(
     mmt,mmsb1,mmst1)*log(mmt/mmu)*pow2(mmsb1)*pow2(mmt))/(mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb1)) - (8*mmgl*DeltaInv(mmt,mmsb1,mmst2)*log(mmt/mmu)*
     pow2(mmsb1)*pow2(mmt))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*
     mmgl*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*pow2(mmsb1)*pow2(mmt))/(mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*DeltaInv(mmt,mmst2,mmgl)*
     log(mmt/mmu)*pow2(mmsb1)*pow2(mmt))/(mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) - (8*mmgl*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*
     pow2(mmsb1)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*
     mmgl*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmsb1)*pow2(
     mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*DeltaInv(mmt,
     mmsb1,mmst1)*log(mmsb1/mmu)*log(mmt/mmu)*pow2(mmsb1)*pow2(mmt))/(3.*mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*DeltaInv(mmt,mmsb1,mmst2)*
     log(mmsb1/mmu)*log(mmt/mmu)*pow2(mmsb1)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)
     *pow2(-mmgl + mmsb1)) - (112*mmgl*mmsb2*mmst1*DeltaInv(mmt,mmsb2,mmst1)*
     pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (16*mmgl*mmsb2*
     mmst1*zt2*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) + (112*mmgl*mmsb2*mmst2*DeltaInv(mmt,mmsb2,mmst2)*
     pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (16*mmgl*mmsb2*
     mmst2*zt2*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) - (112*mmgl*mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (16*mmgl*mmsb2*
     mmst1*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) + (112*mmgl*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (16*mmgl*mmsb2*
     mmst2*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) - (8*mmgl*mmsb2*DeltaInv(mmt,mmsb2,mmst1)*Fin3(mmt,
     mmsb2,mmst1,mmu)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) +
     (8*mmgl*mmsb2*DeltaInv(mmt,mmsb2,mmst2)*Fin3(mmt,mmsb2,mmst2,mmu)*pow2(
     mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmsb2*
     DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu)*pow2(mmt))/(3.*mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmsb2*DeltaInv(mmt,mmst2,
     mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-
     mmgl + mmsb2)) + (8*mmgl*mmsb2*mmst1*DeltaInv(mmt,mmsb2,mmst1)*log(mmst1/
     mmu)*pow2(mmt))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmsb2*
     mmst1*DeltaInv(mmt,mmst1,mmgl)*log(mmst1/mmu)*pow2(mmt))/(mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)
     *log(mmgl/mmu)*log(mmst1/mmu)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-
     mmgl + mmsb2)) + (8*mmgl*mmsb2*mmst1*DeltaInv(mmt,mmsb2,mmst1)*log(mmsb2/
     mmu)*log(mmst1/mmu)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)
     ) - (8*mmgl*mmsb2*mmst2*DeltaInv(mmt,mmsb2,mmst2)*log(mmst2/mmu)*pow2(mmt)
     )/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmsb2*mmst2*
     DeltaInv(mmt,mmst2,mmgl)*log(mmst2/mmu)*pow2(mmt))/(mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) - (8*mmgl*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(
     mmgl/mmu)*log(mmst2/mmu)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb2)) - (8*mmgl*mmsb2*mmst2*DeltaInv(mmt,mmsb2,mmst2)*log(mmsb2/mmu)*
     log(mmst2/mmu)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (
     24*mmgl*mmsb2*mmst1*DeltaInv(mmt,mmsb2,mmst1)*log(mmt/mmu)*pow2(mmt))/(
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (24*mmgl*mmsb2*mmst2*DeltaInv(
     mmt,mmsb2,mmst2)*log(mmt/mmu)*pow2(mmt))/(mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb2)) + (24*mmgl*mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*pow2(
     mmt))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (24*mmgl*mmsb2*mmst2*
     DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*pow2(mmt))/(mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) - (8*mmgl*mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)*log(
     mmgl/mmu)*log(mmt/mmu)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb2)) + (8*mmgl*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(
     mmt/mmu)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*
     mmgl*mmsb2*mmst1*DeltaInv(mmt,mmsb2,mmst1)*log(mmsb2/mmu)*log(mmt/mmu)*
     pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmsb2*
     mmst2*DeltaInv(mmt,mmsb2,mmst2)*log(mmsb2/mmu)*log(mmt/mmu)*pow2(mmt))/(3.
     *mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (16*mmgl*mmsb2*mmst1*DeltaInv(
     mmt,mmsb2,mmst1)*log(mmst1/mmu)*log(mmt/mmu)*pow2(mmt))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb2)) - (16*mmgl*mmsb2*mmst1*DeltaInv(mmt,mmst1,
     mmgl)*log(mmst1/mmu)*log(mmt/mmu)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(
     -mmgl + mmsb2)) + (16*mmgl*mmsb2*mmst2*DeltaInv(mmt,mmsb2,mmst2)*log(
     mmst2/mmu)*log(mmt/mmu)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb2)) + (16*mmgl*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(mmst2/mmu)*
     log(mmt/mmu)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (
     56*mmgl*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2)*pow2(mmt))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*zt2*DeltaInv(mmt,mmsb2,mmst1)*pow2(
     mmsb2)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (56*mmgl*
     DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) + (8*mmgl*zt2*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*
     pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (56*mmgl*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) - (8*mmgl*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*
     pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (56*mmgl*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) + (8*mmgl*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*
     pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*
     DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*pow2(mmsb2)*pow2(mmt))/(mgl*(mmst1
     - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/
     mmu)*pow2(mmsb2)*pow2(mmt))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (
     8*mmgl*DeltaInv(mmt,mmsb2,mmst1)*log(mmsb2/mmu)*pow2(mmsb2)*pow2(mmt))/(
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*DeltaInv(mmt,mmsb2,
     mmst2)*log(mmsb2/mmu)*pow2(mmsb2)*pow2(mmt))/(mgl*(mmst1 - mmst2)*pow2(-
     mmgl + mmsb2)) + (8*mmgl*DeltaInv(mmt,mmsb2,mmst1)*log(mmt/mmu)*pow2(
     mmsb2)*pow2(mmt))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*
     DeltaInv(mmt,mmsb2,mmst2)*log(mmt/mmu)*pow2(mmsb2)*pow2(mmt))/(mgl*(mmst1
     - mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*DeltaInv(mmt,mmst1,mmgl)*log(mmt/
     mmu)*pow2(mmsb2)*pow2(mmt))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (
     8*mmgl*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*pow2(mmsb2)*pow2(mmt))/(mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*DeltaInv(mmt,mmst1,mmgl)*
     log(mmgl/mmu)*log(mmt/mmu)*pow2(mmsb2)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) + (8*mmgl*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(
     mmt/mmu)*pow2(mmsb2)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb2)) - (8*mmgl*DeltaInv(mmt,mmsb2,mmst1)*log(mmsb2/mmu)*log(mmt/mmu)*
     pow2(mmsb2)*pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*
     mmgl*DeltaInv(mmt,mmsb2,mmst2)*log(mmsb2/mmu)*log(mmt/mmu)*pow2(mmsb2)*
     pow2(mmt))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmsb1*
     mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmgl/mmu)))/(3.*mgl*(mmst1 - mmst2))
     + (8*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmgl/mmu)))/(3.*mgl*
     (mmst1 - mmst2)) - (8*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(log(
     mmgl/mmu)))/(3.*mgl*(mmst1 - mmst2)) + (8*mmgl*mmsb1*mmst1*mmt*DeltaInv(
     mmt,mmst1,mmgl)*pow2(log(mmgl/mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 -
     mmst2)) + (8*mmgl*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmgl/
     mmu)))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (8*mmgl*mmsb1*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmgl/mmu)))/(3.*mgl*(mmst1 - mmst2)) - (
     8*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmgl/mmu)))/(3.*mgl*(
     mmst1 - mmst2)) + (8*mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(
     mmgl/mmu)))/(3.*mgl*(mmst1 - mmst2)) - (8*mmgl*mmsb1*mmst2*mmt*DeltaInv(
     mmt,mmst2,mmgl)*pow2(log(mmgl/mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 -
     mmst2)) - (8*mmgl*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmgl/
     mmu)))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (8*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmgl)*pow2(log(mmgl/mmu)))/(3.*mgl*(mmst1 - mmst2)) - (8*
     mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmgl)*pow2(log(mmgl/mmu)))/(3.*mgl*(
     mmst1 - mmst2)) - (4*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(
     log(mmgl/mmu)))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (4*mmgl*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(log(mmgl/mmu)))/(mgl*(-mmgl +
     mmsb1)*(mmst1 - mmst2)) - (4*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmsb1)*pow2(log(mmgl/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) +
     (4*mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(log(mmgl/mmu))
     )/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (4*mmgl*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmsb2)*pow2(log(mmgl/mmu)))/(mgl*(-mmgl + mmsb2)*(mmst1 -
     mmst2)) + (4*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(log(mmgl/
     mmu)))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (4*mmgl*mmst1*mmt*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(log(mmgl/mmu)))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) + (4*mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmsb2)*pow2(log(mmgl/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) -
     (8*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow2(log(mmgl/mmu)))/(3.*mgl*(
     mmst1 - mmst2)) + (8*mmgl*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow2(
     log(mmgl/mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (8*mmgl*mmsb2*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow2(log(mmgl/mmu)))/(3.*mgl*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) + (8*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow2(
     log(mmgl/mmu)))/(3.*mgl*(mmst1 - mmst2)) - (8*mmgl*mmsb1*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmt)*pow2(log(mmgl/mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmst1
     - mmst2)) - (8*mmgl*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow2(log(
     mmgl/mmu)))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (4*mmgl*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(mmt)*pow2(log(mmgl/mmu)))/(3.*mgl*(mmst1
     - mmst2)*pow2(-mmgl + mmsb1)) + (4*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmsb1)*pow2(mmt)*pow2(log(mmgl/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) - (4*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(mmt)*pow2(log(
     mmgl/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (4*mmgl*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(mmt)*pow2(log(mmgl/mmu)))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (4*mmgl*mmst1*mmt*DeltaInv(mmt,
     mmsb1,mmst1)*pow2(mmsb1)*pow2(log(mmsb1/mmu)))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb1)) + (4*mmgl*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(
     mmsb1)*pow2(log(mmsb1/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1))
     - (4*mmgl*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow2(mmt)*pow2(log(mmsb1/
     mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (4*mmgl*DeltaInv(
     mmt,mmsb1,mmst2)*pow2(mmsb1)*pow2(mmt)*pow2(log(mmsb1/mmu)))/(3.*mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (4*mmgl*mmst1*mmt*DeltaInv(mmt,
     mmsb2,mmst1)*pow2(mmsb2)*pow2(log(mmsb2/mmu)))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) + (4*mmgl*mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(
     mmsb2)*pow2(log(mmsb2/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2))
     - (4*mmgl*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2)*pow2(mmt)*pow2(log(mmsb2/
     mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (4*mmgl*DeltaInv(
     mmt,mmsb2,mmst2)*pow2(mmsb2)*pow2(mmt)*pow2(log(mmsb2/mmu)))/(3.*mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmst1*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow2(log(mmst1/mmu)))/(3.*mgl*(mmst1 - mmst2)) - (8*mmgl*
     mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmst1/mmu)))/(3.*mgl*(-
     mmgl + mmsb1)*(mmst1 - mmst2)) - (8*mmgl*mmsb2*mmst1*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow2(log(mmst1/mmu)))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2))
     + (4*mmgl*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow2(log(mmst1/
     mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (4*mmgl*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(log(mmst1/mmu)))/(3.*mgl*(mmst1
     - mmst2)*pow2(-mmgl + mmsb1)) + (4*mmgl*mmst1*mmt*DeltaInv(mmt,mmsb2,
     mmst1)*pow2(mmsb2)*pow2(log(mmst1/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-
     mmgl + mmsb2)) + (4*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*
     pow2(log(mmst1/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (4*
     mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow2(log(mmst1/mmu)))/(3.*
     mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (4*mmgl*mmt*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmst1)*pow2(log(mmst1/mmu)))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 -
     mmst2)) - (4*mmgl*mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1)*pow2(
     log(mmst1/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (4*mmgl*
     mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow2(log(mmst1/mmu)))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (4*mmgl*mmsb2*mmt*DeltaInv(mmt,
     mmsb2,mmst1)*pow2(mmst1)*pow2(log(mmst1/mmu)))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) - (4*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmst1)*pow2(log(mmst1/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2))
     + (4*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow2(log(mmst1/mmu)))/(
     3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (4*mmgl*mmst1*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmt)*pow2(log(mmst1/mmu)))/(3.*mgl*(-mmgl + mmsb2)*(mmst1
     - mmst2)) - (4*mmgl*mmsb1*mmst1*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmt)*pow2(
     log(mmst1/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (4*mmgl*
     mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow2(log(mmst1/mmu)))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (4*mmgl*mmsb2*mmst1*DeltaInv(
     mmt,mmsb2,mmst1)*pow2(mmt)*pow2(log(mmst1/mmu)))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) - (4*mmgl*mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmt)*pow2(log(mmst1/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) -
     (8*mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmst2/mmu)))/(3.*mgl*(
     mmst1 - mmst2)) + (8*mmgl*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(
     log(mmst2/mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (8*mmgl*mmsb2*
     mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmst2/mmu)))/(3.*mgl*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) - (4*mmgl*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2)*
     pow2(mmsb1)*pow2(log(mmst2/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) - (4*mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(log(
     mmst2/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (4*mmgl*mmst2*
     mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow2(log(mmst2/mmu)))/(3.*mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (4*mmgl*mmst2*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmsb2)*pow2(log(mmst2/mmu)))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) - (4*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*
     pow2(log(mmst2/mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (4*mmgl*
     mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow2(log(mmst2/mmu)))/(3.*mgl*(-
     mmgl + mmsb2)*(mmst1 - mmst2)) + (4*mmgl*mmsb1*mmt*DeltaInv(mmt,mmsb1,
     mmst2)*pow2(mmst2)*pow2(log(mmst2/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-
     mmgl + mmsb1)) + (4*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*
     pow2(log(mmst2/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (4*
     mmgl*mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2)*pow2(log(mmst2/mmu)))
     /(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (4*mmgl*mmsb2*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow2(log(mmst2/mmu)))/(3.*mgl*(mmst1
     - mmst2)*pow2(-mmgl + mmsb2)) - (4*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmt)*pow2(log(mmst2/mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) -
     (4*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow2(log(mmst2/mmu)))/(3.
     *mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (4*mmgl*mmsb1*mmst2*DeltaInv(mmt,
     mmsb1,mmst2)*pow2(mmt)*pow2(log(mmst2/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(
     -mmgl + mmsb1)) + (4*mmgl*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*
     pow2(log(mmst2/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (4*
     mmgl*mmsb2*mmst2*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmt)*pow2(log(mmst2/mmu)))
     /(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (4*mmgl*mmsb2*mmst2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow2(log(mmst2/mmu)))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow2(log(mmt/mmu)))/(3.*mgl*(mmst1 - mmst2)) + (8*mmgl*mmsb2*mmt*DeltaInv(
     mmt,mmst1,mmgl)*pow2(log(mmt/mmu)))/(3.*mgl*(mmst1 - mmst2)) - (16*mmgl*
     mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmt/mmu)))/(3.*mgl*(mmst1 -
     mmst2)) + (16*mmgl*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmt/
     mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (16*mmgl*mmsb2*mmst1*
     mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmt/mmu)))/(3.*mgl*(-mmgl + mmsb2)*(
     mmst1 - mmst2)) - (8*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmt/
     mmu)))/(3.*mgl*(mmst1 - mmst2)) - (8*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,
     mmgl)*pow2(log(mmt/mmu)))/(3.*mgl*(mmst1 - mmst2)) + (16*mmgl*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmt/mmu)))/(3.*mgl*(mmst1 - mmst2)) - (
     16*mmgl*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmt/mmu)))/(3.*
     mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (16*mmgl*mmsb2*mmst2*mmt*DeltaInv(
     mmt,mmst2,mmgl)*pow2(log(mmt/mmu)))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 -
     mmst2)) + (8*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmgl)*pow2(log(mmt/mmu)))/(
     3.*mgl*(mmst1 - mmst2)) - (8*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmgl)*pow2(
     log(mmt/mmu)))/(3.*mgl*(mmst1 - mmst2)) - (4*mmgl*mmt*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmsb1)*pow2(log(mmt/mmu)))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)
     ) + (4*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(log(mmt/mmu)))/(
     mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (8*mmgl*mmst1*mmt*DeltaInv(mmt,
     mmsb1,mmst1)*pow2(mmsb1)*pow2(log(mmt/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(
     -mmgl + mmsb1)) + (8*mmgl*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*
     pow2(log(mmt/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*
     mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(log(mmt/mmu)))/(
     3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmst2*mmt*DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(log(mmt/mmu)))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb1)) - (4*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*
     pow2(log(mmt/mmu)))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (4*mmgl*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(log(mmt/mmu)))/(mgl*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) - (8*mmgl*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1)*
     pow2(mmsb2)*pow2(log(mmt/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb2)) + (8*mmgl*mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow2(
     log(mmt/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*
     mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(log(mmt/mmu)))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmst2*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmsb2)*pow2(log(mmt/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-
     mmgl + mmsb2)) - (4*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow2(
     log(mmt/mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) - (4*mmgl*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow2(log(mmt/mmu)))/(3.*mgl*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) + (4*mmgl*mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst1)*
     pow2(mmst1)*pow2(log(mmt/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) + (4*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow2(log(
     mmt/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (4*mmgl*mmsb2*
     mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1)*pow2(log(mmt/mmu)))/(3.*mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (4*mmgl*mmsb2*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmst1)*pow2(log(mmt/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-
     mmgl + mmsb2)) + (4*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow2(
     log(mmt/mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (4*mmgl*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow2(log(mmt/mmu)))/(3.*mgl*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) - (4*mmgl*mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst2)*
     pow2(mmst2)*pow2(log(mmt/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) - (4*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow2(log(
     mmt/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (4*mmgl*mmsb2*
     mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2)*pow2(log(mmt/mmu)))/(3.*mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (4*mmgl*mmsb2*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmst2)*pow2(log(mmt/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-
     mmgl + mmsb2)) - (8*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow2(log(mmt/
     mmu)))/(3.*mgl*(mmst1 - mmst2)) + (8*mmgl*mmsb1*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmt)*pow2(log(mmt/mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmst1 - mmst2)) + (
     8*mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow2(log(mmt/mmu)))/(3.*
     mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (4*mmgl*mmst1*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmt)*pow2(log(mmt/mmu)))/(mgl*(-mmgl + mmsb1)*(mmst1 - mmst2))
     + (4*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow2(log(mmt/mmu)))/(
     mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (8*mmgl*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmt)*pow2(log(mmt/mmu)))/(3.*mgl*(mmst1 - mmst2)) - (8*mmgl*mmsb1*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow2(log(mmt/mmu)))/(3.*mgl*(-mmgl +
     mmsb1)*(mmst1 - mmst2)) - (8*mmgl*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)
     *pow2(log(mmt/mmu)))/(3.*mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (4*mmgl*
     mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow2(log(mmt/mmu)))/(mgl*(-mmgl +
     mmsb1)*(mmst1 - mmst2)) - (4*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)
     *pow2(log(mmt/mmu)))/(mgl*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (4*mmgl*
     mmsb1*mmst1*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmt)*pow2(log(mmt/mmu)))/(mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (4*mmgl*mmsb1*mmst2*DeltaInv(mmt,
     mmsb1,mmst2)*pow2(mmt)*pow2(log(mmt/mmu)))/(mgl*(mmst1 - mmst2)*pow2(-mmgl
     + mmsb1)) - (4*mmgl*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow2(
     log(mmt/mmu)))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (4*mmgl*mmsb1*
     mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow2(log(mmt/mmu)))/(mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb1)) - (4*mmgl*DeltaInv(mmt,mmsb1,mmst1)*pow2(
     mmsb1)*pow2(mmt)*pow2(log(mmt/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) + (4*mmgl*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*pow2(mmt)*pow2(
     log(mmt/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (4*mmgl*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(mmt)*pow2(log(mmt/mmu)))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (4*mmgl*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmsb1)*pow2(mmt)*pow2(log(mmt/mmu)))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb1)) - (4*mmgl*mmsb2*mmst1*DeltaInv(mmt,mmsb2,mmst1)*pow2(
     mmt)*pow2(log(mmt/mmu)))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (4*
     mmgl*mmsb2*mmst2*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmt)*pow2(log(mmt/mmu)))/(
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (4*mmgl*mmsb2*mmst1*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmt)*pow2(log(mmt/mmu)))/(mgl*(mmst1 - mmst2)*pow2(-
     mmgl + mmsb2)) + (4*mmgl*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*
     pow2(log(mmt/mmu)))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (4*mmgl*
     DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2)*pow2(mmt)*pow2(log(mmt/mmu)))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (4*mmgl*DeltaInv(mmt,mmsb2,
     mmst2)*pow2(mmsb2)*pow2(mmt)*pow2(log(mmt/mmu)))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) - (4*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(
     mmt)*pow2(log(mmt/mmu)))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (
     4*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(mmt)*pow2(log(mmt/mmu)))/
     (3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (56*mmgl*mmt*DeltaInv(mmt,
     mmsb1,mmst1)*pow3(mmsb1))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (
     8*mmgl*mmt*zt2*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb1)) - (56*mmgl*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow3(
     mmsb1))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmt*zt2*
     DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl
     + mmsb1)) + (56*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1))/(3.*mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmt*zt2*DeltaInv(mmt,mmst1,
     mmgl)*pow3(mmsb1))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (56*
     mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb1)) - (8*mmgl*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow3(
     mmsb1))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmt*
     DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*pow3(mmsb1))/(mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb1)) + (8*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*
     pow3(mmsb1))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmt*
     DeltaInv(mmt,mmsb1,mmst1)*log(mmsb1/mmu)*pow3(mmsb1))/(mgl*(mmst1 - mmst2)
     *pow2(-mmgl + mmsb1)) + (8*mmgl*mmt*DeltaInv(mmt,mmsb1,mmst2)*log(mmsb1/
     mmu)*pow3(mmsb1))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmt*
     DeltaInv(mmt,mmsb1,mmst1)*log(mmt/mmu)*pow3(mmsb1))/(mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb1)) + (8*mmgl*mmt*DeltaInv(mmt,mmsb1,mmst2)*log(mmt/mmu)*
     pow3(mmsb1))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmt*
     DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*pow3(mmsb1))/(mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb1)) + (8*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*
     pow3(mmsb1))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmt*
     DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow3(mmsb1))/(3.*mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)
     *log(mmgl/mmu)*log(mmt/mmu)*pow3(mmsb1))/(3.*mgl*(mmst1 - mmst2)*pow2(-
     mmgl + mmsb1)) + (8*mmgl*mmt*DeltaInv(mmt,mmsb1,mmst1)*log(mmsb1/mmu)*log(
     mmt/mmu)*pow3(mmsb1))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (8*
     mmgl*mmt*DeltaInv(mmt,mmsb1,mmst2)*log(mmsb1/mmu)*log(mmt/mmu)*pow3(mmsb1)
     )/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (4*mmgl*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow2(log(mmgl/mmu))*pow3(mmsb1))/(3.*mgl*(mmst1 - mmst2)*pow2(
     -mmgl + mmsb1)) - (4*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmgl/mmu))
     *pow3(mmsb1))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (4*mmgl*mmt*
     DeltaInv(mmt,mmsb1,mmst1)*pow2(log(mmsb1/mmu))*pow3(mmsb1))/(3.*mgl*(mmst1
     - mmst2)*pow2(-mmgl + mmsb1)) - (4*mmgl*mmt*DeltaInv(mmt,mmsb1,mmst2)*
     pow2(log(mmsb1/mmu))*pow3(mmsb1))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) + (4*mmgl*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(log(mmt/mmu))*pow3(
     mmsb1))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (4*mmgl*mmt*
     DeltaInv(mmt,mmsb1,mmst2)*pow2(log(mmt/mmu))*pow3(mmsb1))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb1)) + (4*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     log(mmt/mmu))*pow3(mmsb1))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) -
     (4*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmt/mmu))*pow3(mmsb1))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (16*mmgl*mmsb1*mmt*Fin3(mmt,
     mmsb1,mmst1,mmu))/(3.*mgl*(mmst1 - mmst2)*pow3(-mmgl + mmsb1)) + (16*mmgl*
     mmsb1*mmt*Fin3(mmt,mmsb1,mmst2,mmu))/(3.*mgl*(mmst1 - mmst2)*pow3(-mmgl +
     mmsb1)) + (16*mmgl*mmsb1*mmt*Fin3(mmt,mmst1,mmgl,mmu))/(3.*mgl*(mmst1 -
     mmst2)*pow3(-mmgl + mmsb1)) - (16*mmgl*mmsb1*mmt*Fin3(mmt,mmst2,mmgl,mmu))
     /(3.*mgl*(mmst1 - mmst2)*pow3(-mmgl + mmsb1)) + (56*mmgl*mmt*DeltaInv(mmt,
     mmsb2,mmst1)*pow3(mmsb2))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (
     8*mmgl*mmt*zt2*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmsb2))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb2)) - (56*mmgl*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow3(
     mmsb2))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmt*zt2*
     DeltaInv(mmt,mmsb2,mmst2)*pow3(mmsb2))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl
     + mmsb2)) + (56*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2))/(3.*mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmt*zt2*DeltaInv(mmt,mmst1,
     mmgl)*pow3(mmsb2))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (56*
     mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2))/(3.*mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) - (8*mmgl*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow3(
     mmsb2))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmt*
     DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*pow3(mmsb2))/(mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) + (8*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*
     pow3(mmsb2))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmt*
     DeltaInv(mmt,mmsb2,mmst1)*log(mmsb2/mmu)*pow3(mmsb2))/(mgl*(mmst1 - mmst2)
     *pow2(-mmgl + mmsb2)) + (8*mmgl*mmt*DeltaInv(mmt,mmsb2,mmst2)*log(mmsb2/
     mmu)*pow3(mmsb2))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmt*
     DeltaInv(mmt,mmsb2,mmst1)*log(mmt/mmu)*pow3(mmsb2))/(mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) + (8*mmgl*mmt*DeltaInv(mmt,mmsb2,mmst2)*log(mmt/mmu)*
     pow3(mmsb2))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmt*
     DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*pow3(mmsb2))/(mgl*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) + (8*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*
     pow3(mmsb2))/(mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmt*
     DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow3(mmsb2))/(3.*mgl*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)
     *log(mmgl/mmu)*log(mmt/mmu)*pow3(mmsb2))/(3.*mgl*(mmst1 - mmst2)*pow2(-
     mmgl + mmsb2)) + (8*mmgl*mmt*DeltaInv(mmt,mmsb2,mmst1)*log(mmsb2/mmu)*log(
     mmt/mmu)*pow3(mmsb2))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (8*
     mmgl*mmt*DeltaInv(mmt,mmsb2,mmst2)*log(mmsb2/mmu)*log(mmt/mmu)*pow3(mmsb2)
     )/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (4*mmgl*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow2(log(mmgl/mmu))*pow3(mmsb2))/(3.*mgl*(mmst1 - mmst2)*pow2(
     -mmgl + mmsb2)) - (4*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmgl/mmu))
     *pow3(mmsb2))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (4*mmgl*mmt*
     DeltaInv(mmt,mmsb2,mmst1)*pow2(log(mmsb2/mmu))*pow3(mmsb2))/(3.*mgl*(mmst1
     - mmst2)*pow2(-mmgl + mmsb2)) - (4*mmgl*mmt*DeltaInv(mmt,mmsb2,mmst2)*
     pow2(log(mmsb2/mmu))*pow3(mmsb2))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb2)) + (4*mmgl*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(log(mmt/mmu))*pow3(
     mmsb2))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (4*mmgl*mmt*
     DeltaInv(mmt,mmsb2,mmst2)*pow2(log(mmt/mmu))*pow3(mmsb2))/(3.*mgl*(mmst1 -
     mmst2)*pow2(-mmgl + mmsb2)) + (4*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     log(mmt/mmu))*pow3(mmsb2))/(3.*mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) -
     (4*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmt/mmu))*pow3(mmsb2))/(3.*
     mgl*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + Xb*((-224*mmgl*mmsb1*mmt*
     DeltaInv(mmt,mmst1,mmgl))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (224*
     mmgl*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl))/(3.*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) + (224*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl))/(3.*(mmsb1 -
     mmsb2)*(mmst1 - mmst2)) - (224*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl))/(
     3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmgl*mmsb1*mmt*zt2*DeltaInv(mmt,
     mmst1,mmgl))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*mmgl*mmsb2*mmt*
     zt2*DeltaInv(mmt,mmst1,mmgl))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*
     mmsb1*mmst1*mmt*zt2*DeltaInv(mmt,mmst1,mmgl))/(3.*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) - (32*mmsb2*mmst1*mmt*zt2*DeltaInv(mmt,mmst1,mmgl))/(3.*(mmsb1 -
     mmsb2)*(mmst1 - mmst2)) + (224*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl))/(
     3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (224*mmgl*mmsb2*mmt*DeltaInv(mmt,
     mmst2,mmgl))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (224*mmsb1*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (224*
     mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl))/(3.*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) + (32*mmgl*mmsb1*mmt*zt2*DeltaInv(mmt,mmst2,mmgl))/(3.*(mmsb1 -
     mmsb2)*(mmst1 - mmst2)) - (32*mmgl*mmsb2*mmt*zt2*DeltaInv(mmt,mmst2,mmgl))
     /(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmsb1*mmst2*mmt*zt2*DeltaInv(
     mmt,mmst2,mmgl))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*mmsb2*mmst2*
     mmt*zt2*DeltaInv(mmt,mmst2,mmgl))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (
     16*mmt*Fin3(mmt,mmsb1,mmst1,mmu))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) - (16*mmt*Fin3(mmt,mmsb1,mmst2,mmu))/(3.*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)*(mmst1 - mmst2)) - (16*mmt*Fin3(mmt,mmsb2,mmst1,mmu))/(3.*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (16*mmt*Fin3(mmt,mmsb2,
     mmst2,mmu))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (16*
     mmt*Fin3(mmt,mmst1,mmgl,mmu))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) + (16*mmt*Fin3(mmt,mmst1,mmgl,mmu))/(3.*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) - (32*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,
     mmst1,mmgl,mmu))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*mmsb2*mmt*
     DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu))/(3.*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) - (32*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,
     mmst1,mmgl,mmu))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (
     32*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu))/(3.*
     (mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (16*mmt*Fin3(mmt,mmst2,
     mmgl,mmu))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (16*mmt*
     Fin3(mmt,mmst2,mmgl,mmu))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 -
     mmst2)) + (32*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu))
     /(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmsb2*mmt*DeltaInv(mmt,mmst2,
     mmgl)*Fin3(mmt,mmst2,mmgl,mmu))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (
     32*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu))/(3.*
     (-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmsb2*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu))/(3.*(mmsb1 - mmsb2)*(-
     mmgl + mmsb2)*(mmst1 - mmst2)) + (32*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,
     mmgl)*log(mmgl/mmu))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmgl*mmsb2*
     mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu))/((mmsb1 - mmsb2)*(mmst1 -
     mmst2)) - (32*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu))/((
     mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,
     mmgl)*log(mmgl/mmu))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmgl*mmsb1*
     mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu))/((mmsb1 - mmsb2)*(mmst1 -
     mmst2)) + (32*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu))/((
     mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,
     mmgl)*log(mmgl/mmu))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmsb2*mmst2*
     mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu))/((mmsb1 - mmsb2)*(mmst1 -
     mmst2)) + (32*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmst1/mmu))/((
     mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,
     mmgl)*log(mmst1/mmu))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmsb1*mmst1*
     mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmst1/mmu))/(3.*(mmsb1 -
     mmsb2)*(mmst1 - mmst2)) + (32*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*
     log(mmgl/mmu)*log(mmst1/mmu))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*
     mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmst2/mmu))/((mmsb1 - mmsb2)*
     (mmst1 - mmst2)) + (32*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmst2/
     mmu))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*mmsb1*mmst2*mmt*DeltaInv(
     mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmst2/mmu))/(3.*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) - (32*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(
     mmst2/mmu))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*mmgl*mmsb1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) -
     (32*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu))/((mmsb1 - mmsb2)
     *(mmst1 - mmst2)) - (64*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmt/
     mmu))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) + (64*mmsb2*mmst1*mmt*DeltaInv(
     mmt,mmst1,mmgl)*log(mmt/mmu))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*
     mmgl*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu))/((mmsb1 - mmsb2)*(
     mmst1 - mmst2)) + (32*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)
     )/((mmsb1 - mmsb2)*(mmst1 - mmst2)) + (64*mmsb1*mmst2*mmt*DeltaInv(mmt,
     mmst2,mmgl)*log(mmt/mmu))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) - (64*mmsb2*
     mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu))/((mmsb1 - mmsb2)*(mmst1 -
     mmst2)) - (32*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(
     mmt/mmu))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*mmgl*mmsb2*mmt*
     DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu))/(3.*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) + (64*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/
     mmu)*log(mmt/mmu))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (64*mmsb2*mmst1*
     mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu))/(3.*(mmsb1 -
     mmsb2)*(mmst1 - mmst2)) + (32*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*log(
     mmgl/mmu)*log(mmt/mmu))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmgl*
     mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu))/(3.*(mmsb1
     - mmsb2)*(mmst1 - mmst2)) - (64*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*
     log(mmgl/mmu)*log(mmt/mmu))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (64*
     mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu))/(3.*(
     mmsb1 - mmsb2)*(mmst1 - mmst2)) - (224*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmsb1))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (224*mmst1*mmt*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmsb1))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) - (32*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1))/(3.*(mmsb1 -
     mmsb2)*(mmst1 - mmst2)) - (32*mmst1*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmsb1))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (224*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2))
     + (224*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/(3.*(-mmgl + mmsb1)
     *(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmsb1))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*mmst2*mmt*zt2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*
     (mmst1 - mmst2)) + (32*mmt*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,
     mmu)*pow2(mmsb1))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (
     32*mmt*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*pow2(mmsb1))/(3.*
     (-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*mmt*DeltaInv(mmt,
     mmst1,mmgl)*log(mmgl/mmu)*pow2(mmsb1))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) +
     (32*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*pow2(mmsb1))/((-mmgl
     + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmt*DeltaInv(mmt,mmst2,
     mmgl)*log(mmgl/mmu)*pow2(mmsb1))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*
     mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*pow2(mmsb1))/((-mmgl +
     mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmst1*mmt*DeltaInv(mmt,
     mmst1,mmgl)*log(mmst1/mmu)*pow2(mmsb1))/((-mmgl + mmsb1)*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) + (32*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*
     log(mmst1/mmu)*pow2(mmsb1))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) + (32*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmst2/mmu)*pow2(
     mmsb1))/((-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmst2/mmu)*pow2(mmsb1))/(3.*(-
     mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*mmt*DeltaInv(mmt,
     mmst1,mmgl)*log(mmt/mmu)*pow2(mmsb1))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) +
     (64*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*pow2(mmsb1))/((-mmgl +
     mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmt*DeltaInv(mmt,mmst2,mmgl)
     *log(mmt/mmu)*pow2(mmsb1))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) - (64*mmst2*
     mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*pow2(mmsb1))/((-mmgl + mmsb1)*(
     mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmt*DeltaInv(mmt,mmst1,mmgl)*log(
     mmgl/mmu)*log(mmt/mmu)*pow2(mmsb1))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) -
     (64*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(
     mmsb1))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*mmt*
     DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmsb1))/(3.*(
     mmsb1 - mmsb2)*(mmst1 - mmst2)) + (64*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*
     log(mmgl/mmu)*log(mmt/mmu)*pow2(mmsb1))/(3.*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)*(mmst1 - mmst2)) - (32*mmsb1*mmt*Fin3(mmt,mmsb1,mmst1,mmu))/(3.*(
     mmsb1 - mmsb2)*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (32*mmsb1*mmt*Fin3(
     mmt,mmsb1,mmst2,mmu))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)*pow2(-mmgl +
     mmsb1)) + (32*mmsb1*mmt*Fin3(mmt,mmst1,mmgl,mmu))/(3.*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb1)) - (32*mmsb1*mmt*Fin3(mmt,mmst2,mmgl,
     mmu))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)*pow2(-mmgl + mmsb1)) + (224*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2))
     + (224*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2))/(3.*(mmsb1 - mmsb2)
     *(-mmgl + mmsb2)*(mmst1 - mmst2)) + (32*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmsb2))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*mmst1*mmt*zt2*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*
     (mmst1 - mmst2)) - (224*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/(3.*(
     mmsb1 - mmsb2)*(mmst1 - mmst2)) - (224*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmsb2))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (32*
     mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/(3.*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) - (32*mmst2*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/(3.*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (32*mmt*DeltaInv(mmt,
     mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu)*pow2(mmsb2))/(3.*(mmsb1 - mmsb2)*(-
     mmgl + mmsb2)*(mmst1 - mmst2)) + (32*mmt*DeltaInv(mmt,mmst2,mmgl)*Fin3(
     mmt,mmst2,mmgl,mmu)*pow2(mmsb2))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(
     mmst1 - mmst2)) - (32*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*pow2(
     mmsb2))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmst1*mmt*DeltaInv(mmt,
     mmst1,mmgl)*log(mmgl/mmu)*pow2(mmsb2))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(
     mmst1 - mmst2)) + (32*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*pow2(
     mmsb2))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*mmst2*mmt*DeltaInv(mmt,
     mmst2,mmgl)*log(mmgl/mmu)*pow2(mmsb2))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(
     mmst1 - mmst2)) + (32*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmst1/mmu)*
     pow2(mmsb2))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (32*
     mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmst1/mmu)*pow2(
     mmsb2))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (32*mmst2*
     mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmst2/mmu)*pow2(mmsb2))/((mmsb1 - mmsb2)*
     (-mmgl + mmsb2)*(mmst1 - mmst2)) + (32*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*
     log(mmgl/mmu)*log(mmst2/mmu)*pow2(mmsb2))/(3.*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) - (32*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*
     pow2(mmsb2))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) - (64*mmst1*mmt*DeltaInv(
     mmt,mmst1,mmgl)*log(mmt/mmu)*pow2(mmsb2))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)
     *(mmst1 - mmst2)) + (32*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*pow2(
     mmsb2))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) + (64*mmst2*mmt*DeltaInv(mmt,
     mmst2,mmgl)*log(mmt/mmu)*pow2(mmsb2))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(
     mmst1 - mmst2)) + (32*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/
     mmu)*pow2(mmsb2))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (64*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmsb2))/(3.*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (32*mmt*DeltaInv(mmt,
     mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmsb2))/(3.*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) - (64*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*
     log(mmt/mmu)*pow2(mmsb2))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 -
     mmst2)) + (32*mmsb2*mmt*Fin3(mmt,mmsb2,mmst1,mmu))/(3.*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (32*mmsb2*mmt*Fin3(mmt,mmsb2,mmst2,
     mmu))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) - (32*
     mmsb2*mmt*Fin3(mmt,mmst1,mmgl,mmu))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)*
     pow2(-mmgl + mmsb2)) + (32*mmsb2*mmt*Fin3(mmt,mmst2,mmgl,mmu))/(3.*(mmsb1
     - mmsb2)*(mmst1 - mmst2)*pow2(-mmgl + mmsb2)) + (32*mmsb1*mmt*DeltaInv(
     mmt,mmst1,mmgl)*log(mmst1/mmu)*pow2(mmst1))/((-mmgl + mmsb1)*(mmsb1 -
     mmsb2)*(mmst1 - mmst2)) - (32*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*log(
     mmst1/mmu)*pow2(mmst1))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2))
     - (32*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmst1/mmu)*
     pow2(mmst1))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*
     mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmst1/mmu)*pow2(
     mmst1))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (32*mmsb1*
     mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*pow2(mmst1))/((-mmgl + mmsb1)*(
     mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*
     log(mmt/mmu)*pow2(mmst1))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)
     ) + (32*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*
     pow2(mmst1))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*
     mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmst1))
     /(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (32*mmsb1*mmt*
     DeltaInv(mmt,mmst2,mmgl)*log(mmst2/mmu)*pow2(mmst2))/((-mmgl + mmsb1)*(
     mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*
     log(mmst2/mmu)*pow2(mmst2))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 -
     mmst2)) + (32*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmst2/
     mmu)*pow2(mmst2))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (
     32*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmst2/mmu)*pow2(
     mmst2))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (32*mmsb1*
     mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*pow2(mmst2))/((-mmgl + mmsb1)*(
     mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*
     log(mmt/mmu)*pow2(mmst2))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)
     ) - (32*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*
     pow2(mmst2))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*
     mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmst2))
     /(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (224*mmsb1*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmt))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) -
     (224*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt))/(3.*(mmsb1 - mmsb2)*(mmst1
     - mmst2)) - (448*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt))/(3.*(-
     mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (448*mmsb2*mmst1*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmt))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(
     mmst1 - mmst2)) + (32*mmsb1*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt))/(3.*(
     mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmsb2*zt2*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmt))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (64*mmsb1*mmst1*zt2*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmt))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) + (64*mmsb2*mmst1*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt))/
     (3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (224*mmsb1*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmt))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) +
     (224*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt))/(3.*(mmsb1 - mmsb2)*(mmst1
     - mmst2)) + (448*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt))/(3.*(-
     mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (448*mmsb2*mmst2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmt))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(
     mmst1 - mmst2)) - (32*mmsb1*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt))/(3.*(
     mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*mmsb2*zt2*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmt))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (64*mmsb1*mmst2*zt2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmt))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) - (64*mmsb2*mmst2*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt))/
     (3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (32*mmsb1*DeltaInv(
     mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu)*pow2(mmt))/(3.*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*mmsb2*DeltaInv(mmt,mmst1,mmgl)*Fin3(
     mmt,mmst1,mmgl,mmu)*pow2(mmt))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1
     - mmst2)) + (32*mmsb1*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*
     pow2(mmt))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*
     mmsb2*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*pow2(mmt))/(3.*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (32*mmsb1*DeltaInv(mmt,
     mmst1,mmgl)*log(mmgl/mmu)*pow2(mmt))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) + (
     32*mmsb2*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*pow2(mmt))/((mmsb1 -
     mmsb2)*(mmst1 - mmst2)) + (32*mmsb1*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)
     *pow2(mmt))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmsb2*DeltaInv(mmt,
     mmst2,mmgl)*log(mmgl/mmu)*pow2(mmt))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) + (
     32*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*log(mmst1/mmu)*pow2(mmt))/((-mmgl
     + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmsb2*mmst1*DeltaInv(mmt,
     mmst1,mmgl)*log(mmst1/mmu)*pow2(mmt))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(
     mmst1 - mmst2)) + (32*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*
     log(mmst1/mmu)*pow2(mmt))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) - (32*mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(
     mmst1/mmu)*pow2(mmt))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2))
     - (32*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(mmst2/mmu)*pow2(mmt))/((-
     mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*mmsb2*mmst2*DeltaInv(
     mmt,mmst2,mmgl)*log(mmst2/mmu)*pow2(mmt))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)
     *(mmst1 - mmst2)) - (32*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)
     *log(mmst2/mmu)*pow2(mmt))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) + (32*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(
     mmst2/mmu)*pow2(mmt))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2))
     - (32*mmsb1*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*pow2(mmt))/((mmsb1 -
     mmsb2)*(mmst1 - mmst2)) + (32*mmsb2*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*
     pow2(mmt))/((mmsb1 - mmsb2)*(mmst1 - mmst2)) + (96*mmsb1*mmst1*DeltaInv(
     mmt,mmst1,mmgl)*log(mmt/mmu)*pow2(mmt))/((-mmgl + mmsb1)*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) - (96*mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*
     pow2(mmt))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (32*mmsb1*
     DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*pow2(mmt))/((mmsb1 - mmsb2)*(mmst1 -
     mmst2)) - (32*mmsb2*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*pow2(mmt))/((
     mmsb1 - mmsb2)*(mmst1 - mmst2)) - (96*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)
     *log(mmt/mmu)*pow2(mmt))/((-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2))
     + (96*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*pow2(mmt))/((mmsb1
     - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (32*mmsb1*DeltaInv(mmt,mmst1,
     mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmt))/(3.*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) - (32*mmsb2*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*
     pow2(mmt))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmsb1*mmst1*
     DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmt))/(3.*(-mmgl
     + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*mmsb2*mmst1*DeltaInv(mmt,
     mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmt))/(3.*(mmsb1 - mmsb2)*(-
     mmgl + mmsb2)*(mmst1 - mmst2)) - (32*mmsb1*DeltaInv(mmt,mmst2,mmgl)*log(
     mmgl/mmu)*log(mmt/mmu)*pow2(mmt))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (
     32*mmsb2*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmt))/(
     3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*mmsb1*mmst2*DeltaInv(mmt,mmst2,
     mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmt))/(3.*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)*(mmst1 - mmst2)) - (32*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(
     mmgl/mmu)*log(mmt/mmu)*pow2(mmt))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(
     mmst1 - mmst2)) - (64*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*log(mmst1/mmu)*
     log(mmt/mmu)*pow2(mmt))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) + (64*mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)*log(mmst1/mmu)*log(mmt/
     mmu)*pow2(mmt))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (
     64*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(mmst2/mmu)*log(mmt/mmu)*pow2(
     mmt))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (64*mmsb2*
     mmst2*DeltaInv(mmt,mmst2,mmgl)*log(mmst2/mmu)*log(mmt/mmu)*pow2(mmt))/(3.*
     (mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (224*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmsb1)*pow2(mmt))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) - (32*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(mmt))/
     (3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (224*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmsb1)*pow2(mmt))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) + (32*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(mmt))/
     (3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*DeltaInv(mmt,
     mmst1,mmgl)*log(mmgl/mmu)*pow2(mmsb1)*pow2(mmt))/((-mmgl + mmsb1)*(mmsb1 -
     mmsb2)*(mmst1 - mmst2)) - (32*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*pow2(
     mmsb1)*pow2(mmt))/((-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*
     DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*pow2(mmsb1)*pow2(mmt))/((-mmgl +
     mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*DeltaInv(mmt,mmst2,mmgl)*
     log(mmt/mmu)*pow2(mmsb1)*pow2(mmt))/((-mmgl + mmsb1)*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) - (32*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*
     pow2(mmsb1)*pow2(mmt))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)
     ) + (32*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmsb1)*
     pow2(mmt))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (224*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(mmt))/(3.*(mmsb1 - mmsb2)*(-mmgl
     + mmsb2)*(mmst1 - mmst2)) + (32*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*
     pow2(mmt))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (224*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(mmt))/(3.*(mmsb1 - mmsb2)*(-mmgl
     + mmsb2)*(mmst1 - mmst2)) - (32*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*
     pow2(mmt))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (32*
     DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*pow2(mmsb2)*pow2(mmt))/((mmsb1 -
     mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (32*DeltaInv(mmt,mmst2,mmgl)*
     log(mmgl/mmu)*pow2(mmsb2)*pow2(mmt))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(
     mmst1 - mmst2)) - (32*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*pow2(mmsb2)*
     pow2(mmt))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (32*
     DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*pow2(mmsb2)*pow2(mmt))/((mmsb1 -
     mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (32*DeltaInv(mmt,mmst1,mmgl)*
     log(mmgl/mmu)*log(mmt/mmu)*pow2(mmsb2)*pow2(mmt))/(3.*(mmsb1 - mmsb2)*(-
     mmgl + mmsb2)*(mmst1 - mmst2)) - (32*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/
     mmu)*log(mmt/mmu)*pow2(mmsb2)*pow2(mmt))/(3.*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) - (16*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow2(log(mmgl/mmu)))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*mmgl*
     mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmgl/mmu)))/(3.*(mmsb1 -
     mmsb2)*(mmst1 - mmst2)) + (16*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow2(log(mmgl/mmu)))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (16*mmsb2*
     mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmgl/mmu)))/(3.*(mmsb1 -
     mmsb2)*(mmst1 - mmst2)) + (16*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*
     pow2(log(mmgl/mmu)))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (16*mmgl*
     mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmgl/mmu)))/(3.*(mmsb1 -
     mmsb2)*(mmst1 - mmst2)) - (16*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*
     pow2(log(mmgl/mmu)))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*mmsb2*
     mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmgl/mmu)))/(3.*(mmsb1 -
     mmsb2)*(mmst1 - mmst2)) - (16*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*
     pow2(log(mmgl/mmu)))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (16*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(log(mmgl/mmu)))/(3.*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*mmt*DeltaInv(mmt,mmst2,mmgl)
     *pow2(mmsb1)*pow2(log(mmgl/mmu)))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (
     16*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(log(mmgl/mmu)))/(3.
     *(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmsb2)*pow2(log(mmgl/mmu)))/(3.*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) + (16*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(log(
     mmgl/mmu)))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (16*
     mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(log(mmgl/mmu)))/(3.*(mmsb1 -
     mmsb2)*(mmst1 - mmst2)) - (16*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmsb2)*pow2(log(mmgl/mmu)))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 -
     mmst2)) + (16*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow2(log(mmgl/mmu))
     )/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (16*mmsb2*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmt)*pow2(log(mmgl/mmu)))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2))
     - (16*mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow2(log(mmgl/mmu)))/(3.*(
     mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmt)*pow2(log(mmgl/mmu)))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (16*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(mmt)*pow2(log(mmgl/mmu)))/(3.*(-
     mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmsb1)*pow2(mmt)*pow2(log(mmgl/mmu)))/(3.*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)
     *pow2(mmt)*pow2(log(mmgl/mmu)))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1
     - mmst2)) - (16*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(mmt)*pow2(log(
     mmgl/mmu)))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (16*
     mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmst1/mmu)))/(3.*(mmsb1
     - mmsb2)*(mmst1 - mmst2)) + (16*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow2(log(mmst1/mmu)))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*mmst1*
     mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(log(mmst1/mmu)))/(3.*(-mmgl
     + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (16*mmst1*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmsb2)*pow2(log(mmst1/mmu)))/(3.*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) - (16*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmst1)*pow2(log(mmst1/mmu)))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) + (16*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow2(log(
     mmst1/mmu)))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (16*
     mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow2(log(mmst1/mmu)))/(3.*(
     -mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*mmsb2*mmst1*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow2(log(mmst1/mmu)))/(3.*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (16*mmsb1*mmst2*mmt*DeltaInv(
     mmt,mmst2,mmgl)*pow2(log(mmst2/mmu)))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2))
     - (16*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmst2/mmu)))/(3.*(
     mmsb1 - mmsb2)*(mmst1 - mmst2)) - (16*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmsb1)*pow2(log(mmst2/mmu)))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) + (16*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(
     log(mmst2/mmu)))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (
     16*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow2(log(mmst2/mmu)))/(
     3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (16*mmsb2*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow2(log(mmst2/mmu)))/(3.*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (16*mmsb1*mmst2*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmt)*pow2(log(mmst2/mmu)))/(3.*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)*(mmst1 - mmst2)) - (16*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmt)*pow2(log(mmst2/mmu)))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 -
     mmst2)) - (16*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmt/mmu)))/
     (3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*mmgl*mmsb2*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow2(log(mmt/mmu)))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (
     32*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmt/mmu)))/(3.*(mmsb1
     - mmsb2)*(mmst1 - mmst2)) - (32*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow2(log(mmt/mmu)))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*mmgl*mmsb1*
     mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmt/mmu)))/(3.*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) - (16*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(
     mmt/mmu)))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmsb1*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmt/mmu)))/(3.*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) + (32*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmt/mmu)))
     /(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (16*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmsb1)*pow2(log(mmt/mmu)))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (
     32*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(log(mmt/mmu)))/(3.*
     (-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmsb1)*pow2(log(mmt/mmu)))/(3.*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) + (32*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(log(mmt/
     mmu)))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(log(mmt/mmu)))/(3.*(mmsb1 -
     mmsb2)*(mmst1 - mmst2)) + (32*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmsb2)*pow2(log(mmt/mmu)))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 -
     mmst2)) - (16*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(log(mmt/mmu)))
     /(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmst2*mmt*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmsb2)*pow2(log(mmt/mmu)))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*
     (mmst1 - mmst2)) + (16*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*
     pow2(log(mmt/mmu)))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) -
     (16*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow2(log(mmt/mmu)))/(3.
     *(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (16*mmsb1*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow2(log(mmt/mmu)))/(3.*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*mmsb2*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmst2)*pow2(log(mmt/mmu)))/(3.*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)*(mmst1 - mmst2)) + (16*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*
     pow2(log(mmt/mmu)))/(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (16*mmsb2*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow2(log(mmt/mmu)))/(3.*(mmsb1 - mmsb2)
     *(mmst1 - mmst2)) - (16*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*
     pow2(log(mmt/mmu)))/((-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (
     16*mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmt)*pow2(log(mmt/mmu)))/((
     mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (16*mmsb1*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmt)*pow2(log(mmt/mmu)))/(3.*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) + (16*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*pow2(log(mmt/mmu)))
     /(3.*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*mmsb1*mmst2*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmt)*pow2(log(mmt/mmu)))/((-mmgl + mmsb1)*(mmsb1 - mmsb2)
     *(mmst1 - mmst2)) - (16*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmt)*
     pow2(log(mmt/mmu)))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (
     16*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(mmt)*pow2(log(mmt/mmu)))/(3.*
     (-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmsb1)*pow2(mmt)*pow2(log(mmt/mmu)))/(3.*(-mmgl + mmsb1)*(mmsb1
     - mmsb2)*(mmst1 - mmst2)) + (16*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(
     mmt)*pow2(log(mmt/mmu)))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 -
     mmst2)) - (16*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(mmt)*pow2(log(mmt/
     mmu)))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (224*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*
     (mmst1 - mmst2)) + (32*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1))/(3.*(
     -mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (224*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow3(mmsb1))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) - (32*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1))/(3.*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmt*DeltaInv(mmt,mmst1,mmgl)
     *log(mmgl/mmu)*pow3(mmsb1))/((-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) + (32*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*pow3(mmsb1))/((-
     mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmt*DeltaInv(mmt,
     mmst1,mmgl)*log(mmt/mmu)*pow3(mmsb1))/((-mmgl + mmsb1)*(mmsb1 - mmsb2)*(
     mmst1 - mmst2)) + (32*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*pow3(
     mmsb1))/((-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (32*mmt*
     DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow3(mmsb1))/(3.*(-
     mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) - (32*mmt*DeltaInv(mmt,
     mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow3(mmsb1))/(3.*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     log(mmgl/mmu))*pow3(mmsb1))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 -
     mmst2)) - (16*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmgl/mmu))*pow3(mmsb1)
     )/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) + (16*mmt*DeltaInv(
     mmt,mmst1,mmgl)*pow2(log(mmt/mmu))*pow3(mmsb1))/(3.*(-mmgl + mmsb1)*(mmsb1
     - mmsb2)*(mmst1 - mmst2)) - (16*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmt/
     mmu))*pow3(mmsb1))/(3.*(-mmgl + mmsb1)*(mmsb1 - mmsb2)*(mmst1 - mmst2)) -
     (224*mmt*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2))/(3.*(mmsb1 - mmsb2)*(-mmgl
     + mmsb2)*(mmst1 - mmst2)) - (32*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow3(
     mmsb2))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (224*mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*
     (mmst1 - mmst2)) + (32*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2))/(3.*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (32*mmt*DeltaInv(mmt,
     mmst1,mmgl)*log(mmgl/mmu)*pow3(mmsb2))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(
     mmst1 - mmst2)) - (32*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*pow3(
     mmsb2))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (32*mmt*
     DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*pow3(mmsb2))/((mmsb1 - mmsb2)*(-mmgl
     + mmsb2)*(mmst1 - mmst2)) - (32*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*
     pow3(mmsb2))/((mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (32*mmt*
     DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow3(mmsb2))/(3.*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) + (32*mmt*DeltaInv(mmt,
     mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow3(mmsb2))/(3.*(mmsb1 - mmsb2)*(-
     mmgl + mmsb2)*(mmst1 - mmst2)) - (16*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     log(mmgl/mmu))*pow3(mmsb2))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 -
     mmst2)) + (16*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmgl/mmu))*pow3(mmsb2)
     )/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2)) - (16*mmt*DeltaInv(
     mmt,mmst1,mmgl)*pow2(log(mmt/mmu))*pow3(mmsb2))/(3.*(mmsb1 - mmsb2)*(-mmgl
     + mmsb2)*(mmst1 - mmst2)) + (16*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmt/
     mmu))*pow3(mmsb2))/(3.*(mmsb1 - mmsb2)*(-mmgl + mmsb2)*(mmst1 - mmst2))) -
     (16*mmgl*mmsb2*mmt*Fin3(mmt,mmsb2,mmst1,mmu))/(3.*mgl*(mmst1 - mmst2)*
     pow3(-mmgl + mmsb2)) + (16*mmgl*mmsb2*mmt*Fin3(mmt,mmsb2,mmst2,mmu))/(3.*
     mgl*(mmst1 - mmst2)*pow3(-mmgl + mmsb2)) + (16*mmgl*mmsb2*mmt*Fin3(mmt,
     mmst1,mmgl,mmu))/(3.*mgl*(mmst1 - mmst2)*pow3(-mmgl + mmsb2)) - (16*mmgl*
     mmsb2*mmt*Fin3(mmt,mmst2,mmgl,mmu))/(3.*mgl*(mmst1 - mmst2)*pow3(-mmgl +
     mmsb2))) + (4*Fin20(mmgl,mmsusy,mmu)*(2/(-mmgl + mmsb1) + 2/(-mmgl +
     mmsb2) + (6*mmsb1)/pow2(-mmgl + mmsb1) - (6*mmsusy)/pow2(-mmgl + mmsb1) +
     (6*mmsb2)/pow2(-mmgl + mmsb2) - (6*mmsusy)/pow2(-mmgl + mmsb2) + (8*mmsb1*
     mmsusy)/pow3(-mmgl + mmsb1) - (8*pow2(mmsb1))/pow3(-mmgl + mmsb1) + (8*
     mmsb2*mmsusy)/pow3(-mmgl + mmsb2) - (8*pow2(mmsb2))/pow3(-mmgl + mmsb2)))/
     3. + (4*Fin3(mmt,mmst1,mmgl,mmu)*(-3/(4.*(-mmgl + mmsb1)) - 3/(4.*(-mmgl +
     mmsb2)) - 2*mmgl*DeltaInv(mmt,mmst1,mmgl) - 2*mmsb1*DeltaInv(mmt,mmst1,
     mmgl) - 2*mmsb2*DeltaInv(mmt,mmst1,mmgl) + 4*mmst1*DeltaInv(mmt,mmst1,
     mmgl) - (4*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl))/(-mmgl + mmsb1) - (4*
     mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl))/(-mmgl + mmsb2) + 2*mmt*DeltaInv(
     mmt,mmst1,mmgl) - (2*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl))/(-mmgl + mmsb1) -
     (2*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl))/(-mmgl + mmsb2) - (mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl))/(-mmgl + mmsb1) - (mmst1*mmt*DeltaInv(mmt,mmst1,
     mmgl))/(-mmgl + mmsb2) + (3*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1))/(-mmgl +
     mmsb1) + (7*mmsb1)/(4.*pow2(-mmgl + mmsb1)) - (3*mmst1)/(4.*pow2(-mmgl +
     mmsb1)) + (3*mmt)/(4.*pow2(-mmgl + mmsb1)) + (mmsb1*mmst1*mmt*DeltaInv(
     mmt,mmst1,mmgl))/pow2(-mmgl + mmsb1) + (2*mmst1*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmsb1))/pow2(-mmgl + mmsb1) + (mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmsb1))/pow2(-mmgl + mmsb1) + (3*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2))/(-
     mmgl + mmsb2) + (7*mmsb2)/(4.*pow2(-mmgl + mmsb2)) - (3*mmst1)/(4.*pow2(-
     mmgl + mmsb2)) + (3*mmt)/(4.*pow2(-mmgl + mmsb2)) + (mmsb2*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl))/pow2(-mmgl + mmsb2) + (2*mmst1*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmsb2))/pow2(-mmgl + mmsb2) + (mmt*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmsb2))/pow2(-mmgl + mmsb2) + (DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmst1))/(-mmgl + mmsb1) + (DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(-mmgl +
     mmsb2) - (mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/pow2(-mmgl + mmsb1)
     - (mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/pow2(-mmgl + mmsb2) - (
     DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1))/pow2(-mmgl + mmsb1) + (mmsb1*mmst1)/
     pow3(-mmgl + mmsb1) - (mmsb1*mmt)/pow3(-mmgl + mmsb1) - pow2(mmsb1)/pow3(-
     mmgl + mmsb1) - (DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2))/pow2(-mmgl + mmsb2)
     + (mmsb2*mmst1)/pow3(-mmgl + mmsb2) - (mmsb2*mmt)/pow3(-mmgl + mmsb2) -
     pow2(mmsb2)/pow3(-mmgl + mmsb2)))/3. + (4*Fin3(mmt,mmst2,mmgl,mmu)*(-3/(4.
     *(-mmgl + mmsb1)) - 3/(4.*(-mmgl + mmsb2)) - 2*mmgl*DeltaInv(mmt,mmst2,
     mmgl) - 2*mmsb1*DeltaInv(mmt,mmst2,mmgl) - 2*mmsb2*DeltaInv(mmt,mmst2,
     mmgl) + 4*mmst2*DeltaInv(mmt,mmst2,mmgl) - (4*mmsb1*mmst2*DeltaInv(mmt,
     mmst2,mmgl))/(-mmgl + mmsb1) - (4*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl))/(-
     mmgl + mmsb2) + 2*mmt*DeltaInv(mmt,mmst2,mmgl) - (2*mmsb1*mmt*DeltaInv(
     mmt,mmst2,mmgl))/(-mmgl + mmsb1) - (2*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl))/
     (-mmgl + mmsb2) - (mmst2*mmt*DeltaInv(mmt,mmst2,mmgl))/(-mmgl + mmsb1) - (
     mmst2*mmt*DeltaInv(mmt,mmst2,mmgl))/(-mmgl + mmsb2) + (3*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmsb1))/(-mmgl + mmsb1) + (7*mmsb1)/(4.*pow2(-mmgl +
     mmsb1)) - (3*mmst2)/(4.*pow2(-mmgl + mmsb1)) + (3*mmt)/(4.*pow2(-mmgl +
     mmsb1)) + (mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl))/pow2(-mmgl + mmsb1) +
     (2*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/pow2(-mmgl + mmsb1) + (mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/pow2(-mmgl + mmsb1) + (3*DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmsb2))/(-mmgl + mmsb2) + (7*mmsb2)/(4.*pow2(-mmgl +
     mmsb2)) - (3*mmst2)/(4.*pow2(-mmgl + mmsb2)) + (3*mmt)/(4.*pow2(-mmgl +
     mmsb2)) + (mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl))/pow2(-mmgl + mmsb2) +
     (2*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/pow2(-mmgl + mmsb2) + (mmt*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/pow2(-mmgl + mmsb2) + (DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmst2))/(-mmgl + mmsb1) + (DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmst2))/(-mmgl + mmsb2) - (mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/
     pow2(-mmgl + mmsb1) - (mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/pow2(-
     mmgl + mmsb2) - (DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1))/pow2(-mmgl + mmsb1)
     + (mmsb1*mmst2)/pow3(-mmgl + mmsb1) - (mmsb1*mmt)/pow3(-mmgl + mmsb1) -
     pow2(mmsb1)/pow3(-mmgl + mmsb1) - (DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2))/
     pow2(-mmgl + mmsb2) + (mmsb2*mmst2)/pow3(-mmgl + mmsb2) - (mmsb2*mmt)/
     pow3(-mmgl + mmsb2) - pow2(mmsb2)/pow3(-mmgl + mmsb2)))/3. + (4*Fin3(mmt,
     mmsb2,mmst1,mmu)*(-((mmsb2*mmst1*DeltaInv(mmt,mmsb2,mmst1))/(-mmgl +
     mmsb2)) - (mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst1))/(2.*(-mmgl + mmsb2)) - (
     mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1))/(2.*(-mmgl + mmsb2)) + (DeltaInv(mmt,
     mmsb2,mmst1)*pow2(mmsb2))/(2.*(-mmgl + mmsb2)) - (3*mmsb2)/(4.*pow2(-mmgl
     + mmsb2)) + (3*mmst1)/(4.*pow2(-mmgl + mmsb2)) - (3*mmt)/(4.*pow2(-mmgl +
     mmsb2)) + (mmsb2*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1))/pow2(-mmgl + mmsb2)
     + (2*mmst1*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2))/pow2(-mmgl + mmsb2) + (
     mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2))/pow2(-mmgl + mmsb2) + (
     DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/(2.*(-mmgl + mmsb2)) - (mmsb2*
     DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/pow2(-mmgl + mmsb2) - (DeltaInv(
     mmt,mmsb2,mmst1)*pow3(mmsb2))/pow2(-mmgl + mmsb2) - (mmsb2*mmst1)/pow3(-
     mmgl + mmsb2) + (mmsb2*mmt)/pow3(-mmgl + mmsb2) + pow2(mmsb2)/pow3(-mmgl +
     mmsb2)))/3. + (4*Fin3(mmt,mmsb2,mmst2,mmu)*(-((mmsb2*mmst2*DeltaInv(mmt,
     mmsb2,mmst2))/(-mmgl + mmsb2)) - (mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst2))/(2.
     *(-mmgl + mmsb2)) - (mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2))/(2.*(-mmgl +
     mmsb2)) + (DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2))/(2.*(-mmgl + mmsb2)) - (
     3*mmsb2)/(4.*pow2(-mmgl + mmsb2)) + (3*mmst2)/(4.*pow2(-mmgl + mmsb2)) - (
     3*mmt)/(4.*pow2(-mmgl + mmsb2)) + (mmsb2*mmst2*mmt*DeltaInv(mmt,mmsb2,
     mmst2))/pow2(-mmgl + mmsb2) + (2*mmst2*DeltaInv(mmt,mmsb2,mmst2)*pow2(
     mmsb2))/pow2(-mmgl + mmsb2) + (mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2))/
     pow2(-mmgl + mmsb2) + (DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2))/(2.*(-mmgl +
     mmsb2)) - (mmsb2*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2))/pow2(-mmgl +
     mmsb2) - (DeltaInv(mmt,mmsb2,mmst2)*pow3(mmsb2))/pow2(-mmgl + mmsb2) - (
     mmsb2*mmst2)/pow3(-mmgl + mmsb2) + (mmsb2*mmt)/pow3(-mmgl + mmsb2) + pow2(
     mmsb2)/pow3(-mmgl + mmsb2)))/3. + (4*Fin20(mmsb2,mmsusy,mmu)*(4/(-mmgl +
     mmsb2) - (14*mmsb2)/pow2(-mmgl + mmsb2) + (6*mmsusy)/pow2(-mmgl + mmsb2) -
     (8*mmsb2*mmsusy)/pow3(-mmgl + mmsb2) + (8*pow2(mmsb2))/pow3(-mmgl + mmsb2)
     ))/3. + (4*log(mmsb1/mmu)*log(mmsb2/mmu)*((3*mmsb1*mmsb2)/(4.*pow2(-mmgl +
     mmsb1)) + (3*mmsb1*mmsb2)/(4.*pow2(-mmgl + mmsb2)) - (mmsb2*pow2(mmsb1))/
     pow3(-mmgl + mmsb1) - (mmsb1*pow2(mmsb2))/pow3(-mmgl + mmsb2)))/3. + (4*
     log(mmsb2/mmu)*log(mmst1/mmu)*((mmsb2*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1))
     /(-mmgl + mmsb2) + (3*mmsb2*mmst1)/(4.*pow2(-mmgl + mmsb2)) - (2*mmst1*
     mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2))/pow2(-mmgl + mmsb2) - (mmst1*
     pow2(mmsb2))/pow3(-mmgl + mmsb2)))/3. + (4*log(mmgl/mmu)*log(mmst1/mmu)*((
     3*mmst1)/(4.*(-mmgl + mmsb1)) + (3*mmst1)/(4.*(-mmgl + mmsb2)) - 4*mmst1*
     mmt*DeltaInv(mmt,mmst1,mmgl) + (4*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)
     )/(-mmgl + mmsb1) + (4*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl))/(-mmgl +
     mmsb2) - (7*mmsb1*mmst1)/(4.*pow2(-mmgl + mmsb1)) - (2*mmst1*mmt*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmsb1))/pow2(-mmgl + mmsb1) - (7*mmsb2*mmst1)/(4.*
     pow2(-mmgl + mmsb2)) - (2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2))/
     pow2(-mmgl + mmsb2) + (mmst1*pow2(mmsb1))/pow3(-mmgl + mmsb1) + (mmst1*
     pow2(mmsb2))/pow3(-mmgl + mmsb2)))/3. + (4*log(mmsb2/mmu)*log(mmst2/mmu)*(
     (mmsb2*mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2))/(-mmgl + mmsb2) + (3*mmsb2*
     mmst2)/(4.*pow2(-mmgl + mmsb2)) - (2*mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2)*
     pow2(mmsb2))/pow2(-mmgl + mmsb2) - (mmst2*pow2(mmsb2))/pow3(-mmgl + mmsb2)
     ))/3. + (4*log(mmgl/mmu)*log(mmst2/mmu)*((3*mmst2)/(4.*(-mmgl + mmsb1)) +
     (3*mmst2)/(4.*(-mmgl + mmsb2)) - 4*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl) + (
     4*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl))/(-mmgl + mmsb1) + (4*mmsb2*
     mmst2*mmt*DeltaInv(mmt,mmst2,mmgl))/(-mmgl + mmsb2) - (7*mmsb1*mmst2)/(4.*
     pow2(-mmgl + mmsb1)) - (2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/
     pow2(-mmgl + mmsb1) - (7*mmsb2*mmst2)/(4.*pow2(-mmgl + mmsb2)) - (2*mmst2*
     mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/pow2(-mmgl + mmsb2) + (mmst2*
     pow2(mmsb1))/pow3(-mmgl + mmsb1) + (mmst2*pow2(mmsb2))/pow3(-mmgl + mmsb2)
     ))/3. + (4*log(mmsb2/mmu)*log(mmsusy/mmu)*((6*mmsb2*mmsusy)/pow2(-mmgl +
     mmsb2) - (8*mmsusy*pow2(mmsb2))/pow3(-mmgl + mmsb2)))/3. + (4*log(mmgl/
     mmu)*log(mmsusy/mmu)*((6*mmsusy)/(-mmgl + mmsb1) + (6*mmsusy)/(-mmgl +
     mmsb2) - (14*mmsb1*mmsusy)/pow2(-mmgl + mmsb1) - (14*mmsb2*mmsusy)/pow2(-
     mmgl + mmsb2) + (8*mmsusy*pow2(mmsb1))/pow3(-mmgl + mmsb1) + (8*mmsusy*
     pow2(mmsb2))/pow3(-mmgl + mmsb2)))/3. + (4*mmsb1*mmsb2*Fin20(mmsb1,mmgl,
     mmu))/(3.*pow3(-mmgl + mmsb2)) - (4*mmsb1*mmsb2*Fin20(mmsb1,mmsb2,mmu))/(
     3.*pow3(-mmgl + mmsb2)) - (4*Fin20(mmsb1,mmgl,mmu)*pow2(mmsb2))/(3.*pow3(-
     mmgl + mmsb2)) + (4*Fin20(mmsb1,mmsb2,mmu)*pow2(mmsb2))/(3.*pow3(-mmgl +
     mmsb2)) - (16*Fin20(mmsb2,mmgl,mmu)*pow2(mmsb2))/(3.*pow3(-mmgl + mmsb2))
     - (8*mmsb1*log(mmgl/mmu)*pow2(mmsb2))/(3.*pow3(-mmgl + mmsb2)) - (8*mmst1*
     log(mmgl/mmu)*pow2(mmsb2))/(3.*pow3(-mmgl + mmsb2)) - (8*mmst2*log(mmgl/
     mmu)*pow2(mmsb2))/(3.*pow3(-mmgl + mmsb2)) - (64*mmsusy*log(mmgl/mmu)*
     pow2(mmsb2))/(3.*pow3(-mmgl + mmsb2)) + (16*mmt*log(mmgl/mmu)*pow2(mmsb2))
     /(3.*pow3(-mmgl + mmsb2)) + (4*mmsb1*log(mmgl/mmu)*log(mmsb1/mmu)*pow2(
     mmsb2))/(3.*pow3(-mmgl + mmsb2)) + (8*mmsb1*log(mmsb2/mmu)*pow2(mmsb2))/(
     3.*pow3(-mmgl + mmsb2)) + (8*mmst1*log(mmsb2/mmu)*pow2(mmsb2))/(3.*pow3(-
     mmgl + mmsb2)) + (8*mmst2*log(mmsb2/mmu)*pow2(mmsb2))/(3.*pow3(-mmgl +
     mmsb2)) + (64*mmsusy*log(mmsb2/mmu)*pow2(mmsb2))/(3.*pow3(-mmgl + mmsb2))
     - (16*mmt*log(mmsb2/mmu)*pow2(mmsb2))/(3.*pow3(-mmgl + mmsb2)) + (2*mmsb1*
     pow2(mmsb2)*pow2(log(mmgl/mmu)))/(3.*pow3(-mmgl + mmsb2)) + (2*mmst1*pow2(
     mmsb2)*pow2(log(mmgl/mmu)))/(3.*pow3(-mmgl + mmsb2)) + (2*mmst2*pow2(
     mmsb2)*pow2(log(mmgl/mmu)))/(3.*pow3(-mmgl + mmsb2)) + (16*mmsusy*pow2(
     mmsb2)*pow2(log(mmgl/mmu)))/(3.*pow3(-mmgl + mmsb2)) - (4*mmt*pow2(mmsb2)*
     pow2(log(mmgl/mmu)))/(3.*pow3(-mmgl + mmsb2)) - (2*mmsb1*pow2(mmsb2)*pow2(
     log(mmsb2/mmu)))/(3.*pow3(-mmgl + mmsb2)) - (2*mmst1*pow2(mmsb2)*pow2(log(
     mmsb2/mmu)))/(3.*pow3(-mmgl + mmsb2)) - (2*mmst2*pow2(mmsb2)*pow2(log(
     mmsb2/mmu)))/(3.*pow3(-mmgl + mmsb2)) - (16*mmsusy*pow2(mmsb2)*pow2(log(
     mmsb2/mmu)))/(3.*pow3(-mmgl + mmsb2)) + (4*mmt*pow2(mmsb2)*pow2(log(mmsb2/
     mmu)))/(3.*pow3(-mmgl + mmsb2)) - (112*pow3(mmsb2))/(3.*pow3(-mmgl +
     mmsb2)) - (16*zt2*pow3(mmsb2))/(3.*pow3(-mmgl + mmsb2)) + (508*log(mmgl/
     mmu)*pow3(mmsb2))/(9.*pow3(-mmgl + mmsb2)) - (220*log(mmsb2/mmu)*pow3(
     mmsb2))/(9.*pow3(-mmgl + mmsb2)) + (116*log(mmgl/mmu)*log(mmsb2/mmu)*pow3(
     mmsb2))/(9.*pow3(-mmgl + mmsb2)) - (242*pow2(log(mmgl/mmu))*pow3(mmsb2))/(
     9.*pow3(-mmgl + mmsb2)) + (10*pow2(log(mmsb2/mmu))*pow3(mmsb2))/(3.*pow3(-
     mmgl + mmsb2)) + (14*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmst1))/(3.*(-mmgl +
     mmsb1)) + (2*zt2*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmst1))/(3.*(-mmgl +
     mmsb1)) + (14*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmst1))/(3.*(-mmgl + mmsb2))
     + (2*zt2*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmst1))/(3.*(-mmgl + mmsb2)) + (
     28*DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1))/(3.*(-mmgl + mmsb1)) + (28*
     DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1))/(3.*(-mmgl + mmsb2)) + (4*zt2*
     DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1))/(3.*(-mmgl + mmsb1)) + (4*zt2*
     DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1))/(3.*(-mmgl + mmsb2)) - (28*mmsb1*
     DeltaInv(mmt,mmsb1,mmst1)*pow3(mmst1))/(3.*pow2(-mmgl + mmsb1)) - (4*
     mmsb1*zt2*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmst1))/(3.*pow2(-mmgl + mmsb1))
     - (28*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1))/(3.*pow2(-mmgl + mmsb1))
     - (4*mmsb1*zt2*DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1))/(3.*pow2(-mmgl +
     mmsb1)) - (28*mmsb2*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmst1))/(3.*pow2(-mmgl
     + mmsb2)) - (4*mmsb2*zt2*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmst1))/(3.*pow2(-
     mmgl + mmsb2)) - (28*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1))/(3.*pow2(
     -mmgl + mmsb2)) - (4*mmsb2*zt2*DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1))/(3.*
     pow2(-mmgl + mmsb2)) + (4*log(mmst1/mmu)*log(mmt/mmu)*((-3*mmst1)/(2.*(-
     mmgl + mmsb1)) - (3*mmst1)/(2.*(-mmgl + mmsb2)) - (3*mmsb1*mmst1*mmt*
     DeltaInv(mmt,mmsb1,mmst1))/(2.*(-mmgl + mmsb1)) - (3*mmsb2*mmst1*mmt*
     DeltaInv(mmt,mmsb2,mmst1))/(2.*(-mmgl + mmsb2)) - 2*mmgl*mmst1*DeltaInv(
     mmt,mmst1,mmgl) - 2*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl) - 2*mmsb2*mmst1*
     DeltaInv(mmt,mmst1,mmgl) + 6*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl) - (6*
     mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl))/(-mmgl + mmsb1) - (6*mmsb2*
     mmst1*mmt*DeltaInv(mmt,mmst1,mmgl))/(-mmgl + mmsb2) + (mmst1*DeltaInv(mmt,
     mmsb1,mmst1)*pow2(mmsb1))/(2.*(-mmgl + mmsb1)) + (3*mmst1*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmsb1))/(-mmgl + mmsb1) + (2*mmsb1*mmst1)/pow2(-mmgl +
     mmsb1) + (3*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1))/pow2(-mmgl +
     mmsb1) + (3*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1))/pow2(-mmgl +
     mmsb1) + (mmst1*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2))/(2.*(-mmgl + mmsb2)
     ) + (3*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2))/(-mmgl + mmsb2) + (2*
     mmsb2*mmst1)/pow2(-mmgl + mmsb2) + (3*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1)*
     pow2(mmsb2))/pow2(-mmgl + mmsb2) + (3*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmsb2))/pow2(-mmgl + mmsb2) - (mmsb1*DeltaInv(mmt,mmsb1,mmst1)*pow2(
     mmst1))/(-mmgl + mmsb1) - (mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1))/(2.*
     (-mmgl + mmsb1)) - (mmsb2*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/(-mmgl +
     mmsb2) - (mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/(2.*(-mmgl + mmsb2))
     + 4*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1) - (4*mmsb1*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmst1))/(-mmgl + mmsb1) - (4*mmsb2*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmst1))/(-mmgl + mmsb2) - (mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/
     (-mmgl + mmsb1) - (mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(-mmgl +
     mmsb2) + (mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1))/pow2(-mmgl +
     mmsb1) + (mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/pow2(-mmgl +
     mmsb1) + (2*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow2(mmst1))/pow2(-mmgl
     + mmsb1) + (2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(mmst1))/pow2(-mmgl
     + mmsb1) + (mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/pow2(-mmgl +
     mmsb2) + (mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/pow2(-mmgl +
     mmsb2) + (2*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2)*pow2(mmst1))/pow2(-mmgl
     + mmsb2) + (2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(mmst1))/pow2(-mmgl
     + mmsb2) - (mmst1*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1))/pow2(-mmgl +
     mmsb1) - (mmst1*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1))/pow2(-mmgl + mmsb1)
     - (mmst1*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmsb2))/pow2(-mmgl + mmsb2) - (
     mmst1*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2))/pow2(-mmgl + mmsb2) + (
     DeltaInv(mmt,mmsb1,mmst1)*pow3(mmst1))/(2.*(-mmgl + mmsb1)) + (DeltaInv(
     mmt,mmsb2,mmst1)*pow3(mmst1))/(2.*(-mmgl + mmsb2)) + (DeltaInv(mmt,mmst1,
     mmgl)*pow3(mmst1))/(-mmgl + mmsb1) + (DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1)
     )/(-mmgl + mmsb2) - (mmsb1*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmst1))/pow2(-
     mmgl + mmsb1) - (mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1))/pow2(-mmgl +
     mmsb1) - (mmsb2*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmst1))/pow2(-mmgl + mmsb2)
     - (mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1))/pow2(-mmgl + mmsb2)))/3. +
     (4*pow2(log(mmst1/mmu))*(0.25 - (3*mmst1)/(8.*(-mmgl + mmsb1)) - (3*mmst1)
     /(8.*(-mmgl + mmsb2)) - (mmsb1*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1))/(4.*(-
     mmgl + mmsb1)) - (mmsb2*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1))/(4.*(-mmgl +
     mmsb2)) - mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl) - mmsb1*mmst1*DeltaInv(mmt,
     mmst1,mmgl) - mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl) + mmst1*mmt*DeltaInv(
     mmt,mmst1,mmgl) - (mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl))/(-mmgl +
     mmsb1) - (mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl))/(-mmgl + mmsb2) + (
     mmst1*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1))/(4.*(-mmgl + mmsb1)) + (3*
     mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1))/(2.*(-mmgl + mmsb1)) + (mmsb1*
     mmst1)/(2.*pow2(-mmgl + mmsb1)) + (mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1)*
     pow2(mmsb1))/(2.*pow2(-mmgl + mmsb1)) + (mmst1*mmt*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmsb1))/(2.*pow2(-mmgl + mmsb1)) + (mmst1*DeltaInv(mmt,mmsb2,
     mmst1)*pow2(mmsb2))/(4.*(-mmgl + mmsb2)) + (3*mmst1*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmsb2))/(2.*(-mmgl + mmsb2)) + (mmsb2*mmst1)/(2.*pow2(-mmgl +
     mmsb2)) + (mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2))/(2.*pow2(-mmgl
     + mmsb2)) + (mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2))/(2.*pow2(-
     mmgl + mmsb2)) - (mmsb1*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1))/(2.*(-mmgl
     + mmsb1)) - (mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1))/(4.*(-mmgl +
     mmsb1)) - (mmsb2*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/(2.*(-mmgl +
     mmsb2)) - (mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/(4.*(-mmgl + mmsb2))
     + 2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1) - (2*mmsb1*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmst1))/(-mmgl + mmsb1) - (2*mmsb2*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmst1))/(-mmgl + mmsb2) - (mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/
     (2.*(-mmgl + mmsb1)) - (mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(2.*(-
     mmgl + mmsb2)) + (mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1))/(2.*
     pow2(-mmgl + mmsb1)) + (mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(
     2.*pow2(-mmgl + mmsb1)) + (DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow2(
     mmst1))/pow2(-mmgl + mmsb1) + (DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(
     mmst1))/pow2(-mmgl + mmsb1) + (mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(
     mmst1))/(2.*pow2(-mmgl + mmsb2)) + (mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmst1))/(2.*pow2(-mmgl + mmsb2)) + (DeltaInv(mmt,mmsb2,mmst1)*pow2(
     mmsb2)*pow2(mmst1))/pow2(-mmgl + mmsb2) + (DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmsb2)*pow2(mmst1))/pow2(-mmgl + mmsb2) - (mmst1*DeltaInv(mmt,mmsb1,mmst1)
     *pow3(mmsb1))/(2.*pow2(-mmgl + mmsb1)) - (mmst1*DeltaInv(mmt,mmst1,mmgl)*
     pow3(mmsb1))/(2.*pow2(-mmgl + mmsb1)) - (mmst1*DeltaInv(mmt,mmsb2,mmst1)*
     pow3(mmsb2))/(2.*pow2(-mmgl + mmsb2)) - (mmst1*DeltaInv(mmt,mmst1,mmgl)*
     pow3(mmsb2))/(2.*pow2(-mmgl + mmsb2)) + (DeltaInv(mmt,mmsb1,mmst1)*pow3(
     mmst1))/(4.*(-mmgl + mmsb1)) + (DeltaInv(mmt,mmsb2,mmst1)*pow3(mmst1))/(4.
     *(-mmgl + mmsb2)) + (DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1))/(2.*(-mmgl +
     mmsb1)) + (DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1))/(2.*(-mmgl + mmsb2)) - (
     mmsb1*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmst1))/(2.*pow2(-mmgl + mmsb1)) - (
     mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1))/(2.*pow2(-mmgl + mmsb1)) - (
     mmsb2*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmst1))/(2.*pow2(-mmgl + mmsb2)) - (
     mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1))/(2.*pow2(-mmgl + mmsb2))))/3.
      + (4*log(mmst1/mmu)*(
     0.9166666666666666 + (3*mmst1)/(2.*(-mmgl + mmsb1)) + (3*mmst1)/(2.*(-mmgl
     + mmsb2)) + (3*mmsb1*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1))/(2.*(-mmgl +
     mmsb1)) + (3*mmsb2*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1))/(2.*(-mmgl +
     mmsb2)) + 6*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl) + 6*mmsb1*mmst1*DeltaInv(
     mmt,mmst1,mmgl) + 6*mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl) - 6*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl) + (6*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl))/(-
     mmgl + mmsb1) + (6*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl))/(-mmgl +
     mmsb2) - (3*mmst1*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1))/(2.*(-mmgl +
     mmsb1)) - (9*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1))/(-mmgl + mmsb1) -
     (2*mmsb1*mmst1)/pow2(-mmgl + mmsb1) - (3*mmst1*mmt*DeltaInv(mmt,mmsb1,
     mmst1)*pow2(mmsb1))/pow2(-mmgl + mmsb1) - (3*mmst1*mmt*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmsb1))/pow2(-mmgl + mmsb1) - (3*mmst1*DeltaInv(mmt,mmsb2,
     mmst1)*pow2(mmsb2))/(2.*(-mmgl + mmsb2)) - (9*mmst1*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmsb2))/(-mmgl + mmsb2) - (2*mmsb2*mmst1)/pow2(-mmgl + mmsb2) -
     (3*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2))/pow2(-mmgl + mmsb2) -
     (3*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2))/pow2(-mmgl + mmsb2) + (
     3*mmsb1*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1))/(-mmgl + mmsb1) + (3*mmt*
     DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1))/(2.*(-mmgl + mmsb1)) + (3*mmsb2*
     DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/(-mmgl + mmsb2) + (3*mmt*DeltaInv(
     mmt,mmsb2,mmst1)*pow2(mmst1))/(2.*(-mmgl + mmsb2)) - 12*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmst1) + (12*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/
     (-mmgl + mmsb1) + (12*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(-mmgl +
     mmsb2) + (3*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(-mmgl + mmsb1) + (
     3*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(-mmgl + mmsb2) - (3*mmsb1*
     mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1))/pow2(-mmgl + mmsb1) - (3*mmsb1*
     mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/pow2(-mmgl + mmsb1) - (6*
     DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow2(mmst1))/pow2(-mmgl + mmsb1) - (
     6*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(mmst1))/pow2(-mmgl + mmsb1) -
     (3*mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/pow2(-mmgl + mmsb2) -
     (3*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/pow2(-mmgl + mmsb2) - (
     6*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2)*pow2(mmst1))/pow2(-mmgl + mmsb2) -
     (6*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(mmst1))/pow2(-mmgl + mmsb2) +
     (3*mmst1*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1))/pow2(-mmgl + mmsb1) + (3*
     mmst1*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1))/pow2(-mmgl + mmsb1) + (3*
     mmst1*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmsb2))/pow2(-mmgl + mmsb2) + (3*
     mmst1*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2))/pow2(-mmgl + mmsb2) - (3*
     DeltaInv(mmt,mmsb1,mmst1)*pow3(mmst1))/(2.*(-mmgl + mmsb1)) - (3*DeltaInv(
     mmt,mmsb2,mmst1)*pow3(mmst1))/(2.*(-mmgl + mmsb2)) - (3*DeltaInv(mmt,
     mmst1,mmgl)*pow3(mmst1))/(-mmgl + mmsb1) - (3*DeltaInv(mmt,mmst1,mmgl)*
     pow3(mmst1))/(-mmgl + mmsb2) + (3*mmsb1*DeltaInv(mmt,mmsb1,mmst1)*pow3(
     mmst1))/pow2(-mmgl + mmsb1) + (3*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow3(
     mmst1))/pow2(-mmgl + mmsb1) + (3*mmsb2*DeltaInv(mmt,mmsb2,mmst1)*pow3(
     mmst1))/pow2(-mmgl + mmsb2) + (3*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow3(
     mmst1))/pow2(-mmgl + mmsb2)))/3. + (14*DeltaInv(mmt,mmsb1,mmst2)*pow3(
     mmst2))/(3.*(-mmgl + mmsb1)) + (2*zt2*DeltaInv(mmt,mmsb1,mmst2)*pow3(
     mmst2))/(3.*(-mmgl + mmsb1)) + (14*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmst2))/
     (3.*(-mmgl + mmsb2)) + (2*zt2*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmst2))/(3.*(
     -mmgl + mmsb2)) + (28*DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2))/(3.*(-mmgl +
     mmsb1)) + (28*DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2))/(3.*(-mmgl + mmsb2)) +
     (4*zt2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2))/(3.*(-mmgl + mmsb1)) + (4*
     zt2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2))/(3.*(-mmgl + mmsb2)) - (28*
     mmsb1*DeltaInv(mmt,mmsb1,mmst2)*pow3(mmst2))/(3.*pow2(-mmgl + mmsb1)) - (
     4*mmsb1*zt2*DeltaInv(mmt,mmsb1,mmst2)*pow3(mmst2))/(3.*pow2(-mmgl + mmsb1)
     ) - (28*mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2))/(3.*pow2(-mmgl +
     mmsb1)) - (4*mmsb1*zt2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2))/(3.*pow2(-
     mmgl + mmsb1)) - (28*mmsb2*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmst2))/(3.*
     pow2(-mmgl + mmsb2)) - (4*mmsb2*zt2*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmst2))
     /(3.*pow2(-mmgl + mmsb2)) - (28*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2)
     )/(3.*pow2(-mmgl + mmsb2)) - (4*mmsb2*zt2*DeltaInv(mmt,mmst2,mmgl)*pow3(
     mmst2))/(3.*pow2(-mmgl + mmsb2)) + (4*log(mmst2/mmu)*log(mmt/mmu)*((-3*
     mmst2)/(2.*(-mmgl + mmsb1)) - (3*mmst2)/(2.*(-mmgl + mmsb2)) - (3*mmsb1*
     mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2))/(2.*(-mmgl + mmsb1)) - (3*mmsb2*
     mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2))/(2.*(-mmgl + mmsb2)) - 2*mmgl*mmst2*
     DeltaInv(mmt,mmst2,mmgl) - 2*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl) - 2*
     mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl) + 6*mmst2*mmt*DeltaInv(mmt,mmst2,
     mmgl) - (6*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl))/(-mmgl + mmsb1) - (6*
     mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl))/(-mmgl + mmsb2) + (mmst2*
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
     mmst2))/(2.*(-mmgl + mmsb1)) - (mmsb2*DeltaInv(mmt,mmsb2,mmst2)*pow2(
     mmst2))/(-mmgl + mmsb2) - (mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2))/(2.*
     (-mmgl + mmsb2)) + 4*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2) - (4*mmsb1*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(-mmgl + mmsb1) - (4*mmsb2*DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmst2))/(-mmgl + mmsb2) - (mmt*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmst2))/(-mmgl + mmsb1) - (mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmst2))/(-mmgl + mmsb2) + (mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2)
     )/pow2(-mmgl + mmsb1) + (mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/
     pow2(-mmgl + mmsb1) + (2*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*pow2(mmst2)
     )/pow2(-mmgl + mmsb1) + (2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(
     mmst2))/pow2(-mmgl + mmsb1) + (mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(
     mmst2))/pow2(-mmgl + mmsb2) + (mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmst2))/pow2(-mmgl + mmsb2) + (2*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*
     pow2(mmst2))/pow2(-mmgl + mmsb2) + (2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)
     *pow2(mmst2))/pow2(-mmgl + mmsb2) - (mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow3(
     mmsb1))/pow2(-mmgl + mmsb1) - (mmst2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1))
     /pow2(-mmgl + mmsb1) - (mmst2*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmsb2))/pow2(
     -mmgl + mmsb2) - (mmst2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2))/pow2(-mmgl +
     mmsb2) + (DeltaInv(mmt,mmsb1,mmst2)*pow3(mmst2))/(2.*(-mmgl + mmsb1)) + (
     DeltaInv(mmt,mmsb2,mmst2)*pow3(mmst2))/(2.*(-mmgl + mmsb2)) + (DeltaInv(
     mmt,mmst2,mmgl)*pow3(mmst2))/(-mmgl + mmsb1) + (DeltaInv(mmt,mmst2,mmgl)*
     pow3(mmst2))/(-mmgl + mmsb2) - (mmsb1*DeltaInv(mmt,mmsb1,mmst2)*pow3(
     mmst2))/pow2(-mmgl + mmsb1) - (mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2))
     /pow2(-mmgl + mmsb1) - (mmsb2*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmst2))/pow2(
     -mmgl + mmsb2) - (mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2))/pow2(-mmgl +
     mmsb2)))/3. + (4*pow2(log(mmst2/mmu))*(0.25 - (3*mmst2)/(8.*(-mmgl +
     mmsb1)) - (3*mmst2)/(8.*(-mmgl + mmsb2)) - (mmsb1*mmst2*mmt*DeltaInv(mmt,
     mmsb1,mmst2))/(4.*(-mmgl + mmsb1)) - (mmsb2*mmst2*mmt*DeltaInv(mmt,mmsb2,
     mmst2))/(4.*(-mmgl + mmsb2)) - mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl) -
     mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl) - mmsb2*mmst2*DeltaInv(mmt,mmst2,
     mmgl) + mmst2*mmt*DeltaInv(mmt,mmst2,mmgl) - (mmsb1*mmst2*mmt*DeltaInv(
     mmt,mmst2,mmgl))/(-mmgl + mmsb1) - (mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,
     mmgl))/(-mmgl + mmsb2) + (mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))/(4.
     *(-mmgl + mmsb1)) + (3*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/(2.*(-
     mmgl + mmsb1)) + (mmsb1*mmst2)/(2.*pow2(-mmgl + mmsb1)) + (mmst2*mmt*
     DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))/(2.*pow2(-mmgl + mmsb1)) + (mmst2*
     mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/(2.*pow2(-mmgl + mmsb1)) + (
     mmst2*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2))/(4.*(-mmgl + mmsb2)) + (3*
     mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/(2.*(-mmgl + mmsb2)) + (mmsb2*
     mmst2)/(2.*pow2(-mmgl + mmsb2)) + (mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2)*
     pow2(mmsb2))/(2.*pow2(-mmgl + mmsb2)) + (mmst2*mmt*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmsb2))/(2.*pow2(-mmgl + mmsb2)) - (mmsb1*DeltaInv(mmt,mmsb1,
     mmst2)*pow2(mmst2))/(2.*(-mmgl + mmsb1)) - (mmt*DeltaInv(mmt,mmsb1,mmst2)*
     pow2(mmst2))/(4.*(-mmgl + mmsb1)) - (mmsb2*DeltaInv(mmt,mmsb2,mmst2)*pow2(
     mmst2))/(2.*(-mmgl + mmsb2)) - (mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2))
     /(4.*(-mmgl + mmsb2)) + 2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2) - (2*mmsb1*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(-mmgl + mmsb1) - (2*mmsb2*DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmst2))/(-mmgl + mmsb2) - (mmt*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmst2))/(2.*(-mmgl + mmsb1)) - (mmt*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmst2))/(2.*(-mmgl + mmsb2)) + (mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst2)*
     pow2(mmst2))/(2.*pow2(-mmgl + mmsb1)) + (mmsb1*mmt*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmst2))/(2.*pow2(-mmgl + mmsb1)) + (DeltaInv(mmt,mmsb1,mmst2)*
     pow2(mmsb1)*pow2(mmst2))/pow2(-mmgl + mmsb1) + (DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmsb1)*pow2(mmst2))/pow2(-mmgl + mmsb1) + (mmsb2*mmt*DeltaInv(mmt,
     mmsb2,mmst2)*pow2(mmst2))/(2.*pow2(-mmgl + mmsb2)) + (mmsb2*mmt*DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmst2))/(2.*pow2(-mmgl + mmsb2)) + (DeltaInv(mmt,
     mmsb2,mmst2)*pow2(mmsb2)*pow2(mmst2))/pow2(-mmgl + mmsb2) + (DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmsb2)*pow2(mmst2))/pow2(-mmgl + mmsb2) - (mmst2*
     DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1))/(2.*pow2(-mmgl + mmsb1)) - (mmst2*
     DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1))/(2.*pow2(-mmgl + mmsb1)) - (mmst2*
     DeltaInv(mmt,mmsb2,mmst2)*pow3(mmsb2))/(2.*pow2(-mmgl + mmsb2)) - (mmst2*
     DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2))/(2.*pow2(-mmgl + mmsb2)) + (
     DeltaInv(mmt,mmsb1,mmst2)*pow3(mmst2))/(4.*(-mmgl + mmsb1)) + (DeltaInv(
     mmt,mmsb2,mmst2)*pow3(mmst2))/(4.*(-mmgl + mmsb2)) + (DeltaInv(mmt,mmst2,
     mmgl)*pow3(mmst2))/(2.*(-mmgl + mmsb1)) + (DeltaInv(mmt,mmst2,mmgl)*pow3(
     mmst2))/(2.*(-mmgl + mmsb2)) - (mmsb1*DeltaInv(mmt,mmsb1,mmst2)*pow3(
     mmst2))/(2.*pow2(-mmgl + mmsb1)) - (mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow3(
     mmst2))/(2.*pow2(-mmgl + mmsb1)) - (mmsb2*DeltaInv(mmt,mmsb2,mmst2)*pow3(
     mmst2))/(2.*pow2(-mmgl + mmsb2)) - (mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow3(
     mmst2))/(2.*pow2(-mmgl + mmsb2))))/3. + (4*log(mmst2/mmu)*(
     0.9166666666666666 + (3*mmst2)/(2.*(-mmgl + mmsb1)) + (3*mmst2)/(2.*(-mmgl
     + mmsb2)) + (3*mmsb1*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2))/(2.*(-mmgl +
     mmsb1)) + (3*mmsb2*mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2))/(2.*(-mmgl +
     mmsb2)) + 6*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl) + 6*mmsb1*mmst2*DeltaInv(
     mmt,mmst2,mmgl) + 6*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl) - 6*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl) + (6*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl))/(-
     mmgl + mmsb1) + (6*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl))/(-mmgl +
     mmsb2) - (3*mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))/(2.*(-mmgl +
     mmsb1)) - (9*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/(-mmgl + mmsb1) -
     (2*mmsb1*mmst2)/pow2(-mmgl + mmsb1) - (3*mmst2*mmt*DeltaInv(mmt,mmsb1,
     mmst2)*pow2(mmsb1))/pow2(-mmgl + mmsb1) - (3*mmst2*mmt*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmsb1))/pow2(-mmgl + mmsb1) - (3*mmst2*DeltaInv(mmt,mmsb2,
     mmst2)*pow2(mmsb2))/(2.*(-mmgl + mmsb2)) - (9*mmst2*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmsb2))/(-mmgl + mmsb2) - (2*mmsb2*mmst2)/pow2(-mmgl + mmsb2) -
     (3*mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2))/pow2(-mmgl + mmsb2) -
     (3*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/pow2(-mmgl + mmsb2) + (
     3*mmsb1*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2))/(-mmgl + mmsb1) + (3*mmt*
     DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2))/(2.*(-mmgl + mmsb1)) + (3*mmsb2*
     DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2))/(-mmgl + mmsb2) + (3*mmt*DeltaInv(
     mmt,mmsb2,mmst2)*pow2(mmst2))/(2.*(-mmgl + mmsb2)) - 12*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmst2) + (12*mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/
     (-mmgl + mmsb1) + (12*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(-mmgl +
     mmsb2) + (3*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(-mmgl + mmsb1) + (
     3*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(-mmgl + mmsb2) - (3*mmsb1*
     mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2))/pow2(-mmgl + mmsb1) - (3*mmsb1*
     mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/pow2(-mmgl + mmsb1) - (6*
     DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*pow2(mmst2))/pow2(-mmgl + mmsb1) - (
     6*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(mmst2))/pow2(-mmgl + mmsb1) -
     (3*mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2))/pow2(-mmgl + mmsb2) -
     (3*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/pow2(-mmgl + mmsb2) - (
     6*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow2(mmst2))/pow2(-mmgl + mmsb2) -
     (6*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(mmst2))/pow2(-mmgl + mmsb2) +
     (3*mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1))/pow2(-mmgl + mmsb1) + (3*
     mmst2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1))/pow2(-mmgl + mmsb1) + (3*
     mmst2*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmsb2))/pow2(-mmgl + mmsb2) + (3*
     mmst2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2))/pow2(-mmgl + mmsb2) - (3*
     DeltaInv(mmt,mmsb1,mmst2)*pow3(mmst2))/(2.*(-mmgl + mmsb1)) - (3*DeltaInv(
     mmt,mmsb2,mmst2)*pow3(mmst2))/(2.*(-mmgl + mmsb2)) - (3*DeltaInv(mmt,
     mmst2,mmgl)*pow3(mmst2))/(-mmgl + mmsb1) - (3*DeltaInv(mmt,mmst2,mmgl)*
     pow3(mmst2))/(-mmgl + mmsb2) + (3*mmsb1*DeltaInv(mmt,mmsb1,mmst2)*pow3(
     mmst2))/pow2(-mmgl + mmsb1) + (3*mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow3(
     mmst2))/pow2(-mmgl + mmsb1) + (3*mmsb2*DeltaInv(mmt,mmsb2,mmst2)*pow3(
     mmst2))/pow2(-mmgl + mmsb2) + (3*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow3(
     mmst2))/pow2(-mmgl + mmsb2)))/3. + Xb*((-72*mmgl*mmsb1)/(mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*mmsb2)/(mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) + (16*mmgl*mmsb1)/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (72*
     mmgl*mmsb2)/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (64*mmgl*mmst1)/(3.*
     mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (64*mmgl*mmst1)/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) + (64*mmgl*mmst2)/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) - (64*mmgl*mmst2)/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (128*
     mmgl*mmsusy)/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (128*mmgl*mmsusy)/(
     mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (32*mmgl*mmt)/(mgl*(-mmgl + mmsb1)*
     (mmsb1 - mmsb2)) - (32*mmgl*mmt)/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (
     80*mmgl*mmsb1*zt2)/(9.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (8*mmgl*
     mmsb2*zt2)/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (8*mmgl*mmsb1*zt2)/(
     3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (80*mmgl*mmsb2*zt2)/(9.*mgl*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (8*mmgl*mmst1*zt2)/(3.*mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) - (8*mmgl*mmst1*zt2)/(3.*mgl*(mmsb1 - mmsb2)*(-
     mmgl + mmsb2)) + (8*mmgl*mmst2*zt2)/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) - (8*mmgl*mmst2*zt2)/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (
     64*mmgl*mmsusy*zt2)/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (64*mmgl*
     mmsusy*zt2)/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*mmt*zt2)/(
     3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*mmt*zt2)/(3.*mgl*(mmsb1
     - mmsb2)*(-mmgl + mmsb2)) - (112*mmgl*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)
     )/(3.*mgl*(mmsb1 - mmsb2)) + (112*mmgl*mmsb2*mmst1*DeltaInv(mmt,mmst1,
     mmgl))/(3.*mgl*(mmsb1 - mmsb2)) - (112*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,
     mmgl))/(3.*mgl*(mmsb1 - mmsb2)) + (112*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst1,
     mmgl))/(3.*mgl*(mmsb1 - mmsb2)) + (448*mmgl*mmsb1*mmst1*mmt*DeltaInv(mmt,
     mmst1,mmgl))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (448*mmgl*mmsb2*
     mmst1*mmt*DeltaInv(mmt,mmst1,mmgl))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) - (16*mmgl*mmsb1*mmst1*zt2*DeltaInv(mmt,mmst1,mmgl))/(3.*mgl*(
     mmsb1 - mmsb2)) + (16*mmgl*mmsb2*mmst1*zt2*DeltaInv(mmt,mmst1,mmgl))/(3.*
     mgl*(mmsb1 - mmsb2)) - (16*mmgl*mmsb1*mmt*zt2*DeltaInv(mmt,mmst1,mmgl))/(
     3.*mgl*(mmsb1 - mmsb2)) + (16*mmgl*mmsb2*mmt*zt2*DeltaInv(mmt,mmst1,mmgl))
     /(3.*mgl*(mmsb1 - mmsb2)) + (64*mmgl*mmsb1*mmst1*mmt*zt2*DeltaInv(mmt,
     mmst1,mmgl))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (64*mmgl*mmsb2*
     mmst1*mmt*zt2*DeltaInv(mmt,mmst1,mmgl))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) - (112*mmgl*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl))/(3.*mgl*(mmsb1 -
     mmsb2)) + (112*mmgl*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl))/(3.*mgl*(mmsb1 -
     mmsb2)) - (112*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl))/(3.*mgl*(mmsb1 -
     mmsb2)) + (112*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl))/(3.*mgl*(mmsb1 -
     mmsb2)) + (448*mmgl*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl))/(3.*mgl*(-
     mmgl + mmsb1)*(mmsb1 - mmsb2)) - (448*mmgl*mmsb2*mmst2*mmt*DeltaInv(mmt,
     mmst2,mmgl))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*mmsb1*
     mmst2*zt2*DeltaInv(mmt,mmst2,mmgl))/(3.*mgl*(mmsb1 - mmsb2)) + (16*mmgl*
     mmsb2*mmst2*zt2*DeltaInv(mmt,mmst2,mmgl))/(3.*mgl*(mmsb1 - mmsb2)) - (16*
     mmgl*mmsb1*mmt*zt2*DeltaInv(mmt,mmst2,mmgl))/(3.*mgl*(mmsb1 - mmsb2)) + (
     16*mmgl*mmsb2*mmt*zt2*DeltaInv(mmt,mmst2,mmgl))/(3.*mgl*(mmsb1 - mmsb2)) +
     (64*mmgl*mmsb1*mmst2*mmt*zt2*DeltaInv(mmt,mmst2,mmgl))/(3.*mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) - (64*mmgl*mmsb2*mmst2*mmt*zt2*DeltaInv(mmt,mmst2,
     mmgl))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (64*mmgl*Fin20(mmgl,
     mmsusy,mmu))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (64*mmgl*Fin20(
     mmgl,mmsusy,mmu))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (176*mmgl*
     Fin20(mmsb1,mmgl,mmu))/(9.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (32*
     mmgl*Fin20(mmsb1,mmgl,mmu))/(9.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (8*
     mmgl*Fin20(mmsb1,mmsb2,mmu))/(9.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (
     8*mmgl*Fin20(mmsb1,mmsb2,mmu))/(9.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) -
     (32*mmgl*Fin20(mmsb2,mmgl,mmu))/(9.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) +
     (176*mmgl*Fin20(mmsb2,mmgl,mmu))/(9.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2))
     + (8*mmgl*Fin3(mmt,mmst1,mmgl,mmu))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) - (8*mmgl*Fin3(mmt,mmst1,mmgl,mmu))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl
     + mmsb2)) + (16*mmgl*mmsb1*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,
     mmu))/(3.*mgl*(mmsb1 - mmsb2)) - (16*mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl)*
     Fin3(mmt,mmst1,mmgl,mmu))/(3.*mgl*(mmsb1 - mmsb2)) + (32*mmgl*mmsb1*mmst1*
     DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu))/(3.*mgl*(-mmgl + mmsb1)
     *(mmsb1 - mmsb2)) - (32*mmgl*mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)*Fin3(
     mmt,mmst1,mmgl,mmu))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*
     mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu))/(3.*mgl*(-
     mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst1,
     mmgl)*Fin3(mmt,mmst1,mmgl,mmu))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) +
     (16*mmgl*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu))/(3.*
     mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*mmst1*mmt*DeltaInv(mmt,
     mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) + (8*mmgl*Fin3(mmt,mmst2,mmgl,mmu))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1
     - mmsb2)) - (8*mmgl*Fin3(mmt,mmst2,mmgl,mmu))/(3.*mgl*(mmsb1 - mmsb2)*(-
     mmgl + mmsb2)) + (16*mmgl*mmsb1*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,
     mmgl,mmu))/(3.*mgl*(mmsb1 - mmsb2)) - (16*mmgl*mmsb2*DeltaInv(mmt,mmst2,
     mmgl)*Fin3(mmt,mmst2,mmgl,mmu))/(3.*mgl*(mmsb1 - mmsb2)) + (32*mmgl*mmsb1*
     mmst2*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu))/(3.*mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) - (32*mmgl*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl)*
     Fin3(mmt,mmst2,mmgl,mmu))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*
     mmgl*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu))/(3.*mgl*
     (-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,
     mmgl)*Fin3(mmt,mmst2,mmgl,mmu))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) +
     (16*mmgl*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu))/(3.*
     mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*mmst2*mmt*DeltaInv(mmt,
     mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) + (96*mmgl*mmsb1*log(mmgl/mmu))/(mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) + (16*mmgl*mmsb2*log(mmgl/mmu))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) - (16*mmgl*mmsb1*log(mmgl/mmu))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) - (96*mmgl*mmsb2*log(mmgl/mmu))/(mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) + (16*mmgl*mmst1*log(mmgl/mmu))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) - (16*mmgl*mmst1*log(mmgl/mmu))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) + (16*mmgl*mmst2*log(mmgl/mmu))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) - (16*mmgl*mmst2*log(mmgl/mmu))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) + (128*mmgl*mmsusy*log(mmgl/mmu))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) - (128*mmgl*mmsusy*log(mmgl/mmu))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) - (32*mmgl*mmt*log(mmgl/mmu))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) + (32*mmgl*mmt*log(mmgl/mmu))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) + (32*mmgl*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu))/(
     mgl*(mmsb1 - mmsb2)) - (32*mmgl*mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)*log(
     mmgl/mmu))/(mgl*(mmsb1 - mmsb2)) + (16*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,
     mmgl)*log(mmgl/mmu))/(mgl*(mmsb1 - mmsb2)) - (16*mmgl*mmsb2*mmt*DeltaInv(
     mmt,mmst1,mmgl)*log(mmgl/mmu))/(mgl*(mmsb1 - mmsb2)) - (16*mmgl*mmsb1*
     mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu))/(mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) + (16*mmgl*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(
     mmgl/mmu))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (32*mmgl*mmsb1*mmst2*
     DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu))/(mgl*(mmsb1 - mmsb2)) - (32*mmgl*
     mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu))/(mgl*(mmsb1 - mmsb2))
     + (16*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu))/(mgl*(mmsb1 -
     mmsb2)) - (16*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu))/(mgl*
     (mmsb1 - mmsb2)) - (16*mmgl*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(
     mmgl/mmu))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*mmgl*mmsb2*mmst2*
     mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu))/(mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) + (128*mmgl*log(mmsb1/mmu))/(9.*mgl*(mmsb1 - mmsb2)) - (112*mmgl*
     mmsb1*log(mmsb1/mmu))/(9.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*mmgl*
     mmsb2*log(mmsb1/mmu))/(9.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (32*mmgl*
     mmsb1*log(mmsb1/mmu))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*
     mmsb2*log(mmsb1/mmu))/(9.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (128*
     mmgl*log(mmgl/mmu)*log(mmsb1/mmu))/(9.*mgl*(mmsb1 - mmsb2)) - (200*mmgl*
     mmsb1*log(mmgl/mmu)*log(mmsb1/mmu))/(9.*mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) - (80*mmgl*mmsb2*log(mmgl/mmu)*log(mmsb1/mmu))/(9.*mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) + (8*mmgl*mmsb1*log(mmgl/mmu)*log(mmsb1/mmu))/(3.*
     mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (80*mmgl*mmsb2*log(mmgl/mmu)*log(
     mmsb1/mmu))/(9.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (128*mmgl*log(
     mmsb2/mmu))/(9.*mgl*(mmsb1 - mmsb2)) + (80*mmgl*mmsb2*log(mmsb2/mmu))/(9.*
     mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (128*mmgl*mmsb2*log(mmsb2/mmu))/(9.
     *mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (128*mmgl*log(mmgl/mmu)*log(mmsb2/
     mmu))/(9.*mgl*(mmsb1 - mmsb2)) + (8*mmgl*mmsb2*log(mmgl/mmu)*log(mmsb2/
     mmu))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (152*mmgl*mmsb2*log(mmgl/
     mmu)*log(mmsb2/mmu))/(9.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*
     mmsb2*log(mmsb1/mmu)*log(mmsb2/mmu))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) - (80*mmgl*mmsb2*log(mmsb1/mmu)*log(mmsb2/mmu))/(9.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*mmst1*log(mmst1/mmu))/(3.*mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) + (16*mmgl*mmst1*log(mmst1/mmu))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*
     log(mmst1/mmu))/(mgl*(mmsb1 - mmsb2)) + (16*mmgl*mmsb2*mmst1*DeltaInv(mmt,
     mmst1,mmgl)*log(mmst1/mmu))/(mgl*(mmsb1 - mmsb2)) - (16*mmgl*mmsb1*mmst1*
     mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmst1/mmu))/(mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) + (16*mmgl*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmst1/mmu)
     )/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (8*mmgl*mmst1*log(mmgl/mmu)*log(
     mmst1/mmu))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (8*mmgl*mmst1*log(
     mmgl/mmu)*log(mmst1/mmu))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (32*
     mmgl*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmst1/mmu)
     )/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (32*mmgl*mmsb2*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmst1/mmu))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*mmst2*log(mmst2/mmu))/(3.*mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) + (16*mmgl*mmst2*log(mmst2/mmu))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*
     log(mmst2/mmu))/(mgl*(mmsb1 - mmsb2)) + (16*mmgl*mmsb2*mmst2*DeltaInv(mmt,
     mmst2,mmgl)*log(mmst2/mmu))/(mgl*(mmsb1 - mmsb2)) - (16*mmgl*mmsb1*mmst2*
     mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmst2/mmu))/(mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) + (16*mmgl*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmst2/mmu)
     )/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (8*mmgl*mmst2*log(mmgl/mmu)*log(
     mmst2/mmu))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (8*mmgl*mmst2*log(
     mmgl/mmu)*log(mmst2/mmu))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (32*
     mmgl*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmst2/mmu)
     )/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (32*mmgl*mmsb2*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmst2/mmu))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) + (256*mmgl*mmsusy*log(mmsusy/mmu))/(3.*mgl*(-mmgl
     + mmsb1)*(mmsb1 - mmsb2)) - (256*mmgl*mmsusy*log(mmsusy/mmu))/(3.*mgl*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (64*mmgl*mmsusy*log(mmgl/mmu)*log(
     mmsusy/mmu))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (64*mmgl*mmsusy*
     log(mmgl/mmu)*log(mmsusy/mmu))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) -
     (32*mmgl*mmsb1*log(mmt/mmu))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (32*
     mmgl*mmsb2*log(mmt/mmu))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*
     mmst1*log(mmt/mmu))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*mmgl*
     mmst1*log(mmt/mmu))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*
     mmst2*log(mmt/mmu))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*mmgl*
     mmst2*log(mmt/mmu))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (64*mmgl*mmt*
     log(mmt/mmu))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (64*mmgl*mmt*log(
     mmt/mmu))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*mmsb1*mmst1*
     DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu))/(mgl*(mmsb1 - mmsb2)) - (16*mmgl*
     mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu))/(mgl*(mmsb1 - mmsb2)) +
     (16*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu))/(mgl*(mmsb1 -
     mmsb2)) - (16*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu))/(mgl*(
     mmsb1 - mmsb2)) - (96*mmgl*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(
     mmt/mmu))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (96*mmgl*mmsb2*mmst1*
     mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu))/(mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) + (16*mmgl*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu))/(
     mgl*(mmsb1 - mmsb2)) - (16*mmgl*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(
     mmt/mmu))/(mgl*(mmsb1 - mmsb2)) + (16*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst2,
     mmgl)*log(mmt/mmu))/(mgl*(mmsb1 - mmsb2)) - (16*mmgl*mmsb2*mmt*DeltaInv(
     mmt,mmst2,mmgl)*log(mmt/mmu))/(mgl*(mmsb1 - mmsb2)) - (96*mmgl*mmsb1*
     mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu))/(mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) + (96*mmgl*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(
     mmt/mmu))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (32*mmgl*mmsb1*log(mmgl/
     mmu)*log(mmt/mmu))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (32*mmgl*
     mmsb2*log(mmgl/mmu)*log(mmt/mmu))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2))
     + (16*mmgl*mmt*log(mmgl/mmu)*log(mmt/mmu))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1
     - mmsb2)) - (16*mmgl*mmt*log(mmgl/mmu)*log(mmt/mmu))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) - (32*mmgl*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*
     log(mmgl/mmu)*log(mmt/mmu))/(3.*mgl*(mmsb1 - mmsb2)) + (32*mmgl*mmsb2*
     mmst1*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu))/(3.*mgl*(mmsb1
     - mmsb2)) - (16*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(
     mmt/mmu))/(3.*mgl*(mmsb1 - mmsb2)) + (16*mmgl*mmsb2*mmt*DeltaInv(mmt,
     mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu))/(3.*mgl*(mmsb1 - mmsb2)) + (16*
     mmgl*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu))/
     (mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*mmsb2*mmst1*mmt*DeltaInv(
     mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu))/(mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) - (32*mmgl*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(
     mmt/mmu))/(3.*mgl*(mmsb1 - mmsb2)) + (32*mmgl*mmsb2*mmst2*DeltaInv(mmt,
     mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu))/(3.*mgl*(mmsb1 - mmsb2)) - (16*
     mmgl*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu))/(3.*
     mgl*(mmsb1 - mmsb2)) + (16*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*log(
     mmgl/mmu)*log(mmt/mmu))/(3.*mgl*(mmsb1 - mmsb2)) + (16*mmgl*mmsb1*mmst2*
     mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu))/(mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,
     mmgl)*log(mmgl/mmu)*log(mmt/mmu))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) +
     (16*mmgl*mmst1*log(mmst1/mmu)*log(mmt/mmu))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1
     - mmsb2)) - (16*mmgl*mmst1*log(mmst1/mmu)*log(mmt/mmu))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*
     log(mmst1/mmu)*log(mmt/mmu))/(3.*mgl*(mmsb1 - mmsb2)) - (16*mmgl*mmsb2*
     mmst1*DeltaInv(mmt,mmst1,mmgl)*log(mmst1/mmu)*log(mmt/mmu))/(3.*mgl*(mmsb1
     - mmsb2)) + (16*mmgl*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmst1/
     mmu)*log(mmt/mmu))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*mmsb2*
     mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmst1/mmu)*log(mmt/mmu))/(mgl*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*mmst2*log(mmst2/mmu)*log(mmt/
     mmu))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*mmst2*log(mmst2/
     mmu)*log(mmt/mmu))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*
     mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(mmst2/mmu)*log(mmt/mmu))/(3.*mgl*
     (mmsb1 - mmsb2)) - (16*mmgl*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(
     mmst2/mmu)*log(mmt/mmu))/(3.*mgl*(mmsb1 - mmsb2)) + (16*mmgl*mmsb1*mmst2*
     mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmst2/mmu)*log(mmt/mmu))/(mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,
     mmgl)*log(mmst2/mmu)*log(mmt/mmu))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) +
     (112*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmgl))/(3.*mgl*(mmsb1 - mmsb2)) -
     (112*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmgl))/(3.*mgl*(mmsb1 - mmsb2)) +
     (16*mmsb1*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmgl))/(3.*mgl*(mmsb1 - mmsb2)
     ) - (16*mmsb2*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmgl))/(3.*mgl*(mmsb1 -
     mmsb2)) + (112*mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(mmgl))/(3.*mgl*(mmsb1 -
     mmsb2)) - (112*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmgl))/(3.*mgl*(mmsb1 -
     mmsb2)) + (16*mmsb1*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmgl))/(3.*mgl*(
     mmsb1 - mmsb2)) - (16*mmsb2*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmgl))/(3.*
     mgl*(mmsb1 - mmsb2)) - (16*mmsb1*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*
     pow2(mmgl))/(mgl*(mmsb1 - mmsb2)) + (16*mmsb2*DeltaInv(mmt,mmst1,mmgl)*
     log(mmgl/mmu)*pow2(mmgl))/(mgl*(mmsb1 - mmsb2)) - (16*mmsb1*DeltaInv(mmt,
     mmst2,mmgl)*log(mmgl/mmu)*pow2(mmgl))/(mgl*(mmsb1 - mmsb2)) + (16*mmsb2*
     DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*pow2(mmgl))/(mgl*(mmsb1 - mmsb2)) -
     (16*mmsb1*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*pow2(mmgl))/(mgl*(mmsb1 -
     mmsb2)) + (16*mmsb2*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*pow2(mmgl))/(
     mgl*(mmsb1 - mmsb2)) - (16*mmsb1*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*
     pow2(mmgl))/(mgl*(mmsb1 - mmsb2)) + (16*mmsb2*DeltaInv(mmt,mmst2,mmgl)*
     log(mmt/mmu)*pow2(mmgl))/(mgl*(mmsb1 - mmsb2)) + (16*mmsb1*DeltaInv(mmt,
     mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmgl))/(3.*mgl*(mmsb1 - mmsb2)
     ) - (16*mmsb2*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(
     mmgl))/(3.*mgl*(mmsb1 - mmsb2)) + (16*mmsb1*DeltaInv(mmt,mmst2,mmgl)*log(
     mmgl/mmu)*log(mmt/mmu)*pow2(mmgl))/(3.*mgl*(mmsb1 - mmsb2)) - (16*mmsb2*
     DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmgl))/(3.*mgl*(
     mmsb1 - mmsb2)) + (112*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1))/(3.*mgl*
     (mmsb1 - mmsb2)) + (112*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1))/(
     3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (112*mmgl*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmsb1))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*
     mmgl*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1))/(3.*mgl*(mmsb1 - mmsb2)) +
     (16*mmgl*mmst1*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1))/(3.*mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) + (16*mmgl*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmsb1))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (112*mmgl*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmsb1))/(3.*mgl*(mmsb1 - mmsb2)) + (112*mmgl*mmst2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) + (112*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/(3.*mgl*(-
     mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*mmgl*zt2*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmsb1))/(3.*mgl*(mmsb1 - mmsb2)) + (16*mmgl*mmst2*zt2*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmsb1))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*
     mmgl*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/(3.*mgl*(-mmgl + mmsb1)
     *(mmsb1 - mmsb2)) - (16*mmgl*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,
     mmu)*pow2(mmsb1))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*
     DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*pow2(mmsb1))/(3.*mgl*(-
     mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*DeltaInv(mmt,mmst1,mmgl)*log(
     mmgl/mmu)*pow2(mmsb1))/(mgl*(mmsb1 - mmsb2)) - (32*mmgl*mmst1*DeltaInv(
     mmt,mmst1,mmgl)*log(mmgl/mmu)*pow2(mmsb1))/(mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) - (16*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*pow2(mmsb1))
     /(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*DeltaInv(mmt,mmst2,mmgl)
     *log(mmgl/mmu)*pow2(mmsb1))/(mgl*(mmsb1 - mmsb2)) - (32*mmgl*mmst2*
     DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*pow2(mmsb1))/(mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) - (16*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*
     pow2(mmsb1))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*mmgl*mmst1*
     DeltaInv(mmt,mmst1,mmgl)*log(mmst1/mmu)*pow2(mmsb1))/(mgl*(-mmgl + mmsb1)*
     (mmsb1 - mmsb2)) + (16*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(mmst2/mmu)*
     pow2(mmsb1))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*DeltaInv(
     mmt,mmst1,mmgl)*log(mmt/mmu)*pow2(mmsb1))/(mgl*(mmsb1 - mmsb2)) - (16*
     mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*pow2(mmsb1))/(mgl*(-mmgl
     + mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmt/
     mmu)*pow2(mmsb1))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*
     DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*pow2(mmsb1))/(mgl*(mmsb1 - mmsb2)) -
     (16*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*pow2(mmsb1))/(mgl*(-
     mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*
     log(mmt/mmu)*pow2(mmsb1))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*
     mmgl*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmsb1))/(3.*
     mgl*(mmsb1 - mmsb2)) + (32*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/
     mmu)*log(mmt/mmu)*pow2(mmsb1))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) +
     (16*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(
     mmsb1))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*mmgl*DeltaInv(mmt,
     mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmsb1))/(3.*mgl*(mmsb1 -
     mmsb2)) + (32*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmt/
     mmu)*pow2(mmsb1))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*mmgl*mmt*
     DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmsb1))/(3.*mgl*(
     -mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*
     log(mmst1/mmu)*log(mmt/mmu)*pow2(mmsb1))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) - (16*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(mmst2/mmu)*log(mmt/
     mmu)*pow2(mmsb1))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (64*mmgl*
     mmsb1*Fin20(mmgl,mmsusy,mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1))
     + (64*mmgl*mmsusy*Fin20(mmgl,mmsusy,mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-
     mmgl + mmsb1)) + (8*mmgl*mmsb1*Fin20(mmsb1,mmsb2,mmu))/(3.*mgl*(mmsb1 -
     mmsb2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmsb2*Fin20(mmsb1,mmsb2,mmu))/(3.*
     mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) + (64*mmgl*mmsb1*Fin20(mmsb1,
     mmsusy,mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (64*mmgl*
     mmsusy*Fin20(mmsb1,mmsusy,mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl +
     mmsb1)) - (8*mmgl*mmsb1*Fin20(mmsb2,mmgl,mmu))/(3.*mgl*(mmsb1 - mmsb2)*
     pow2(-mmgl + mmsb1)) + (8*mmgl*mmsb2*Fin20(mmsb2,mmgl,mmu))/(3.*mgl*(mmsb1
     - mmsb2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmsb1*Fin3(mmt,mmsb1,mmst1,mmu))/(
     3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmst1*Fin3(mmt,
     mmsb1,mmst1,mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) + (8*mmgl*
     mmt*Fin3(mmt,mmsb1,mmst1,mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)
     ) + (8*mmgl*mmsb1*Fin3(mmt,mmsb1,mmst2,mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(
     -mmgl + mmsb1)) - (8*mmgl*mmst2*Fin3(mmt,mmsb1,mmst2,mmu))/(3.*mgl*(mmsb1
     - mmsb2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmt*Fin3(mmt,mmsb1,mmst2,mmu))/(3.
     *mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmsb1*Fin3(mmt,mmst1,
     mmgl,mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmst1*
     Fin3(mmt,mmst1,mmgl,mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (
     8*mmgl*mmt*Fin3(mmt,mmst1,mmgl,mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl +
     mmsb1)) - (8*mmgl*mmsb1*Fin3(mmt,mmst2,mmgl,mmu))/(3.*mgl*(mmsb1 - mmsb2)*
     pow2(-mmgl + mmsb1)) + (8*mmgl*mmst2*Fin3(mmt,mmst2,mmgl,mmu))/(3.*mgl*(
     mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmt*Fin3(mmt,mmst2,mmgl,mmu)
     )/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (16*mmgl*mmsb1*mmsb2*log(
     mmgl/mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (16*mmgl*mmsb1*
     mmst1*log(mmgl/mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (16*
     mmgl*mmsb1*mmst2*log(mmgl/mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl +
     mmsb1)) - (128*mmgl*mmsb1*mmsusy*log(mmgl/mmu))/(3.*mgl*(mmsb1 - mmsb2)*
     pow2(-mmgl + mmsb1)) + (32*mmgl*mmsb1*mmt*log(mmgl/mmu))/(3.*mgl*(mmsb1 -
     mmsb2)*pow2(-mmgl + mmsb1)) + (16*mmgl*mmsb1*mmsb2*log(mmsb1/mmu))/(3.*
     mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) + (16*mmgl*mmsb1*mmst1*log(mmsb1/
     mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) + (16*mmgl*mmsb1*mmst2*
     log(mmsb1/mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) + (128*mmgl*
     mmsb1*mmsusy*log(mmsb1/mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1))
     - (32*mmgl*mmsb1*mmt*log(mmsb1/mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl +
     mmsb1)) - (16*mmgl*mmsb1*mmsb2*log(mmgl/mmu)*log(mmsb1/mmu))/(9.*mgl*(
     mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmsb1*mmsb2*log(mmgl/mmu)*
     log(mmsb2/mmu))/(9.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (8*mmgl*
     mmsb1*mmsb2*log(mmsb1/mmu)*log(mmsb2/mmu))/(9.*mgl*(mmsb1 - mmsb2)*pow2(-
     mmgl + mmsb1)) + (8*mmgl*mmsb1*mmst1*log(mmgl/mmu)*log(mmst1/mmu))/(3.*
     mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmsb1*mmst1*log(mmsb1/
     mmu)*log(mmst1/mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) + (8*
     mmgl*mmsb1*mmst2*log(mmgl/mmu)*log(mmst2/mmu))/(3.*mgl*(mmsb1 - mmsb2)*
     pow2(-mmgl + mmsb1)) - (8*mmgl*mmsb1*mmst2*log(mmsb1/mmu)*log(mmst2/mmu))/
     (3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) + (64*mmgl*mmsb1*mmsusy*log(
     mmgl/mmu)*log(mmsusy/mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) -
     (64*mmgl*mmsb1*mmsusy*log(mmsb1/mmu)*log(mmsusy/mmu))/(3.*mgl*(mmsb1 -
     mmsb2)*pow2(-mmgl + mmsb1)) - (16*mmgl*mmsb1*mmt*log(mmgl/mmu)*log(mmt/
     mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) + (16*mmgl*mmsb1*mmt*
     log(mmsb1/mmu)*log(mmt/mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1))
     + (80*mmgl*log(mmgl/mmu)*pow2(mmsb1))/(mgl*(mmsb1 - mmsb2)*pow2(-mmgl +
     mmsb1)) - (80*mmgl*log(mmsb1/mmu)*pow2(mmsb1))/(mgl*(mmsb1 - mmsb2)*pow2(-
     mmgl + mmsb1)) + (296*mmgl*log(mmgl/mmu)*log(mmsb1/mmu)*pow2(mmsb1))/(9.*
     mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (112*mmgl*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmsb2))/(3.*mgl*(mmsb1 - mmsb2)) - (112*mmgl*mmst1*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmsb2))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (
     112*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2))/(3.*mgl*(mmsb1 - mmsb2)
     *(-mmgl + mmsb2)) - (16*mmgl*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2))/(3.
     *mgl*(mmsb1 - mmsb2)) - (16*mmgl*mmst1*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmsb2))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*mmt*zt2*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) - (112*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/(3.*mgl*(mmsb1 -
     mmsb2)) - (112*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/(3.*mgl*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (112*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmsb2))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*zt2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/(3.*mgl*(mmsb1 - mmsb2)) - (16*mmgl*
     mmst2*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/(3.*mgl*(mmsb1 - mmsb2)*(-
     mmgl + mmsb2)) - (16*mmgl*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/(
     3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*DeltaInv(mmt,mmst1,
     mmgl)*Fin3(mmt,mmst1,mmgl,mmu)*pow2(mmsb2))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl
     + mmsb2)) + (16*mmgl*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*
     pow2(mmsb2))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*DeltaInv(
     mmt,mmst1,mmgl)*log(mmgl/mmu)*pow2(mmsb2))/(mgl*(mmsb1 - mmsb2)) + (32*
     mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*pow2(mmsb2))/(mgl*(mmsb1
     - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*log(
     mmgl/mmu)*pow2(mmsb2))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*
     DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*pow2(mmsb2))/(mgl*(mmsb1 - mmsb2))
     + (32*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*pow2(mmsb2))/(mgl*
     (mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*
     log(mmgl/mmu)*pow2(mmsb2))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*
     mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*log(mmst1/mmu)*pow2(mmsb2))/(mgl*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*
     log(mmst2/mmu)*pow2(mmsb2))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*
     mmgl*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*pow2(mmsb2))/(mgl*(mmsb1 -
     mmsb2)) + (16*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*pow2(mmsb2)
     )/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*mmt*DeltaInv(mmt,mmst1,
     mmgl)*log(mmt/mmu)*pow2(mmsb2))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (
     16*mmgl*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*pow2(mmsb2))/(mgl*(mmsb1 -
     mmsb2)) + (16*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*pow2(mmsb2)
     )/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*mmt*DeltaInv(mmt,mmst2,
     mmgl)*log(mmt/mmu)*pow2(mmsb2))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (
     16*mmgl*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmsb2))/(
     3.*mgl*(mmsb1 - mmsb2)) - (32*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*log(
     mmgl/mmu)*log(mmt/mmu)*pow2(mmsb2))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) - (16*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu)
     *pow2(mmsb2))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*
     DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmsb2))/(3.*mgl*(
     mmsb1 - mmsb2)) - (32*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*
     log(mmt/mmu)*pow2(mmsb2))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*
     mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmsb2))/
     (3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*mmst1*DeltaInv(mmt,
     mmst1,mmgl)*log(mmst1/mmu)*log(mmt/mmu)*pow2(mmsb2))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*log(
     mmst2/mmu)*log(mmt/mmu)*pow2(mmsb2))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) - (16*mmgl*log(mmgl/mmu)*log(mmsb1/mmu)*pow2(mmsb2))/(9.*mgl*(
     mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (16*mmgl*log(mmgl/mmu)*log(mmsb2/
     mmu)*pow2(mmsb2))/(9.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) + (16*mmgl*
     log(mmsb1/mmu)*log(mmsb2/mmu)*pow2(mmsb2))/(9.*mgl*(mmsb1 - mmsb2)*pow2(-
     mmgl + mmsb1)) + (16*mmgl*log(mmsb1/mmu)*pow2(mmsb2))/(9.*mgl*(-mmgl +
     mmsb1)*pow2(mmsb1 - mmsb2)) - (16*mmgl*log(mmsb1/mmu)*pow2(mmsb2))/(9.*
     mgl*(-mmgl + mmsb2)*pow2(mmsb1 - mmsb2)) - (112*mmgl*log(mmgl/mmu)*log(
     mmsb1/mmu)*pow2(mmsb2))/(9.*mgl*(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) + (
     112*mmgl*log(mmgl/mmu)*log(mmsb1/mmu)*pow2(mmsb2))/(9.*mgl*(-mmgl + mmsb2)
     *pow2(mmsb1 - mmsb2)) - (16*mmgl*log(mmsb2/mmu)*pow2(mmsb2))/(9.*mgl*(-
     mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) + (16*mmgl*log(mmsb2/mmu)*pow2(mmsb2))/
     (9.*mgl*(-mmgl + mmsb2)*pow2(mmsb1 - mmsb2)) + (16*mmgl*log(mmgl/mmu)*log(
     mmsb2/mmu)*pow2(mmsb2))/(9.*mgl*(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (
     16*mmgl*log(mmgl/mmu)*log(mmsb2/mmu)*pow2(mmsb2))/(9.*mgl*(-mmgl + mmsb2)*
     pow2(mmsb1 - mmsb2)) - (16*mmgl*log(mmsb1/mmu)*log(mmsb2/mmu)*pow2(mmsb2))
     /(9.*mgl*(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (112*mmgl*log(mmsb1/mmu)*
     log(mmsb2/mmu)*pow2(mmsb2))/(9.*mgl*(-mmgl + mmsb2)*pow2(mmsb1 - mmsb2)) +
     (64*mmgl*mmsb2*Fin20(mmgl,mmsusy,mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl
     + mmsb2)) - (64*mmgl*mmsusy*Fin20(mmgl,mmsusy,mmu))/(3.*mgl*(mmsb1 -
     mmsb2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmsb1*Fin20(mmsb1,mmgl,mmu))/(3.*
     mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmsb2*Fin20(mmsb1,mmgl,
     mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmsb1*Fin20(
     mmsb1,mmsb2,mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (8*mmgl*
     mmsb2*Fin20(mmsb1,mmsb2,mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2))
     - (64*mmgl*mmsb2*Fin20(mmsb2,mmsusy,mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-
     mmgl + mmsb2)) + (64*mmgl*mmsusy*Fin20(mmsb2,mmsusy,mmu))/(3.*mgl*(mmsb1 -
     mmsb2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmsb2*Fin3(mmt,mmsb2,mmst1,mmu))/(3.
     *mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmst1*Fin3(mmt,mmsb2,
     mmst1,mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmt*
     Fin3(mmt,mmsb2,mmst1,mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) -
     (8*mmgl*mmsb2*Fin3(mmt,mmsb2,mmst2,mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-
     mmgl + mmsb2)) + (8*mmgl*mmst2*Fin3(mmt,mmsb2,mmst2,mmu))/(3.*mgl*(mmsb1 -
     mmsb2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmt*Fin3(mmt,mmsb2,mmst2,mmu))/(3.*
     mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmsb2*Fin3(mmt,mmst1,
     mmgl,mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmst1*
     Fin3(mmt,mmst1,mmgl,mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (
     8*mmgl*mmt*Fin3(mmt,mmst1,mmgl,mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl +
     mmsb2)) + (8*mmgl*mmsb2*Fin3(mmt,mmst2,mmgl,mmu))/(3.*mgl*(mmsb1 - mmsb2)*
     pow2(-mmgl + mmsb2)) - (8*mmgl*mmst2*Fin3(mmt,mmst2,mmgl,mmu))/(3.*mgl*(
     mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmt*Fin3(mmt,mmst2,mmgl,mmu)
     )/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (16*mmgl*mmsb1*mmsb2*log(
     mmgl/mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (16*mmgl*mmsb2*
     mmst1*log(mmgl/mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (16*
     mmgl*mmsb2*mmst2*log(mmgl/mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl +
     mmsb2)) + (128*mmgl*mmsb2*mmsusy*log(mmgl/mmu))/(3.*mgl*(mmsb1 - mmsb2)*
     pow2(-mmgl + mmsb2)) - (32*mmgl*mmsb2*mmt*log(mmgl/mmu))/(3.*mgl*(mmsb1 -
     mmsb2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmsb1*mmsb2*log(mmgl/mmu)*log(mmsb1/
     mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (16*mmgl*mmsb1*mmsb2*
     log(mmsb2/mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (16*mmgl*
     mmsb2*mmst1*log(mmsb2/mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) -
     (16*mmgl*mmsb2*mmst2*log(mmsb2/mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl +
     mmsb2)) - (128*mmgl*mmsb2*mmsusy*log(mmsb2/mmu))/(3.*mgl*(mmsb1 - mmsb2)*
     pow2(-mmgl + mmsb2)) + (32*mmgl*mmsb2*mmt*log(mmsb2/mmu))/(3.*mgl*(mmsb1 -
     mmsb2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmsb1*mmsb2*log(mmsb1/mmu)*log(
     mmsb2/mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmsb2*
     mmst1*log(mmgl/mmu)*log(mmst1/mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl +
     mmsb2)) + (8*mmgl*mmsb2*mmst1*log(mmsb2/mmu)*log(mmst1/mmu))/(3.*mgl*(
     mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmsb2*mmst2*log(mmgl/mmu)*
     log(mmst2/mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (8*mmgl*
     mmsb2*mmst2*log(mmsb2/mmu)*log(mmst2/mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-
     mmgl + mmsb2)) - (64*mmgl*mmsb2*mmsusy*log(mmgl/mmu)*log(mmsusy/mmu))/(3.*
     mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (64*mmgl*mmsb2*mmsusy*log(
     mmsb2/mmu)*log(mmsusy/mmu))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) +
     (16*mmgl*mmsb2*mmt*log(mmgl/mmu)*log(mmt/mmu))/(3.*mgl*(mmsb1 - mmsb2)*
     pow2(-mmgl + mmsb2)) - (16*mmgl*mmsb2*mmt*log(mmsb2/mmu)*log(mmt/mmu))/(3.
     *mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (80*mmgl*log(mmgl/mmu)*pow2(
     mmsb2))/(mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (16*mmgl*log(mmgl/mmu)
     *log(mmsb1/mmu)*pow2(mmsb2))/(9.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2))
     + (80*mmgl*log(mmsb2/mmu)*pow2(mmsb2))/(mgl*(mmsb1 - mmsb2)*pow2(-mmgl +
     mmsb2)) - (104*mmgl*log(mmgl/mmu)*log(mmsb2/mmu)*pow2(mmsb2))/(3.*mgl*(
     mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (16*mmgl*log(mmsb1/mmu)*log(mmsb2/
     mmu)*pow2(mmsb2))/(9.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (112*
     mmgl*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(3.*mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) - (112*mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(
     3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (112*mmgl*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmst1))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (112*
     mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(3.*mgl*(mmsb1 - mmsb2)*(-
     mmgl + mmsb2)) + (16*mmgl*mmsb1*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/
     (3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*mmsb2*zt2*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmst1))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (
     16*mmgl*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(3.*mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*mmt*zt2*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmst1))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*DeltaInv(mmt,
     mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,mmu)*pow2(mmst1))/(3.*mgl*(-mmgl + mmsb1)*
     (mmsb1 - mmsb2)) + (16*mmgl*DeltaInv(mmt,mmst1,mmgl)*Fin3(mmt,mmst1,mmgl,
     mmu)*pow2(mmst1))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*
     mmsb1*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*pow2(mmst1))/(mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl)*log(
     mmgl/mmu)*pow2(mmst1))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (32*mmgl*
     mmsb1*DeltaInv(mmt,mmst1,mmgl)*log(mmst1/mmu)*pow2(mmst1))/(mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) + (32*mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl)*log(
     mmst1/mmu)*pow2(mmst1))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*
     mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmst1/mmu)*pow2(mmst1))/(mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) + (16*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmst1/
     mmu)*pow2(mmst1))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*mmsb1*
     DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*pow2(mmst1))/(mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) + (16*mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*
     pow2(mmst1))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*mmt*
     DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*pow2(mmst1))/(mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) + (16*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*pow2(
     mmst1))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*mmsb1*DeltaInv(
     mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmst1))/(3.*mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) + (16*mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl)*log(
     mmgl/mmu)*log(mmt/mmu)*pow2(mmst1))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) + (32*mmgl*mmsb1*DeltaInv(mmt,mmst1,mmgl)*log(mmst1/mmu)*log(mmt/
     mmu)*pow2(mmst1))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (32*mmgl*
     mmsb2*DeltaInv(mmt,mmst1,mmgl)*log(mmst1/mmu)*log(mmt/mmu)*pow2(mmst1))/(
     3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*mmt*DeltaInv(mmt,mmst1,
     mmgl)*log(mmst1/mmu)*log(mmt/mmu)*pow2(mmst1))/(3.*mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) - (16*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*log(mmst1/mmu)*
     log(mmt/mmu)*pow2(mmst1))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (112*
     mmgl*mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(3.*mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) - (112*mmgl*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(
     3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (112*mmgl*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmst2))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (112*
     mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(3.*mgl*(mmsb1 - mmsb2)*(-
     mmgl + mmsb2)) + (16*mmgl*mmsb1*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/
     (3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*mmsb2*zt2*DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmst2))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (
     16*mmgl*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(3.*mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*mmt*zt2*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmst2))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*DeltaInv(mmt,
     mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,mmu)*pow2(mmst2))/(3.*mgl*(-mmgl + mmsb1)*
     (mmsb1 - mmsb2)) + (16*mmgl*DeltaInv(mmt,mmst2,mmgl)*Fin3(mmt,mmst2,mmgl,
     mmu)*pow2(mmst2))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*
     mmsb1*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*pow2(mmst2))/(mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*mmsb2*DeltaInv(mmt,mmst2,mmgl)*log(
     mmgl/mmu)*pow2(mmst2))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (32*mmgl*
     mmsb1*DeltaInv(mmt,mmst2,mmgl)*log(mmst2/mmu)*pow2(mmst2))/(mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) + (32*mmgl*mmsb2*DeltaInv(mmt,mmst2,mmgl)*log(
     mmst2/mmu)*pow2(mmst2))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*
     mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmst2/mmu)*pow2(mmst2))/(mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) + (16*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmst2/
     mmu)*pow2(mmst2))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*mmsb1*
     DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*pow2(mmst2))/(mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) + (16*mmgl*mmsb2*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*
     pow2(mmst2))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*mmt*
     DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*pow2(mmst2))/(mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) + (16*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*pow2(
     mmst2))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*mmsb1*DeltaInv(
     mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow2(mmst2))/(3.*mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) + (16*mmgl*mmsb2*DeltaInv(mmt,mmst2,mmgl)*log(
     mmgl/mmu)*log(mmt/mmu)*pow2(mmst2))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) + (32*mmgl*mmsb1*DeltaInv(mmt,mmst2,mmgl)*log(mmst2/mmu)*log(mmt/
     mmu)*pow2(mmst2))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (32*mmgl*
     mmsb2*DeltaInv(mmt,mmst2,mmgl)*log(mmst2/mmu)*log(mmt/mmu)*pow2(mmst2))/(
     3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*mmt*DeltaInv(mmt,mmst2,
     mmgl)*log(mmst2/mmu)*log(mmt/mmu)*pow2(mmst2))/(3.*mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) - (16*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*log(mmst2/mmu)*
     log(mmt/mmu)*pow2(mmst2))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (148*
     mmgl*mmsb1*pow2(log(mmgl/mmu)))/(9.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) +
     (4*mmgl*mmsb2*pow2(log(mmgl/mmu)))/(9.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)
     ) + (4*mmgl*mmsb1*pow2(log(mmgl/mmu)))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) + (44*mmgl*mmsb2*pow2(log(mmgl/mmu)))/(3.*mgl*(mmsb1 - mmsb2)*(-
     mmgl + mmsb2)) - (4*mmgl*mmst1*pow2(log(mmgl/mmu)))/(3.*mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) + (4*mmgl*mmst1*pow2(log(mmgl/mmu)))/(3.*mgl*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (4*mmgl*mmst2*pow2(log(mmgl/mmu)))/(3.*
     mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (4*mmgl*mmst2*pow2(log(mmgl/mmu)))/
     (3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (32*mmgl*mmsusy*pow2(log(mmgl/
     mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (32*mmgl*mmsusy*pow2(
     log(mmgl/mmu)))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (8*mmgl*mmt*
     pow2(log(mmgl/mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (8*mmgl*
     mmt*pow2(log(mmgl/mmu)))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*
     mmgl*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmgl/mmu)))/(3.*mgl*(
     mmsb1 - mmsb2)) + (16*mmgl*mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(log(
     mmgl/mmu)))/(3.*mgl*(mmsb1 - mmsb2)) - (8*mmgl*mmsb1*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow2(log(mmgl/mmu)))/(3.*mgl*(mmsb1 - mmsb2)) + (8*mmgl*mmsb2*
     mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmgl/mmu)))/(3.*mgl*(mmsb1 - mmsb2))
     + (8*mmgl*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmgl/mmu)))/(
     3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (8*mmgl*mmsb2*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmgl/mmu)))/(3.*mgl*(mmsb1 - mmsb2)*(-
     mmgl + mmsb2)) - (16*mmgl*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(log(
     mmgl/mmu)))/(3.*mgl*(mmsb1 - mmsb2)) + (16*mmgl*mmsb2*mmst2*DeltaInv(mmt,
     mmst2,mmgl)*pow2(log(mmgl/mmu)))/(3.*mgl*(mmsb1 - mmsb2)) - (8*mmgl*mmsb1*
     mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmgl/mmu)))/(3.*mgl*(mmsb1 - mmsb2))
     + (8*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmgl/mmu)))/(3.*mgl*
     (mmsb1 - mmsb2)) + (8*mmgl*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(
     log(mmgl/mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (8*mmgl*mmsb2*
     mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmgl/mmu)))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) + (8*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmgl)*
     pow2(log(mmgl/mmu)))/(3.*mgl*(mmsb1 - mmsb2)) - (8*mmsb2*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmgl)*pow2(log(mmgl/mmu)))/(3.*mgl*(mmsb1 - mmsb2)) + (8*
     mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(mmgl)*pow2(log(mmgl/mmu)))/(3.*mgl*(
     mmsb1 - mmsb2)) - (8*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmgl)*pow2(log(
     mmgl/mmu)))/(3.*mgl*(mmsb1 - mmsb2)) + (8*mmgl*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmsb1)*pow2(log(mmgl/mmu)))/(3.*mgl*(mmsb1 - mmsb2)) + (16*mmgl*
     mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(log(mmgl/mmu)))/(3.*mgl*(-
     mmgl + mmsb1)*(mmsb1 - mmsb2)) + (8*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmsb1)*pow2(log(mmgl/mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2))
     + (8*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(log(mmgl/mmu)))/(3.*
     mgl*(mmsb1 - mmsb2)) + (16*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)
     *pow2(log(mmgl/mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (8*mmgl*
     mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(log(mmgl/mmu)))/(3.*mgl*(-
     mmgl + mmsb1)*(mmsb1 - mmsb2)) + (28*mmgl*mmsb1*mmsb2*pow2(log(mmgl/mmu)))
     /(9.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) + (4*mmgl*mmsb1*mmst1*pow2(
     log(mmgl/mmu)))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) + (4*mmgl*
     mmsb1*mmst2*pow2(log(mmgl/mmu)))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl +
     mmsb1)) + (32*mmgl*mmsb1*mmsusy*pow2(log(mmgl/mmu)))/(3.*mgl*(mmsb1 -
     mmsb2)*pow2(-mmgl + mmsb1)) - (8*mmgl*mmsb1*mmt*pow2(log(mmgl/mmu)))/(3.*
     mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (140*mmgl*pow2(mmsb1)*pow2(log(
     mmgl/mmu)))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (8*mmgl*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(log(mmgl/mmu)))/(3.*mgl*(mmsb1 -
     mmsb2)) - (16*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(log(
     mmgl/mmu)))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (8*mmgl*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(log(mmgl/mmu)))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) - (8*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*
     pow2(log(mmgl/mmu)))/(3.*mgl*(mmsb1 - mmsb2)) - (16*mmgl*mmst2*DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(log(mmgl/mmu)))/(3.*mgl*(mmsb1 - mmsb2)*(
     -mmgl + mmsb2)) - (8*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(
     log(mmgl/mmu)))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*pow2(
     mmsb2)*pow2(log(mmgl/mmu)))/(9.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) +
     (16*mmgl*pow2(mmsb2)*pow2(log(mmgl/mmu)))/(3.*mgl*(-mmgl + mmsb1)*pow2(
     mmsb1 - mmsb2)) - (16*mmgl*pow2(mmsb2)*pow2(log(mmgl/mmu)))/(3.*mgl*(-mmgl
     + mmsb2)*pow2(mmsb1 - mmsb2)) - (4*mmgl*mmsb1*mmsb2*pow2(log(mmgl/mmu)))/(
     3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (4*mmgl*mmsb2*mmst1*pow2(
     log(mmgl/mmu)))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (4*mmgl*
     mmsb2*mmst2*pow2(log(mmgl/mmu)))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl +
     mmsb2)) - (32*mmgl*mmsb2*mmsusy*pow2(log(mmgl/mmu)))/(3.*mgl*(mmsb1 -
     mmsb2)*pow2(-mmgl + mmsb2)) + (8*mmgl*mmsb2*mmt*pow2(log(mmgl/mmu)))/(3.*
     mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (436*mmgl*pow2(mmsb2)*pow2(log(
     mmgl/mmu)))/(9.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmsb1*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow2(log(mmgl/mmu)))/(3.*mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) + (8*mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmst1)*pow2(log(mmgl/mmu)))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (8*
     mmgl*mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow2(log(mmgl/mmu)))/(3.*
     mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (8*mmgl*mmsb2*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmst2)*pow2(log(mmgl/mmu)))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) + (44*mmgl*mmsb1*pow2(log(mmsb1/mmu)))/(9.*mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) + (64*mmgl*mmsb2*pow2(log(mmsb1/mmu)))/(9.*mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) + (4*mmgl*mmsb1*pow2(log(mmsb1/mmu)))/(3.*mgl*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (4*mmgl*mmsb1*mmsb2*pow2(log(mmsb1/mmu))
     )/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (4*mmgl*mmsb1*mmst1*pow2(
     log(mmsb1/mmu)))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) - (4*mmgl*
     mmsb1*mmst2*pow2(log(mmsb1/mmu)))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl +
     mmsb1)) - (32*mmgl*mmsb1*mmsusy*pow2(log(mmsb1/mmu)))/(3.*mgl*(mmsb1 -
     mmsb2)*pow2(-mmgl + mmsb1)) + (8*mmgl*mmsb1*mmt*pow2(log(mmsb1/mmu)))/(3.*
     mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) + (124*mmgl*pow2(mmsb1)*pow2(log(
     mmsb1/mmu)))/(9.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb1)) + (64*mmgl*pow2(
     mmsb2)*pow2(log(mmsb1/mmu)))/(9.*mgl*(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2))
     - (4*mmgl*mmsb2*pow2(log(mmsb2/mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) + (20*mmgl*mmsb2*pow2(log(mmsb2/mmu)))/(9.*mgl*(mmsb1 - mmsb2)*(-
     mmgl + mmsb2)) + (64*mmgl*pow2(mmsb2)*pow2(log(mmsb2/mmu)))/(9.*mgl*(-mmgl
     + mmsb2)*pow2(mmsb1 - mmsb2)) + (4*mmgl*mmsb1*mmsb2*pow2(log(mmsb2/mmu)))/
     (3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (4*mmgl*mmsb2*mmst1*pow2(
     log(mmsb2/mmu)))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (4*mmgl*
     mmsb2*mmst2*pow2(log(mmsb2/mmu)))/(3.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl +
     mmsb2)) + (32*mmgl*mmsb2*mmsusy*pow2(log(mmsb2/mmu)))/(3.*mgl*(mmsb1 -
     mmsb2)*pow2(-mmgl + mmsb2)) - (8*mmgl*mmsb2*mmt*pow2(log(mmsb2/mmu)))/(3.*
     mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (124*mmgl*pow2(mmsb2)*pow2(log(
     mmsb2/mmu)))/(9.*mgl*(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (4*mmgl*mmst1*
     pow2(log(mmst1/mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (4*mmgl*
     mmst1*pow2(log(mmst1/mmu)))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (8*
     mmgl*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmst1/mmu)))/(3.*mgl*(
     mmsb1 - mmsb2)) - (8*mmgl*mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(log(
     mmst1/mmu)))/(3.*mgl*(mmsb1 - mmsb2)) + (8*mmgl*mmsb1*mmst1*mmt*DeltaInv(
     mmt,mmst1,mmgl)*pow2(log(mmst1/mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) - (8*mmgl*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmst1/
     mmu)))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (8*mmgl*mmst1*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(log(mmst1/mmu)))/(3.*mgl*(-mmgl + mmsb1)*
     (mmsb1 - mmsb2)) + (8*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*
     pow2(log(mmst1/mmu)))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*
     mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow2(log(mmst1/mmu)))/(3.*mgl*(
     -mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmst1)*pow2(log(mmst1/mmu)))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2))
     + (8*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow2(log(mmst1/mmu)))/(
     3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (8*mmgl*mmt*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmst1)*pow2(log(mmst1/mmu)))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) + (4*mmgl*mmst2*pow2(log(mmst2/mmu)))/(3.*mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) - (4*mmgl*mmst2*pow2(log(mmst2/mmu)))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) + (8*mmgl*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl)*
     pow2(log(mmst2/mmu)))/(3.*mgl*(mmsb1 - mmsb2)) - (8*mmgl*mmsb2*mmst2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmst2/mmu)))/(3.*mgl*(mmsb1 - mmsb2)) +
     (8*mmgl*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmst2/mmu)))/(3.
     *mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (8*mmgl*mmsb2*mmst2*mmt*DeltaInv(
     mmt,mmst2,mmgl)*pow2(log(mmst2/mmu)))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) - (8*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(log(
     mmst2/mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (8*mmgl*mmst2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(log(mmst2/mmu)))/(3.*mgl*(mmsb1
     - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmst2)*pow2(log(mmst2/mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (
     16*mmgl*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow2(log(mmst2/mmu)))/(
     3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (8*mmgl*mmt*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmst2)*pow2(log(mmst2/mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) - (8*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow2(log(mmst2/
     mmu)))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (32*mmgl*mmsusy*pow2(
     log(mmsusy/mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (32*mmgl*
     mmsusy*pow2(log(mmsusy/mmu)))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (
     16*mmgl*mmsb1*pow2(log(mmt/mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2))
     - (16*mmgl*mmsb2*pow2(log(mmt/mmu)))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) + (8*mmgl*mmst1*pow2(log(mmt/mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1
     - mmsb2)) - (8*mmgl*mmst1*pow2(log(mmt/mmu)))/(3.*mgl*(mmsb1 - mmsb2)*(-
     mmgl + mmsb2)) + (8*mmgl*mmst2*pow2(log(mmt/mmu)))/(3.*mgl*(-mmgl + mmsb1)
     *(mmsb1 - mmsb2)) - (8*mmgl*mmst2*pow2(log(mmt/mmu)))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) + (8*mmgl*mmt*pow2(log(mmt/mmu)))/(3.*mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) - (8*mmgl*mmt*pow2(log(mmt/mmu)))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) - (8*mmgl*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl)*
     pow2(log(mmt/mmu)))/(3.*mgl*(mmsb1 - mmsb2)) + (8*mmgl*mmsb2*mmst1*
     DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmt/mmu)))/(3.*mgl*(mmsb1 - mmsb2)) - (
     8*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmt/mmu)))/(3.*mgl*(
     mmsb1 - mmsb2)) + (8*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmt/
     mmu)))/(3.*mgl*(mmsb1 - mmsb2)) + (16*mmgl*mmsb1*mmst1*mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow2(log(mmt/mmu)))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (
     16*mmgl*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmt/mmu)))/(mgl*
     (mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (8*mmgl*mmsb1*mmst2*DeltaInv(mmt,mmst2,
     mmgl)*pow2(log(mmt/mmu)))/(3.*mgl*(mmsb1 - mmsb2)) + (8*mmgl*mmsb2*mmst2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmt/mmu)))/(3.*mgl*(mmsb1 - mmsb2)) - (
     8*mmgl*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmt/mmu)))/(3.*mgl*(
     mmsb1 - mmsb2)) + (8*mmgl*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmt/
     mmu)))/(3.*mgl*(mmsb1 - mmsb2)) + (16*mmgl*mmsb1*mmst2*mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow2(log(mmt/mmu)))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (
     16*mmgl*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmt/mmu)))/(mgl*
     (mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (8*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmgl)*pow2(log(mmt/mmu)))/(3.*mgl*(mmsb1 - mmsb2)) - (8*mmsb2*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmgl)*pow2(log(mmt/mmu)))/(3.*mgl*(mmsb1 - mmsb2)) +
     (8*mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(mmgl)*pow2(log(mmt/mmu)))/(3.*mgl*(
     mmsb1 - mmsb2)) - (8*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmgl)*pow2(log(
     mmt/mmu)))/(3.*mgl*(mmsb1 - mmsb2)) + (8*mmgl*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmsb1)*pow2(log(mmt/mmu)))/(3.*mgl*(mmsb1 - mmsb2)) + (8*mmgl*mmst1*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(log(mmt/mmu)))/(3.*mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) + (8*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)
     *pow2(log(mmt/mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (8*mmgl*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(log(mmt/mmu)))/(3.*mgl*(mmsb1 -
     mmsb2)) + (8*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(log(mmt/
     mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (8*mmgl*mmt*DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(log(mmt/mmu)))/(3.*mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) - (8*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(log(
     mmt/mmu)))/(3.*mgl*(mmsb1 - mmsb2)) - (8*mmgl*mmst1*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmsb2)*pow2(log(mmt/mmu)))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) - (8*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(log(mmt/
     mmu)))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (8*mmgl*DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmsb2)*pow2(log(mmt/mmu)))/(3.*mgl*(mmsb1 - mmsb2)) - (8*
     mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(log(mmt/mmu)))/(3.*
     mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (8*mmgl*mmt*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmsb2)*pow2(log(mmt/mmu)))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) + (8*mmgl*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow2(log(mmt/
     mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (8*mmgl*mmsb2*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmst1)*pow2(log(mmt/mmu)))/(3.*mgl*(mmsb1 - mmsb2)*(-
     mmgl + mmsb2)) + (8*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow2(
     log(mmt/mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (8*mmgl*mmt*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1)*pow2(log(mmt/mmu)))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) + (8*mmgl*mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmst2)*pow2(log(mmt/mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (8*
     mmgl*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow2(log(mmt/mmu)))/(3.*
     mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (8*mmgl*mmt*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmst2)*pow2(log(mmt/mmu)))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) - (8*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)*pow2(log(mmt/
     mmu)))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (112*mmgl*DeltaInv(mmt,
     mmst1,mmgl)*pow3(mmsb1))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*
     mmgl*zt2*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1))/(3.*mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) - (112*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1))/(3.*mgl*
     (-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*zt2*DeltaInv(mmt,mmst2,mmgl)*
     pow3(mmsb1))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*mmgl*DeltaInv(
     mmt,mmst1,mmgl)*log(mmgl/mmu)*pow3(mmsb1))/(mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) + (16*mmgl*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*pow3(mmsb1))/(
     mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*mmgl*DeltaInv(mmt,mmst1,mmgl)*
     log(mmt/mmu)*pow3(mmsb1))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*
     mmgl*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*pow3(mmsb1))/(mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*
     log(mmt/mmu)*pow3(mmsb1))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (16*
     mmgl*DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow3(mmsb1))/(3.*
     mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (8*mmgl*DeltaInv(mmt,mmst1,mmgl)*
     pow2(log(mmgl/mmu))*pow3(mmsb1))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2))
     - (8*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmgl/mmu))*pow3(mmsb1))/(3.*
     mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) - (8*mmgl*DeltaInv(mmt,mmst1,mmgl)*
     pow2(log(mmt/mmu))*pow3(mmsb1))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) -
     (8*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmt/mmu))*pow3(mmsb1))/(3.*mgl*(
     -mmgl + mmsb1)*(mmsb1 - mmsb2)) - (32*mmgl*log(mmgl/mmu)*log(mmsb1/mmu)*
     pow3(mmsb1))/(9.*mgl*(mmsb1 - mmsb2)*pow3(-mmgl + mmsb1)) + (16*mmgl*pow2(
     log(mmgl/mmu))*pow3(mmsb1))/(9.*mgl*(mmsb1 - mmsb2)*pow3(-mmgl + mmsb1)) +
     (16*mmgl*pow2(log(mmsb1/mmu))*pow3(mmsb1))/(9.*mgl*(mmsb1 - mmsb2)*pow3(-
     mmgl + mmsb1)) + (112*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2))/(3.*mgl*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*zt2*DeltaInv(mmt,mmst1,mmgl)*
     pow3(mmsb2))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (112*mmgl*
     DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl +
     mmsb2)) + (16*mmgl*zt2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2))/(3.*mgl*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*DeltaInv(mmt,mmst1,mmgl)*log(
     mmgl/mmu)*pow3(mmsb2))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*
     DeltaInv(mmt,mmst2,mmgl)*log(mmgl/mmu)*pow3(mmsb2))/(mgl*(mmsb1 - mmsb2)*(
     -mmgl + mmsb2)) - (16*mmgl*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*pow3(
     mmsb2))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*DeltaInv(mmt,
     mmst2,mmgl)*log(mmt/mmu)*pow3(mmsb2))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)
     ) + (16*mmgl*DeltaInv(mmt,mmst1,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow3(
     mmsb2))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*DeltaInv(mmt,
     mmst2,mmgl)*log(mmgl/mmu)*log(mmt/mmu)*pow3(mmsb2))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*log(mmgl/mmu)*log(mmsb1/mmu)*pow3(
     mmsb2))/(9.*mgl*pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (16*mmgl*log(
     mmgl/mmu)*log(mmsb2/mmu)*pow3(mmsb2))/(9.*mgl*pow2(-mmgl + mmsb1)*pow2(
     mmsb1 - mmsb2)) + (16*mmgl*log(mmsb1/mmu)*log(mmsb2/mmu)*pow3(mmsb2))/(9.*
     mgl*pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (16*mmgl*log(mmgl/mmu)*log(
     mmsb1/mmu)*pow3(mmsb2))/(9.*mgl*pow2(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) -
     (16*mmgl*log(mmgl/mmu)*log(mmsb2/mmu)*pow3(mmsb2))/(9.*mgl*pow2(mmsb1 -
     mmsb2)*pow2(-mmgl + mmsb2)) + (16*mmgl*log(mmsb1/mmu)*log(mmsb2/mmu)*pow3(
     mmsb2))/(9.*mgl*pow2(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (8*mmgl*
     DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmgl/mmu))*pow3(mmsb2))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) + (8*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmgl/
     mmu))*pow3(mmsb2))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*
     pow2(log(mmgl/mmu))*pow3(mmsb2))/(9.*mgl*pow2(-mmgl + mmsb1)*pow2(mmsb1 -
     mmsb2)) + (16*mmgl*pow2(log(mmgl/mmu))*pow3(mmsb2))/(9.*mgl*pow2(mmsb1 -
     mmsb2)*pow2(-mmgl + mmsb2)) + (8*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(log(
     mmt/mmu))*pow3(mmsb2))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (8*mmgl*
     DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmt/mmu))*pow3(mmsb2))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) - (32*mmgl*log(mmgl/mmu)*log(mmsb1/mmu)*pow3(
     mmsb2))/(9.*mgl*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) + (32*mmgl*log(mmgl/
     mmu)*log(mmsb1/mmu)*pow3(mmsb2))/(9.*mgl*(-mmgl + mmsb2)*pow3(mmsb1 -
     mmsb2)) - (32*mmgl*log(mmgl/mmu)*log(mmsb2/mmu)*pow3(mmsb2))/(9.*mgl*(-
     mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) + (32*mmgl*log(mmgl/mmu)*log(mmsb2/mmu)
     *pow3(mmsb2))/(9.*mgl*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) + (32*mmgl*log(
     mmsb1/mmu)*log(mmsb2/mmu)*pow3(mmsb2))/(9.*mgl*(-mmgl + mmsb1)*pow3(mmsb1
     - mmsb2)) - (32*mmgl*log(mmsb1/mmu)*log(mmsb2/mmu)*pow3(mmsb2))/(9.*mgl*(-
     mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) + (32*mmgl*pow2(log(mmgl/mmu))*pow3(
     mmsb2))/(9.*mgl*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) - (32*mmgl*pow2(log(
     mmgl/mmu))*pow3(mmsb2))/(9.*mgl*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) + (
     32*mmgl*log(mmgl/mmu)*log(mmsb2/mmu)*pow3(mmsb2))/(9.*mgl*(mmsb1 - mmsb2)*
     pow3(-mmgl + mmsb2)) - (16*mmgl*pow2(log(mmgl/mmu))*pow3(mmsb2))/(9.*mgl*(
     mmsb1 - mmsb2)*pow3(-mmgl + mmsb2)) - (16*mmgl*pow2(log(mmsb2/mmu))*pow3(
     mmsb2))/(9.*mgl*(mmsb1 - mmsb2)*pow3(-mmgl + mmsb2)) - (112*mmgl*DeltaInv(
     mmt,mmst1,mmgl)*pow3(mmst1))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (
     112*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1))/(3.*mgl*(mmsb1 - mmsb2)*(-
     mmgl + mmsb2)) - (16*mmgl*zt2*DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1))/(3.*
     mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*mmgl*zt2*DeltaInv(mmt,mmst1,
     mmgl)*pow3(mmst1))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*
     DeltaInv(mmt,mmst1,mmgl)*log(mmst1/mmu)*pow3(mmst1))/(mgl*(-mmgl + mmsb1)*
     (mmsb1 - mmsb2)) - (16*mmgl*DeltaInv(mmt,mmst1,mmgl)*log(mmst1/mmu)*pow3(
     mmst1))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*DeltaInv(mmt,
     mmst1,mmgl)*log(mmt/mmu)*pow3(mmst1))/(mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)
     ) - (16*mmgl*DeltaInv(mmt,mmst1,mmgl)*log(mmt/mmu)*pow3(mmst1))/(mgl*(
     mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*DeltaInv(mmt,mmst1,mmgl)*log(
     mmst1/mmu)*log(mmt/mmu)*pow3(mmst1))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) + (16*mmgl*DeltaInv(mmt,mmst1,mmgl)*log(mmst1/mmu)*log(mmt/mmu)*
     pow3(mmst1))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (8*mmgl*DeltaInv(
     mmt,mmst1,mmgl)*pow2(log(mmst1/mmu))*pow3(mmst1))/(3.*mgl*(-mmgl + mmsb1)*
     (mmsb1 - mmsb2)) + (8*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmst1/mmu))*
     pow3(mmst1))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (8*mmgl*DeltaInv(
     mmt,mmst1,mmgl)*pow2(log(mmt/mmu))*pow3(mmst1))/(3.*mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) + (8*mmgl*DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmt/mmu))*
     pow3(mmst1))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (112*mmgl*
     DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) + (112*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2))/(3.*mgl*(mmsb1 -
     mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*zt2*DeltaInv(mmt,mmst2,mmgl)*pow3(
     mmst2))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 - mmsb2)) + (16*mmgl*zt2*DeltaInv(
     mmt,mmst2,mmgl)*pow3(mmst2))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (
     16*mmgl*DeltaInv(mmt,mmst2,mmgl)*log(mmst2/mmu)*pow3(mmst2))/(mgl*(-mmgl +
     mmsb1)*(mmsb1 - mmsb2)) - (16*mmgl*DeltaInv(mmt,mmst2,mmgl)*log(mmst2/mmu)
     *pow3(mmst2))/(mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) + (16*mmgl*DeltaInv(
     mmt,mmst2,mmgl)*log(mmt/mmu)*pow3(mmst2))/(mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) - (16*mmgl*DeltaInv(mmt,mmst2,mmgl)*log(mmt/mmu)*pow3(mmst2))/(
     mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (16*mmgl*DeltaInv(mmt,mmst2,mmgl)*
     log(mmst2/mmu)*log(mmt/mmu)*pow3(mmst2))/(3.*mgl*(-mmgl + mmsb1)*(mmsb1 -
     mmsb2)) + (16*mmgl*DeltaInv(mmt,mmst2,mmgl)*log(mmst2/mmu)*log(mmt/mmu)*
     pow3(mmst2))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (8*mmgl*DeltaInv(
     mmt,mmst2,mmgl)*pow2(log(mmst2/mmu))*pow3(mmst2))/(3.*mgl*(-mmgl + mmsb1)*
     (mmsb1 - mmsb2)) + (8*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmst2/mmu))*
     pow3(mmst2))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2)) - (8*mmgl*DeltaInv(
     mmt,mmst2,mmgl)*pow2(log(mmt/mmu))*pow3(mmst2))/(3.*mgl*(-mmgl + mmsb1)*(
     mmsb1 - mmsb2)) + (8*mmgl*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmt/mmu))*
     pow3(mmst2))/(3.*mgl*(mmsb1 - mmsb2)*(-mmgl + mmsb2))) - (28*DeltaInv(mmt,
     mmsb1,mmst1)*pow4(mmsb1))/(3.*pow2(-mmgl + mmsb1)) - (4*zt2*DeltaInv(mmt,
     mmsb1,mmst1)*pow4(mmsb1))/(3.*pow2(-mmgl + mmsb1)) - (28*DeltaInv(mmt,
     mmsb1,mmst2)*pow4(mmsb1))/(3.*pow2(-mmgl + mmsb1)) - (4*zt2*DeltaInv(mmt,
     mmsb1,mmst2)*pow4(mmsb1))/(3.*pow2(-mmgl + mmsb1)) - (28*DeltaInv(mmt,
     mmst1,mmgl)*pow4(mmsb1))/(3.*pow2(-mmgl + mmsb1)) - (4*zt2*DeltaInv(mmt,
     mmst1,mmgl)*pow4(mmsb1))/(3.*pow2(-mmgl + mmsb1)) - (28*DeltaInv(mmt,
     mmst2,mmgl)*pow4(mmsb1))/(3.*pow2(-mmgl + mmsb1)) - (4*zt2*DeltaInv(mmt,
     mmst2,mmgl)*pow4(mmsb1))/(3.*pow2(-mmgl + mmsb1)) + (4*DeltaInv(mmt,mmst1,
     mmgl)*log(mmgl/mmu)*pow4(mmsb1))/pow2(-mmgl + mmsb1) + (4*DeltaInv(mmt,
     mmst2,mmgl)*log(mmgl/mmu)*pow4(mmsb1))/pow2(-mmgl + mmsb1) + (4*DeltaInv(
     mmt,mmsb1,mmst1)*log(mmsb1/mmu)*pow4(mmsb1))/pow2(-mmgl + mmsb1) + (4*
     DeltaInv(mmt,mmsb1,mmst2)*log(mmsb1/mmu)*pow4(mmsb1))/pow2(-mmgl + mmsb1)
     - (2*DeltaInv(mmt,mmst1,mmgl)*pow2(log(mmgl/mmu))*pow4(mmsb1))/(3.*pow2(-
     mmgl + mmsb1)) - (2*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmgl/mmu))*pow4(
     mmsb1))/(3.*pow2(-mmgl + mmsb1)) - (2*DeltaInv(mmt,mmsb1,mmst1)*pow2(log(
     mmsb1/mmu))*pow4(mmsb1))/(3.*pow2(-mmgl + mmsb1)) - (2*DeltaInv(mmt,mmsb1,
     mmst2)*pow2(log(mmsb1/mmu))*pow4(mmsb1))/(3.*pow2(-mmgl + mmsb1)) + (4*
     log(mmsb1/mmu)*log(mmt/mmu)*(-(mmsb1/(-mmgl + mmsb1)) - (3*mmsb1*mmst1*
     mmt*DeltaInv(mmt,mmsb1,mmst1))/(2.*(-mmgl + mmsb1)) - (3*mmsb1*mmst2*mmt*
     DeltaInv(mmt,mmsb1,mmst2))/(2.*(-mmgl + mmsb1)) - (mmst1*DeltaInv(mmt,
     mmsb1,mmst1)*pow2(mmsb1))/(-mmgl + mmsb1) - (mmt*DeltaInv(mmt,mmsb1,mmst1)
     *pow2(mmsb1))/(2.*(-mmgl + mmsb1)) - (mmst2*DeltaInv(mmt,mmsb1,mmst2)*
     pow2(mmsb1))/(-mmgl + mmsb1) - (mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))
     /(2.*(-mmgl + mmsb1)) - (3*mmsb1*mmt)/(2.*pow2(-mmgl + mmsb1)) + (2*pow2(
     mmsb1))/pow2(-mmgl + mmsb1) + (3*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(
     mmsb1))/pow2(-mmgl + mmsb1) + (3*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(
     mmsb1))/pow2(-mmgl + mmsb1) + (mmsb1*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1)
     )/(2.*(-mmgl + mmsb1)) - (DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow2(
     mmst1))/pow2(-mmgl + mmsb1) + (mmsb1*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2)
     )/(2.*(-mmgl + mmsb1)) - (DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*pow2(
     mmst2))/pow2(-mmgl + mmsb1) + (DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1))/(2.*
     (-mmgl + mmsb1)) + (DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1))/(2.*(-mmgl +
     mmsb1)) + (2*mmst1*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1))/pow2(-mmgl +
     mmsb1) + (mmt*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1))/pow2(-mmgl + mmsb1) +
     (2*mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1))/pow2(-mmgl + mmsb1) + (
     mmt*DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1))/pow2(-mmgl + mmsb1) + (2*mmt*
     pow2(mmsb1))/pow3(-mmgl + mmsb1) - (DeltaInv(mmt,mmsb1,mmst1)*pow4(mmsb1))
     /pow2(-mmgl + mmsb1) - (DeltaInv(mmt,mmsb1,mmst2)*pow4(mmsb1))/pow2(-mmgl
     + mmsb1)))/3. - (4*log(mmgl/mmu)*log(mmsb1/mmu)*pow4(mmsb1))/(3.*pow4(-
     mmgl + mmsb1)) + (2*pow2(log(mmgl/mmu))*pow4(mmsb1))/(3.*pow4(-mmgl +
     mmsb1)) + (2*pow2(log(mmsb1/mmu))*pow4(mmsb1))/(3.*pow4(-mmgl + mmsb1)) +
     pow3(Xb)*((-128*mmb*mmgl*mmsb1)/(9.*mgl*(-mmgl + mmsb1)*pow3(mmsb1 -
     mmsb2)) + (128*mmb*mmgl*mmsb2)/(9.*mgl*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)
     ) - (128*mmb*mmgl*mmsb1)/(9.*mgl*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) + (
     128*mmb*mmgl*mmsb2)/(9.*mgl*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) - (128*
     mmb*mmgl*mmsb1*log(mmgl/mmu))/(9.*mgl*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2))
     - (128*mmb*mmgl*mmsb2*log(mmgl/mmu))/(9.*mgl*(-mmgl + mmsb1)*pow3(mmsb1 -
     mmsb2)) + (128*mmb*mmgl*mmsb1*log(mmgl/mmu))/(9.*mgl*(-mmgl + mmsb2)*pow3(
     mmsb1 - mmsb2)) + (128*mmb*mmgl*mmsb2*log(mmgl/mmu))/(9.*mgl*(-mmgl +
     mmsb2)*pow3(mmsb1 - mmsb2)) + (256*mmb*mmgl*mmsb1*log(mmsb1/mmu))/(9.*mgl*
     (-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) + (128*mmb*mmgl*mmsb2*log(mmsb1/mmu))
     /(9.*mgl*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) + (128*mmb*mmgl*mmsb1*log(
     mmsb1/mmu))/(9.*mgl*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) + (128*mmb*mmgl*
     mmsb1*log(mmgl/mmu)*log(mmsb1/mmu))/(9.*mgl*(-mmgl + mmsb1)*pow3(mmsb1 -
     mmsb2)) + (256*mmb*mmgl*mmsb2*log(mmgl/mmu)*log(mmsb1/mmu))/(9.*mgl*(-mmgl
     + mmsb1)*pow3(mmsb1 - mmsb2)) - (128*mmb*mmgl*mmsb1*log(mmgl/mmu)*log(
     mmsb1/mmu))/(9.*mgl*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) - (256*mmb*mmgl*
     mmsb2*log(mmgl/mmu)*log(mmsb1/mmu))/(9.*mgl*(-mmgl + mmsb2)*pow3(mmsb1 -
     mmsb2)) - (128*mmb*mmgl*mmsb2*log(mmsb2/mmu))/(9.*mgl*(-mmgl + mmsb1)*
     pow3(mmsb1 - mmsb2)) - (128*mmb*mmgl*mmsb1*log(mmsb2/mmu))/(9.*mgl*(-mmgl
     + mmsb2)*pow3(mmsb1 - mmsb2)) - (256*mmb*mmgl*mmsb2*log(mmsb2/mmu))/(9.*
     mgl*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) - (128*mmb*mmgl*mmsb2*log(mmgl/
     mmu)*log(mmsb2/mmu))/(9.*mgl*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) + (128*
     mmb*mmgl*mmsb2*log(mmgl/mmu)*log(mmsb2/mmu))/(9.*mgl*(-mmgl + mmsb2)*pow3(
     mmsb1 - mmsb2)) + (128*mmb*mmgl*mmsb2*log(mmsb1/mmu)*log(mmsb2/mmu))/(9.*
     mgl*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) + (128*mmb*mmgl*mmsb1*log(mmsb1/
     mmu)*log(mmsb2/mmu))/(9.*mgl*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) + (256*
     mmb*mmgl*mmsb2*log(mmsb1/mmu)*log(mmsb2/mmu))/(9.*mgl*(-mmgl + mmsb2)*
     pow3(mmsb1 - mmsb2)) + (128*mmb*mmgl*mmsb1*mmsb2*log(mmgl/mmu))/(9.*mgl*
     pow2(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) - (128*mmb*mmgl*mmsb1*mmsb2*log(
     mmsb1/mmu))/(9.*mgl*pow2(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) - (128*mmb*
     mmgl*mmsb1*mmsb2*log(mmgl/mmu)*log(mmsb2/mmu))/(9.*mgl*pow2(-mmgl + mmsb1)
     *pow3(mmsb1 - mmsb2)) + (128*mmb*mmgl*mmsb1*mmsb2*log(mmsb1/mmu)*log(
     mmsb2/mmu))/(9.*mgl*pow2(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) - (128*mmb*
     mmgl*log(mmgl/mmu)*pow2(mmsb1))/(9.*mgl*pow2(-mmgl + mmsb1)*pow3(mmsb1 -
     mmsb2)) + (128*mmb*mmgl*log(mmsb1/mmu)*pow2(mmsb1))/(9.*mgl*pow2(-mmgl +
     mmsb1)*pow3(mmsb1 - mmsb2)) + (128*mmb*mmgl*log(mmgl/mmu)*log(mmsb1/mmu)*
     pow2(mmsb1))/(9.*mgl*pow2(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) - (128*mmb*
     mmgl*mmsb1*mmsb2*log(mmgl/mmu))/(9.*mgl*pow2(-mmgl + mmsb2)*pow3(mmsb1 -
     mmsb2)) + (128*mmb*mmgl*mmsb1*mmsb2*log(mmgl/mmu)*log(mmsb1/mmu))/(9.*mgl*
     pow2(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) + (128*mmb*mmgl*mmsb1*mmsb2*log(
     mmsb2/mmu))/(9.*mgl*pow2(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) - (128*mmb*
     mmgl*mmsb1*mmsb2*log(mmsb1/mmu)*log(mmsb2/mmu))/(9.*mgl*pow2(-mmgl +
     mmsb2)*pow3(mmsb1 - mmsb2)) + (128*mmb*mmgl*log(mmgl/mmu)*pow2(mmsb2))/(9.
     *mgl*pow2(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) - (128*mmb*mmgl*log(mmsb2/
     mmu)*pow2(mmsb2))/(9.*mgl*pow2(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) - (128*
     mmb*mmgl*log(mmgl/mmu)*log(mmsb2/mmu)*pow2(mmsb2))/(9.*mgl*pow2(-mmgl +
     mmsb2)*pow3(mmsb1 - mmsb2)) - (128*mmb*mmgl*mmsb1*pow2(log(mmsb1/mmu)))/(
     9.*mgl*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) - (256*mmb*mmgl*mmsb2*pow2(
     log(mmsb1/mmu)))/(9.*mgl*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) - (128*mmb*
     mmgl*pow2(mmsb1)*pow2(log(mmsb1/mmu)))/(9.*mgl*pow2(-mmgl + mmsb1)*pow3(
     mmsb1 - mmsb2)) - (128*mmb*mmgl*mmsb2*pow2(log(mmsb2/mmu)))/(9.*mgl*(-mmgl
     + mmsb2)*pow3(mmsb1 - mmsb2)) + (128*mmb*mmgl*pow2(mmsb2)*pow2(log(mmsb2/
     mmu)))/(9.*mgl*pow2(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) + (256*mmb*mmgl*
     log(mmgl/mmu)*log(mmsb1/mmu)*pow2(mmsb2))/(9.*mgl*(-mmgl + mmsb1)*pow4(
     mmsb1 - mmsb2)) - (256*mmb*mmgl*log(mmgl/mmu)*log(mmsb1/mmu)*pow2(mmsb2))/
     (9.*mgl*(-mmgl + mmsb2)*pow4(mmsb1 - mmsb2)) - (256*mmb*mmgl*log(mmgl/mmu)
     *log(mmsb2/mmu)*pow2(mmsb2))/(9.*mgl*(-mmgl + mmsb1)*pow4(mmsb1 - mmsb2))
     + (256*mmb*mmgl*log(mmgl/mmu)*log(mmsb2/mmu)*pow2(mmsb2))/(9.*mgl*(-mmgl +
     mmsb2)*pow4(mmsb1 - mmsb2)) + (256*mmb*mmgl*log(mmsb1/mmu)*log(mmsb2/mmu)*
     pow2(mmsb2))/(9.*mgl*(-mmgl + mmsb1)*pow4(mmsb1 - mmsb2)) + (256*mmb*mmgl*
     log(mmsb1/mmu)*log(mmsb2/mmu)*pow2(mmsb2))/(9.*mgl*(-mmgl + mmsb2)*pow4(
     mmsb1 - mmsb2)) - (256*mmb*mmgl*pow2(mmsb2)*pow2(log(mmsb1/mmu)))/(9.*mgl*
     (-mmgl + mmsb1)*pow4(mmsb1 - mmsb2)) - (256*mmb*mmgl*pow2(mmsb2)*pow2(log(
     mmsb2/mmu)))/(9.*mgl*(-mmgl + mmsb2)*pow4(mmsb1 - mmsb2))) - (4*log(mmgl/
     mmu)*log(mmsb1/mmu)*pow4(mmsb2))/(9.*pow2(-mmgl + mmsb1)*pow2(mmsb1 -
     mmsb2)) - (4*log(mmgl/mmu)*log(mmsb2/mmu)*pow4(mmsb2))/(9.*pow2(-mmgl +
     mmsb1)*pow2(mmsb1 - mmsb2)) + (4*log(mmsb1/mmu)*log(mmsb2/mmu)*pow4(mmsb2)
     )/(9.*pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (28*DeltaInv(mmt,mmsb2,
     mmst1)*pow4(mmsb2))/(3.*pow2(-mmgl + mmsb2)) - (4*zt2*DeltaInv(mmt,mmsb2,
     mmst1)*pow4(mmsb2))/(3.*pow2(-mmgl + mmsb2)) - (28*DeltaInv(mmt,mmsb2,
     mmst2)*pow4(mmsb2))/(3.*pow2(-mmgl + mmsb2)) - (4*zt2*DeltaInv(mmt,mmsb2,
     mmst2)*pow4(mmsb2))/(3.*pow2(-mmgl + mmsb2)) - (28*DeltaInv(mmt,mmst1,
     mmgl)*pow4(mmsb2))/(3.*pow2(-mmgl + mmsb2)) - (4*zt2*DeltaInv(mmt,mmst1,
     mmgl)*pow4(mmsb2))/(3.*pow2(-mmgl + mmsb2)) - (28*DeltaInv(mmt,mmst2,mmgl)
     *pow4(mmsb2))/(3.*pow2(-mmgl + mmsb2)) - (4*zt2*DeltaInv(mmt,mmst2,mmgl)*
     pow4(mmsb2))/(3.*pow2(-mmgl + mmsb2)) + (4*DeltaInv(mmt,mmst1,mmgl)*log(
     mmgl/mmu)*pow4(mmsb2))/pow2(-mmgl + mmsb2) + (4*DeltaInv(mmt,mmst2,mmgl)*
     log(mmgl/mmu)*pow4(mmsb2))/pow2(-mmgl + mmsb2) + (4*DeltaInv(mmt,mmsb2,
     mmst1)*log(mmsb2/mmu)*pow4(mmsb2))/pow2(-mmgl + mmsb2) + (4*DeltaInv(mmt,
     mmsb2,mmst2)*log(mmsb2/mmu)*pow4(mmsb2))/pow2(-mmgl + mmsb2) - (4*log(
     mmgl/mmu)*log(mmsb1/mmu)*pow4(mmsb2))/(9.*pow2(mmsb1 - mmsb2)*pow2(-mmgl +
     mmsb2)) - (4*log(mmgl/mmu)*log(mmsb2/mmu)*pow4(mmsb2))/(9.*pow2(mmsb1 -
     mmsb2)*pow2(-mmgl + mmsb2)) + (4*log(mmsb1/mmu)*log(mmsb2/mmu)*pow4(mmsb2)
     )/(9.*pow2(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (4*pow2(log(mmgl/mmu))*
     pow4(mmsb2))/(9.*pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (2*DeltaInv(
     mmt,mmst1,mmgl)*pow2(log(mmgl/mmu))*pow4(mmsb2))/(3.*pow2(-mmgl + mmsb2))
     - (2*DeltaInv(mmt,mmst2,mmgl)*pow2(log(mmgl/mmu))*pow4(mmsb2))/(3.*pow2(-
     mmgl + mmsb2)) + (4*pow2(log(mmgl/mmu))*pow4(mmsb2))/(9.*pow2(mmsb1 -
     mmsb2)*pow2(-mmgl + mmsb2)) - (2*DeltaInv(mmt,mmsb2,mmst1)*pow2(log(mmsb2/
     mmu))*pow4(mmsb2))/(3.*pow2(-mmgl + mmsb2)) - (2*DeltaInv(mmt,mmsb2,mmst2)
     *pow2(log(mmsb2/mmu))*pow4(mmsb2))/(3.*pow2(-mmgl + mmsb2)) - (8*log(mmgl/
     mmu)*log(mmsb1/mmu)*pow4(mmsb2))/(9.*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2))
     + (8*log(mmgl/mmu)*log(mmsb1/mmu)*pow4(mmsb2))/(9.*(-mmgl + mmsb2)*pow3(
     mmsb1 - mmsb2)) - (8*log(mmgl/mmu)*log(mmsb2/mmu)*pow4(mmsb2))/(9.*(-mmgl
     + mmsb1)*pow3(mmsb1 - mmsb2)) + (8*log(mmgl/mmu)*log(mmsb2/mmu)*pow4(
     mmsb2))/(9.*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) + (8*log(mmsb1/mmu)*log(
     mmsb2/mmu)*pow4(mmsb2))/(9.*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) - (8*log(
     mmsb1/mmu)*log(mmsb2/mmu)*pow4(mmsb2))/(9.*(-mmgl + mmsb2)*pow3(mmsb1 -
     mmsb2)) + (8*pow2(log(mmgl/mmu))*pow4(mmsb2))/(9.*(-mmgl + mmsb1)*pow3(
     mmsb1 - mmsb2)) - (8*pow2(log(mmgl/mmu))*pow4(mmsb2))/(9.*(-mmgl + mmsb2)*
     pow3(mmsb1 - mmsb2)) + (4*log(mmsb2/mmu)*log(mmt/mmu)*(-(mmsb2/(-mmgl +
     mmsb2)) - (3*mmsb2*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1))/(2.*(-mmgl +
     mmsb2)) - (3*mmsb2*mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2))/(2.*(-mmgl +
     mmsb2)) - (mmst1*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2))/(-mmgl + mmsb2) -
     (mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2))/(2.*(-mmgl + mmsb2)) - (mmst2*
     DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2))/(-mmgl + mmsb2) - (mmt*DeltaInv(
     mmt,mmsb2,mmst2)*pow2(mmsb2))/(2.*(-mmgl + mmsb2)) - (3*mmsb2*mmt)/(2.*
     pow2(-mmgl + mmsb2)) + (2*pow2(mmsb2))/pow2(-mmgl + mmsb2) + (3*mmst1*mmt*
     DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2))/pow2(-mmgl + mmsb2) + (3*mmst2*mmt*
     DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2))/pow2(-mmgl + mmsb2) + (mmsb2*
     DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/(2.*(-mmgl + mmsb2)) - (DeltaInv(
     mmt,mmsb2,mmst1)*pow2(mmsb2)*pow2(mmst1))/pow2(-mmgl + mmsb2) + (mmsb2*
     DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2))/(2.*(-mmgl + mmsb2)) - (DeltaInv(
     mmt,mmsb2,mmst2)*pow2(mmsb2)*pow2(mmst2))/pow2(-mmgl + mmsb2) + (DeltaInv(
     mmt,mmsb2,mmst1)*pow3(mmsb2))/(2.*(-mmgl + mmsb2)) + (DeltaInv(mmt,mmsb2,
     mmst2)*pow3(mmsb2))/(2.*(-mmgl + mmsb2)) + (2*mmst1*DeltaInv(mmt,mmsb2,
     mmst1)*pow3(mmsb2))/pow2(-mmgl + mmsb2) + (mmt*DeltaInv(mmt,mmsb2,mmst1)*
     pow3(mmsb2))/pow2(-mmgl + mmsb2) + (2*mmst2*DeltaInv(mmt,mmsb2,mmst2)*
     pow3(mmsb2))/pow2(-mmgl + mmsb2) + (mmt*DeltaInv(mmt,mmsb2,mmst2)*pow3(
     mmsb2))/pow2(-mmgl + mmsb2) + (2*mmt*pow2(mmsb2))/pow3(-mmgl + mmsb2) - (
     DeltaInv(mmt,mmsb2,mmst1)*pow4(mmsb2))/pow2(-mmgl + mmsb2) - (DeltaInv(
     mmt,mmsb2,mmst2)*pow4(mmsb2))/pow2(-mmgl + mmsb2)))/3. + (4*log(mmgl/mmu)*
     log(mmt/mmu)*(4 - (4*mmsb1)/(-mmgl + mmsb1) - (4*mmsb2)/(-mmgl + mmsb2) -
     (3*mmt)/(2.*(-mmgl + mmsb1)) - (3*mmt)/(2.*(-mmgl + mmsb2)) - 2*mmgl*
     mmsb1*DeltaInv(mmt,mmst1,mmgl) - 2*mmgl*mmsb2*DeltaInv(mmt,mmst1,mmgl) +
     4*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl) + 4*mmsb1*mmst1*DeltaInv(mmt,mmst1,
     mmgl) + 4*mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl) + 2*mmgl*mmt*DeltaInv(mmt,
     mmst1,mmgl) + 2*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl) + 2*mmsb2*mmt*DeltaInv(
     mmt,mmst1,mmgl) + 6*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl) - (6*mmsb1*mmst1*
     mmt*DeltaInv(mmt,mmst1,mmgl))/(-mmgl + mmsb1) - (6*mmsb2*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl))/(-mmgl + mmsb2) - 2*mmgl*mmsb1*DeltaInv(mmt,
     mmst2,mmgl) - 2*mmgl*mmsb2*DeltaInv(mmt,mmst2,mmgl) + 4*mmgl*mmst2*
     DeltaInv(mmt,mmst2,mmgl) + 4*mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl) + 4*
     mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl) + 2*mmgl*mmt*DeltaInv(mmt,mmst2,mmgl)
     + 2*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl) + 2*mmsb2*mmt*DeltaInv(mmt,mmst2,
     mmgl) + 6*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl) - (6*mmsb1*mmst2*mmt*
     DeltaInv(mmt,mmst2,mmgl))/(-mmgl + mmsb1) - (6*mmsb2*mmst2*mmt*DeltaInv(
     mmt,mmst2,mmgl))/(-mmgl + mmsb2) - 2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmgl) -
     2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmgl) - 3*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmsb1) - (6*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1))/(-mmgl + mmsb1) -
     (3*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1))/(-mmgl + mmsb1) - 3*DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmsb1) - (6*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmsb1))/(-mmgl + mmsb1) - (3*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/(-
     mmgl + mmsb1) + (7*mmsb1*mmt)/(2.*pow2(-mmgl + mmsb1)) + (2*pow2(mmsb1))/
     pow2(-mmgl + mmsb1) + (3*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1))/
     pow2(-mmgl + mmsb1) + (3*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/
     pow2(-mmgl + mmsb1) - 3*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2) - (6*mmst1*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2))/(-mmgl + mmsb2) - (3*mmt*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmsb2))/(-mmgl + mmsb2) - 3*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmsb2) - (6*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/(-mmgl +
     mmsb2) - (3*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/(-mmgl + mmsb2) + (
     7*mmsb2*mmt)/(2.*pow2(-mmgl + mmsb2)) + (2*pow2(mmsb2))/pow2(-mmgl +
     mmsb2) + (3*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2))/pow2(-mmgl +
     mmsb2) + (3*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/pow2(-mmgl +
     mmsb2) - 2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1) + (2*mmsb1*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmst1))/(-mmgl + mmsb1) + (2*mmsb2*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmst1))/(-mmgl + mmsb2) - (DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)
     *pow2(mmst1))/pow2(-mmgl + mmsb1) - (DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*
     pow2(mmst1))/pow2(-mmgl + mmsb2) - 2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2)
     + (2*mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(-mmgl + mmsb1) + (2*
     mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(-mmgl + mmsb2) - (DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(mmst2))/pow2(-mmgl + mmsb1) - (DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(mmst2))/pow2(-mmgl + mmsb2) + (4*
     DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1))/(-mmgl + mmsb1) + (4*DeltaInv(mmt,
     mmst2,mmgl)*pow3(mmsb1))/(-mmgl + mmsb1) + (2*mmst1*DeltaInv(mmt,mmst1,
     mmgl)*pow3(mmsb1))/pow2(-mmgl + mmsb1) + (mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow3(mmsb1))/pow2(-mmgl + mmsb1) + (2*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow3(
     mmsb1))/pow2(-mmgl + mmsb1) + (mmt*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1))/
     pow2(-mmgl + mmsb1) - (2*mmt*pow2(mmsb1))/pow3(-mmgl + mmsb1) + (4*
     DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2))/(-mmgl + mmsb2) + (4*DeltaInv(mmt,
     mmst2,mmgl)*pow3(mmsb2))/(-mmgl + mmsb2) + (2*mmst1*DeltaInv(mmt,mmst1,
     mmgl)*pow3(mmsb2))/pow2(-mmgl + mmsb2) + (mmt*DeltaInv(mmt,mmst1,mmgl)*
     pow3(mmsb2))/pow2(-mmgl + mmsb2) + (2*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow3(
     mmsb2))/pow2(-mmgl + mmsb2) + (mmt*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2))/
     pow2(-mmgl + mmsb2) - (2*mmt*pow2(mmsb2))/pow3(-mmgl + mmsb2) - (DeltaInv(
     mmt,mmst1,mmgl)*pow4(mmsb1))/pow2(-mmgl + mmsb1) - (DeltaInv(mmt,mmst2,
     mmgl)*pow4(mmsb1))/pow2(-mmgl + mmsb1) - (DeltaInv(mmt,mmst1,mmgl)*pow4(
     mmsb2))/pow2(-mmgl + mmsb2) - (DeltaInv(mmt,mmst2,mmgl)*pow4(mmsb2))/pow2(
     -mmgl + mmsb2)))/3. + (4*pow2(log(mmt/mmu))*(3 - (5*mmsb1)/(2.*(-mmgl +
     mmsb1)) - (5*mmsb2)/(2.*(-mmgl + mmsb2)) - (3*mmst1)/(4.*(-mmgl + mmsb1))
     - (3*mmst1)/(4.*(-mmgl + mmsb2)) - (3*mmst2)/(4.*(-mmgl + mmsb1)) - (3*
     mmst2)/(4.*(-mmgl + mmsb2)) - (3*mmt)/(4.*(-mmgl + mmsb1)) - (3*mmt)/(4.*(
     -mmgl + mmsb2)) - (3*mmsb1*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1))/(2.*(-mmgl
     + mmsb1)) - (3*mmsb1*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2))/(2.*(-mmgl +
     mmsb1)) - (3*mmsb2*mmst1*mmt*DeltaInv(mmt,mmsb2,mmst1))/(2.*(-mmgl +
     mmsb2)) - (3*mmsb2*mmst2*mmt*DeltaInv(mmt,mmsb2,mmst2))/(2.*(-mmgl +
     mmsb2)) - mmgl*mmsb1*DeltaInv(mmt,mmst1,mmgl) - mmgl*mmsb2*DeltaInv(mmt,
     mmst1,mmgl) + mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl) + mmsb1*mmst1*DeltaInv(
     mmt,mmst1,mmgl) + mmsb2*mmst1*DeltaInv(mmt,mmst1,mmgl) + mmgl*mmt*
     DeltaInv(mmt,mmst1,mmgl) + mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl) + mmsb2*mmt*
     DeltaInv(mmt,mmst1,mmgl) + 6*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl) - (6*
     mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl))/(-mmgl + mmsb1) - (6*mmsb2*
     mmst1*mmt*DeltaInv(mmt,mmst1,mmgl))/(-mmgl + mmsb2) - mmgl*mmsb1*DeltaInv(
     mmt,mmst2,mmgl) - mmgl*mmsb2*DeltaInv(mmt,mmst2,mmgl) + mmgl*mmst2*
     DeltaInv(mmt,mmst2,mmgl) + mmsb1*mmst2*DeltaInv(mmt,mmst2,mmgl) + mmsb2*
     mmst2*DeltaInv(mmt,mmst2,mmgl) + mmgl*mmt*DeltaInv(mmt,mmst2,mmgl) +
     mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl) + mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl) +
     6*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl) - (6*mmsb1*mmst2*mmt*DeltaInv(mmt,
     mmst2,mmgl))/(-mmgl + mmsb1) - (6*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)
     )/(-mmgl + mmsb2) - DeltaInv(mmt,mmst1,mmgl)*pow2(mmgl) - DeltaInv(mmt,
     mmst2,mmgl)*pow2(mmgl) - (mmst1*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1))/(4.
     *(-mmgl + mmsb1)) - (mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1))/(4.*(-mmgl
     + mmsb1)) - (mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))/(4.*(-mmgl +
     mmsb1)) - (mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))/(4.*(-mmgl + mmsb1))
     - (3*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1))/2. - (3*mmst1*DeltaInv(mmt,
     mmst1,mmgl)*pow2(mmsb1))/(2.*(-mmgl + mmsb1)) - (3*mmt*DeltaInv(mmt,mmst1,
     mmgl)*pow2(mmsb1))/(2.*(-mmgl + mmsb1)) - (3*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmsb1))/2. - (3*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/(2.*(-
     mmgl + mmsb1)) - (3*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/(2.*(-mmgl +
     mmsb1)) + (mmsb1*mmst1)/pow2(-mmgl + mmsb1) + (mmsb1*mmst2)/pow2(-mmgl +
     mmsb1) + (mmsb1*mmt)/pow2(-mmgl + mmsb1) + (2*pow2(mmsb1))/pow2(-mmgl +
     mmsb1) + (3*mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1))/pow2(-mmgl +
     mmsb1) + (3*mmst2*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))/pow2(-mmgl +
     mmsb1) + (3*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1))/pow2(-mmgl +
     mmsb1) + (3*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/pow2(-mmgl +
     mmsb1) - (mmst1*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2))/(4.*(-mmgl + mmsb2)
     ) - (mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2))/(4.*(-mmgl + mmsb2)) - (
     mmst2*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2))/(4.*(-mmgl + mmsb2)) - (mmt*
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
     DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/(4.*(-mmgl + mmsb2)) + DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmst1) - (mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))
     /(-mmgl + mmsb1) - (mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(-mmgl +
     mmsb2) - (mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(2.*(-mmgl + mmsb1)) -
     (mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(2.*(-mmgl + mmsb2)) + (mmsb1*
     mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1))/(2.*pow2(-mmgl + mmsb1)) + (
     mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(2.*pow2(-mmgl + mmsb1)) +
     (DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow2(mmst1))/(2.*pow2(-mmgl +
     mmsb1)) + (DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(mmst1))/(2.*pow2(-
     mmgl + mmsb1)) + (mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/(2.*
     pow2(-mmgl + mmsb2)) + (mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(
     2.*pow2(-mmgl + mmsb2)) + (DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2)*pow2(
     mmst1))/(2.*pow2(-mmgl + mmsb2)) + (DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*
     pow2(mmst1))/(2.*pow2(-mmgl + mmsb2)) - (mmsb1*DeltaInv(mmt,mmsb1,mmst2)*
     pow2(mmst2))/(4.*(-mmgl + mmsb1)) - (mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(
     mmst2))/(4.*(-mmgl + mmsb1)) - (mmsb2*DeltaInv(mmt,mmsb2,mmst2)*pow2(
     mmst2))/(4.*(-mmgl + mmsb2)) - (mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2))
     /(4.*(-mmgl + mmsb2)) + DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2) - (mmsb1*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(-mmgl + mmsb1) - (mmsb2*DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmst2))/(-mmgl + mmsb2) - (mmt*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmst2))/(2.*(-mmgl + mmsb1)) - (mmt*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmst2))/(2.*(-mmgl + mmsb2)) + (mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst2)*
     pow2(mmst2))/(2.*pow2(-mmgl + mmsb1)) + (mmsb1*mmt*DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmst2))/(2.*pow2(-mmgl + mmsb1)) + (DeltaInv(mmt,mmsb1,mmst2)*
     pow2(mmsb1)*pow2(mmst2))/(2.*pow2(-mmgl + mmsb1)) + (DeltaInv(mmt,mmst2,
     mmgl)*pow2(mmsb1)*pow2(mmst2))/(2.*pow2(-mmgl + mmsb1)) + (mmsb2*mmt*
     DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2))/(2.*pow2(-mmgl + mmsb2)) + (mmsb2*
     mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(2.*pow2(-mmgl + mmsb2)) + (
     DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow2(mmst2))/(2.*pow2(-mmgl + mmsb2)
     ) + (DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(mmst2))/(2.*pow2(-mmgl +
     mmsb2)) + (DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1))/(4.*(-mmgl + mmsb1)) + (
     DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1))/(4.*(-mmgl + mmsb1)) + (2*DeltaInv(
     mmt,mmst1,mmgl)*pow3(mmsb1))/(-mmgl + mmsb1) + (2*DeltaInv(mmt,mmst2,mmgl)
     *pow3(mmsb1))/(-mmgl + mmsb1) + (mmst1*DeltaInv(mmt,mmsb1,mmst1)*pow3(
     mmsb1))/(2.*pow2(-mmgl + mmsb1)) + (mmt*DeltaInv(mmt,mmsb1,mmst1)*pow3(
     mmsb1))/(2.*pow2(-mmgl + mmsb1)) + (mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow3(
     mmsb1))/(2.*pow2(-mmgl + mmsb1)) + (mmt*DeltaInv(mmt,mmsb1,mmst2)*pow3(
     mmsb1))/(2.*pow2(-mmgl + mmsb1)) + (mmst1*DeltaInv(mmt,mmst1,mmgl)*pow3(
     mmsb1))/(2.*pow2(-mmgl + mmsb1)) + (mmt*DeltaInv(mmt,mmst1,mmgl)*pow3(
     mmsb1))/(2.*pow2(-mmgl + mmsb1)) + (mmst2*DeltaInv(mmt,mmst2,mmgl)*pow3(
     mmsb1))/(2.*pow2(-mmgl + mmsb1)) + (mmt*DeltaInv(mmt,mmst2,mmgl)*pow3(
     mmsb1))/(2.*pow2(-mmgl + mmsb1)) + (DeltaInv(mmt,mmsb2,mmst1)*pow3(mmsb2))
     /(4.*(-mmgl + mmsb2)) + (DeltaInv(mmt,mmsb2,mmst2)*pow3(mmsb2))/(4.*(-mmgl
     + mmsb2)) + (2*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2))/(-mmgl + mmsb2) + (2*
     DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2))/(-mmgl + mmsb2) + (mmst1*DeltaInv(
     mmt,mmsb2,mmst1)*pow3(mmsb2))/(2.*pow2(-mmgl + mmsb2)) + (mmt*DeltaInv(
     mmt,mmsb2,mmst1)*pow3(mmsb2))/(2.*pow2(-mmgl + mmsb2)) + (mmst2*DeltaInv(
     mmt,mmsb2,mmst2)*pow3(mmsb2))/(2.*pow2(-mmgl + mmsb2)) + (mmt*DeltaInv(
     mmt,mmsb2,mmst2)*pow3(mmsb2))/(2.*pow2(-mmgl + mmsb2)) + (mmst1*DeltaInv(
     mmt,mmst1,mmgl)*pow3(mmsb2))/(2.*pow2(-mmgl + mmsb2)) + (mmt*DeltaInv(mmt,
     mmst1,mmgl)*pow3(mmsb2))/(2.*pow2(-mmgl + mmsb2)) + (mmst2*DeltaInv(mmt,
     mmst2,mmgl)*pow3(mmsb2))/(2.*pow2(-mmgl + mmsb2)) + (mmt*DeltaInv(mmt,
     mmst2,mmgl)*pow3(mmsb2))/(2.*pow2(-mmgl + mmsb2)) + (DeltaInv(mmt,mmsb1,
     mmst1)*pow3(mmst1))/(4.*(-mmgl + mmsb1)) + (DeltaInv(mmt,mmsb2,mmst1)*
     pow3(mmst1))/(4.*(-mmgl + mmsb2)) + (DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1))
     /(2.*(-mmgl + mmsb1)) + (DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1))/(2.*(-mmgl
     + mmsb2)) - (mmsb1*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmst1))/(2.*pow2(-mmgl +
     mmsb1)) - (mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1))/(2.*pow2(-mmgl +
     mmsb1)) - (mmsb2*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmst1))/(2.*pow2(-mmgl +
     mmsb2)) - (mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1))/(2.*pow2(-mmgl +
     mmsb2)) + (DeltaInv(mmt,mmsb1,mmst2)*pow3(mmst2))/(4.*(-mmgl + mmsb1)) + (
     DeltaInv(mmt,mmsb2,mmst2)*pow3(mmst2))/(4.*(-mmgl + mmsb2)) + (DeltaInv(
     mmt,mmst2,mmgl)*pow3(mmst2))/(2.*(-mmgl + mmsb1)) + (DeltaInv(mmt,mmst2,
     mmgl)*pow3(mmst2))/(2.*(-mmgl + mmsb2)) - (mmsb1*DeltaInv(mmt,mmsb1,mmst2)
     *pow3(mmst2))/(2.*pow2(-mmgl + mmsb1)) - (mmsb1*DeltaInv(mmt,mmst2,mmgl)*
     pow3(mmst2))/(2.*pow2(-mmgl + mmsb1)) - (mmsb2*DeltaInv(mmt,mmsb2,mmst2)*
     pow3(mmst2))/(2.*pow2(-mmgl + mmsb2)) - (mmsb2*DeltaInv(mmt,mmst2,mmgl)*
     pow3(mmst2))/(2.*pow2(-mmgl + mmsb2)) - (DeltaInv(mmt,mmsb1,mmst1)*pow4(
     mmsb1))/(2.*pow2(-mmgl + mmsb1)) - (DeltaInv(mmt,mmsb1,mmst2)*pow4(mmsb1))
     /(2.*pow2(-mmgl + mmsb1)) - (DeltaInv(mmt,mmst1,mmgl)*pow4(mmsb1))/(2.*
     pow2(-mmgl + mmsb1)) - (DeltaInv(mmt,mmst2,mmgl)*pow4(mmsb1))/(2.*pow2(-
     mmgl + mmsb1)) - (DeltaInv(mmt,mmsb2,mmst1)*pow4(mmsb2))/(2.*pow2(-mmgl +
     mmsb2)) - (DeltaInv(mmt,mmsb2,mmst2)*pow4(mmsb2))/(2.*pow2(-mmgl + mmsb2))
     - (DeltaInv(mmt,mmst1,mmgl)*pow4(mmsb2))/(2.*pow2(-mmgl + mmsb2)) - (
     DeltaInv(mmt,mmst2,mmgl)*pow4(mmsb2))/(2.*pow2(-mmgl + mmsb2))))/3. + (4*
     log(mmt/mmu)*(-10.333333333333334 + (15*mmsb1)/(-mmgl + mmsb1) + (15*
     mmsb2)/(-mmgl + mmsb2) + (9*mmst1)/(2.*(-mmgl + mmsb1)) + (9*mmst1)/(2.*(-
     mmgl + mmsb2)) + (9*mmst2)/(2.*(-mmgl + mmsb1)) + (9*mmst2)/(2.*(-mmgl +
     mmsb2)) + (6*mmt)/(-mmgl + mmsb1) + (6*mmt)/(-mmgl + mmsb2) + (9*mmsb1*
     mmst1*mmt*DeltaInv(mmt,mmsb1,mmst1))/(-mmgl + mmsb1) + (9*mmsb1*mmst2*mmt*
     DeltaInv(mmt,mmsb1,mmst2))/(-mmgl + mmsb1) + (9*mmsb2*mmst1*mmt*DeltaInv(
     mmt,mmsb2,mmst1))/(-mmgl + mmsb2) + (9*mmsb2*mmst2*mmt*DeltaInv(mmt,mmsb2,
     mmst2))/(-mmgl + mmsb2) + 6*mmgl*mmsb1*DeltaInv(mmt,mmst1,mmgl) + 6*mmgl*
     mmsb2*DeltaInv(mmt,mmst1,mmgl) - 6*mmgl*mmst1*DeltaInv(mmt,mmst1,mmgl) -
     6*mmsb1*mmst1*DeltaInv(mmt,mmst1,mmgl) - 6*mmsb2*mmst1*DeltaInv(mmt,mmst1,
     mmgl) - 6*mmgl*mmt*DeltaInv(mmt,mmst1,mmgl) - 6*mmsb1*mmt*DeltaInv(mmt,
     mmst1,mmgl) - 6*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl) - 36*mmst1*mmt*
     DeltaInv(mmt,mmst1,mmgl) + (36*mmsb1*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl))/(
     -mmgl + mmsb1) + (36*mmsb2*mmst1*mmt*DeltaInv(mmt,mmst1,mmgl))/(-mmgl +
     mmsb2) + 6*mmgl*mmsb1*DeltaInv(mmt,mmst2,mmgl) + 6*mmgl*mmsb2*DeltaInv(
     mmt,mmst2,mmgl) - 6*mmgl*mmst2*DeltaInv(mmt,mmst2,mmgl) - 6*mmsb1*mmst2*
     DeltaInv(mmt,mmst2,mmgl) - 6*mmsb2*mmst2*DeltaInv(mmt,mmst2,mmgl) - 6*
     mmgl*mmt*DeltaInv(mmt,mmst2,mmgl) - 6*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl) -
     6*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl) - 36*mmst2*mmt*DeltaInv(mmt,mmst2,
     mmgl) + (36*mmsb1*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl))/(-mmgl + mmsb1) + (
     36*mmsb2*mmst2*mmt*DeltaInv(mmt,mmst2,mmgl))/(-mmgl + mmsb2) + 6*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmgl) + 6*DeltaInv(mmt,mmst2,mmgl)*pow2(mmgl) + (3*
     mmst1*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1))/(2.*(-mmgl + mmsb1)) + (3*
     mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1))/(2.*(-mmgl + mmsb1)) + (3*
     mmst2*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))/(2.*(-mmgl + mmsb1)) + (3*
     mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))/(2.*(-mmgl + mmsb1)) + 9*
     DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1) + (9*mmst1*DeltaInv(mmt,mmst1,mmgl)*
     pow2(mmsb1))/(-mmgl + mmsb1) + (9*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)
     )/(-mmgl + mmsb1) + 9*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1) + (9*mmst2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/(-mmgl + mmsb1) + (9*mmt*DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmsb1))/(-mmgl + mmsb1) - (6*mmsb1*mmst1)/pow2(-mmgl
     + mmsb1) - (6*mmsb1*mmst2)/pow2(-mmgl + mmsb1) - (8*mmsb1*mmt)/pow2(-mmgl
     + mmsb1) - (12*pow2(mmsb1))/pow2(-mmgl + mmsb1) - (18*mmst1*mmt*DeltaInv(
     mmt,mmsb1,mmst1)*pow2(mmsb1))/pow2(-mmgl + mmsb1) - (18*mmst2*mmt*
     DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1))/pow2(-mmgl + mmsb1) - (18*mmst1*
     mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1))/pow2(-mmgl + mmsb1) - (18*mmst2*
     mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1))/pow2(-mmgl + mmsb1) + (3*mmst1*
     DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2))/(2.*(-mmgl + mmsb2)) + (3*mmt*
     DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2))/(2.*(-mmgl + mmsb2)) + (3*mmst2*
     DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2))/(2.*(-mmgl + mmsb2)) + (3*mmt*
     DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2))/(2.*(-mmgl + mmsb2)) + 9*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmsb2) + (9*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmsb2))/(-mmgl + mmsb2) + (9*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2))/(-
     mmgl + mmsb2) + 9*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2) + (9*mmst2*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/(-mmgl + mmsb2) + (9*mmt*DeltaInv(
     mmt,mmst2,mmgl)*pow2(mmsb2))/(-mmgl + mmsb2) - (6*mmsb2*mmst1)/pow2(-mmgl
     + mmsb2) - (6*mmsb2*mmst2)/pow2(-mmgl + mmsb2) - (8*mmsb2*mmt)/pow2(-mmgl
     + mmsb2) - (12*pow2(mmsb2))/pow2(-mmgl + mmsb2) - (18*mmst1*mmt*DeltaInv(
     mmt,mmsb2,mmst1)*pow2(mmsb2))/pow2(-mmgl + mmsb2) - (18*mmst2*mmt*
     DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2))/pow2(-mmgl + mmsb2) - (18*mmst1*
     mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2))/pow2(-mmgl + mmsb2) - (18*mmst2*
     mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2))/pow2(-mmgl + mmsb2) + (3*mmsb1*
     DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1))/(2.*(-mmgl + mmsb1)) + (3*mmt*
     DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1))/(2.*(-mmgl + mmsb1)) + (3*mmsb2*
     DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/(2.*(-mmgl + mmsb2)) + (3*mmt*
     DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/(2.*(-mmgl + mmsb2)) - 6*DeltaInv(
     mmt,mmst1,mmgl)*pow2(mmst1) + (6*mmsb1*DeltaInv(mmt,mmst1,mmgl)*pow2(
     mmst1))/(-mmgl + mmsb1) + (6*mmsb2*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(
     -mmgl + mmsb2) + (3*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(-mmgl +
     mmsb1) + (3*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/(-mmgl + mmsb2) - (
     3*mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmst1))/pow2(-mmgl + mmsb1) - (
     3*mmsb1*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/pow2(-mmgl + mmsb1) - (
     3*DeltaInv(mmt,mmsb1,mmst1)*pow2(mmsb1)*pow2(mmst1))/pow2(-mmgl + mmsb1) -
     (3*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb1)*pow2(mmst1))/pow2(-mmgl + mmsb1) -
     (3*mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmst1))/pow2(-mmgl + mmsb2) -
     (3*mmsb2*mmt*DeltaInv(mmt,mmst1,mmgl)*pow2(mmst1))/pow2(-mmgl + mmsb2) - (
     3*DeltaInv(mmt,mmsb2,mmst1)*pow2(mmsb2)*pow2(mmst1))/pow2(-mmgl + mmsb2) -
     (3*DeltaInv(mmt,mmst1,mmgl)*pow2(mmsb2)*pow2(mmst1))/pow2(-mmgl + mmsb2) +
     (3*mmsb1*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2))/(2.*(-mmgl + mmsb1)) + (3*
     mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2))/(2.*(-mmgl + mmsb1)) + (3*
     mmsb2*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2))/(2.*(-mmgl + mmsb2)) + (3*
     mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2))/(2.*(-mmgl + mmsb2)) - 6*
     DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2) + (6*mmsb1*DeltaInv(mmt,mmst2,mmgl)*
     pow2(mmst2))/(-mmgl + mmsb1) + (6*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow2(
     mmst2))/(-mmgl + mmsb2) + (3*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(-
     mmgl + mmsb1) + (3*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/(-mmgl +
     mmsb2) - (3*mmsb1*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmst2))/pow2(-mmgl +
     mmsb1) - (3*mmsb1*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/pow2(-mmgl +
     mmsb1) - (3*DeltaInv(mmt,mmsb1,mmst2)*pow2(mmsb1)*pow2(mmst2))/pow2(-mmgl
     + mmsb1) - (3*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb1)*pow2(mmst2))/pow2(-mmgl
     + mmsb1) - (3*mmsb2*mmt*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmst2))/pow2(-mmgl
     + mmsb2) - (3*mmsb2*mmt*DeltaInv(mmt,mmst2,mmgl)*pow2(mmst2))/pow2(-mmgl +
     mmsb2) - (3*DeltaInv(mmt,mmsb2,mmst2)*pow2(mmsb2)*pow2(mmst2))/pow2(-mmgl
     + mmsb2) - (3*DeltaInv(mmt,mmst2,mmgl)*pow2(mmsb2)*pow2(mmst2))/pow2(-mmgl
     + mmsb2) - (3*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmsb1))/(2.*(-mmgl + mmsb1))
     - (3*DeltaInv(mmt,mmsb1,mmst2)*pow3(mmsb1))/(2.*(-mmgl + mmsb1)) - (12*
     DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1))/(-mmgl + mmsb1) - (12*DeltaInv(mmt,
     mmst2,mmgl)*pow3(mmsb1))/(-mmgl + mmsb1) - (3*mmst1*DeltaInv(mmt,mmsb1,
     mmst1)*pow3(mmsb1))/pow2(-mmgl + mmsb1) - (3*mmt*DeltaInv(mmt,mmsb1,mmst1)
     *pow3(mmsb1))/pow2(-mmgl + mmsb1) - (3*mmst2*DeltaInv(mmt,mmsb1,mmst2)*
     pow3(mmsb1))/pow2(-mmgl + mmsb1) - (3*mmt*DeltaInv(mmt,mmsb1,mmst2)*pow3(
     mmsb1))/pow2(-mmgl + mmsb1) - (3*mmst1*DeltaInv(mmt,mmst1,mmgl)*pow3(
     mmsb1))/pow2(-mmgl + mmsb1) - (3*mmt*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb1))
     /pow2(-mmgl + mmsb1) - (3*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1))/
     pow2(-mmgl + mmsb1) - (3*mmt*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb1))/pow2(-
     mmgl + mmsb1) - (3*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmsb2))/(2.*(-mmgl +
     mmsb2)) - (3*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmsb2))/(2.*(-mmgl + mmsb2)) -
     (12*DeltaInv(mmt,mmst1,mmgl)*pow3(mmsb2))/(-mmgl + mmsb2) - (12*DeltaInv(
     mmt,mmst2,mmgl)*pow3(mmsb2))/(-mmgl + mmsb2) - (3*mmst1*DeltaInv(mmt,
     mmsb2,mmst1)*pow3(mmsb2))/pow2(-mmgl + mmsb2) - (3*mmt*DeltaInv(mmt,mmsb2,
     mmst1)*pow3(mmsb2))/pow2(-mmgl + mmsb2) - (3*mmst2*DeltaInv(mmt,mmsb2,
     mmst2)*pow3(mmsb2))/pow2(-mmgl + mmsb2) - (3*mmt*DeltaInv(mmt,mmsb2,mmst2)
     *pow3(mmsb2))/pow2(-mmgl + mmsb2) - (3*mmst1*DeltaInv(mmt,mmst1,mmgl)*
     pow3(mmsb2))/pow2(-mmgl + mmsb2) - (3*mmt*DeltaInv(mmt,mmst1,mmgl)*pow3(
     mmsb2))/pow2(-mmgl + mmsb2) - (3*mmst2*DeltaInv(mmt,mmst2,mmgl)*pow3(
     mmsb2))/pow2(-mmgl + mmsb2) - (3*mmt*DeltaInv(mmt,mmst2,mmgl)*pow3(mmsb2))
     /pow2(-mmgl + mmsb2) - (3*DeltaInv(mmt,mmsb1,mmst1)*pow3(mmst1))/(2.*(-
     mmgl + mmsb1)) - (3*DeltaInv(mmt,mmsb2,mmst1)*pow3(mmst1))/(2.*(-mmgl +
     mmsb2)) - (3*DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1))/(-mmgl + mmsb1) - (3*
     DeltaInv(mmt,mmst1,mmgl)*pow3(mmst1))/(-mmgl + mmsb2) + (3*mmsb1*DeltaInv(
     mmt,mmsb1,mmst1)*pow3(mmst1))/pow2(-mmgl + mmsb1) + (3*mmsb1*DeltaInv(mmt,
     mmst1,mmgl)*pow3(mmst1))/pow2(-mmgl + mmsb1) + (3*mmsb2*DeltaInv(mmt,
     mmsb2,mmst1)*pow3(mmst1))/pow2(-mmgl + mmsb2) + (3*mmsb2*DeltaInv(mmt,
     mmst1,mmgl)*pow3(mmst1))/pow2(-mmgl + mmsb2) - (3*DeltaInv(mmt,mmsb1,
     mmst2)*pow3(mmst2))/(2.*(-mmgl + mmsb1)) - (3*DeltaInv(mmt,mmsb2,mmst2)*
     pow3(mmst2))/(2.*(-mmgl + mmsb2)) - (3*DeltaInv(mmt,mmst2,mmgl)*pow3(
     mmst2))/(-mmgl + mmsb1) - (3*DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2))/(-mmgl
     + mmsb2) + (3*mmsb1*DeltaInv(mmt,mmsb1,mmst2)*pow3(mmst2))/pow2(-mmgl +
     mmsb1) + (3*mmsb1*DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2))/pow2(-mmgl +
     mmsb1) + (3*mmsb2*DeltaInv(mmt,mmsb2,mmst2)*pow3(mmst2))/pow2(-mmgl +
     mmsb2) + (3*mmsb2*DeltaInv(mmt,mmst2,mmgl)*pow3(mmst2))/pow2(-mmgl +
     mmsb2) + (3*DeltaInv(mmt,mmsb1,mmst1)*pow4(mmsb1))/pow2(-mmgl + mmsb1) + (
     3*DeltaInv(mmt,mmsb1,mmst2)*pow4(mmsb1))/pow2(-mmgl + mmsb1) + (3*
     DeltaInv(mmt,mmst1,mmgl)*pow4(mmsb1))/pow2(-mmgl + mmsb1) + (3*DeltaInv(
     mmt,mmst2,mmgl)*pow4(mmsb1))/pow2(-mmgl + mmsb1) + (3*DeltaInv(mmt,mmsb2,
     mmst1)*pow4(mmsb2))/pow2(-mmgl + mmsb2) + (3*DeltaInv(mmt,mmsb2,mmst2)*
     pow4(mmsb2))/pow2(-mmgl + mmsb2) + (3*DeltaInv(mmt,mmst1,mmgl)*pow4(mmsb2)
     )/pow2(-mmgl + mmsb2) + (3*DeltaInv(mmt,mmst2,mmgl)*pow4(mmsb2))/pow2(-
     mmgl + mmsb2)))/3. - (4*log(mmgl/mmu)*log(mmsb2/mmu)*pow4(mmsb2))/(3.*
     pow4(-mmgl + mmsb2)) + (2*pow2(log(mmgl/mmu))*pow4(mmsb2))/(3.*pow4(-mmgl
     + mmsb2)) + (2*pow2(log(mmsb2/mmu))*pow4(mmsb2))/(3.*pow4(-mmgl + mmsb2))
     + pow2(Xb)*((2560*mmb)/(9.*pow2(mmsb1 - mmsb2)) - (3284*mmb*mmsb1)/(9.*(-
     mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (5152*mmb*mmsb2)/(9.*(-mmgl + mmsb1)*
     pow2(mmsb1 - mmsb2)) - (32*mmb*mmsb1)/(3.*(-mmgl + mmsb2)*pow2(mmsb1 -
     mmsb2)) + (1772*mmb*mmsb2)/(9.*(-mmgl + mmsb2)*pow2(mmsb1 - mmsb2)) + (
     512*mmb*zt2)/(9.*pow2(mmsb1 - mmsb2)) - (556*mmb*mmsb1*zt2)/(9.*(-mmgl +
     mmsb1)*pow2(mmsb1 - mmsb2)) - (80*mmb*mmsb2*zt2)/((-mmgl + mmsb1)*pow2(
     mmsb1 - mmsb2)) + (164*mmb*mmsb2*zt2)/(9.*(-mmgl + mmsb2)*pow2(mmsb1 -
     mmsb2)) - (556*mmb*Fin20(mmsb1,mmgl,mmu))/(9.*(-mmgl + mmsb1)*pow2(mmsb1 -
     mmsb2)) + (44*mmb*Fin20(mmsb1,mmgl,mmu))/(9.*(-mmgl + mmsb2)*pow2(mmsb1 -
     mmsb2)) - (44*mmb*Fin20(mmsb1,mmsb2,mmu))/(9.*(-mmgl + mmsb1)*pow2(mmsb1 -
     mmsb2)) - (44*mmb*Fin20(mmsb1,mmsb2,mmu))/(9.*(-mmgl + mmsb2)*pow2(mmsb1 -
     mmsb2)) - (20*mmb*Fin20(mmsb2,mmgl,mmu))/(3.*(-mmgl + mmsb1)*pow2(mmsb1 -
     mmsb2)) + (20*mmb*Fin20(mmsb2,mmgl,mmu))/(3.*(-mmgl + mmsb2)*pow2(mmsb1 -
     mmsb2)) - (512*mmb*log(mmgl/mmu))/(3.*pow2(mmsb1 - mmsb2)) + (1472*mmb*
     mmsb1*log(mmgl/mmu))/(9.*(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) + (3296*mmb*
     mmsb2*log(mmgl/mmu))/(9.*(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) + (64*mmb*
     mmsb1*log(mmgl/mmu))/(9.*(-mmgl + mmsb2)*pow2(mmsb1 - mmsb2)) - (1760*mmb*
     mmsb2*log(mmgl/mmu))/(9.*(-mmgl + mmsb2)*pow2(mmsb1 - mmsb2)) - (256*mmb*
     log(mmsb1/mmu))/(3.*pow2(mmsb1 - mmsb2)) + (168*mmb*mmsb1*log(mmsb1/mmu))/
     ((-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) + (248*mmb*mmsb2*log(mmsb1/mmu))/(3.
     *(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) + (32*mmb*mmsb1*log(mmsb1/mmu))/(3.*
     (-mmgl + mmsb2)*pow2(mmsb1 - mmsb2)) - (296*mmb*mmsb2*log(mmsb1/mmu))/(9.*
     (-mmgl + mmsb2)*pow2(mmsb1 - mmsb2)) + (256*mmb*log(mmgl/mmu)*log(mmsb1/
     mmu))/(9.*pow2(mmsb1 - mmsb2)) - (448*mmb*mmsb1*log(mmgl/mmu)*log(mmsb1/
     mmu))/(9.*(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (52*mmb*mmsb2*log(mmgl/
     mmu)*log(mmsb1/mmu))/(9.*(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (64*mmb*
     mmsb1*log(mmgl/mmu)*log(mmsb1/mmu))/(9.*(-mmgl + mmsb2)*pow2(mmsb1 -
     mmsb2)) + (52*mmb*mmsb2*log(mmgl/mmu)*log(mmsb1/mmu))/(9.*(-mmgl + mmsb2)*
     pow2(mmsb1 - mmsb2)) - (256*mmb*log(mmsb2/mmu))/(9.*pow2(mmsb1 - mmsb2)) +
     (376*mmb*mmsb2*log(mmsb2/mmu))/(9.*(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) -
     (64*mmb*mmsb1*log(mmsb2/mmu))/(9.*(-mmgl + mmsb2)*pow2(mmsb1 - mmsb2)) + (
     80*mmb*mmsb2*log(mmsb2/mmu))/((-mmgl + mmsb2)*pow2(mmsb1 - mmsb2)) + (256*
     mmb*log(mmgl/mmu)*log(mmsb2/mmu))/(9.*pow2(mmsb1 - mmsb2)) - (116*mmb*
     mmsb2*log(mmgl/mmu)*log(mmsb2/mmu))/(9.*(-mmgl + mmsb1)*pow2(mmsb1 -
     mmsb2)) - (44*mmb*mmsb2*log(mmgl/mmu)*log(mmsb2/mmu))/((-mmgl + mmsb2)*
     pow2(mmsb1 - mmsb2)) + (4*mmb*mmsb2*log(mmsb1/mmu)*log(mmsb2/mmu))/(3.*(-
     mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) + (64*mmb*mmsb1*log(mmsb1/mmu)*log(
     mmsb2/mmu))/(9.*(-mmgl + mmsb2)*pow2(mmsb1 - mmsb2)) + (52*mmb*mmsb2*log(
     mmsb1/mmu)*log(mmsb2/mmu))/(9.*(-mmgl + mmsb2)*pow2(mmsb1 - mmsb2)) + (64*
     mmb*mmsb1*mmsb2)/(9.*pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) + (104*mmb*
     mmsb1*Fin20(mmsb1,mmsb2,mmu))/(9.*pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2))
     - (104*mmb*mmsb1*Fin20(mmsb2,mmgl,mmu))/(9.*pow2(-mmgl + mmsb1)*pow2(mmsb1
     - mmsb2)) - (16*mmb*mmsb1*mmsb2*log(mmgl/mmu))/(pow2(-mmgl + mmsb1)*pow2(
     mmsb1 - mmsb2)) + (16*mmb*mmsb1*mmsb2*log(mmsb1/mmu))/(pow2(-mmgl + mmsb1)
     *pow2(mmsb1 - mmsb2)) - (32*mmb*mmsb1*mmsb2*log(mmgl/mmu)*log(mmsb1/mmu))/
     (9.*pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (64*mmb*mmsb1*mmsb2*log(
     mmsb2/mmu))/(9.*pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) + (32*mmb*mmsb1*
     mmsb2*log(mmgl/mmu)*log(mmsb2/mmu))/(3.*pow2(-mmgl + mmsb1)*pow2(mmsb1 -
     mmsb2)) - (32*mmb*mmsb1*mmsb2*log(mmsb1/mmu)*log(mmsb2/mmu))/(3.*pow2(-
     mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (8*mmb*pow2(mmsb1))/(pow2(-mmgl +
     mmsb1)*pow2(mmsb1 - mmsb2)) + (328*mmb*log(mmgl/mmu)*pow2(mmsb1))/(3.*
     pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (920*mmb*log(mmsb1/mmu)*pow2(
     mmsb1))/(9.*pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) + (12*mmb*log(mmgl/
     mmu)*log(mmsb1/mmu)*pow2(mmsb1))/(pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2))
     - (16*mmb*log(mmgl/mmu)*pow2(mmsb2))/(9.*pow2(-mmgl + mmsb1)*pow2(mmsb1 -
     mmsb2)) + (16*mmb*log(mmsb1/mmu)*pow2(mmsb2))/(9.*pow2(-mmgl + mmsb1)*
     pow2(mmsb1 - mmsb2)) - (16*mmb*log(mmgl/mmu)*log(mmsb1/mmu)*pow2(mmsb2))/(
     3.*pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (16*mmb*log(mmgl/mmu)*log(
     mmsb2/mmu)*pow2(mmsb2))/(3.*pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) + (
     16*mmb*log(mmsb1/mmu)*log(mmsb2/mmu)*pow2(mmsb2))/(3.*pow2(-mmgl + mmsb1)*
     pow2(mmsb1 - mmsb2)) + (64*mmb*mmsb1*mmsb2)/(9.*pow2(mmsb1 - mmsb2)*pow2(-
     mmgl + mmsb2)) - (104*mmb*mmsb2*Fin20(mmsb1,mmgl,mmu))/(9.*pow2(mmsb1 -
     mmsb2)*pow2(-mmgl + mmsb2)) + (104*mmb*mmsb2*Fin20(mmsb1,mmsb2,mmu))/(9.*
     pow2(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (128*mmb*mmsb1*mmsb2*log(mmgl/
     mmu))/(9.*pow2(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (64*mmb*mmsb1*mmsb2*
     log(mmsb1/mmu))/(9.*pow2(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (128*mmb*
     mmsb1*mmsb2*log(mmgl/mmu)*log(mmsb1/mmu))/(9.*pow2(mmsb1 - mmsb2)*pow2(-
     mmgl + mmsb2)) + (128*mmb*mmsb1*mmsb2*log(mmsb2/mmu))/(9.*pow2(mmsb1 -
     mmsb2)*pow2(-mmgl + mmsb2)) - (128*mmb*mmsb1*mmsb2*log(mmsb1/mmu)*log(
     mmsb2/mmu))/(9.*pow2(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (8*mmb*pow2(
     mmsb2))/(pow2(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (1000*mmb*log(mmgl/
     mmu)*pow2(mmsb2))/(9.*pow2(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (16*mmb*
     log(mmgl/mmu)*log(mmsb1/mmu)*pow2(mmsb2))/(9.*pow2(mmsb1 - mmsb2)*pow2(-
     mmgl + mmsb2)) - (104*mmb*log(mmsb2/mmu)*pow2(mmsb2))/(pow2(mmsb1 - mmsb2)
     *pow2(-mmgl + mmsb2)) + (124*mmb*log(mmgl/mmu)*log(mmsb2/mmu)*pow2(mmsb2))
     /(9.*pow2(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) - (16*mmb*log(mmsb1/mmu)*
     log(mmsb2/mmu)*pow2(mmsb2))/(9.*pow2(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) +
     (256*mmb*pow2(log(mmgl/mmu)))/(9.*pow2(mmsb1 - mmsb2)) - (256*mmb*mmsb1*
     pow2(log(mmgl/mmu)))/(9.*(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (308*mmb*
     mmsb2*pow2(log(mmgl/mmu)))/(3.*(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) + (
     668*mmb*mmsb2*pow2(log(mmgl/mmu)))/(9.*(-mmgl + mmsb2)*pow2(mmsb1 - mmsb2)
     ) + (32*mmb*mmsb1*mmsb2*pow2(log(mmgl/mmu)))/(9.*pow2(-mmgl + mmsb1)*pow2(
     mmsb1 - mmsb2)) - (350*mmb*pow2(mmsb1)*pow2(log(mmgl/mmu)))/(9.*pow2(-mmgl
     + mmsb1)*pow2(mmsb1 - mmsb2)) + (16*mmb*pow2(mmsb2)*pow2(log(mmgl/mmu)))/(
     3.*pow2(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (122*mmb*pow2(mmsb2)*pow2(
     log(mmgl/mmu)))/(3.*pow2(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (256*mmb*
     pow2(log(mmsb1/mmu)))/(9.*pow2(mmsb1 - mmsb2)) - (136*mmb*mmsb1*pow2(log(
     mmsb1/mmu)))/(3.*(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) - (308*mmb*mmsb2*
     pow2(log(mmsb1/mmu)))/(9.*(-mmgl + mmsb1)*pow2(mmsb1 - mmsb2)) + (52*mmb*
     mmsb2*pow2(log(mmsb1/mmu)))/(9.*(-mmgl + mmsb2)*pow2(mmsb1 - mmsb2)) + (
     242*mmb*pow2(mmsb1)*pow2(log(mmsb1/mmu)))/(9.*pow2(-mmgl + mmsb1)*pow2(
     mmsb1 - mmsb2)) - (52*mmb*mmsb2*pow2(log(mmsb2/mmu)))/(9.*(-mmgl + mmsb1)*
     pow2(mmsb1 - mmsb2)) - (100*mmb*mmsb2*pow2(log(mmsb2/mmu)))/(9.*(-mmgl +
     mmsb2)*pow2(mmsb1 - mmsb2)) + (242*mmb*pow2(mmsb2)*pow2(log(mmsb2/mmu)))/(
     9.*pow2(mmsb1 - mmsb2)*pow2(-mmgl + mmsb2)) + (64*mmb*mmsb2*log(mmgl/mmu)*
     pow2(mmsb1))/(9.*pow2(mmsb1 - mmsb2)*pow3(-mmgl + mmsb1)) - (64*mmb*mmsb2*
     log(mmsb1/mmu)*pow2(mmsb1))/(9.*pow2(mmsb1 - mmsb2)*pow3(-mmgl + mmsb1)) -
     (64*mmb*mmsb2*log(mmgl/mmu)*log(mmsb2/mmu)*pow2(mmsb1))/(9.*pow2(mmsb1 -
     mmsb2)*pow3(-mmgl + mmsb1)) + (64*mmb*mmsb2*log(mmsb1/mmu)*log(mmsb2/mmu)*
     pow2(mmsb1))/(9.*pow2(mmsb1 - mmsb2)*pow3(-mmgl + mmsb1)) - (80*mmb*log(
     mmgl/mmu)*pow3(mmsb1))/(9.*pow2(mmsb1 - mmsb2)*pow3(-mmgl + mmsb1)) + (80*
     mmb*log(mmsb1/mmu)*pow3(mmsb1))/(9.*pow2(mmsb1 - mmsb2)*pow3(-mmgl +
     mmsb1)) + (32*mmb*pow2(log(mmgl/mmu))*pow3(mmsb1))/(9.*pow2(mmsb1 - mmsb2)
     *pow3(-mmgl + mmsb1)) - (32*mmb*pow2(log(mmsb1/mmu))*pow3(mmsb1))/(9.*
     pow2(mmsb1 - mmsb2)*pow3(-mmgl + mmsb1)) + (512*mmb*Fin20(mmsb1,mmgl,mmu))
     /(9.*pow3(mmsb1 - mmsb2)) - (616*mmb*mmsb2*Fin20(mmsb1,mmgl,mmu))/(9.*(-
     mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) + (104*mmb*mmsb2*Fin20(mmsb1,mmgl,mmu))
     /(9.*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) - (512*mmb*Fin20(mmsb2,mmgl,mmu)
     )/(9.*pow3(mmsb1 - mmsb2)) - (104*mmb*mmsb2*Fin20(mmsb2,mmgl,mmu))/(9.*(-
     mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) + (616*mmb*mmsb2*Fin20(mmsb2,mmgl,mmu))
     /(9.*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) - (512*mmb*mmsb2*log(mmsb1/mmu))
     /(9.*pow3(mmsb1 - mmsb2)) + (512*mmb*mmsb2*log(mmsb2/mmu))/(9.*pow3(mmsb1
     - mmsb2)) - (5056*mmb*pow2(mmsb2))/(9.*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)
     ) + (5056*mmb*pow2(mmsb2))/(9.*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) - (80*
     mmb*zt2*pow2(mmsb2))/((-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) + (80*mmb*zt2*
     pow2(mmsb2))/((-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) + (3232*mmb*log(mmgl/
     mmu)*pow2(mmsb2))/(9.*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) - (3232*mmb*
     log(mmgl/mmu)*pow2(mmsb2))/(9.*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) + (
     824*mmb*log(mmsb1/mmu)*pow2(mmsb2))/(9.*(-mmgl + mmsb1)*pow3(mmsb1 -
     mmsb2)) - (104*mmb*log(mmsb1/mmu)*pow2(mmsb2))/(3.*(-mmgl + mmsb2)*pow3(
     mmsb1 - mmsb2)) - (28*mmb*log(mmgl/mmu)*log(mmsb1/mmu)*pow2(mmsb2))/(3.*(-
     mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) + (28*mmb*log(mmgl/mmu)*log(mmsb1/mmu)*
     pow2(mmsb2))/(3.*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) + (88*mmb*log(mmsb2/
     mmu)*pow2(mmsb2))/(3.*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) - (776*mmb*log(
     mmsb2/mmu)*pow2(mmsb2))/(9.*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) - (28*
     mmb*log(mmgl/mmu)*log(mmsb2/mmu)*pow2(mmsb2))/(3.*(-mmgl + mmsb1)*pow3(
     mmsb1 - mmsb2)) + (28*mmb*log(mmgl/mmu)*log(mmsb2/mmu)*pow2(mmsb2))/(3.*(-
     mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) - (20*mmb*log(mmsb1/mmu)*log(mmsb2/mmu)
     *pow2(mmsb2))/(9.*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) + (20*mmb*log(
     mmsb1/mmu)*log(mmsb2/mmu)*pow2(mmsb2))/(9.*(-mmgl + mmsb2)*pow3(mmsb1 -
     mmsb2)) - (892*mmb*pow2(mmsb2)*pow2(log(mmgl/mmu)))/(9.*(-mmgl + mmsb1)*
     pow3(mmsb1 - mmsb2)) + (892*mmb*pow2(mmsb2)*pow2(log(mmgl/mmu)))/(9.*(-
     mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) + (256*mmb*mmsb2*pow2(log(mmsb1/mmu)))/
     (9.*pow3(mmsb1 - mmsb2)) - (308*mmb*pow2(mmsb2)*pow2(log(mmsb1/mmu)))/(9.*
     (-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) + (52*mmb*pow2(mmsb2)*pow2(log(mmsb1/
     mmu)))/(9.*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) - (256*mmb*mmsb2*pow2(log(
     mmsb2/mmu)))/(9.*pow3(mmsb1 - mmsb2)) - (52*mmb*pow2(mmsb2)*pow2(log(
     mmsb2/mmu)))/(9.*(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) + (308*mmb*pow2(
     mmsb2)*pow2(log(mmsb2/mmu)))/(9.*(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) - (
     16*mmb*log(mmgl/mmu)*pow3(mmsb2))/(9.*pow2(-mmgl + mmsb1)*pow3(mmsb1 -
     mmsb2)) + (16*mmb*log(mmsb1/mmu)*pow3(mmsb2))/(9.*pow2(-mmgl + mmsb1)*
     pow3(mmsb1 - mmsb2)) - (64*mmb*log(mmgl/mmu)*log(mmsb1/mmu)*pow3(mmsb2))/(
     9.*pow2(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) - (64*mmb*log(mmgl/mmu)*log(
     mmsb2/mmu)*pow3(mmsb2))/(9.*pow2(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) + (
     64*mmb*log(mmsb1/mmu)*log(mmsb2/mmu)*pow3(mmsb2))/(9.*pow2(-mmgl + mmsb1)*
     pow3(mmsb1 - mmsb2)) + (16*mmb*log(mmgl/mmu)*pow3(mmsb2))/(9.*pow2(-mmgl +
     mmsb2)*pow3(mmsb1 - mmsb2)) - (16*mmb*log(mmsb2/mmu)*pow3(mmsb2))/(9.*
     pow2(-mmgl + mmsb2)*pow3(mmsb1 - mmsb2)) + (64*mmb*pow2(log(mmgl/mmu))*
     pow3(mmsb2))/(9.*pow2(-mmgl + mmsb1)*pow3(mmsb1 - mmsb2)) + (64*mmb*mmsb1*
     log(mmgl/mmu)*pow2(mmsb2))/(9.*pow2(mmsb1 - mmsb2)*pow3(-mmgl + mmsb2)) -
     (64*mmb*mmsb1*log(mmgl/mmu)*log(mmsb1/mmu)*pow2(mmsb2))/(9.*pow2(mmsb1 -
     mmsb2)*pow3(-mmgl + mmsb2)) - (64*mmb*mmsb1*log(mmsb2/mmu)*pow2(mmsb2))/(
     9.*pow2(mmsb1 - mmsb2)*pow3(-mmgl + mmsb2)) + (64*mmb*mmsb1*log(mmsb1/mmu)
     *log(mmsb2/mmu)*pow2(mmsb2))/(9.*pow2(mmsb1 - mmsb2)*pow3(-mmgl + mmsb2))
     - (80*mmb*log(mmgl/mmu)*pow3(mmsb2))/(9.*pow2(mmsb1 - mmsb2)*pow3(-mmgl +
     mmsb2)) + (80*mmb*log(mmsb2/mmu)*pow3(mmsb2))/(9.*pow2(mmsb1 - mmsb2)*
     pow3(-mmgl + mmsb2)) + (32*mmb*pow2(log(mmgl/mmu))*pow3(mmsb2))/(9.*pow2(
     mmsb1 - mmsb2)*pow3(-mmgl + mmsb2)) - (32*mmb*pow2(log(mmsb2/mmu))*pow3(
     mmsb2))/(9.*pow2(mmsb1 - mmsb2)*pow3(-mmgl + mmsb2)) + (16*mmb*log(mmgl/
     mmu)*log(mmsb1/mmu)*pow4(mmsb1))/(9.*pow2(mmsb1 - mmsb2)*pow4(-mmgl +
     mmsb1)) - (8*mmb*pow2(log(mmgl/mmu))*pow4(mmsb1))/(9.*pow2(mmsb1 - mmsb2)*
     pow4(-mmgl + mmsb1)) - (8*mmb*pow2(log(mmsb1/mmu))*pow4(mmsb1))/(9.*pow2(
     mmsb1 - mmsb2)*pow4(-mmgl + mmsb1)) + (16*mmb*log(mmsb1/mmu)*pow3(mmsb2))/
     (9.*(-mmgl + mmsb1)*pow4(mmsb1 - mmsb2)) - (16*mmb*log(mmsb1/mmu)*pow3(
     mmsb2))/(9.*(-mmgl + mmsb2)*pow4(mmsb1 - mmsb2)) - (64*mmb*log(mmgl/mmu)*
     log(mmsb1/mmu)*pow3(mmsb2))/(9.*(-mmgl + mmsb1)*pow4(mmsb1 - mmsb2)) + (
     64*mmb*log(mmgl/mmu)*log(mmsb1/mmu)*pow3(mmsb2))/(9.*(-mmgl + mmsb2)*pow4(
     mmsb1 - mmsb2)) - (16*mmb*log(mmsb2/mmu)*pow3(mmsb2))/(9.*(-mmgl + mmsb1)*
     pow4(mmsb1 - mmsb2)) + (16*mmb*log(mmsb2/mmu)*pow3(mmsb2))/(9.*(-mmgl +
     mmsb2)*pow4(mmsb1 - mmsb2)) - (64*mmb*log(mmgl/mmu)*log(mmsb2/mmu)*pow3(
     mmsb2))/(9.*(-mmgl + mmsb1)*pow4(mmsb1 - mmsb2)) + (64*mmb*log(mmgl/mmu)*
     log(mmsb2/mmu)*pow3(mmsb2))/(9.*(-mmgl + mmsb2)*pow4(mmsb1 - mmsb2)) + (
     64*mmb*log(mmsb1/mmu)*log(mmsb2/mmu)*pow3(mmsb2))/(9.*(-mmgl + mmsb1)*
     pow4(mmsb1 - mmsb2)) - (64*mmb*log(mmsb1/mmu)*log(mmsb2/mmu)*pow3(mmsb2))/
     (9.*(-mmgl + mmsb2)*pow4(mmsb1 - mmsb2)) + (64*mmb*pow2(log(mmgl/mmu))*
     pow3(mmsb2))/(9.*(-mmgl + mmsb1)*pow4(mmsb1 - mmsb2)) - (64*mmb*pow2(log(
     mmgl/mmu))*pow3(mmsb2))/(9.*(-mmgl + mmsb2)*pow4(mmsb1 - mmsb2)) - (16*
     mmb*log(mmgl/mmu)*log(mmsb1/mmu)*pow4(mmsb2))/(9.*pow2(-mmgl + mmsb1)*
     pow4(mmsb1 - mmsb2)) - (16*mmb*log(mmgl/mmu)*log(mmsb2/mmu)*pow4(mmsb2))/(
     9.*pow2(-mmgl + mmsb1)*pow4(mmsb1 - mmsb2)) + (16*mmb*log(mmsb1/mmu)*log(
     mmsb2/mmu)*pow4(mmsb2))/(9.*pow2(-mmgl + mmsb1)*pow4(mmsb1 - mmsb2)) - (
     16*mmb*log(mmgl/mmu)*log(mmsb1/mmu)*pow4(mmsb2))/(9.*pow2(-mmgl + mmsb2)*
     pow4(mmsb1 - mmsb2)) - (16*mmb*log(mmgl/mmu)*log(mmsb2/mmu)*pow4(mmsb2))/(
     9.*pow2(-mmgl + mmsb2)*pow4(mmsb1 - mmsb2)) + (16*mmb*log(mmsb1/mmu)*log(
     mmsb2/mmu)*pow4(mmsb2))/(9.*pow2(-mmgl + mmsb2)*pow4(mmsb1 - mmsb2)) + (
     16*mmb*pow2(log(mmgl/mmu))*pow4(mmsb2))/(9.*pow2(-mmgl + mmsb1)*pow4(mmsb1
     - mmsb2)) + (16*mmb*pow2(log(mmgl/mmu))*pow4(mmsb2))/(9.*pow2(-mmgl +
     mmsb2)*pow4(mmsb1 - mmsb2)) + (16*mmb*log(mmgl/mmu)*log(mmsb2/mmu)*pow4(
     mmsb2))/(9.*pow2(mmsb1 - mmsb2)*pow4(-mmgl + mmsb2)) - (8*mmb*pow2(log(
     mmgl/mmu))*pow4(mmsb2))/(9.*pow2(mmsb1 - mmsb2)*pow4(-mmgl + mmsb2)) - (8*
     mmb*pow2(log(mmsb2/mmu))*pow4(mmsb2))/(9.*pow2(mmsb1 - mmsb2)*pow4(-mmgl +
     mmsb2)) - (32*mmb*log(mmgl/mmu)*log(mmsb1/mmu)*pow4(mmsb2))/(9.*(-mmgl +
     mmsb1)*pow5(mmsb1 - mmsb2)) + (32*mmb*log(mmgl/mmu)*log(mmsb1/mmu)*pow4(
     mmsb2))/(9.*(-mmgl + mmsb2)*pow5(mmsb1 - mmsb2)) - (32*mmb*log(mmgl/mmu)*
     log(mmsb2/mmu)*pow4(mmsb2))/(9.*(-mmgl + mmsb1)*pow5(mmsb1 - mmsb2)) + (
     32*mmb*log(mmgl/mmu)*log(mmsb2/mmu)*pow4(mmsb2))/(9.*(-mmgl + mmsb2)*pow5(
     mmsb1 - mmsb2)) + (32*mmb*log(mmsb1/mmu)*log(mmsb2/mmu)*pow4(mmsb2))/(9.*(
     -mmgl + mmsb1)*pow5(mmsb1 - mmsb2)) - (32*mmb*log(mmsb1/mmu)*log(mmsb2/
     mmu)*pow4(mmsb2))/(9.*(-mmgl + mmsb2)*pow5(mmsb1 - mmsb2)) + (32*mmb*pow2(
     log(mmgl/mmu))*pow4(mmsb2))/(9.*(-mmgl + mmsb1)*pow5(mmsb1 - mmsb2)) - (
     32*mmb*pow2(log(mmgl/mmu))*pow4(mmsb2))/(9.*(-mmgl + mmsb2)*pow5(mmsb1 -
     mmsb2))));

   return result * twoLoop;
}

} // namespace mssm_twoloop_mb
} // namespace flexiblesusy
