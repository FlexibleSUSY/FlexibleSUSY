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

// This file has been generated at Sun 26 Feb 2023 12:30:46
// with the script "tquark_to_cpp.m".

#include "mssm_twoloop_mt.hpp"
#include "Li2.hpp"
#include <cmath>
#include <limits>

namespace flexiblesusy {
namespace mssm_twoloop_mt {

namespace {
   constexpr double Pi      = 3.1415926535897932384626433832795;
   constexpr double zt2     = 1.6449340668482264364724151666460;    // Zeta[2]
   constexpr double zt3     = 1.2020569031595942853997381615114;    // Zeta[3]
   constexpr double log2    = 6.9314718055994530941723212145818e-1; // Log[2]
   constexpr double oneLoop = 6.3325739776461107152424664506080e-3; // 1/(4Pi)^2
   constexpr double twoLoop = 4.0101493182360684332628059637182e-5; // 1/(4Pi)^4

   constexpr double pow2(double x) noexcept { return x*x; }
   constexpr double pow3(double x) noexcept { return x*x*x; }
   constexpr double pow4(double x) noexcept { return pow2(pow2(x)); }
   constexpr double pow5(double x) noexcept { return x*pow4(x); }
   constexpr double pow6(double x) noexcept { return pow2(pow3(x)); }
   constexpr double pow7(double x) noexcept { return x*pow6(x); }
   constexpr double pow8(double x) noexcept { return pow2(pow4(x)); }

   bool is_zero(double a, double prec) noexcept
   {
      return std::abs(a) < prec;
   }

   bool is_equal(double a, double b, double prec) noexcept
   {
      return is_zero(a - b, prec);
   }

   bool is_equal_rel(double a, double b, double prec) noexcept
   {
      if (is_equal(a, b, std::numeric_limits<double>::epsilon()))
         return true;

      const double min = std::min(std::abs(a), std::abs(b));

      if (min < std::numeric_limits<double>::epsilon())
         return is_equal(a, b, prec);

      const double max = std::max(std::abs(a), std::abs(b));

      return is_equal(a, b, prec*max);
   }

   /**
    * fin[] function from arXiv:hep-ph/0507139 .
    *
    * @param mm1 squared mass \f$m_1^2\f$
    * @param mm2 squared mass \f$m_2^2\f$
    * @param mmu squared renormalization scale
    *
    * @return fin(m12, m22)
    */
   double fin(double mm1, double mm2, double mmu) noexcept
   {
      const double log1u = std::log(mm1/mmu);
      const double log2u = std::log(mm2/mmu);
      const double log12 = std::log(mm1/mm2);

      return (6*(mm1*log1u + mm2*log2u) +
         (-mm1 - mm2)*(7 + pow2(Pi)/6.) +
         (mm1 - mm2)*(2*Li2(1 - mm1/mm2) + pow2(log12)/2.) +
         ((mm1 + mm2)*pow2(log12))/2. -
         2*(mm1*pow2(log1u) + mm2*pow2(log2u)))/2.;
   }

   /// shift gluino mass away from mst1 and mst2 if too close
   double shift_mg(double mg, double mst1, double mst2) noexcept
   {
      if (is_equal_rel(std::min(mst1, mst2), mg, 0.0003))
         return mg * 0.9995;

      if (is_equal_rel(std::max(mst1, mst2), mg, 0.0003))
         return mg * 1.0005;

      return mg;
   }

} // anonymous namespace

/// 1-loop QCD contributions to Delta Mt over mt [hep-ph/0210258]
double dMt_over_mt_1loop_qcd(const Parameters& pars)
{
   const double g32 = pow2(pars.g3);
   const double mmt = pow2(pars.mt);
   const double mmu = pow2(pars.Q);
   const double logmmtmmu = std::log(mmt/mmu);

   const double result = 6.666666666666667 - 4*logmmtmmu;

   return result * g32 * oneLoop;
}

/// 1-loop SUSY contributions to Delta Mt over mt [hep-ph/0210258]
double dMt_over_mt_1loop_susy(const Parameters& pars)
{
   const double g32    = pow2(pars.g3);
   const double mt     = pars.mt;
   const double mgl    = pars.mg;
   const double mmu    = pow2(pars.Q);
   const double mmgl   = pow2(pars.mg);
   const double mmst1  = pow2(pars.mst1);
   const double mmst2  = pow2(pars.mst2);
   const double SX     = 2*mt*pars.xt;
   const double s2t    = SX / (mmst1 - mmst2);

   if (is_equal(mmst1, mmst2, 1e-6) && is_equal(mmst1, mmgl, 1e-6)) {
      const double logmmglmmu = std::log(mmgl/mmu);
      const double result = (4*logmmglmmu)/3.;

      return result * g32 * oneLoop;
   }

   if (is_equal(mmst1, mmst2, 1e-6)) {
      const double logmmglmmu  = std::log(mmgl/mmu);
      const double logmmst2mmu = std::log(mmst2/mmu);

      const double result =
      (4*logmmglmmu*pow2(mmgl) - 2*(-4*mmgl*mmst2 + 2*logmmst2mmu*(2*mmgl -
         mmst2)*mmst2 + 3*pow2(mmgl) + pow2(mmst2)))/(3.*pow2(mmgl - mmst2));

      return result * g32 * oneLoop;
   }

   if (is_equal(mmgl, mmst1, 1e-6)) {
      const double logmmglmmu  = std::log(mmgl/mmu);
      const double logmmst2mmu = std::log(mmst2/mmu);

      const double result =
      (2*logmmst2mmu*mmst2*(-2*mmgl*(mt + mgl*s2t) + mmst2*(mt + 2*mgl*s2t)) -
         (mmgl - mmst2)*(-(mmst2*(mt + 4*mgl*s2t)) + mmgl*(3*mt + 4*mgl*s2t)) +
         2*logmmglmmu*(-2*mmgl*mmst2*(mt - mgl*s2t) + 2*mt*pow2(mmgl) + (mt - 2
         *mgl*s2t)*pow2(mmst2)))/(3.*mt*pow2(mmgl - mmst2));

      return result * g32 * oneLoop;
   }

   if (is_equal(mmgl, mmst2, 1e-6)) {
      const double logmmglmmu = std::log(mmgl/mmu);
      const double logmmst1mmu = std::log(mmst1/mmu);

      const double result =
      (2*logmmst1mmu*mmst1*(-2*mmgl*mt + mmst1*mt + 2*mgl*mmgl*s2t - 2*mgl*
         mmst1*s2t) - (mmgl - mmst1)*(3*mmgl*mt - mmst1*mt - 4*mgl*mmgl*s2t + 4
         *mgl*mmst1*s2t) + 2*logmmglmmu*(-2*mmgl*mmst1*(mt + mgl*s2t) + 2*mt*
         pow2(mmgl) + (mt + 2*mgl*s2t)*pow2(mmst1)))/(3.*mt*pow2(mmgl - mmst1))
         ;
      return result * g32 * oneLoop;

   }

   const double logmmglmmu  = std::log(mmgl/mmu);
   const double logmmst1mmu = std::log(mmst1/mmu);
   const double logmmst2mmu = std::log(mmst2/mmu);

   const double result =
   (2*(-3 + mmst1/(-mmgl + mmst1) + mmst2/(-mmgl + mmst2) - (2*logmmglmmu*mgl*
      mmgl*(mmst1 - mmst2)*s2t)/((mmgl - mmst1)*(mmgl - mmst2)*mt) + (2*
      logmmst1mmu*mgl*mmst1*s2t)/(mmgl*mt - mmst1*mt) - (2*logmmst2mmu*mgl*
      mmst2*s2t)/(mmgl*mt - mmst2*mt) + (logmmst1mmu*mmst1*(-2*mmgl + mmst1))/
      pow2(mmgl - mmst1) + (logmmst2mmu*mmst2*(-2*mmgl + mmst2))/pow2(mmgl -
      mmst2) + (logmmglmmu*pow2(mmgl)*(-2*mmgl*(mmst1 + mmst2) + 2*pow2(mmgl) +
      pow2(mmst1) + pow2(mmst2)))/(pow2(mmgl - mmst1)*pow2(mmgl - mmst2))))/3.;

   return result * g32 * oneLoop;
}

/// 1-loop full SQCD contributions to Delta Mt over mt [hep-ph/0210258]
double dMt_over_mt_1loop(const Parameters& pars)
{
   return dMt_over_mt_1loop_qcd(pars) + dMt_over_mt_1loop_susy(pars);
}

/// 2-loop QCD contributions to Delta Mt over mt [hep-ph/0507139]
double dMt_over_mt_2loop_qcd(const Parameters& pars)
{
   const double g34 = pow4(pars.g3);
   const double mmt = pow2(pars.mt);
   const double mmu = pow2(pars.Q);
   const double logmmtmmu = std::log(mmt/mmu);

   const double result =
   -82*logmmtmmu + (2011 + 96*(1 + 2*log2)*zt2 - 48*zt3)/18. + 22*pow2(
      logmmtmmu);

   return result * g34 * twoLoop;
}

/// 2-loop SUSY contributions to Delta Mt over mt [hep-ph/0507139]
double dMt_over_mt_2loop_susy(const Parameters& pars)
{
   const double g34    = pow4(pars.g3);
   const double mt     = pars.mt;
   const double mgl    = shift_mg(pars.mg, pars.mst1, pars.mst2);
   const double mmu    = pow2(pars.Q);
   const double mmt    = pow2(pars.mt);
   const double mmgl   = pow2(mgl);
   const double mmst1  = pow2(pars.mst1);
   const double mmst2  = pow2(pars.mst2);
   const double mmsusy = pow2(pars.msusy);
   const double SX     = 2*mt*pars.xt;
   const double s2t    = SX / (mmst1 - mmst2);

   if (is_equal(mmst1, mmst2, mmt) && is_equal(mmst1, pow2(pars.mg), mmt) &&
       is_equal(mmst1, mmsusy, mmt)) {
      const double logmmsusymmu = std::log(mmsusy/mmu);
      const double logmmtmmu    = std::log(mmt/mmu);

      const double result =
      (1745 - 4*logmmsusymmu*(-677 + 288*logmmtmmu) + 372*pow2(logmmsusymmu))/
         54.;

      return result * g34 * twoLoop;
   }

   const double logmmsusymmu = std::log(mmsusy/mmu);
   const double logmmtmmu    = std::log(mmt/mmu);
   const double logmmst1mmu  = std::log(mmst1/mmu);
   const double logmmst2mmu  = std::log(mmst2/mmu);
   const double logmmglmmu   = std::log(mmgl/mmu );

   const double result =
   -10.055555555555555 + 82*logmmtmmu - (20*logmmsusymmu*logmmtmmu)/3. - (56*
      mmst1)/(-mmgl + mmst1) - (896*mmst1)/(9.*(mmst1 - mmst2)) + (6*mmst2)/(-
      mmgl + mmst1) + (896*mmst2)/(9.*(mmst1 - mmst2)) - (2*mmst1)/(3.*(-mmgl +
      mmst2)) - (188*mmst2)/(3.*(-mmgl + mmst2)) + logmmtmmu*(-74 - (8*mmst1)/(
      3.*(-mmgl + mmst1)) - (8*mmst2)/(3.*(-mmgl + mmst2))) + (60*mmsusy)/(-
      mmgl + mmst1) + (60*mmsusy)/(-mmgl + mmst2) + (16*logmmst1mmu*logmmtmmu*
      mgl*mmst1*s2t)/(3.*(-mmgl + mmst1)*mt) - (16*logmmst2mmu*logmmtmmu*mgl*
      mmst2*s2t)/(3.*(-mmgl + mmst2)*mt) + (logmmglmmu*logmmtmmu*((-16*mgl*
      mmst1*s2t)/(3.*(-mmgl + mmst1)) + (16*mgl*mmst2*s2t)/(3.*(-mmgl + mmst2))
      ))/mt + (logmmsusymmu*((160*mgl*mmsusy*s2t)/(3.*(-mmgl + mmst1)) - (160*
      mgl*mmsusy*s2t)/(3.*(-mmgl + mmst2))))/mt + (8*zt2)/9. - (91*mmst1*zt2)/(
      9.*(-mmgl + mmst1)) - (128*mmst1*zt2)/(9.*(mmst1 - mmst2)) + (mmst2*zt2)/
      (-mmgl + mmst1) + (128*mmst2*zt2)/(9.*(mmst1 - mmst2)) + (mmst1*zt2)/(9.*
      (-mmgl + mmst2)) - (11*mmst2*zt2)/(-mmgl + mmst2) + (10*mmsusy*zt2)/(-
      mmgl + mmst1) + (10*mmsusy*zt2)/(-mmgl + mmst2) + (((-20*mgl*mmsusy*s2t)/
      (3.*(-mmgl + mmst1)) + (20*mgl*mmsusy*s2t)/(3.*(-mmgl + mmst2)))*pow2(
      logmmsusymmu))/mt + (896*pow2(mmst1))/(9.*(-mmgl + mmst1)*(mmst1 - mmst2)
      ) + (20*pow2(mmst1))/(3.*(-mmgl + mmst1)*(-mmgl + mmst2)) + (128*zt2*pow2
      (mmst1))/(9.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (8*zt2*pow2(mmst1))/(9.*(
      -mmgl + mmst1)*(-mmgl + mmst2)) + (fin(mmst2,mmgl,mmu)*((-16*mgl*s2t)/(9.
      *(-mmgl + mmst1)) + (88*mgl*s2t)/(9.*(-mmgl + mmst2)) - (4*mgl*mmst1*s2t)
      /(3.*pow2(-mmgl + mmst1)) + (4*mgl*mmst2*s2t)/(3.*pow2(-mmgl + mmst1))))/
      mt + (fin(mmst1,mmsusy,mmu)*((40*mgl*mmst1*s2t)/(3.*pow2(-mmgl + mmst1))
      - (40*mgl*mmsusy*s2t)/(3.*pow2(-mmgl + mmst1))))/mt + logmmst1mmu*
      logmmtmmu*(-0.6666666666666666 - (16*mmst1)/(3.*(-mmgl + mmst1)) + (8*
      pow2(mmst1))/(3.*pow2(-mmgl + mmst1))) - (2825*mmst1*mmst2)/(324.*pow2(-
      mmgl + mmst1)) - (27275*mmst1*mmsusy)/(324.*pow2(-mmgl + mmst1)) - (40*
      logmmst1mmu*logmmsusymmu*mgl*mmst1*mmsusy*s2t)/(3.*mt*pow2(-mmgl + mmst1)
      ) - (4*mmst1*mmst2*zt2)/(3.*pow2(-mmgl + mmst1)) - (130*mmst1*mmsusy*zt2)
      /(9.*pow2(-mmgl + mmst1)) + (886*pow2(mmst1))/(9.*pow2(-mmgl + mmst1)) +
      (233*mmst2*pow2(mmst1))/(324.*(mmst1 - mmst2)*pow2(-mmgl + mmst1)) + (
      1355*mmsusy*pow2(mmst1))/(324.*(mmst1 - mmsusy)*pow2(-mmgl + mmst1)) + (
      44*zt2*pow2(mmst1))/(3.*pow2(-mmgl + mmst1)) + (10*mmsusy*zt2*pow2(mmst1)
      )/(9.*(mmst1 - mmsusy)*pow2(-mmgl + mmst1)) - (896*pow2(mmst2))/(9.*(
      mmst1 - mmst2)*(-mmgl + mmst2)) - (128*zt2*pow2(mmst2))/(9.*(mmst1 -
      mmst2)*(-mmgl + mmst2)) - (233*mmst1*pow2(mmst2))/(216.*(mmst1 - mmst2)*
      pow2(-mmgl + mmst1)) + (233*pow2(mmst1)*pow2(mmst2))/(648.*pow2(-mmgl +
      mmst1)*pow2(mmst1 - mmst2)) + (fin(mmst1,mmst2,mmu)*((4*mgl*s2t)/(9.*(-
      mmgl + mmst1)) - (4*mgl*s2t)/(9.*(-mmgl + mmst2)) + (4*mgl*mmst1*s2t)/(3.
      *pow2(-mmgl + mmst1)) - (4*mgl*mmst2*s2t)/(3.*pow2(-mmgl + mmst1)) + (4*
      mgl*mmst1*s2t)/(3.*pow2(-mmgl + mmst2)) - (4*mgl*mmst2*s2t)/(3.*pow2(-
      mmgl + mmst2))))/mt + (fin(mmst1,mmgl,mmu)*((-88*mgl*s2t)/(9.*(-mmgl +
      mmst1)) + (16*mgl*s2t)/(9.*(-mmgl + mmst2)) - (4*mgl*mmst1*s2t)/(3.*pow2(
      -mmgl + mmst2)) + (4*mgl*mmst2*s2t)/(3.*pow2(-mmgl + mmst2))))/mt + (fin(
      mmgl,mmsusy,mmu)*((-40*mgl*s2t)/(3.*(-mmgl + mmst1)) + (40*mgl*s2t)/(3.*(
      -mmgl + mmst2)) - (40*mgl*mmst1*s2t)/(3.*pow2(-mmgl + mmst1)) + (40*mgl*
      mmsusy*s2t)/(3.*pow2(-mmgl + mmst1)) + (40*mgl*mmst2*s2t)/(3.*pow2(-mmgl
      + mmst2)) - (40*mgl*mmsusy*s2t)/(3.*pow2(-mmgl + mmst2))))/mt + (fin(
      mmst2,mmsusy,mmu)*((-40*mgl*mmst2*s2t)/(3.*pow2(-mmgl + mmst2)) + (40*mgl
      *mmsusy*s2t)/(3.*pow2(-mmgl + mmst2))))/mt + (logmmglmmu*logmmsusymmu*((-
      40*mgl*mmsusy*s2t)/(3.*(-mmgl + mmst1)) + (40*mgl*mmsusy*s2t)/(3.*(-mmgl
      + mmst2)) + (40*mgl*mmst1*mmsusy*s2t)/(3.*pow2(-mmgl + mmst1)) - (40*mgl*
      mmst2*mmsusy*s2t)/(3.*pow2(-mmgl + mmst2))))/mt + logmmglmmu*logmmtmmu*(-
      13.333333333333334 + (16*mmst1)/(3.*(-mmgl + mmst1)) + (16*mmst2)/(3.*(-
      mmgl + mmst2)) - (8*pow2(mmst1))/(3.*pow2(-mmgl + mmst1)) - (8*pow2(mmst2
      ))/(3.*pow2(-mmgl + mmst2))) + logmmst2mmu*logmmtmmu*(-0.6666666666666666
       - (16*mmst2)/(3.*(-mmgl + mmst2)) + (8*pow2(mmst2))/(3.*pow2(-mmgl +
      mmst2))) - (8*mmst1*mmst2)/pow2(-mmgl + mmst2) - (27275*mmst2*mmsusy)/(
      324.*pow2(-mmgl + mmst2)) + (40*logmmst2mmu*logmmsusymmu*mgl*mmst2*mmsusy
      *s2t)/(3.*mt*pow2(-mmgl + mmst2)) - (4*mmst1*mmst2*zt2)/(3.*pow2(-mmgl +
      mmst2)) - (130*mmst2*mmsusy*zt2)/(9.*pow2(-mmgl + mmst2)) + (886*pow2(
      mmst2))/(9.*pow2(-mmgl + mmst2)) + (1355*mmsusy*pow2(mmst2))/(324.*(mmst2
       - mmsusy)*pow2(-mmgl + mmst2)) + (44*zt2*pow2(mmst2))/(3.*pow2(-mmgl +
      mmst2)) + (10*mmsusy*zt2*pow2(mmst2))/(9.*(mmst2 - mmsusy)*pow2(-mmgl +
      mmst2)) - (515*mmst1*pow2(mmsusy))/(108.*(mmst1 - mmsusy)*pow2(-mmgl +
      mmst1)) - (20*mmst1*zt2*pow2(mmsusy))/(9.*(mmst1 - mmsusy)*pow2(-mmgl +
      mmst1)) - (515*mmst2*pow2(mmsusy))/(108.*(mmst2 - mmsusy)*pow2(-mmgl +
      mmst2)) - (20*mmst2*zt2*pow2(mmsusy))/(9.*(mmst2 - mmsusy)*pow2(-mmgl +
      mmst2)) + (95*pow2(mmst1)*pow2(mmsusy))/(162.*pow2(-mmgl + mmst1)*pow2(
      mmst1 - mmsusy)) + (10*zt2*pow2(mmst1)*pow2(mmsusy))/(9.*pow2(-mmgl +
      mmst1)*pow2(mmst1 - mmsusy)) + (95*pow2(mmst2)*pow2(mmsusy))/(162.*pow2(-
      mmgl + mmst2)*pow2(mmst2 - mmsusy)) + (10*zt2*pow2(mmst2)*pow2(mmsusy))/(
      9.*pow2(-mmgl + mmst2)*pow2(mmst2 - mmsusy)) - (256*pow2(s2t))/9. + (
      63625*mmst1*pow2(s2t))/(1296.*(-mmgl + mmst1)) + (896*mmst1*pow2(s2t))/(
      9.*(mmst1 - mmst2)) - (10535*mmst2*pow2(s2t))/(3888.*(-mmgl + mmst1)) - (
      896*mmst2*pow2(s2t))/(9.*(mmst1 - mmst2)) + (mmst1*mmst2*pow2(s2t))/(648.
      *(-mmgl + mmst1)*(mmst1 - mmst2)) - (56281*mmst1*pow2(s2t))/(1296.*(-mmgl
       + mmst2)) + (10967*mmst2*pow2(s2t))/(1296.*(-mmgl + mmst2)) - (85*mmst1*
      mmst2*pow2(s2t))/(972.*(mmst1 - mmst2)*(-mmgl + mmst2)) + (41*mmst1*zt2*
      pow2(s2t))/(9.*(-mmgl + mmst1)) + (128*mmst1*zt2*pow2(s2t))/(9.*(mmst1 -
      mmst2)) - (128*mmst2*zt2*pow2(s2t))/(9.*(mmst1 - mmst2)) + (mmst1*mmst2*
      zt2*pow2(s2t))/(27.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (52*mmst1*zt2*pow2
      (s2t))/(9.*(-mmgl + mmst2)) - (11*mmst2*zt2*pow2(s2t))/(9.*(-mmgl + mmst2
      )) - (mmst1*mmst2*zt2*pow2(s2t))/(27.*(mmst1 - mmst2)*(-mmgl + mmst2)) -
      (96683*pow2(mmst1)*pow2(s2t))/(972.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (
      52825*pow2(mmst1)*pow2(s2t))/(1296.*(-mmgl + mmst1)*(-mmgl + mmst2)) - (
      85*pow2(mmst1)*pow2(s2t))/(972.*(mmst1 - mmst2)*(-mmgl + mmst2)) - (383*
      zt2*pow2(mmst1)*pow2(s2t))/(27.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (52*
      zt2*pow2(mmst1)*pow2(s2t))/(9.*(-mmgl + mmst1)*(-mmgl + mmst2)) - (zt2*
      pow2(mmst1)*pow2(s2t))/(27.*(mmst1 - mmst2)*(-mmgl + mmst2)) + (16*mmst1*
      mmst2*pow2(s2t))/(9.*pow2(-mmgl + mmst1)) - (2*pow2(mmst1)*pow2(s2t))/
      pow2(-mmgl + mmst1) + (85*pow2(mmst2)*pow2(s2t))/(972.*(-mmgl + mmst1)*(
      mmst1 - mmst2)) + (387233*pow2(mmst2)*pow2(s2t))/(3888.*(mmst1 - mmst2)*(
      -mmgl + mmst2)) + (zt2*pow2(mmst2)*pow2(s2t))/(27.*(-mmgl + mmst1)*(mmst1
       - mmst2)) + (383*zt2*pow2(mmst2)*pow2(s2t))/(27.*(mmst1 - mmst2)*(-mmgl
      + mmst2)) - (167*mmst1*pow2(mmst2)*pow2(s2t))/(3888.*(-mmgl + mmst1)*pow2
      (mmst1 - mmst2)) + (16*mmst1*mmst2*pow2(s2t))/(9.*pow2(-mmgl + mmst2)) -
      (2*pow2(mmst2)*pow2(s2t))/pow2(-mmgl + mmst2) + (85*pow2(s2t)*pow3(mmst1)
      )/(972.*(-mmgl + mmst1)*(mmst1 - mmst2)*(-mmgl + mmst2)) + (zt2*pow2(s2t)
      *pow3(mmst1))/(27.*(-mmgl + mmst1)*(mmst1 - mmst2)*(-mmgl + mmst2)) + fin
      (mmst1,mmsusy,mmu)*(20/(3.*(-mmgl + mmst1)) - (220*mmst1)/(9.*pow2(-mmgl
      + mmst1)) + (10*mmsusy)/pow2(-mmgl + mmst1) - (20*mmst1*mmsusy)/(9.*(
      mmst1 - mmsusy)*pow2(-mmgl + mmst1)) + (10*pow2(mmst1))/(9.*(mmst1 -
      mmsusy)*pow2(-mmgl + mmst1)) + (10*mmsusy*pow2(mmst1))/(9.*pow2(-mmgl +
      mmst1)*pow2(mmst1 - mmsusy)) - (10*mmst1*pow2(mmsusy))/(9.*pow2(-mmgl +
      mmst1)*pow2(mmst1 - mmsusy)) - (40*mmst1*mmsusy)/(3.*pow3(-mmgl + mmst1))
      + (40*pow2(mmst1))/(3.*pow3(-mmgl + mmst1))) - (112*pow3(mmst1))/(3.*pow3
      (-mmgl + mmst1)) - (16*zt2*pow3(mmst1))/(3.*pow3(-mmgl + mmst1)) + (271*
      pow3(mmst2))/(648.*(mmst1 - mmst2)*pow2(-mmgl + mmst1)) + (zt2*pow3(mmst2
      ))/(9.*(mmst1 - mmst2)*pow2(-mmgl + mmst1)) - (7*mmst1*pow3(mmst2))/(9.*
      pow2(-mmgl + mmst1)*pow2(mmst1 - mmst2)) - (mmst1*zt2*pow3(mmst2))/(9.*
      pow2(-mmgl + mmst1)*pow2(mmst1 - mmst2)) + (169*pow2(s2t)*pow3(mmst2))/(
      1296.*(-mmgl + mmst1)*pow2(mmst1 - mmst2)) - (85*pow2(s2t)*pow3(mmst2))/(
      972.*(-mmgl + mmst2)*pow2(mmst1 - mmst2)) + (zt2*pow2(s2t)*pow3(mmst2))/(
      27.*(-mmgl + mmst1)*pow2(mmst1 - mmst2)) - (zt2*pow2(s2t)*pow3(mmst2))/(
      27.*(-mmgl + mmst2)*pow2(mmst1 - mmst2)) + fin(mmgl,mmsusy,mmu)*(10/(3.*(
      -mmgl + mmst1)) + 10/(3.*(-mmgl + mmst2)) + (10*mmst1)/pow2(-mmgl + mmst1
      ) - (10*mmsusy)/pow2(-mmgl + mmst1) + (10*mmst2)/pow2(-mmgl + mmst2) - (
      10*mmsusy)/pow2(-mmgl + mmst2) + (40*mmst1*mmsusy)/(3.*pow3(-mmgl + mmst1
      )) - (40*pow2(mmst1))/(3.*pow3(-mmgl + mmst1)) + (40*mmst2*mmsusy)/(3.*
      pow3(-mmgl + mmst2)) - (40*pow2(mmst2))/(3.*pow3(-mmgl + mmst2))) + fin(
      mmst2,mmgl,mmu)*(7/(9.*(-mmgl + mmst1)) + 128/(9.*(mmst1 - mmst2)) - 5/(-
      mmgl + mmst2) + (8*mmst1)/(9.*(-mmgl + mmst1)*(-mmgl + mmst2)) - (128*
      mmst2)/(9.*(mmst1 - mmst2)*(-mmgl + mmst2)) + mmst1/pow2(-mmgl + mmst1) -
      mmst2/pow2(-mmgl + mmst1) + (12*mmst2)/pow2(-mmgl + mmst2) + (11*pow2(s2t
      ))/(9.*(-mmgl + mmst1)) - (128*pow2(s2t))/(9.*(mmst1 - mmst2)) - (11*pow2
      (s2t))/(9.*(-mmgl + mmst2)) + (26*mmst1*pow2(s2t))/(9.*(-mmgl + mmst1)*(-
      mmgl + mmst2)) + (128*mmst2*pow2(s2t))/(9.*(mmst1 - mmst2)*(-mmgl + mmst2
      )) - (26*mmst1*pow2(s2t))/(9.*pow2(-mmgl + mmst1)) + (4*mmst1*mmst2)/(3.*
      pow3(-mmgl + mmst1)) - (4*pow2(mmst1))/(3.*pow3(-mmgl + mmst1)) - (16*
      pow2(mmst2))/(3.*pow3(-mmgl + mmst2))) + fin(mmst1,mmgl,mmu)*(-37/(9.*(-
      mmgl + mmst1)) - 128/(9.*(mmst1 - mmst2)) + (128*mmst1)/(9.*(-mmgl +
      mmst1)*(mmst1 - mmst2)) - 1/(9.*(-mmgl + mmst2)) + (8*mmst1)/(9.*(-mmgl +
      mmst1)*(-mmgl + mmst2)) + (12*mmst1)/pow2(-mmgl + mmst1) - mmst1/pow2(-
      mmgl + mmst2) + mmst2/pow2(-mmgl + mmst2) + (5*pow2(s2t))/(3.*(-mmgl +
      mmst1)) + (128*pow2(s2t))/(9.*(mmst1 - mmst2)) - (128*mmst1*pow2(s2t))/(
      9.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (5*pow2(s2t))/(3.*(-mmgl + mmst2))
      + (26*mmst1*pow2(s2t))/(9.*(-mmgl + mmst1)*(-mmgl + mmst2)) - (26*mmst2*
      pow2(s2t))/(9.*pow2(-mmgl + mmst2)) - (16*pow2(mmst1))/(3.*pow3(-mmgl +
      mmst1)) + (4*mmst1*mmst2)/(3.*pow3(-mmgl + mmst2)) - (4*pow2(mmst2))/(3.*
      pow3(-mmgl + mmst2))) + fin(mmst1,mmst2,mmu)*(2/(9.*(-mmgl + mmst1)) + 10
      /(9.*(-mmgl + mmst2)) - (8*mmst1)/(9.*(-mmgl + mmst1)*(-mmgl + mmst2)) -
      (7*mmst1)/(3.*pow2(-mmgl + mmst1)) + mmst2/pow2(-mmgl + mmst1) + pow2(
      mmst2)/(9.*(mmst1 - mmst2)*pow2(-mmgl + mmst1)) - (mmst1*pow2(mmst2))/(9.
      *pow2(-mmgl + mmst1)*pow2(mmst1 - mmst2)) + mmst1/pow2(-mmgl + mmst2) - (
      7*mmst2)/(3.*pow2(-mmgl + mmst2)) - (67*pow2(s2t))/(54.*(-mmgl + mmst1))
      - (mmst2*pow2(s2t))/(54.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (65*pow2(s2t)
      )/(54.*(-mmgl + mmst2)) - (mmst1*pow2(s2t))/(54.*(-mmgl + mmst1)*(-mmgl +
      mmst2)) + (mmst2*pow2(s2t))/(54.*(mmst1 - mmst2)*(-mmgl + mmst2)) + (26*
      mmst1*pow2(s2t))/(9.*pow2(-mmgl + mmst1)) + (26*mmst2*pow2(s2t))/(9.*pow2
      (-mmgl + mmst2)) - (4*mmst1*mmst2)/(3.*pow3(-mmgl + mmst1)) + (4*pow2(
      mmst1))/(3.*pow3(-mmgl + mmst1)) + pow3(mmst2)/(9.*pow2(-mmgl + mmst1)*
      pow2(mmst1 - mmst2)) - (4*mmst1*mmst2)/(3.*pow3(-mmgl + mmst2)) + (4*pow2
      (mmst2))/(3.*pow3(-mmgl + mmst2))) + fin(mmst2,mmsusy,mmu)*(20/(3.*(-mmgl
       + mmst2)) - (220*mmst2)/(9.*pow2(-mmgl + mmst2)) + (10*mmsusy)/pow2(-
      mmgl + mmst2) - (20*mmst2*mmsusy)/(9.*(mmst2 - mmsusy)*pow2(-mmgl + mmst2
      )) + (10*pow2(mmst2))/(9.*(mmst2 - mmsusy)*pow2(-mmgl + mmst2)) + (10*
      mmsusy*pow2(mmst2))/(9.*pow2(-mmgl + mmst2)*pow2(mmst2 - mmsusy)) - (10*
      mmst2*pow2(mmsusy))/(9.*pow2(-mmgl + mmst2)*pow2(mmst2 - mmsusy)) - (40*
      mmst2*mmsusy)/(3.*pow3(-mmgl + mmst2)) + (40*pow2(mmst2))/(3.*pow3(-mmgl
      + mmst2))) + logmmglmmu*logmmsusymmu*((10*mmsusy)/(-mmgl + mmst1) + (10*
      mmsusy)/(-mmgl + mmst2) - (70*mmst1*mmsusy)/(3.*pow2(-mmgl + mmst1)) - (
      70*mmst2*mmsusy)/(3.*pow2(-mmgl + mmst2)) + (40*mmsusy*pow2(mmst1))/(3.*
      pow3(-mmgl + mmst1)) + (40*mmsusy*pow2(mmst2))/(3.*pow3(-mmgl + mmst2)))
      + (pow2(logmmglmmu)*((-10*mgl*mmst1*s2t)/(-mmgl + mmst1) - (2*mgl*mmst2*
      s2t)/(3.*(-mmgl + mmst1)) - (2*mgl*mmst1*s2t)/(9.*(-mmgl + mmst2)) + (98*
      mgl*mmst2*s2t)/(9.*(-mmgl + mmst2)) - (20*mgl*mmsusy*s2t)/(3.*(-mmgl +
      mmst1)) + (20*mgl*mmsusy*s2t)/(3.*(-mmgl + mmst2)) + (8*mgl*s2t*pow2(
      mmst1))/(3.*(-mmgl + mmst1)*(-mmgl + mmst2)) + (2*mgl*mmst1*mmst2*s2t)/(
      3.*pow2(-mmgl + mmst1)) + (20*mgl*mmst1*mmsusy*s2t)/(3.*pow2(-mmgl +
      mmst1)) - (218*mgl*s2t*pow2(mmst1))/(9.*pow2(-mmgl + mmst1)) - (14*mgl*
      mmst1*mmst2*s2t)/(9.*pow2(-mmgl + mmst2)) - (20*mgl*mmst2*mmsusy*s2t)/(3.
      *pow2(-mmgl + mmst2)) - (8*mgl*s2t*pow2(mmst1))/(9.*pow2(-mmgl + mmst2))
      + (70*mgl*s2t*pow2(mmst2))/(3.*pow2(-mmgl + mmst2)) - (8*mgl*s2t*pow3(
      mmst1))/(9.*(-mmgl + mmst2)*pow2(-mmgl + mmst1)) + (8*mgl*s2t*pow3(mmst1)
      )/(9.*(-mmgl + mmst1)*pow2(-mmgl + mmst2)) + (8*mgl*s2t*pow3(mmst1))/(9.*
      pow3(-mmgl + mmst1)) - (8*mgl*s2t*pow3(mmst2))/(9.*pow3(-mmgl + mmst2))))
      /mt + logmmglmmu*(64.88888888888889 + (98*mmst1)/(3.*(-mmgl + mmst1)) + (
      640*mmst1)/(9.*(mmst1 - mmst2)) - (2*mmst2)/(-mmgl + mmst1) - (640*mmst2)
      /(9.*(mmst1 - mmst2)) + (14*mmst1)/(3.*(-mmgl + mmst2)) + (118*mmst2)/(3.
      *(-mmgl + mmst2)) - (20*mmsusy)/(-mmgl + mmst1) - (20*mmsusy)/(-mmgl +
      mmst2) - (640*pow2(mmst1))/(9.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (20*
      pow2(mmst1))/(3.*(-mmgl + mmst1)*(-mmgl + mmst2)) + (14*mmst1*mmst2)/(3.*
      pow2(-mmgl + mmst1)) + (140*mmst1*mmsusy)/(3.*pow2(-mmgl + mmst1)) - (116
      *pow2(mmst1))/pow2(-mmgl + mmst1) + (640*pow2(mmst2))/(9.*(mmst1 - mmst2)
      *(-mmgl + mmst2)) + (38*mmst1*mmst2)/(9.*pow2(-mmgl + mmst2)) + (140*
      mmst2*mmsusy)/(3.*pow2(-mmgl + mmst2)) - (4*pow2(mmst1))/(9.*pow2(-mmgl +
      mmst2)) - (1048*pow2(mmst2))/(9.*pow2(-mmgl + mmst2)) + (256*pow2(s2t))/
      9. - (376*mmst1*pow2(s2t))/(9.*(-mmgl + mmst1)) - (640*mmst1*pow2(s2t))/(
      9.*(mmst1 - mmst2)) + (16*mmst2*pow2(s2t))/(9.*(-mmgl + mmst1)) + (640*
      mmst2*pow2(s2t))/(9.*(mmst1 - mmst2)) + (40*mmst1*pow2(s2t))/(3.*(-mmgl +
      mmst2)) - (272*mmst2*pow2(s2t))/(9.*(-mmgl + mmst2)) + (640*pow2(mmst1)*
      pow2(s2t))/(9.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (104*pow2(mmst1)*pow2(
      s2t))/(9.*(-mmgl + mmst1)*(-mmgl + mmst2)) - (32*mmst1*mmst2*pow2(s2t))/(
      9.*pow2(-mmgl + mmst1)) + (170*pow2(mmst1)*pow2(s2t))/(9.*pow2(-mmgl +
      mmst1)) - (640*pow2(mmst2)*pow2(s2t))/(9.*(mmst1 - mmst2)*(-mmgl + mmst2)
      ) - (4*mmst1*mmst2*pow2(s2t))/(9.*pow2(-mmgl + mmst2)) + (28*pow2(mmst1)*
      pow2(s2t))/(9.*pow2(-mmgl + mmst2)) + (22*pow2(mmst2)*pow2(s2t))/pow2(-
      mmgl + mmst2) + (4*pow3(mmst1))/(9.*(-mmgl + mmst2)*pow2(-mmgl + mmst1))
      + (4*pow3(mmst1))/(9.*(-mmgl + mmst1)*pow2(-mmgl + mmst2)) - (28*pow2(s2t
      )*pow3(mmst1))/(9.*(-mmgl + mmst2)*pow2(-mmgl + mmst1)) - (28*pow2(s2t)*
      pow3(mmst1))/(9.*(-mmgl + mmst1)*pow2(-mmgl + mmst2)) - (8*mmst2*pow2(
      mmst1))/(3.*pow3(-mmgl + mmst1)) - (80*mmsusy*pow2(mmst1))/(3.*pow3(-mmgl
       + mmst1)) + (16*mmst2*pow2(mmst1)*pow2(s2t))/(9.*pow3(-mmgl + mmst1)) +
      (508*pow3(mmst1))/(9.*pow3(-mmgl + mmst1)) + (4*pow2(s2t)*pow3(mmst1))/(
      3.*pow3(-mmgl + mmst1)) - (8*mmst1*pow2(mmst2))/(3.*pow3(-mmgl + mmst2))
      - (80*mmsusy*pow2(mmst2))/(3.*pow3(-mmgl + mmst2)) + (16*mmst1*pow2(mmst2
      )*pow2(s2t))/(9.*pow3(-mmgl + mmst2)) + (508*pow3(mmst2))/(9.*pow3(-mmgl
      + mmst2)) + (4*pow2(s2t)*pow3(mmst2))/(3.*pow3(-mmgl + mmst2))) - (112*
      pow3(mmst2))/(3.*pow3(-mmgl + mmst2)) - (16*zt2*pow3(mmst2))/(3.*pow3(-
      mmgl + mmst2)) - (95*mmst1*pow3(mmsusy))/(162.*pow2(-mmgl + mmst1)*pow2(
      mmst1 - mmsusy)) - (10*mmst1*zt2*pow3(mmsusy))/(9.*pow2(-mmgl + mmst1)*
      pow2(mmst1 - mmsusy)) - (95*mmst2*pow3(mmsusy))/(162.*pow2(-mmgl + mmst2)
      *pow2(mmst2 - mmsusy)) - (10*mmst2*zt2*pow3(mmsusy))/(9.*pow2(-mmgl +
      mmst2)*pow2(mmst2 - mmsusy)) + logmmst1mmu*logmmsusymmu*((80*mmst1*mmsusy
      )/(9.*pow2(-mmgl + mmst1)) + (10*mmsusy*pow2(mmst1))/(9.*(mmst1 - mmsusy)
      *pow2(-mmgl + mmst1)) - (20*mmst1*pow2(mmsusy))/(9.*(mmst1 - mmsusy)*pow2
      (-mmgl + mmst1)) + (10*pow2(mmst1)*pow2(mmsusy))/(9.*pow2(-mmgl + mmst1)*
      pow2(mmst1 - mmsusy)) - (40*mmsusy*pow2(mmst1))/(3.*pow3(-mmgl + mmst1))
      - (10*mmst1*pow3(mmsusy))/(9.*pow2(-mmgl + mmst1)*pow2(mmst1 - mmsusy)))
      + logmmst2mmu*logmmsusymmu*((80*mmst2*mmsusy)/(9.*pow2(-mmgl + mmst2)) +
      (10*mmsusy*pow2(mmst2))/(9.*(mmst2 - mmsusy)*pow2(-mmgl + mmst2)) - (20*
      mmst2*pow2(mmsusy))/(9.*(mmst2 - mmsusy)*pow2(-mmgl + mmst2)) + (10*pow2(
      mmst2)*pow2(mmsusy))/(9.*pow2(-mmgl + mmst2)*pow2(mmst2 - mmsusy)) - (40*
      mmsusy*pow2(mmst2))/(3.*pow3(-mmgl + mmst2)) - (10*mmst2*pow3(mmsusy))/(
      9.*pow2(-mmgl + mmst2)*pow2(mmst2 - mmsusy))) + pow2(logmmsusymmu)*(
      3.3333333333333335 + (5*mmsusy)/(-mmgl + mmst1) + (5*mmsusy)/(-mmgl +
      mmst2) - (65*mmst1*mmsusy)/(9.*pow2(-mmgl + mmst1)) + (5*mmsusy*pow2(
      mmst1))/(9.*(mmst1 - mmsusy)*pow2(-mmgl + mmst1)) - (65*mmst2*mmsusy)/(9.
      *pow2(-mmgl + mmst2)) + (5*mmsusy*pow2(mmst2))/(9.*(mmst2 - mmsusy)*pow2(
      -mmgl + mmst2)) - (10*mmst1*pow2(mmsusy))/(9.*(mmst1 - mmsusy)*pow2(-mmgl
       + mmst1)) - (10*mmst2*pow2(mmsusy))/(9.*(mmst2 - mmsusy)*pow2(-mmgl +
      mmst2)) + (5*pow2(mmst1)*pow2(mmsusy))/(9.*pow2(-mmgl + mmst1)*pow2(mmst1
       - mmsusy)) + (5*pow2(mmst2)*pow2(mmsusy))/(9.*pow2(-mmgl + mmst2)*pow2(
      mmst2 - mmsusy)) - (5*mmst1*pow3(mmsusy))/(9.*pow2(-mmgl + mmst1)*pow2(
      mmst1 - mmsusy)) - (5*mmst2*pow3(mmsusy))/(9.*pow2(-mmgl + mmst2)*pow2(
      mmst2 - mmsusy))) + logmmsusymmu*(21.11111111111111 - (40*mmsusy)/(-mmgl
      + mmst1) - (40*mmsusy)/(-mmgl + mmst2) + (170*mmst1*mmsusy)/(3.*pow2(-
      mmgl + mmst1)) - (10*mmsusy*pow2(mmst1))/(3.*(mmst1 - mmsusy)*pow2(-mmgl
      + mmst1)) + (170*mmst2*mmsusy)/(3.*pow2(-mmgl + mmst2)) - (10*mmsusy*pow2
      (mmst2))/(3.*(mmst2 - mmsusy)*pow2(-mmgl + mmst2)) + (20*mmst1*pow2(
      mmsusy))/(3.*(mmst1 - mmsusy)*pow2(-mmgl + mmst1)) + (20*mmst2*pow2(
      mmsusy))/(3.*(mmst2 - mmsusy)*pow2(-mmgl + mmst2)) - (10*pow2(mmst1)*pow2
      (mmsusy))/(3.*pow2(-mmgl + mmst1)*pow2(mmst1 - mmsusy)) - (10*pow2(mmst2)
      *pow2(mmsusy))/(3.*pow2(-mmgl + mmst2)*pow2(mmst2 - mmsusy)) + (10*mmst1*
      pow3(mmsusy))/(3.*pow2(-mmgl + mmst1)*pow2(mmst1 - mmsusy)) + (10*mmst2*
      pow3(mmsusy))/(3.*pow2(-mmgl + mmst2)*pow2(mmst2 - mmsusy))) + (
      logmmst1mmu*((64*mgl*s2t)/9. - (32*mgl*mmst1*s2t)/(9.*(-mmgl + mmst1)) +
      (64*mgl*mmst2*s2t)/(9.*(-mmgl + mmst1)) - (64*mgl*mmst1*mmst2*s2t)/(9.*(-
      mmgl + mmst1)*(mmst1 - mmst2)) - (40*mgl*mmst1*s2t)/(9.*(-mmgl + mmst2))
      - (8*mgl*s2t*pow2(mmst1))/(9.*(-mmgl + mmst1)*(-mmgl + mmst2)) + (8*mgl*
      mmst1*mmst2*s2t)/(3.*pow2(-mmgl + mmst1)) + (80*mgl*mmst1*mmsusy*s2t)/(3.
      *pow2(-mmgl + mmst1)) - (40*mgl*s2t*pow2(mmst1))/pow2(-mmgl + mmst1) + (
      64*mgl*s2t*pow2(mmst2))/(9.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (32*mgl*
      mmst1*pow3(s2t))/(9.*(-mmgl + mmst1)) - (16*mgl*mmst2*pow3(s2t))/(3.*(-
      mmgl + mmst1)) + (64*mgl*mmst1*mmst2*pow3(s2t))/(9.*(-mmgl + mmst1)*(
      mmst1 - mmst2)) + (16*mgl*mmst1*pow3(s2t))/(9.*(-mmgl + mmst2)) - (16*mgl
      *mmst1*mmst2*pow3(s2t))/(9.*pow2(-mmgl + mmst1)) + (16*mgl*pow2(mmst1)*
      pow3(s2t))/(9.*pow2(-mmgl + mmst1)) - (64*mgl*pow2(mmst2)*pow3(s2t))/(9.*
      (-mmgl + mmst1)*(mmst1 - mmst2))))/mt + (pow2(logmmst1mmu)*((22*mgl*mmst1
      *s2t)/(9.*(-mmgl + mmst1)) - (16*mgl*mmst2*s2t)/(9.*(-mmgl + mmst1)) + (
      16*mgl*mmst1*mmst2*s2t)/(3.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (2*mgl*
      mmst1*s2t)/(3.*(-mmgl + mmst2)) - (2*mgl*mmst1*mmst2*s2t)/(3.*pow2(-mmgl
      + mmst1)) - (20*mgl*mmst1*mmsusy*s2t)/(3.*pow2(-mmgl + mmst1)) + (62*mgl*
      s2t*pow2(mmst1))/(9.*pow2(-mmgl + mmst1)) - (16*mgl*s2t*pow2(mmst2))/(9.*
      (-mmgl + mmst1)*(mmst1 - mmst2)) + (8*mgl*s2t*pow3(mmst1))/(9.*pow3(-mmgl
       + mmst1)) - (16*mgl*mmst1*pow3(s2t))/(9.*(-mmgl + mmst1)) + (16*mgl*
      mmst2*pow3(s2t))/(9.*(-mmgl + mmst1)) - (16*mgl*mmst1*mmst2*pow3(s2t))/(
      3.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (16*mgl*pow2(mmst1)*pow3(s2t))/(9.*
      pow2(-mmgl + mmst1)) + (16*mgl*pow2(mmst2)*pow3(s2t))/(9.*(-mmgl + mmst1)
      *(mmst1 - mmst2))))/mt + ((-220*mgl*mmst1*s2t)/(3.*(-mmgl + mmst1)) - (56
      *mgl*mmst2*s2t)/(3.*(-mmgl + mmst1)) + (32*mgl*mmst1*mmst2*s2t)/(3.*(-
      mmgl + mmst1)*(mmst1 - mmst2)) + (8*mgl*mmst1*s2t)/(-mmgl + mmst2) + (220
      *mgl*mmst2*s2t)/(3.*(-mmgl + mmst2)) - (80*mgl*mmsusy*s2t)/(-mmgl + mmst1
      ) + (80*mgl*mmsusy*s2t)/(-mmgl + mmst2) - (88*mgl*mmst1*s2t*zt2)/(9.*(-
      mmgl + mmst1)) - (44*mgl*mmst2*s2t*zt2)/(9.*(-mmgl + mmst1)) + (32*mgl*
      mmst1*mmst2*s2t*zt2)/(9.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (4*mgl*mmst1*
      s2t*zt2)/(3.*(-mmgl + mmst2)) + (88*mgl*mmst2*s2t*zt2)/(9.*(-mmgl + mmst2
      )) - (40*mgl*mmsusy*s2t*zt2)/(3.*(-mmgl + mmst1)) + (40*mgl*mmsusy*s2t*
      zt2)/(3.*(-mmgl + mmst2)) - (32*mgl*s2t*pow2(mmst2))/(3.*(-mmgl + mmst1)*
      (mmst1 - mmst2)) - (32*mgl*s2t*zt2*pow2(mmst2))/(9.*(-mmgl + mmst1)*(
      mmst1 - mmst2)) - (16*mgl*mmst1*pow3(s2t))/(9.*(-mmgl + mmst1)) + (112*
      mgl*mmst2*pow3(s2t))/(9.*(-mmgl + mmst1)) - (32*mgl*mmst1*mmst2*pow3(s2t)
      )/(3.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (16*mgl*mmst1*pow3(s2t))/(9.*(-
      mmgl + mmst2)) + (16*mgl*mmst2*pow3(s2t))/(9.*(-mmgl + mmst2)) + (32*mgl*
      mmst2*zt2*pow3(s2t))/(9.*(-mmgl + mmst1)) - (32*mgl*mmst1*mmst2*zt2*pow3(
      s2t))/(9.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (32*mgl*pow2(mmst2)*pow3(s2t
      ))/(3.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (32*mgl*zt2*pow2(mmst2)*pow3(
      s2t))/(9.*(-mmgl + mmst1)*(mmst1 - mmst2)))/mt + (logmmst1mmu*logmmst2mmu
      *((8*mgl*mmst1*s2t)/(9.*(-mmgl + mmst1)) - (32*mgl*mmst2*s2t)/(9.*(-mmgl
      + mmst1)) - (8*mgl*mmst1*s2t)/(9.*(-mmgl + mmst2)) - (32*mgl*mmst2*s2t)/(
      9.*(-mmgl + mmst2)) + (8*mgl*s2t*pow2(mmst1))/(3.*(-mmgl + mmst1)*(-mmgl
      + mmst2)) - (4*mgl*mmst1*mmst2*s2t)/(3.*pow2(-mmgl + mmst1)) - (8*mgl*s2t
      *pow2(mmst1))/(9.*pow2(-mmgl + mmst1)) - (32*mgl*s2t*pow2(mmst2))/(9.*(-
      mmgl + mmst1)*(mmst1 - mmst2)) - (32*mgl*s2t*pow2(mmst2))/(9.*(mmst1 -
      mmst2)*(-mmgl + mmst2)) + (4*mgl*mmst1*mmst2*s2t)/(9.*pow2(-mmgl + mmst2)
      ) - (8*mgl*s2t*pow2(mmst1))/(9.*pow2(-mmgl + mmst2)) - (8*mgl*s2t*pow3(
      mmst1))/(9.*(-mmgl + mmst2)*pow2(-mmgl + mmst1)) + (8*mgl*s2t*pow3(mmst1)
      )/(9.*(-mmgl + mmst1)*pow2(-mmgl + mmst2)) + (16*mgl*mmst2*pow3(s2t))/(9.
      *(-mmgl + mmst1)) + (16*mgl*mmst1*pow3(s2t))/(9.*(-mmgl + mmst2)) + (32*
      mgl*mmst2*pow3(s2t))/(9.*(-mmgl + mmst2)) + (16*mgl*mmst1*mmst2*pow3(s2t)
      )/(9.*pow2(-mmgl + mmst1)) + (32*mgl*pow2(mmst2)*pow3(s2t))/(9.*(-mmgl +
      mmst1)*(mmst1 - mmst2)) + (32*mgl*pow2(mmst2)*pow3(s2t))/(9.*(mmst1 -
      mmst2)*(-mmgl + mmst2)) - (16*mgl*mmst1*mmst2*pow3(s2t))/(9.*pow2(-mmgl +
      mmst2))))/mt + (logmmglmmu*logmmst1mmu*((-64*mgl*s2t)/9. - (12*mgl*mmst1*
      s2t)/(-mmgl + mmst1) - (32*mgl*mmst1*mmst2*s2t)/(9.*(-mmgl + mmst1)*(
      mmst1 - mmst2)) + (20*mgl*mmst1*s2t)/(9.*(-mmgl + mmst2)) + (32*mgl*mmst2
      *s2t)/(9.*(-mmgl + mmst2)) - (8*mgl*s2t*pow2(mmst1))/(3.*(-mmgl + mmst1)*
      (-mmgl + mmst2)) + (52*mgl*s2t*pow2(mmst1))/(3.*pow2(-mmgl + mmst1)) + (
      32*mgl*s2t*pow2(mmst2))/(9.*(mmst1 - mmst2)*(-mmgl + mmst2)) - (4*mgl*
      mmst1*mmst2*s2t)/(9.*pow2(-mmgl + mmst2)) + (8*mgl*s2t*pow2(mmst1))/(9.*
      pow2(-mmgl + mmst2)) + (8*mgl*s2t*pow3(mmst1))/(9.*(-mmgl + mmst2)*pow2(-
      mmgl + mmst1)) - (8*mgl*s2t*pow3(mmst1))/(9.*(-mmgl + mmst1)*pow2(-mmgl +
      mmst2)) - (16*mgl*s2t*pow3(mmst1))/(9.*pow3(-mmgl + mmst1)) + (16*mgl*
      mmst1*pow3(s2t))/(9.*(-mmgl + mmst1)) + (32*mgl*mmst1*mmst2*pow3(s2t))/(
      9.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (16*mgl*mmst1*pow3(s2t))/(9.*(-mmgl
       + mmst2)) - (32*mgl*mmst2*pow3(s2t))/(9.*(-mmgl + mmst2)) + (16*mgl*pow2
      (mmst1)*pow3(s2t))/(9.*pow2(-mmgl + mmst1)) - (32*mgl*pow2(mmst2)*pow3(
      s2t))/(9.*(mmst1 - mmst2)*(-mmgl + mmst2)) + (16*mgl*mmst1*mmst2*pow3(s2t
      ))/(9.*pow2(-mmgl + mmst2))))/mt + (logmmglmmu*logmmst2mmu*((64*mgl*s2t)/
      9. - (8*mgl*mmst1*s2t)/(9.*(-mmgl + mmst1)) - (4*mgl*mmst2*s2t)/(3.*(-
      mmgl + mmst1)) + (32*mgl*mmst1*mmst2*s2t)/(9.*(-mmgl + mmst1)*(mmst1 -
      mmst2)) + (8*mgl*mmst1*s2t)/(9.*(-mmgl + mmst2)) + (68*mgl*mmst2*s2t)/(9.
      *(-mmgl + mmst2)) - (8*mgl*s2t*pow2(mmst1))/(3.*(-mmgl + mmst1)*(-mmgl +
      mmst2)) + (4*mgl*mmst1*mmst2*s2t)/(3.*pow2(-mmgl + mmst1)) + (8*mgl*s2t*
      pow2(mmst1))/(9.*pow2(-mmgl + mmst1)) - (32*mgl*s2t*pow2(mmst2))/(9.*(
      mmst1 - mmst2)*(-mmgl + mmst2)) + (8*mgl*mmst1*mmst2*s2t)/(9.*pow2(-mmgl
      + mmst2)) + (8*mgl*s2t*pow2(mmst1))/(9.*pow2(-mmgl + mmst2)) - (148*mgl*
      s2t*pow2(mmst2))/(9.*pow2(-mmgl + mmst2)) + (8*mgl*s2t*pow3(mmst1))/(9.*(
      -mmgl + mmst2)*pow2(-mmgl + mmst1)) - (8*mgl*s2t*pow3(mmst1))/(9.*(-mmgl
      + mmst1)*pow2(-mmgl + mmst2)) + (16*mgl*s2t*pow3(mmst2))/(9.*pow3(-mmgl +
      mmst2)) + (16*mgl*mmst2*pow3(s2t))/(9.*(-mmgl + mmst1)) - (32*mgl*mmst1*
      mmst2*pow3(s2t))/(9.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (16*mgl*mmst2*
      pow3(s2t))/(9.*(-mmgl + mmst2)) - (16*mgl*mmst1*mmst2*pow3(s2t))/(9.*pow2
      (-mmgl + mmst1)) + (32*mgl*pow2(mmst2)*pow3(s2t))/(9.*(mmst1 - mmst2)*(-
      mmgl + mmst2)) - (16*mgl*pow2(mmst2)*pow3(s2t))/(9.*pow2(-mmgl + mmst2)))
      )/mt + (logmmst2mmu*((-64*mgl*s2t)/9. + (8*mgl*mmst1*s2t)/(9.*(-mmgl +
      mmst1)) + (112*mgl*mmst2*s2t)/(9.*(-mmgl + mmst1)) - (64*mgl*mmst1*mmst2*
      s2t)/(9.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (8*mgl*mmst1*s2t)/(9.*(-mmgl
      + mmst2)) + (8*mgl*mmst2*s2t)/(3.*(-mmgl + mmst2)) + (8*mgl*s2t*pow2(
      mmst1))/(9.*(-mmgl + mmst1)*(-mmgl + mmst2)) + (64*mgl*s2t*pow2(mmst2))/(
      9.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (8*mgl*mmst1*mmst2*s2t)/(3.*pow2(-
      mmgl + mmst2)) - (80*mgl*mmst2*mmsusy*s2t)/(3.*pow2(-mmgl + mmst2)) + (40
      *mgl*s2t*pow2(mmst2))/pow2(-mmgl + mmst2) - (80*mgl*mmst2*pow3(s2t))/(9.*
      (-mmgl + mmst1)) + (64*mgl*mmst1*mmst2*pow3(s2t))/(9.*(-mmgl + mmst1)*(
      mmst1 - mmst2)) - (16*mgl*mmst1*pow3(s2t))/(9.*(-mmgl + mmst2)) - (32*mgl
      *mmst2*pow3(s2t))/(9.*(-mmgl + mmst2)) - (64*mgl*pow2(mmst2)*pow3(s2t))/(
      9.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (16*mgl*mmst1*mmst2*pow3(s2t))/(9.*
      pow2(-mmgl + mmst2)) - (16*mgl*pow2(mmst2)*pow3(s2t))/(9.*pow2(-mmgl +
      mmst2))))/mt + (pow2(logmmst2mmu)*((-22*mgl*mmst2*s2t)/(9.*(-mmgl + mmst1
      )) + (16*mgl*mmst1*mmst2*s2t)/(9.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (10*
      mgl*mmst2*s2t)/(9.*(-mmgl + mmst2)) - (16*mgl*s2t*pow2(mmst2))/(9.*(-mmgl
       + mmst1)*(mmst1 - mmst2)) + (32*mgl*s2t*pow2(mmst2))/(9.*(mmst1 - mmst2)
      *(-mmgl + mmst2)) + (2*mgl*mmst1*mmst2*s2t)/(3.*pow2(-mmgl + mmst2)) + (
      20*mgl*mmst2*mmsusy*s2t)/(3.*pow2(-mmgl + mmst2)) - (62*mgl*s2t*pow2(
      mmst2))/(9.*pow2(-mmgl + mmst2)) - (8*mgl*s2t*pow3(mmst2))/(9.*pow3(-mmgl
       + mmst2)) + (16*mgl*mmst2*pow3(s2t))/(9.*(-mmgl + mmst1)) - (16*mgl*
      mmst1*mmst2*pow3(s2t))/(9.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (16*mgl*
      mmst2*pow3(s2t))/(9.*(-mmgl + mmst2)) + (16*mgl*pow2(mmst2)*pow3(s2t))/(
      9.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (32*mgl*pow2(mmst2)*pow3(s2t))/(9.*
      (mmst1 - mmst2)*(-mmgl + mmst2)) + (16*mgl*pow2(mmst2)*pow3(s2t))/(9.*
      pow2(-mmgl + mmst2))))/mt + (logmmglmmu*((544*mgl*mmst1*s2t)/(9.*(-mmgl +
      mmst1)) + (8*mgl*mmst2*s2t)/(3.*(-mmgl + mmst1)) - (8*mgl*mmst1*s2t)/(3.*
      (-mmgl + mmst2)) - (544*mgl*mmst2*s2t)/(9.*(-mmgl + mmst2)) + (80*mgl*
      mmsusy*s2t)/(3.*(-mmgl + mmst1)) - (80*mgl*mmsusy*s2t)/(3.*(-mmgl + mmst2
      )) - (8*mgl*mmst1*mmst2*s2t)/(3.*pow2(-mmgl + mmst1)) - (80*mgl*mmst1*
      mmsusy*s2t)/(3.*pow2(-mmgl + mmst1)) + (40*mgl*s2t*pow2(mmst1))/pow2(-
      mmgl + mmst1) + (8*mgl*mmst1*mmst2*s2t)/(3.*pow2(-mmgl + mmst2)) + (80*
      mgl*mmst2*mmsusy*s2t)/(3.*pow2(-mmgl + mmst2)) - (40*mgl*s2t*pow2(mmst2))
      /pow2(-mmgl + mmst2) - (16*mgl*mmst1*pow3(s2t))/(9.*(-mmgl + mmst1)) - (
      16*mgl*mmst2*pow3(s2t))/(9.*(-mmgl + mmst1)) + (16*mgl*mmst1*pow3(s2t))/(
      9.*(-mmgl + mmst2)) + (16*mgl*mmst2*pow3(s2t))/(9.*(-mmgl + mmst2)) + (16
      *mgl*mmst1*mmst2*pow3(s2t))/(9.*pow2(-mmgl + mmst1)) - (16*mgl*pow2(mmst1
      )*pow3(s2t))/(9.*pow2(-mmgl + mmst1)) - (16*mgl*mmst1*mmst2*pow3(s2t))/(
      9.*pow2(-mmgl + mmst2)) + (16*mgl*pow2(mmst2)*pow3(s2t))/(9.*pow2(-mmgl +
      mmst2))))/mt + logmmglmmu*logmmst1mmu*(7.111111111111111 - (59*mmst1)/(9.
      *(-mmgl + mmst1)) + (13*mmst1)/(9.*(-mmgl + mmst2)) + (4*pow2(mmst1))/(9.
      *(-mmgl + mmst1)*(-mmgl + mmst2)) - (109*pow2(mmst1))/(9.*pow2(-mmgl +
      mmst1)) - (29*mmst1*mmst2)/(9.*pow2(-mmgl + mmst2)) - (4*pow2(mmst1))/(3.
      *pow2(-mmgl + mmst2)) + (64*pow2(s2t))/9. - (67*mmst1*pow2(s2t))/(9.*(-
      mmgl + mmst1)) - (61*mmst1*pow2(s2t))/(9.*(-mmgl + mmst2)) + (13*pow2(
      mmst1)*pow2(s2t))/((-mmgl + mmst1)*(-mmgl + mmst2)) - (97*pow2(mmst1)*
      pow2(s2t))/(9.*pow2(-mmgl + mmst1)) + (56*mmst1*mmst2*pow2(s2t))/(9.*pow2
      (-mmgl + mmst2)) - (4*pow2(mmst1)*pow2(s2t))/(3.*pow2(-mmgl + mmst2)) + (
      16*pow3(mmst1))/(9.*(-mmgl + mmst1)*pow2(-mmgl + mmst2)) - (32*pow2(s2t)*
      pow3(mmst1))/(3.*(-mmgl + mmst2)*pow2(-mmgl + mmst1)) + (16*pow2(s2t)*
      pow3(mmst1))/(3.*(-mmgl + mmst1)*pow2(-mmgl + mmst2)) + (116*pow3(mmst1))
      /(9.*pow3(-mmgl + mmst1)) + (160*pow2(s2t)*pow3(mmst1))/(9.*pow3(-mmgl +
      mmst1)) + (4*mmst1*pow2(mmst2))/(3.*pow3(-mmgl + mmst2)) - (32*mmst2*pow2
      (mmst1)*pow2(s2t))/(9.*pow3(-mmgl + mmst2)) - (16*mmst1*pow2(mmst2)*pow2(
      s2t))/(3.*pow3(-mmgl + mmst2)) - (32*pow2(s2t)*pow3(mmst1))/(9.*pow3(-
      mmgl + mmst2)) - (4*pow4(mmst1))/(9.*pow2(-mmgl + mmst1)*pow2(-mmgl +
      mmst2)) - (4*pow2(s2t)*pow4(mmst1))/(9.*pow2(-mmgl + mmst1)*pow2(-mmgl +
      mmst2)) + (32*pow2(s2t)*pow4(mmst1))/(9.*(-mmgl + mmst2)*pow3(-mmgl +
      mmst1)) + (32*pow2(s2t)*pow4(mmst1))/(9.*(-mmgl + mmst1)*pow3(-mmgl +
      mmst2)) - (4*pow4(mmst1))/(3.*pow4(-mmgl + mmst1)) - (20*pow2(s2t)*pow4(
      mmst1))/(3.*pow4(-mmgl + mmst1))) + (271*pow4(mmst2))/(648.*pow2(-mmgl +
      mmst1)*pow2(mmst1 - mmst2)) + (zt2*pow4(mmst2))/(9.*pow2(-mmgl + mmst1)*
      pow2(mmst1 - mmst2)) + logmmst2mmu*(-5 + (4*mmst1)/(9.*(-mmgl + mmst1)) -
      (4*mmst2)/(-mmgl + mmst1) - (128*mmst2)/(9.*(mmst1 - mmst2)) - (4*mmst1)/
      (9.*(-mmgl + mmst2)) + (178*mmst2)/(9.*(-mmgl + mmst2)) + (16*mmst1*mmst2
      )/(3.*pow2(-mmgl + mmst1)) + (128*pow2(mmst2))/(9.*(mmst1 - mmst2)*(-mmgl
       + mmst2)) - (14*mmst1*mmst2)/(9.*pow2(-mmgl + mmst2)) - (515*mmst2*
      mmsusy)/(27.*pow2(-mmgl + mmst2)) + (4*pow2(mmst1))/(9.*pow2(-mmgl +
      mmst2)) + (248*pow2(mmst2))/(9.*pow2(-mmgl + mmst2)) - (25*mmsusy*pow2(
      mmst2))/(27.*(mmst2 - mmsusy)*pow2(-mmgl + mmst2)) - (5*mmst2*pow2(mmsusy
      ))/(9.*(mmst2 - mmsusy)*pow2(-mmgl + mmst2)) + (40*pow2(mmst2)*pow2(
      mmsusy))/(27.*pow2(-mmgl + mmst2)*pow2(mmst2 - mmsusy)) - (64*pow2(s2t))/
      9. - (4385*mmst1*pow2(s2t))/(324.*(-mmgl + mmst1)) + (871*mmst2*pow2(s2t)
      )/(324.*(-mmgl + mmst1)) + (128*mmst2*pow2(s2t))/(9.*(mmst1 - mmst2)) - (
      11*mmst1*mmst2*pow2(s2t))/(162.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (3809*
      mmst1*pow2(s2t))/(324.*(-mmgl + mmst2)) + (8993*mmst2*pow2(s2t))/(324.*(-
      mmgl + mmst2)) + (11*mmst1*mmst2*pow2(s2t))/(162.*(mmst1 - mmst2)*(-mmgl
      + mmst2)) - (11*pow2(mmst1)*pow2(s2t))/(162.*(-mmgl + mmst1)*(mmst1 -
      mmst2)) - (5681*pow2(mmst1)*pow2(s2t))/(324.*(-mmgl + mmst1)*(-mmgl +
      mmst2)) + (11*pow2(mmst1)*pow2(s2t))/(162.*(mmst1 - mmst2)*(-mmgl + mmst2
      )) - (16*mmst1*mmst2*pow2(s2t))/(9.*pow2(-mmgl + mmst1)) + (32*pow2(mmst1
      )*pow2(s2t))/(9.*pow2(-mmgl + mmst1)) - (5*pow2(mmst2)*pow2(s2t))/(108.*(
      -mmgl + mmst1)*(mmst1 - mmst2)) - (1531*pow2(mmst2)*pow2(s2t))/(108.*(
      mmst1 - mmst2)*(-mmgl + mmst2)) + (4*mmst1*mmst2*pow2(s2t))/pow2(-mmgl +
      mmst2) + (4*pow2(mmst1)*pow2(s2t))/(9.*pow2(-mmgl + mmst2)) - (182*pow2(
      mmst2)*pow2(s2t))/(9.*pow2(-mmgl + mmst2)) - (4*pow3(mmst1))/(9.*(-mmgl +
      mmst1)*pow2(-mmgl + mmst2)) - (11*pow2(s2t)*pow3(mmst1))/(162.*(-mmgl +
      mmst1)*(mmst1 - mmst2)*(-mmgl + mmst2)) + (32*pow2(s2t)*pow3(mmst1))/(9.*
      (-mmgl + mmst2)*pow2(-mmgl + mmst1)) - (4*pow2(s2t)*pow3(mmst1))/(9.*(-
      mmgl + mmst1)*pow2(-mmgl + mmst2)) - pow3(mmst2)/(3.*(mmst1 - mmst2)*pow2
      (-mmgl + mmst1)) + (mmst1*pow3(mmst2))/(3.*pow2(-mmgl + mmst1)*pow2(mmst1
       - mmst2)) - (11*pow2(s2t)*pow3(mmst2))/(162.*(-mmgl + mmst1)*pow2(mmst1
      - mmst2)) + (11*pow2(s2t)*pow3(mmst2))/(162.*(-mmgl + mmst2)*pow2(mmst1 -
      mmst2)) + (8*mmst1*pow2(mmst2))/(3.*pow3(-mmgl + mmst2)) + (80*mmsusy*
      pow2(mmst2))/(3.*pow3(-mmgl + mmst2)) - (16*mmst1*pow2(mmst2)*pow2(s2t))/
      (9.*pow3(-mmgl + mmst2)) - (220*pow3(mmst2))/(9.*pow3(-mmgl + mmst2)) - (
      4*pow2(s2t)*pow3(mmst2))/(3.*pow3(-mmgl + mmst2)) - (40*mmst2*pow3(mmsusy
      ))/(27.*pow2(-mmgl + mmst2)*pow2(mmst2 - mmsusy)) - pow4(mmst2)/(3.*pow2(
      -mmgl + mmst1)*pow2(mmst1 - mmst2))) + logmmst1mmu*(-5 + (62*mmst1)/(3.*(
      -mmgl + mmst1)) + (128*mmst1)/(9.*(mmst1 - mmst2)) - (44*mmst1)/(9.*(-
      mmgl + mmst2)) - (128*pow2(mmst1))/(9.*(-mmgl + mmst1)*(mmst1 - mmst2)) +
      (4*pow2(mmst1))/(3.*(-mmgl + mmst1)*(-mmgl + mmst2)) - (41*mmst1*mmst2)/(
      27.*pow2(-mmgl + mmst1)) - (515*mmst1*mmsusy)/(27.*pow2(-mmgl + mmst1)) +
      (244*pow2(mmst1))/(9.*pow2(-mmgl + mmst1)) - (13*mmst2*pow2(mmst1))/(27.*
      (mmst1 - mmst2)*pow2(-mmgl + mmst1)) - (25*mmsusy*pow2(mmst1))/(27.*(
      mmst1 - mmsusy)*pow2(-mmgl + mmst1)) + (13*mmst1*pow2(mmst2))/(18.*(mmst1
       - mmst2)*pow2(-mmgl + mmst1)) - (13*pow2(mmst1)*pow2(mmst2))/(54.*pow2(-
      mmgl + mmst1)*pow2(mmst1 - mmst2)) + (16*mmst1*mmst2)/(3.*pow2(-mmgl +
      mmst2)) - (5*mmst1*pow2(mmsusy))/(9.*(mmst1 - mmsusy)*pow2(-mmgl + mmst1)
      ) + (40*pow2(mmst1)*pow2(mmsusy))/(27.*pow2(-mmgl + mmst1)*pow2(mmst1 -
      mmsusy)) - (64*pow2(s2t))/9. + (2959*mmst1*pow2(s2t))/(162.*(-mmgl +
      mmst1)) - (128*mmst1*pow2(s2t))/(9.*(mmst1 - mmst2)) - (16*mmst2*pow2(s2t
      ))/(9.*(-mmgl + mmst1)) - (2*mmst1*mmst2*pow2(s2t))/(81.*(-mmgl + mmst1)*
      (mmst1 - mmst2)) + (1973*mmst1*pow2(s2t))/(162.*(-mmgl + mmst2)) - (7*
      mmst2*pow2(s2t))/(162.*(-mmgl + mmst2)) + (11*mmst1*mmst2*pow2(s2t))/(
      162.*(mmst1 - mmst2)*(-mmgl + mmst2)) + (2293*pow2(mmst1)*pow2(s2t))/(
      162.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (893*pow2(mmst1)*pow2(s2t))/(162.
      *(-mmgl + mmst1)*(-mmgl + mmst2)) + (11*pow2(mmst1)*pow2(s2t))/(162.*(
      mmst1 - mmst2)*(-mmgl + mmst2)) + (32*mmst1*mmst2*pow2(s2t))/(9.*pow2(-
      mmgl + mmst1)) - (62*pow2(mmst1)*pow2(s2t))/(3.*pow2(-mmgl + mmst1)) - (
      29*pow2(mmst2)*pow2(s2t))/(324.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (2*
      pow2(mmst2)*pow2(s2t))/(81.*(mmst1 - mmst2)*(-mmgl + mmst2)) + (7*mmst1*
      pow2(mmst2)*pow2(s2t))/(324.*(-mmgl + mmst1)*pow2(mmst1 - mmst2)) - (16*
      mmst1*mmst2*pow2(s2t))/(3.*pow2(-mmgl + mmst2)) - (32*pow2(mmst1)*pow2(
      s2t))/(9.*pow2(-mmgl + mmst2)) - (4*pow3(mmst1))/(9.*(-mmgl + mmst2)*pow2
      (-mmgl + mmst1)) - (11*pow2(s2t)*pow3(mmst1))/(162.*(-mmgl + mmst1)*(
      mmst1 - mmst2)*(-mmgl + mmst2)) - (4*pow2(s2t)*pow3(mmst1))/(9.*(-mmgl +
      mmst2)*pow2(-mmgl + mmst1)) + (32*pow2(s2t)*pow3(mmst1))/(9.*(-mmgl +
      mmst1)*pow2(-mmgl + mmst2)) + (8*mmst2*pow2(mmst1))/(3.*pow3(-mmgl +
      mmst1)) + (80*mmsusy*pow2(mmst1))/(3.*pow3(-mmgl + mmst1)) - (16*mmst2*
      pow2(mmst1)*pow2(s2t))/(9.*pow3(-mmgl + mmst1)) - (220*pow3(mmst1))/(9.*
      pow3(-mmgl + mmst1)) - (4*pow2(s2t)*pow3(mmst1))/(3.*pow3(-mmgl + mmst1))
      - (5*pow3(mmst2))/(54.*(mmst1 - mmst2)*pow2(-mmgl + mmst1)) + (mmst1*pow3
      (mmst2))/(3.*pow2(-mmgl + mmst1)*pow2(mmst1 - mmst2)) - (29*pow2(s2t)*
      pow3(mmst2))/(324.*(-mmgl + mmst1)*pow2(mmst1 - mmst2)) + (11*pow2(s2t)*
      pow3(mmst2))/(162.*(-mmgl + mmst2)*pow2(mmst1 - mmst2)) - (40*mmst1*pow3(
      mmsusy))/(27.*pow2(-mmgl + mmst1)*pow2(mmst1 - mmsusy)) - (5*pow4(mmst2))
      /(54.*pow2(-mmgl + mmst1)*pow2(mmst1 - mmst2))) + pow2(logmmst1mmu)*(
      0.3333333333333333 - (35*mmst1)/(18.*(-mmgl + mmst1)) - (64*mmst1)/(9.*(
      mmst1 - mmst2)) + mmst1/(2.*(-mmgl + mmst2)) + (64*pow2(mmst1))/(9.*(-
      mmgl + mmst1)*(mmst1 - mmst2)) + (mmst1*mmst2)/(2.*pow2(-mmgl + mmst1)) +
      (40*mmst1*mmsusy)/(9.*pow2(-mmgl + mmst1)) - (25*pow2(mmst1))/(6.*pow2(-
      mmgl + mmst1)) + (5*mmsusy*pow2(mmst1))/(9.*(mmst1 - mmsusy)*pow2(-mmgl +
      mmst1)) - (2*mmst1*mmst2)/(3.*pow2(-mmgl + mmst2)) - (10*mmst1*pow2(
      mmsusy))/(9.*(mmst1 - mmsusy)*pow2(-mmgl + mmst1)) + (5*pow2(mmst1)*pow2(
      mmsusy))/(9.*pow2(-mmgl + mmst1)*pow2(mmst1 - mmsusy)) - (25*mmst1*pow2(
      s2t))/(9.*(-mmgl + mmst1)) + (64*mmst1*pow2(s2t))/(9.*(mmst1 - mmst2)) +
      (mmst1*mmst2*pow2(s2t))/(54.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (13*mmst1
      *pow2(s2t))/(9.*(-mmgl + mmst2)) - (mmst1*mmst2*pow2(s2t))/(54.*(mmst1 -
      mmst2)*(-mmgl + mmst2)) - (383*pow2(mmst1)*pow2(s2t))/(54.*(-mmgl + mmst1
      )*(mmst1 - mmst2)) + (13*pow2(mmst1)*pow2(s2t))/(9.*(-mmgl + mmst1)*(-
      mmgl + mmst2)) - (pow2(mmst1)*pow2(s2t))/(54.*(mmst1 - mmst2)*(-mmgl +
      mmst2)) + (185*pow2(mmst1)*pow2(s2t))/(18.*pow2(-mmgl + mmst1)) + (pow2(
      mmst2)*pow2(s2t))/(54.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (pow2(mmst2)*
      pow2(s2t))/(54.*(mmst1 - mmst2)*(-mmgl + mmst2)) + (pow2(s2t)*pow3(mmst1)
      )/(54.*(-mmgl + mmst1)*(mmst1 - mmst2)*(-mmgl + mmst2)) - (2*mmst2*pow2(
      mmst1))/(3.*pow3(-mmgl + mmst1)) - (20*mmsusy*pow2(mmst1))/(3.*pow3(-mmgl
       + mmst1)) + (10*pow3(mmst1))/(3.*pow3(-mmgl + mmst1)) - (8*pow2(s2t)*
      pow3(mmst1))/pow3(-mmgl + mmst1) + pow3(mmst2)/(18.*(mmst1 - mmst2)*pow2(
      -mmgl + mmst1)) - (mmst1*pow3(mmst2))/(18.*pow2(-mmgl + mmst1)*pow2(mmst1
       - mmst2)) + (pow2(s2t)*pow3(mmst2))/(54.*(-mmgl + mmst1)*pow2(mmst1 -
      mmst2)) - (pow2(s2t)*pow3(mmst2))/(54.*(-mmgl + mmst2)*pow2(mmst1 - mmst2
      )) - (5*mmst1*pow3(mmsusy))/(9.*pow2(-mmgl + mmst1)*pow2(mmst1 - mmsusy))
      + (2*pow4(mmst1))/(3.*pow4(-mmgl + mmst1)) + (10*pow2(s2t)*pow4(mmst1))/(
      3.*pow4(-mmgl + mmst1)) + pow4(mmst2)/(18.*pow2(-mmgl + mmst1)*pow2(mmst1
       - mmst2))) + logmmst1mmu*logmmst2mmu*((4*mmst1)/(9.*(-mmgl + mmst1)) - (
      4*mmst1)/(9.*(-mmgl + mmst2)) - (4*pow2(mmst1))/(9.*(-mmgl + mmst1)*(-
      mmgl + mmst2)) + (mmst1*mmst2)/pow2(-mmgl + mmst1) - (4*pow2(mmst1))/(9.*
      pow2(-mmgl + mmst1)) + (17*mmst1*mmst2)/(9.*pow2(-mmgl + mmst2)) + (4*
      pow2(mmst1))/(3.*pow2(-mmgl + mmst2)) - (19*mmst1*pow2(s2t))/(9.*(-mmgl +
      mmst1)) + (16*mmst2*pow2(s2t))/(9.*(-mmgl + mmst1)) + (mmst1*mmst2*pow2(
      s2t))/(27.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (35*mmst1*pow2(s2t))/(9.*(-
      mmgl + mmst2)) - (mmst1*mmst2*pow2(s2t))/(27.*(mmst1 - mmst2)*(-mmgl +
      mmst2)) + (pow2(mmst1)*pow2(s2t))/(27.*(-mmgl + mmst1)*(mmst1 - mmst2)) -
      (91*pow2(mmst1)*pow2(s2t))/(9.*(-mmgl + mmst1)*(-mmgl + mmst2)) - (pow2(
      mmst1)*pow2(s2t))/(27.*(mmst1 - mmst2)*(-mmgl + mmst2)) - (32*mmst1*mmst2
      *pow2(s2t))/(9.*pow2(-mmgl + mmst1)) + (20*pow2(mmst1)*pow2(s2t))/(3.*
      pow2(-mmgl + mmst1)) + (pow2(mmst2)*pow2(s2t))/(27.*(-mmgl + mmst1)*(
      mmst1 - mmst2)) - (pow2(mmst2)*pow2(s2t))/(27.*(mmst1 - mmst2)*(-mmgl +
      mmst2)) - (56*mmst1*mmst2*pow2(s2t))/(9.*pow2(-mmgl + mmst2)) + (4*pow2(
      mmst1)*pow2(s2t))/(3.*pow2(-mmgl + mmst2)) - (16*pow3(mmst1))/(9.*(-mmgl
      + mmst1)*pow2(-mmgl + mmst2)) + (pow2(s2t)*pow3(mmst1))/(27.*(-mmgl +
      mmst1)*(mmst1 - mmst2)*(-mmgl + mmst2)) + (32*pow2(s2t)*pow3(mmst1))/(3.*
      (-mmgl + mmst2)*pow2(-mmgl + mmst1)) - (16*pow2(s2t)*pow3(mmst1))/(3.*(-
      mmgl + mmst1)*pow2(-mmgl + mmst2)) - (4*mmst2*pow2(mmst1))/(3.*pow3(-mmgl
       + mmst1)) + (16*mmst2*pow2(mmst1)*pow2(s2t))/(9.*pow3(-mmgl + mmst1)) -
      (32*pow2(s2t)*pow3(mmst1))/(9.*pow3(-mmgl + mmst1)) + pow3(mmst2)/(9.*(
      mmst1 - mmst2)*pow2(-mmgl + mmst1)) - (mmst1*pow3(mmst2))/(9.*pow2(-mmgl
      + mmst1)*pow2(mmst1 - mmst2)) + (pow2(s2t)*pow3(mmst2))/(27.*(-mmgl +
      mmst1)*pow2(mmst1 - mmst2)) - (pow2(s2t)*pow3(mmst2))/(27.*(-mmgl + mmst2
      )*pow2(mmst1 - mmst2)) - (4*mmst1*pow2(mmst2))/(3.*pow3(-mmgl + mmst2)) +
      (32*mmst2*pow2(mmst1)*pow2(s2t))/(9.*pow3(-mmgl + mmst2)) + (16*mmst1*
      pow2(mmst2)*pow2(s2t))/(3.*pow3(-mmgl + mmst2)) + (32*pow2(s2t)*pow3(
      mmst1))/(9.*pow3(-mmgl + mmst2)) + (4*pow4(mmst1))/(9.*pow2(-mmgl + mmst1
      )*pow2(-mmgl + mmst2)) + (4*pow2(s2t)*pow4(mmst1))/(9.*pow2(-mmgl + mmst1
      )*pow2(-mmgl + mmst2)) - (32*pow2(s2t)*pow4(mmst1))/(9.*(-mmgl + mmst2)*
      pow3(-mmgl + mmst1)) - (32*pow2(s2t)*pow4(mmst1))/(9.*(-mmgl + mmst1)*
      pow3(-mmgl + mmst2)) + pow4(mmst2)/(9.*pow2(-mmgl + mmst1)*pow2(mmst1 -
      mmst2))) + logmmglmmu*logmmst2mmu*(7.111111111111111 - (4*mmst1)/(9.*(-
      mmgl + mmst1)) + mmst2/(-mmgl + mmst1) + (4*mmst1)/(9.*(-mmgl + mmst2)) -
      (55*mmst2)/(9.*(-mmgl + mmst2)) + (4*pow2(mmst1))/(9.*(-mmgl + mmst1)*(-
      mmgl + mmst2)) - (7*mmst1*mmst2)/(3.*pow2(-mmgl + mmst1)) + (4*pow2(mmst1
      ))/(9.*pow2(-mmgl + mmst1)) - (8*mmst1*mmst2)/(9.*pow2(-mmgl + mmst2)) -
      (4*pow2(mmst1))/(3.*pow2(-mmgl + mmst2)) - (113*pow2(mmst2))/(9.*pow2(-
      mmgl + mmst2)) + (64*pow2(s2t))/9. + (5*mmst1*pow2(s2t))/(-mmgl + mmst1)
      - (16*mmst2*pow2(s2t))/(9.*(-mmgl + mmst1)) - (5*mmst1*pow2(s2t))/(-mmgl
      + mmst2) - (112*mmst2*pow2(s2t))/(9.*(-mmgl + mmst2)) + (13*pow2(mmst1)*
      pow2(s2t))/((-mmgl + mmst1)*(-mmgl + mmst2)) + (32*mmst1*mmst2*pow2(s2t))
      /(9.*pow2(-mmgl + mmst1)) - (20*pow2(mmst1)*pow2(s2t))/(3.*pow2(-mmgl +
      mmst1)) + (8*mmst1*mmst2*pow2(s2t))/(3.*pow2(-mmgl + mmst2)) - (4*pow2(
      mmst1)*pow2(s2t))/(3.*pow2(-mmgl + mmst2)) - (37*pow2(mmst2)*pow2(s2t))/(
      9.*pow2(-mmgl + mmst2)) + (16*pow3(mmst1))/(9.*(-mmgl + mmst1)*pow2(-mmgl
       + mmst2)) - (32*pow2(s2t)*pow3(mmst1))/(3.*(-mmgl + mmst2)*pow2(-mmgl +
      mmst1)) + (16*pow2(s2t)*pow3(mmst1))/(3.*(-mmgl + mmst1)*pow2(-mmgl +
      mmst2)) + (4*mmst2*pow2(mmst1))/(3.*pow3(-mmgl + mmst1)) - (16*mmst2*pow2
      (mmst1)*pow2(s2t))/(9.*pow3(-mmgl + mmst1)) + (32*pow2(s2t)*pow3(mmst1))/
      (9.*pow3(-mmgl + mmst1)) - (32*mmst2*pow2(mmst1)*pow2(s2t))/(9.*pow3(-
      mmgl + mmst2)) - (32*mmst1*pow2(mmst2)*pow2(s2t))/(9.*pow3(-mmgl + mmst2)
      ) - (32*pow2(s2t)*pow3(mmst1))/(9.*pow3(-mmgl + mmst2)) + (116*pow3(mmst2
      ))/(9.*pow3(-mmgl + mmst2)) + (128*pow2(s2t)*pow3(mmst2))/(9.*pow3(-mmgl
      + mmst2)) - (4*pow4(mmst1))/(9.*pow2(-mmgl + mmst1)*pow2(-mmgl + mmst2))
      - (4*pow2(s2t)*pow4(mmst1))/(9.*pow2(-mmgl + mmst1)*pow2(-mmgl + mmst2))
      + (32*pow2(s2t)*pow4(mmst1))/(9.*(-mmgl + mmst2)*pow3(-mmgl + mmst1)) + (
      32*pow2(s2t)*pow4(mmst1))/(9.*(-mmgl + mmst1)*pow3(-mmgl + mmst2)) - (4*
      pow4(mmst2))/(3.*pow4(-mmgl + mmst2)) - (20*pow2(s2t)*pow4(mmst2))/(3.*
      pow4(-mmgl + mmst2))) + pow2(logmmglmmu)*(-9.555555555555555 - (211*mmst1
      )/(18.*(-mmgl + mmst1)) - (64*mmst1)/(3.*(mmst1 - mmst2)) + mmst2/(2.*(-
      mmgl + mmst1)) + (64*mmst2)/(3.*(mmst1 - mmst2)) - (31*mmst1)/(18.*(-mmgl
       + mmst2)) - (251*mmst2)/(18.*(-mmgl + mmst2)) + (5*mmsusy)/(-mmgl +
      mmst1) + (5*mmsusy)/(-mmgl + mmst2) + (64*pow2(mmst1))/(3.*(-mmgl + mmst1
      )*(mmst1 - mmst2)) + (4*pow2(mmst1))/(3.*(-mmgl + mmst1)*(-mmgl + mmst2))
      - (7*mmst1*mmst2)/(6.*pow2(-mmgl + mmst1)) - (35*mmst1*mmsusy)/(3.*pow2(-
      mmgl + mmst1)) + (821*pow2(mmst1))/(18.*pow2(-mmgl + mmst1)) - (64*pow2(
      mmst2))/(3.*(mmst1 - mmst2)*(-mmgl + mmst2)) - (5*mmst1*mmst2)/(18.*pow2(
      -mmgl + mmst2)) - (35*mmst2*mmsusy)/(3.*pow2(-mmgl + mmst2)) + (4*pow2(
      mmst1))/(3.*pow2(-mmgl + mmst2)) + (829*pow2(mmst2))/(18.*pow2(-mmgl +
      mmst2)) - (128*pow2(s2t))/9. + (15*mmst1*pow2(s2t))/(-mmgl + mmst1) + (64
      *mmst1*pow2(s2t))/(3.*(mmst1 - mmst2)) - (64*mmst2*pow2(s2t))/(3.*(mmst1
      - mmst2)) - (7*mmst1*pow2(s2t))/(9.*(-mmgl + mmst2)) + (128*mmst2*pow2(
      s2t))/(9.*(-mmgl + mmst2)) - (64*pow2(mmst1)*pow2(s2t))/(3.*(-mmgl +
      mmst1)*(mmst1 - mmst2)) - (65*pow2(mmst1)*pow2(s2t))/(9.*(-mmgl + mmst1)*
      (-mmgl + mmst2)) + (pow2(mmst1)*pow2(s2t))/(2.*pow2(-mmgl + mmst1)) + (64
      *pow2(mmst2)*pow2(s2t))/(3.*(mmst1 - mmst2)*(-mmgl + mmst2)) - (8*mmst1*
      mmst2*pow2(s2t))/(3.*pow2(-mmgl + mmst2)) + (4*pow2(mmst1)*pow2(s2t))/(3.
      *pow2(-mmgl + mmst2)) - (37*pow2(mmst2)*pow2(s2t))/(6.*pow2(-mmgl + mmst2
      )) - (16*pow3(mmst1))/(9.*(-mmgl + mmst1)*pow2(-mmgl + mmst2)) + (32*pow2
      (s2t)*pow3(mmst1))/(3.*(-mmgl + mmst2)*pow2(-mmgl + mmst1)) - (16*pow2(
      s2t)*pow3(mmst1))/(3.*(-mmgl + mmst1)*pow2(-mmgl + mmst2)) + (2*mmst2*
      pow2(mmst1))/(3.*pow3(-mmgl + mmst1)) + (20*mmsusy*pow2(mmst1))/(3.*pow3(
      -mmgl + mmst1)) - (242*pow3(mmst1))/(9.*pow3(-mmgl + mmst1)) - (88*pow2(
      s2t)*pow3(mmst1))/(9.*pow3(-mmgl + mmst1)) + (2*mmst1*pow2(mmst2))/(3.*
      pow3(-mmgl + mmst2)) + (20*mmsusy*pow2(mmst2))/(3.*pow3(-mmgl + mmst2)) +
      (32*mmst2*pow2(mmst1)*pow2(s2t))/(9.*pow3(-mmgl + mmst2)) + (32*mmst1*
      pow2(mmst2)*pow2(s2t))/(9.*pow3(-mmgl + mmst2)) + (32*pow2(s2t)*pow3(
      mmst1))/(9.*pow3(-mmgl + mmst2)) - (242*pow3(mmst2))/(9.*pow3(-mmgl +
      mmst2)) - (56*pow2(s2t)*pow3(mmst2))/(9.*pow3(-mmgl + mmst2)) + (4*pow4(
      mmst1))/(9.*pow2(-mmgl + mmst1)*pow2(-mmgl + mmst2)) + (4*pow2(s2t)*pow4(
      mmst1))/(9.*pow2(-mmgl + mmst1)*pow2(-mmgl + mmst2)) - (32*pow2(s2t)*pow4
      (mmst1))/(9.*(-mmgl + mmst2)*pow3(-mmgl + mmst1)) - (32*pow2(s2t)*pow4(
      mmst1))/(9.*(-mmgl + mmst1)*pow3(-mmgl + mmst2)) + (2*pow4(mmst1))/(3.*
      pow4(-mmgl + mmst1)) + (10*pow2(s2t)*pow4(mmst1))/(3.*pow4(-mmgl + mmst1)
      ) + (2*pow4(mmst2))/(3.*pow4(-mmgl + mmst2)) + (10*pow2(s2t)*pow4(mmst2))
      /(3.*pow4(-mmgl + mmst2))) + pow2(logmmst2mmu)*(0.3333333333333333 +
      mmst2/(2.*(-mmgl + mmst1)) + (64*mmst2)/(9.*(mmst1 - mmst2)) - (35*mmst2)
      /(18.*(-mmgl + mmst2)) - (2*mmst1*mmst2)/(3.*pow2(-mmgl + mmst1)) - (64*
      pow2(mmst2))/(9.*(mmst1 - mmst2)*(-mmgl + mmst2)) + (mmst1*mmst2)/(2.*
      pow2(-mmgl + mmst2)) + (40*mmst2*mmsusy)/(9.*pow2(-mmgl + mmst2)) - (25*
      pow2(mmst2))/(6.*pow2(-mmgl + mmst2)) + (5*mmsusy*pow2(mmst2))/(9.*(mmst2
       - mmsusy)*pow2(-mmgl + mmst2)) - (10*mmst2*pow2(mmsusy))/(9.*(mmst2 -
      mmsusy)*pow2(-mmgl + mmst2)) + (5*pow2(mmst2)*pow2(mmsusy))/(9.*pow2(-
      mmgl + mmst2)*pow2(mmst2 - mmsusy)) + (13*mmst1*pow2(s2t))/(9.*(-mmgl +
      mmst1)) - (64*mmst2*pow2(s2t))/(9.*(mmst1 - mmst2)) + (mmst1*mmst2*pow2(
      s2t))/(54.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (13*mmst1*pow2(s2t))/(9.*(-
      mmgl + mmst2)) - (38*mmst2*pow2(s2t))/(9.*(-mmgl + mmst2)) - (mmst1*mmst2
      *pow2(s2t))/(54.*(mmst1 - mmst2)*(-mmgl + mmst2)) + (pow2(mmst1)*pow2(s2t
      ))/(54.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (13*pow2(mmst1)*pow2(s2t))/(9.
      *(-mmgl + mmst1)*(-mmgl + mmst2)) - (pow2(mmst1)*pow2(s2t))/(54.*(mmst1 -
      mmst2)*(-mmgl + mmst2)) + (pow2(mmst2)*pow2(s2t))/(54.*(-mmgl + mmst1)*(
      mmst1 - mmst2)) + (383*pow2(mmst2)*pow2(s2t))/(54.*(mmst1 - mmst2)*(-mmgl
       + mmst2)) + (185*pow2(mmst2)*pow2(s2t))/(18.*pow2(-mmgl + mmst2)) + (
      pow2(s2t)*pow3(mmst1))/(54.*(-mmgl + mmst1)*(mmst1 - mmst2)*(-mmgl +
      mmst2)) + pow3(mmst2)/(18.*(mmst1 - mmst2)*pow2(-mmgl + mmst1)) - (mmst1*
      pow3(mmst2))/(18.*pow2(-mmgl + mmst1)*pow2(mmst1 - mmst2)) + (pow2(s2t)*
      pow3(mmst2))/(54.*(-mmgl + mmst1)*pow2(mmst1 - mmst2)) - (pow2(s2t)*pow3(
      mmst2))/(54.*(-mmgl + mmst2)*pow2(mmst1 - mmst2)) - (2*mmst1*pow2(mmst2))
      /(3.*pow3(-mmgl + mmst2)) - (20*mmsusy*pow2(mmst2))/(3.*pow3(-mmgl +
      mmst2)) + (10*pow3(mmst2))/(3.*pow3(-mmgl + mmst2)) - (8*pow2(s2t)*pow3(
      mmst2))/pow3(-mmgl + mmst2) - (5*mmst2*pow3(mmsusy))/(9.*pow2(-mmgl +
      mmst2)*pow2(mmst2 - mmsusy)) + pow4(mmst2)/(18.*pow2(-mmgl + mmst1)*pow2(
      mmst1 - mmst2)) + (2*pow4(mmst2))/(3.*pow4(-mmgl + mmst2)) + (10*pow2(s2t
      )*pow4(mmst2))/(3.*pow4(-mmgl + mmst2)));

   return result * g34 * twoLoop;
}

/// 2-loop full SQCD contributions to Delta Mt over mt [hep-ph/0507139]
double dMt_over_mt_2loop(const Parameters& pars)
{
   return dMt_over_mt_2loop_qcd(pars) + dMt_over_mt_2loop_susy(pars);
}

} // namespace mssm_twoloop_mt
} // namespace flexiblesusy
