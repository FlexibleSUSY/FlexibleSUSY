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

// This file has been generated at Wed 17 Jun 2020 11:16:10
// with the script "tquark_to_cpp.m".

#include "mssm_twoloop_mt.hpp"
#include "dilog.hpp"
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
         (mm1 - mm2)*(2*dilog(1 - mm1/mm2) + pow2(log12)/2.) +
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
   const double Xt     = pars.xt;
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

   if (is_equal(mmst1, mmst2, mmt) && is_equal(mmst1, mmgl, mmt) &&
       is_equal(mmst1, mmsusy, mmt) && is_equal(std::abs(Xt), 0., 1e-1)) {
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
   -10.055555555555555 + 82*logmmtmmu - (20*logmmsusymmu*logmmtmmu)/3. + (56*
      mmst1)/(mmgl - mmst1) + (2*mmst1)/(3*mmgl - 3*mmst2) - (896*mmst1)/(9.*(
      mmst1 - mmst2)) + (6*mmst2)/(-mmgl + mmst1) + (188*mmst2)/(3*mmgl - 3*
      mmst2) + (896*mmst2)/(9.*(mmst1 - mmst2)) + logmmtmmu*(-74 + (8*mmst1)/(3
      *mmgl - 3*mmst1) + (8*mmst2)/(3*mmgl - 3*mmst2)) + (60*mmsusy)/(-mmgl +
      mmst1) + (60*mmsusy)/(-mmgl + mmst2) + (16*logmmglmmu*logmmtmmu*mgl*mmgl*
      (mmst1 - mmst2)*s2t)/(3.*(mmgl - mmst1)*(mmgl - mmst2)*mt) - (160*
      logmmsusymmu*mgl*(mmst1 - mmst2)*mmsusy*s2t)/(3.*(mmgl - mmst1)*(mmgl -
      mmst2)*mt) - (16*logmmst1mmu*logmmtmmu*mgl*mmst1*s2t)/(3*mmgl*mt - 3*
      mmst1*mt) + (16*logmmst2mmu*logmmtmmu*mgl*mmst2*s2t)/(3*mmgl*mt - 3*mmst2
      *mt) + (8*zt2)/9. + (91*mmst1*zt2)/(9*mmgl - 9*mmst1) - (128*mmst1*zt2)/(
      9.*(mmst1 - mmst2)) + (mmst2*zt2)/(-mmgl + mmst1) + (11*mmst2*zt2)/(mmgl
      - mmst2) + (128*mmst2*zt2)/(9.*(mmst1 - mmst2)) + (mmst1*zt2)/(9.*(-mmgl
      + mmst2)) + (10*mmsusy*zt2)/(-mmgl + mmst1) + (10*mmsusy*zt2)/(-mmgl +
      mmst2) + (20*mgl*(mmst1 - mmst2)*mmsusy*s2t*pow2(logmmsusymmu))/(3.*(mmgl
       - mmst1)*(mmgl - mmst2)*mt) - (2825*mmst1*mmst2)/(324.*pow2(mmgl - mmst1
      )) - (27275*mmst1*mmsusy)/(324.*pow2(mmgl - mmst1)) - (40*logmmst1mmu*
      logmmsusymmu*mgl*mmst1*mmsusy*s2t)/(3.*mt*pow2(mmgl - mmst1)) - (4*mmst1*
      mmst2*zt2)/(3.*pow2(mmgl - mmst1)) - (130*mmst1*mmsusy*zt2)/(9.*pow2(mmgl
       - mmst1)) + (40*mgl*(mmst1 - mmsusy)*s2t*fin(mmst1,mmsusy,mmu))/(3.*mt*
      pow2(mmgl - mmst1)) + (20*pow2(mmst1))/(3.*(mmgl - mmst1)*(mmgl - mmst2))
      + (896*pow2(mmst1))/(9.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (8*zt2*pow2(
      mmst1))/(9.*(mmgl - mmst1)*(mmgl - mmst2)) + (128*zt2*pow2(mmst1))/(9.*(-
      mmgl + mmst1)*(mmst1 - mmst2)) + (886*pow2(mmst1))/(9.*pow2(mmgl - mmst1)
      ) + (233*mmst2*pow2(mmst1))/(324.*(mmst1 - mmst2)*pow2(mmgl - mmst1)) + (
      1355*mmsusy*pow2(mmst1))/(324.*(mmst1 - mmsusy)*pow2(mmgl - mmst1)) + (44
      *zt2*pow2(mmst1))/(3.*pow2(mmgl - mmst1)) + (10*mmsusy*zt2*pow2(mmst1))/(
      9.*(mmst1 - mmsusy)*pow2(mmgl - mmst1)) - (2*logmmst1mmu*logmmtmmu*(-10*
      mmgl*mmst1 + pow2(mmgl) + 5*pow2(mmst1)))/(3.*pow2(mmgl - mmst1)) - (8*
      mmst1*mmst2)/pow2(mmgl - mmst2) - (27275*mmst2*mmsusy)/(324.*pow2(mmgl -
      mmst2)) + (40*logmmst2mmu*logmmsusymmu*mgl*mmst2*mmsusy*s2t)/(3.*mt*pow2(
      mmgl - mmst2)) - (4*mmst1*mmst2*zt2)/(3.*pow2(mmgl - mmst2)) - (130*mmst2
      *mmsusy*zt2)/(9.*pow2(mmgl - mmst2)) - (40*mgl*(mmst2 - mmsusy)*s2t*fin(
      mmst2,mmsusy,mmu))/(3.*mt*pow2(mmgl - mmst2)) + (40*logmmglmmu*
      logmmsusymmu*mgl*mmgl*(2*mmgl - mmst1 - mmst2)*(mmst1 - mmst2)*mmsusy*s2t
      )/(3.*mt*pow2(mmgl - mmst1)*pow2(mmgl - mmst2)) - (40*mgl*(mmst1 - mmst2)
      *(mmgl*(mmst1 + mmst2 - 2*mmsusy) + mmst2*mmsusy + mmst1*(-2*mmst2 +
      mmsusy))*s2t*fin(mmgl,mmsusy,mmu))/(3.*mt*pow2(mmgl - mmst1)*pow2(mmgl -
      mmst2)) + (8*logmmglmmu*logmmtmmu*((2*mmst1)/(-mmgl + mmst1) - pow2(mmst1
      )/pow2(mmgl - mmst1) + (8*mmgl*mmst2 - 5*pow2(mmgl) - 4*pow2(mmst2))/pow2
      (mmgl - mmst2)))/3. - (896*pow2(mmst2))/(9.*(mmst1 - mmst2)*(-mmgl +
      mmst2)) - (128*zt2*pow2(mmst2))/(9.*(mmst1 - mmst2)*(-mmgl + mmst2)) - (
      233*mmst1*pow2(mmst2))/(216.*(mmst1 - mmst2)*pow2(mmgl - mmst1)) + (886*
      pow2(mmst2))/(9.*pow2(mmgl - mmst2)) + (1355*mmsusy*pow2(mmst2))/(324.*(
      mmst2 - mmsusy)*pow2(mmgl - mmst2)) + (44*zt2*pow2(mmst2))/(3.*pow2(mmgl
      - mmst2)) + (10*mmsusy*zt2*pow2(mmst2))/(9.*(mmst2 - mmsusy)*pow2(mmgl -
      mmst2)) + (233*pow2(mmst1)*pow2(mmst2))/(648.*pow2(mmgl - mmst1)*pow2(
      mmst1 - mmst2)) + (4*mgl*(mmst1 - mmst2)*s2t*fin(mmst1,mmst2,mmu)*(-(
      mmst1*mmst2) - 5*mmgl*(mmst1 + mmst2) + 5*pow2(mmgl) + 3*pow2(mmst1) + 3*
      pow2(mmst2)))/(9.*mt*pow2(mmgl - mmst1)*pow2(mmgl - mmst2)) - (4*mgl*s2t*
      fin(mmst2,mmgl,mmu)*(-7*mmst1*mmst2 + mmgl*(-37*mmst1 + mmst2) + 18*pow2(
      mmgl) + 22*pow2(mmst1) + 3*pow2(mmst2)))/(9.*(mmgl - mmst2)*mt*pow2(mmgl
      - mmst1)) - (2*logmmst2mmu*logmmtmmu*(-10*mmgl*mmst2 + pow2(mmgl) + 5*
      pow2(mmst2)))/(3.*pow2(mmgl - mmst2)) + (4*mgl*s2t*fin(mmst1,mmgl,mmu)*(
      mmgl*(mmst1 - 37*mmst2) - 7*mmst1*mmst2 + 18*pow2(mmgl) + 3*pow2(mmst1) +
      22*pow2(mmst2)))/(9.*(mmgl - mmst1)*mt*pow2(mmgl - mmst2)) - (515*mmst1*
      pow2(mmsusy))/(108.*(mmst1 - mmsusy)*pow2(mmgl - mmst1)) - (20*mmst1*zt2*
      pow2(mmsusy))/(9.*(mmst1 - mmsusy)*pow2(mmgl - mmst1)) - (515*mmst2*pow2(
      mmsusy))/(108.*(mmst2 - mmsusy)*pow2(mmgl - mmst2)) - (20*mmst2*zt2*pow2(
      mmsusy))/(9.*(mmst2 - mmsusy)*pow2(mmgl - mmst2)) + (95*pow2(mmst1)*pow2(
      mmsusy))/(162.*pow2(mmgl - mmst1)*pow2(mmst1 - mmsusy)) + (10*zt2*pow2(
      mmst1)*pow2(mmsusy))/(9.*pow2(mmgl - mmst1)*pow2(mmst1 - mmsusy)) + (95*
      pow2(mmst2)*pow2(mmsusy))/(162.*pow2(mmgl - mmst2)*pow2(mmst2 - mmsusy))
      + (10*zt2*pow2(mmst2)*pow2(mmsusy))/(9.*pow2(mmgl - mmst2)*pow2(mmst2 -
      mmsusy)) - (256*pow2(s2t))/9. + (63625*mmst1*pow2(s2t))/(1296.*(-mmgl +
      mmst1)) + (56281*mmst1*pow2(s2t))/(1296.*(mmgl - mmst2)) + (896*mmst1*
      pow2(s2t))/(9.*(mmst1 - mmst2)) + (10535*mmst2*pow2(s2t))/(3888.*(mmgl -
      mmst1)) - (896*mmst2*pow2(s2t))/(9.*(mmst1 - mmst2)) + (mmst1*mmst2*pow2(
      s2t))/(648.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (10967*mmst2*pow2(s2t))/(
      1296.*(-mmgl + mmst2)) - (85*mmst1*mmst2*pow2(s2t))/(972.*(mmst1 - mmst2)
      *(-mmgl + mmst2)) + (41*mmst1*zt2*pow2(s2t))/(9.*(-mmgl + mmst1)) + (52*
      mmst1*zt2*pow2(s2t))/(9*mmgl - 9*mmst2) + (128*mmst1*zt2*pow2(s2t))/(9.*(
      mmst1 - mmst2)) + (11*mmst2*zt2*pow2(s2t))/(9*mmgl - 9*mmst2) - (128*
      mmst2*zt2*pow2(s2t))/(9.*(mmst1 - mmst2)) + (mmst1*mmst2*zt2*pow2(s2t))/(
      27.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (mmst1*mmst2*zt2*pow2(s2t))/(27.*(
      mmst1 - mmst2)*(-mmgl + mmst2)) + (16*mmst1*mmst2*pow2(s2t))/(9.*pow2(
      mmgl - mmst1)) + (52825*pow2(mmst1)*pow2(s2t))/(1296.*(mmgl - mmst1)*(
      mmgl - mmst2)) + (96683*pow2(mmst1)*pow2(s2t))/(972.*(mmgl - mmst1)*(
      mmst1 - mmst2)) - (85*pow2(mmst1)*pow2(s2t))/(972.*(mmst1 - mmst2)*(-mmgl
       + mmst2)) + (52*zt2*pow2(mmst1)*pow2(s2t))/(9.*(mmgl - mmst1)*(mmgl -
      mmst2)) + (383*zt2*pow2(mmst1)*pow2(s2t))/(27.*(mmgl - mmst1)*(mmst1 -
      mmst2)) - (zt2*pow2(mmst1)*pow2(s2t))/(27.*(mmst1 - mmst2)*(-mmgl + mmst2
      )) - (2*pow2(mmst1)*pow2(s2t))/pow2(mmgl - mmst1) + (16*mmst1*mmst2*pow2(
      s2t))/(9.*pow2(mmgl - mmst2)) + (85*pow2(mmst2)*pow2(s2t))/(972.*(-mmgl +
      mmst1)*(mmst1 - mmst2)) + (387233*pow2(mmst2)*pow2(s2t))/(3888.*(mmst1 -
      mmst2)*(-mmgl + mmst2)) + (zt2*pow2(mmst2)*pow2(s2t))/(27.*(-mmgl + mmst1
      )*(mmst1 - mmst2)) + (383*zt2*pow2(mmst2)*pow2(s2t))/(27.*(mmst1 - mmst2)
      *(-mmgl + mmst2)) - (2*pow2(mmst2)*pow2(s2t))/pow2(mmgl - mmst2) + (167*
      mmst1*pow2(mmst2)*pow2(s2t))/(3888.*(mmgl - mmst1)*pow2(mmst1 - mmst2)) +
      (4*mgl*(mmst1 - mmst2)*s2t*(18*mmst2 + 180*mmsusy + 3*mmst2*zt2 + 30*
      mmsusy*zt2 + mmst1*(18 + 3*zt2 - 4*pow2(s2t)) - 4*mmst2*pow2(s2t) + mmgl*
      (147 + 19*zt2 + 8*pow2(s2t))))/(9.*(mmgl - mmst1)*(mmgl - mmst2)*mt) - (8
      *logmmglmmu*mgl*mmgl*(mmst1 - mmst2)*s2t*(mmgl*(-107*mmst1 - 107*mmst2 +
      60*mmsusy) + 65*pow2(mmgl) + pow2(mmst1)*(-3 + 2*pow2(s2t)) + mmst1*(155*
      mmst2 - 30*mmsusy - 4*mmst2*pow2(s2t)) + mmst2*(-3*mmst2 - 30*mmsusy + 2*
      mmst2*pow2(s2t))))/(9.*mt*pow2(mmgl - mmst1)*pow2(mmgl - mmst2)) + (8*
      logmmst2mmu*mgl*s2t*(mmst2*(30*mmst1*mmsusy + 3*pow2(mmst1) + 2*pow2(
      mmst2)*(-3 + pow2(s2t)) + mmst1*mmst2*(-41 + 6*pow2(s2t))) + pow2(mmgl)*(
      2*mmst1*(4 + pow2(s2t)) + mmst2*(7 + 6*pow2(s2t))) - mmgl*(2*pow2(mmst1)*
      pow2(s2t) + mmst1*mmst2*(15 + 4*pow2(s2t)) + 2*mmst2*(15*mmsusy + mmst2*(
      -26 + 5*pow2(s2t)))) - 8*pow3(mmgl)))/(9.*(mmgl - mmst1)*mt*pow2(mmgl -
      mmst2)) + (8*logmmst1mmu*mgl*s2t*(mmgl*(30*mmst1*mmsusy + 2*pow2(mmst2)*
      pow2(s2t) + mmst1*mmst2*(15 + 4*pow2(s2t)) + 2*pow2(mmst1)*(-26 + 5*pow2(
      s2t))) - mmst1*(3*mmst2*(mmst2 + 10*mmsusy) + 2*pow2(mmst1)*(-3 + pow2(
      s2t)) + mmst1*mmst2*(-41 + 6*pow2(s2t))) - pow2(mmgl)*(2*mmst2*(4 + pow2(
      s2t)) + mmst1*(7 + 6*pow2(s2t))) + 8*pow3(mmgl)))/(9.*(mmgl - mmst2)*mt*
      pow2(mmgl - mmst1)) + (10*logmmst1mmu*logmmsusymmu*mmst1*(3*mmgl + mmst1)
      *mmsusy)/(3.*pow3(mmgl - mmst1)) + (10*fin(mmst1,mmsusy,mmu)*(-3*mmgl*(
      mmst1 - mmsusy) + mmst1*(mmst1 + mmsusy) - 2*pow2(mmgl)))/(3.*pow3(mmgl -
      mmst1)) + (85*pow2(s2t)*pow3(mmst1))/(972.*(mmgl - mmst1)*(mmgl - mmst2)*
      (mmst1 - mmst2)) + (zt2*pow2(s2t)*pow3(mmst1))/(27.*(mmgl - mmst1)*(mmgl
      - mmst2)*(mmst1 - mmst2)) + (112*pow3(mmst1))/(3.*pow3(mmgl - mmst1)) + (
      16*zt2*pow3(mmst1))/(3.*pow3(mmgl - mmst1)) - (2*mgl*mmst1*s2t*pow2(
      logmmst1mmu)*(-2*pow2(mmgl)*(15*mmst1*mmsusy - 15*mmst2*mmsusy + 4*pow2(
      mmst2)*(-1 + pow2(s2t)) + mmst1*mmst2*(11 + 8*pow2(s2t)) + pow2(mmst1)*(-
      31 + 12*pow2(s2t))) + 2*(mmst1*(-7 + 4*pow2(s2t)) + mmst2*(-1 + 4*pow2(
      s2t)))*pow3(mmgl) + mmst1*(-3*mmst1*mmst2*(11*mmst2 + 10*mmsusy) + 3*(
      mmst2 + 10*mmsusy)*pow2(mmst2) + mmst2*pow2(mmst1)*(43 - 16*pow2(s2t)) +
      3*pow3(mmst1)) + mmgl*(-3*(mmst2 + 10*mmsusy)*pow2(mmst2) + mmst1*pow2(
      mmst2)*(21 + 8*pow2(s2t)) + pow2(mmst1)*(-11*mmst2 + 30*mmsusy + 24*mmst2
      *pow2(s2t)) + (-55 + 16*pow2(s2t))*pow3(mmst1))))/(9.*(mmgl - mmst2)*(
      mmst1 - mmst2)*mt*pow3(-mmgl + mmst1)) + (10*fin(mmgl,mmsusy,mmu)*(1/(-
      mmgl + mmst1) + 1/(-mmgl + mmst2) - (3*mmsusy)/pow2(mmgl - mmst1) - (3*
      mmsusy)/pow2(mmgl - mmst2) + (4*pow2(mmst1))/pow3(mmgl - mmst1) + (mmst1*
      (-3*mmgl + 3*mmst1 + 4*mmsusy))/pow3(-mmgl + mmst1) + (mmst2*(3*mmgl - 4*
      mmsusy))/pow3(mmgl - mmst2) + pow2(mmst2)/pow3(mmgl - mmst2)))/3. + (fin(
      mmst2,mmgl,mmu)*(7/(-mmgl + mmst1) + 45/(mmgl - mmst2) + 128/(mmst1 -
      mmst2) + (11*pow2(s2t))/(-mmgl + mmst1) + (11*pow2(s2t))/(mmgl - mmst2) -
      (128*pow2(s2t))/(mmst1 - mmst2) + mmst2*(-128/((mmst1 - mmst2)*(-mmgl +
      mmst2)) - 9/pow2(mmgl - mmst1) + 108/pow2(mmgl - mmst2) + (128*pow2(s2t))
      /((mmst1 - mmst2)*(-mmgl + mmst2))) + (mmst1*(17*mmgl - 2*mmst1*(4 + 13*
      pow2(s2t)) + mmst2*(-9 + 26*pow2(s2t))))/((mmgl - mmst2)*pow2(mmgl -
      mmst1)) + (12*pow2(mmst1))/pow3(mmgl - mmst1) + (12*mmst1*mmst2)/pow3(-
      mmgl + mmst1) + (48*pow2(mmst2))/pow3(mmgl - mmst2)))/9. + (10*
      logmmst2mmu*logmmsusymmu*mmst2*(3*mmgl + mmst2)*mmsusy)/(3.*pow3(mmgl -
      mmst2)) + (10*fin(mmst2,mmsusy,mmu)*(-3*mmgl*(mmst2 - mmsusy) + mmst2*(
      mmst2 + mmsusy) - 2*pow2(mmgl)))/(3.*pow3(mmgl - mmst2)) + (271*pow3(
      mmst2))/(648.*(mmst1 - mmst2)*pow2(mmgl - mmst1)) + (zt2*pow3(mmst2))/(9.
      *(mmst1 - mmst2)*pow2(mmgl - mmst1)) - (7*mmst1*pow3(mmst2))/(9.*pow2(
      mmgl - mmst1)*pow2(mmst1 - mmst2)) - (mmst1*zt2*pow3(mmst2))/(9.*pow2(
      mmgl - mmst1)*pow2(mmst1 - mmst2)) + (169*pow2(s2t)*pow3(mmst2))/(1296.*(
      -mmgl + mmst1)*pow2(mmst1 - mmst2)) - (85*pow2(s2t)*pow3(mmst2))/(972.*(-
      mmgl + mmst2)*pow2(mmst1 - mmst2)) + (zt2*pow2(s2t)*pow3(mmst2))/(27.*(-
      mmgl + mmst1)*pow2(mmst1 - mmst2)) - (zt2*pow2(s2t)*pow3(mmst2))/(27.*(-
      mmgl + mmst2)*pow2(mmst1 - mmst2)) + (112*pow3(mmst2))/(3.*pow3(mmgl -
      mmst2)) + (16*zt2*pow3(mmst2))/(3.*pow3(mmgl - mmst2)) + (2*mgl*mmst2*s2t
      *pow2(logmmst2mmu)*(-2*pow2(mmgl)*(-15*mmst1*mmsusy + 15*mmst2*mmsusy + 4
      *pow2(mmst1)*(-1 + pow2(s2t)) + mmst1*mmst2*(11 + 8*pow2(s2t)) + pow2(
      mmst2)*(-31 + 12*pow2(s2t))) + 2*(mmst2*(-7 + 4*pow2(s2t)) + mmst1*(-1 +
      4*pow2(s2t)))*pow3(mmgl) + mmgl*(mmst1*pow2(mmst2)*(-11 + 24*pow2(s2t)) +
      pow2(mmst1)*(-30*mmsusy + mmst2*(21 + 8*pow2(s2t))) + pow2(mmst2)*(30*
      mmsusy + mmst2*(-55 + 16*pow2(s2t))) - 3*pow3(mmst1)) + mmst2*((-33*mmst2
       + 30*mmsusy)*pow2(mmst1) + mmst1*mmst2*(-30*mmsusy + mmst2*(43 - 16*pow2
      (s2t))) + 3*pow3(mmst1) + 3*pow3(mmst2))))/(9.*(mmgl - mmst1)*(mmst1 -
      mmst2)*mt*pow3(mmgl - mmst2)) + (fin(mmst1,mmgl,mmu)*(37/(mmgl - mmst1) +
      1/(mmgl - mmst2) - 128/(mmst1 - mmst2) + (mmst2*(9 - 26*pow2(s2t)))/pow2(
      mmgl - mmst2) + (15*pow2(s2t))/(-mmgl + mmst1) + (15*pow2(s2t))/(mmgl -
      mmst2) + (128*pow2(s2t))/(mmst1 - mmst2) + (48*pow2(mmst1))/pow3(mmgl -
      mmst1) + (12*pow2(mmst2))/pow3(mmgl - mmst2) + mmst1*(8/((mmgl - mmst1)*(
      mmgl - mmst2)) + 128/((-mmgl + mmst1)*(mmst1 - mmst2)) + 108/pow2(mmgl -
      mmst1) - 9/pow2(mmgl - mmst2) + (2*(64*mmgl + 13*mmst1 - 77*mmst2)*pow2(
      s2t))/((mmgl - mmst1)*(mmgl - mmst2)*(mmst1 - mmst2)) + (12*mmst2)/pow3(-
      mmgl + mmst2))))/9. - (95*mmst1*pow3(mmsusy))/(162.*pow2(mmgl - mmst1)*
      pow2(mmst1 - mmsusy)) - (10*mmst1*zt2*pow3(mmsusy))/(9.*pow2(mmgl - mmst1
      )*pow2(mmst1 - mmsusy)) - (95*mmst2*pow3(mmsusy))/(162.*pow2(mmgl - mmst2
      )*pow2(mmst2 - mmsusy)) - (10*mmst2*zt2*pow3(mmsusy))/(9.*pow2(mmgl -
      mmst2)*pow2(mmst2 - mmsusy)) + (10*logmmglmmu*logmmsusymmu*mmgl*mmsusy*(3
      *mmgl*(mmst1 + mmst2)*pow2(mmst1 - mmst2) + pow2(mmgl)*(6*mmst1*mmst2 - 9
      *pow2(mmst1) - 9*pow2(mmst2)) + mmst1*mmst2*(pow2(mmst1) + pow2(mmst2)) +
      8*(mmst1 + mmst2)*pow3(mmgl) - 6*pow4(mmgl)))/(3.*pow3(mmgl - mmst1)*pow3
      (mmgl - mmst2)) + (5*pow2(logmmsusymmu)*(mmst1*mmst2*(2*mmst1*mmst2 -
      mmst1*mmsusy - mmst2*mmsusy) + pow2(mmgl)*(8*mmst1*mmst2 + 5*mmst1*mmsusy
       + 5*mmst2*mmsusy + 2*pow2(mmst1) + 2*pow2(mmst2)) - mmgl*(4*mmst1*mmst2*
      (mmst2 - mmsusy) + (4*mmst2 + 3*mmsusy)*pow2(mmst1) + 3*mmsusy*pow2(mmst2
      )) - 2*(2*mmst1 + 2*mmst2 + 3*mmsusy)*pow3(mmgl) + 2*pow4(mmgl)))/(3.*
      pow2(mmgl - mmst1)*pow2(mmgl - mmst2)) - (4*logmmglmmu*logmmst1mmu*mgl*
      mmgl*s2t*(2*pow2(mmgl)*(69*mmst1*mmst2 + 26*pow2(mmst1) + 8*pow2(mmst2))
      + pow2(mmst1)*(mmst1*mmst2*(6 - 8*pow2(s2t)) + pow2(mmst1)*(-3 + 4*pow2(
      s2t)) + pow2(mmst2)*(47 + 4*pow2(s2t))) + mmgl*mmst1*(pow2(mmst1)*(1 - 4*
      pow2(s2t)) + 8*mmst1*mmst2*(-13 + pow2(s2t)) - pow2(mmst2)*(67 + 4*pow2(
      s2t))) - 2*(35*mmst1 + 16*mmst2)*pow3(mmgl) + 16*pow4(mmgl)))/(9.*mt*pow2
      (mmgl - mmst2)*pow3(mmgl - mmst1)) + (4*logmmglmmu*logmmst2mmu*mgl*mmgl*
      s2t*(2*pow2(mmgl)*(69*mmst1*mmst2 + 8*pow2(mmst1) + 26*pow2(mmst2)) +
      pow2(mmst2)*(mmst1*mmst2*(6 - 8*pow2(s2t)) + pow2(mmst2)*(-3 + 4*pow2(s2t
      )) + pow2(mmst1)*(47 + 4*pow2(s2t))) + mmgl*mmst2*(pow2(mmst2)*(1 - 4*
      pow2(s2t)) + 8*mmst1*mmst2*(-13 + pow2(s2t)) - pow2(mmst1)*(67 + 4*pow2(
      s2t))) - 2*(16*mmst1 + 35*mmst2)*pow3(mmgl) + 16*pow4(mmgl)))/(9.*mt*pow2
      (mmgl - mmst1)*pow3(mmgl - mmst2)) + (10*logmmsusymmu*(mmst1*mmst2*(19*
      mmst1*mmst2 + 12*mmst1*mmsusy + 12*mmst2*mmsusy) + pow2(mmgl)*(76*mmst1*
      mmst2 - 60*mmst1*mmsusy - 60*mmst2*mmsusy + 19*pow2(mmst1) + 19*pow2(
      mmst2)) + mmgl*(-2*mmst1*mmst2*(19*mmst2 + 24*mmsusy) + (-38*mmst2 + 36*
      mmsusy)*pow2(mmst1) + 36*mmsusy*pow2(mmst2)) + (-38*mmst1 - 38*mmst2 + 72
      *mmsusy)*pow3(mmgl) + 19*pow4(mmgl)))/(9.*pow2(mmgl - mmst1)*pow2(mmgl -
      mmst2)) + (2*mgl*mmgl*(mmst1 - mmst2)*s2t*pow2(logmmglmmu)*(-(mmst1*mmst2
      *(3*mmst2*(mmst2 + 10*mmsusy) + mmst1*(-248*mmst2 + 30*mmsusy) + 3*pow2(
      mmst1))) + pow2(mmgl)*(598*mmst1*mmst2 - 90*mmst1*mmsusy - 90*mmst2*
      mmsusy + 141*pow2(mmst1) + 141*pow2(mmst2)) - 2*(97*mmst1 + 97*mmst2 - 30
      *mmsusy)*pow3(mmgl) + 3*mmgl*(mmst1*mmst2*(-131*mmst2 + 40*mmsusy) + (-
      131*mmst2 + 10*mmsusy)*pow2(mmst1) + (mmst2 + 10*mmsusy)*pow2(mmst2) +
      pow3(mmst1)) + 46*pow4(mmgl)))/(9.*mt*pow3(mmgl - mmst1)*pow3(mmgl -
      mmst2)) + (271*pow4(mmst2))/(648.*pow2(mmgl - mmst1)*pow2(mmst1 - mmst2))
      + (zt2*pow4(mmst2))/(9.*pow2(mmgl - mmst1)*pow2(mmst1 - mmst2)) - (4*
      logmmst1mmu*logmmst2mmu*mgl*s2t*(-8*(mmst1 + mmst2)*pow2(mmgl)*(mmst1*
      mmst2*(-3 + pow2(s2t)) + pow2(mmst1)*pow2(s2t) + pow2(mmst2)*pow2(s2t)) -
      mmst1*mmst2*(mmst1 + mmst2)*(3*pow2(mmst1) + 3*pow2(mmst2) + 2*mmst1*
      mmst2*(-7 + 4*pow2(s2t))) + 4*(2*mmst1*mmst2*(-2 + pow2(s2t)) + pow2(
      mmst1)*pow2(s2t) + pow2(mmst2)*pow2(s2t))*pow3(mmgl) + 4*mmgl*(2*pow2(
      mmst1)*pow2(mmst2)*(-5 + 4*pow2(s2t)) + mmst2*(-1 + pow2(s2t))*pow3(mmst1
      ) + mmst1*(-1 + pow2(s2t))*pow3(mmst2) + pow2(s2t)*pow4(mmst1) + pow2(s2t
      )*pow4(mmst2))))/(9.*(mmst1 - mmst2)*mt*pow2(mmgl - mmst1)*pow2(mmgl -
      mmst2)) - (logmmst1mmu*logmmst2mmu*(3*pow2(mmst1)*pow2(mmst2)*(pow2(mmst1
      ) + pow2(mmst2) - 3*mmst1*mmst2*pow2(s2t)) + mmgl*mmst1*mmst2*(mmst1 +
      mmst2)*(9*pow2(mmst1) + 9*pow2(mmst2) + 2*mmst1*mmst2*(-7 + 11*pow2(s2t))
      ) + 2*(mmst1 + mmst2)*(mmst1*mmst2*(20 - 39*pow2(s2t)) + 24*pow2(mmst1)*
      pow2(s2t) + 24*pow2(mmst2)*pow2(s2t))*pow3(mmgl) - (mmst1*mmst2*(30 - 43*
      pow2(s2t)) + 48*pow2(mmst1)*pow2(s2t) + 48*pow2(mmst2)*pow2(s2t))*pow4(
      mmgl) - pow2(mmgl)*(2*pow2(mmst1)*pow2(mmst2)*(1 + 28*pow2(s2t)) + mmst2*
      (31 - 19*pow2(s2t))*pow3(mmst1) + mmst1*(31 - 19*pow2(s2t))*pow3(mmst2) +
      16*pow2(s2t)*pow4(mmst1) + 16*pow2(s2t)*pow4(mmst2)) + 16*(mmst1 + mmst2)
      *pow2(s2t)*pow5(mmgl)))/(9.*pow3(mmgl - mmst1)*pow3(mmgl - mmst2)) + (fin
      (mmst1,mmst2,mmu)*(-(mmst1*mmst2*(mmst1 + mmst2)*(3*pow2(mmst1) + 3*pow2(
      mmst2) + mmst1*mmst2*(4 - 15*pow2(s2t)))) - (mmst1 + mmst2)*pow2(mmgl)*(
      mmst1*mmst2*(94 - 68*pow2(s2t)) + pow2(mmst1)*(-29 + 11*pow2(s2t)) + pow2
      (mmst2)*(-29 + 11*pow2(s2t))) + 2*(2*mmst1*mmst2*(19 - 6*pow2(s2t)) +
      pow2(mmst1)*(-17 + 9*pow2(s2t)) + pow2(mmst2)*(-17 + 9*pow2(s2t)))*pow3(
      mmgl) - (mmst1 + mmst2)*(-14 + 29*pow2(s2t))*pow4(mmgl) - mmgl*(2*pow2(
      mmst1)*pow2(mmst2)*(-17 + 45*pow2(s2t)) + 2*mmst2*(-13 + 2*pow2(s2t))*
      pow3(mmst1) + 2*mmst1*(-13 + 2*pow2(s2t))*pow3(mmst2) + 9*pow4(mmst1) + 9
      *pow4(mmst2)) + 2*(-6 + 11*pow2(s2t))*pow5(mmgl)))/(9.*pow3(mmgl - mmst1)
      *pow3(mmgl - mmst2)) - (logmmst2mmu*(mmst1*(mmst1 - mmst2)*pow2(mmst2)*(
      60*mmst1*mmsusy + 6*pow2(mmst1) + mmst1*mmst2*(165 - 98*pow2(s2t)) + 4*
      pow2(mmst2)*(3 + 2*pow2(s2t))) + pow3(mmgl)*(mmst1*mmst2*(180*mmsusy +
      mmst2*(239 - 704*pow2(s2t))) - 5*pow2(mmst1)*(mmst2 + 46*mmst2*pow2(s2t))
      + pow2(mmst2)*(-180*mmsusy + mmst2*(489 + 70*pow2(s2t))) + (45 + 96*pow2(
      s2t))*pow3(mmst1)) - (mmst1*mmst2*(159 - 482*pow2(s2t)) + pow2(mmst2)*(
      263 - 174*pow2(s2t)) + 18*pow2(mmst1)*(5 + 8*pow2(s2t)))*pow4(mmgl) -
      pow2(mmgl)*(mmst1*pow2(mmst2)*(-420*mmsusy + mmst2*(1195 - 284*pow2(s2t))
      ) + mmst2*pow2(mmst1)*(360*mmsusy - mmst2*(721 + 320*pow2(s2t))) + 3*
      mmst2*(-5 + 4*pow2(s2t))*pow3(mmst1) + (60*mmsusy + mmst2*(53 + 64*pow2(
      s2t)))*pow3(mmst2) + 16*pow2(s2t)*pow4(mmst1)) + mmgl*mmst2*(2*mmst1*pow2
      (mmst2)*(60*mmsusy + mmst2*(199 - 69*pow2(s2t))) + mmst2*pow2(mmst1)*(-
      300*mmsusy + mmst2*(259 - 64*pow2(s2t))) + (180*mmsusy + mmst2*(-511 + 50
      *pow2(s2t)))*pow3(mmst1) + 18*pow4(mmst1) + 12*(-3 + 2*pow2(s2t))*pow4(
      mmst2)) + (mmst2*(83 - 192*pow2(s2t)) + mmst1*(45 + 64*pow2(s2t)))*pow5(
      mmgl)))/(9.*(mmst1 - mmst2)*pow2(mmgl - mmst1)*pow3(mmgl - mmst2)) - (
      logmmst1mmu*((mmst1 - mmst2)*mmst2*pow2(mmst1)*(6*mmst2*(mmst2 + 10*
      mmsusy) + mmst1*mmst2*(165 - 98*pow2(s2t)) + 4*pow2(mmst1)*(3 + 2*pow2(
      s2t))) - pow3(mmgl)*(pow2(mmst1)*(239*mmst2 - 180*mmsusy - 704*mmst2*pow2
      (s2t)) - 5*mmst1*mmst2*(mmst2 - 36*mmsusy + 46*mmst2*pow2(s2t)) + (489 +
      70*pow2(s2t))*pow3(mmst1) + 3*(15 + 32*pow2(s2t))*pow3(mmst2)) + (mmst1*
      mmst2*(159 - 482*pow2(s2t)) + pow2(mmst1)*(263 - 174*pow2(s2t)) + 18*pow2
      (mmst2)*(5 + 8*pow2(s2t)))*pow4(mmgl) - mmgl*mmst1*(mmst2*pow2(mmst1)*(
      259*mmst2 + 120*mmsusy - 64*mmst2*pow2(s2t)) + mmst1*pow2(mmst2)*(-511*
      mmst2 - 300*mmsusy + 50*mmst2*pow2(s2t)) + 2*mmst2*(199 - 69*pow2(s2t))*
      pow3(mmst1) + 18*(mmst2 + 10*mmsusy)*pow3(mmst2) + 12*(-3 + 2*pow2(s2t))*
      pow4(mmst1)) + pow2(mmgl)*(3*mmst1*pow2(mmst2)*(-5*mmst2 + 120*mmsusy + 4
      *mmst2*pow2(s2t)) - mmst2*pow2(mmst1)*(721*mmst2 + 420*mmsusy + 320*mmst2
      *pow2(s2t)) + (1195*mmst2 + 60*mmsusy - 284*mmst2*pow2(s2t))*pow3(mmst1)
      + (53 + 64*pow2(s2t))*pow4(mmst1) + 16*pow2(s2t)*pow4(mmst2)) + (-(mmst2*
      (45 + 64*pow2(s2t))) + mmst1*(-83 + 192*pow2(s2t)))*pow5(mmgl)))/(9.*(
      mmst1 - mmst2)*pow2(mmgl - mmst2)*pow3(mmgl - mmst1)) + (pow2(logmmst1mmu
      )*(-((mmst1 - mmst2)*mmst2*(3*mmst2*(mmst2 + 10*mmsusy) + 3*pow2(mmst1) +
      mmst1*mmst2*(32 - 51*pow2(s2t)))*pow3(mmst1)) + mmst1*pow3(mmgl)*(pow2(
      mmst1)*(945*mmst2 - 60*mmsusy - 762*mmst2*pow2(s2t)) + mmst1*mmst2*(83*
      mmst2 - 120*mmsusy - 342*mmst2*pow2(s2t)) + pow2(mmst2)*(7*mmst2 + 180*
      mmsusy - 50*mmst2*pow2(s2t)) - 7*(-35 + 18*pow2(s2t))*pow3(mmst1)) + (2*
      mmst1*mmst2*(-53*mmst2 - 45*mmsusy + 127*mmst2*pow2(s2t)) + pow2(mmst1)*(
      -676*mmst2 + 90*mmsusy + 685*mmst2*pow2(s2t)) + (-492 + 341*pow2(s2t))*
      pow3(mmst1) - 6*pow3(mmst2))*pow4(mmgl) + mmgl*pow2(mmst1)*(mmst2*pow2(
      mmst1)*(229*mmst2 + 60*mmsusy - 128*mmst2*pow2(s2t)) + mmst1*pow2(mmst2)*
      (-183*mmst2 - 120*mmsusy + 76*mmst2*pow2(s2t)) + mmst2*(85 - 76*pow2(s2t)
      )*pow3(mmst1) + 6*(mmst2 + 10*mmsusy)*pow3(mmst2) - 9*pow4(mmst1)) +
      mmst1*pow2(mmgl)*(mmst1*pow2(mmst2)*(141*mmst2 - 30*mmsusy - 35*mmst2*
      pow2(s2t)) + 3*mmst2*pow2(mmst1)*(-51*mmst2 + 50*mmsusy + 115*mmst2*pow2(
      s2t)) + (-623*mmst2 - 30*mmsusy + 305*mmst2*pow2(s2t))*pow3(mmst1) - 9*(
      mmst2 + 10*mmsusy)*pow3(mmst2) + (4 + 25*pow2(s2t))*pow4(mmst1)) + (12*
      pow2(mmst2) + mmst1*mmst2*(242 - 332*pow2(s2t)) + pow2(mmst1)*(386 - 308*
      pow2(s2t)))*pow5(mmgl) - 2*(3*mmst2 + mmst1*(61 - 64*pow2(s2t)))*pow6(
      mmgl)))/(18.*(mmst1 - mmst2)*pow2(mmgl - mmst2)*pow4(mmgl - mmst1)) + (
      pow2(logmmst2mmu)*(mmst2*pow3(mmgl)*(pow2(mmst2)*(60*mmsusy + 7*mmst2*(-
      35 + 18*pow2(s2t))) + 3*mmst1*mmst2*(40*mmsusy + mmst2*(-315 + 254*pow2(
      s2t))) + pow2(mmst1)*(-180*mmsusy + mmst2*(-83 + 342*pow2(s2t))) + (-7 +
      50*pow2(s2t))*pow3(mmst1)) - mmst1*(mmst1 - mmst2)*(30*mmst1*mmsusy + 3*
      pow2(mmst1) + 3*pow2(mmst2) + mmst1*mmst2*(32 - 51*pow2(s2t)))*pow3(mmst2
      ) + (mmst1*mmst2*(90*mmsusy + mmst2*(676 - 685*pow2(s2t))) + pow2(mmst2)*
      (-90*mmsusy + mmst2*(492 - 341*pow2(s2t))) + 2*mmst2*pow2(mmst1)*(53 -
      127*pow2(s2t)) + 6*pow3(mmst1))*pow4(mmgl) + mmst2*pow2(mmgl)*(mmst1*pow2
      (mmst2)*(-150*mmsusy + mmst2*(623 - 305*pow2(s2t))) + 3*mmst2*pow2(mmst1)
      *(10*mmsusy + mmst2*(51 - 115*pow2(s2t))) + (90*mmsusy + mmst2*(-141 + 35
      *pow2(s2t)))*pow3(mmst1) + (30*mmsusy - mmst2*(4 + 25*pow2(s2t)))*pow3(
      mmst2) + 9*pow4(mmst1)) + mmgl*pow2(mmst2)*(mmst1*pow2(mmst2)*(-60*mmsusy
       + mmst2*(-85 + 76*pow2(s2t))) + mmst2*pow2(mmst1)*(120*mmsusy + mmst2*(-
      229 + 128*pow2(s2t))) + (-60*mmsusy + mmst2*(183 - 76*pow2(s2t)))*pow3(
      mmst1) - 6*pow4(mmst1) + 9*pow4(mmst2)) - 2*(6*pow2(mmst1) + mmst1*mmst2*
      (121 - 166*pow2(s2t)) + pow2(mmst2)*(193 - 154*pow2(s2t)))*pow5(mmgl) + 2
      *(3*mmst1 + mmst2*(61 - 64*pow2(s2t)))*pow6(mmgl)))/(18.*(mmst1 - mmst2)*
      pow2(mmgl - mmst1)*pow4(mmgl - mmst2)) + (logmmglmmu*logmmst1mmu*mmgl*(-(
      mmst2*(3*pow2(mmst1) + pow2(mmst2)*(23 - 21*pow2(s2t)) + mmst1*mmst2*(4 +
      21*pow2(s2t)))*pow3(mmst1)) + pow3(mmgl)*(-2*mmst2*pow2(mmst1)*(139 + 92*
      pow2(s2t)) - mmst1*pow2(mmst2)*(587 + 514*pow2(s2t)) + (-35 + 122*pow2(
      s2t))*pow3(mmst1) - 64*(1 + pow2(s2t))*pow3(mmst2)) + mmgl*pow2(mmst1)*(
      81*mmst1*pow2(mmst2) + 2*mmst2*pow2(mmst1)*(8 + 3*pow2(s2t)) + (-9 + 16*
      pow2(s2t))*pow3(mmst1) - 2*(49 + 43*pow2(s2t))*pow3(mmst2)) + mmst1*pow2(
      mmgl)*(mmst2*pow2(mmst1)*(-95 + 49*pow2(s2t)) + mmst1*pow2(mmst2)*(282 +
      163*pow2(s2t)) - 9*(-4 + 9*pow2(s2t))*pow3(mmst1) + (197 + 189*pow2(s2t))
      *pow3(mmst2)) + (pow2(mmst1)*(142 - 53*pow2(s2t)) + 192*pow2(mmst2)*(1 +
      pow2(s2t)) + 3*mmst1*mmst2*(196 + 167*pow2(s2t)))*pow4(mmgl) - 2*(96*
      mmst2*(1 + pow2(s2t)) + mmst1*(105 + 64*pow2(s2t)))*pow5(mmgl) + 64*(1 +
      pow2(s2t))*pow6(mmgl)))/(9.*pow3(mmgl - mmst2)*pow4(mmgl - mmst1)) + (
      logmmglmmu*logmmst2mmu*mmgl*(mmst1*(-3*pow2(mmst2) + pow2(mmst1)*(-23 +
      21*pow2(s2t)) - mmst1*mmst2*(4 + 21*pow2(s2t)))*pow3(mmst2) - pow3(mmgl)*
      (2*mmst1*pow2(mmst2)*(139 + 92*pow2(s2t)) + mmst2*pow2(mmst1)*(587 + 514*
      pow2(s2t)) + 64*(1 + pow2(s2t))*pow3(mmst1) + (35 - 122*pow2(s2t))*pow3(
      mmst2)) + mmst2*pow2(mmgl)*(mmst1*pow2(mmst2)*(-95 + 49*pow2(s2t)) +
      mmst2*pow2(mmst1)*(282 + 163*pow2(s2t)) + (197 + 189*pow2(s2t))*pow3(
      mmst1) + 9*(4 - 9*pow2(s2t))*pow3(mmst2)) + mmgl*pow2(mmst2)*(81*mmst2*
      pow2(mmst1) + 2*mmst1*pow2(mmst2)*(8 + 3*pow2(s2t)) - 2*(49 + 43*pow2(s2t
      ))*pow3(mmst1) + (-9 + 16*pow2(s2t))*pow3(mmst2)) + (pow2(mmst2)*(142 -
      53*pow2(s2t)) + 192*pow2(mmst1)*(1 + pow2(s2t)) + 3*mmst1*mmst2*(196 +
      167*pow2(s2t)))*pow4(mmgl) - 2*(96*mmst1*(1 + pow2(s2t)) + mmst2*(105 +
      64*pow2(s2t)))*pow5(mmgl) + 64*(1 + pow2(s2t))*pow6(mmgl)))/(9.*pow3(mmgl
       - mmst1)*pow4(mmgl - mmst2)) - (2*logmmglmmu*(3*mmgl*mmst1*mmst2*((123*
      mmst2 + 10*mmsusy)*pow2(mmst1) + 123*mmst1*pow2(mmst2) + (mmst2 + 10*
      mmsusy)*pow2(mmst2) + pow3(mmst1)) + pow3(mmgl)*(pow2(mmst1)*(-270*mmsusy
       + mmst2*(3236 - 505*pow2(s2t))) + mmst1*mmst2*(180*mmsusy + mmst2*(3236
      - 505*pow2(s2t))) + pow2(mmst2)*(-270*mmsusy + mmst2*(172 + 121*pow2(s2t)
      )) + (172 + 121*pow2(s2t))*pow3(mmst1)) - 78*pow3(mmst1)*pow3(mmst2) + (
      240*mmst1*mmsusy + 240*mmst2*mmsusy + pow2(mmst1)*(-999 + 85*pow2(s2t)) +
      pow2(mmst2)*(-999 + 85*pow2(s2t)) + mmst1*mmst2*(-4944 + 982*pow2(s2t)))*
      pow4(mmgl) - pow2(mmgl)*(2*mmst2*pow2(mmst1)*(45*mmsusy + mmst2*(904 -
      207*pow2(s2t))) + mmst1*pow2(mmst2)*(90*mmsusy + mmst2*(753 + 103*pow2(
      s2t))) + (-90*mmsusy + mmst2*(753 + 103*pow2(s2t)))*pow3(mmst1) + (-90*
      mmsusy + mmst2*(-9 + 8*pow2(s2t)))*pow3(mmst2) + (-9 + 8*pow2(s2t))*pow4(
      mmst1)) - 4*(45*mmsusy + mmst1*(-421 + 96*pow2(s2t)) + mmst2*(-421 + 96*
      pow2(s2t)))*pow5(mmgl) + 12*(-51 + 16*pow2(s2t))*pow6(mmgl)))/(9.*pow3(
      mmgl - mmst1)*pow3(mmgl - mmst2)) + (pow2(logmmglmmu)*(3*mmgl*pow2(mmst1)
      *pow2(mmst2)*((-87*mmst2 + 10*mmsusy)*pow2(mmst1) - 87*mmst1*pow2(mmst2)
      + (mmst2 + 10*mmsusy)*pow2(mmst2) + pow3(mmst1)) + mmst1*mmst2*pow2(mmgl)
      *(-(mmst1*pow2(mmst2)*(220*mmst2 + 120*mmsusy + 159*mmst2*pow2(s2t))) + 2
      *mmst2*pow2(mmst1)*(288*mmst2 - 60*mmsusy + 223*mmst2*pow2(s2t)) + (-220*
      mmst2 + 60*mmsusy - 159*mmst2*pow2(s2t))*pow3(mmst1) + 6*(mmst2 + 10*
      mmsusy)*pow3(mmst2) + 6*pow4(mmst1)) - pow4(mmgl)*(2*mmst1*pow2(mmst2)*(-
      120*mmsusy + 7*mmst2*(218 + 9*pow2(s2t))) - 2*mmst2*pow2(mmst1)*(120*
      mmsusy + mmst2*(-4172 + 1249*pow2(s2t))) + 2*(-180*mmsusy + 7*mmst2*(218
      + 9*pow2(s2t)))*pow3(mmst1) + (-360*mmsusy + mmst2*(30 + 163*pow2(s2t)))*
      pow3(mmst2) + (30 + 163*pow2(s2t))*pow4(mmst1)) + 72*pow4(mmst1)*pow4(
      mmst2) + (pow2(mmst1)*(8613*mmst2 - 510*mmsusy - 1486*mmst2*pow2(s2t)) +
      mmst1*mmst2*(8613*mmst2 - 480*mmsusy - 1486*mmst2*pow2(s2t)) + pow2(mmst2
      )*(597*mmst2 - 510*mmsusy + 206*mmst2*pow2(s2t)) + (597 + 206*pow2(s2t))*
      pow3(mmst1))*pow5(mmgl) + pow3(mmgl)*(2*pow2(mmst1)*pow2(mmst2)*(1046*
      mmst2 + 180*mmsusy - 383*mmst2*pow2(s2t)) - 2*mmst2*(-1046*mmst2 + 120*
      mmsusy + 383*mmst2*pow2(s2t))*pow3(mmst1) + mmst1*(451*mmst2 - 240*mmsusy
       + 382*mmst2*pow2(s2t))*pow3(mmst2) + (451*mmst2 - 90*mmsusy + 382*mmst2*
      pow2(s2t))*pow4(mmst1) - 9*(mmst2 + 10*mmsusy)*pow4(mmst2) - 9*pow5(mmst1
      )) + (7*pow2(mmst1)*(-296 + 39*pow2(s2t)) + 7*mmst2*(-296*mmst2 + 60*
      mmsusy + 39*mmst2*pow2(s2t)) + 6*mmst1*(-1358*mmst2 + 70*mmsusy + 229*
      mmst2*pow2(s2t)))*pow6(mmgl) - 6*(-347*mmst2 + 30*mmsusy + 64*mmst2*pow2(
      s2t) + mmst1*(-347 + 64*pow2(s2t)))*pow7(mmgl) + 4*(-139 + 32*pow2(s2t))*
      pow8(mmgl)))/(18.*pow4(mmgl - mmst1)*pow4(mmgl - mmst2));

   return result * g34 * twoLoop;
}

/// 2-loop full SQCD contributions to Delta Mt over mt [hep-ph/0507139]
double dMt_over_mt_2loop(const Parameters& pars)
{
   return dMt_over_mt_2loop_qcd(pars) + dMt_over_mt_2loop_susy(pars);
}

} // namespace mssm_twoloop_mt
} // namespace flexiblesusy
