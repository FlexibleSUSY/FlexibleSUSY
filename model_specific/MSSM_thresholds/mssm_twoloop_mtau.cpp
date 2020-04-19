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

// This file has been generated at Sun 19 Apr 2020 19:46:31
// with the script "tau_to_cpp.m".

#include "mssm_twoloop_mtau.hpp"
#include "dilog.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <ostream>

namespace flexiblesusy {
namespace mssm_twoloop_mtau {

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
         (-mm1 - mm2)*(7 + power2(Pi)/6) +
         (mm1 - mm2)*(2*dilog(1 - mm1/mm2) +
            power2(log12)/2) +
         ((mm1 + mm2)*power2(log12))/2 -
         2*(mm1*power2(log1u) + mm2*power2(log2u)))/2;
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

} // anonymous namespace

double delta_mtau_2loop_atau_atau(const Parameters& pars)
{
   const Real ytau  = pars.ytau;
   const Real xtau  = pars.xtau;
   const Real mstau1 = pars.mstau1;
   const Real mstau12 = power2(mstau1);
   const Real mstau2  = pars.mstau2;
   const Real mstau22 = power2(pars.mstau2);
   const Real mstau24 = power4(pars.mstau2);
   const Real msntau2 = power2(pars.msntau);
   const Real mw2   = power2(pars.mw);
   const Real mz2   = power2(pars.mz);
   const Real mh2   = power2(pars.mh);
   const Real mH2   = power2(pars.mH);
   const Real mC2   = power2(pars.mC);
   const Real mA2   = power2(pars.mA);
   const Real mu    = pars.mu;
   const Real mu2   = power2(pars.mu);
   const Real tb    = pars.tb;
   const Real sb    = tb / std::sqrt(1 + power2(tb));
   const Real cb    = 1  / std::sqrt(1 + power2(tb));
   const Real Q2    = power2(pars.Q);
   const Real alpha = calc_alpha(mh2, mH2, tb);
   const Real sa    = std::sin(alpha);
   const Real ca    = std::cos(alpha);
   const Real Al    = xtau + mu*tb;

   const Real invdmstau     = 1/(mstau12 - mstau22);
   const Real invdmstau1mu  = 1/(-mstau12 + mu2);
   const Real invdmstau2mu  = 1/(-mstau22 + mu2);
   const Real invdmsntau2mu = 1/(-msntau2 + mu2);
   const Real invdmhH       = 1/(-mh2 + mH2);
   const Real invdmAh       = 1/(-mA2 + mh2);
   const Real invdmAH       = 1/(-mA2 + mH2);
   const Real invdmAC       = 1/(-mA2 + mC2);
   const Real invdmCh       = 1/(mC2 - mh2);
   const Real invdmCH       = 1/(mC2 - mH2);

   const Real logmstau12Q2  = std::log(mstau12/Q2);
   const Real logmstau22Q2  = std::log(mstau22/Q2);
   const Real logmH2Q2      = std::log(mH2/Q2);
   const Real logmA2Q2      = std::log(mA2/Q2);
   const Real logmh2Q2      = std::log(mh2/Q2);
   const Real logmu2Q2      = std::log(mu2/Q2);
   const Real logmC2Q2      = std::log(mC2/Q2);
   const Real logmw2Q2      = std::log(mw2/Q2);
   const Real logmz2Q2      = std::log(mz2/Q2);
   const Real logmsntau2Q2  = std::log(msntau2/Q2);

   const double result =
   Fin3(mA2,msntau2,mu2,Q2)*((3*mA2*msntau2 - 2*mA2*mu2 + 5*msntau2*mu2 - 5*
     power2(msntau2))/(16.*msntau2*power2(mA2)) + (mu2*DeltaInv(mA2,msntau2,mu2
     )*(-3*mA2*msntau2 - 2*mA2*mu2 - msntau2*mu2 + power2(mA2) + power2(mu2)))/
     (8.*mA2*msntau2)) + Fin3(mC2,msntau2,mu2,Q2)*(-((msntau2 + mu2 + 2*
     invdmstau2mu*msntau2*mu2)*power2(Al)*power2(invdmstau2mu))/(4.*msntau2) +
     (mu2*DeltaInv(mC2,msntau2,mu2)*power2(Al)*power2(invdmstau2mu)*(-3*mC2*
     msntau2 - 2*mC2*mu2 - msntau2*mu2 + power2(mC2) + power2(mu2)))/(4.*
     msntau2)) + Fin3(mstau22,mC2,mu2,Q2)*(-((mstau22 + mu2 + 2*invdmsntau2mu*
     mstau22*mu2)*power2(Al)*power2(invdmsntau2mu))/(4.*mstau22) - (mC2*mstau22
      + 2*mC2*mu2 - mstau22*mu2 + power2(mstau22))/(8.*mstau22*power2(mC2)) +
     DeltaInv(mstau22,mC2,mu2)*((mu2*(-3*mC2*mstau22 - 2*mC2*mu2 - mstau22*mu2
     + power2(mC2) + power2(mu2)))/(4.*mC2*mstau22) + (mu2*power2(Al)*power2(
     invdmsntau2mu)*(-3*mC2*mstau22 - 2*mC2*mu2 - mstau22*mu2 + power2(mC2) +
     power2(mu2)))/(4.*mstau22))) + Fin3(msntau2,mw2,mu2,Q2)*(-(mu2*(msntau2 +
     mu2 + 2*invdmstau2mu*msntau2*mu2)*power2(invdmstau2mu))/(4.*msntau2) + (
     DeltaInv(msntau2,mw2,mu2)*power2(invdmstau2mu)*power2(mu2)*(-(msntau2*(mu2
      + 3*mw2)) + power2(mu2 - mw2)))/(4.*msntau2)) + Fin3(mstau22,mw2,mu2,Q2)*
     (-(mu2*(mstau22 + mu2 + 2*invdmsntau2mu*mstau22*mu2)*power2(invdmsntau2mu)
     )/(4.*mstau22) + (DeltaInv(mstau22,mw2,mu2)*power2(invdmsntau2mu)*power2(
     mu2)*(-(mstau22*(mu2 + 3*mw2)) + power2(mu2 - mw2)))/(4.*mstau22)) + (
     DeltaInv(msntau2,mw2,mu2)*power2(invdmstau2mu)*power2(mu2)*(power2(mu2 -
     mw2)*(mu2*(42 + 6*logmsntau2Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*power2(
     logmsntau2Q2) + 3*power2(logmu2Q2) + power2(Pi)) + mw2*(-6*(3 + logmu2Q2)*
     logmw2Q2 + 6*logmsntau2Q2*(-9 + logmu2Q2 + 2*logmw2Q2) + 9*power2(
     logmsntau2Q2) + 3*power2(logmw2Q2) + 2*(42 + power2(Pi)))) - msntau2*(
     power2(mu2)*(42 + 6*logmsntau2Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*power2(
     logmsntau2Q2) + 3*power2(logmu2Q2) + power2(Pi)) + mu2*mw2*(294 - 18*
     logmw2Q2 - 18*logmu2Q2*(3 + logmw2Q2) + 12*logmsntau2Q2*(-15 + 3*logmu2Q2
     + 2*logmw2Q2) + 30*power2(logmsntau2Q2) + 9*power2(logmu2Q2) + 3*power2(
     logmw2Q2) + 7*power2(Pi)) + power2(mw2)*(-6*(9 + logmu2Q2)*logmw2Q2 + 6*
     logmsntau2Q2*(-15 + logmu2Q2 + 4*logmw2Q2) + 15*power2(logmsntau2Q2) + 9*
     power2(logmw2Q2) + 4*(42 + power2(Pi))))))/(24.*msntau2) + (DeltaInv(
     mstau22,mw2,mu2)*power2(invdmsntau2mu)*power2(mu2)*(power2(mu2 - mw2)*(mu2
     *(42 + 6*logmstau22Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*power2(
     logmstau22Q2) + 3*power2(logmu2Q2) + power2(Pi)) + mw2*(-6*(3 + logmu2Q2)*
     logmw2Q2 + 6*logmstau22Q2*(-9 + logmu2Q2 + 2*logmw2Q2) + 9*power2(
     logmstau22Q2) + 3*power2(logmw2Q2) + 2*(42 + power2(Pi)))) - mstau22*(
     power2(mu2)*(42 + 6*logmstau22Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*power2(
     logmstau22Q2) + 3*power2(logmu2Q2) + power2(Pi)) + mu2*mw2*(294 - 18*
     logmw2Q2 - 18*logmu2Q2*(3 + logmw2Q2) + 12*logmstau22Q2*(-15 + 3*logmu2Q2
     + 2*logmw2Q2) + 30*power2(logmstau22Q2) + 9*power2(logmu2Q2) + 3*power2(
     logmw2Q2) + 7*power2(Pi)) + power2(mw2)*(-6*(9 + logmu2Q2)*logmw2Q2 + 6*
     logmstau22Q2*(-15 + logmu2Q2 + 4*logmw2Q2) + 15*power2(logmstau22Q2) + 9*
     power2(logmw2Q2) + 4*(42 + power2(Pi))))))/(24.*mstau22) + (DeltaInv(
     mstau12,mz2,mu2)*power2(invdmstau2mu)*power2(mu2)*(power2(mu2 - mz2)*(mu2*
     (42 + 6*logmstau12Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*power2(logmstau12Q2
     ) + 3*power2(logmu2Q2) + power2(Pi)) + mz2*(-6*(3 + logmu2Q2)*logmz2Q2 + 6
     *logmstau12Q2*(-9 + logmu2Q2 + 2*logmz2Q2) + 9*power2(logmstau12Q2) + 3*
     power2(logmz2Q2) + 2*(42 + power2(Pi)))) - mstau12*(power2(mu2)*(42 + 6*
     logmstau12Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*power2(logmstau12Q2) + 3*
     power2(logmu2Q2) + power2(Pi)) + mu2*mz2*(294 - 18*logmz2Q2 - 18*logmu2Q2*
     (3 + logmz2Q2) + 12*logmstau12Q2*(-15 + 3*logmu2Q2 + 2*logmz2Q2) + 30*
     power2(logmstau12Q2) + 9*power2(logmu2Q2) + 3*power2(logmz2Q2) + 7*power2(
     Pi)) + power2(mz2)*(-6*(9 + logmu2Q2)*logmz2Q2 + 6*logmstau12Q2*(-15 +
     logmu2Q2 + 4*logmz2Q2) + 15*power2(logmstau12Q2) + 9*power2(logmz2Q2) + 4*
     (42 + power2(Pi))))))/(48.*mstau12) + (DeltaInv(mstau22,mz2,mu2)*power2(
     invdmstau1mu)*power2(mu2)*(power2(mu2 - mz2)*(mu2*(42 + 6*logmstau22Q2*(-3
      + logmu2Q2) - 18*logmu2Q2 + 3*power2(logmstau22Q2) + 3*power2(logmu2Q2) +
     power2(Pi)) + mz2*(-6*(3 + logmu2Q2)*logmz2Q2 + 6*logmstau22Q2*(-9 +
     logmu2Q2 + 2*logmz2Q2) + 9*power2(logmstau22Q2) + 3*power2(logmz2Q2) + 2*(
     42 + power2(Pi)))) - mstau22*(power2(mu2)*(42 + 6*logmstau22Q2*(-3 +
     logmu2Q2) - 18*logmu2Q2 + 3*power2(logmstau22Q2) + 3*power2(logmu2Q2) +
     power2(Pi)) + mu2*mz2*(294 - 18*logmz2Q2 - 18*logmu2Q2*(3 + logmz2Q2) + 12
     *logmstau22Q2*(-15 + 3*logmu2Q2 + 2*logmz2Q2) + 30*power2(logmstau22Q2) +
     9*power2(logmu2Q2) + 3*power2(logmz2Q2) + 7*power2(Pi)) + power2(mz2)*(-6*
     (9 + logmu2Q2)*logmz2Q2 + 6*logmstau22Q2*(-15 + logmu2Q2 + 4*logmz2Q2) +
     15*power2(logmstau22Q2) + 9*power2(logmz2Q2) + 4*(42 + power2(Pi))))))/(
     48.*mstau22) + Fin3(mH2,msntau2,mu2,Q2)*(((3*msntau2*(-msntau2 + mu2) +
     mH2*(5*msntau2 + 2*mu2))*(-1 + power2(sa)))/(16.*msntau2*power2(mH2)) - (
     mu2*DeltaInv(mH2,msntau2,mu2)*(-3*mH2*msntau2 - 2*mH2*mu2 - msntau2*mu2 +
     power2(mH2) + power2(mu2))*(-1 + power2(sa)))/(8.*mH2*msntau2)) + (Fin20(
     mA2,mh2,Q2)*(-invdmAh + 1/mA2 + 10/mh2 - (5*mh2)/power2(mA2) + (4*mA2)/
     power2(mh2) + ((-mA2 + mh2)*power2(invdmAh)*power2(mA2))/power2(mh2))*
     power2(sa))/16. + Fin3(mstau12,mstau12,mh2,Q2)*((invdmstau1mu*sa*(cb*sa -
     ca*invdmstau*mu2*sb))/(2.*cb) + (Al*invdmstau*invdmstau1mu*mu*sa*(ca*cb -
     sa*sb))/(2.*cb) + (invdmstau*invdmstau1mu*power2(Al)*power2(sa))/2.) +
     Fin3(mstau22,mstau22,mh2,Q2)*((invdmstau2mu*sa*(cb*sa + ca*invdmstau*mu2*
     sb))/(2.*cb) + (Al*invdmstau*invdmstau2mu*mu*sa*(-(ca*cb) + sa*sb))/(2.*cb
     ) - (invdmstau*invdmstau2mu*power2(Al)*power2(sa))/2.) + (Fin20(msntau2,
     mu2,Q2)*(2*invdmsntau2mu + (8*(-1 + invdmsntau2mu*mu2))/mA2 + (5*(msntau2
     - mu2))/power2(mA2) + (8*(-1 + invdmsntau2mu*mu2)*(-1 + power2(sa)))/mH2 +
     (3*(msntau2 - mu2)*(-1 + power2(sa)))/power2(mH2) - (8*(-1 + invdmsntau2mu
     *mu2)*power2(sa))/mh2 + (3*(-msntau2 + mu2)*power2(sa))/power2(mh2)))/16.
     + Fin3(msntau2,mh2,mu2,Q2)*(-((3*msntau2*(-msntau2 + mu2) + mh2*(5*msntau2
      + 2*mu2))*power2(sa))/(16.*msntau2*power2(mh2)) + (mu2*DeltaInv(msntau2,
     mh2,mu2)*(-3*mh2*msntau2 - 2*mh2*mu2 - msntau2*mu2 + power2(mh2) + power2(
     mu2))*power2(sa))/(8.*mh2*msntau2)) + Fin3(mH2,mstau12,mstau12,Q2)*(-(
     invdmstau*invdmstau1mu*power2(Al)*(-1 + power2(sa)))/2. + (invdmstau1mu*(
     cb + ca*invdmstau*mu2*sa*sb - cb*power2(sa)))/(2.*cb) - (Al*invdmstau*
     invdmstau1mu*mu*(ca*cb*sa + sb - sb*power2(sa)))/(2.*cb)) + Fin3(mstau22,
     mstau22,mH2,Q2)*(-(invdmstau2mu*(ca*invdmstau*mu2*sa*sb + cb*(-1 + power2(
     sa))))/(2.*cb) + (invdmstau*invdmstau2mu*power2(Al)*(-1 + power2(sa)))/2.
     + (Al*invdmstau*invdmstau2mu*mu*(ca*cb*sa + sb - sb*power2(sa)))/(2.*cb))
     + Fin3(mstau22,mH2,mu2,Q2)*((Al*mu*(ca*cb*sa*(-2*invdmstau*(mstau22 +
     invdmstau2mu*mH2*mstau22) + mH2*(mstau22 + mu2 + 2*invdmstau1mu*mstau22*
     mu2)*power2(invdmstau1mu)) + 2*invdmstau*(1 + invdmstau2mu*mH2)*mstau22*sb
     *(-1 + power2(sa))))/(4.*cb*mH2*mstau22) - (power2(Al)*(4*invdmstau*(
     mstau22 + invdmstau2mu*mH2*mstau22) - mH2*(mstau22 + mu2 + 2*invdmstau1mu*
     mstau22*mu2)*power2(invdmstau1mu))*(-1 + power2(sa)))/(8.*mH2*mstau22) +
     DeltaInv(mstau22,mH2,mu2)*(-(Al*ca*mu*mu2*sa*power2(invdmstau1mu)*(-3*mH2*
     mstau22 - 2*mH2*mu2 - mstau22*mu2 + power2(mH2) + power2(mu2)))/(4.*
     mstau22) - (mu2*power2(Al)*power2(invdmstau1mu)*(-3*mH2*mstau22 - 2*mH2*
     mu2 - mstau22*mu2 + power2(mH2) + power2(mu2))*(-1 + power2(sa)))/(8.*
     mstau22) + (mu2*(-3*mH2*mstau22 - 2*mH2*mu2 - mstau22*mu2 + power2(mH2) +
     power2(mu2))*(1 + (-1 + mH2*mu2*power2(invdmstau1mu))*power2(sa)))/(8.*mH2
     *mstau22)) + (8*ca*invdmstau*mH2*(1 + invdmstau2mu*mH2)*mstau22*mu2*sa*sb
     + cb*(-3*mstau22*(mstau22 - mu2)*(-1 + power2(sa)) + mH2*(13*mstau22 + 2*
     mu2)*(-1 + power2(sa)) + power2(mH2)*(8*invdmstau2mu*mstau22*(-1 + power2(
     sa)) - 2*mu2*(mstau22 + mu2 + 2*invdmstau1mu*mstau22*mu2)*power2(
     invdmstau1mu)*power2(sa))))/(16.*cb*mstau22*power2(mH2))) + Fin3(mH2,
     mstau12,mu2,Q2)*((Al*mu*(ca*cb*sa*(2*invdmstau*(mstau12 + invdmstau1mu*mH2
     *mstau12) + mH2*(mstau12 + mu2 + 2*invdmstau2mu*mstau12*mu2)*power2(
     invdmstau2mu)) - 2*invdmstau*(1 + invdmstau1mu*mH2)*mstau12*sb*(-1 +
     power2(sa))))/(4.*cb*mH2*mstau12) + (power2(Al)*(4*invdmstau*(mstau12 +
     invdmstau1mu*mH2*mstau12) + mH2*(mstau12 + mu2 + 2*invdmstau2mu*mstau12*
     mu2)*power2(invdmstau2mu))*(-1 + power2(sa)))/(8.*mH2*mstau12) + DeltaInv(
     mH2,mstau12,mu2)*(-(Al*ca*mu*mu2*sa*power2(invdmstau2mu)*(-3*mH2*mstau12 -
     2*mH2*mu2 - mstau12*mu2 + power2(mH2) + power2(mu2)))/(4.*mstau12) - (mu2*
     power2(Al)*power2(invdmstau2mu)*(-3*mH2*mstau12 - 2*mH2*mu2 - mstau12*mu2
     + power2(mH2) + power2(mu2))*(-1 + power2(sa)))/(8.*mstau12) + (mu2*(-3*
     mH2*mstau12 - 2*mH2*mu2 - mstau12*mu2 + power2(mH2) + power2(mu2))*(1 + (-
     1 + mH2*mu2*power2(invdmstau2mu))*power2(sa)))/(8.*mH2*mstau12)) + (-8*ca*
     invdmstau*mH2*(1 + invdmstau1mu*mH2)*mstau12*mu2*sa*sb + cb*(-3*mstau12*(
     mstau12 - mu2)*(-1 + power2(sa)) + mH2*(13*mstau12 + 2*mu2)*(-1 + power2(
     sa)) + power2(mH2)*(8*invdmstau1mu*mstau12*(-1 + power2(sa)) - 2*mu2*(
     mstau12 + mu2 + 2*invdmstau2mu*mstau12*mu2)*power2(invdmstau2mu)*power2(sa
     ))))/(16.*cb*mstau12*power2(mH2))) + Fin3(mstau22,mz2,mu2,Q2)*((DeltaInv(
     mstau22,mz2,mu2)*power2(invdmstau1mu)*power2(mu2)*(-(mstau22*(mu2 + 3*mz2)
     ) + power2(mu2 - mz2)))/(8.*mstau22) - (invdmstau1mu*mu2*(invdmstau1mu*(
     mstau22 + mu2)*mz2 + 2*mstau22*mu2*mz2*power2(invdmstau1mu) + 4*invdmstau*
     mstau22*(mstau22 - mu2 - mz2)*power2(sb)))/(8.*mstau22*mz2)) + Fin3(
     mstau12,mz2,mu2,Q2)*((DeltaInv(mstau12,mz2,mu2)*power2(invdmstau2mu)*
     power2(mu2)*(-(mstau12*(mu2 + 3*mz2)) + power2(mu2 - mz2)))/(8.*mstau12) -
     (invdmstau2mu*mu2*(invdmstau2mu*(mstau12 + mu2)*mz2 + 2*mstau12*mu2*mz2*
     power2(invdmstau2mu) + 4*invdmstau*mstau12*(-mstau12 + mu2 + mz2)*power2(
     sb)))/(8.*mstau12*mz2)) + Fin3(mA2,mstau22,mu2,Q2)*(-(Al*invdmstau*
     invdmstau1mu*mu*(mA2 - mstau22 + mu2)*sb)/(2.*cb*mA2) - (invdmstau1mu*(-4*
     invdmstau*mstau22*(mA2 - mstau22 + mu2) + invdmstau1mu*mA2*(mstau22 + mu2
     + 2*invdmstau1mu*mstau22*mu2))*power2(Al))/(8.*mA2*mstau22) + DeltaInv(mA2
     ,mstau22,mu2)*((mu2*(-3*mA2*mstau22 - 2*mA2*mu2 - mstau22*mu2 + power2(mA2
     ) + power2(mu2)))/(8.*mA2*mstau22) + (mu2*power2(Al)*power2(invdmstau1mu)*
     (-3*mA2*mstau22 - 2*mA2*mu2 - mstau22*mu2 + power2(mA2) + power2(mu2)))/(
     8.*mstau22)) + (5*mstau22*(-mstau22 + mu2) - 8*invdmstau*invdmstau1mu*
     mstau22*mu2*power2(mA2)*power2(sb) + mA2*(-2*mu2 + 8*invdmstau*
     invdmstau1mu*mu2*power2(mstau22)*power2(sb) + mstau22*(3 - 8*invdmstau*
     invdmstau1mu*power2(mu2)*power2(sb))))/(16.*mstau22*power2(mA2))) + Fin3(
     mA2,mstau12,mu2,Q2)*((Al*invdmstau*invdmstau2mu*mu*(mA2 - mstau12 + mu2)*
     sb)/(2.*cb*mA2) - (invdmstau2mu*(4*invdmstau*mstau12*(mA2 - mstau12 + mu2)
     + invdmstau2mu*mA2*(mstau12 + mu2 + 2*invdmstau2mu*mstau12*mu2))*power2(Al
     ))/(8.*mA2*mstau12) + DeltaInv(mA2,mstau12,mu2)*((mu2*(-3*mA2*mstau12 - 2*
     mA2*mu2 - mstau12*mu2 + power2(mA2) + power2(mu2)))/(8.*mA2*mstau12) + (
     mu2*power2(Al)*power2(invdmstau2mu)*(-3*mA2*mstau12 - 2*mA2*mu2 - mstau12*
     mu2 + power2(mA2) + power2(mu2)))/(8.*mstau12)) + (5*mstau12*(-mstau12 +
     mu2) + 8*invdmstau*invdmstau2mu*mstau12*mu2*power2(mA2)*power2(sb) + mA2*(
     -2*mu2 - 8*invdmstau*invdmstau2mu*mu2*power2(mstau12)*power2(sb) + mstau12
     *(3 + 8*invdmstau*invdmstau2mu*power2(mu2)*power2(sb))))/(16.*mstau12*
     power2(mA2))) + Fin20(mstau12,mu2,Q2)*((invdmstau*power2(Al)*(invdmstau2mu
     *mh2*mH2*(-mstau12 + mu2) + mA2*(mh2 - mh2*power2(sa) + mH2*power2(sa))))/
     (2.*mA2*mh2*mH2) - (Al*invdmstau*mu*(ca*cb*mA2*(mh2 - mH2)*sa +
     invdmstau2mu*mh2*mH2*(-mstau12 + mu2)*sb + mA2*sb*(mh2 - mh2*power2(sa) +
     mH2*power2(sa))))/(2.*cb*mA2*mh2*mH2) + (8*ca*invdmstau*mh2*(mh2 - mH2)*
     mH2*mu2*mz2*sa*sb*power2(mA2) + cb*(5*(mstau12 - mu2)*mz2*power2(mh2)*
     power2(mH2) + 8*mA2*mz2*power2(mh2)*power2(mH2)*(-1 + invdmstau1mu*mu2 +
     invdmstau*invdmstau2mu*(mstau12 - mu2)*mu2*power2(sb)) + power2(mA2)*(3*(-
     mstau12 + mu2)*mz2*power2(mH2)*power2(sa) - 8*mh2*(-2 + invdmstau1mu*mu2)*
     mz2*power2(mH2)*power2(sa) + power2(mh2)*(3*(mstau12 - mu2)*mz2*(-1 +
     power2(sa)) + 8*mH2*(-2 + invdmstau1mu*mu2)*mz2*(-1 + power2(sa)) + 2*
     power2(mH2)*(invdmstau1mu*(mz2 - 8*invdmstau2mu*mu2*mz2) + invdmstau2mu*
     mz2*(4 + invdmstau2mu*(mstau12 - 3*mu2) + 2*mu2*(-mstau12 + mu2)*power2(
     invdmstau2mu)) + 4*invdmstau*invdmstau2mu*mu2*(-mstau12 + mu2)*power2(sb))
     ))))/(16.*cb*mz2*power2(mA2)*power2(mh2)*power2(mH2))) + Fin3(mstau22,mh2,
     mu2,Q2)*(-(Al*mu*sa*(2*invdmstau*(1 + invdmstau2mu*mh2)*mstau22*sa*sb + ca
     *cb*(-2*invdmstau*(mstau22 + invdmstau2mu*mh2*mstau22) + mh2*(mstau22 +
     mu2 + 2*invdmstau1mu*mstau22*mu2)*power2(invdmstau1mu))))/(4.*cb*mh2*
     mstau22) + (power2(Al)*(4*invdmstau*(mstau22 + invdmstau2mu*mh2*mstau22) -
     mh2*(mstau22 + mu2 + 2*invdmstau1mu*mstau22*mu2)*power2(invdmstau1mu))*
     power2(sa))/(8.*mh2*mstau22) + DeltaInv(mstau22,mh2,mu2)*((Al*ca*mu*mu2*sa
     *power2(invdmstau1mu)*(-3*mh2*mstau22 - 2*mh2*mu2 - mstau22*mu2 + power2(
     mh2) + power2(mu2)))/(4.*mstau22) - (mu2*(-3*mh2*mstau22 - 2*mh2*mu2 -
     mstau22*mu2 + power2(mh2) + power2(mu2))*(mh2*mu2*power2(invdmstau1mu)*(-1
      + power2(sa)) - power2(sa)))/(8.*mh2*mstau22) + (mu2*power2(Al)*power2(
     invdmstau1mu)*(-3*mh2*mstau22 - 2*mh2*mu2 - mstau22*mu2 + power2(mh2) +
     power2(mu2))*power2(sa))/(8.*mstau22)) + (-8*ca*invdmstau*mh2*(1 +
     invdmstau2mu*mh2)*mstau22*mu2*sa*sb + cb*(2*mu2*(mstau22 + mu2)*power2(
     invdmstau1mu)*power2(mh2)*(-1 + power2(sa)) - (13*mh2*mstau22 + 2*mh2*mu2
     + 3*mstau22*mu2 + 8*invdmstau2mu*mstau22*power2(mh2) - 3*power2(mstau22))*
     power2(sa) + 4*mstau22*power2(mh2)*power2(mu2)*(-1 + power2(sa))*power3(
     invdmstau1mu)))/(16.*cb*mstau22*power2(mh2))) + Fin3(mstau22,msntau2,mw2,
     Q2)*((mu2*DeltaInv(mstau22,msntau2,mw2)*(invdmstau2mu*(1 + invdmstau2mu*
     mu2)*(mstau24 + mw2*(-2*mstau22 + mw2)) + msntau2*(-3 + invdmsntau2mu*(-3*
     mstau22 + 2*mu2 - 3*mw2) - invdmstau2mu*(2*mstau22 + 3*mw2) + (-mstau24 +
     (mu2 - mw2)*mw2 + mstau22*(mu2 + 2*mw2))*power2(invdmsntau2mu) + mu2*(-2*
     mstau22 + mu2 - 3*mw2)*power2(invdmstau2mu))))/(4.*msntau2) + (mu2*(-
     invdmstau2mu + 2*msntau2*(1 + invdmsntau2mu*mu2)*power2(invdmsntau2mu) + (
     msntau2 - mu2)*power2(invdmstau2mu) + 2*msntau2*mu2*power3(invdmstau2mu)))
     /(4.*msntau2)) + Fin3(mstau22,mstau12,mz2,Q2)*((mu2*DeltaInv(mstau22,
     mstau12,mz2)*(invdmstau2mu*(1 + invdmstau2mu*mu2)*(mstau24 + mz2*(-2*
     mstau22 + mz2)) + mstau12*(-3 + invdmstau1mu*(-3*mstau22 + 2*mu2 - 3*mz2)
     - invdmstau2mu*(2*mstau22 + 3*mz2) + (-mstau24 + (mu2 - mz2)*mz2 + mstau22
     *(mu2 + 2*mz2))*power2(invdmstau1mu) + mu2*(-2*mstau22 + mu2 - 3*mz2)*
     power2(invdmstau2mu))))/(8.*mstau12) + (mu2*((mstau12 - mu2)*mz2*power2(
     invdmstau2mu) + 2*invdmstau1mu*mstau12*(invdmstau1mu*mz2 + mu2*mz2*power2(
     invdmstau1mu) + 2*invdmstau*(mstau22 - mu2 - mz2)*power2(sb)) +
     invdmstau2mu*(4*invdmstau*mstau12*(-mstau12 + mu2)*power2(sb) + mz2*(-1 +
     4*invdmstau*mstau12*power2(sb))) + 2*mstau12*mu2*mz2*power3(invdmstau2mu))
     )/(8.*mstau12*mz2)) + Fin3(mstau12,mh2,mu2,Q2)*((Al*mu*sa*(2*invdmstau*(1
     + invdmstau1mu*mh2)*mstau12*sa*sb - ca*cb*(2*invdmstau*(mstau12 +
     invdmstau1mu*mh2*mstau12) + mh2*(mstau12 + mu2 + 2*invdmstau2mu*mstau12*
     mu2)*power2(invdmstau2mu))))/(4.*cb*mh2*mstau12) - (power2(Al)*(4*
     invdmstau*(mstau12 + invdmstau1mu*mh2*mstau12) + mh2*(mstau12 + mu2 + 2*
     invdmstau2mu*mstau12*mu2)*power2(invdmstau2mu))*power2(sa))/(8.*mh2*
     mstau12) + DeltaInv(mstau12,mh2,mu2)*((Al*ca*mu*mu2*sa*power2(invdmstau2mu
     )*(-3*mh2*mstau12 - 2*mh2*mu2 - mstau12*mu2 + power2(mh2) + power2(mu2)))/
     (4.*mstau12) - (mu2*(-3*mh2*mstau12 - 2*mh2*mu2 - mstau12*mu2 + power2(mh2
     ) + power2(mu2))*(mh2*mu2*power2(invdmstau2mu)*(-1 + power2(sa)) - power2(
     sa)))/(8.*mh2*mstau12) + (mu2*power2(Al)*power2(invdmstau2mu)*(-3*mh2*
     mstau12 - 2*mh2*mu2 - mstau12*mu2 + power2(mh2) + power2(mu2))*power2(sa))
     /(8.*mstau12)) + (8*ca*invdmstau*mh2*(1 + invdmstau1mu*mh2)*mstau12*mu2*sa
     *sb + cb*(2*mu2*(mstau12 + mu2)*power2(invdmstau2mu)*power2(mh2)*(-1 +
     power2(sa)) - (13*mh2*mstau12 + 2*mh2*mu2 + 3*mstau12*mu2 + 8*invdmstau1mu
     *mstau12*power2(mh2) - 3*power2(mstau12))*power2(sa) + 4*mstau12*power2(
     mh2)*power2(mu2)*(-1 + power2(sa))*power3(invdmstau2mu)))/(16.*cb*mstau12*
     power2(mh2))) + (mu2*DeltaInv(mA2,msntau2,mu2)*(-((msntau2 - mu2)*power2(
     mu2)*(42 + 6*logmsntau2Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*power2(
     logmsntau2Q2) + 3*power2(logmu2Q2) + power2(Pi))) - mA2*mu2*(-3*mu2*(4*
     logmA2Q2*logmsntau2Q2 - 2*(-6 + logmu2Q2)*logmu2Q2 - 2*logmA2Q2*(3 +
     logmu2Q2) - 2*logmsntau2Q2*(3 + logmu2Q2) + power2(logmA2Q2) + power2(
     logmsntau2Q2)) + msntau2*(294 + 36*logmsntau2Q2*(-5 + logmu2Q2) - 54*
     logmu2Q2 + 6*logmA2Q2*(4*logmsntau2Q2 - 3*(1 + logmu2Q2)) + 3*power2(
     logmA2Q2) + 30*power2(logmsntau2Q2) + 9*power2(logmu2Q2) + 7*power2(Pi)))
     - power2(mA2)*(3*mu2*(42 + 4*logmA2Q2*(-3 + 2*logmsntau2Q2 - logmu2Q2) + 2
     *logmsntau2Q2*(-15 + logmu2Q2) + 6*logmu2Q2 + 2*power2(logmA2Q2) + 5*
     power2(logmsntau2Q2) - power2(logmu2Q2) + power2(Pi)) + msntau2*(6*
     logmA2Q2*(-9 + 4*logmsntau2Q2 - logmu2Q2) + 6*logmsntau2Q2*(-15 + logmu2Q2
     ) + 9*power2(logmA2Q2) + 15*power2(logmsntau2Q2) + 4*(42 + power2(Pi)))) +
     (6*logmA2Q2*(-3 + 2*logmsntau2Q2 - logmu2Q2) + 6*logmsntau2Q2*(-9 +
     logmu2Q2) + 3*power2(logmA2Q2) + 9*power2(logmsntau2Q2) + 2*(42 + power2(
     Pi)))*power3(mA2)))/(48.*mA2*msntau2) + Fin3(mA2,mstau22,mstau12,Q2)*((Al*
     invdmstau*mu*(-(invdmstau2mu*(mA2 - mstau12 + mu2)) + invdmstau1mu*(mA2 -
     mstau22 + mu2))*sb)/(2.*cb*mA2) + (invdmstau*mu2*(-(invdmstau2mu*(mA2 -
     mstau12 + mu2)) + invdmstau1mu*(mA2 - mstau22 + mu2))*power2(sb))/(2.*mA2)
     + (power2(Al)*(-2*mstau12 + mA2*(-2 - 4*invdmstau*mstau12*(invdmstau2mu*
     mstau12 - invdmstau1mu*mstau22 + invdmstau1mu*mu2 - invdmstau2mu*mu2)) +
     power2(mA2)*(invdmstau2mu*(-1 + 4*invdmstau*mstau12) + 2*invdmstau1mu*
     mstau12*(-2*invdmstau + invdmstau1mu + mu2*power2(invdmstau1mu)) + (
     mstau12 - mu2)*power2(invdmstau2mu) + 2*mstau12*mu2*power3(invdmstau2mu)))
     )/(8.*mstau12*power2(mA2)) + (DeltaInv(mA2,mstau22,mstau12)*power2(Al)*(-6
     *mstau12*mstau22 + 2*mstau24 + mA2*(-4*mstau22 + invdmstau2mu*mstau24 +
     mstau24*mu2*power2(invdmstau2mu) + mstau12*(-5 - 2*invdmstau2mu*mstau22 +
     invdmstau1mu*(-3*mstau22 + 2*mu2) + (-mstau24 + mstau22*mu2)*power2(
     invdmstau1mu) + mu2*(-2*mstau22 + mu2)*power2(invdmstau2mu))) + (2 - 3*
     invdmstau1mu*mstau12 - invdmstau2mu*(3*mstau12 + 2*mstau22) + mstau12*(2*
     mstau22 + mu2)*power2(invdmstau1mu) - (3*mstau12 + 2*mstau22)*mu2*power2(
     invdmstau2mu))*power2(mA2) + (invdmstau2mu - mstau12*power2(invdmstau1mu)
     + mu2*power2(invdmstau2mu))*power3(mA2)))/(8.*mA2*mstau12)) + DeltaInv(mA2
     ,mstau12,mu2)*(-(mu2*power2(Al)*power2(invdmstau2mu)*((mstau12 - mu2)*
     power2(mu2)*(42 + 6*logmstau12Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*power2(
     logmstau12Q2) + 3*power2(logmu2Q2) + power2(Pi)) + mA2*mu2*(-3*mu2*(4*
     logmA2Q2*logmstau12Q2 - 2*(-6 + logmu2Q2)*logmu2Q2 - 2*logmA2Q2*(3 +
     logmu2Q2) - 2*logmstau12Q2*(3 + logmu2Q2) + power2(logmA2Q2) + power2(
     logmstau12Q2)) + mstau12*(294 + 36*logmstau12Q2*(-5 + logmu2Q2) - 54*
     logmu2Q2 + 6*logmA2Q2*(4*logmstau12Q2 - 3*(1 + logmu2Q2)) + 3*power2(
     logmA2Q2) + 30*power2(logmstau12Q2) + 9*power2(logmu2Q2) + 7*power2(Pi)))
     + power2(mA2)*(3*mu2*(42 + 4*logmA2Q2*(-3 + 2*logmstau12Q2 - logmu2Q2) + 2
     *logmstau12Q2*(-15 + logmu2Q2) + 6*logmu2Q2 + 2*power2(logmA2Q2) + 5*
     power2(logmstau12Q2) - power2(logmu2Q2) + power2(Pi)) + mstau12*(6*
     logmA2Q2*(-9 + 4*logmstau12Q2 - logmu2Q2) + 6*logmstau12Q2*(-15 + logmu2Q2
     ) + 9*power2(logmA2Q2) + 15*power2(logmstau12Q2) + 4*(42 + power2(Pi)))) -
     (6*logmA2Q2*(-3 + 2*logmstau12Q2 - logmu2Q2) + 6*logmstau12Q2*(-9 +
     logmu2Q2) + 3*power2(logmA2Q2) + 9*power2(logmstau12Q2) + 2*(42 + power2(
     Pi)))*power3(mA2)))/(48.*mstau12) + (mu2*(-((mstau12 - mu2)*power2(mu2)*(
     42 + 6*logmstau12Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*power2(logmstau12Q2)
     + 3*power2(logmu2Q2) + power2(Pi))) - mA2*mu2*(-3*mu2*(4*logmA2Q2*
     logmstau12Q2 - 2*(-6 + logmu2Q2)*logmu2Q2 - 2*logmA2Q2*(3 + logmu2Q2) - 2*
     logmstau12Q2*(3 + logmu2Q2) + power2(logmA2Q2) + power2(logmstau12Q2)) +
     mstau12*(294 + 36*logmstau12Q2*(-5 + logmu2Q2) - 54*logmu2Q2 + 6*logmA2Q2*
     (4*logmstau12Q2 - 3*(1 + logmu2Q2)) + 3*power2(logmA2Q2) + 30*power2(
     logmstau12Q2) + 9*power2(logmu2Q2) + 7*power2(Pi))) - power2(mA2)*(3*mu2*(
     42 + 4*logmA2Q2*(-3 + 2*logmstau12Q2 - logmu2Q2) + 2*logmstau12Q2*(-15 +
     logmu2Q2) + 6*logmu2Q2 + 2*power2(logmA2Q2) + 5*power2(logmstau12Q2) -
     power2(logmu2Q2) + power2(Pi)) + mstau12*(6*logmA2Q2*(-9 + 4*logmstau12Q2
     - logmu2Q2) + 6*logmstau12Q2*(-15 + logmu2Q2) + 9*power2(logmA2Q2) + 15*
     power2(logmstau12Q2) + 4*(42 + power2(Pi)))) + (6*logmA2Q2*(-3 + 2*
     logmstau12Q2 - logmu2Q2) + 6*logmstau12Q2*(-9 + logmu2Q2) + 3*power2(
     logmA2Q2) + 9*power2(logmstau12Q2) + 2*(42 + power2(Pi)))*power3(mA2)))/(
     48.*mA2*mstau12)) + DeltaInv(mA2,mstau22,mu2)*(-(mu2*power2(Al)*power2(
     invdmstau1mu)*((mstau22 - mu2)*power2(mu2)*(42 + 6*logmstau22Q2*(-3 +
     logmu2Q2) - 18*logmu2Q2 + 3*power2(logmstau22Q2) + 3*power2(logmu2Q2) +
     power2(Pi)) + mA2*mu2*(-3*mu2*(4*logmA2Q2*logmstau22Q2 - 2*(-6 + logmu2Q2)
     *logmu2Q2 - 2*logmA2Q2*(3 + logmu2Q2) - 2*logmstau22Q2*(3 + logmu2Q2) +
     power2(logmA2Q2) + power2(logmstau22Q2)) + mstau22*(294 + 36*logmstau22Q2*
     (-5 + logmu2Q2) - 54*logmu2Q2 + 6*logmA2Q2*(4*logmstau22Q2 - 3*(1 +
     logmu2Q2)) + 3*power2(logmA2Q2) + 30*power2(logmstau22Q2) + 9*power2(
     logmu2Q2) + 7*power2(Pi))) + power2(mA2)*(3*mu2*(42 + 4*logmA2Q2*(-3 + 2*
     logmstau22Q2 - logmu2Q2) + 2*logmstau22Q2*(-15 + logmu2Q2) + 6*logmu2Q2 +
     2*power2(logmA2Q2) + 5*power2(logmstau22Q2) - power2(logmu2Q2) + power2(Pi
     )) + mstau22*(6*logmA2Q2*(-9 + 4*logmstau22Q2 - logmu2Q2) + 6*logmstau22Q2
     *(-15 + logmu2Q2) + 9*power2(logmA2Q2) + 15*power2(logmstau22Q2) + 4*(42 +
     power2(Pi)))) - (6*logmA2Q2*(-3 + 2*logmstau22Q2 - logmu2Q2) + 6*
     logmstau22Q2*(-9 + logmu2Q2) + 3*power2(logmA2Q2) + 9*power2(logmstau22Q2)
     + 2*(42 + power2(Pi)))*power3(mA2)))/(48.*mstau22) + (mu2*(-((mstau22 -
     mu2)*power2(mu2)*(42 + 6*logmstau22Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*
     power2(logmstau22Q2) + 3*power2(logmu2Q2) + power2(Pi))) - mA2*mu2*(-3*mu2
     *(4*logmA2Q2*logmstau22Q2 - 2*(-6 + logmu2Q2)*logmu2Q2 - 2*logmA2Q2*(3 +
     logmu2Q2) - 2*logmstau22Q2*(3 + logmu2Q2) + power2(logmA2Q2) + power2(
     logmstau22Q2)) + mstau22*(294 + 36*logmstau22Q2*(-5 + logmu2Q2) - 54*
     logmu2Q2 + 6*logmA2Q2*(4*logmstau22Q2 - 3*(1 + logmu2Q2)) + 3*power2(
     logmA2Q2) + 30*power2(logmstau22Q2) + 9*power2(logmu2Q2) + 7*power2(Pi)))
     - power2(mA2)*(3*mu2*(42 + 4*logmA2Q2*(-3 + 2*logmstau22Q2 - logmu2Q2) + 2
     *logmstau22Q2*(-15 + logmu2Q2) + 6*logmu2Q2 + 2*power2(logmA2Q2) + 5*
     power2(logmstau22Q2) - power2(logmu2Q2) + power2(Pi)) + mstau22*(6*
     logmA2Q2*(-9 + 4*logmstau22Q2 - logmu2Q2) + 6*logmstau22Q2*(-15 + logmu2Q2
     ) + 9*power2(logmA2Q2) + 15*power2(logmstau22Q2) + 4*(42 + power2(Pi)))) +
     (6*logmA2Q2*(-3 + 2*logmstau22Q2 - logmu2Q2) + 6*logmstau22Q2*(-9 +
     logmu2Q2) + 3*power2(logmA2Q2) + 9*power2(logmstau22Q2) + 2*(42 + power2(
     Pi)))*power3(mA2)))/(48.*mA2*mstau22)) - (Fin20(mA2,mC2,Q2)*(-3*mA2 + 4*
     mC2 + mC2*power2(invdmAC)*power2(mA2) - 2*mA2*power2(invdmAC)*power2(mC2)
     + power2(invdmAC)*power3(mC2)))/(16.*power2(mA2)) - (Fin20(mC2,mh2,Q2)*
     power2(sa)*(5*mh2 - 2*mh2*power2(invdmCh)*power2(mC2) + mC2*(-4 + power2(
     invdmCh)*power2(mh2)) + power2(invdmCh)*power3(mC2)))/(16.*power2(mh2)) +
     (Fin20(mC2,mH2,Q2)*(-1 + power2(sa))*(5*mH2 - 2*mH2*power2(invdmCH)*power2
     (mC2) + mC2*(-4 + power2(invdmCH)*power2(mH2)) + power2(invdmCH)*power3(
     mC2)))/(16.*power2(mH2)) - (mu2*DeltaInv(mC2,msntau2,mu2)*power2(Al)*
     power2(invdmstau2mu)*((msntau2 - mu2)*power2(mu2)*(42 + 6*logmsntau2Q2*(-3
      + logmu2Q2) - 18*logmu2Q2 + 3*power2(logmsntau2Q2) + 3*power2(logmu2Q2) +
     power2(Pi)) + mC2*mu2*(-3*mu2*(4*logmC2Q2*logmsntau2Q2 - 2*(-6 + logmu2Q2)
     *logmu2Q2 - 2*logmC2Q2*(3 + logmu2Q2) - 2*logmsntau2Q2*(3 + logmu2Q2) +
     power2(logmC2Q2) + power2(logmsntau2Q2)) + msntau2*(294 + 36*logmsntau2Q2*
     (-5 + logmu2Q2) - 54*logmu2Q2 + 6*logmC2Q2*(4*logmsntau2Q2 - 3*(1 +
     logmu2Q2)) + 3*power2(logmC2Q2) + 30*power2(logmsntau2Q2) + 9*power2(
     logmu2Q2) + 7*power2(Pi))) + power2(mC2)*(3*mu2*(42 + 4*logmC2Q2*(-3 + 2*
     logmsntau2Q2 - logmu2Q2) + 2*logmsntau2Q2*(-15 + logmu2Q2) + 6*logmu2Q2 +
     2*power2(logmC2Q2) + 5*power2(logmsntau2Q2) - power2(logmu2Q2) + power2(Pi
     )) + msntau2*(6*logmC2Q2*(-9 + 4*logmsntau2Q2 - logmu2Q2) + 6*logmsntau2Q2
     *(-15 + logmu2Q2) + 9*power2(logmC2Q2) + 15*power2(logmsntau2Q2) + 4*(42 +
     power2(Pi)))) - (6*logmC2Q2*(-3 + 2*logmsntau2Q2 - logmu2Q2) + 6*
     logmsntau2Q2*(-9 + logmu2Q2) + 3*power2(logmC2Q2) + 9*power2(logmsntau2Q2)
     + 2*(42 + power2(Pi)))*power3(mC2)))/(24.*msntau2) + Fin3(mstau22,mC2,
     msntau2,Q2)*((power2(Al)*(mC2 + msntau2 + power2(mC2)*(-invdmstau2mu + 2*
     msntau2*(1 + invdmsntau2mu*mu2)*power2(invdmsntau2mu) + (msntau2 - mu2)*
     power2(invdmstau2mu) + 2*msntau2*mu2*power3(invdmstau2mu))))/(4.*msntau2*
     power2(mC2)) + (DeltaInv(mstau22,mC2,msntau2)*power2(Al)*(3*msntau2*
     mstau22 - mstau24 + mC2*(2*mstau22 + invdmstau2mu*mstau24 + mstau24*mu2*
     power2(invdmstau2mu) + msntau2*(-2 - 2*invdmstau2mu*mstau22 +
     invdmsntau2mu*(-3*mstau22 + 2*mu2) + (-mstau24 + mstau22*mu2)*power2(
     invdmsntau2mu) + mu2*(-2*mstau22 + mu2)*power2(invdmstau2mu))) + (-1 - 3*
     invdmsntau2mu*msntau2 - invdmstau2mu*(3*msntau2 + 2*mstau22) + msntau2*(2*
     mstau22 + mu2)*power2(invdmsntau2mu) - (3*msntau2 + 2*mstau22)*mu2*power2(
     invdmstau2mu))*power2(mC2) + (invdmstau2mu - msntau2*power2(invdmsntau2mu)
     + mu2*power2(invdmstau2mu))*power3(mC2)))/(4.*mC2*msntau2)) + DeltaInv(
     mstau22,mC2,mu2)*(-(mu2*power2(Al)*power2(invdmsntau2mu)*((mstau22 - mu2)*
     power2(mu2)*(42 + 6*logmstau22Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*power2(
     logmstau22Q2) + 3*power2(logmu2Q2) + power2(Pi)) + mC2*mu2*(-3*mu2*(4*
     logmC2Q2*logmstau22Q2 - 2*(-6 + logmu2Q2)*logmu2Q2 - 2*logmC2Q2*(3 +
     logmu2Q2) - 2*logmstau22Q2*(3 + logmu2Q2) + power2(logmC2Q2) + power2(
     logmstau22Q2)) + mstau22*(294 + 36*logmstau22Q2*(-5 + logmu2Q2) - 54*
     logmu2Q2 + 6*logmC2Q2*(4*logmstau22Q2 - 3*(1 + logmu2Q2)) + 3*power2(
     logmC2Q2) + 30*power2(logmstau22Q2) + 9*power2(logmu2Q2) + 7*power2(Pi)))
     + power2(mC2)*(3*mu2*(42 + 4*logmC2Q2*(-3 + 2*logmstau22Q2 - logmu2Q2) + 2
     *logmstau22Q2*(-15 + logmu2Q2) + 6*logmu2Q2 + 2*power2(logmC2Q2) + 5*
     power2(logmstau22Q2) - power2(logmu2Q2) + power2(Pi)) + mstau22*(6*
     logmC2Q2*(-9 + 4*logmstau22Q2 - logmu2Q2) + 6*logmstau22Q2*(-15 + logmu2Q2
     ) + 9*power2(logmC2Q2) + 15*power2(logmstau22Q2) + 4*(42 + power2(Pi)))) -
     (6*logmC2Q2*(-3 + 2*logmstau22Q2 - logmu2Q2) + 6*logmstau22Q2*(-9 +
     logmu2Q2) + 3*power2(logmC2Q2) + 9*power2(logmstau22Q2) + 2*(42 + power2(
     Pi)))*power3(mC2)))/(24.*mstau22) + (mu2*(-((mstau22 - mu2)*power2(mu2)*(
     42 + 6*logmstau22Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*power2(logmstau22Q2)
     + 3*power2(logmu2Q2) + power2(Pi))) - mC2*mu2*(-3*mu2*(4*logmC2Q2*
     logmstau22Q2 - 2*(-6 + logmu2Q2)*logmu2Q2 - 2*logmC2Q2*(3 + logmu2Q2) - 2*
     logmstau22Q2*(3 + logmu2Q2) + power2(logmC2Q2) + power2(logmstau22Q2)) +
     mstau22*(294 + 36*logmstau22Q2*(-5 + logmu2Q2) - 54*logmu2Q2 + 6*logmC2Q2*
     (4*logmstau22Q2 - 3*(1 + logmu2Q2)) + 3*power2(logmC2Q2) + 30*power2(
     logmstau22Q2) + 9*power2(logmu2Q2) + 7*power2(Pi))) - power2(mC2)*(3*mu2*(
     42 + 4*logmC2Q2*(-3 + 2*logmstau22Q2 - logmu2Q2) + 2*logmstau22Q2*(-15 +
     logmu2Q2) + 6*logmu2Q2 + 2*power2(logmC2Q2) + 5*power2(logmstau22Q2) -
     power2(logmu2Q2) + power2(Pi)) + mstau22*(6*logmC2Q2*(-9 + 4*logmstau22Q2
     - logmu2Q2) + 6*logmstau22Q2*(-15 + logmu2Q2) + 9*power2(logmC2Q2) + 15*
     power2(logmstau22Q2) + 4*(42 + power2(Pi)))) + (6*logmC2Q2*(-3 + 2*
     logmstau22Q2 - logmu2Q2) + 6*logmstau22Q2*(-9 + logmu2Q2) + 3*power2(
     logmC2Q2) + 9*power2(logmstau22Q2) + 2*(42 + power2(Pi)))*power3(mC2)))/(
     24.*mC2*mstau22)) + (mu2*DeltaInv(msntau2,mh2,mu2)*power2(sa)*(-((msntau2
     - mu2)*power2(mu2)*(42 + 6*logmsntau2Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*
     power2(logmsntau2Q2) + 3*power2(logmu2Q2) + power2(Pi))) - mh2*mu2*(-3*mu2
     *(4*logmh2Q2*logmsntau2Q2 - 2*(-6 + logmu2Q2)*logmu2Q2 - 2*logmh2Q2*(3 +
     logmu2Q2) - 2*logmsntau2Q2*(3 + logmu2Q2) + power2(logmh2Q2) + power2(
     logmsntau2Q2)) + msntau2*(294 + 36*logmsntau2Q2*(-5 + logmu2Q2) - 54*
     logmu2Q2 + 6*logmh2Q2*(4*logmsntau2Q2 - 3*(1 + logmu2Q2)) + 3*power2(
     logmh2Q2) + 30*power2(logmsntau2Q2) + 9*power2(logmu2Q2) + 7*power2(Pi)))
     - power2(mh2)*(3*mu2*(42 + 4*logmh2Q2*(-3 + 2*logmsntau2Q2 - logmu2Q2) + 2
     *logmsntau2Q2*(-15 + logmu2Q2) + 6*logmu2Q2 + 2*power2(logmh2Q2) + 5*
     power2(logmsntau2Q2) - power2(logmu2Q2) + power2(Pi)) + msntau2*(6*
     logmh2Q2*(-9 + 4*logmsntau2Q2 - logmu2Q2) + 6*logmsntau2Q2*(-15 + logmu2Q2
     ) + 9*power2(logmh2Q2) + 15*power2(logmsntau2Q2) + 4*(42 + power2(Pi)))) +
     (6*logmh2Q2*(-3 + 2*logmsntau2Q2 - logmu2Q2) + 6*logmsntau2Q2*(-9 +
     logmu2Q2) + 3*power2(logmh2Q2) + 9*power2(logmsntau2Q2) + 2*(42 + power2(
     Pi)))*power3(mh2)))/(48.*mh2*msntau2) + DeltaInv(mstau12,mh2,mu2)*(-(Al*ca
     *mu*mu2*sa*power2(invdmstau2mu)*((mstau12 - mu2)*power2(mu2)*(42 + 6*
     logmstau12Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*power2(logmstau12Q2) + 3*
     power2(logmu2Q2) + power2(Pi)) + mh2*mu2*(-3*mu2*(4*logmh2Q2*logmstau12Q2
     - 2*(-6 + logmu2Q2)*logmu2Q2 - 2*logmh2Q2*(3 + logmu2Q2) - 2*logmstau12Q2*
     (3 + logmu2Q2) + power2(logmh2Q2) + power2(logmstau12Q2)) + mstau12*(294 +
     36*logmstau12Q2*(-5 + logmu2Q2) - 54*logmu2Q2 + 6*logmh2Q2*(4*logmstau12Q2
      - 3*(1 + logmu2Q2)) + 3*power2(logmh2Q2) + 30*power2(logmstau12Q2) + 9*
     power2(logmu2Q2) + 7*power2(Pi))) + power2(mh2)*(3*mu2*(42 + 4*logmh2Q2*(-
     3 + 2*logmstau12Q2 - logmu2Q2) + 2*logmstau12Q2*(-15 + logmu2Q2) + 6*
     logmu2Q2 + 2*power2(logmh2Q2) + 5*power2(logmstau12Q2) - power2(logmu2Q2)
     + power2(Pi)) + mstau12*(6*logmh2Q2*(-9 + 4*logmstau12Q2 - logmu2Q2) + 6*
     logmstau12Q2*(-15 + logmu2Q2) + 9*power2(logmh2Q2) + 15*power2(
     logmstau12Q2) + 4*(42 + power2(Pi)))) - (6*logmh2Q2*(-3 + 2*logmstau12Q2 -
     logmu2Q2) + 6*logmstau12Q2*(-9 + logmu2Q2) + 3*power2(logmh2Q2) + 9*power2
     (logmstau12Q2) + 2*(42 + power2(Pi)))*power3(mh2)))/(24.*mstau12) - (mu2*
     power2(Al)*power2(invdmstau2mu)*power2(sa)*((mstau12 - mu2)*power2(mu2)*(
     42 + 6*logmstau12Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*power2(logmstau12Q2)
     + 3*power2(logmu2Q2) + power2(Pi)) + mh2*mu2*(-3*mu2*(4*logmh2Q2*
     logmstau12Q2 - 2*(-6 + logmu2Q2)*logmu2Q2 - 2*logmh2Q2*(3 + logmu2Q2) - 2*
     logmstau12Q2*(3 + logmu2Q2) + power2(logmh2Q2) + power2(logmstau12Q2)) +
     mstau12*(294 + 36*logmstau12Q2*(-5 + logmu2Q2) - 54*logmu2Q2 + 6*logmh2Q2*
     (4*logmstau12Q2 - 3*(1 + logmu2Q2)) + 3*power2(logmh2Q2) + 30*power2(
     logmstau12Q2) + 9*power2(logmu2Q2) + 7*power2(Pi))) + power2(mh2)*(3*mu2*(
     42 + 4*logmh2Q2*(-3 + 2*logmstau12Q2 - logmu2Q2) + 2*logmstau12Q2*(-15 +
     logmu2Q2) + 6*logmu2Q2 + 2*power2(logmh2Q2) + 5*power2(logmstau12Q2) -
     power2(logmu2Q2) + power2(Pi)) + mstau12*(6*logmh2Q2*(-9 + 4*logmstau12Q2
     - logmu2Q2) + 6*logmstau12Q2*(-15 + logmu2Q2) + 9*power2(logmh2Q2) + 15*
     power2(logmstau12Q2) + 4*(42 + power2(Pi)))) - (6*logmh2Q2*(-3 + 2*
     logmstau12Q2 - logmu2Q2) + 6*logmstau12Q2*(-9 + logmu2Q2) + 3*power2(
     logmh2Q2) + 9*power2(logmstau12Q2) + 2*(42 + power2(Pi)))*power3(mh2)))/(
     48.*mstau12) - (mu2*(mh2*mu2*power2(invdmstau2mu)*(-1 + power2(sa)) -
     power2(sa))*(-((mstau12 - mu2)*power2(mu2)*(42 + 6*logmstau12Q2*(-3 +
     logmu2Q2) - 18*logmu2Q2 + 3*power2(logmstau12Q2) + 3*power2(logmu2Q2) +
     power2(Pi))) - mh2*mu2*(-3*mu2*(4*logmh2Q2*logmstau12Q2 - 2*(-6 + logmu2Q2
     )*logmu2Q2 - 2*logmh2Q2*(3 + logmu2Q2) - 2*logmstau12Q2*(3 + logmu2Q2) +
     power2(logmh2Q2) + power2(logmstau12Q2)) + mstau12*(294 + 36*logmstau12Q2*
     (-5 + logmu2Q2) - 54*logmu2Q2 + 6*logmh2Q2*(4*logmstau12Q2 - 3*(1 +
     logmu2Q2)) + 3*power2(logmh2Q2) + 30*power2(logmstau12Q2) + 9*power2(
     logmu2Q2) + 7*power2(Pi))) - power2(mh2)*(3*mu2*(42 + 4*logmh2Q2*(-3 + 2*
     logmstau12Q2 - logmu2Q2) + 2*logmstau12Q2*(-15 + logmu2Q2) + 6*logmu2Q2 +
     2*power2(logmh2Q2) + 5*power2(logmstau12Q2) - power2(logmu2Q2) + power2(Pi
     )) + mstau12*(6*logmh2Q2*(-9 + 4*logmstau12Q2 - logmu2Q2) + 6*logmstau12Q2
     *(-15 + logmu2Q2) + 9*power2(logmh2Q2) + 15*power2(logmstau12Q2) + 4*(42 +
     power2(Pi)))) + (6*logmh2Q2*(-3 + 2*logmstau12Q2 - logmu2Q2) + 6*
     logmstau12Q2*(-9 + logmu2Q2) + 3*power2(logmh2Q2) + 9*power2(logmstau12Q2)
     + 2*(42 + power2(Pi)))*power3(mh2)))/(48.*mh2*mstau12)) + DeltaInv(mstau22
     ,mh2,mu2)*(-(Al*ca*mu*mu2*sa*power2(invdmstau1mu)*((mstau22 - mu2)*power2(
     mu2)*(42 + 6*logmstau22Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*power2(
     logmstau22Q2) + 3*power2(logmu2Q2) + power2(Pi)) + mh2*mu2*(-3*mu2*(4*
     logmh2Q2*logmstau22Q2 - 2*(-6 + logmu2Q2)*logmu2Q2 - 2*logmh2Q2*(3 +
     logmu2Q2) - 2*logmstau22Q2*(3 + logmu2Q2) + power2(logmh2Q2) + power2(
     logmstau22Q2)) + mstau22*(294 + 36*logmstau22Q2*(-5 + logmu2Q2) - 54*
     logmu2Q2 + 6*logmh2Q2*(4*logmstau22Q2 - 3*(1 + logmu2Q2)) + 3*power2(
     logmh2Q2) + 30*power2(logmstau22Q2) + 9*power2(logmu2Q2) + 7*power2(Pi)))
     + power2(mh2)*(3*mu2*(42 + 4*logmh2Q2*(-3 + 2*logmstau22Q2 - logmu2Q2) + 2
     *logmstau22Q2*(-15 + logmu2Q2) + 6*logmu2Q2 + 2*power2(logmh2Q2) + 5*
     power2(logmstau22Q2) - power2(logmu2Q2) + power2(Pi)) + mstau22*(6*
     logmh2Q2*(-9 + 4*logmstau22Q2 - logmu2Q2) + 6*logmstau22Q2*(-15 + logmu2Q2
     ) + 9*power2(logmh2Q2) + 15*power2(logmstau22Q2) + 4*(42 + power2(Pi)))) -
     (6*logmh2Q2*(-3 + 2*logmstau22Q2 - logmu2Q2) + 6*logmstau22Q2*(-9 +
     logmu2Q2) + 3*power2(logmh2Q2) + 9*power2(logmstau22Q2) + 2*(42 + power2(
     Pi)))*power3(mh2)))/(24.*mstau22) - (mu2*power2(Al)*power2(invdmstau1mu)*
     power2(sa)*((mstau22 - mu2)*power2(mu2)*(42 + 6*logmstau22Q2*(-3 +
     logmu2Q2) - 18*logmu2Q2 + 3*power2(logmstau22Q2) + 3*power2(logmu2Q2) +
     power2(Pi)) + mh2*mu2*(-3*mu2*(4*logmh2Q2*logmstau22Q2 - 2*(-6 + logmu2Q2)
     *logmu2Q2 - 2*logmh2Q2*(3 + logmu2Q2) - 2*logmstau22Q2*(3 + logmu2Q2) +
     power2(logmh2Q2) + power2(logmstau22Q2)) + mstau22*(294 + 36*logmstau22Q2*
     (-5 + logmu2Q2) - 54*logmu2Q2 + 6*logmh2Q2*(4*logmstau22Q2 - 3*(1 +
     logmu2Q2)) + 3*power2(logmh2Q2) + 30*power2(logmstau22Q2) + 9*power2(
     logmu2Q2) + 7*power2(Pi))) + power2(mh2)*(3*mu2*(42 + 4*logmh2Q2*(-3 + 2*
     logmstau22Q2 - logmu2Q2) + 2*logmstau22Q2*(-15 + logmu2Q2) + 6*logmu2Q2 +
     2*power2(logmh2Q2) + 5*power2(logmstau22Q2) - power2(logmu2Q2) + power2(Pi
     )) + mstau22*(6*logmh2Q2*(-9 + 4*logmstau22Q2 - logmu2Q2) + 6*logmstau22Q2
     *(-15 + logmu2Q2) + 9*power2(logmh2Q2) + 15*power2(logmstau22Q2) + 4*(42 +
     power2(Pi)))) - (6*logmh2Q2*(-3 + 2*logmstau22Q2 - logmu2Q2) + 6*
     logmstau22Q2*(-9 + logmu2Q2) + 3*power2(logmh2Q2) + 9*power2(logmstau22Q2)
     + 2*(42 + power2(Pi)))*power3(mh2)))/(48.*mstau22) - (mu2*(mh2*mu2*power2(
     invdmstau1mu)*(-1 + power2(sa)) - power2(sa))*(-((mstau22 - mu2)*power2(
     mu2)*(42 + 6*logmstau22Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*power2(
     logmstau22Q2) + 3*power2(logmu2Q2) + power2(Pi))) - mh2*mu2*(-3*mu2*(4*
     logmh2Q2*logmstau22Q2 - 2*(-6 + logmu2Q2)*logmu2Q2 - 2*logmh2Q2*(3 +
     logmu2Q2) - 2*logmstau22Q2*(3 + logmu2Q2) + power2(logmh2Q2) + power2(
     logmstau22Q2)) + mstau22*(294 + 36*logmstau22Q2*(-5 + logmu2Q2) - 54*
     logmu2Q2 + 6*logmh2Q2*(4*logmstau22Q2 - 3*(1 + logmu2Q2)) + 3*power2(
     logmh2Q2) + 30*power2(logmstau22Q2) + 9*power2(logmu2Q2) + 7*power2(Pi)))
     - power2(mh2)*(3*mu2*(42 + 4*logmh2Q2*(-3 + 2*logmstau22Q2 - logmu2Q2) + 2
     *logmstau22Q2*(-15 + logmu2Q2) + 6*logmu2Q2 + 2*power2(logmh2Q2) + 5*
     power2(logmstau22Q2) - power2(logmu2Q2) + power2(Pi)) + mstau22*(6*
     logmh2Q2*(-9 + 4*logmstau22Q2 - logmu2Q2) + 6*logmstau22Q2*(-15 + logmu2Q2
     ) + 9*power2(logmh2Q2) + 15*power2(logmstau22Q2) + 4*(42 + power2(Pi)))) +
     (6*logmh2Q2*(-3 + 2*logmstau22Q2 - logmu2Q2) + 6*logmstau22Q2*(-9 +
     logmu2Q2) + 3*power2(logmh2Q2) + 9*power2(logmstau22Q2) + 2*(42 + power2(
     Pi)))*power3(mh2)))/(48.*mh2*mstau22)) + Fin3(mstau22,mstau12,mh2,Q2)*((
     mu2*(-1 + power2(sa))*(invdmstau2mu*power2(mh2) + (-mstau12 + mu2)*power2(
     invdmstau2mu)*power2(mh2) - 2*(mstau12*power2(invdmstau1mu)*power2(mh2) +
     3*(mh2 + mstau12 - 2*invdmhH*mh2*mstau12)*power2(sa) + mstau12*mu2*power2(
     mh2)*power3(invdmstau1mu)) - 2*mstau12*mu2*power2(mh2)*power3(invdmstau2mu
     )))/(8.*mstau12*power2(mh2)) + (Al*ca*mu*sa*(-(invdmstau2mu*power2(mh2)) +
     2*mstau12*power2(invdmstau1mu)*power2(mh2) + (mstau12 - mu2)*power2(
     invdmstau2mu)*power2(mh2) + 6*invdmhH*mh2*mstau12*(1 - 2*power2(sa)) + 6*
     mh2*power2(sa) + 6*mstau12*power2(sa) + 2*mstau12*mu2*power2(mh2)*power3(
     invdmstau1mu) + 2*mstau12*mu2*power2(mh2)*power3(invdmstau2mu)))/(4.*
     mstau12*power2(mh2)) + (power2(Al)*power2(sa)*(-(invdmstau2mu*power2(mh2))
     + (mstau12 - mu2)*power2(invdmstau2mu)*power2(mh2) + 2*(mstau12*power2(
     invdmstau1mu)*power2(mh2) - 6*invdmhH*mh2*mstau12*(-1 + power2(sa)) + 3*(
     mh2 + mstau12)*power2(sa) + mstau12*mu2*power2(mh2)*power3(invdmstau1mu))
     + 2*mstau12*mu2*power2(mh2)*power3(invdmstau2mu)))/(8.*mstau12*power2(mh2)
     ) + DeltaInv(mstau22,mstau12,mh2)*((Al*ca*mu*sa*(invdmstau2mu*mh2*(-3*mh2*
     mstau12 - 2*mh2*mstau22 - 2*mstau12*mstau22 + mstau24 + power2(mh2)) + mh2
     *mu2*power2(invdmstau2mu)*(-3*mh2*mstau12 - 2*mh2*mstau22 - 2*mstau12*
     mstau22 + mstau24 + mstau12*mu2 + power2(mh2)) + power2(mh2)*(-3*
     invdmstau1mu*mstau12 + mstau12*(2*mstau22 + mu2)*power2(invdmstau1mu) - 6*
     power2(sa)) + 6*(3*mstau12*mstau22 - mstau24)*power2(sa) + mh2*(12*mstau22
     *power2(sa) + mstau12*(-3 - 3*invdmstau1mu*mstau22 + 2*invdmstau1mu*mu2 -
     mstau24*power2(invdmstau1mu) + mstau22*mu2*power2(invdmstau1mu) + 6*power2
     (sa))) - mstau12*power2(invdmstau1mu)*power3(mh2)))/(4.*mh2*mstau12) - (
     mu2*(-1 + power2(sa))*(invdmstau2mu*mh2*(-3*mh2*mstau12 - 2*mh2*mstau22 -
     2*mstau12*mstau22 + mstau24 + power2(mh2)) + mh2*mu2*power2(invdmstau2mu)*
     (-3*mh2*mstau12 - 2*mh2*mstau22 - 2*mstau12*mstau22 + mstau24 + mstau12*
     mu2 + power2(mh2)) + power2(mh2)*(-3*invdmstau1mu*mstau12 + mstau12*(2*
     mstau22 + mu2)*power2(invdmstau1mu) - 6*power2(sa)) + 6*(3*mstau12*mstau22
      - mstau24)*power2(sa) + mh2*(12*mstau22*power2(sa) + mstau12*(-3 - 3*
     invdmstau1mu*mstau22 + 2*invdmstau1mu*mu2 - mstau24*power2(invdmstau1mu) +
     mstau22*mu2*power2(invdmstau1mu) + 6*power2(sa))) - mstau12*power2(
     invdmstau1mu)*power3(mh2)))/(8.*mh2*mstau12) + (power2(Al)*power2(sa)*(
     invdmstau2mu*mh2*(-3*mh2*mstau12 - 2*mh2*mstau22 - 2*mstau12*mstau22 +
     mstau24 + power2(mh2)) + mh2*mu2*power2(invdmstau2mu)*(-3*mh2*mstau12 - 2*
     mh2*mstau22 - 2*mstau12*mstau22 + mstau24 + mstau12*mu2 + power2(mh2)) +
     power2(mh2)*(-3*invdmstau1mu*mstau12 + mstau12*(2*mstau22 + mu2)*power2(
     invdmstau1mu) - 6*power2(sa)) + 6*(3*mstau12*mstau22 - mstau24)*power2(sa)
     + mh2*(12*mstau22*power2(sa) + mstau12*(-3 - 3*invdmstau1mu*mstau22 + 2*
     invdmstau1mu*mu2 - mstau24*power2(invdmstau1mu) + mstau22*mu2*power2(
     invdmstau1mu) + 6*power2(sa))) - mstau12*power2(invdmstau1mu)*power3(mh2))
     )/(8.*mh2*mstau12))) - (mu2*DeltaInv(mH2,msntau2,mu2)*(-1 + power2(sa))*(-
     ((msntau2 - mu2)*power2(mu2)*(42 + 6*logmsntau2Q2*(-3 + logmu2Q2) - 18*
     logmu2Q2 + 3*power2(logmsntau2Q2) + 3*power2(logmu2Q2) + power2(Pi))) -
     mH2*mu2*(-3*mu2*(4*logmH2Q2*logmsntau2Q2 - 2*(-6 + logmu2Q2)*logmu2Q2 - 2*
     logmH2Q2*(3 + logmu2Q2) - 2*logmsntau2Q2*(3 + logmu2Q2) + power2(logmH2Q2)
     + power2(logmsntau2Q2)) + msntau2*(294 + 36*logmsntau2Q2*(-5 + logmu2Q2) -
     54*logmu2Q2 + 6*logmH2Q2*(4*logmsntau2Q2 - 3*(1 + logmu2Q2)) + 3*power2(
     logmH2Q2) + 30*power2(logmsntau2Q2) + 9*power2(logmu2Q2) + 7*power2(Pi)))
     - power2(mH2)*(3*mu2*(42 + 4*logmH2Q2*(-3 + 2*logmsntau2Q2 - logmu2Q2) + 2
     *logmsntau2Q2*(-15 + logmu2Q2) + 6*logmu2Q2 + 2*power2(logmH2Q2) + 5*
     power2(logmsntau2Q2) - power2(logmu2Q2) + power2(Pi)) + msntau2*(6*
     logmH2Q2*(-9 + 4*logmsntau2Q2 - logmu2Q2) + 6*logmsntau2Q2*(-15 + logmu2Q2
     ) + 9*power2(logmH2Q2) + 15*power2(logmsntau2Q2) + 4*(42 + power2(Pi)))) +
     (6*logmH2Q2*(-3 + 2*logmsntau2Q2 - logmu2Q2) + 6*logmsntau2Q2*(-9 +
     logmu2Q2) + 3*power2(logmH2Q2) + 9*power2(logmsntau2Q2) + 2*(42 + power2(
     Pi)))*power3(mH2)))/(48.*mH2*msntau2) + DeltaInv(mH2,mstau12,mu2)*((Al*ca*
     mu*mu2*sa*power2(invdmstau2mu)*((mstau12 - mu2)*power2(mu2)*(42 + 6*
     logmstau12Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*power2(logmstau12Q2) + 3*
     power2(logmu2Q2) + power2(Pi)) + mH2*mu2*(-3*mu2*(4*logmH2Q2*logmstau12Q2
     - 2*(-6 + logmu2Q2)*logmu2Q2 - 2*logmH2Q2*(3 + logmu2Q2) - 2*logmstau12Q2*
     (3 + logmu2Q2) + power2(logmH2Q2) + power2(logmstau12Q2)) + mstau12*(294 +
     36*logmstau12Q2*(-5 + logmu2Q2) - 54*logmu2Q2 + 6*logmH2Q2*(4*logmstau12Q2
      - 3*(1 + logmu2Q2)) + 3*power2(logmH2Q2) + 30*power2(logmstau12Q2) + 9*
     power2(logmu2Q2) + 7*power2(Pi))) + power2(mH2)*(3*mu2*(42 + 4*logmH2Q2*(-
     3 + 2*logmstau12Q2 - logmu2Q2) + 2*logmstau12Q2*(-15 + logmu2Q2) + 6*
     logmu2Q2 + 2*power2(logmH2Q2) + 5*power2(logmstau12Q2) - power2(logmu2Q2)
     + power2(Pi)) + mstau12*(6*logmH2Q2*(-9 + 4*logmstau12Q2 - logmu2Q2) + 6*
     logmstau12Q2*(-15 + logmu2Q2) + 9*power2(logmH2Q2) + 15*power2(
     logmstau12Q2) + 4*(42 + power2(Pi)))) - (6*logmH2Q2*(-3 + 2*logmstau12Q2 -
     logmu2Q2) + 6*logmstau12Q2*(-9 + logmu2Q2) + 3*power2(logmH2Q2) + 9*power2
     (logmstau12Q2) + 2*(42 + power2(Pi)))*power3(mH2)))/(24.*mstau12) + (mu2*
     power2(Al)*power2(invdmstau2mu)*(-1 + power2(sa))*((mstau12 - mu2)*power2(
     mu2)*(42 + 6*logmstau12Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*power2(
     logmstau12Q2) + 3*power2(logmu2Q2) + power2(Pi)) + mH2*mu2*(-3*mu2*(4*
     logmH2Q2*logmstau12Q2 - 2*(-6 + logmu2Q2)*logmu2Q2 - 2*logmH2Q2*(3 +
     logmu2Q2) - 2*logmstau12Q2*(3 + logmu2Q2) + power2(logmH2Q2) + power2(
     logmstau12Q2)) + mstau12*(294 + 36*logmstau12Q2*(-5 + logmu2Q2) - 54*
     logmu2Q2 + 6*logmH2Q2*(4*logmstau12Q2 - 3*(1 + logmu2Q2)) + 3*power2(
     logmH2Q2) + 30*power2(logmstau12Q2) + 9*power2(logmu2Q2) + 7*power2(Pi)))
     + power2(mH2)*(3*mu2*(42 + 4*logmH2Q2*(-3 + 2*logmstau12Q2 - logmu2Q2) + 2
     *logmstau12Q2*(-15 + logmu2Q2) + 6*logmu2Q2 + 2*power2(logmH2Q2) + 5*
     power2(logmstau12Q2) - power2(logmu2Q2) + power2(Pi)) + mstau12*(6*
     logmH2Q2*(-9 + 4*logmstau12Q2 - logmu2Q2) + 6*logmstau12Q2*(-15 + logmu2Q2
     ) + 9*power2(logmH2Q2) + 15*power2(logmstau12Q2) + 4*(42 + power2(Pi)))) -
     (6*logmH2Q2*(-3 + 2*logmstau12Q2 - logmu2Q2) + 6*logmstau12Q2*(-9 +
     logmu2Q2) + 3*power2(logmH2Q2) + 9*power2(logmstau12Q2) + 2*(42 + power2(
     Pi)))*power3(mH2)))/(48.*mstau12) + (mu2*(1 + (-1 + mH2*mu2*power2(
     invdmstau2mu))*power2(sa))*(-((mstau12 - mu2)*power2(mu2)*(42 + 6*
     logmstau12Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*power2(logmstau12Q2) + 3*
     power2(logmu2Q2) + power2(Pi))) - mH2*mu2*(-3*mu2*(4*logmH2Q2*logmstau12Q2
      - 2*(-6 + logmu2Q2)*logmu2Q2 - 2*logmH2Q2*(3 + logmu2Q2) - 2*logmstau12Q2
     *(3 + logmu2Q2) + power2(logmH2Q2) + power2(logmstau12Q2)) + mstau12*(294
     + 36*logmstau12Q2*(-5 + logmu2Q2) - 54*logmu2Q2 + 6*logmH2Q2*(4*
     logmstau12Q2 - 3*(1 + logmu2Q2)) + 3*power2(logmH2Q2) + 30*power2(
     logmstau12Q2) + 9*power2(logmu2Q2) + 7*power2(Pi))) - power2(mH2)*(3*mu2*(
     42 + 4*logmH2Q2*(-3 + 2*logmstau12Q2 - logmu2Q2) + 2*logmstau12Q2*(-15 +
     logmu2Q2) + 6*logmu2Q2 + 2*power2(logmH2Q2) + 5*power2(logmstau12Q2) -
     power2(logmu2Q2) + power2(Pi)) + mstau12*(6*logmH2Q2*(-9 + 4*logmstau12Q2
     - logmu2Q2) + 6*logmstau12Q2*(-15 + logmu2Q2) + 9*power2(logmH2Q2) + 15*
     power2(logmstau12Q2) + 4*(42 + power2(Pi)))) + (6*logmH2Q2*(-3 + 2*
     logmstau12Q2 - logmu2Q2) + 6*logmstau12Q2*(-9 + logmu2Q2) + 3*power2(
     logmH2Q2) + 9*power2(logmstau12Q2) + 2*(42 + power2(Pi)))*power3(mH2)))/(
     48.*mH2*mstau12)) + DeltaInv(mstau22,mH2,mu2)*((Al*ca*mu*mu2*sa*power2(
     invdmstau1mu)*((mstau22 - mu2)*power2(mu2)*(42 + 6*logmstau22Q2*(-3 +
     logmu2Q2) - 18*logmu2Q2 + 3*power2(logmstau22Q2) + 3*power2(logmu2Q2) +
     power2(Pi)) + mH2*mu2*(-3*mu2*(4*logmH2Q2*logmstau22Q2 - 2*(-6 + logmu2Q2)
     *logmu2Q2 - 2*logmH2Q2*(3 + logmu2Q2) - 2*logmstau22Q2*(3 + logmu2Q2) +
     power2(logmH2Q2) + power2(logmstau22Q2)) + mstau22*(294 + 36*logmstau22Q2*
     (-5 + logmu2Q2) - 54*logmu2Q2 + 6*logmH2Q2*(4*logmstau22Q2 - 3*(1 +
     logmu2Q2)) + 3*power2(logmH2Q2) + 30*power2(logmstau22Q2) + 9*power2(
     logmu2Q2) + 7*power2(Pi))) + power2(mH2)*(3*mu2*(42 + 4*logmH2Q2*(-3 + 2*
     logmstau22Q2 - logmu2Q2) + 2*logmstau22Q2*(-15 + logmu2Q2) + 6*logmu2Q2 +
     2*power2(logmH2Q2) + 5*power2(logmstau22Q2) - power2(logmu2Q2) + power2(Pi
     )) + mstau22*(6*logmH2Q2*(-9 + 4*logmstau22Q2 - logmu2Q2) + 6*logmstau22Q2
     *(-15 + logmu2Q2) + 9*power2(logmH2Q2) + 15*power2(logmstau22Q2) + 4*(42 +
     power2(Pi)))) - (6*logmH2Q2*(-3 + 2*logmstau22Q2 - logmu2Q2) + 6*
     logmstau22Q2*(-9 + logmu2Q2) + 3*power2(logmH2Q2) + 9*power2(logmstau22Q2)
     + 2*(42 + power2(Pi)))*power3(mH2)))/(24.*mstau22) + (mu2*power2(Al)*
     power2(invdmstau1mu)*(-1 + power2(sa))*((mstau22 - mu2)*power2(mu2)*(42 +
     6*logmstau22Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*power2(logmstau22Q2) + 3*
     power2(logmu2Q2) + power2(Pi)) + mH2*mu2*(-3*mu2*(4*logmH2Q2*logmstau22Q2
     - 2*(-6 + logmu2Q2)*logmu2Q2 - 2*logmH2Q2*(3 + logmu2Q2) - 2*logmstau22Q2*
     (3 + logmu2Q2) + power2(logmH2Q2) + power2(logmstau22Q2)) + mstau22*(294 +
     36*logmstau22Q2*(-5 + logmu2Q2) - 54*logmu2Q2 + 6*logmH2Q2*(4*logmstau22Q2
      - 3*(1 + logmu2Q2)) + 3*power2(logmH2Q2) + 30*power2(logmstau22Q2) + 9*
     power2(logmu2Q2) + 7*power2(Pi))) + power2(mH2)*(3*mu2*(42 + 4*logmH2Q2*(-
     3 + 2*logmstau22Q2 - logmu2Q2) + 2*logmstau22Q2*(-15 + logmu2Q2) + 6*
     logmu2Q2 + 2*power2(logmH2Q2) + 5*power2(logmstau22Q2) - power2(logmu2Q2)
     + power2(Pi)) + mstau22*(6*logmH2Q2*(-9 + 4*logmstau22Q2 - logmu2Q2) + 6*
     logmstau22Q2*(-15 + logmu2Q2) + 9*power2(logmH2Q2) + 15*power2(
     logmstau22Q2) + 4*(42 + power2(Pi)))) - (6*logmH2Q2*(-3 + 2*logmstau22Q2 -
     logmu2Q2) + 6*logmstau22Q2*(-9 + logmu2Q2) + 3*power2(logmH2Q2) + 9*power2
     (logmstau22Q2) + 2*(42 + power2(Pi)))*power3(mH2)))/(48.*mstau22) + (mu2*(
     1 + (-1 + mH2*mu2*power2(invdmstau1mu))*power2(sa))*(-((mstau22 - mu2)*
     power2(mu2)*(42 + 6*logmstau22Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*power2(
     logmstau22Q2) + 3*power2(logmu2Q2) + power2(Pi))) - mH2*mu2*(-3*mu2*(4*
     logmH2Q2*logmstau22Q2 - 2*(-6 + logmu2Q2)*logmu2Q2 - 2*logmH2Q2*(3 +
     logmu2Q2) - 2*logmstau22Q2*(3 + logmu2Q2) + power2(logmH2Q2) + power2(
     logmstau22Q2)) + mstau22*(294 + 36*logmstau22Q2*(-5 + logmu2Q2) - 54*
     logmu2Q2 + 6*logmH2Q2*(4*logmstau22Q2 - 3*(1 + logmu2Q2)) + 3*power2(
     logmH2Q2) + 30*power2(logmstau22Q2) + 9*power2(logmu2Q2) + 7*power2(Pi)))
     - power2(mH2)*(3*mu2*(42 + 4*logmH2Q2*(-3 + 2*logmstau22Q2 - logmu2Q2) + 2
     *logmstau22Q2*(-15 + logmu2Q2) + 6*logmu2Q2 + 2*power2(logmH2Q2) + 5*
     power2(logmstau22Q2) - power2(logmu2Q2) + power2(Pi)) + mstau22*(6*
     logmH2Q2*(-9 + 4*logmstau22Q2 - logmu2Q2) + 6*logmstau22Q2*(-15 + logmu2Q2
     ) + 9*power2(logmH2Q2) + 15*power2(logmstau22Q2) + 4*(42 + power2(Pi)))) +
     (6*logmH2Q2*(-3 + 2*logmstau22Q2 - logmu2Q2) + 6*logmstau22Q2*(-9 +
     logmu2Q2) + 3*power2(logmH2Q2) + 9*power2(logmstau22Q2) + 2*(42 + power2(
     Pi)))*power3(mH2)))/(48.*mH2*mstau22)) + Fin3(mstau22,mH2,mstau12,Q2)*((
     power2(Al)*(-1 + power2(sa))*(4*mh2*mH2*(-3 + 4*invdmhH*mH2)*mstau12*
     power2(sa) - 4*(-1 + invdmhH*mH2)*mstau12*power2(mH2)*power2(sa) + power2(
     mh2)*(6*mH2*(-1 + power2(sa)) + 6*mstau12*(-1 + power2(sa)) + power2(mH2)*
     (invdmstau2mu - 2*mstau12*(1 + invdmstau1mu*mu2)*power2(invdmstau1mu) + (-
     mstau12 + mu2)*power2(invdmstau2mu) - 2*mstau12*mu2*power3(invdmstau2mu)))
     ))/(8.*mstau12*power2(mh2)*power2(mH2)) + (Al*ca*mu*sa*(2*mh2*mH2*(-3 + 4*
     invdmhH*mH2)*mstau12*(-1 + 2*power2(sa)) - 2*(-1 + invdmhH*mH2)*mstau12*
     power2(mH2)*(-1 + 2*power2(sa)) + power2(mh2)*(6*mH2*(-1 + power2(sa)) + 6
     *mstau12*(-1 + power2(sa)) + power2(mH2)*(invdmstau2mu - 2*mstau12*(1 +
     invdmstau1mu*mu2)*power2(invdmstau1mu) + (-mstau12 + mu2)*power2(
     invdmstau2mu) - 2*mstau12*mu2*power3(invdmstau2mu)))))/(4.*mstau12*power2(
     mh2)*power2(mH2)) + (mu2*power2(sa)*(-4*mh2*mH2*(-3 + 4*invdmhH*mH2)*
     mstau12*(-1 + power2(sa)) + 4*(-1 + invdmhH*mH2)*mstau12*power2(mH2)*(-1 +
     power2(sa)) + power2(mh2)*(-6*mH2*(-1 + power2(sa)) - 6*mstau12*(-1 +
     power2(sa)) + power2(mH2)*(-invdmstau2mu + 2*mstau12*(1 + invdmstau1mu*mu2
     )*power2(invdmstau1mu) + (mstau12 - mu2)*power2(invdmstau2mu) + 2*mstau12*
     mu2*power3(invdmstau2mu)))))/(8.*mstau12*power2(mh2)*power2(mH2)) +
     DeltaInv(mstau22,mH2,mstau12)*(-(Al*ca*mu*sa*(power2(mH2)*(-3*invdmstau1mu
     *mstau12 - invdmstau2mu*(3*mstau12 + 2*mstau22) + mstau12*(2*mstau22 + mu2
     )*power2(invdmstau1mu) - (3*mstau12 + 2*mstau22)*mu2*power2(invdmstau2mu)
     + 6*(-1 + power2(sa))) - 6*(3*mstau12*mstau22 - mstau24)*(-1 + power2(sa))
     + mH2*(12*mstau22 + invdmstau2mu*mstau24 + mstau24*mu2*power2(invdmstau2mu
     ) + mstau12*(3 - 2*invdmstau2mu*mstau22 + invdmstau1mu*(-3*mstau22 + 2*mu2
     ) + (-mstau24 + mstau22*mu2)*power2(invdmstau1mu) + mu2*(-2*mstau22 + mu2)
     *power2(invdmstau2mu) - 6*power2(sa)) - 12*mstau22*power2(sa)) + (
     invdmstau2mu - mstau12*power2(invdmstau1mu) + mu2*power2(invdmstau2mu))*
     power3(mH2)))/(4.*mH2*mstau12) - (power2(Al)*(-1 + power2(sa))*(power2(mH2
     )*(-3*invdmstau1mu*mstau12 - invdmstau2mu*(3*mstau12 + 2*mstau22) +
     mstau12*(2*mstau22 + mu2)*power2(invdmstau1mu) - (3*mstau12 + 2*mstau22)*
     mu2*power2(invdmstau2mu) + 6*(-1 + power2(sa))) - 6*(3*mstau12*mstau22 -
     mstau24)*(-1 + power2(sa)) + mH2*(12*mstau22 + invdmstau2mu*mstau24 +
     mstau24*mu2*power2(invdmstau2mu) + mstau12*(3 - 2*invdmstau2mu*mstau22 +
     invdmstau1mu*(-3*mstau22 + 2*mu2) + (-mstau24 + mstau22*mu2)*power2(
     invdmstau1mu) + mu2*(-2*mstau22 + mu2)*power2(invdmstau2mu) - 6*power2(sa)
     ) - 12*mstau22*power2(sa)) + (invdmstau2mu - mstau12*power2(invdmstau1mu)
     + mu2*power2(invdmstau2mu))*power3(mH2)))/(8.*mH2*mstau12) + (mu2*power2(
     sa)*(power2(mH2)*(-3*invdmstau1mu*mstau12 - invdmstau2mu*(3*mstau12 + 2*
     mstau22) + mstau12*(2*mstau22 + mu2)*power2(invdmstau1mu) - (3*mstau12 + 2
     *mstau22)*mu2*power2(invdmstau2mu) + 6*(-1 + power2(sa))) - 6*(3*mstau12*
     mstau22 - mstau24)*(-1 + power2(sa)) + mH2*(12*mstau22 + invdmstau2mu*
     mstau24 + mstau24*mu2*power2(invdmstau2mu) + mstau12*(3 - 2*invdmstau2mu*
     mstau22 + invdmstau1mu*(-3*mstau22 + 2*mu2) + (-mstau24 + mstau22*mu2)*
     power2(invdmstau1mu) + mu2*(-2*mstau22 + mu2)*power2(invdmstau2mu) - 6*
     power2(sa)) - 12*mstau22*power2(sa)) + (invdmstau2mu - mstau12*power2(
     invdmstau1mu) + mu2*power2(invdmstau2mu))*power3(mH2)))/(8.*mH2*mstau12)))
     + (Al*mu*(ca*cb*mA2*mC2*sa*(-(mH2*(-1 + invdmhH*mH2)*mstau12*mstau22*(3*(-
     7 + 4*logmH2Q2)*mH2 + 3*(-7 + 4*logmstau12Q2)*mstau12 + 3*(-7 + 4*
     logmstau22Q2)*mstau22 + 2*invdmhH*power2(mH2)*(21 - 18*logmH2Q2 + 6*power2
     (logmH2Q2) + power2(Pi)))*(-1 + 2*power2(sa))) + mh2*mH2*mstau22*(12*
     mstau22*(-6*(3 + logmh2Q2)*logmstau22Q2 + 6*logmstau12Q2*(-9 + logmh2Q2 +
     2*logmstau22Q2) + 9*power2(logmstau12Q2) + 3*power2(logmstau22Q2) + 2*(42
     + power2(Pi)))*power2(sa) + 3*invdmhH*(-7 + 4*logmstau12Q2)*power2(mstau12
     )*(-1 + 2*power2(sa)) + mstau12*(21 + 6*(53 + 12*logmh2Q2*(-2 +
     logmstau12Q2) - 36*logmstau12Q2 + 6*power2(logmh2Q2) + 6*power2(
     logmstau12Q2) + 2*power2(Pi))*power2(sa) + 12*invdmhH*mH2*(-4 + 5*invdmhH*
     mH2)*power2(logmH2Q2)*(-1 + 2*power2(sa)) - 12*logmH2Q2*(-1 - 8*invdmhH*
     mH2 + 11*power2(invdmhH)*power2(mH2))*(-1 + 2*power2(sa)) + 2*power2(
     invdmhH)*power2(mH2)*(69 + 5*power2(Pi))*(-1 + 2*power2(sa)) - invdmhH*(3*
     (7 - 4*logmstau22Q2)*mstau22 + 8*mH2*(12 + power2(Pi)))*(-1 + 2*power2(sa)
     ))) + power2(mh2)*(12*mstau22*(mstau12*(30 + 6*logmH2Q2*(-2 + logmstau12Q2
     ) - 18*logmstau12Q2 + 3*power2(logmH2Q2) + 3*power2(logmstau12Q2) + power2
     (Pi)) + mstau22*(-6*(3 + logmH2Q2)*logmstau22Q2 + 6*logmstau12Q2*(-9 +
     logmH2Q2 + 2*logmstau22Q2) + 9*power2(logmstau12Q2) + 3*power2(
     logmstau22Q2) + 2*(42 + power2(Pi))))*(-1 + power2(sa)) + 3*mH2*mstau22*(-
     168 + 7*invdmhH*mstau12 - 24*invdmhH*logmh2Q2*mstau12 + 8*invdmstau*
     logmh2Q2*logmstau22Q2*mstau12 - 8*invdmstau*invdmstau2mu*logmh2Q2*
     logmstau22Q2*mstau12*mu2 - 8*invdmstau*invdmstau1mu*logmh2Q2*logmu2Q2*
     mstau12*mu2 + 8*invdmstau*invdmstau2mu*logmh2Q2*logmu2Q2*mstau12*mu2 + 12*
     invdmhH*mstau12*power2(logmh2Q2) - 4*power2(Pi) + 336*power2(sa) - 72*
     logmh2Q2*power2(sa) - 14*invdmhH*mstau12*power2(sa) + 48*invdmhH*logmh2Q2*
     mstau12*power2(sa) + 12*power2(logmh2Q2)*power2(sa) - 24*invdmhH*mstau12*
     power2(logmh2Q2)*power2(sa) + 8*power2(Pi)*power2(sa) + 12*power2(
     logmstau12Q2)*(-1 + 2*power2(sa)) + 8*logmstau12Q2*(9 + invdmstau*logmh2Q2
     *mstau12*(-1 + invdmstau1mu*mu2) + 3*(-6 + logmh2Q2)*power2(sa)) + 12*
     power2(logmH2Q2)*(-1 + power2(sa) + invdmhH*mstau12*(-1 + 2*power2(sa))) -
     4*logmH2Q2*(-2*(9 + invdmstau*mstau12*((invdmstau1mu - invdmstau2mu)*
     logmu2Q2*mu2 + logmstau22Q2*(-1 + invdmstau2mu*mu2)) - 9*power2(sa)) + 2*
     logmstau12Q2*(3 + invdmstau*mstau12*(-1 + invdmstau1mu*mu2) - 3*power2(sa)
     ) + 5*invdmhH*mstau12*(-1 + 2*power2(sa)))) + 2*power2(mH2)*(invdmstau2mu*
     mstau22*(6*logmH2Q2*(-3 + 2*logmstau12Q2 - logmstau22Q2) + 6*logmstau12Q2*
     (-9 + logmstau22Q2) + 3*power2(logmH2Q2) + 9*power2(logmstau12Q2) + 2*(42
     + power2(Pi))) + 2*mstau22*mu2*power2(invdmstau2mu)*(3*logmH2Q2*(-6 + 4*
     logmstau12Q2 - logmstau22Q2 - logmu2Q2) + 3*logmstau12Q2*(-18 +
     logmstau22Q2 + logmu2Q2) + 3*power2(logmH2Q2) + 9*power2(logmstau12Q2) + 2
     *(42 + power2(Pi))) + mstau12*(power2(invdmstau1mu)*(-(mstau22*(42 + 6*
     logmH2Q2*(-3 + logmstau12Q2) - 18*logmstau12Q2 + 3*power2(logmH2Q2) + 3*
     power2(logmstau12Q2) + power2(Pi))) + mu2*(6*logmH2Q2*(-3 + 2*logmstau22Q2
      - logmu2Q2) + 6*logmstau22Q2*(-9 + logmu2Q2) + 3*power2(logmH2Q2) + 9*
     power2(logmstau22Q2) + 2*(42 + power2(Pi)))) - 4*mstau22*power2(invdmhH)*(
     12 - 12*logmH2Q2 + 6*power2(logmH2Q2) + power2(Pi))*(-1 + 2*power2(sa)))))
     - 2*mH2*(invdmstau2mu*mstau22*(6*logmh2Q2*(-3 + 2*logmstau12Q2 -
     logmstau22Q2) + 6*logmstau12Q2*(-9 + logmstau22Q2) + 3*power2(logmh2Q2) +
     9*power2(logmstau12Q2) + 2*(42 + power2(Pi))) + 2*mstau22*mu2*power2(
     invdmstau2mu)*(3*logmh2Q2*(-6 + 4*logmstau12Q2 - logmstau22Q2 - logmu2Q2)
     + 3*logmstau12Q2*(-18 + logmstau22Q2 + logmu2Q2) + 3*power2(logmh2Q2) + 9*
     power2(logmstau12Q2) + 2*(42 + power2(Pi))) + mstau12*power2(invdmstau1mu)
     *(-(mstau22*(42 + 6*logmh2Q2*(-3 + logmstau12Q2) - 18*logmstau12Q2 + 3*
     power2(logmh2Q2) + 3*power2(logmstau12Q2) + power2(Pi))) + mu2*(6*logmh2Q2
     *(-3 + 2*logmstau22Q2 - logmu2Q2) + 6*logmstau22Q2*(-9 + logmu2Q2) + 3*
     power2(logmh2Q2) + 9*power2(logmstau22Q2) + 2*(42 + power2(Pi)))))*power3(
     mh2)) + 12*mstau12*mstau22*sb*(mA2*mH2*(-1 + logmstau12Q2 + invdmstau*
     logmstau12Q2*mstau22 - invdmstau*logmstau22Q2*mstau22)*power2(mh2) - mC2*
     mH2*(-1 + logmstau12Q2 + invdmstau*logmstau12Q2*mstau22 - invdmstau*
     logmstau22Q2*mstau22)*power2(mh2) + mA2*mC2*((-3 + 2*logmH2Q2)*(-1 +
     invdmhH*mH2)*(-1 + logmstau12Q2 + invdmstau*logmstau12Q2*mstau22 -
     invdmstau*logmstau22Q2*mstau22)*power2(mH2)*(-1 + power2(sa))*power2(sa) -
     mh2*mH2*(-1 + logmstau12Q2 + invdmstau*logmstau12Q2*mstau22 - invdmstau*
     logmstau22Q2*mstau22)*power2(sa)*(-5 - 8*invdmhH*mH2*(-1 + power2(sa)) + 2
     *logmH2Q2*(-3 + 4*invdmhH*mH2)*(-1 + power2(sa)) + 2*power2(sa)) + power2(
     mh2)*(-((-1 + power2(sa))*(-3 + (3 + invdmhH*(-5 + 6*logmh2Q2)*mH2)*power2
     (sa))) + invdmstau*(-3*logmstau22Q2*mstau22 - 2*logmA2Q2*mH2*((
     invdmstau1mu - invdmstau2mu)*logmu2Q2*mu2 + logmstau22Q2*(-1 +
     invdmstau2mu*mu2)) - 2*logmH2Q2*mH2*((invdmstau1mu - invdmstau2mu)*
     logmu2Q2*mu2 + logmstau22Q2*(-1 + invdmstau2mu*mu2))*(-1 + power2(sa)) - 2
     *logmh2Q2*logmstau22Q2*mH2*power2(sa) + 6*logmstau22Q2*mstau22*power2(sa)
     - 5*invdmhH*logmstau22Q2*mH2*mstau22*power2(sa) + 6*invdmhH*logmh2Q2*
     logmstau22Q2*mH2*mstau22*power2(sa) + 2*invdmstau2mu*logmh2Q2*logmstau22Q2
     *mH2*mu2*power2(sa) + 2*invdmstau1mu*logmh2Q2*logmu2Q2*mH2*mu2*power2(sa)
     - 2*invdmstau2mu*logmh2Q2*logmu2Q2*mH2*mu2*power2(sa) - 3*logmstau22Q2*
     mstau22*power4(sa) + 5*invdmhH*logmstau22Q2*mH2*mstau22*power4(sa) - 6*
     invdmhH*logmh2Q2*logmstau22Q2*mH2*mstau22*power4(sa)) + logmstau12Q2*((-1
     + power2(sa))*(-3 + (3 + invdmhH*(-5 + 6*logmh2Q2)*mH2)*power2(sa)) +
     invdmstau*(3*mstau22 + 2*logmA2Q2*mH2*(-1 + invdmstau1mu*mu2) + 2*logmH2Q2
     *mH2*(-1 + invdmstau1mu*mu2)*(-1 + power2(sa)) + 2*logmh2Q2*mH2*power2(sa)
     - 6*mstau22*power2(sa) + 5*invdmhH*mH2*mstau22*power2(sa) - 6*invdmhH*
     logmh2Q2*mH2*mstau22*power2(sa) - 2*invdmstau1mu*logmh2Q2*mH2*mu2*power2(
     sa) + 3*mstau22*power4(sa) - 5*invdmhH*mH2*mstau22*power4(sa) + 6*invdmhH*
     logmh2Q2*mH2*mstau22*power4(sa))))))))/(48.*cb*mA2*mC2*mH2*mstau12*mstau22
     *power2(mh2)) - (power2(Al)*(-2*mC2*mH2*msntau2*mstau22*power2(mh2)*(-(
     mstau12*(36 + 6*logmA2Q2*(-2 + logmstau12Q2) + 6*invdmstau*logmstau22Q2*
     mstau22 - 6*logmstau12Q2*(4 + invdmstau*mstau22) + 3*power2(logmA2Q2) + 3*
     power2(logmstau12Q2) + power2(Pi))) - mstau22*(-6*(3 + logmA2Q2)*
     logmstau22Q2 + 6*logmstau12Q2*(-9 + logmA2Q2 + 2*logmstau22Q2) + 9*power2(
     logmstau12Q2) + 3*power2(logmstau22Q2) + 2*(42 + power2(Pi)))) + mC2*mH2*
     msntau2*power2(mA2)*power2(mh2)*(invdmstau2mu*mstau22*(6*logmA2Q2*(-3 + 2*
     logmstau12Q2 - logmstau22Q2) + 6*logmstau12Q2*(-9 + logmstau22Q2) + 3*
     power2(logmA2Q2) + 9*power2(logmstau12Q2) + 2*(42 + power2(Pi))) + 2*
     mstau22*mu2*power2(invdmstau2mu)*(3*logmA2Q2*(-6 + 4*logmstau12Q2 -
     logmstau22Q2 - logmu2Q2) + 3*logmstau12Q2*(-18 + logmstau22Q2 + logmu2Q2)
     + 3*power2(logmA2Q2) + 9*power2(logmstau12Q2) + 2*(42 + power2(Pi))) +
     mstau12*power2(invdmstau1mu)*(-(mstau22*(42 + 6*logmA2Q2*(-3 +
     logmstau12Q2) - 18*logmstau12Q2 + 3*power2(logmA2Q2) + 3*power2(
     logmstau12Q2) + power2(Pi))) + mu2*(6*logmA2Q2*(-3 + 2*logmstau22Q2 -
     logmu2Q2) + 6*logmstau22Q2*(-9 + logmu2Q2) + 3*power2(logmA2Q2) + 9*power2
     (logmstau22Q2) + 2*(42 + power2(Pi))))) + mA2*(-2*mH2*mstau12*mstau22*
     power2(mh2)*(msntau2*(36 + 6*logmC2Q2*(-2 + logmsntau2Q2) - 18*
     logmsntau2Q2 - 6*logmstau12Q2 - 6*invdmstau*logmstau12Q2*mstau22 + 6*
     invdmstau*logmstau22Q2*mstau22 + 3*power2(logmC2Q2) + 3*power2(
     logmsntau2Q2) + power2(Pi)) + mstau22*(-6*(3 + logmC2Q2)*logmstau22Q2 + 6*
     logmsntau2Q2*(-9 + logmC2Q2 + 2*logmstau22Q2) + 9*power2(logmsntau2Q2) + 3
     *power2(logmstau22Q2) + 2*(42 + power2(Pi)))) + 2*mH2*mstau12*power2(mC2)*
     power2(mh2)*(invdmstau2mu*mstau22*(6*logmC2Q2*(-3 + 2*logmsntau2Q2 -
     logmstau22Q2) + 6*logmsntau2Q2*(-9 + logmstau22Q2) + 3*power2(logmC2Q2) +
     9*power2(logmsntau2Q2) + 2*(42 + power2(Pi))) + 2*mstau22*mu2*power2(
     invdmstau2mu)*(3*logmC2Q2*(-6 + 4*logmsntau2Q2 - logmstau22Q2 - logmu2Q2)
     + 3*logmsntau2Q2*(-18 + logmstau22Q2 + logmu2Q2) + 3*power2(logmC2Q2) + 9*
     power2(logmsntau2Q2) + 2*(42 + power2(Pi))) + msntau2*power2(invdmsntau2mu
     )*(-(mstau22*(42 + 6*logmC2Q2*(-3 + logmsntau2Q2) - 18*logmsntau2Q2 + 3*
     power2(logmC2Q2) + 3*power2(logmsntau2Q2) + power2(Pi))) + mu2*(6*logmC2Q2
     *(-3 + 2*logmstau22Q2 - logmu2Q2) + 6*logmstau22Q2*(-9 + logmu2Q2) + 3*
     power2(logmC2Q2) + 9*power2(logmstau22Q2) + 2*(42 + power2(Pi))))) + mC2*(
     mH2*(-1 + invdmhH*mH2)*msntau2*mstau12*mstau22*(3*((-7 + 4*logmstau12Q2)*
     mstau12 + (-7 + 4*logmstau22Q2)*mstau22) + 3*mH2*(5 + 12*invdmstau*
     logmstau22Q2*mstau22 - 12*logmstau12Q2*(1 + invdmstau*mstau22) + logmH2Q2*
     (-4 - 8*invdmstau*logmstau22Q2*mstau22 + 8*logmstau12Q2*(1 + invdmstau*
     mstau22))) + 2*invdmhH*power2(mH2)*(21 - 18*logmH2Q2 + 6*power2(logmH2Q2)
     + power2(Pi)))*(-1 + power2(sa))*power2(sa) - mh2*mH2*msntau2*mstau22*
     power2(sa)*(3*invdmhH*(-7 + 4*logmstau12Q2)*power2(mstau12)*(-1 + power2(
     sa)) + 6*mstau22*(-6*(3 + logmh2Q2)*logmstau22Q2 + 6*logmstau12Q2*(-9 +
     logmh2Q2 + 2*logmstau22Q2) + 9*power2(logmstau12Q2) + 3*power2(
     logmstau22Q2) + 2*(42 + power2(Pi)))*power2(sa) + mstau12*(81 + 21*invdmhH
     *mstau22 - 12*invdmhH*logmstau22Q2*mstau22 + 60*invdmstau*logmstau22Q2*
     mstau22 - 96*invdmhH*invdmstau*logmstau22Q2*mH2*mstau22 - 138*power2(
     invdmhH)*power2(mH2) + 8*invdmhH*mH2*power2(Pi) - 10*power2(invdmhH)*
     power2(mH2)*power2(Pi) + 12*invdmhH*mH2*(-4 + 5*invdmhH*mH2)*power2(
     logmH2Q2)*(-1 + power2(sa)) + 12*logmH2Q2*(7 + 6*invdmstau*logmstau22Q2*
     mstau22 - 8*invdmhH*invdmstau*logmstau22Q2*mH2*mstau22 + 2*logmstau12Q2*(-
     3 + 4*invdmhH*mH2)*(1 + invdmstau*mstau22) - 11*power2(invdmhH)*power2(mH2
     ))*(-1 + power2(sa)) + 135*power2(sa) - 72*logmh2Q2*power2(sa) - 21*
     invdmhH*mstau22*power2(sa) + 12*invdmhH*logmstau22Q2*mstau22*power2(sa) -
     24*invdmstau*logmstau22Q2*mstau22*power2(sa) + 96*invdmhH*invdmstau*
     logmstau22Q2*mH2*mstau22*power2(sa) + 18*power2(logmh2Q2)*power2(sa) + 18*
     power2(logmstau12Q2)*power2(sa) + 138*power2(invdmhH)*power2(mH2)*power2(
     sa) + 6*power2(Pi)*power2(sa) - 8*invdmhH*mH2*power2(Pi)*power2(sa) + 10*
     power2(invdmhH)*power2(mH2)*power2(Pi)*power2(sa) - 12*logmstau12Q2*(5 +
     invdmstau*mstau22*(5 - 2*power2(sa)) + 8*invdmhH*mH2*(1 + invdmstau*
     mstau22)*(-1 + power2(sa)) + 7*power2(sa) - 3*logmh2Q2*power2(sa)))) + mH2
     *msntau2*(invdmstau2mu*mstau22*(6*logmh2Q2*(-3 + 2*logmstau12Q2 -
     logmstau22Q2) + 6*logmstau12Q2*(-9 + logmstau22Q2) + 3*power2(logmh2Q2) +
     9*power2(logmstau12Q2) + 2*(42 + power2(Pi))) + 2*mstau22*mu2*power2(
     invdmstau2mu)*(3*logmh2Q2*(-6 + 4*logmstau12Q2 - logmstau22Q2 - logmu2Q2)
     + 3*logmstau12Q2*(-18 + logmstau22Q2 + logmu2Q2) + 3*power2(logmh2Q2) + 9*
     power2(logmstau12Q2) + 2*(42 + power2(Pi))) + mstau12*power2(invdmstau1mu)
     *(-(mstau22*(42 + 6*logmh2Q2*(-3 + logmstau12Q2) - 18*logmstau12Q2 + 3*
     power2(logmh2Q2) + 3*power2(logmstau12Q2) + power2(Pi))) + mu2*(6*logmh2Q2
     *(-3 + 2*logmstau22Q2 - logmu2Q2) + 6*logmstau22Q2*(-9 + logmu2Q2) + 3*
     power2(logmh2Q2) + 9*power2(logmstau22Q2) + 2*(42 + power2(Pi)))))*power2(
     sa)*power3(mh2) + power2(mh2)*(-(msntau2*power2(mH2)*(-1 + power2(sa))*(
     invdmstau2mu*mstau22*(6*logmH2Q2*(-3 + 2*logmstau12Q2 - logmstau22Q2) + 6*
     logmstau12Q2*(-9 + logmstau22Q2) + 3*power2(logmH2Q2) + 9*power2(
     logmstau12Q2) + 2*(42 + power2(Pi))) + 2*mstau22*mu2*power2(invdmstau2mu)*
     (3*logmH2Q2*(-6 + 4*logmstau12Q2 - logmstau22Q2 - logmu2Q2) + 3*
     logmstau12Q2*(-18 + logmstau22Q2 + logmu2Q2) + 3*power2(logmH2Q2) + 9*
     power2(logmstau12Q2) + 2*(42 + power2(Pi))) - mstau12*(power2(invdmstau1mu
     )*(mstau22*(42 + 6*logmH2Q2*(-3 + logmstau12Q2) - 18*logmstau12Q2 + 3*
     power2(logmH2Q2) + 3*power2(logmstau12Q2) + power2(Pi)) - mu2*(6*logmH2Q2*
     (-3 + 2*logmstau22Q2 - logmu2Q2) + 6*logmstau22Q2*(-9 + logmu2Q2) + 3*
     power2(logmH2Q2) + 9*power2(logmstau22Q2) + 2*(42 + power2(Pi)))) + 8*
     mstau22*power2(invdmhH)*(12 - 12*logmH2Q2 + 6*power2(logmH2Q2) + power2(Pi
     ))*power2(sa)))) - 6*msntau2*mstau22*(mstau12*(36 + 6*logmH2Q2*(-2 +
     logmstau12Q2) + 6*invdmstau*logmstau22Q2*mstau22 - 6*logmstau12Q2*(4 +
     invdmstau*mstau22) + 3*power2(logmH2Q2) + 3*power2(logmstau12Q2) + power2(
     Pi)) + mstau22*(-6*(3 + logmH2Q2)*logmstau22Q2 + 6*logmstau12Q2*(-9 +
     logmH2Q2 + 2*logmstau22Q2) + 9*power2(logmstau12Q2) + 3*power2(
     logmstau22Q2) + 2*(42 + power2(Pi))))*power2(-1 + power2(sa)) + mH2*(-2*
     mstau12*mstau22*(87 + 6*logmC2Q2*(-3 + logmsntau2Q2) - 18*logmstau22Q2 + 6
     *invdmsntau2mu*mu2 + 3*power2(logmC2Q2) + 3*power2(logmstau22Q2) - 9*
     power2(invdmsntau2mu)*power2(mu2) - 84*power2(invdmstau2mu)*power2(mu2) +
     18*logmstau22Q2*power2(invdmstau2mu)*power2(mu2) + 18*logmu2Q2*power2(
     invdmstau2mu)*power2(mu2) - 3*power2(invdmstau2mu)*power2(logmstau22Q2)*
     power2(mu2) - 3*power2(invdmstau2mu)*power2(logmu2Q2)*power2(mu2) + power2
     (logmsntau2Q2)*(6 - 6*power2(invdmstau2mu)*power2(mu2)) + 6*logmsntau2Q2*(
     -7 + logmstau22Q2 + power2(invdmsntau2mu)*power2(mu2) + 6*power2(
     invdmstau2mu)*power2(mu2) - logmstau22Q2*power2(invdmstau2mu)*power2(mu2)
     - logmu2Q2*power2(invdmstau2mu)*power2(mu2)) + 2*power2(Pi) - 2*power2(
     invdmstau2mu)*power2(mu2)*power2(Pi)) + msntau2*(mstau12*(-12 - 30*
     invdmstau1mu*mstau22 + 108*invdmstau2mu*mstau22 - 36*invdmstau2mu*
     logmsntau2Q2*mstau22 + 6*invdmsntau2mu*(-5 + 4*logmsntau2Q2)*mstau22 + 24*
     invdmstau1mu*logmstau12Q2*mstau22 - 36*invdmstau2mu*logmstau12Q2*mstau22 -
     24*invdmstau*logmA2Q2*logmstau12Q2*mstau22 + 24*invdmstau*logmH2Q2*
     logmstau12Q2*mstau22 - 24*invdmstau2mu*mu2 + 24*invdmstau*invdmstau1mu*
     logmA2Q2*logmstau12Q2*mstau22*mu2 - 24*invdmstau*invdmstau1mu*logmH2Q2*
     logmstau12Q2*mstau22*mu2 - 24*invdmstau*invdmstau1mu*logmA2Q2*logmu2Q2*
     mstau22*mu2 + 24*invdmstau*invdmstau2mu*logmA2Q2*logmu2Q2*mstau22*mu2 + 24
     *invdmstau*invdmstau1mu*logmH2Q2*logmu2Q2*mstau22*mu2 - 24*invdmstau*
     invdmstau2mu*logmH2Q2*logmu2Q2*mstau22*mu2 + 18*mstau22*mu2*power2(
     invdmstau1mu) + 24*logmstau12Q2*mstau22*mu2*power2(invdmstau1mu) - 24*
     logmu2Q2*mstau22*mu2*power2(invdmstau1mu) + 204*mstau22*mu2*power2(
     invdmstau2mu) - 72*logmsntau2Q2*mstau22*mu2*power2(invdmstau2mu) - 72*
     logmstau12Q2*mstau22*mu2*power2(invdmstau2mu) - 48*logmu2Q2*mstau22*mu2*
     power2(invdmstau2mu) + 12*logmsntau2Q2*logmu2Q2*mstau22*mu2*power2(
     invdmstau2mu) + 12*logmstau12Q2*logmu2Q2*mstau22*mu2*power2(invdmstau2mu)
     + 6*invdmstau2mu*mstau22*power2(logmsntau2Q2) + 12*mstau22*mu2*power2(
     invdmstau2mu)*power2(logmsntau2Q2) + 6*invdmstau2mu*mstau22*power2(
     logmstau12Q2) + 12*mstau22*mu2*power2(invdmstau2mu)*power2(logmstau12Q2) +
     6*mstau22*mu2*power2(invdmstau1mu)*power2(logmu2Q2) + 12*mstau22*mu2*
     power2(invdmstau2mu)*power2(logmu2Q2) - 84*power2(invdmstau1mu)*power2(
     mstau22) + 36*logmstau12Q2*power2(invdmstau1mu)*power2(mstau22) - 6*power2
     (invdmstau1mu)*power2(logmstau12Q2)*power2(mstau22) + 84*power2(
     invdmstau1mu)*power2(mu2) - 36*logmu2Q2*power2(invdmstau1mu)*power2(mu2) +
     36*power2(invdmstau2mu)*power2(mu2) + 6*power2(invdmstau1mu)*power2(
     logmu2Q2)*power2(mu2) + 4*invdmstau2mu*mstau22*power2(Pi) + 2*mstau22*mu2*
     power2(invdmstau1mu)*power2(Pi) + 8*mstau22*mu2*power2(invdmstau2mu)*
     power2(Pi) - 2*power2(invdmstau1mu)*power2(mstau22)*power2(Pi) + 2*power2(
     invdmstau1mu)*power2(mu2)*power2(Pi) + 2*power2(invdmsntau2mu)*(-(power2(
     mstau22)*(42 - 18*logmsntau2Q2 + 3*power2(logmsntau2Q2) + power2(Pi))) +
     power2(mu2)*(42 - 18*logmu2Q2 + 3*power2(logmu2Q2) + power2(Pi)) + mstau22
     *mu2*(9 + 12*logmsntau2Q2 - 12*logmu2Q2 + 3*power2(logmu2Q2) + power2(Pi))
     ) - 81*invdmhH*mstau22*power2(sa) + 144*invdmhH*logmh2Q2*mstau22*power2(sa
     ) - 60*invdmhH*logmH2Q2*mstau22*power2(sa) + 60*invdmhH*logmstau12Q2*
     mstau22*power2(sa) - 72*invdmhH*logmh2Q2*logmstau12Q2*mstau22*power2(sa) +
     24*invdmstau*logmh2Q2*logmstau12Q2*mstau22*power2(sa) - 24*invdmstau*
     logmH2Q2*logmstau12Q2*mstau22*power2(sa) - 24*invdmstau*invdmstau1mu*
     logmh2Q2*logmstau12Q2*mstau22*mu2*power2(sa) + 24*invdmstau*invdmstau1mu*
     logmH2Q2*logmstau12Q2*mstau22*mu2*power2(sa) + 24*invdmstau*invdmstau1mu*
     logmh2Q2*logmu2Q2*mstau22*mu2*power2(sa) - 24*invdmstau*invdmstau2mu*
     logmh2Q2*logmu2Q2*mstau22*mu2*power2(sa) - 24*invdmstau*invdmstau1mu*
     logmH2Q2*logmu2Q2*mstau22*mu2*power2(sa) + 24*invdmstau*invdmstau2mu*
     logmH2Q2*logmu2Q2*mstau22*mu2*power2(sa) - 36*invdmhH*mstau22*power2(
     logmh2Q2)*power2(sa) + 36*invdmhH*mstau22*power2(logmH2Q2)*power2(sa) + 60
     *invdmhH*invdmstau*logmstau12Q2*power2(mstau22)*power2(sa) - 72*invdmhH*
     invdmstau*logmh2Q2*logmstau12Q2*power2(mstau22)*power2(sa) - 12*mstau22*(-
     2*logmsntau2Q2 - (-2 + logmu2Q2)*logmu2Q2 + power2(logmsntau2Q2))*power2(
     mu2)*power3(invdmsntau2mu) + 24*logmstau12Q2*mstau22*power2(mu2)*power3(
     invdmstau1mu) - 24*logmu2Q2*mstau22*power2(mu2)*power3(invdmstau1mu) - 12*
     mstau22*power2(logmstau12Q2)*power2(mu2)*power3(invdmstau1mu) + 12*mstau22
     *power2(logmu2Q2)*power2(mu2)*power3(invdmstau1mu) - 48*logmu2Q2*mstau22*
     power2(mu2)*power3(invdmstau2mu) + 24*mstau22*power2(logmu2Q2)*power2(mu2)
     *power3(invdmstau2mu) - 6*power2(logmstau22Q2)*(-2*invdmstau2mu*mstau22 -
     2*mstau22*mu2*power2(invdmstau2mu) + (power2(invdmsntau2mu) + power2(
     invdmstau1mu))*(-(mstau22*mu2) + power2(mstau22) - power2(mu2)) + 4*
     mstau22*power2(mu2)*power3(invdmstau2mu)) + 12*logmstau22Q2*(2 - 3*mstau22
     *mu2*power2(invdmsntau2mu) + logmu2Q2*mstau22*mu2*power2(invdmsntau2mu) -
     3*mstau22*mu2*power2(invdmstau1mu) + logmu2Q2*mstau22*mu2*power2(
     invdmstau1mu) + ((-2 + logmsntau2Q2 + logmstau12Q2)*mstau22 - 2*mu2)*mu2*
     power2(invdmstau2mu) + 3*power2(invdmsntau2mu)*power2(mstau22) -
     logmsntau2Q2*power2(invdmsntau2mu)*power2(mstau22) + 3*power2(invdmstau1mu
     )*power2(mstau22) - logmstau12Q2*power2(invdmstau1mu)*power2(mstau22) - 3*
     power2(invdmsntau2mu)*power2(mu2) + logmu2Q2*power2(invdmsntau2mu)*power2(
     mu2) - 3*power2(invdmstau1mu)*power2(mu2) + logmu2Q2*power2(invdmstau1mu)*
     power2(mu2) + invdmstau2mu*mstau22*(-2 + logmsntau2Q2 + logmstau12Q2 - 2*
     invdmstau*logmA2Q2*mu2 + 2*invdmstau*logmH2Q2*mu2 + 2*invdmstau*logmh2Q2*
     mu2*power2(sa) - 2*invdmstau*logmH2Q2*mu2*power2(sa)) + invdmstau*mstau22*
     (2*logmA2Q2 + 2*logmH2Q2*(-1 + power2(sa)) + (logmh2Q2*(-2 - 6*invdmhH*
     mstau22*(-1 + power2(sa))) + 5*invdmhH*mstau22*(-1 + power2(sa)))*power2(
     sa)) + 4*mstau22*power2(mu2)*power3(invdmstau2mu)) + 81*invdmhH*mstau22*
     power4(sa) - 144*invdmhH*logmh2Q2*mstau22*power4(sa) + 60*invdmhH*logmH2Q2
     *mstau22*power4(sa) - 60*invdmhH*logmstau12Q2*mstau22*power4(sa) + 72*
     invdmhH*logmh2Q2*logmstau12Q2*mstau22*power4(sa) + 36*invdmhH*mstau22*
     power2(logmh2Q2)*power4(sa) - 36*invdmhH*mstau22*power2(logmH2Q2)*power4(
     sa) - 60*invdmhH*invdmstau*logmstau12Q2*power2(mstau22)*power4(sa) + 72*
     invdmhH*invdmstau*logmh2Q2*logmstau12Q2*power2(mstau22)*power4(sa)) + 2*
     mstau22*(-129 + 6*logmA2Q2*(-3 + logmstau12Q2) + 60*logmstau12Q2 + 18*
     logmstau22Q2 - 6*logmstau12Q2*logmstau22Q2 - 6*invdmstau1mu*mu2 + 3*power2
     (logmA2Q2) - 9*power2(logmstau12Q2) - 3*power2(logmstau22Q2) + 9*power2(
     invdmstau1mu)*power2(mu2) - 6*logmstau12Q2*power2(invdmstau1mu)*power2(mu2
     ) + 84*power2(invdmstau2mu)*power2(mu2) - 36*logmstau12Q2*power2(
     invdmstau2mu)*power2(mu2) - 18*logmstau22Q2*power2(invdmstau2mu)*power2(
     mu2) + 6*logmstau12Q2*logmstau22Q2*power2(invdmstau2mu)*power2(mu2) - 18*
     logmu2Q2*power2(invdmstau2mu)*power2(mu2) + 6*logmstau12Q2*logmu2Q2*power2
     (invdmstau2mu)*power2(mu2) + 6*power2(invdmstau2mu)*power2(logmstau12Q2)*
     power2(mu2) + 3*power2(invdmstau2mu)*power2(logmstau22Q2)*power2(mu2) + 3*
     power2(invdmstau2mu)*power2(logmu2Q2)*power2(mu2) - 3*power2(Pi) + 2*
     power2(invdmstau2mu)*power2(mu2)*power2(Pi) + 252*power2(sa) - 108*
     logmstau12Q2*power2(sa) + 18*power2(logmstau12Q2)*power2(sa) + 6*power2(Pi
     )*power2(sa) - 18*logmH2Q2*(-3 + logmstau12Q2)*power2(-1 + power2(sa)) - 9
     *power2(logmH2Q2)*power2(-1 + power2(sa)) - 252*power4(sa) + 54*logmh2Q2*
     power4(sa) + 108*logmstau12Q2*power4(sa) - 18*logmh2Q2*logmstau12Q2*power4
     (sa) - 9*power2(logmh2Q2)*power4(sa) - 18*power2(logmstau12Q2)*power4(sa)
     - 6*power2(Pi)*power4(sa)))))))))/(48.*mA2*mC2*mH2*msntau2*mstau12*mstau22
     *power2(mh2)) + (Fin20(mA2,mH2,Q2)*(-1 + power2(sa))*(mH2*(-10 + invdmAH*
     mH2)*power2(mA2) - mA2*power2(mH2) - 4*power3(mA2) + 5*power3(mH2) - mH2*
     power2(invdmAH)*power4(mA2) + power2(invdmAH)*power5(mA2)))/(16.*power2(
     mA2)*power2(mH2)) - (Fin20(mH2,mh2,Q2)*(-1 + power2(sa))*power2(sa)*(-(mH2
     *power2(mh2)*(3 + invdmhH*mH2 + 2*power2(invdmhH)*power2(mH2))) + mh2*
     power2(mH2)*(-4 + 3*power2(invdmhH)*power2(mH2)) + 3*power3(mh2) + 4*
     power3(mH2) - power2(invdmhH)*power5(mH2)))/(16.*power2(mh2)*power2(mH2))
     + (-1548 + 72*logmH2Q2 - 72*logmsntau2Q2 + 96*logmH2Q2*logmsntau2Q2 - 72*
     logmstau12Q2 - 96*logmH2Q2*logmstau12Q2 + 72*logmsntau2Q2*logmstau12Q2 +
     168*logmstau22Q2 - 96*logmH2Q2*logmstau22Q2 + 24*logmsntau2Q2*logmstau22Q2
      + 216*logmstau12Q2*logmstau22Q2 - 168*invdmAH*mA2 + 72*invdmstau1mu*mA2 +
     72*invdmstau2mu*mA2 + 72*invdmAH*logmH2Q2*mA2 - 48*invdmstau1mu*
     logmstau12Q2*mA2 - 48*invdmstau2mu*logmstau22Q2*mA2 + 105*invdmAC*mC2 +
     105*invdmCH*mC2 + 144*invdmsntau2mu*mC2 + 144*invdmstau2mu*mC2 - 72*
     invdmCH*logmH2Q2*mC2 - 96*invdmsntau2mu*logmsntau2Q2*mC2 - 96*invdmstau2mu
     *logmstau22Q2*mC2 + (57*mC2)/mA2 + (1104*mA2)/mH2 - (216*logmH2Q2*mA2)/mH2
      - (135*mC2)/mH2 - (216*logmH2Q2*mC2)/mH2 + 72*invdmstau1mu*mH2 + 72*
     invdmstau2mu*mH2 - 48*invdmstau1mu*logmH2Q2*mH2 - 48*invdmstau2mu*logmH2Q2
     *mH2 - 48*invdmstau1mu*logmstau12Q2*mH2 - 48*invdmstau2mu*logmstau22Q2*mH2
      - (216*mH2)/mA2 + (48*logmH2Q2*mH2)/mA2 + (48*mC2)/msntau2 - (96*
     logmsntau2Q2*mC2)/msntau2 + 144*invdmstau2mu*msntau2 - 96*invdmstau2mu*
     logmsntau2Q2*msntau2 - 96*invdmstau2mu*logmstau22Q2*msntau2 - (720*msntau2
     )/mA2 + (480*logmsntau2Q2*msntau2)/mA2 - (96*msntau2)/mC2 + (96*
     logmsntau2Q2*msntau2)/mC2 + (432*msntau2)/mH2 - (144*logmH2Q2*msntau2)/mH2
      - (288*logmsntau2Q2*msntau2)/mH2 + (72*logmH2Q2*logmsntau2Q2*msntau2)/mH2
      + (24*mA2)/mstau12 - (48*logmstau12Q2*mA2)/mstau12 + (24*mH2)/mstau12 - (
     48*logmstau12Q2*mH2)/mstau12 + 432*invdmstau2mu*mstau12 - 288*invdmstau2mu
     *logmstau12Q2*mstau12 - 192*invdmstau2mu*logmstau22Q2*mstau12 + 48*
     invdmstau2mu*logmstau12Q2*logmstau22Q2*mstau12 - (720*mstau12)/mA2 + (480*
     logmstau12Q2*mstau12)/mA2 + (96*mstau12)/mC2 - (96*logmstau12Q2*mstau12)/
     mC2 + (432*mstau12)/mH2 - (144*logmH2Q2*mstau12)/mH2 - (288*logmstau12Q2*
     mstau12)/mH2 + (72*logmH2Q2*logmstau12Q2*mstau12)/mH2 + (24*mA2)/mstau22 -
     (48*logmstau22Q2*mA2)/mstau22 + (48*mC2)/mstau22 - (96*logmstau22Q2*mC2)/
     mstau22 + (24*mH2)/mstau22 - (48*logmstau22Q2*mH2)/mstau22 + (48*msntau2)/
     mstau22 - (96*logmstau22Q2*msntau2)/mstau22 + (48*mstau12)/mstau22 - (96*
     logmstau22Q2*mstau12)/mstau22 + 306*invdmsntau2mu*mstau22 + 210*invdmstau*
     mstau22 + 306*invdmstau1mu*mstau22 - 192*invdmsntau2mu*logmsntau2Q2*
     mstau22 - 72*invdmstau*logmstau12Q2*mstau22 - 192*invdmstau1mu*
     logmstau12Q2*mstau22 - 216*invdmsntau2mu*logmstau22Q2*mstau22 - 144*
     invdmstau*logmstau22Q2*mstau22 - 216*invdmstau1mu*logmstau22Q2*mstau22 +
     48*invdmsntau2mu*logmsntau2Q2*logmstau22Q2*mstau22 + 48*invdmstau*
     logmstau12Q2*logmstau22Q2*mstau22 + 48*invdmstau1mu*logmstau12Q2*
     logmstau22Q2*mstau22 - (720*mstau22)/mA2 + (480*logmstau22Q2*mstau22)/mA2
     - (288*mstau22)/mC2 + (192*logmstau22Q2*mstau22)/mC2 + (432*mstau22)/mH2 -
     (144*logmH2Q2*mstau22)/mH2 - (288*logmstau22Q2*mstau22)/mH2 + (72*logmH2Q2
     *logmstau22Q2*mstau22)/mH2 + (48*mstau22)/msntau2 - (96*logmsntau2Q2*
     mstau22)/msntau2 + (48*mstau22)/mstau12 - (96*logmstau12Q2*mstau22)/
     mstau12 + 126*invdmstau*invdmstau1mu*mstau24 + 210*invdmstau*invdmstau2mu*
     mstau24 - 72*invdmstau*invdmstau2mu*logmstau12Q2*mstau24 - 72*invdmstau*
     invdmstau1mu*logmstau22Q2*mstau24 - 144*invdmstau*invdmstau2mu*
     logmstau22Q2*mstau24 + 48*invdmstau*invdmstau2mu*logmstau12Q2*logmstau22Q2
     *mstau24 + (126*invdmsntau2mu*mstau24)/(msntau2 - mstau22) - (72*
     invdmsntau2mu*logmstau22Q2*mstau24)/(msntau2 - mstau22) + 330*
     invdmsntau2mu*mu2 - 168*invdmstau*mu2 - 246*invdmstau1mu*mu2 - 552*
     invdmstau2mu*mu2 - 96*invdmsntau2mu*logmH2Q2*mu2 - 96*invdmstau1mu*
     logmH2Q2*mu2 - 96*invdmstau2mu*logmH2Q2*mu2 - 168*invdmsntau2mu*
     logmsntau2Q2*mu2 - 72*invdmstau1mu*logmsntau2Q2*mu2 + 264*invdmstau2mu*
     logmsntau2Q2*mu2 - 72*invdmsntau2mu*logmstau12Q2*mu2 - 168*invdmstau1mu*
     logmstau12Q2*mu2 + 264*invdmstau2mu*logmstau12Q2*mu2 + 192*invdmstau1mu*
     logmH2Q2*logmstau12Q2*mu2 - 24*invdmsntau2mu*logmstau22Q2*mu2 + 288*
     invdmstau*logmstau22Q2*mu2 + 360*invdmstau1mu*logmstau22Q2*mu2 + 312*
     invdmstau2mu*logmstau22Q2*mu2 + 192*invdmstau2mu*logmH2Q2*logmstau22Q2*mu2
      - 96*invdmstau2mu*logmsntau2Q2*logmstau22Q2*mu2 - 96*invdmstau*
     logmstau12Q2*logmstau22Q2*mu2 - 192*invdmstau1mu*logmstau12Q2*logmstau22Q2
     *mu2 - 288*invdmstau2mu*logmstau12Q2*logmstau22Q2*mu2 - 24*invdmsntau2mu*
     logmu2Q2*mu2 + 360*invdmstau1mu*logmu2Q2*mu2 - 288*invdmstau2mu*logmu2Q2*
     mu2 - 192*invdmstau1mu*logmH2Q2*logmu2Q2*mu2 - 192*invdmstau2mu*logmH2Q2*
     logmu2Q2*mu2 + 48*invdmsntau2mu*logmsntau2Q2*logmu2Q2*mu2 + 48*
     invdmstau1mu*logmstau12Q2*logmu2Q2*mu2 + 192*invdmstau2mu*logmstau12Q2*
     logmu2Q2*mu2 + 144*invdmstau2mu*logmstau22Q2*logmu2Q2*mu2 - (2880*mu2)/mA2
      + (720*logmsntau2Q2*mu2)/mA2 + (720*logmstau12Q2*mu2)/mA2 + (720*
     logmstau22Q2*mu2)/mA2 + (720*logmu2Q2*mu2)/mA2 - (240*logmsntau2Q2*
     logmu2Q2*mu2)/mA2 - (240*logmstau12Q2*logmu2Q2*mu2)/mA2 - (240*
     logmstau22Q2*logmu2Q2*mu2)/mA2 - (384*mu2)/mC2 + (288*logmstau22Q2*mu2)/
     mC2 + (96*logmu2Q2*mu2)/mC2 - (96*logmstau22Q2*logmu2Q2*mu2)/mC2 + (1728*
     mu2)/mH2 + (432*logmH2Q2*mu2)/mH2 - (432*logmsntau2Q2*mu2)/mH2 - (432*
     logmstau12Q2*mu2)/mH2 - (432*logmstau22Q2*mu2)/mH2 - (432*logmu2Q2*mu2)/
     mH2 - (216*logmH2Q2*logmu2Q2*mu2)/mH2 + (144*logmsntau2Q2*logmu2Q2*mu2)/
     mH2 + (144*logmstau12Q2*logmu2Q2*mu2)/mH2 + (144*logmstau22Q2*logmu2Q2*mu2
     )/mH2 - (720*mu2)/msntau2 + (144*logmH2Q2*mu2)/msntau2 + (672*logmsntau2Q2
     *mu2)/msntau2 - (96*logmH2Q2*logmsntau2Q2*mu2)/msntau2 - (288*logmstau22Q2
     *mu2)/msntau2 + (96*logmsntau2Q2*logmstau22Q2*mu2)/msntau2 + (48*logmH2Q2*
     logmu2Q2*mu2)/msntau2 - (96*logmsntau2Q2*logmu2Q2*mu2)/msntau2 + (96*
     invdmsntau2mu*mC2*mu2)/msntau2 - (720*mu2)/mstau12 + (144*logmH2Q2*mu2)/
     mstau12 + (672*logmstau12Q2*mu2)/mstau12 - (96*logmH2Q2*logmstau12Q2*mu2)/
     mstau12 - (288*logmstau22Q2*mu2)/mstau12 + (96*logmstau12Q2*logmstau22Q2*
     mu2)/mstau12 + (48*logmH2Q2*logmu2Q2*mu2)/mstau12 - (96*logmstau12Q2*
     logmu2Q2*mu2)/mstau12 + (48*invdmstau1mu*mA2*mu2)/mstau12 - (672*
     invdmstau2mu*mh2*mu2)/mstau12 + (144*invdmstau2mu*logmh2Q2*mh2*mu2)/
     mstau12 + (432*invdmstau2mu*logmstau12Q2*mh2*mu2)/mstau12 - (96*
     invdmstau2mu*logmh2Q2*logmstau12Q2*mh2*mu2)/mstau12 + (48*invdmstau2mu*
     logmh2Q2*logmstau22Q2*mh2*mu2)/mstau12 - (48*invdmstau2mu*logmstau12Q2*
     logmstau22Q2*mh2*mu2)/mstau12 + (48*invdmstau1mu*mH2*mu2)/mstau12 - (2784*
     mu2)/mstau22 + (144*logmH2Q2*mu2)/mstau22 + (1920*logmstau22Q2*mu2)/
     mstau22 - (96*logmH2Q2*logmstau22Q2*mu2)/mstau22 + (48*logmH2Q2*logmu2Q2*
     mu2)/mstau22 - (192*logmstau22Q2*logmu2Q2*mu2)/mstau22 + (48*invdmstau2mu*
     mA2*mu2)/mstau22 + (96*invdmstau2mu*mC2*mu2)/mstau22 + (48*invdmstau2mu*
     mH2*mu2)/mstau22 + (96*invdmstau2mu*msntau2*mu2)/mstau22 + (96*
     invdmstau2mu*mstau12*mu2)/mstau22 - 588*invdmstau*invdmstau2mu*mstau22*mu2
      - 576*invdmstau1mu*invdmstau2mu*mstau22*mu2 + 144*invdmstau*invdmstau2mu*
     logmstau12Q2*mstau22*mu2 + 576*invdmstau*invdmstau2mu*logmstau22Q2*mstau22
     *mu2 + 384*invdmstau1mu*invdmstau2mu*logmstau22Q2*mstau22*mu2 - 192*
     invdmstau*invdmstau2mu*logmstau12Q2*logmstau22Q2*mstau22*mu2 + 384*
     invdmstau1mu*invdmstau2mu*logmu2Q2*mstau22*mu2 - 192*invdmstau1mu*
     invdmstau2mu*logmstau22Q2*logmu2Q2*mstau22*mu2 + (96*invdmsntau2mu*mstau22
     *mu2)/msntau2 + (96*invdmstau1mu*mstau22*mu2)/mstau12 - (1344*invdmstau2mu
     *mu2*mw2)/msntau2 + (864*invdmstau2mu*logmsntau2Q2*mu2*mw2)/msntau2 - (96*
     invdmstau2mu*logmsntau2Q2*logmstau22Q2*mu2*mw2)/msntau2 + (288*
     invdmstau2mu*logmw2Q2*mu2*mw2)/msntau2 - (192*invdmstau2mu*logmsntau2Q2*
     logmw2Q2*mu2*mw2)/msntau2 + (96*invdmstau2mu*logmstau22Q2*logmw2Q2*mu2*mw2
     )/msntau2 - (672*invdmstau2mu*mu2*mz2)/mstau12 + (432*invdmstau2mu*
     logmstau12Q2*mu2*mz2)/mstau12 - (48*invdmstau2mu*logmstau12Q2*logmstau22Q2
     *mu2*mz2)/mstau12 + (144*invdmstau2mu*logmz2Q2*mu2*mz2)/mstau12 - (96*
     invdmstau2mu*logmstau12Q2*logmz2Q2*mu2*mz2)/mstau12 + (48*invdmstau2mu*
     logmstau22Q2*logmz2Q2*mu2*mz2)/mstau12 + (192*ca*invdmstau*logmh2Q2*
     logmstau12Q2*mu2*sa*sb)/cb - (192*ca*invdmstau*logmH2Q2*logmstau12Q2*mu2*
     sa*sb)/cb - (192*ca*invdmstau*logmh2Q2*logmstau22Q2*mu2*sa*sb)/cb + (192*
     ca*invdmstau*logmH2Q2*logmstau22Q2*mu2*sa*sb)/cb + 336*mC2*mu2*power2(
     invdmsntau2mu) - 96*logmsntau2Q2*mC2*mu2*power2(invdmsntau2mu) - 144*
     msntau2*mu2*power2(invdmsntau2mu) + 96*logmsntau2Q2*msntau2*mu2*power2(
     invdmsntau2mu) + 642*mstau22*mu2*power2(invdmsntau2mu) - 96*logmsntau2Q2*
     mstau22*mu2*power2(invdmsntau2mu) - 312*logmstau22Q2*mstau22*mu2*power2(
     invdmsntau2mu) - 48*logmsntau2Q2*logmstau22Q2*mstau22*mu2*power2(
     invdmsntau2mu) - 96*logmu2Q2*mstau22*mu2*power2(invdmsntau2mu) + 96*
     logmstau22Q2*logmu2Q2*mstau22*mu2*power2(invdmsntau2mu) - (126*mstau24*mu2
     *power2(invdmsntau2mu))/(msntau2 - mstau22) + (72*logmstau22Q2*mstau24*mu2
     *power2(invdmsntau2mu))/(msntau2 - mstau22) + 672*mu2*mw2*power2(
     invdmsntau2mu) - 288*logmsntau2Q2*mu2*mw2*power2(invdmsntau2mu) - 288*
     logmw2Q2*mu2*mw2*power2(invdmsntau2mu) + 96*logmsntau2Q2*logmw2Q2*mu2*mw2*
     power2(invdmsntau2mu) + 210*mstau24*power2(invdmstau) - 72*logmstau12Q2*
     mstau24*power2(invdmstau) - 144*logmstau22Q2*mstau24*power2(invdmstau) +
     48*logmstau12Q2*logmstau22Q2*mstau24*power2(invdmstau) - 168*mstau22*mu2*
     power2(invdmstau) + 288*logmstau22Q2*mstau22*mu2*power2(invdmstau) - 96*
     logmstau12Q2*logmstau22Q2*mstau22*mu2*power2(invdmstau) - 588*invdmstau2mu
     *mstau24*mu2*power2(invdmstau) + 144*invdmstau2mu*logmstau12Q2*mstau24*mu2
     *power2(invdmstau) + 576*invdmstau2mu*logmstau22Q2*mstau24*mu2*power2(
     invdmstau) - 192*invdmstau2mu*logmstau12Q2*logmstau22Q2*mstau24*mu2*power2
     (invdmstau) + 168*mA2*mu2*power2(invdmstau1mu) - 48*logmstau12Q2*mA2*mu2*
     power2(invdmstau1mu) + 336*mh2*mu2*power2(invdmstau1mu) - 144*logmh2Q2*mh2
     *mu2*power2(invdmstau1mu) - 144*logmstau12Q2*mh2*mu2*power2(invdmstau1mu)
     + 48*logmh2Q2*logmstau12Q2*mh2*mu2*power2(invdmstau1mu) + 168*mH2*mu2*
     power2(invdmstau1mu) - 96*logmH2Q2*mH2*mu2*power2(invdmstau1mu) - 48*
     logmstau12Q2*mH2*mu2*power2(invdmstau1mu) - 144*mstau12*mu2*power2(
     invdmstau1mu) + 96*logmstau12Q2*mstau12*mu2*power2(invdmstau1mu) + 642*
     mstau22*mu2*power2(invdmstau1mu) - 96*logmstau12Q2*mstau22*mu2*power2(
     invdmstau1mu) - 312*logmstau22Q2*mstau22*mu2*power2(invdmstau1mu) - 48*
     logmstau12Q2*logmstau22Q2*mstau22*mu2*power2(invdmstau1mu) - 96*logmu2Q2*
     mstau22*mu2*power2(invdmstau1mu) + 96*logmstau22Q2*logmu2Q2*mstau22*mu2*
     power2(invdmstau1mu) - 126*invdmstau*mstau24*mu2*power2(invdmstau1mu) + 72
     *invdmstau*logmstau22Q2*mstau24*mu2*power2(invdmstau1mu) + 336*mu2*mz2*
     power2(invdmstau1mu) - 144*logmstau12Q2*mu2*mz2*power2(invdmstau1mu) - 144
     *logmz2Q2*mu2*mz2*power2(invdmstau1mu) + 48*logmstau12Q2*logmz2Q2*mu2*mz2*
     power2(invdmstau1mu) + 168*mA2*mu2*power2(invdmstau2mu) - 48*logmstau22Q2*
     mA2*mu2*power2(invdmstau2mu) + 336*mC2*mu2*power2(invdmstau2mu) - 96*
     logmstau22Q2*mC2*mu2*power2(invdmstau2mu) + 168*mH2*mu2*power2(
     invdmstau2mu) - 96*logmH2Q2*mH2*mu2*power2(invdmstau2mu) - 48*logmstau22Q2
     *mH2*mu2*power2(invdmstau2mu) + 336*msntau2*mu2*power2(invdmstau2mu) - 192
     *logmsntau2Q2*msntau2*mu2*power2(invdmstau2mu) - 96*logmstau22Q2*msntau2*
     mu2*power2(invdmstau2mu) - 240*mstau12*mu2*power2(invdmstau2mu) + 192*
     logmstau12Q2*mstau12*mu2*power2(invdmstau2mu) + 192*logmstau22Q2*mstau12*
     mu2*power2(invdmstau2mu) - 144*logmstau12Q2*logmstau22Q2*mstau12*mu2*
     power2(invdmstau2mu) - 96*logmu2Q2*mstau12*mu2*power2(invdmstau2mu) + 48*
     logmstau12Q2*logmu2Q2*mstau12*mu2*power2(invdmstau2mu) - 288*mstau22*mu2*
     power2(invdmstau2mu) + 192*logmstau22Q2*mstau22*mu2*power2(invdmstau2mu) -
     342*invdmsntau2mu*mstau24*mu2*power2(invdmstau2mu) - 210*invdmstau*mstau24
     *mu2*power2(invdmstau2mu) - 342*invdmstau1mu*mstau24*mu2*power2(
     invdmstau2mu) + 72*invdmstau*logmstau12Q2*mstau24*mu2*power2(invdmstau2mu)
     + 144*invdmsntau2mu*logmstau22Q2*mstau24*mu2*power2(invdmstau2mu) + 144*
     invdmstau*logmstau22Q2*mstau24*mu2*power2(invdmstau2mu) + 144*invdmstau1mu
     *logmstau22Q2*mstau24*mu2*power2(invdmstau2mu) - 48*invdmstau*logmstau12Q2
     *logmstau22Q2*mstau24*mu2*power2(invdmstau2mu) + 216*invdmsntau2mu*
     logmu2Q2*mstau24*mu2*power2(invdmstau2mu) + 216*invdmstau1mu*logmu2Q2*
     mstau24*mu2*power2(invdmstau2mu) - (24*invdmstau2mu*mh2*mu2*power2(
     logmh2Q2))/mstau12 + 24*mh2*mu2*power2(invdmstau1mu)*power2(logmh2Q2) +
     144*power2(logmH2Q2) - 12*invdmAH*mA2*power2(logmH2Q2) + 12*invdmCH*mC2*
     power2(logmH2Q2) + (48*mA2*power2(logmH2Q2))/mH2 + (48*mC2*power2(logmH2Q2
     ))/mH2 + (84*mH2*power2(logmH2Q2))/mA2 + (36*msntau2*power2(logmH2Q2))/mH2
      + (36*mstau12*power2(logmH2Q2))/mH2 + (36*mstau22*power2(logmH2Q2))/mH2 -
     (108*mu2*power2(logmH2Q2))/mH2 - (24*mu2*power2(logmH2Q2))/msntau2 - (24*
     mu2*power2(logmH2Q2))/mstau12 - (24*mu2*power2(logmH2Q2))/mstau22 + 36*
     power2(logmsntau2Q2) - (60*msntau2*power2(logmsntau2Q2))/mA2 + (36*msntau2
     *power2(logmsntau2Q2))/mH2 + 24*invdmsntau2mu*mstau22*power2(logmsntau2Q2)
     + 24*invdmsntau2mu*mu2*power2(logmsntau2Q2) - 48*invdmstau2mu*mu2*power2(
     logmsntau2Q2) - (120*mu2*power2(logmsntau2Q2))/mA2 + (72*mu2*power2(
     logmsntau2Q2))/mH2 - (96*mu2*power2(logmsntau2Q2))/msntau2 - (144*
     invdmstau2mu*mu2*mw2*power2(logmsntau2Q2))/msntau2 - 24*mstau22*mu2*power2
     (invdmsntau2mu)*power2(logmsntau2Q2) + 48*mu2*mw2*power2(invdmsntau2mu)*
     power2(logmsntau2Q2) + 36*power2(logmstau12Q2) + 24*invdmstau2mu*mstau12*
     power2(logmstau12Q2) - (60*mstau12*power2(logmstau12Q2))/mA2 + (36*mstau12
     *power2(logmstau12Q2))/mH2 + 24*invdmstau*mstau22*power2(logmstau12Q2) +
     24*invdmstau1mu*mstau22*power2(logmstau12Q2) + 24*invdmstau*invdmstau2mu*
     mstau24*power2(logmstau12Q2) - 48*invdmstau*mu2*power2(logmstau12Q2) + 24*
     invdmstau1mu*mu2*power2(logmstau12Q2) - 48*invdmstau2mu*mu2*power2(
     logmstau12Q2) - (120*mu2*power2(logmstau12Q2))/mA2 + (72*mu2*power2(
     logmstau12Q2))/mH2 - (96*mu2*power2(logmstau12Q2))/mstau12 - (72*
     invdmstau2mu*mh2*mu2*power2(logmstau12Q2))/mstau12 - 96*invdmstau*
     invdmstau2mu*mstau22*mu2*power2(logmstau12Q2) - (72*invdmstau2mu*mu2*mz2*
     power2(logmstau12Q2))/mstau12 + 24*mstau24*power2(invdmstau)*power2(
     logmstau12Q2) - 48*mstau22*mu2*power2(invdmstau)*power2(logmstau12Q2) - 96
     *invdmstau2mu*mstau24*mu2*power2(invdmstau)*power2(logmstau12Q2) + 24*mh2*
     mu2*power2(invdmstau1mu)*power2(logmstau12Q2) - 24*mstau22*mu2*power2(
     invdmstau1mu)*power2(logmstau12Q2) + 24*mu2*mz2*power2(invdmstau1mu)*
     power2(logmstau12Q2) - 48*mstau12*mu2*power2(invdmstau2mu)*power2(
     logmstau12Q2) - 24*invdmstau*mstau24*mu2*power2(invdmstau2mu)*power2(
     logmstau12Q2) - 12*power2(logmstau22Q2) + 24*invdmstau2mu*mstau12*power2(
     logmstau22Q2) + 24*invdmsntau2mu*mstau22*power2(logmstau22Q2) + 24*
     invdmstau*mstau22*power2(logmstau22Q2) + 24*invdmstau1mu*mstau22*power2(
     logmstau22Q2) - (60*mstau22*power2(logmstau22Q2))/mA2 - (24*mstau22*power2
     (logmstau22Q2))/mC2 + (36*mstau22*power2(logmstau22Q2))/mH2 + 24*invdmstau
     *invdmstau2mu*mstau24*power2(logmstau22Q2) - 48*invdmstau*mu2*power2(
     logmstau22Q2) - 96*invdmstau1mu*mu2*power2(logmstau22Q2) - 24*invdmstau2mu
     *mu2*power2(logmstau22Q2) - (120*mu2*power2(logmstau22Q2))/mA2 - (48*mu2*
     power2(logmstau22Q2))/mC2 + (72*mu2*power2(logmstau22Q2))/mH2 + (48*mu2*
     power2(logmstau22Q2))/msntau2 + (48*mu2*power2(logmstau22Q2))/mstau12 - (
     288*mu2*power2(logmstau22Q2))/mstau22 - 96*invdmstau*invdmstau2mu*mstau22*
     mu2*power2(logmstau22Q2) - 96*invdmstau1mu*invdmstau2mu*mstau22*mu2*power2
     (logmstau22Q2) + 24*mstau22*mu2*power2(invdmsntau2mu)*power2(logmstau22Q2)
     + 24*mstau24*power2(invdmstau)*power2(logmstau22Q2) - 48*mstau22*mu2*
     power2(invdmstau)*power2(logmstau22Q2) - 96*invdmstau2mu*mstau24*mu2*
     power2(invdmstau)*power2(logmstau22Q2) + 24*mstau22*mu2*power2(
     invdmstau1mu)*power2(logmstau22Q2) - 72*mstau12*mu2*power2(invdmstau2mu)*
     power2(logmstau22Q2) - 24*invdmstau*mstau24*mu2*power2(invdmstau2mu)*
     power2(logmstau22Q2) + 24*invdmsntau2mu*mu2*power2(logmu2Q2) - 72*
     invdmstau1mu*mu2*power2(logmu2Q2) + 72*invdmstau2mu*mu2*power2(logmu2Q2) -
     (180*mu2*power2(logmu2Q2))/mA2 - (24*mu2*power2(logmu2Q2))/mC2 + (108*mu2*
     power2(logmu2Q2))/mH2 - 96*invdmstau1mu*invdmstau2mu*mstau22*mu2*power2(
     logmu2Q2) + 48*mstau22*mu2*power2(invdmsntau2mu)*power2(logmu2Q2) + 48*
     mstau22*mu2*power2(invdmstau1mu)*power2(logmu2Q2) + 24*mstau12*mu2*power2(
     invdmstau2mu)*power2(logmu2Q2) - (48*invdmstau2mu*mu2*mw2*power2(logmw2Q2)
     )/msntau2 + 48*mu2*mw2*power2(invdmsntau2mu)*power2(logmw2Q2) - (24*
     invdmstau2mu*mu2*mz2*power2(logmz2Q2))/mstau12 + 24*mu2*mz2*power2(
     invdmstau1mu)*power2(logmz2Q2) + 105*power2(invdmAH)*power2(mA2) - 72*
     logmH2Q2*power2(invdmAH)*power2(mA2) + 12*power2(invdmAH)*power2(logmH2Q2)
     *power2(mA2) + 42*power2(invdmAC)*power2(mC2) + 42*power2(invdmCH)*power2(
     mC2) - 72*logmH2Q2*power2(invdmCH)*power2(mC2) + 12*power2(invdmCH)*power2
     (logmH2Q2)*power2(mC2) - (357*power2(mC2))/power2(mA2) + (36*logmC2Q2*
     power2(mC2)*(9 + power2(invdmAC)*power2(mC2)))/power2(mA2) - (12*logmC2Q2*
     mC2*(11 + 3*power2(invdmAC)*power2(mC2)))/mA2 + (315*power2(mA2))/power2(
     mH2) + (315*power2(mC2))/power2(mH2) - (420*power2(mH2))/power2(mA2) + (
     360*logmH2Q2*power2(mH2))/power2(mA2) - (120*power2(logmH2Q2)*power2(mH2))
     /power2(mA2) + 72*invdmsntau2mu*invdmstau1mu*power2(mu2) + 24*
     invdmsntau2mu*invdmstau2mu*power2(mu2) + 252*invdmstau*invdmstau2mu*power2
     (mu2) - 2088*invdmstau1mu*invdmstau2mu*power2(mu2) - 432*invdmstau*
     invdmstau2mu*logmstau22Q2*power2(mu2) - 384*invdmstau1mu*invdmstau2mu*
     logmstau22Q2*power2(mu2) + 144*invdmstau*invdmstau2mu*logmstau12Q2*
     logmstau22Q2*power2(mu2) + 192*invdmstau1mu*invdmstau2mu*logmstau12Q2*
     logmstau22Q2*power2(mu2) + 1920*invdmstau1mu*invdmstau2mu*logmu2Q2*power2(
     mu2) - 192*invdmstau1mu*invdmstau2mu*logmstau12Q2*logmu2Q2*power2(mu2) + (
     1344*invdmsntau2mu*power2(mu2))/mA2 + (1344*invdmstau1mu*power2(mu2))/mA2
     + (1344*invdmstau2mu*power2(mu2))/mA2 - (576*invdmsntau2mu*logmsntau2Q2*
     power2(mu2))/mA2 - (576*invdmstau1mu*logmstau12Q2*power2(mu2))/mA2 - (576*
     invdmstau2mu*logmstau22Q2*power2(mu2))/mA2 - (576*invdmsntau2mu*logmu2Q2*
     power2(mu2))/mA2 - (576*invdmstau1mu*logmu2Q2*power2(mu2))/mA2 - (576*
     invdmstau2mu*logmu2Q2*power2(mu2))/mA2 + (192*invdmsntau2mu*logmsntau2Q2*
     logmu2Q2*power2(mu2))/mA2 + (192*invdmstau1mu*logmstau12Q2*logmu2Q2*power2
     (mu2))/mA2 + (192*invdmstau2mu*logmstau22Q2*logmu2Q2*power2(mu2))/mA2 - (
     1344*invdmsntau2mu*power2(mu2))/mH2 - (1344*invdmstau1mu*power2(mu2))/mH2
     - (1344*invdmstau2mu*power2(mu2))/mH2 + (576*invdmsntau2mu*logmsntau2Q2*
     power2(mu2))/mH2 + (576*invdmstau1mu*logmstau12Q2*power2(mu2))/mH2 + (576*
     invdmstau2mu*logmstau22Q2*power2(mu2))/mH2 + (576*invdmsntau2mu*logmu2Q2*
     power2(mu2))/mH2 + (576*invdmstau1mu*logmu2Q2*power2(mu2))/mH2 + (576*
     invdmstau2mu*logmu2Q2*power2(mu2))/mH2 - (192*invdmsntau2mu*logmsntau2Q2*
     logmu2Q2*power2(mu2))/mH2 - (192*invdmstau1mu*logmstau12Q2*logmu2Q2*power2
     (mu2))/mH2 - (192*invdmstau2mu*logmstau22Q2*logmu2Q2*power2(mu2))/mH2 - (
     96*invdmsntau2mu*power2(mu2))/msntau2 - (336*power2(mu2))/(mA2*msntau2) +
     (144*logmsntau2Q2*power2(mu2))/(mA2*msntau2) + (144*logmu2Q2*power2(mu2))/
     (mA2*msntau2) - (48*logmsntau2Q2*logmu2Q2*power2(mu2))/(mA2*msntau2) - (
     336*power2(mu2))/(mH2*msntau2) + (144*logmsntau2Q2*power2(mu2))/(mH2*
     msntau2) + (144*logmu2Q2*power2(mu2))/(mH2*msntau2) - (48*logmsntau2Q2*
     logmu2Q2*power2(mu2))/(mH2*msntau2) - (96*invdmstau1mu*power2(mu2))/
     mstau12 - (336*power2(mu2))/(mA2*mstau12) + (144*logmstau12Q2*power2(mu2))
     /(mA2*mstau12) + (144*logmu2Q2*power2(mu2))/(mA2*mstau12) - (48*
     logmstau12Q2*logmu2Q2*power2(mu2))/(mA2*mstau12) - (336*power2(mu2))/(mH2*
     mstau12) + (144*logmstau12Q2*power2(mu2))/(mH2*mstau12) + (144*logmu2Q2*
     power2(mu2))/(mH2*mstau12) - (48*logmstau12Q2*logmu2Q2*power2(mu2))/(mH2*
     mstau12) - (192*invdmstau2mu*power2(mu2))/mstau22 - (336*power2(mu2))/(mA2
     *mstau22) + (144*logmstau22Q2*power2(mu2))/(mA2*mstau22) + (144*logmu2Q2*
     power2(mu2))/(mA2*mstau22) - (48*logmstau22Q2*logmu2Q2*power2(mu2))/(mA2*
     mstau22) - (672*power2(mu2))/(mC2*mstau22) + (288*logmstau22Q2*power2(mu2)
     )/(mC2*mstau22) + (288*logmu2Q2*power2(mu2))/(mC2*mstau22) - (96*
     logmstau22Q2*logmu2Q2*power2(mu2))/(mC2*mstau22) - (336*power2(mu2))/(mH2*
     mstau22) + (144*logmstau22Q2*power2(mu2))/(mH2*mstau22) + (144*logmu2Q2*
     power2(mu2))/(mH2*mstau22) - (48*logmstau22Q2*logmu2Q2*power2(mu2))/(mH2*
     mstau22) - (192*ca*invdmstau*invdmstau1mu*logmh2Q2*logmstau12Q2*sa*sb*
     power2(mu2))/cb + (192*ca*invdmstau*invdmstau1mu*logmH2Q2*logmstau12Q2*sa*
     sb*power2(mu2))/cb + (192*ca*invdmstau*invdmstau2mu*logmh2Q2*logmstau22Q2*
     sa*sb*power2(mu2))/cb - (192*ca*invdmstau*invdmstau2mu*logmH2Q2*
     logmstau22Q2*sa*sb*power2(mu2))/cb + (192*ca*invdmstau*invdmstau1mu*
     logmh2Q2*logmu2Q2*sa*sb*power2(mu2))/cb - (192*ca*invdmstau*invdmstau2mu*
     logmh2Q2*logmu2Q2*sa*sb*power2(mu2))/cb - (192*ca*invdmstau*invdmstau1mu*
     logmH2Q2*logmu2Q2*sa*sb*power2(mu2))/cb + (192*ca*invdmstau*invdmstau2mu*
     logmH2Q2*logmu2Q2*sa*sb*power2(mu2))/cb - 480*power2(invdmsntau2mu)*power2
     (mu2) - 312*logmsntau2Q2*power2(invdmsntau2mu)*power2(mu2) - 96*logmH2Q2*
     logmsntau2Q2*power2(invdmsntau2mu)*power2(mu2) - 72*logmsntau2Q2*
     logmstau12Q2*power2(invdmsntau2mu)*power2(mu2) + 288*logmstau22Q2*power2(
     invdmsntau2mu)*power2(mu2) - 24*logmsntau2Q2*logmstau22Q2*power2(
     invdmsntau2mu)*power2(mu2) + 648*logmu2Q2*power2(invdmsntau2mu)*power2(mu2
     ) + 96*logmH2Q2*logmu2Q2*power2(invdmsntau2mu)*power2(mu2) + 312*
     logmsntau2Q2*logmu2Q2*power2(invdmsntau2mu)*power2(mu2) + 72*logmstau12Q2*
     logmu2Q2*power2(invdmsntau2mu)*power2(mu2) - 72*logmstau22Q2*logmu2Q2*
     power2(invdmsntau2mu)*power2(mu2) - (144*mC2*power2(invdmsntau2mu)*power2(
     mu2))/msntau2 + (96*logmsntau2Q2*mC2*power2(invdmsntau2mu)*power2(mu2))/
     msntau2 - (144*mstau22*power2(invdmsntau2mu)*power2(mu2))/msntau2 + (96*
     logmsntau2Q2*mstau22*power2(invdmsntau2mu)*power2(mu2))/msntau2 - (1344*
     mw2*power2(invdmsntau2mu)*power2(mu2))/mstau22 + (864*logmstau22Q2*mw2*
     power2(invdmsntau2mu)*power2(mu2))/mstau22 - (96*logmstau22Q2*logmu2Q2*mw2
     *power2(invdmsntau2mu)*power2(mu2))/mstau22 + (288*logmw2Q2*mw2*power2(
     invdmsntau2mu)*power2(mu2))/mstau22 - (192*logmstau22Q2*logmw2Q2*mw2*
     power2(invdmsntau2mu)*power2(mu2))/mstau22 + (96*logmu2Q2*logmw2Q2*mw2*
     power2(invdmsntau2mu)*power2(mu2))/mstau22 + 252*invdmstau2mu*mstau22*
     power2(invdmstau)*power2(mu2) - 432*invdmstau2mu*logmstau22Q2*mstau22*
     power2(invdmstau)*power2(mu2) + 144*invdmstau2mu*logmstau12Q2*logmstau22Q2
     *mstau22*power2(invdmstau)*power2(mu2) - 480*power2(invdmstau1mu)*power2(
     mu2) - 312*logmstau12Q2*power2(invdmstau1mu)*power2(mu2) - 96*logmH2Q2*
     logmstau12Q2*power2(invdmstau1mu)*power2(mu2) - 72*logmsntau2Q2*
     logmstau12Q2*power2(invdmstau1mu)*power2(mu2) + 288*logmstau22Q2*power2(
     invdmstau1mu)*power2(mu2) - 24*logmstau12Q2*logmstau22Q2*power2(
     invdmstau1mu)*power2(mu2) + 648*logmu2Q2*power2(invdmstau1mu)*power2(mu2)
     + 96*logmH2Q2*logmu2Q2*power2(invdmstau1mu)*power2(mu2) + 72*logmsntau2Q2*
     logmu2Q2*power2(invdmstau1mu)*power2(mu2) + 312*logmstau12Q2*logmu2Q2*
     power2(invdmstau1mu)*power2(mu2) - 72*logmstau22Q2*logmu2Q2*power2(
     invdmstau1mu)*power2(mu2) - (72*mA2*power2(invdmstau1mu)*power2(mu2))/
     mstau12 + (48*logmstau12Q2*mA2*power2(invdmstau1mu)*power2(mu2))/mstau12 -
     (72*mH2*power2(invdmstau1mu)*power2(mu2))/mstau12 + (48*logmstau12Q2*mH2*
     power2(invdmstau1mu)*power2(mu2))/mstau12 - (672*mh2*power2(invdmstau1mu)*
     power2(mu2))/mstau22 + (144*logmh2Q2*mh2*power2(invdmstau1mu)*power2(mu2))
     /mstau22 + (432*logmstau22Q2*mh2*power2(invdmstau1mu)*power2(mu2))/mstau22
      - (96*logmh2Q2*logmstau22Q2*mh2*power2(invdmstau1mu)*power2(mu2))/mstau22
      + (48*logmh2Q2*logmu2Q2*mh2*power2(invdmstau1mu)*power2(mu2))/mstau22 - (
     48*logmstau22Q2*logmu2Q2*mh2*power2(invdmstau1mu)*power2(mu2))/mstau22 - (
     144*mstau22*power2(invdmstau1mu)*power2(mu2))/mstau12 + (96*logmstau12Q2*
     mstau22*power2(invdmstau1mu)*power2(mu2))/mstau12 - (672*mz2*power2(
     invdmstau1mu)*power2(mu2))/mstau22 + (432*logmstau22Q2*mz2*power2(
     invdmstau1mu)*power2(mu2))/mstau22 - (48*logmstau22Q2*logmu2Q2*mz2*power2(
     invdmstau1mu)*power2(mu2))/mstau22 + (144*logmz2Q2*mz2*power2(invdmstau1mu
     )*power2(mu2))/mstau22 - (96*logmstau22Q2*logmz2Q2*mz2*power2(invdmstau1mu
     )*power2(mu2))/mstau22 + (48*logmu2Q2*logmz2Q2*mz2*power2(invdmstau1mu)*
     power2(mu2))/mstau22 - 2556*power2(invdmstau2mu)*power2(mu2) + 576*
     logmsntau2Q2*power2(invdmstau2mu)*power2(mu2) + 576*logmstau12Q2*power2(
     invdmstau2mu)*power2(mu2) - 360*logmstau22Q2*power2(invdmstau2mu)*power2(
     mu2) - 96*logmH2Q2*logmstau22Q2*power2(invdmstau2mu)*power2(mu2) - 120*
     logmsntau2Q2*logmstau22Q2*power2(invdmstau2mu)*power2(mu2) - 120*
     logmstau12Q2*logmstau22Q2*power2(invdmstau2mu)*power2(mu2) + 1704*logmu2Q2
     *power2(invdmstau2mu)*power2(mu2) + 96*logmH2Q2*logmu2Q2*power2(
     invdmstau2mu)*power2(mu2) - 72*logmsntau2Q2*logmu2Q2*power2(invdmstau2mu)*
     power2(mu2) - 72*logmstau12Q2*logmu2Q2*power2(invdmstau2mu)*power2(mu2) +
     408*logmstau22Q2*logmu2Q2*power2(invdmstau2mu)*power2(mu2) - (1344*mh2*
     power2(invdmstau2mu)*power2(mu2))/mstau12 + (288*logmh2Q2*mh2*power2(
     invdmstau2mu)*power2(mu2))/mstau12 + (864*logmstau12Q2*mh2*power2(
     invdmstau2mu)*power2(mu2))/mstau12 - (192*logmh2Q2*logmstau12Q2*mh2*power2
     (invdmstau2mu)*power2(mu2))/mstau12 + (48*logmh2Q2*logmstau22Q2*mh2*power2
     (invdmstau2mu)*power2(mu2))/mstau12 - (48*logmstau12Q2*logmstau22Q2*mh2*
     power2(invdmstau2mu)*power2(mu2))/mstau12 + (48*logmh2Q2*logmu2Q2*mh2*
     power2(invdmstau2mu)*power2(mu2))/mstau12 - (48*logmstau12Q2*logmu2Q2*mh2*
     power2(invdmstau2mu)*power2(mu2))/mstau12 - (72*mA2*power2(invdmstau2mu)*
     power2(mu2))/mstau22 + (48*logmstau22Q2*mA2*power2(invdmstau2mu)*power2(
     mu2))/mstau22 - (144*mC2*power2(invdmstau2mu)*power2(mu2))/mstau22 + (96*
     logmstau22Q2*mC2*power2(invdmstau2mu)*power2(mu2))/mstau22 - (72*mH2*
     power2(invdmstau2mu)*power2(mu2))/mstau22 + (48*logmstau22Q2*mH2*power2(
     invdmstau2mu)*power2(mu2))/mstau22 - (144*msntau2*power2(invdmstau2mu)*
     power2(mu2))/mstau22 + (96*logmstau22Q2*msntau2*power2(invdmstau2mu)*
     power2(mu2))/mstau22 - (144*mstau12*power2(invdmstau2mu)*power2(mu2))/
     mstau22 + (96*logmstau22Q2*mstau12*power2(invdmstau2mu)*power2(mu2))/
     mstau22 + 6*invdmsntau2mu*mstau22*power2(invdmstau2mu)*power2(mu2) + 294*
     invdmstau*mstau22*power2(invdmstau2mu)*power2(mu2) + 6*invdmstau1mu*
     mstau22*power2(invdmstau2mu)*power2(mu2) - 72*invdmstau*logmstau12Q2*
     mstau22*power2(invdmstau2mu)*power2(mu2) - 288*invdmstau*logmstau22Q2*
     mstau22*power2(invdmstau2mu)*power2(mu2) + 96*invdmstau*logmstau12Q2*
     logmstau22Q2*mstau22*power2(invdmstau2mu)*power2(mu2) - 72*invdmsntau2mu*
     logmu2Q2*mstau22*power2(invdmstau2mu)*power2(mu2) - 72*invdmstau1mu*
     logmu2Q2*mstau22*power2(invdmstau2mu)*power2(mu2) - 48*invdmsntau2mu*
     logmstau22Q2*logmu2Q2*mstau22*power2(invdmstau2mu)*power2(mu2) - 48*
     invdmstau1mu*logmstau22Q2*logmu2Q2*mstau22*power2(invdmstau2mu)*power2(mu2
     ) - (2688*mw2*power2(invdmstau2mu)*power2(mu2))/msntau2 + (1728*
     logmsntau2Q2*mw2*power2(invdmstau2mu)*power2(mu2))/msntau2 - (96*
     logmsntau2Q2*logmstau22Q2*mw2*power2(invdmstau2mu)*power2(mu2))/msntau2 -
     (96*logmsntau2Q2*logmu2Q2*mw2*power2(invdmstau2mu)*power2(mu2))/msntau2 +
     (576*logmw2Q2*mw2*power2(invdmstau2mu)*power2(mu2))/msntau2 - (384*
     logmsntau2Q2*logmw2Q2*mw2*power2(invdmstau2mu)*power2(mu2))/msntau2 + (96*
     logmstau22Q2*logmw2Q2*mw2*power2(invdmstau2mu)*power2(mu2))/msntau2 + (96*
     logmu2Q2*logmw2Q2*mw2*power2(invdmstau2mu)*power2(mu2))/msntau2 - (1344*
     mz2*power2(invdmstau2mu)*power2(mu2))/mstau12 + (864*logmstau12Q2*mz2*
     power2(invdmstau2mu)*power2(mu2))/mstau12 - (48*logmstau12Q2*logmstau22Q2*
     mz2*power2(invdmstau2mu)*power2(mu2))/mstau12 - (48*logmstau12Q2*logmu2Q2*
     mz2*power2(invdmstau2mu)*power2(mu2))/mstau12 + (288*logmz2Q2*mz2*power2(
     invdmstau2mu)*power2(mu2))/mstau12 - (192*logmstau12Q2*logmz2Q2*mz2*power2
     (invdmstau2mu)*power2(mu2))/mstau12 + (48*logmstau22Q2*logmz2Q2*mz2*power2
     (invdmstau2mu)*power2(mu2))/mstau12 + (48*logmu2Q2*logmz2Q2*mz2*power2(
     invdmstau2mu)*power2(mu2))/mstau12 + 294*mstau24*power2(invdmsntau2mu)*
     power2(invdmstau2mu)*power2(mu2) - 504*logmstau22Q2*mstau24*power2(
     invdmsntau2mu)*power2(invdmstau2mu)*power2(mu2) + 144*logmu2Q2*mstau24*
     power2(invdmsntau2mu)*power2(invdmstau2mu)*power2(mu2) + 96*logmstau22Q2*
     logmu2Q2*mstau24*power2(invdmsntau2mu)*power2(invdmstau2mu)*power2(mu2) +
     294*mstau24*power2(invdmstau)*power2(invdmstau2mu)*power2(mu2) - 72*
     logmstau12Q2*mstau24*power2(invdmstau)*power2(invdmstau2mu)*power2(mu2) -
     288*logmstau22Q2*mstau24*power2(invdmstau)*power2(invdmstau2mu)*power2(mu2
     ) + 96*logmstau12Q2*logmstau22Q2*mstau24*power2(invdmstau)*power2(
     invdmstau2mu)*power2(mu2) + 294*mstau24*power2(invdmstau1mu)*power2(
     invdmstau2mu)*power2(mu2) - 504*logmstau22Q2*mstau24*power2(invdmstau1mu)*
     power2(invdmstau2mu)*power2(mu2) + 144*logmu2Q2*mstau24*power2(
     invdmstau1mu)*power2(invdmstau2mu)*power2(mu2) + 96*logmstau22Q2*logmu2Q2*
     mstau24*power2(invdmstau1mu)*power2(invdmstau2mu)*power2(mu2) - (24*mh2*
     power2(invdmstau1mu)*power2(logmh2Q2)*power2(mu2))/mstau22 - (48*mh2*
     power2(invdmstau2mu)*power2(logmh2Q2)*power2(mu2))/mstau12 + (96*
     invdmsntau2mu*power2(logmsntau2Q2)*power2(mu2))/mA2 - (96*invdmsntau2mu*
     power2(logmsntau2Q2)*power2(mu2))/mH2 - (24*power2(logmsntau2Q2)*power2(
     mu2))/(mA2*msntau2) - (24*power2(logmsntau2Q2)*power2(mu2))/(mH2*msntau2)
     - 48*power2(invdmsntau2mu)*power2(logmsntau2Q2)*power2(mu2) - 96*power2(
     invdmstau2mu)*power2(logmsntau2Q2)*power2(mu2) - (288*mw2*power2(
     invdmstau2mu)*power2(logmsntau2Q2)*power2(mu2))/msntau2 + 72*invdmstau*
     invdmstau2mu*power2(logmstau12Q2)*power2(mu2) + (96*invdmstau1mu*power2(
     logmstau12Q2)*power2(mu2))/mA2 - (96*invdmstau1mu*power2(logmstau12Q2)*
     power2(mu2))/mH2 - (24*power2(logmstau12Q2)*power2(mu2))/(mA2*mstau12) - (
     24*power2(logmstau12Q2)*power2(mu2))/(mH2*mstau12) + 72*invdmstau2mu*
     mstau22*power2(invdmstau)*power2(logmstau12Q2)*power2(mu2) - 48*power2(
     invdmstau1mu)*power2(logmstau12Q2)*power2(mu2) - 96*power2(invdmstau2mu)*
     power2(logmstau12Q2)*power2(mu2) - (144*mh2*power2(invdmstau2mu)*power2(
     logmstau12Q2)*power2(mu2))/mstau12 + 48*invdmstau*mstau22*power2(
     invdmstau2mu)*power2(logmstau12Q2)*power2(mu2) - (144*mz2*power2(
     invdmstau2mu)*power2(logmstau12Q2)*power2(mu2))/mstau12 + 48*mstau24*
     power2(invdmstau)*power2(invdmstau2mu)*power2(logmstau12Q2)*power2(mu2) +
     72*invdmstau*invdmstau2mu*power2(logmstau22Q2)*power2(mu2) + 96*
     invdmstau1mu*invdmstau2mu*power2(logmstau22Q2)*power2(mu2) + (96*
     invdmstau2mu*power2(logmstau22Q2)*power2(mu2))/mA2 - (96*invdmstau2mu*
     power2(logmstau22Q2)*power2(mu2))/mH2 - (24*power2(logmstau22Q2)*power2(
     mu2))/(mA2*mstau22) - (48*power2(logmstau22Q2)*power2(mu2))/(mC2*mstau22)
     - (24*power2(logmstau22Q2)*power2(mu2))/(mH2*mstau22) - 48*power2(
     invdmsntau2mu)*power2(logmstau22Q2)*power2(mu2) - (144*mw2*power2(
     invdmsntau2mu)*power2(logmstau22Q2)*power2(mu2))/mstau22 + 72*invdmstau2mu
     *mstau22*power2(invdmstau)*power2(logmstau22Q2)*power2(mu2) - 48*power2(
     invdmstau1mu)*power2(logmstau22Q2)*power2(mu2) - (72*mh2*power2(
     invdmstau1mu)*power2(logmstau22Q2)*power2(mu2))/mstau22 - (72*mz2*power2(
     invdmstau1mu)*power2(logmstau22Q2)*power2(mu2))/mstau22 - 48*power2(
     invdmstau2mu)*power2(logmstau22Q2)*power2(mu2) - 24*invdmsntau2mu*mstau22*
     power2(invdmstau2mu)*power2(logmstau22Q2)*power2(mu2) + 48*invdmstau*
     mstau22*power2(invdmstau2mu)*power2(logmstau22Q2)*power2(mu2) - 24*
     invdmstau1mu*mstau22*power2(invdmstau2mu)*power2(logmstau22Q2)*power2(mu2)
     + 48*mstau24*power2(invdmsntau2mu)*power2(invdmstau2mu)*power2(
     logmstau22Q2)*power2(mu2) + 48*mstau24*power2(invdmstau)*power2(
     invdmstau2mu)*power2(logmstau22Q2)*power2(mu2) + 48*mstau24*power2(
     invdmstau1mu)*power2(invdmstau2mu)*power2(logmstau22Q2)*power2(mu2) - 480*
     invdmstau1mu*invdmstau2mu*power2(logmu2Q2)*power2(mu2) + (96*invdmsntau2mu
     *power2(logmu2Q2)*power2(mu2))/mA2 + (96*invdmstau1mu*power2(logmu2Q2)*
     power2(mu2))/mA2 + (96*invdmstau2mu*power2(logmu2Q2)*power2(mu2))/mA2 - (
     96*invdmsntau2mu*power2(logmu2Q2)*power2(mu2))/mH2 - (96*invdmstau1mu*
     power2(logmu2Q2)*power2(mu2))/mH2 - (96*invdmstau2mu*power2(logmu2Q2)*
     power2(mu2))/mH2 - (24*power2(logmu2Q2)*power2(mu2))/(mA2*msntau2) - (24*
     power2(logmu2Q2)*power2(mu2))/(mH2*msntau2) - (24*power2(logmu2Q2)*power2(
     mu2))/(mA2*mstau12) - (24*power2(logmu2Q2)*power2(mu2))/(mH2*mstau12) - (
     24*power2(logmu2Q2)*power2(mu2))/(mA2*mstau22) - (48*power2(logmu2Q2)*
     power2(mu2))/(mC2*mstau22) - (24*power2(logmu2Q2)*power2(mu2))/(mH2*
     mstau22) - 312*power2(invdmsntau2mu)*power2(logmu2Q2)*power2(mu2) - 312*
     power2(invdmstau1mu)*power2(logmu2Q2)*power2(mu2) - 552*power2(
     invdmstau2mu)*power2(logmu2Q2)*power2(mu2) - 24*invdmsntau2mu*mstau22*
     power2(invdmstau2mu)*power2(logmu2Q2)*power2(mu2) - 24*invdmstau1mu*
     mstau22*power2(invdmstau2mu)*power2(logmu2Q2)*power2(mu2) + 48*mstau24*
     power2(invdmsntau2mu)*power2(invdmstau2mu)*power2(logmu2Q2)*power2(mu2) +
     48*mstau24*power2(invdmstau1mu)*power2(invdmstau2mu)*power2(logmu2Q2)*
     power2(mu2) - (48*mw2*power2(invdmsntau2mu)*power2(logmw2Q2)*power2(mu2))/
     mstau22 - (96*mw2*power2(invdmstau2mu)*power2(logmw2Q2)*power2(mu2))/
     msntau2 - (24*mz2*power2(invdmstau1mu)*power2(logmz2Q2)*power2(mu2))/
     mstau22 - (48*mz2*power2(invdmstau2mu)*power2(logmz2Q2)*power2(mu2))/
     mstau12 - 16*power2(Pi) - 4*invdmAH*mA2*power2(Pi) + 4*invdmAC*mC2*power2(
     Pi) + 4*invdmCH*mC2*power2(Pi) + (32*mC2*power2(Pi))/mA2 + (52*mA2*power2(
     Pi))/mH2 - (32*mC2*power2(Pi))/mH2 + (16*mH2*power2(Pi))/mA2 - (20*msntau2
     *power2(Pi))/mA2 + (12*msntau2*power2(Pi))/mH2 + 8*invdmstau2mu*mstau12*
     power2(Pi) - (20*mstau12*power2(Pi))/mA2 + (12*mstau12*power2(Pi))/mH2 + 8
     *invdmsntau2mu*mstau22*power2(Pi) + 8*invdmstau*mstau22*power2(Pi) + 8*
     invdmstau1mu*mstau22*power2(Pi) - (20*mstau22*power2(Pi))/mA2 - (8*mstau22
     *power2(Pi))/mC2 + (12*mstau22*power2(Pi))/mH2 + 8*invdmstau*invdmstau2mu*
     mstau24*power2(Pi) + 8*invdmsntau2mu*mu2*power2(Pi) - 16*invdmstau*mu2*
     power2(Pi) - 24*invdmstau1mu*mu2*power2(Pi) - 8*invdmstau2mu*mu2*power2(Pi
     ) - (60*mu2*power2(Pi))/mA2 - (8*mu2*power2(Pi))/mC2 + (36*mu2*power2(Pi))
     /mH2 - (16*mu2*power2(Pi))/msntau2 - (16*mu2*power2(Pi))/mstau12 - (16*
     invdmstau2mu*mh2*mu2*power2(Pi))/mstau12 - (64*mu2*power2(Pi))/mstau22 -
     32*invdmstau*invdmstau2mu*mstau22*mu2*power2(Pi) - 32*invdmstau1mu*
     invdmstau2mu*mstau22*mu2*power2(Pi) - (32*invdmstau2mu*mu2*mw2*power2(Pi))
     /msntau2 - (16*invdmstau2mu*mu2*mz2*power2(Pi))/mstau12 + 8*mstau22*mu2*
     power2(invdmsntau2mu)*power2(Pi) + 16*mu2*mw2*power2(invdmsntau2mu)*power2
     (Pi) + 8*mstau24*power2(invdmstau)*power2(Pi) - 16*mstau22*mu2*power2(
     invdmstau)*power2(Pi) - 32*invdmstau2mu*mstau24*mu2*power2(invdmstau)*
     power2(Pi) + 8*mh2*mu2*power2(invdmstau1mu)*power2(Pi) + 8*mstau22*mu2*
     power2(invdmstau1mu)*power2(Pi) + 8*mu2*mz2*power2(invdmstau1mu)*power2(Pi
     ) - 16*mstau12*mu2*power2(invdmstau2mu)*power2(Pi) - 8*invdmstau*mstau24*
     mu2*power2(invdmstau2mu)*power2(Pi) + 4*power2(invdmAH)*power2(mA2)*power2
     (Pi) + 4*power2(invdmAC)*power2(mC2)*power2(Pi) + 4*power2(invdmCH)*power2
     (mC2)*power2(Pi) - (30*power2(mC2)*power2(Pi))/power2(mA2) + (18*power2(
     mA2)*power2(Pi))/power2(mH2) + (18*power2(mC2)*power2(Pi))/power2(mH2) - (
     30*power2(mH2)*power2(Pi))/power2(mA2) + 24*invdmstau*invdmstau2mu*power2(
     mu2)*power2(Pi) - 32*invdmstau1mu*invdmstau2mu*power2(mu2)*power2(Pi) + (
     32*invdmsntau2mu*power2(mu2)*power2(Pi))/mA2 + (32*invdmstau1mu*power2(mu2
     )*power2(Pi))/mA2 + (32*invdmstau2mu*power2(mu2)*power2(Pi))/mA2 - (32*
     invdmsntau2mu*power2(mu2)*power2(Pi))/mH2 - (32*invdmstau1mu*power2(mu2)*
     power2(Pi))/mH2 - (32*invdmstau2mu*power2(mu2)*power2(Pi))/mH2 - (8*power2
     (mu2)*power2(Pi))/(mA2*msntau2) - (8*power2(mu2)*power2(Pi))/(mH2*msntau2)
     - (8*power2(mu2)*power2(Pi))/(mA2*mstau12) - (8*power2(mu2)*power2(Pi))/(
     mH2*mstau12) - (8*power2(mu2)*power2(Pi))/(mA2*mstau22) - (16*power2(mu2)*
     power2(Pi))/(mC2*mstau22) - (8*power2(mu2)*power2(Pi))/(mH2*mstau22) - 16*
     power2(invdmsntau2mu)*power2(mu2)*power2(Pi) - (32*mw2*power2(
     invdmsntau2mu)*power2(mu2)*power2(Pi))/mstau22 + 24*invdmstau2mu*mstau22*
     power2(invdmstau)*power2(mu2)*power2(Pi) - 16*power2(invdmstau1mu)*power2(
     mu2)*power2(Pi) - (16*mh2*power2(invdmstau1mu)*power2(mu2)*power2(Pi))/
     mstau22 - (16*mz2*power2(invdmstau1mu)*power2(mu2)*power2(Pi))/mstau22 -
     64*power2(invdmstau2mu)*power2(mu2)*power2(Pi) - (32*mh2*power2(
     invdmstau2mu)*power2(mu2)*power2(Pi))/mstau12 - 8*invdmsntau2mu*mstau22*
     power2(invdmstau2mu)*power2(mu2)*power2(Pi) + 16*invdmstau*mstau22*power2(
     invdmstau2mu)*power2(mu2)*power2(Pi) - 8*invdmstau1mu*mstau22*power2(
     invdmstau2mu)*power2(mu2)*power2(Pi) - (64*mw2*power2(invdmstau2mu)*power2
     (mu2)*power2(Pi))/msntau2 - (32*mz2*power2(invdmstau2mu)*power2(mu2)*
     power2(Pi))/mstau12 + 16*mstau24*power2(invdmsntau2mu)*power2(invdmstau2mu
     )*power2(mu2)*power2(Pi) + 16*mstau24*power2(invdmstau)*power2(
     invdmstau2mu)*power2(mu2)*power2(Pi) + 16*mstau24*power2(invdmstau1mu)*
     power2(invdmstau2mu)*power2(mu2)*power2(Pi) + 792*power2(sa) - 576*
     logmh2Q2*power2(sa) + 432*logmH2Q2*power2(sa) + 576*logmh2Q2*logmH2Q2*
     power2(sa) + 96*logmh2Q2*logmsntau2Q2*power2(sa) - 96*logmH2Q2*
     logmsntau2Q2*power2(sa) - 96*logmh2Q2*logmstau12Q2*power2(sa) + 96*
     logmH2Q2*logmstau12Q2*power2(sa) - 96*logmh2Q2*logmstau22Q2*power2(sa) +
     96*logmH2Q2*logmstau22Q2*power2(sa) - 168*invdmAh*mA2*power2(sa) + 168*
     invdmAH*mA2*power2(sa) + 72*invdmAh*logmh2Q2*mA2*power2(sa) - 72*invdmAH*
     logmH2Q2*mA2*power2(sa) + 105*invdmCh*mC2*power2(sa) - 105*invdmCH*mC2*
     power2(sa) - 72*invdmCh*logmh2Q2*mC2*power2(sa) + 72*invdmCH*logmH2Q2*mC2*
     power2(sa) + (1104*mA2*power2(sa))/mh2 - (216*logmh2Q2*mA2*power2(sa))/mh2
      - (135*mC2*power2(sa))/mh2 - (216*logmh2Q2*mC2*power2(sa))/mh2 + 72*
     invdmstau1mu*mh2*power2(sa) + 72*invdmstau2mu*mh2*power2(sa) - 48*
     invdmstau1mu*logmh2Q2*mh2*power2(sa) - 48*invdmstau2mu*logmh2Q2*mh2*power2
     (sa) - 48*invdmstau1mu*logmstau12Q2*mh2*power2(sa) - 48*invdmstau2mu*
     logmstau22Q2*mh2*power2(sa) - (216*mh2*power2(sa))/mA2 + (48*logmh2Q2*mh2*
     power2(sa))/mA2 - (1104*mA2*power2(sa))/mH2 + (216*logmH2Q2*mA2*power2(sa)
     )/mH2 + (135*mC2*power2(sa))/mH2 + (216*logmH2Q2*mC2*power2(sa))/mH2 - (72
     *mh2*power2(sa))/mH2 + (144*logmh2Q2*mh2*power2(sa))/mH2 - (144*logmH2Q2*
     mh2*power2(sa))/mH2 + (72*logmh2Q2*logmH2Q2*mh2*power2(sa))/mH2 - 1758*
     invdmhH*mH2*power2(sa) - 72*invdmstau1mu*mH2*power2(sa) - 72*invdmstau2mu*
     mH2*power2(sa) + 1080*invdmhH*logmh2Q2*mH2*power2(sa) + 384*invdmhH*
     logmH2Q2*mH2*power2(sa) + 48*invdmstau1mu*logmH2Q2*mH2*power2(sa) + 48*
     invdmstau2mu*logmH2Q2*mH2*power2(sa) + 24*invdmhH*logmh2Q2*logmH2Q2*mH2*
     power2(sa) + 48*invdmstau1mu*logmstau12Q2*mH2*power2(sa) + 48*invdmstau2mu
     *logmstau22Q2*mH2*power2(sa) + (216*mH2*power2(sa))/mA2 - (48*logmH2Q2*mH2
     *power2(sa))/mA2 - (1830*mH2*power2(sa))/mh2 - (216*logmh2Q2*mH2*power2(sa
     ))/mh2 + (1680*logmH2Q2*mH2*power2(sa))/mh2 + (96*logmh2Q2*logmH2Q2*mH2*
     power2(sa))/mh2 + (432*msntau2*power2(sa))/mh2 - (144*logmh2Q2*msntau2*
     power2(sa))/mh2 - (288*logmsntau2Q2*msntau2*power2(sa))/mh2 + (72*logmh2Q2
     *logmsntau2Q2*msntau2*power2(sa))/mh2 - (432*msntau2*power2(sa))/mH2 + (
     144*logmH2Q2*msntau2*power2(sa))/mH2 + (288*logmsntau2Q2*msntau2*power2(sa
     ))/mH2 - (72*logmH2Q2*logmsntau2Q2*msntau2*power2(sa))/mH2 + (24*mh2*
     power2(sa))/mstau12 - (48*logmstau12Q2*mh2*power2(sa))/mstau12 - (24*mH2*
     power2(sa))/mstau12 + (48*logmstau12Q2*mH2*power2(sa))/mstau12 - 1008*
     invdmhH*mstau12*power2(sa) + 480*invdmhH*logmH2Q2*mstau12*power2(sa) + 576
     *invdmhH*logmstau12Q2*mstau12*power2(sa) - (576*mstau12*power2(sa))/mh2 -
     (144*logmh2Q2*mstau12*power2(sa))/mh2 + (480*logmH2Q2*mstau12*power2(sa))/
     mh2 + (288*logmstau12Q2*mstau12*power2(sa))/mh2 + (72*logmh2Q2*
     logmstau12Q2*mstau12*power2(sa))/mh2 - (432*mstau12*power2(sa))/mH2 + (144
     *logmH2Q2*mstau12*power2(sa))/mH2 + (288*logmstau12Q2*mstau12*power2(sa))/
     mH2 - (72*logmH2Q2*logmstau12Q2*mstau12*power2(sa))/mH2 + (2304*invdmhH*
     mH2*mstau12*power2(sa))/mh2 - (1536*invdmhH*logmH2Q2*mH2*mstau12*power2(sa
     ))/mh2 - (768*invdmhH*logmstau12Q2*mH2*mstau12*power2(sa))/mh2 + (24*mh2*
     power2(sa))/mstau22 - (48*logmstau22Q2*mh2*power2(sa))/mstau22 - (24*mH2*
     power2(sa))/mstau22 + (48*logmstau22Q2*mH2*power2(sa))/mstau22 - 1008*
     invdmhH*mstau22*power2(sa) + 480*invdmhH*logmH2Q2*mstau22*power2(sa) + 576
     *invdmhH*logmstau22Q2*mstau22*power2(sa) - (576*mstau22*power2(sa))/mh2 -
     (144*logmh2Q2*mstau22*power2(sa))/mh2 + (480*logmH2Q2*mstau22*power2(sa))/
     mh2 + (288*logmstau22Q2*mstau22*power2(sa))/mh2 + (72*logmh2Q2*
     logmstau22Q2*mstau22*power2(sa))/mh2 - (432*mstau22*power2(sa))/mH2 + (144
     *logmH2Q2*mstau22*power2(sa))/mH2 + (288*logmstau22Q2*mstau22*power2(sa))/
     mH2 - (72*logmH2Q2*logmstau22Q2*mstau22*power2(sa))/mH2 + (2304*invdmhH*
     mH2*mstau22*power2(sa))/mh2 - (1536*invdmhH*logmH2Q2*mH2*mstau22*power2(sa
     ))/mh2 - (768*invdmhH*logmstau22Q2*mH2*mstau22*power2(sa))/mh2 - 648*
     invdmhH*mu2*power2(sa) + 1152*invdmhH*logmh2Q2*mu2*power2(sa) - 96*
     invdmsntau2mu*logmh2Q2*mu2*power2(sa) - 96*invdmstau1mu*logmh2Q2*mu2*
     power2(sa) - 96*invdmstau2mu*logmh2Q2*mu2*power2(sa) - 480*invdmhH*
     logmH2Q2*mu2*power2(sa) + 96*invdmsntau2mu*logmH2Q2*mu2*power2(sa) + 96*
     invdmstau1mu*logmH2Q2*mu2*power2(sa) + 96*invdmstau2mu*logmH2Q2*mu2*power2
     (sa) + 480*invdmhH*logmstau12Q2*mu2*power2(sa) - 576*invdmhH*logmh2Q2*
     logmstau12Q2*mu2*power2(sa) + 192*invdmstau1mu*logmh2Q2*logmstau12Q2*mu2*
     power2(sa) - 192*invdmstau1mu*logmH2Q2*logmstau12Q2*mu2*power2(sa) + 192*
     invdmstau2mu*logmh2Q2*logmstau22Q2*mu2*power2(sa) - 192*invdmstau2mu*
     logmH2Q2*logmstau22Q2*mu2*power2(sa) - 192*invdmstau1mu*logmh2Q2*logmu2Q2*
     mu2*power2(sa) - 192*invdmstau2mu*logmh2Q2*logmu2Q2*mu2*power2(sa) + 192*
     invdmstau1mu*logmH2Q2*logmu2Q2*mu2*power2(sa) + 192*invdmstau2mu*logmH2Q2*
     logmu2Q2*mu2*power2(sa) + (2808*mu2*power2(sa))/mh2 - (144*logmh2Q2*mu2*
     power2(sa))/mh2 + (672*logmH2Q2*mu2*power2(sa))/mh2 - (432*logmsntau2Q2*
     mu2*power2(sa))/mh2 - (1104*logmstau12Q2*mu2*power2(sa))/mh2 + (288*
     logmh2Q2*logmstau12Q2*mu2*power2(sa))/mh2 - (576*logmH2Q2*logmstau12Q2*mu2
     *power2(sa))/mh2 - (432*logmstau22Q2*mu2*power2(sa))/mh2 - (432*logmu2Q2*
     mu2*power2(sa))/mh2 - (216*logmh2Q2*logmu2Q2*mu2*power2(sa))/mh2 + (144*
     logmsntau2Q2*logmu2Q2*mu2*power2(sa))/mh2 + (144*logmstau12Q2*logmu2Q2*mu2
     *power2(sa))/mh2 + (144*logmstau22Q2*logmu2Q2*mu2*power2(sa))/mh2 - (1008*
     logmH2Q2*mu2*power2(sa))/mH2 + (432*logmsntau2Q2*mu2*power2(sa))/mH2 - (
     720*logmstau12Q2*mu2*power2(sa))/mH2 + (288*logmH2Q2*logmstau12Q2*mu2*
     power2(sa))/mH2 + (432*logmstau22Q2*mu2*power2(sa))/mH2 + (432*logmu2Q2*
     mu2*power2(sa))/mH2 + (216*logmH2Q2*logmu2Q2*mu2*power2(sa))/mH2 - (144*
     logmsntau2Q2*logmu2Q2*mu2*power2(sa))/mH2 - (144*logmstau12Q2*logmu2Q2*mu2
     *power2(sa))/mH2 - (144*logmstau22Q2*logmu2Q2*mu2*power2(sa))/mH2 - (768*
     invdmhH*logmstau12Q2*mH2*mu2*power2(sa))/mh2 + (768*invdmhH*logmH2Q2*
     logmstau12Q2*mH2*mu2*power2(sa))/mh2 + (144*logmh2Q2*mu2*power2(sa))/
     msntau2 - (144*logmH2Q2*mu2*power2(sa))/msntau2 - (96*logmh2Q2*
     logmsntau2Q2*mu2*power2(sa))/msntau2 + (96*logmH2Q2*logmsntau2Q2*mu2*
     power2(sa))/msntau2 + (48*logmh2Q2*logmu2Q2*mu2*power2(sa))/msntau2 - (48*
     logmH2Q2*logmu2Q2*mu2*power2(sa))/msntau2 + (4032*mu2*power2(sa))/mstau12
     - (720*logmh2Q2*mu2*power2(sa))/mstau12 - (1008*logmH2Q2*mu2*power2(sa))/
     mstau12 - (1728*logmstau12Q2*mu2*power2(sa))/mstau12 + (192*logmh2Q2*
     logmstau12Q2*mu2*power2(sa))/mstau12 + (384*logmH2Q2*logmstau12Q2*mu2*
     power2(sa))/mstau12 + (48*logmh2Q2*logmu2Q2*mu2*power2(sa))/mstau12 - (48*
     logmH2Q2*logmu2Q2*mu2*power2(sa))/mstau12 + (48*invdmstau1mu*mh2*mu2*
     power2(sa))/mstau12 + (672*invdmstau2mu*mh2*mu2*power2(sa))/mstau12 - (144
     *invdmstau2mu*logmh2Q2*mh2*mu2*power2(sa))/mstau12 - (432*invdmstau2mu*
     logmstau12Q2*mh2*mu2*power2(sa))/mstau12 + (96*invdmstau2mu*logmh2Q2*
     logmstau12Q2*mh2*mu2*power2(sa))/mstau12 - (48*invdmstau2mu*logmh2Q2*
     logmstau22Q2*mh2*mu2*power2(sa))/mstau12 + (48*invdmstau2mu*logmstau12Q2*
     logmstau22Q2*mh2*mu2*power2(sa))/mstau12 - (48*invdmstau1mu*mH2*mu2*power2
     (sa))/mstau12 - (672*invdmstau2mu*mH2*mu2*power2(sa))/mstau12 + (144*
     invdmstau2mu*logmH2Q2*mH2*mu2*power2(sa))/mstau12 + (432*invdmstau2mu*
     logmstau12Q2*mH2*mu2*power2(sa))/mstau12 - (96*invdmstau2mu*logmH2Q2*
     logmstau12Q2*mH2*mu2*power2(sa))/mstau12 + (48*invdmstau2mu*logmH2Q2*
     logmstau22Q2*mH2*mu2*power2(sa))/mstau12 - (48*invdmstau2mu*logmstau12Q2*
     logmstau22Q2*mH2*mu2*power2(sa))/mstau12 - (168*invdmhH*mstau12*mu2*power2
     (sa))/mh2 + (96*invdmhH*logmstau12Q2*mstau12*mu2*power2(sa))/mh2 + (144*
     logmh2Q2*mu2*power2(sa))/mstau22 - (144*logmH2Q2*mu2*power2(sa))/mstau22 -
     (96*logmh2Q2*logmstau22Q2*mu2*power2(sa))/mstau22 + (96*logmH2Q2*
     logmstau22Q2*mu2*power2(sa))/mstau22 + (48*logmh2Q2*logmu2Q2*mu2*power2(sa
     ))/mstau22 - (48*logmH2Q2*logmu2Q2*mu2*power2(sa))/mstau22 + (48*
     invdmstau2mu*mh2*mu2*power2(sa))/mstau22 - (48*invdmstau2mu*mH2*mu2*power2
     (sa))/mstau22 + 480*invdmhH*invdmstau*logmstau12Q2*mstau22*mu2*power2(sa)
     - 576*invdmhH*invdmstau*logmh2Q2*logmstau12Q2*mstau22*mu2*power2(sa) - 480
     *invdmhH*invdmstau*logmstau22Q2*mstau22*mu2*power2(sa) + 576*invdmhH*
     invdmstau*logmh2Q2*logmstau22Q2*mstau22*mu2*power2(sa) - (168*invdmhH*
     mstau22*mu2*power2(sa))/mh2 + (192*invdmstau*logmstau12Q2*mstau22*mu2*
     power2(sa))/mh2 - (576*invdmstau*logmH2Q2*logmstau12Q2*mstau22*mu2*power2(
     sa))/mh2 + (96*invdmhH*logmstau22Q2*mstau22*mu2*power2(sa))/mh2 - (192*
     invdmstau*logmstau22Q2*mstau22*mu2*power2(sa))/mh2 + (576*invdmstau*
     logmH2Q2*logmstau22Q2*mstau22*mu2*power2(sa))/mh2 - (288*invdmstau*
     logmstau12Q2*mstau22*mu2*power2(sa))/mH2 + (288*invdmstau*logmstau22Q2*
     mstau22*mu2*power2(sa))/mH2 - (768*invdmhH*invdmstau*logmstau12Q2*mH2*
     mstau22*mu2*power2(sa))/mh2 + (768*invdmhH*invdmstau*logmH2Q2*logmstau12Q2
     *mH2*mstau22*mu2*power2(sa))/mh2 + (768*invdmhH*invdmstau*logmstau22Q2*mH2
     *mstau22*mu2*power2(sa))/mh2 - (768*invdmhH*invdmstau*logmH2Q2*
     logmstau22Q2*mH2*mstau22*mu2*power2(sa))/mh2 + (4032*mstau22*mu2*power2(sa
     ))/(mh2*mstau12) - (2592*logmstau12Q2*mstau22*mu2*power2(sa))/(mh2*mstau12
     ) + (288*logmh2Q2*logmstau12Q2*mstau22*mu2*power2(sa))/(mh2*mstau12) - (
     864*logmstau22Q2*mstau22*mu2*power2(sa))/(mh2*mstau12) - (288*logmh2Q2*
     logmstau22Q2*mstau22*mu2*power2(sa))/(mh2*mstau12) + (576*logmstau12Q2*
     logmstau22Q2*mstau22*mu2*power2(sa))/(mh2*mstau12) + (4032*mstau22*mu2*
     power2(sa))/(mH2*mstau12) - (2592*logmstau12Q2*mstau22*mu2*power2(sa))/(
     mH2*mstau12) + (288*logmH2Q2*logmstau12Q2*mstau22*mu2*power2(sa))/(mH2*
     mstau12) - (864*logmstau22Q2*mstau22*mu2*power2(sa))/(mH2*mstau12) - (288*
     logmH2Q2*logmstau22Q2*mstau22*mu2*power2(sa))/(mH2*mstau12) + (576*
     logmstau12Q2*logmstau22Q2*mstau22*mu2*power2(sa))/(mH2*mstau12) + 768*mH2*
     mstau12*power2(invdmhH)*power2(sa) - 768*logmH2Q2*mH2*mstau12*power2(
     invdmhH)*power2(sa) + 768*mH2*mstau22*power2(invdmhH)*power2(sa) - 768*
     logmH2Q2*mH2*mstau22*power2(invdmhH)*power2(sa) - 768*mH2*mu2*power2(
     invdmhH)*power2(sa) + 768*logmH2Q2*mH2*mu2*power2(invdmhH)*power2(sa) -
     168*mh2*mu2*power2(invdmstau1mu)*power2(sa) + 48*logmh2Q2*mh2*mu2*power2(
     invdmstau1mu)*power2(sa) + 96*logmstau12Q2*mh2*mu2*power2(invdmstau1mu)*
     power2(sa) - 48*logmh2Q2*logmstau12Q2*mh2*mu2*power2(invdmstau1mu)*power2(
     sa) + 168*mH2*mu2*power2(invdmstau1mu)*power2(sa) - 48*logmH2Q2*mH2*mu2*
     power2(invdmstau1mu)*power2(sa) - 96*logmstau12Q2*mH2*mu2*power2(
     invdmstau1mu)*power2(sa) + 48*logmH2Q2*logmstau12Q2*mH2*mu2*power2(
     invdmstau1mu)*power2(sa) + 168*mh2*mu2*power2(invdmstau2mu)*power2(sa) -
     96*logmh2Q2*mh2*mu2*power2(invdmstau2mu)*power2(sa) - 48*logmstau22Q2*mh2*
     mu2*power2(invdmstau2mu)*power2(sa) - 168*mH2*mu2*power2(invdmstau2mu)*
     power2(sa) + 96*logmH2Q2*mH2*mu2*power2(invdmstau2mu)*power2(sa) + 48*
     logmstau22Q2*mH2*mu2*power2(invdmstau2mu)*power2(sa) - 72*power2(logmh2Q2)
     *power2(sa) - 12*invdmAh*mA2*power2(logmh2Q2)*power2(sa) + 12*invdmCh*mC2*
     power2(logmh2Q2)*power2(sa) + (48*mA2*power2(logmh2Q2)*power2(sa))/mh2 + (
     48*mC2*power2(logmh2Q2)*power2(sa))/mh2 + (84*mh2*power2(logmh2Q2)*power2(
     sa))/mA2 - (108*mh2*power2(logmh2Q2)*power2(sa))/mH2 - 276*invdmhH*mH2*
     power2(logmh2Q2)*power2(sa) + (48*mH2*power2(logmh2Q2)*power2(sa))/mh2 + (
     36*msntau2*power2(logmh2Q2)*power2(sa))/mh2 + (36*mstau12*power2(logmh2Q2)
     *power2(sa))/mh2 + (36*mstau22*power2(logmh2Q2)*power2(sa))/mh2 - 288*
     invdmhH*mu2*power2(logmh2Q2)*power2(sa) + (36*mu2*power2(logmh2Q2)*power2(
     sa))/mh2 - (24*mu2*power2(logmh2Q2)*power2(sa))/msntau2 + (120*mu2*power2(
     logmh2Q2)*power2(sa))/mstau12 + (24*invdmstau2mu*mh2*mu2*power2(logmh2Q2)*
     power2(sa))/mstau12 - (24*mu2*power2(logmh2Q2)*power2(sa))/mstau22 - 24*
     mh2*mu2*power2(invdmstau1mu)*power2(logmh2Q2)*power2(sa) - 648*power2(
     logmH2Q2)*power2(sa) + 12*invdmAH*mA2*power2(logmH2Q2)*power2(sa) - 12*
     invdmCH*mC2*power2(logmH2Q2)*power2(sa) - (48*mA2*power2(logmH2Q2)*power2(
     sa))/mH2 - (48*mC2*power2(logmH2Q2)*power2(sa))/mH2 + (36*mh2*power2(
     logmH2Q2)*power2(sa))/mH2 - 180*invdmhH*mH2*power2(logmH2Q2)*power2(sa) -
     (84*mH2*power2(logmH2Q2)*power2(sa))/mA2 - (576*mH2*power2(logmH2Q2)*
     power2(sa))/mh2 - (36*msntau2*power2(logmH2Q2)*power2(sa))/mH2 - 288*
     invdmhH*mstau12*power2(logmH2Q2)*power2(sa) - (288*mstau12*power2(logmH2Q2
     )*power2(sa))/mh2 - (36*mstau12*power2(logmH2Q2)*power2(sa))/mH2 + (768*
     invdmhH*mH2*mstau12*power2(logmH2Q2)*power2(sa))/mh2 - 288*invdmhH*mstau22
     *power2(logmH2Q2)*power2(sa) - (288*mstau22*power2(logmH2Q2)*power2(sa))/
     mh2 - (36*mstau22*power2(logmH2Q2)*power2(sa))/mH2 + (768*invdmhH*mH2*
     mstau22*power2(logmH2Q2)*power2(sa))/mh2 + 288*invdmhH*mu2*power2(logmH2Q2
     )*power2(sa) + (252*mu2*power2(logmH2Q2)*power2(sa))/mH2 - (384*invdmhH*
     mH2*mu2*power2(logmH2Q2)*power2(sa))/mh2 + (24*mu2*power2(logmH2Q2)*power2
     (sa))/msntau2 + (168*mu2*power2(logmH2Q2)*power2(sa))/mstau12 - (24*
     invdmstau2mu*mH2*mu2*power2(logmH2Q2)*power2(sa))/mstau12 + (24*mu2*power2
     (logmH2Q2)*power2(sa))/mstau22 + 384*mH2*mstau12*power2(invdmhH)*power2(
     logmH2Q2)*power2(sa) + 384*mH2*mstau22*power2(invdmhH)*power2(logmH2Q2)*
     power2(sa) - 384*mH2*mu2*power2(invdmhH)*power2(logmH2Q2)*power2(sa) + 24*
     mH2*mu2*power2(invdmstau1mu)*power2(logmH2Q2)*power2(sa) + (36*msntau2*
     power2(logmsntau2Q2)*power2(sa))/mh2 - (36*msntau2*power2(logmsntau2Q2)*
     power2(sa))/mH2 + (72*mu2*power2(logmsntau2Q2)*power2(sa))/mh2 - (72*mu2*
     power2(logmsntau2Q2)*power2(sa))/mH2 - 288*invdmhH*mstau12*power2(
     logmstau12Q2)*power2(sa) - (252*mstau12*power2(logmstau12Q2)*power2(sa))/
     mh2 - (36*mstau12*power2(logmstau12Q2)*power2(sa))/mH2 + (384*invdmhH*mH2*
     mstau12*power2(logmstau12Q2)*power2(sa))/mh2 + (216*mu2*power2(
     logmstau12Q2)*power2(sa))/mh2 + (72*mu2*power2(logmstau12Q2)*power2(sa))/
     mH2 + (288*mu2*power2(logmstau12Q2)*power2(sa))/mstau12 + (72*invdmstau2mu
     *mh2*mu2*power2(logmstau12Q2)*power2(sa))/mstau12 - (72*invdmstau2mu*mH2*
     mu2*power2(logmstau12Q2)*power2(sa))/mstau12 + (432*mstau22*mu2*power2(
     logmstau12Q2)*power2(sa))/(mh2*mstau12) + (432*mstau22*mu2*power2(
     logmstau12Q2)*power2(sa))/(mH2*mstau12) - 24*mh2*mu2*power2(invdmstau1mu)*
     power2(logmstau12Q2)*power2(sa) + 24*mH2*mu2*power2(invdmstau1mu)*power2(
     logmstau12Q2)*power2(sa) - 288*invdmhH*mstau22*power2(logmstau22Q2)*power2
     (sa) - (252*mstau22*power2(logmstau22Q2)*power2(sa))/mh2 - (36*mstau22*
     power2(logmstau22Q2)*power2(sa))/mH2 + (384*invdmhH*mH2*mstau22*power2(
     logmstau22Q2)*power2(sa))/mh2 + (72*mu2*power2(logmstau22Q2)*power2(sa))/
     mh2 - (72*mu2*power2(logmstau22Q2)*power2(sa))/mH2 + (144*mstau22*mu2*
     power2(logmstau22Q2)*power2(sa))/(mh2*mstau12) + (144*mstau22*mu2*power2(
     logmstau22Q2)*power2(sa))/(mH2*mstau12) + (108*mu2*power2(logmu2Q2)*power2
     (sa))/mh2 - (108*mu2*power2(logmu2Q2)*power2(sa))/mH2 + 105*power2(invdmAh
     )*power2(mA2)*power2(sa) - 72*logmh2Q2*power2(invdmAh)*power2(mA2)*power2(
     sa) - 105*power2(invdmAH)*power2(mA2)*power2(sa) + 72*logmH2Q2*power2(
     invdmAH)*power2(mA2)*power2(sa) + 12*power2(invdmAh)*power2(logmh2Q2)*
     power2(mA2)*power2(sa) - 12*power2(invdmAH)*power2(logmH2Q2)*power2(mA2)*
     power2(sa) + 42*power2(invdmCh)*power2(mC2)*power2(sa) - 72*logmh2Q2*
     power2(invdmCh)*power2(mC2)*power2(sa) - 42*power2(invdmCH)*power2(mC2)*
     power2(sa) + 72*logmH2Q2*power2(invdmCH)*power2(mC2)*power2(sa) + 12*
     power2(invdmCh)*power2(logmh2Q2)*power2(mC2)*power2(sa) - 12*power2(
     invdmCH)*power2(logmH2Q2)*power2(mC2)*power2(sa) + (528*mH2*mstau12*power2
     (sa))/power2(mh2) - (288*logmH2Q2*mH2*mstau12*power2(sa))/power2(mh2) - (
     192*logmstau12Q2*mH2*mstau12*power2(sa))/power2(mh2) + (528*mH2*mstau22*
     power2(sa))/power2(mh2) - (288*logmH2Q2*mH2*mstau22*power2(sa))/power2(mh2
     ) - (192*logmstau22Q2*mH2*mstau22*power2(sa))/power2(mh2) + (120*mH2*mu2*
     power2(sa))/power2(mh2) - (96*logmH2Q2*mH2*mu2*power2(sa))/power2(mh2) - (
     288*logmstau12Q2*mH2*mu2*power2(sa))/power2(mh2) + (192*logmH2Q2*
     logmstau12Q2*mH2*mu2*power2(sa))/power2(mh2) - (168*mstau12*mu2*power2(sa)
     )/power2(mh2) + (96*logmstau12Q2*mstau12*mu2*power2(sa))/power2(mh2) + (
     168*invdmhH*mH2*mstau12*mu2*power2(sa))/power2(mh2) - (96*invdmhH*
     logmstau12Q2*mH2*mstau12*mu2*power2(sa))/power2(mh2) - (168*mstau22*mu2*
     power2(sa))/power2(mh2) + (96*logmstau22Q2*mstau22*mu2*power2(sa))/power2(
     mh2) + (168*invdmhH*mH2*mstau22*mu2*power2(sa))/power2(mh2) - (288*
     invdmstau*logmstau12Q2*mH2*mstau22*mu2*power2(sa))/power2(mh2) + (192*
     invdmstau*logmH2Q2*logmstau12Q2*mH2*mstau22*mu2*power2(sa))/power2(mh2) -
     (96*invdmhH*logmstau22Q2*mH2*mstau22*mu2*power2(sa))/power2(mh2) + (288*
     invdmstau*logmstau22Q2*mH2*mstau22*mu2*power2(sa))/power2(mh2) - (192*
     invdmstau*logmH2Q2*logmstau22Q2*mH2*mstau22*mu2*power2(sa))/power2(mh2) +
     (96*mH2*mstau12*power2(logmH2Q2)*power2(sa))/power2(mh2) + (96*mH2*mstau22
     *power2(logmH2Q2)*power2(sa))/power2(mh2) + (96*mH2*mstau12*power2(
     logmstau12Q2)*power2(sa))/power2(mh2) + (96*mH2*mstau22*power2(
     logmstau22Q2)*power2(sa))/power2(mh2) + (315*power2(mA2)*power2(sa))/
     power2(mh2) + (315*power2(mC2)*power2(sa))/power2(mh2) - (420*power2(mh2)*
     power2(sa))/power2(mA2) + (360*logmh2Q2*power2(mh2)*power2(sa))/power2(mA2
     ) - (120*power2(logmh2Q2)*power2(mh2)*power2(sa))/power2(mA2) - (315*
     power2(mA2)*power2(sa))/power2(mH2) - (315*power2(mC2)*power2(sa))/power2(
     mH2) + (252*power2(mh2)*power2(sa))/power2(mH2) - (216*logmh2Q2*power2(mh2
     )*power2(sa))/power2(mH2) + (72*power2(logmh2Q2)*power2(mh2)*power2(sa))/
     power2(mH2) + (2472*invdmhH*power2(mH2)*power2(sa))/mh2 - (2064*invdmhH*
     logmH2Q2*power2(mH2)*power2(sa))/mh2 - 21*power2(invdmhH)*power2(mH2)*
     power2(sa) - 72*logmh2Q2*power2(invdmhH)*power2(mH2)*power2(sa) + 36*
     logmH2Q2*power2(invdmhH)*power2(mH2)*power2(sa) + 24*logmh2Q2*logmH2Q2*
     power2(invdmhH)*power2(mH2)*power2(sa) - (1104*mstau12*power2(invdmhH)*
     power2(mH2)*power2(sa))/mh2 + (1056*logmH2Q2*mstau12*power2(invdmhH)*
     power2(mH2)*power2(sa))/mh2 - (1104*mstau22*power2(invdmhH)*power2(mH2)*
     power2(sa))/mh2 + (1056*logmH2Q2*mstau22*power2(invdmhH)*power2(mH2)*
     power2(sa))/mh2 + (1104*mu2*power2(invdmhH)*power2(mH2)*power2(sa))/mh2 -
     (1056*logmH2Q2*mu2*power2(invdmhH)*power2(mH2)*power2(sa))/mh2 + 12*power2
     (invdmhH)*power2(logmh2Q2)*power2(mH2)*power2(sa) + (624*invdmhH*power2(
     logmH2Q2)*power2(mH2)*power2(sa))/mh2 + 12*power2(invdmhH)*power2(logmH2Q2
     )*power2(mH2)*power2(sa) - (480*mstau12*power2(invdmhH)*power2(logmH2Q2)*
     power2(mH2)*power2(sa))/mh2 - (480*mstau22*power2(invdmhH)*power2(logmH2Q2
     )*power2(mH2)*power2(sa))/mh2 + (480*mu2*power2(invdmhH)*power2(logmH2Q2)*
     power2(mH2)*power2(sa))/mh2 + (420*power2(mH2)*power2(sa))/power2(mA2) - (
     360*logmH2Q2*power2(mH2)*power2(sa))/power2(mA2) + (120*power2(logmH2Q2)*
     power2(mH2)*power2(sa))/power2(mA2) + (987*power2(mH2)*power2(sa))/power2(
     mh2) - (780*logmH2Q2*power2(mH2)*power2(sa))/power2(mh2) - (864*invdmhH*
     mstau12*power2(mH2)*power2(sa))/power2(mh2) + (576*invdmhH*logmH2Q2*
     mstau12*power2(mH2)*power2(sa))/power2(mh2) + (192*invdmhH*logmstau12Q2*
     mstau12*power2(mH2)*power2(sa))/power2(mh2) - (864*invdmhH*mstau22*power2(
     mH2)*power2(sa))/power2(mh2) + (576*invdmhH*logmH2Q2*mstau22*power2(mH2)*
     power2(sa))/power2(mh2) + (192*invdmhH*logmstau22Q2*mstau22*power2(mH2)*
     power2(sa))/power2(mh2) + (216*invdmhH*mu2*power2(mH2)*power2(sa))/power2(
     mh2) - (192*invdmhH*logmH2Q2*mu2*power2(mH2)*power2(sa))/power2(mh2) + (
     288*invdmhH*logmstau12Q2*mu2*power2(mH2)*power2(sa))/power2(mh2) - (192*
     invdmhH*logmH2Q2*logmstau12Q2*mu2*power2(mH2)*power2(sa))/power2(mh2) + (
     288*invdmhH*invdmstau*logmstau12Q2*mstau22*mu2*power2(mH2)*power2(sa))/
     power2(mh2) - (192*invdmhH*invdmstau*logmH2Q2*logmstau12Q2*mstau22*mu2*
     power2(mH2)*power2(sa))/power2(mh2) - (288*invdmhH*invdmstau*logmstau22Q2*
     mstau22*mu2*power2(mH2)*power2(sa))/power2(mh2) + (192*invdmhH*invdmstau*
     logmH2Q2*logmstau22Q2*mstau22*mu2*power2(mH2)*power2(sa))/power2(mh2) + (
     216*power2(logmH2Q2)*power2(mH2)*power2(sa))/power2(mh2) - (192*invdmhH*
     mstau12*power2(logmH2Q2)*power2(mH2)*power2(sa))/power2(mh2) - (192*
     invdmhH*mstau22*power2(logmH2Q2)*power2(mH2)*power2(sa))/power2(mh2) + (96
     *invdmhH*mu2*power2(logmH2Q2)*power2(mH2)*power2(sa))/power2(mh2) - (96*
     invdmhH*mstau12*power2(logmstau12Q2)*power2(mH2)*power2(sa))/power2(mh2) -
     (96*invdmhH*mstau22*power2(logmstau22Q2)*power2(mH2)*power2(sa))/power2(
     mh2) - (1344*invdmsntau2mu*power2(mu2)*power2(sa))/mh2 - (1344*
     invdmstau1mu*power2(mu2)*power2(sa))/mh2 - (1344*invdmstau2mu*power2(mu2)*
     power2(sa))/mh2 + (576*invdmsntau2mu*logmsntau2Q2*power2(mu2)*power2(sa))/
     mh2 + (576*invdmstau1mu*logmstau12Q2*power2(mu2)*power2(sa))/mh2 + (576*
     invdmstau2mu*logmstau22Q2*power2(mu2)*power2(sa))/mh2 + (576*invdmsntau2mu
     *logmu2Q2*power2(mu2)*power2(sa))/mh2 + (576*invdmstau1mu*logmu2Q2*power2(
     mu2)*power2(sa))/mh2 + (576*invdmstau2mu*logmu2Q2*power2(mu2)*power2(sa))/
     mh2 - (192*invdmsntau2mu*logmsntau2Q2*logmu2Q2*power2(mu2)*power2(sa))/mh2
      - (192*invdmstau1mu*logmstau12Q2*logmu2Q2*power2(mu2)*power2(sa))/mh2 - (
     192*invdmstau2mu*logmstau22Q2*logmu2Q2*power2(mu2)*power2(sa))/mh2 + (1344
     *invdmsntau2mu*power2(mu2)*power2(sa))/mH2 + (1344*invdmstau1mu*power2(mu2
     )*power2(sa))/mH2 + (1344*invdmstau2mu*power2(mu2)*power2(sa))/mH2 - (576*
     invdmsntau2mu*logmsntau2Q2*power2(mu2)*power2(sa))/mH2 - (576*invdmstau1mu
     *logmstau12Q2*power2(mu2)*power2(sa))/mH2 - (576*invdmstau2mu*logmstau22Q2
     *power2(mu2)*power2(sa))/mH2 - (576*invdmsntau2mu*logmu2Q2*power2(mu2)*
     power2(sa))/mH2 - (576*invdmstau1mu*logmu2Q2*power2(mu2)*power2(sa))/mH2 -
     (576*invdmstau2mu*logmu2Q2*power2(mu2)*power2(sa))/mH2 + (192*
     invdmsntau2mu*logmsntau2Q2*logmu2Q2*power2(mu2)*power2(sa))/mH2 + (192*
     invdmstau1mu*logmstau12Q2*logmu2Q2*power2(mu2)*power2(sa))/mH2 + (192*
     invdmstau2mu*logmstau22Q2*logmu2Q2*power2(mu2)*power2(sa))/mH2 - (336*
     power2(mu2)*power2(sa))/(mh2*msntau2) + (144*logmsntau2Q2*power2(mu2)*
     power2(sa))/(mh2*msntau2) + (144*logmu2Q2*power2(mu2)*power2(sa))/(mh2*
     msntau2) - (48*logmsntau2Q2*logmu2Q2*power2(mu2)*power2(sa))/(mh2*msntau2)
     + (336*power2(mu2)*power2(sa))/(mH2*msntau2) - (144*logmsntau2Q2*power2(
     mu2)*power2(sa))/(mH2*msntau2) - (144*logmu2Q2*power2(mu2)*power2(sa))/(
     mH2*msntau2) + (48*logmsntau2Q2*logmu2Q2*power2(mu2)*power2(sa))/(mH2*
     msntau2) - (336*power2(mu2)*power2(sa))/(mh2*mstau12) + (144*logmstau12Q2*
     power2(mu2)*power2(sa))/(mh2*mstau12) + (144*logmu2Q2*power2(mu2)*power2(
     sa))/(mh2*mstau12) - (48*logmstau12Q2*logmu2Q2*power2(mu2)*power2(sa))/(
     mh2*mstau12) + (336*power2(mu2)*power2(sa))/(mH2*mstau12) - (144*
     logmstau12Q2*power2(mu2)*power2(sa))/(mH2*mstau12) - (144*logmu2Q2*power2(
     mu2)*power2(sa))/(mH2*mstau12) + (48*logmstau12Q2*logmu2Q2*power2(mu2)*
     power2(sa))/(mH2*mstau12) - (336*power2(mu2)*power2(sa))/(mh2*mstau22) + (
     144*logmstau22Q2*power2(mu2)*power2(sa))/(mh2*mstau22) + (144*logmu2Q2*
     power2(mu2)*power2(sa))/(mh2*mstau22) - (48*logmstau22Q2*logmu2Q2*power2(
     mu2)*power2(sa))/(mh2*mstau22) + (336*power2(mu2)*power2(sa))/(mH2*mstau22
     ) - (144*logmstau22Q2*power2(mu2)*power2(sa))/(mH2*mstau22) - (144*
     logmu2Q2*power2(mu2)*power2(sa))/(mH2*mstau22) + (48*logmstau22Q2*logmu2Q2
     *power2(mu2)*power2(sa))/(mH2*mstau22) - 96*logmh2Q2*logmsntau2Q2*power2(
     invdmsntau2mu)*power2(mu2)*power2(sa) + 96*logmH2Q2*logmsntau2Q2*power2(
     invdmsntau2mu)*power2(mu2)*power2(sa) + 96*logmh2Q2*logmu2Q2*power2(
     invdmsntau2mu)*power2(mu2)*power2(sa) - 96*logmH2Q2*logmu2Q2*power2(
     invdmsntau2mu)*power2(mu2)*power2(sa) - 96*logmh2Q2*logmstau12Q2*power2(
     invdmstau1mu)*power2(mu2)*power2(sa) + 96*logmH2Q2*logmstau12Q2*power2(
     invdmstau1mu)*power2(mu2)*power2(sa) + 96*logmh2Q2*logmu2Q2*power2(
     invdmstau1mu)*power2(mu2)*power2(sa) - 96*logmH2Q2*logmu2Q2*power2(
     invdmstau1mu)*power2(mu2)*power2(sa) - (72*mh2*power2(invdmstau1mu)*power2
     (mu2)*power2(sa))/mstau12 + (48*logmstau12Q2*mh2*power2(invdmstau1mu)*
     power2(mu2)*power2(sa))/mstau12 + (72*mH2*power2(invdmstau1mu)*power2(mu2)
     *power2(sa))/mstau12 - (48*logmstau12Q2*mH2*power2(invdmstau1mu)*power2(
     mu2)*power2(sa))/mstau12 + (672*mh2*power2(invdmstau1mu)*power2(mu2)*
     power2(sa))/mstau22 - (144*logmh2Q2*mh2*power2(invdmstau1mu)*power2(mu2)*
     power2(sa))/mstau22 - (432*logmstau22Q2*mh2*power2(invdmstau1mu)*power2(
     mu2)*power2(sa))/mstau22 + (96*logmh2Q2*logmstau22Q2*mh2*power2(
     invdmstau1mu)*power2(mu2)*power2(sa))/mstau22 - (48*logmh2Q2*logmu2Q2*mh2*
     power2(invdmstau1mu)*power2(mu2)*power2(sa))/mstau22 + (48*logmstau22Q2*
     logmu2Q2*mh2*power2(invdmstau1mu)*power2(mu2)*power2(sa))/mstau22 - (672*
     mH2*power2(invdmstau1mu)*power2(mu2)*power2(sa))/mstau22 + (144*logmH2Q2*
     mH2*power2(invdmstau1mu)*power2(mu2)*power2(sa))/mstau22 + (432*
     logmstau22Q2*mH2*power2(invdmstau1mu)*power2(mu2)*power2(sa))/mstau22 - (
     96*logmH2Q2*logmstau22Q2*mH2*power2(invdmstau1mu)*power2(mu2)*power2(sa))/
     mstau22 + (48*logmH2Q2*logmu2Q2*mH2*power2(invdmstau1mu)*power2(mu2)*
     power2(sa))/mstau22 - (48*logmstau22Q2*logmu2Q2*mH2*power2(invdmstau1mu)*
     power2(mu2)*power2(sa))/mstau22 - 96*logmh2Q2*logmstau22Q2*power2(
     invdmstau2mu)*power2(mu2)*power2(sa) + 96*logmH2Q2*logmstau22Q2*power2(
     invdmstau2mu)*power2(mu2)*power2(sa) + 96*logmh2Q2*logmu2Q2*power2(
     invdmstau2mu)*power2(mu2)*power2(sa) - 96*logmH2Q2*logmu2Q2*power2(
     invdmstau2mu)*power2(mu2)*power2(sa) + (1344*mh2*power2(invdmstau2mu)*
     power2(mu2)*power2(sa))/mstau12 - (288*logmh2Q2*mh2*power2(invdmstau2mu)*
     power2(mu2)*power2(sa))/mstau12 - (864*logmstau12Q2*mh2*power2(
     invdmstau2mu)*power2(mu2)*power2(sa))/mstau12 + (192*logmh2Q2*logmstau12Q2
     *mh2*power2(invdmstau2mu)*power2(mu2)*power2(sa))/mstau12 - (48*logmh2Q2*
     logmstau22Q2*mh2*power2(invdmstau2mu)*power2(mu2)*power2(sa))/mstau12 + (
     48*logmstau12Q2*logmstau22Q2*mh2*power2(invdmstau2mu)*power2(mu2)*power2(
     sa))/mstau12 - (48*logmh2Q2*logmu2Q2*mh2*power2(invdmstau2mu)*power2(mu2)*
     power2(sa))/mstau12 + (48*logmstau12Q2*logmu2Q2*mh2*power2(invdmstau2mu)*
     power2(mu2)*power2(sa))/mstau12 - (1344*mH2*power2(invdmstau2mu)*power2(
     mu2)*power2(sa))/mstau12 + (288*logmH2Q2*mH2*power2(invdmstau2mu)*power2(
     mu2)*power2(sa))/mstau12 + (864*logmstau12Q2*mH2*power2(invdmstau2mu)*
     power2(mu2)*power2(sa))/mstau12 - (192*logmH2Q2*logmstau12Q2*mH2*power2(
     invdmstau2mu)*power2(mu2)*power2(sa))/mstau12 + (48*logmH2Q2*logmstau22Q2*
     mH2*power2(invdmstau2mu)*power2(mu2)*power2(sa))/mstau12 - (48*
     logmstau12Q2*logmstau22Q2*mH2*power2(invdmstau2mu)*power2(mu2)*power2(sa))
     /mstau12 + (48*logmH2Q2*logmu2Q2*mH2*power2(invdmstau2mu)*power2(mu2)*
     power2(sa))/mstau12 - (48*logmstau12Q2*logmu2Q2*mH2*power2(invdmstau2mu)*
     power2(mu2)*power2(sa))/mstau12 - (72*mh2*power2(invdmstau2mu)*power2(mu2)
     *power2(sa))/mstau22 + (48*logmstau22Q2*mh2*power2(invdmstau2mu)*power2(
     mu2)*power2(sa))/mstau22 + (72*mH2*power2(invdmstau2mu)*power2(mu2)*power2
     (sa))/mstau22 - (48*logmstau22Q2*mH2*power2(invdmstau2mu)*power2(mu2)*
     power2(sa))/mstau22 + (24*mh2*power2(invdmstau1mu)*power2(logmh2Q2)*power2
     (mu2)*power2(sa))/mstau22 + (48*mh2*power2(invdmstau2mu)*power2(logmh2Q2)*
     power2(mu2)*power2(sa))/mstau12 - (24*mH2*power2(invdmstau1mu)*power2(
     logmH2Q2)*power2(mu2)*power2(sa))/mstau22 - (48*mH2*power2(invdmstau2mu)*
     power2(logmH2Q2)*power2(mu2)*power2(sa))/mstau12 - (96*invdmsntau2mu*
     power2(logmsntau2Q2)*power2(mu2)*power2(sa))/mh2 + (96*invdmsntau2mu*
     power2(logmsntau2Q2)*power2(mu2)*power2(sa))/mH2 - (24*power2(logmsntau2Q2
     )*power2(mu2)*power2(sa))/(mh2*msntau2) + (24*power2(logmsntau2Q2)*power2(
     mu2)*power2(sa))/(mH2*msntau2) - (96*invdmstau1mu*power2(logmstau12Q2)*
     power2(mu2)*power2(sa))/mh2 + (96*invdmstau1mu*power2(logmstau12Q2)*power2
     (mu2)*power2(sa))/mH2 - (24*power2(logmstau12Q2)*power2(mu2)*power2(sa))/(
     mh2*mstau12) + (24*power2(logmstau12Q2)*power2(mu2)*power2(sa))/(mH2*
     mstau12) + (144*mh2*power2(invdmstau2mu)*power2(logmstau12Q2)*power2(mu2)*
     power2(sa))/mstau12 - (144*mH2*power2(invdmstau2mu)*power2(logmstau12Q2)*
     power2(mu2)*power2(sa))/mstau12 - (96*invdmstau2mu*power2(logmstau22Q2)*
     power2(mu2)*power2(sa))/mh2 + (96*invdmstau2mu*power2(logmstau22Q2)*power2
     (mu2)*power2(sa))/mH2 - (24*power2(logmstau22Q2)*power2(mu2)*power2(sa))/(
     mh2*mstau22) + (24*power2(logmstau22Q2)*power2(mu2)*power2(sa))/(mH2*
     mstau22) + (72*mh2*power2(invdmstau1mu)*power2(logmstau22Q2)*power2(mu2)*
     power2(sa))/mstau22 - (72*mH2*power2(invdmstau1mu)*power2(logmstau22Q2)*
     power2(mu2)*power2(sa))/mstau22 - (96*invdmsntau2mu*power2(logmu2Q2)*
     power2(mu2)*power2(sa))/mh2 - (96*invdmstau1mu*power2(logmu2Q2)*power2(mu2
     )*power2(sa))/mh2 - (96*invdmstau2mu*power2(logmu2Q2)*power2(mu2)*power2(
     sa))/mh2 + (96*invdmsntau2mu*power2(logmu2Q2)*power2(mu2)*power2(sa))/mH2
     + (96*invdmstau1mu*power2(logmu2Q2)*power2(mu2)*power2(sa))/mH2 + (96*
     invdmstau2mu*power2(logmu2Q2)*power2(mu2)*power2(sa))/mH2 - (24*power2(
     logmu2Q2)*power2(mu2)*power2(sa))/(mh2*msntau2) + (24*power2(logmu2Q2)*
     power2(mu2)*power2(sa))/(mH2*msntau2) - (24*power2(logmu2Q2)*power2(mu2)*
     power2(sa))/(mh2*mstau12) + (24*power2(logmu2Q2)*power2(mu2)*power2(sa))/(
     mH2*mstau12) - (24*power2(logmu2Q2)*power2(mu2)*power2(sa))/(mh2*mstau22)
     + (24*power2(logmu2Q2)*power2(mu2)*power2(sa))/(mH2*mstau22) + 12*power2(
     Pi)*power2(sa) - 4*invdmAh*mA2*power2(Pi)*power2(sa) + 4*invdmAH*mA2*
     power2(Pi)*power2(sa) + 4*invdmCh*mC2*power2(Pi)*power2(sa) - 4*invdmCH*
     mC2*power2(Pi)*power2(sa) + (52*mA2*power2(Pi)*power2(sa))/mh2 - (32*mC2*
     power2(Pi)*power2(sa))/mh2 + (16*mh2*power2(Pi)*power2(sa))/mA2 - (52*mA2*
     power2(Pi)*power2(sa))/mH2 + (32*mC2*power2(Pi)*power2(sa))/mH2 - (24*mh2*
     power2(Pi)*power2(sa))/mH2 - 124*invdmhH*mH2*power2(Pi)*power2(sa) - (16*
     mH2*power2(Pi)*power2(sa))/mA2 - (148*mH2*power2(Pi)*power2(sa))/mh2 + (12
     *msntau2*power2(Pi)*power2(sa))/mh2 - (12*msntau2*power2(Pi)*power2(sa))/
     mH2 - 96*invdmhH*mstau12*power2(Pi)*power2(sa) - (84*mstau12*power2(Pi)*
     power2(sa))/mh2 - (12*mstau12*power2(Pi)*power2(sa))/mH2 + (192*invdmhH*
     mH2*mstau12*power2(Pi)*power2(sa))/mh2 - 96*invdmhH*mstau22*power2(Pi)*
     power2(sa) - (84*mstau22*power2(Pi)*power2(sa))/mh2 - (12*mstau22*power2(
     Pi)*power2(sa))/mH2 + (192*invdmhH*mH2*mstau22*power2(Pi)*power2(sa))/mh2
     + (84*mu2*power2(Pi)*power2(sa))/mh2 + (12*mu2*power2(Pi)*power2(sa))/mH2
     - (64*invdmhH*mH2*mu2*power2(Pi)*power2(sa))/mh2 + (96*mu2*power2(Pi)*
     power2(sa))/mstau12 + (16*invdmstau2mu*mh2*mu2*power2(Pi)*power2(sa))/
     mstau12 - (16*invdmstau2mu*mH2*mu2*power2(Pi)*power2(sa))/mstau12 + (96*
     mstau22*mu2*power2(Pi)*power2(sa))/(mh2*mstau12) + (96*mstau22*mu2*power2(
     Pi)*power2(sa))/(mH2*mstau12) + 64*mH2*mstau12*power2(invdmhH)*power2(Pi)*
     power2(sa) + 64*mH2*mstau22*power2(invdmhH)*power2(Pi)*power2(sa) - 64*mH2
     *mu2*power2(invdmhH)*power2(Pi)*power2(sa) - 8*mh2*mu2*power2(invdmstau1mu
     )*power2(Pi)*power2(sa) + 8*mH2*mu2*power2(invdmstau1mu)*power2(Pi)*power2
     (sa) + 4*power2(invdmAh)*power2(mA2)*power2(Pi)*power2(sa) - 4*power2(
     invdmAH)*power2(mA2)*power2(Pi)*power2(sa) + 4*power2(invdmCh)*power2(mC2)
     *power2(Pi)*power2(sa) - 4*power2(invdmCH)*power2(mC2)*power2(Pi)*power2(
     sa) + (32*mH2*mstau12*power2(Pi)*power2(sa))/power2(mh2) + (32*mH2*mstau22
     *power2(Pi)*power2(sa))/power2(mh2) + (18*power2(mA2)*power2(Pi)*power2(sa
     ))/power2(mh2) + (18*power2(mC2)*power2(Pi)*power2(sa))/power2(mh2) - (30*
     power2(mh2)*power2(Pi)*power2(sa))/power2(mA2) - (18*power2(mA2)*power2(Pi
     )*power2(sa))/power2(mH2) - (18*power2(mC2)*power2(Pi)*power2(sa))/power2(
     mH2) + (18*power2(mh2)*power2(Pi)*power2(sa))/power2(mH2) + (168*invdmhH*
     power2(mH2)*power2(Pi)*power2(sa))/mh2 + 4*power2(invdmhH)*power2(mH2)*
     power2(Pi)*power2(sa) - (80*mstau12*power2(invdmhH)*power2(mH2)*power2(Pi)
     *power2(sa))/mh2 - (80*mstau22*power2(invdmhH)*power2(mH2)*power2(Pi)*
     power2(sa))/mh2 + (80*mu2*power2(invdmhH)*power2(mH2)*power2(Pi)*power2(sa
     ))/mh2 + (30*power2(mH2)*power2(Pi)*power2(sa))/power2(mA2) + (58*power2(
     mH2)*power2(Pi)*power2(sa))/power2(mh2) - (48*invdmhH*mstau12*power2(mH2)*
     power2(Pi)*power2(sa))/power2(mh2) - (48*invdmhH*mstau22*power2(mH2)*
     power2(Pi)*power2(sa))/power2(mh2) + (16*invdmhH*mu2*power2(mH2)*power2(Pi
     )*power2(sa))/power2(mh2) - (32*invdmsntau2mu*power2(mu2)*power2(Pi)*
     power2(sa))/mh2 - (32*invdmstau1mu*power2(mu2)*power2(Pi)*power2(sa))/mh2
     - (32*invdmstau2mu*power2(mu2)*power2(Pi)*power2(sa))/mh2 + (32*
     invdmsntau2mu*power2(mu2)*power2(Pi)*power2(sa))/mH2 + (32*invdmstau1mu*
     power2(mu2)*power2(Pi)*power2(sa))/mH2 + (32*invdmstau2mu*power2(mu2)*
     power2(Pi)*power2(sa))/mH2 - (8*power2(mu2)*power2(Pi)*power2(sa))/(mh2*
     msntau2) + (8*power2(mu2)*power2(Pi)*power2(sa))/(mH2*msntau2) - (8*power2
     (mu2)*power2(Pi)*power2(sa))/(mh2*mstau12) + (8*power2(mu2)*power2(Pi)*
     power2(sa))/(mH2*mstau12) - (8*power2(mu2)*power2(Pi)*power2(sa))/(mh2*
     mstau22) + (8*power2(mu2)*power2(Pi)*power2(sa))/(mH2*mstau22) + (16*mh2*
     power2(invdmstau1mu)*power2(mu2)*power2(Pi)*power2(sa))/mstau22 - (16*mH2*
     power2(invdmstau1mu)*power2(mu2)*power2(Pi)*power2(sa))/mstau22 + (32*mh2*
     power2(invdmstau2mu)*power2(mu2)*power2(Pi)*power2(sa))/mstau12 - (32*mH2*
     power2(invdmstau2mu)*power2(mu2)*power2(Pi)*power2(sa))/mstau12 + 192*
     invdmstau*logmstau12Q2*logmz2Q2*mu2*power2(sb) - 192*invdmstau*
     logmstau22Q2*logmz2Q2*mu2*power2(sb) - 192*invdmstau*invdmstau1mu*
     logmstau12Q2*logmz2Q2*power2(mu2)*power2(sb) + 192*invdmstau*invdmstau2mu*
     logmstau22Q2*logmz2Q2*power2(mu2)*power2(sb) + 192*invdmstau*invdmstau1mu*
     logmu2Q2*logmz2Q2*power2(mu2)*power2(sb) - 192*invdmstau*invdmstau2mu*
     logmu2Q2*logmz2Q2*power2(mu2)*power2(sb) + 192*logmsntau2Q2*mC2*power2(mu2
     )*power3(invdmsntau2mu) - 192*logmu2Q2*mC2*power2(mu2)*power3(
     invdmsntau2mu) - 288*logmsntau2Q2*msntau2*power2(mu2)*power3(invdmsntau2mu
     ) + 288*logmu2Q2*msntau2*power2(mu2)*power3(invdmsntau2mu) - 96*
     logmsntau2Q2*logmstau22Q2*mstau22*power2(mu2)*power3(invdmsntau2mu) + 96*
     logmstau22Q2*logmu2Q2*mstau22*power2(mu2)*power3(invdmsntau2mu) + 96*
     msntau2*power2(logmsntau2Q2)*power2(mu2)*power3(invdmsntau2mu) + 48*
     mstau22*power2(logmsntau2Q2)*power2(mu2)*power3(invdmsntau2mu) - 96*
     msntau2*power2(logmu2Q2)*power2(mu2)*power3(invdmsntau2mu) - 48*mstau22*
     power2(logmu2Q2)*power2(mu2)*power3(invdmsntau2mu) + 96*logmstau12Q2*mA2*
     power2(mu2)*power3(invdmstau1mu) - 96*logmu2Q2*mA2*power2(mu2)*power3(
     invdmstau1mu) + 96*logmstau12Q2*mH2*power2(mu2)*power3(invdmstau1mu) - 96*
     logmH2Q2*logmstau12Q2*mH2*power2(mu2)*power3(invdmstau1mu) - 96*logmu2Q2*
     mH2*power2(mu2)*power3(invdmstau1mu) + 96*logmH2Q2*logmu2Q2*mH2*power2(mu2
     )*power3(invdmstau1mu) - 288*logmstau12Q2*mstau12*power2(mu2)*power3(
     invdmstau1mu) + 288*logmu2Q2*mstau12*power2(mu2)*power3(invdmstau1mu) - 96
     *logmstau12Q2*logmstau22Q2*mstau22*power2(mu2)*power3(invdmstau1mu) + 96*
     logmstau22Q2*logmu2Q2*mstau22*power2(mu2)*power3(invdmstau1mu) + 96*
     mstau12*power2(logmstau12Q2)*power2(mu2)*power3(invdmstau1mu) + 48*mstau22
     *power2(logmstau12Q2)*power2(mu2)*power3(invdmstau1mu) - 96*mstau12*power2
     (logmu2Q2)*power2(mu2)*power3(invdmstau1mu) - 48*mstau22*power2(logmu2Q2)*
     power2(mu2)*power3(invdmstau1mu) + 96*logmstau12Q2*mh2*power2(mu2)*power2(
     sa)*power3(invdmstau1mu) - 96*logmh2Q2*logmstau12Q2*mh2*power2(mu2)*power2
     (sa)*power3(invdmstau1mu) - 96*logmu2Q2*mh2*power2(mu2)*power2(sa)*power3(
     invdmstau1mu) + 96*logmh2Q2*logmu2Q2*mh2*power2(mu2)*power2(sa)*power3(
     invdmstau1mu) - 96*logmstau12Q2*mH2*power2(mu2)*power2(sa)*power3(
     invdmstau1mu) + 96*logmH2Q2*logmstau12Q2*mH2*power2(mu2)*power2(sa)*power3
     (invdmstau1mu) + 96*logmu2Q2*mH2*power2(mu2)*power2(sa)*power3(
     invdmstau1mu) - 96*logmH2Q2*logmu2Q2*mH2*power2(mu2)*power2(sa)*power3(
     invdmstau1mu) + 96*logmstau22Q2*mA2*power2(mu2)*power3(invdmstau2mu) - 96*
     logmu2Q2*mA2*power2(mu2)*power3(invdmstau2mu) + 192*logmstau22Q2*mC2*
     power2(mu2)*power3(invdmstau2mu) - 192*logmu2Q2*mC2*power2(mu2)*power3(
     invdmstau2mu) + 96*logmstau22Q2*mH2*power2(mu2)*power3(invdmstau2mu) - 96*
     logmH2Q2*logmstau22Q2*mH2*power2(mu2)*power3(invdmstau2mu) - 96*logmu2Q2*
     mH2*power2(mu2)*power3(invdmstau2mu) + 96*logmH2Q2*logmu2Q2*mH2*power2(mu2
     )*power3(invdmstau2mu) + 192*logmstau22Q2*msntau2*power2(mu2)*power3(
     invdmstau2mu) - 192*logmsntau2Q2*logmstau22Q2*msntau2*power2(mu2)*power3(
     invdmstau2mu) - 192*logmu2Q2*msntau2*power2(mu2)*power3(invdmstau2mu) +
     192*logmsntau2Q2*logmu2Q2*msntau2*power2(mu2)*power3(invdmstau2mu) - 96*
     logmstau12Q2*logmstau22Q2*mstau12*power2(mu2)*power3(invdmstau2mu) + 96*
     logmstau12Q2*logmu2Q2*mstau12*power2(mu2)*power3(invdmstau2mu) - 576*
     logmstau22Q2*mstau22*power2(mu2)*power3(invdmstau2mu) + 576*logmu2Q2*
     mstau22*power2(mu2)*power3(invdmstau2mu) + 48*mstau12*power2(logmstau22Q2)
     *power2(mu2)*power3(invdmstau2mu) + 192*mstau22*power2(logmstau22Q2)*
     power2(mu2)*power3(invdmstau2mu) - 48*mstau12*power2(logmu2Q2)*power2(mu2)
     *power3(invdmstau2mu) - 192*mstau22*power2(logmu2Q2)*power2(mu2)*power3(
     invdmstau2mu) + 96*logmstau22Q2*mh2*power2(mu2)*power2(sa)*power3(
     invdmstau2mu) - 96*logmh2Q2*logmstau22Q2*mh2*power2(mu2)*power2(sa)*power3
     (invdmstau2mu) - 96*logmu2Q2*mh2*power2(mu2)*power2(sa)*power3(
     invdmstau2mu) + 96*logmh2Q2*logmu2Q2*mh2*power2(mu2)*power2(sa)*power3(
     invdmstau2mu) - 96*logmstau22Q2*mH2*power2(mu2)*power2(sa)*power3(
     invdmstau2mu) + 96*logmH2Q2*logmstau22Q2*mH2*power2(mu2)*power2(sa)*power3
     (invdmstau2mu) + 96*logmu2Q2*mH2*power2(mu2)*power2(sa)*power3(
     invdmstau2mu) - 96*logmH2Q2*logmu2Q2*mH2*power2(mu2)*power2(sa)*power3(
     invdmstau2mu) - (42*power2(invdmAH)*power3(mA2))/mH2 + (72*logmH2Q2*power2
     (invdmAH)*power3(mA2))/mH2 - (12*power2(invdmAH)*power2(logmH2Q2)*power3(
     mA2))/mH2 - (4*power2(invdmAH)*power2(Pi)*power3(mA2))/mH2 - (42*power2(
     invdmAh)*power2(sa)*power3(mA2))/mh2 + (72*logmh2Q2*power2(invdmAh)*power2
     (sa)*power3(mA2))/mh2 + (42*power2(invdmAH)*power2(sa)*power3(mA2))/mH2 -
     (72*logmH2Q2*power2(invdmAH)*power2(sa)*power3(mA2))/mH2 - (12*power2(
     invdmAh)*power2(logmh2Q2)*power2(sa)*power3(mA2))/mh2 + (12*power2(invdmAH
     )*power2(logmH2Q2)*power2(sa)*power3(mA2))/mH2 - (4*power2(invdmAh)*power2
     (Pi)*power2(sa)*power3(mA2))/mh2 + (4*power2(invdmAH)*power2(Pi)*power2(sa
     )*power3(mA2))/mH2 + (21*power2(invdmAC)*power3(mC2))/mA2 + (21*power2(
     invdmCH)*power3(mC2))/mH2 + (72*logmH2Q2*power2(invdmCH)*power3(mC2))/mH2
     - (12*power2(invdmCH)*power2(logmH2Q2)*power3(mC2))/mH2 - (4*power2(
     invdmAC)*power2(Pi)*power3(mC2))/mA2 - (4*power2(invdmCH)*power2(Pi)*
     power3(mC2))/mH2 + (21*power2(invdmCh)*power2(sa)*power3(mC2))/mh2 + (72*
     logmh2Q2*power2(invdmCh)*power2(sa)*power3(mC2))/mh2 - (21*power2(invdmCH)
     *power2(sa)*power3(mC2))/mH2 - (72*logmH2Q2*power2(invdmCH)*power2(sa)*
     power3(mC2))/mH2 - (12*power2(invdmCh)*power2(logmh2Q2)*power2(sa)*power3(
     mC2))/mh2 + (12*power2(invdmCH)*power2(logmH2Q2)*power2(sa)*power3(mC2))/
     mH2 - (4*power2(invdmCh)*power2(Pi)*power2(sa)*power3(mC2))/mh2 + (4*
     power2(invdmCH)*power2(Pi)*power2(sa)*power3(mC2))/mH2 + 12*power2(
     logmC2Q2)*(7 + invdmAC*mC2 + (12*mC2)/mA2 - (12*mC2)/mH2 - (2*mstau22)/mC2
      + (2*mu2)/mC2 - (4*mu2)/mstau22 + ((mA2 - mC2)*power2(invdmAC)*power2(mC2
     ))/mA2 - (10*power2(mC2))/power2(mA2) + (6*power2(mC2))/power2(mH2) + ((
     mC2 - mH2)*power2(invdmCH)*power2(mC2)*(-1 + power2(sa)))/mH2 + invdmCh*
     mC2*power2(sa) - (12*mC2*power2(sa))/mh2 + (12*mC2*power2(sa))/mH2 +
     power2(invdmCh)*power2(mC2)*power2(sa) + (6*power2(mC2)*power2(sa))/power2
     (mh2) - (6*power2(mC2)*power2(sa))/power2(mH2) + invdmCH*(mC2 - mC2*power2
     (sa)) - (power2(invdmCh)*power2(sa)*power3(mC2))/mh2) + (84*power2(invdmhH
     )*power2(sa)*power3(mH2))/mh2 + (72*logmh2Q2*power2(invdmhH)*power2(sa)*
     power3(mH2))/mh2 - (72*logmH2Q2*power2(invdmhH)*power2(sa)*power3(mH2))/
     mh2 - (24*logmh2Q2*logmH2Q2*power2(invdmhH)*power2(sa)*power3(mH2))/mh2 -
     (12*power2(invdmhH)*power2(logmh2Q2)*power2(sa)*power3(mH2))/mh2 - (12*
     power2(invdmhH)*power2(logmH2Q2)*power2(sa)*power3(mH2))/mh2 - (672*
     invdmhH*power2(sa)*power3(mH2))/power2(mh2) + (528*invdmhH*logmH2Q2*power2
     (sa)*power3(mH2))/power2(mh2) + (336*mstau12*power2(invdmhH)*power2(sa)*
     power3(mH2))/power2(mh2) - (288*logmH2Q2*mstau12*power2(invdmhH)*power2(sa
     )*power3(mH2))/power2(mh2) + (336*mstau22*power2(invdmhH)*power2(sa)*
     power3(mH2))/power2(mh2) - (288*logmH2Q2*mstau22*power2(invdmhH)*power2(sa
     )*power3(mH2))/power2(mh2) - (336*mu2*power2(invdmhH)*power2(sa)*power3(
     mH2))/power2(mh2) + (288*logmH2Q2*mu2*power2(invdmhH)*power2(sa)*power3(
     mH2))/power2(mh2) - (144*invdmhH*power2(logmH2Q2)*power2(sa)*power3(mH2))/
     power2(mh2) + (96*mstau12*power2(invdmhH)*power2(logmH2Q2)*power2(sa)*
     power3(mH2))/power2(mh2) + (96*mstau22*power2(invdmhH)*power2(logmH2Q2)*
     power2(sa)*power3(mH2))/power2(mh2) - (96*mu2*power2(invdmhH)*power2(
     logmH2Q2)*power2(sa)*power3(mH2))/power2(mh2) - (4*power2(invdmhH)*power2(
     Pi)*power2(sa)*power3(mH2))/mh2 - (40*invdmhH*power2(Pi)*power2(sa)*power3
     (mH2))/power2(mh2) + (16*mstau12*power2(invdmhH)*power2(Pi)*power2(sa)*
     power3(mH2))/power2(mh2) + (16*mstau22*power2(invdmhH)*power2(Pi)*power2(
     sa)*power3(mH2))/power2(mh2) - (16*mu2*power2(invdmhH)*power2(Pi)*power2(
     sa)*power3(mH2))/power2(mh2) + (42*power2(invdmsntau2mu)*power3(mstau22))/
     (msntau2 - mstau22) + (144*logmsntau2Q2*power2(invdmsntau2mu)*power3(
     mstau22))/(msntau2 - mstau22) - (72*logmstau22Q2*power2(invdmsntau2mu)*
     power3(mstau22))/(msntau2 - mstau22) - (48*logmsntau2Q2*logmstau22Q2*
     power2(invdmsntau2mu)*power3(mstau22))/(msntau2 - mstau22) - 84*
     invdmstau1mu*power2(invdmstau)*power3(mstau22) + 210*invdmstau2mu*power2(
     invdmstau)*power3(mstau22) + 144*invdmstau1mu*logmstau12Q2*power2(
     invdmstau)*power3(mstau22) - 72*invdmstau2mu*logmstau12Q2*power2(invdmstau
     )*power3(mstau22) - 144*invdmstau2mu*logmstau22Q2*power2(invdmstau)*power3
     (mstau22) - 48*invdmstau1mu*logmstau12Q2*logmstau22Q2*power2(invdmstau)*
     power3(mstau22) + 48*invdmstau2mu*logmstau12Q2*logmstau22Q2*power2(
     invdmstau)*power3(mstau22) + 42*invdmstau*power2(invdmstau1mu)*power3(
     mstau22) + 144*invdmstau*logmstau12Q2*power2(invdmstau1mu)*power3(mstau22)
     - 72*invdmstau*logmstau22Q2*power2(invdmstau1mu)*power3(mstau22) - 48*
     invdmstau*logmstau12Q2*logmstau22Q2*power2(invdmstau1mu)*power3(mstau22) +
     84*mu2*power2(invdmstau)*power2(invdmstau1mu)*power3(mstau22) - 144*
     logmstau12Q2*mu2*power2(invdmstau)*power2(invdmstau1mu)*power3(mstau22) +
     48*logmstau12Q2*logmstau22Q2*mu2*power2(invdmstau)*power2(invdmstau1mu)*
     power3(mstau22) + 282*invdmsntau2mu*power2(invdmstau2mu)*power3(mstau22) +
     282*invdmstau1mu*power2(invdmstau2mu)*power3(mstau22) - 144*invdmsntau2mu*
     logmstau22Q2*power2(invdmstau2mu)*power3(mstau22) - 144*invdmstau1mu*
     logmstau22Q2*power2(invdmstau2mu)*power3(mstau22) - 168*invdmsntau2mu*
     logmu2Q2*power2(invdmstau2mu)*power3(mstau22) - 168*invdmstau1mu*logmu2Q2*
     power2(invdmstau2mu)*power3(mstau22) + 48*invdmsntau2mu*logmstau22Q2*
     logmu2Q2*power2(invdmstau2mu)*power3(mstau22) + 48*invdmstau1mu*
     logmstau22Q2*logmu2Q2*power2(invdmstau2mu)*power3(mstau22) - 462*mu2*
     power2(invdmsntau2mu)*power2(invdmstau2mu)*power3(mstau22) + 360*
     logmstau22Q2*mu2*power2(invdmsntau2mu)*power2(invdmstau2mu)*power3(mstau22
     ) - 48*logmstau22Q2*logmu2Q2*mu2*power2(invdmsntau2mu)*power2(invdmstau2mu
     )*power3(mstau22) - 210*mu2*power2(invdmstau)*power2(invdmstau2mu)*power3(
     mstau22) + 72*logmstau12Q2*mu2*power2(invdmstau)*power2(invdmstau2mu)*
     power3(mstau22) + 144*logmstau22Q2*mu2*power2(invdmstau)*power2(
     invdmstau2mu)*power3(mstau22) - 48*logmstau12Q2*logmstau22Q2*mu2*power2(
     invdmstau)*power2(invdmstau2mu)*power3(mstau22) - 462*mu2*power2(
     invdmstau1mu)*power2(invdmstau2mu)*power3(mstau22) + 360*logmstau22Q2*mu2*
     power2(invdmstau1mu)*power2(invdmstau2mu)*power3(mstau22) - 48*
     logmstau22Q2*logmu2Q2*mu2*power2(invdmstau1mu)*power2(invdmstau2mu)*power3
     (mstau22) - (24*power2(invdmsntau2mu)*power2(logmsntau2Q2)*power3(mstau22)
     )/(msntau2 - mstau22) - 24*invdmstau1mu*power2(invdmstau)*power2(
     logmstau12Q2)*power3(mstau22) + 24*invdmstau2mu*power2(invdmstau)*power2(
     logmstau12Q2)*power3(mstau22) - 24*invdmstau*power2(invdmstau1mu)*power2(
     logmstau12Q2)*power3(mstau22) + 24*mu2*power2(invdmstau)*power2(
     invdmstau1mu)*power2(logmstau12Q2)*power3(mstau22) - 24*mu2*power2(
     invdmstau)*power2(invdmstau2mu)*power2(logmstau12Q2)*power3(mstau22) - (24
     *power2(invdmsntau2mu)*power2(logmstau22Q2)*power3(mstau22))/(msntau2 -
     mstau22) - 24*invdmstau1mu*power2(invdmstau)*power2(logmstau22Q2)*power3(
     mstau22) + 24*invdmstau2mu*power2(invdmstau)*power2(logmstau22Q2)*power3(
     mstau22) - 24*invdmstau*power2(invdmstau1mu)*power2(logmstau22Q2)*power3(
     mstau22) + 24*mu2*power2(invdmstau)*power2(invdmstau1mu)*power2(
     logmstau22Q2)*power3(mstau22) + 24*invdmsntau2mu*power2(invdmstau2mu)*
     power2(logmstau22Q2)*power3(mstau22) + 24*invdmstau1mu*power2(invdmstau2mu
     )*power2(logmstau22Q2)*power3(mstau22) - 24*mu2*power2(invdmsntau2mu)*
     power2(invdmstau2mu)*power2(logmstau22Q2)*power3(mstau22) - 24*mu2*power2(
     invdmstau)*power2(invdmstau2mu)*power2(logmstau22Q2)*power3(mstau22) - 24*
     mu2*power2(invdmstau1mu)*power2(invdmstau2mu)*power2(logmstau22Q2)*power3(
     mstau22) + 24*invdmsntau2mu*power2(invdmstau2mu)*power2(logmu2Q2)*power3(
     mstau22) + 24*invdmstau1mu*power2(invdmstau2mu)*power2(logmu2Q2)*power3(
     mstau22) - 24*mu2*power2(invdmsntau2mu)*power2(invdmstau2mu)*power2(
     logmu2Q2)*power3(mstau22) - 24*mu2*power2(invdmstau1mu)*power2(
     invdmstau2mu)*power2(logmu2Q2)*power3(mstau22) - (84*invdmsntau2mu*power3(
     mstau22))/power2(msntau2 - mstau22) + (144*invdmsntau2mu*logmsntau2Q2*
     power3(mstau22))/power2(msntau2 - mstau22) - (48*invdmsntau2mu*
     logmsntau2Q2*logmstau22Q2*power3(mstau22))/power2(msntau2 - mstau22) + (84
     *mu2*power2(invdmsntau2mu)*power3(mstau22))/power2(msntau2 - mstau22) - (
     144*logmsntau2Q2*mu2*power2(invdmsntau2mu)*power3(mstau22))/power2(msntau2
      - mstau22) + (48*logmsntau2Q2*logmstau22Q2*mu2*power2(invdmsntau2mu)*
     power3(mstau22))/power2(msntau2 - mstau22) - (24*invdmsntau2mu*power2(
     logmsntau2Q2)*power3(mstau22))/power2(msntau2 - mstau22) + (24*mu2*power2(
     invdmsntau2mu)*power2(logmsntau2Q2)*power3(mstau22))/power2(msntau2 -
     mstau22) - (24*invdmsntau2mu*power2(logmstau22Q2)*power3(mstau22))/power2(
     msntau2 - mstau22) + (24*mu2*power2(invdmsntau2mu)*power2(logmstau22Q2)*
     power3(mstau22))/power2(msntau2 - mstau22) - (8*power2(invdmsntau2mu)*
     power2(Pi)*power3(mstau22))/(msntau2 - mstau22) - 8*invdmstau1mu*power2(
     invdmstau)*power2(Pi)*power3(mstau22) + 8*invdmstau2mu*power2(invdmstau)*
     power2(Pi)*power3(mstau22) - 8*invdmstau*power2(invdmstau1mu)*power2(Pi)*
     power3(mstau22) + 8*mu2*power2(invdmstau)*power2(invdmstau1mu)*power2(Pi)*
     power3(mstau22) + 8*invdmsntau2mu*power2(invdmstau2mu)*power2(Pi)*power3(
     mstau22) + 8*invdmstau1mu*power2(invdmstau2mu)*power2(Pi)*power3(mstau22)
     - 8*mu2*power2(invdmsntau2mu)*power2(invdmstau2mu)*power2(Pi)*power3(
     mstau22) - 8*mu2*power2(invdmstau)*power2(invdmstau2mu)*power2(Pi)*power3(
     mstau22) - 8*mu2*power2(invdmstau1mu)*power2(invdmstau2mu)*power2(Pi)*
     power3(mstau22) - (8*invdmsntau2mu*power2(Pi)*power3(mstau22))/power2(
     msntau2 - mstau22) + (8*mu2*power2(invdmsntau2mu)*power2(Pi)*power3(
     mstau22))/power2(msntau2 - mstau22) + 72*invdmstau1mu*logmsntau2Q2*power2(
     invdmsntau2mu)*power3(mu2) + 24*invdmstau2mu*logmsntau2Q2*power2(
     invdmsntau2mu)*power3(mu2) - 72*invdmstau1mu*logmu2Q2*power2(invdmsntau2mu
     )*power3(mu2) - 24*invdmstau2mu*logmu2Q2*power2(invdmsntau2mu)*power3(mu2)
     + (144*power2(invdmsntau2mu)*power3(mu2))/msntau2 - (96*logmsntau2Q2*
     power2(invdmsntau2mu)*power3(mu2))/msntau2 - (672*power2(invdmsntau2mu)*
     power3(mu2))/mstau22 + (288*logmstau22Q2*power2(invdmsntau2mu)*power3(mu2)
     )/mstau22 + (288*logmu2Q2*power2(invdmsntau2mu)*power3(mu2))/mstau22 - (96
     *logmstau22Q2*logmu2Q2*power2(invdmsntau2mu)*power3(mu2))/mstau22 + 72*
     invdmsntau2mu*logmstau12Q2*power2(invdmstau1mu)*power3(mu2) + 24*
     invdmstau2mu*logmstau12Q2*power2(invdmstau1mu)*power3(mu2) - 72*
     invdmsntau2mu*logmu2Q2*power2(invdmstau1mu)*power3(mu2) - 24*invdmstau2mu*
     logmu2Q2*power2(invdmstau1mu)*power3(mu2) + (144*power2(invdmstau1mu)*
     power3(mu2))/mstau12 - (96*logmstau12Q2*power2(invdmstau1mu)*power3(mu2))/
     mstau12 - (672*power2(invdmstau1mu)*power3(mu2))/mstau22 + (288*
     logmstau22Q2*power2(invdmstau1mu)*power3(mu2))/mstau22 + (288*logmu2Q2*
     power2(invdmstau1mu)*power3(mu2))/mstau22 - (96*logmstau22Q2*logmu2Q2*
     power2(invdmstau1mu)*power3(mu2))/mstau22 + 54*invdmsntau2mu*power2(
     invdmstau2mu)*power3(mu2) - 84*invdmstau*power2(invdmstau2mu)*power3(mu2)
     + 54*invdmstau1mu*power2(invdmstau2mu)*power3(mu2) + 24*invdmsntau2mu*
     logmstau22Q2*power2(invdmstau2mu)*power3(mu2) + 144*invdmstau*logmstau22Q2
     *power2(invdmstau2mu)*power3(mu2) + 24*invdmstau1mu*logmstau22Q2*power2(
     invdmstau2mu)*power3(mu2) - 48*invdmstau*logmstau12Q2*logmstau22Q2*power2(
     invdmstau2mu)*power3(mu2) - (1344*power2(invdmstau2mu)*power3(mu2))/
     msntau2 + (576*logmsntau2Q2*power2(invdmstau2mu)*power3(mu2))/msntau2 + (
     288*logmstau22Q2*power2(invdmstau2mu)*power3(mu2))/msntau2 - (96*
     logmsntau2Q2*logmstau22Q2*power2(invdmstau2mu)*power3(mu2))/msntau2 + (288
     *logmu2Q2*power2(invdmstau2mu)*power3(mu2))/msntau2 - (96*logmsntau2Q2*
     logmu2Q2*power2(invdmstau2mu)*power3(mu2))/msntau2 - (1344*power2(
     invdmstau2mu)*power3(mu2))/mstau12 + (576*logmstau12Q2*power2(invdmstau2mu
     )*power3(mu2))/mstau12 + (288*logmstau22Q2*power2(invdmstau2mu)*power3(mu2
     ))/mstau12 - (96*logmstau12Q2*logmstau22Q2*power2(invdmstau2mu)*power3(mu2
     ))/mstau12 + (288*logmu2Q2*power2(invdmstau2mu)*power3(mu2))/mstau12 - (96
     *logmstau12Q2*logmu2Q2*power2(invdmstau2mu)*power3(mu2))/mstau12 + (288*
     power2(invdmstau2mu)*power3(mu2))/mstau22 - (192*logmstau22Q2*power2(
     invdmstau2mu)*power3(mu2))/mstau22 + 294*mstau22*power2(invdmsntau2mu)*
     power2(invdmstau2mu)*power3(mu2) + 216*logmstau22Q2*mstau22*power2(
     invdmsntau2mu)*power2(invdmstau2mu)*power3(mu2) - 288*logmu2Q2*mstau22*
     power2(invdmsntau2mu)*power2(invdmstau2mu)*power3(mu2) - 48*logmstau22Q2*
     logmu2Q2*mstau22*power2(invdmsntau2mu)*power2(invdmstau2mu)*power3(mu2) -
     84*mstau22*power2(invdmstau)*power2(invdmstau2mu)*power3(mu2) + 144*
     logmstau22Q2*mstau22*power2(invdmstau)*power2(invdmstau2mu)*power3(mu2) -
     48*logmstau12Q2*logmstau22Q2*mstau22*power2(invdmstau)*power2(invdmstau2mu
     )*power3(mu2) + 294*mstau22*power2(invdmstau1mu)*power2(invdmstau2mu)*
     power3(mu2) + 216*logmstau22Q2*mstau22*power2(invdmstau1mu)*power2(
     invdmstau2mu)*power3(mu2) - 288*logmu2Q2*mstau22*power2(invdmstau1mu)*
     power2(invdmstau2mu)*power3(mu2) - 48*logmstau22Q2*logmu2Q2*mstau22*power2
     (invdmstau1mu)*power2(invdmstau2mu)*power3(mu2) - (96*power2(invdmstau2mu)
     *power2(logmsntau2Q2)*power3(mu2))/msntau2 - 24*invdmstau*power2(
     invdmstau2mu)*power2(logmstau12Q2)*power3(mu2) - (96*power2(invdmstau2mu)*
     power2(logmstau12Q2)*power3(mu2))/mstau12 - 24*mstau22*power2(invdmstau)*
     power2(invdmstau2mu)*power2(logmstau12Q2)*power3(mu2) - (48*power2(
     invdmsntau2mu)*power2(logmstau22Q2)*power3(mu2))/mstau22 - (48*power2(
     invdmstau1mu)*power2(logmstau22Q2)*power3(mu2))/mstau22 - 24*invdmstau*
     power2(invdmstau2mu)*power2(logmstau22Q2)*power3(mu2) - (48*power2(
     invdmstau2mu)*power2(logmstau22Q2)*power3(mu2))/msntau2 - (48*power2(
     invdmstau2mu)*power2(logmstau22Q2)*power3(mu2))/mstau12 - 24*mstau22*
     power2(invdmsntau2mu)*power2(invdmstau2mu)*power2(logmstau22Q2)*power3(mu2
     ) - 24*mstau22*power2(invdmstau)*power2(invdmstau2mu)*power2(logmstau22Q2)
     *power3(mu2) - 24*mstau22*power2(invdmstau1mu)*power2(invdmstau2mu)*power2
     (logmstau22Q2)*power3(mu2) - (48*power2(invdmsntau2mu)*power2(logmu2Q2)*
     power3(mu2))/mstau22 - (48*power2(invdmstau1mu)*power2(logmu2Q2)*power3(
     mu2))/mstau22 - (48*power2(invdmstau2mu)*power2(logmu2Q2)*power3(mu2))/
     msntau2 - (48*power2(invdmstau2mu)*power2(logmu2Q2)*power3(mu2))/mstau12 -
     24*mstau22*power2(invdmsntau2mu)*power2(invdmstau2mu)*power2(logmu2Q2)*
     power3(mu2) - 24*mstau22*power2(invdmstau1mu)*power2(invdmstau2mu)*power2(
     logmu2Q2)*power3(mu2) - (16*power2(invdmsntau2mu)*power2(Pi)*power3(mu2))/
     mstau22 - (16*power2(invdmstau1mu)*power2(Pi)*power3(mu2))/mstau22 - 8*
     invdmstau*power2(invdmstau2mu)*power2(Pi)*power3(mu2) - (32*power2(
     invdmstau2mu)*power2(Pi)*power3(mu2))/msntau2 - (32*power2(invdmstau2mu)*
     power2(Pi)*power3(mu2))/mstau12 - 8*mstau22*power2(invdmsntau2mu)*power2(
     invdmstau2mu)*power2(Pi)*power3(mu2) - 8*mstau22*power2(invdmstau)*power2(
     invdmstau2mu)*power2(Pi)*power3(mu2) - 8*mstau22*power2(invdmstau1mu)*
     power2(invdmstau2mu)*power2(Pi)*power3(mu2) - 24*logmsntau2Q2*power3(
     invdmsntau2mu)*power3(mu2) + 24*logmu2Q2*power3(invdmsntau2mu)*power3(mu2)
     + 96*logmsntau2Q2*logmu2Q2*power3(invdmsntau2mu)*power3(mu2) + 48*power2(
     logmsntau2Q2)*power3(invdmsntau2mu)*power3(mu2) - 144*power2(logmu2Q2)*
     power3(invdmsntau2mu)*power3(mu2) - 24*logmstau12Q2*power3(invdmstau1mu)*
     power3(mu2) + 24*logmu2Q2*power3(invdmstau1mu)*power3(mu2) + 96*
     logmstau12Q2*logmu2Q2*power3(invdmstau1mu)*power3(mu2) + 48*power2(
     logmstau12Q2)*power3(invdmstau1mu)*power3(mu2) - 144*power2(logmu2Q2)*
     power3(invdmstau1mu)*power3(mu2) - 120*logmstau22Q2*power3(invdmstau2mu)*
     power3(mu2) + 120*logmu2Q2*power3(invdmstau2mu)*power3(mu2) + 288*
     logmstau22Q2*logmu2Q2*power3(invdmstau2mu)*power3(mu2) + 48*power2(
     logmstau22Q2)*power3(invdmstau2mu)*power3(mu2) - 336*power2(logmu2Q2)*
     power3(invdmstau2mu)*power3(mu2) - 720*invdmsntau2mu*power2(mu2)*power3(
     mstau22)*power4(invdmstau2mu) - 720*invdmstau1mu*power2(mu2)*power3(
     mstau22)*power4(invdmstau2mu) + 960*invdmsntau2mu*logmu2Q2*power2(mu2)*
     power3(mstau22)*power4(invdmstau2mu) + 960*invdmstau1mu*logmu2Q2*power2(
     mu2)*power3(mstau22)*power4(invdmstau2mu) + 552*invdmsntau2mu*mstau24*
     power3(mu2)*power4(invdmstau2mu) + 552*invdmstau1mu*mstau24*power3(mu2)*
     power4(invdmstau2mu) + 144*invdmsntau2mu*logmstau22Q2*mstau24*power3(mu2)*
     power4(invdmstau2mu) + 144*invdmstau1mu*logmstau22Q2*mstau24*power3(mu2)*
     power4(invdmstau2mu) - 816*invdmsntau2mu*logmu2Q2*mstau24*power3(mu2)*
     power4(invdmstau2mu) - 816*invdmstau1mu*logmu2Q2*mstau24*power3(mu2)*
     power4(invdmstau2mu) - 96*invdmsntau2mu*logmstau22Q2*logmu2Q2*mstau24*
     power3(mu2)*power4(invdmstau2mu) - 96*invdmstau1mu*logmstau22Q2*logmu2Q2*
     mstau24*power3(mu2)*power4(invdmstau2mu) - 48*invdmsntau2mu*mstau24*power2
     (logmstau22Q2)*power3(mu2)*power4(invdmstau2mu) - 48*invdmstau1mu*mstau24*
     power2(logmstau22Q2)*power3(mu2)*power4(invdmstau2mu) - 48*invdmsntau2mu*
     mstau24*power2(logmu2Q2)*power3(mu2)*power4(invdmstau2mu) - 48*
     invdmstau1mu*mstau24*power2(logmu2Q2)*power3(mu2)*power4(invdmstau2mu) -
     16*invdmsntau2mu*mstau24*power2(Pi)*power3(mu2)*power4(invdmstau2mu) - 16*
     invdmstau1mu*mstau24*power2(Pi)*power3(mu2)*power4(invdmstau2mu) - (63*
     power2(invdmAH)*power4(mA2))/power2(mH2) - (63*power2(invdmAh)*power2(sa)*
     power4(mA2))/power2(mh2) + (63*power2(invdmAH)*power2(sa)*power4(mA2))/
     power2(mH2) + (12*power2(logmA2Q2)*(mA2*(-2*mstau12*mstau22*mu2 + msntau2*
     (-2*mstau12*mu2 - 2*mstau22*mu2 + mstau12*mstau22*(4 + invdmAC*mC2 +
     power2(invdmAC)*power2(mC2))))*power2(mh2)*power2(mH2) + mh2*mH2*msntau2*
     mstau12*mstau22*power2(mA2)*(16*mH2*power2(sa) + mh2*(16 + invdmAH*mH2*(-1
      + power2(sa)) - (16 + invdmAh*mH2)*power2(sa))) + msntau2*mstau12*mstau22
     *(6*power2(mH2)*power2(sa) + power2(mh2)*(6 - power2(invdmAH)*power2(mH2)*
     (-1 + power2(sa)) + (-6 + power2(invdmAh)*power2(mH2))*power2(sa)))*power3
     (mA2) - msntau2*mstau12*mstau22*power2(mh2)*power2(mH2)*(4*mC2 + 5*(mH2 +
     msntau2 + mstau12 + mstau22 - 3*mu2 + mh2*power2(sa) - mH2*power2(sa)) +
     power2(invdmAC)*power3(mC2)) + mh2*mH2*msntau2*mstau12*mstau22*(mh2*power2
     (invdmAH)*(-1 + power2(sa)) - mH2*power2(invdmAh)*power2(sa))*power4(mA2))
     )/(mA2*msntau2*mstau12*mstau22*power2(mh2)*power2(mH2)) - (63*power2(
     invdmAC)*power4(mC2))/power2(mA2) - (63*power2(invdmCH)*power4(mC2))/
     power2(mH2) - (63*power2(invdmCh)*power2(sa)*power4(mC2))/power2(mh2) + (
     63*power2(invdmCH)*power2(sa)*power4(mC2))/power2(mH2) - (63*power2(
     invdmhH)*power2(sa)*power4(mH2))/power2(mh2) + (36*logmH2Q2*power2(invdmhH
     )*power2(sa)*power4(mH2))/power2(mh2) - 84*power2(invdmstau)*power2(
     invdmstau1mu)*power4(mstau22) + 144*logmstau12Q2*power2(invdmstau)*power2(
     invdmstau1mu)*power4(mstau22) - 48*logmstau12Q2*logmstau22Q2*power2(
     invdmstau)*power2(invdmstau1mu)*power4(mstau22) + 126*power2(invdmsntau2mu
     )*power2(invdmstau2mu)*power4(mstau22) - 72*logmstau22Q2*power2(
     invdmsntau2mu)*power2(invdmstau2mu)*power4(mstau22) + 126*power2(
     invdmstau1mu)*power2(invdmstau2mu)*power4(mstau22) - 72*logmstau22Q2*
     power2(invdmstau1mu)*power2(invdmstau2mu)*power4(mstau22) - 24*power2(
     invdmstau)*power2(invdmstau1mu)*power2(logmstau12Q2)*power4(mstau22) - 24*
     power2(invdmstau)*power2(invdmstau1mu)*power2(logmstau22Q2)*power4(mstau22
     ) - (84*power2(invdmsntau2mu)*power4(mstau22))/power2(msntau2 - mstau22) +
     (144*logmsntau2Q2*power2(invdmsntau2mu)*power4(mstau22))/power2(msntau2 -
     mstau22) - (48*logmsntau2Q2*logmstau22Q2*power2(invdmsntau2mu)*power4(
     mstau22))/power2(msntau2 - mstau22) - (24*power2(invdmsntau2mu)*power2(
     logmsntau2Q2)*power4(mstau22))/power2(msntau2 - mstau22) - (24*power2(
     invdmsntau2mu)*power2(logmstau22Q2)*power4(mstau22))/power2(msntau2 -
     mstau22) - 8*power2(invdmstau)*power2(invdmstau1mu)*power2(Pi)*power4(
     mstau22) - (8*power2(invdmsntau2mu)*power2(Pi)*power4(mstau22))/power2(
     msntau2 - mstau22) + 528*invdmsntau2mu*mu2*power4(invdmstau2mu)*power4(
     mstau22) + 528*invdmstau1mu*mu2*power4(invdmstau2mu)*power4(mstau22) - 144
     *invdmsntau2mu*logmstau22Q2*mu2*power4(invdmstau2mu)*power4(mstau22) - 144
     *invdmstau1mu*logmstau22Q2*mu2*power4(invdmstau2mu)*power4(mstau22) - 624*
     invdmsntau2mu*logmu2Q2*mu2*power4(invdmstau2mu)*power4(mstau22) - 624*
     invdmstau1mu*logmu2Q2*mu2*power4(invdmstau2mu)*power4(mstau22) + 96*
     invdmsntau2mu*logmstau22Q2*logmu2Q2*mu2*power4(invdmstau2mu)*power4(
     mstau22) + 96*invdmstau1mu*logmstau22Q2*logmu2Q2*mu2*power4(invdmstau2mu)*
     power4(mstau22) + 48*invdmsntau2mu*mu2*power2(logmstau22Q2)*power4(
     invdmstau2mu)*power4(mstau22) + 48*invdmstau1mu*mu2*power2(logmstau22Q2)*
     power4(invdmstau2mu)*power4(mstau22) + 48*invdmsntau2mu*mu2*power2(
     logmu2Q2)*power4(invdmstau2mu)*power4(mstau22) + 48*invdmstau1mu*mu2*
     power2(logmu2Q2)*power4(invdmstau2mu)*power4(mstau22) + 16*invdmsntau2mu*
     mu2*power2(Pi)*power4(invdmstau2mu)*power4(mstau22) + 16*invdmstau1mu*mu2*
     power2(Pi)*power4(invdmstau2mu)*power4(mstau22) + 72*logmsntau2Q2*
     logmstau12Q2*power2(invdmsntau2mu)*power2(invdmstau1mu)*power4(mu2) - 72*
     logmsntau2Q2*logmu2Q2*power2(invdmsntau2mu)*power2(invdmstau1mu)*power4(
     mu2) - 72*logmstau12Q2*logmu2Q2*power2(invdmsntau2mu)*power2(invdmstau1mu)
     *power4(mu2) - 252*power2(invdmsntau2mu)*power2(invdmstau2mu)*power4(mu2)
     + 24*logmsntau2Q2*logmstau22Q2*power2(invdmsntau2mu)*power2(invdmstau2mu)*
     power4(mu2) + 144*logmu2Q2*power2(invdmsntau2mu)*power2(invdmstau2mu)*
     power4(mu2) - 24*logmsntau2Q2*logmu2Q2*power2(invdmsntau2mu)*power2(
     invdmstau2mu)*power4(mu2) - 24*logmstau22Q2*logmu2Q2*power2(invdmsntau2mu)
     *power2(invdmstau2mu)*power4(mu2) - 252*power2(invdmstau1mu)*power2(
     invdmstau2mu)*power4(mu2) + 24*logmstau12Q2*logmstau22Q2*power2(
     invdmstau1mu)*power2(invdmstau2mu)*power4(mu2) + 144*logmu2Q2*power2(
     invdmstau1mu)*power2(invdmstau2mu)*power4(mu2) - 24*logmstau12Q2*logmu2Q2*
     power2(invdmstau1mu)*power2(invdmstau2mu)*power4(mu2) - 24*logmstau22Q2*
     logmu2Q2*power2(invdmstau1mu)*power2(invdmstau2mu)*power4(mu2) + 72*power2
     (invdmsntau2mu)*power2(invdmstau1mu)*power2(logmu2Q2)*power4(mu2) + 24*
     power2(invdmsntau2mu)*power2(invdmstau2mu)*power2(logmu2Q2)*power4(mu2) +
     24*power2(invdmstau1mu)*power2(invdmstau2mu)*power2(logmu2Q2)*power4(mu2)
     - 72*logmsntau2Q2*logmu2Q2*power4(invdmsntau2mu)*power4(mu2) + 36*power2(
     logmsntau2Q2)*power4(invdmsntau2mu)*power4(mu2) + 36*power2(logmu2Q2)*
     power4(invdmsntau2mu)*power4(mu2) - 72*logmstau12Q2*logmu2Q2*power4(
     invdmstau1mu)*power4(mu2) + 36*power2(logmstau12Q2)*power4(invdmstau1mu)*
     power4(mu2) + 36*power2(logmu2Q2)*power4(invdmstau1mu)*power4(mu2) - 72*
     logmstau22Q2*logmu2Q2*power4(invdmstau2mu)*power4(mu2) - 276*invdmsntau2mu
     *mstau22*power4(invdmstau2mu)*power4(mu2) - 276*invdmstau1mu*mstau22*
     power4(invdmstau2mu)*power4(mu2) - 72*invdmsntau2mu*logmstau22Q2*mstau22*
     power4(invdmstau2mu)*power4(mu2) - 72*invdmstau1mu*logmstau22Q2*mstau22*
     power4(invdmstau2mu)*power4(mu2) + 408*invdmsntau2mu*logmu2Q2*mstau22*
     power4(invdmstau2mu)*power4(mu2) + 408*invdmstau1mu*logmu2Q2*mstau22*
     power4(invdmstau2mu)*power4(mu2) + 48*invdmsntau2mu*logmstau22Q2*logmu2Q2*
     mstau22*power4(invdmstau2mu)*power4(mu2) + 48*invdmstau1mu*logmstau22Q2*
     logmu2Q2*mstau22*power4(invdmstau2mu)*power4(mu2) + 36*power2(logmstau22Q2
     )*power4(invdmstau2mu)*power4(mu2) + 24*invdmsntau2mu*mstau22*power2(
     logmstau22Q2)*power4(invdmstau2mu)*power4(mu2) + 24*invdmstau1mu*mstau22*
     power2(logmstau22Q2)*power4(invdmstau2mu)*power4(mu2) + 36*power2(logmu2Q2
     )*power4(invdmstau2mu)*power4(mu2) + 24*invdmsntau2mu*mstau22*power2(
     logmu2Q2)*power4(invdmstau2mu)*power4(mu2) + 24*invdmstau1mu*mstau22*
     power2(logmu2Q2)*power4(invdmstau2mu)*power4(mu2) + 8*invdmsntau2mu*
     mstau22*power2(Pi)*power4(invdmstau2mu)*power4(mu2) + 8*invdmstau1mu*
     mstau22*power2(Pi)*power4(invdmstau2mu)*power4(mu2) - 792*power4(sa) + 648
     *logmh2Q2*power4(sa) - 504*logmH2Q2*power4(sa) - 576*logmh2Q2*logmH2Q2*
     power4(sa) + (72*mh2*power4(sa))/mH2 - (144*logmh2Q2*mh2*power4(sa))/mH2 +
     (144*logmH2Q2*mh2*power4(sa))/mH2 - (72*logmh2Q2*logmH2Q2*mh2*power4(sa))/
     mH2 + 1758*invdmhH*mH2*power4(sa) - 1080*invdmhH*logmh2Q2*mH2*power4(sa) -
     384*invdmhH*logmH2Q2*mH2*power4(sa) - 24*invdmhH*logmh2Q2*logmH2Q2*mH2*
     power4(sa) + (1830*mH2*power4(sa))/mh2 + (216*logmh2Q2*mH2*power4(sa))/mh2
      - (1680*logmH2Q2*mH2*power4(sa))/mh2 - (96*logmh2Q2*logmH2Q2*mH2*power4(
     sa))/mh2 + 1008*invdmhH*mstau12*power4(sa) - 480*invdmhH*logmH2Q2*mstau12*
     power4(sa) - 576*invdmhH*logmstau12Q2*mstau12*power4(sa) + (1008*mstau12*
     power4(sa))/mh2 - (480*logmH2Q2*mstau12*power4(sa))/mh2 - (576*
     logmstau12Q2*mstau12*power4(sa))/mh2 - (2304*invdmhH*mH2*mstau12*power4(sa
     ))/mh2 + (1536*invdmhH*logmH2Q2*mH2*mstau12*power4(sa))/mh2 + (768*invdmhH
     *logmstau12Q2*mH2*mstau12*power4(sa))/mh2 + 1008*invdmhH*mstau22*power4(sa
     ) - 480*invdmhH*logmH2Q2*mstau22*power4(sa) - 576*invdmhH*logmstau22Q2*
     mstau22*power4(sa) + (1008*mstau22*power4(sa))/mh2 - (480*logmH2Q2*mstau22
     *power4(sa))/mh2 - (576*logmstau22Q2*mstau22*power4(sa))/mh2 - (2304*
     invdmhH*mH2*mstau22*power4(sa))/mh2 + (1536*invdmhH*logmH2Q2*mH2*mstau22*
     power4(sa))/mh2 + (768*invdmhH*logmstau22Q2*mH2*mstau22*power4(sa))/mh2 +
     648*invdmhH*mu2*power4(sa) - 1152*invdmhH*logmh2Q2*mu2*power4(sa) + 480*
     invdmhH*logmH2Q2*mu2*power4(sa) - 480*invdmhH*logmstau12Q2*mu2*power4(sa)
     + 576*invdmhH*logmh2Q2*logmstau12Q2*mu2*power4(sa) - (1080*mu2*power4(sa))
     /mh2 + (576*logmh2Q2*mu2*power4(sa))/mh2 - (672*logmH2Q2*mu2*power4(sa))/
     mh2 + (672*logmstau12Q2*mu2*power4(sa))/mh2 - (288*logmh2Q2*logmstau12Q2*
     mu2*power4(sa))/mh2 + (576*logmH2Q2*logmstau12Q2*mu2*power4(sa))/mh2 - (
     1728*mu2*power4(sa))/mH2 + (576*logmH2Q2*mu2*power4(sa))/mH2 + (1152*
     logmstau12Q2*mu2*power4(sa))/mH2 - (288*logmH2Q2*logmstau12Q2*mu2*power4(
     sa))/mH2 + (768*invdmhH*logmstau12Q2*mH2*mu2*power4(sa))/mh2 - (768*
     invdmhH*logmH2Q2*logmstau12Q2*mH2*mu2*power4(sa))/mh2 - (4032*mu2*power4(
     sa))/mstau12 + (864*logmh2Q2*mu2*power4(sa))/mstau12 + (864*logmH2Q2*mu2*
     power4(sa))/mstau12 + (1728*logmstau12Q2*mu2*power4(sa))/mstau12 - (288*
     logmh2Q2*logmstau12Q2*mu2*power4(sa))/mstau12 - (288*logmH2Q2*logmstau12Q2
     *mu2*power4(sa))/mstau12 + (168*invdmhH*mstau12*mu2*power4(sa))/mh2 - (96*
     invdmhH*logmstau12Q2*mstau12*mu2*power4(sa))/mh2 - 480*invdmhH*invdmstau*
     logmstau12Q2*mstau22*mu2*power4(sa) + 576*invdmhH*invdmstau*logmh2Q2*
     logmstau12Q2*mstau22*mu2*power4(sa) + 480*invdmhH*invdmstau*logmstau22Q2*
     mstau22*mu2*power4(sa) - 576*invdmhH*invdmstau*logmh2Q2*logmstau22Q2*
     mstau22*mu2*power4(sa) + (168*invdmhH*mstau22*mu2*power4(sa))/mh2 - (192*
     invdmstau*logmstau12Q2*mstau22*mu2*power4(sa))/mh2 + (576*invdmstau*
     logmH2Q2*logmstau12Q2*mstau22*mu2*power4(sa))/mh2 - (96*invdmhH*
     logmstau22Q2*mstau22*mu2*power4(sa))/mh2 + (192*invdmstau*logmstau22Q2*
     mstau22*mu2*power4(sa))/mh2 - (576*invdmstau*logmH2Q2*logmstau22Q2*mstau22
     *mu2*power4(sa))/mh2 + (288*invdmstau*logmstau12Q2*mstau22*mu2*power4(sa))
     /mH2 - (288*invdmstau*logmstau22Q2*mstau22*mu2*power4(sa))/mH2 + (768*
     invdmhH*invdmstau*logmstau12Q2*mH2*mstau22*mu2*power4(sa))/mh2 - (768*
     invdmhH*invdmstau*logmH2Q2*logmstau12Q2*mH2*mstau22*mu2*power4(sa))/mh2 -
     (768*invdmhH*invdmstau*logmstau22Q2*mH2*mstau22*mu2*power4(sa))/mh2 + (768
     *invdmhH*invdmstau*logmH2Q2*logmstau22Q2*mH2*mstau22*mu2*power4(sa))/mh2 -
     (4032*mstau22*mu2*power4(sa))/(mh2*mstau12) + (2592*logmstau12Q2*mstau22*
     mu2*power4(sa))/(mh2*mstau12) - (288*logmh2Q2*logmstau12Q2*mstau22*mu2*
     power4(sa))/(mh2*mstau12) + (864*logmstau22Q2*mstau22*mu2*power4(sa))/(mh2
     *mstau12) + (288*logmh2Q2*logmstau22Q2*mstau22*mu2*power4(sa))/(mh2*
     mstau12) - (576*logmstau12Q2*logmstau22Q2*mstau22*mu2*power4(sa))/(mh2*
     mstau12) - (4032*mstau22*mu2*power4(sa))/(mH2*mstau12) + (2592*
     logmstau12Q2*mstau22*mu2*power4(sa))/(mH2*mstau12) - (288*logmH2Q2*
     logmstau12Q2*mstau22*mu2*power4(sa))/(mH2*mstau12) + (864*logmstau22Q2*
     mstau22*mu2*power4(sa))/(mH2*mstau12) + (288*logmH2Q2*logmstau22Q2*mstau22
     *mu2*power4(sa))/(mH2*mstau12) - (576*logmstau12Q2*logmstau22Q2*mstau22*
     mu2*power4(sa))/(mH2*mstau12) - 768*mH2*mstau12*power2(invdmhH)*power4(sa)
     + 768*logmH2Q2*mH2*mstau12*power2(invdmhH)*power4(sa) - 768*mH2*mstau22*
     power2(invdmhH)*power4(sa) + 768*logmH2Q2*mH2*mstau22*power2(invdmhH)*
     power4(sa) + 768*mH2*mu2*power2(invdmhH)*power4(sa) - 768*logmH2Q2*mH2*mu2
     *power2(invdmhH)*power4(sa) + 216*power2(logmh2Q2)*power4(sa) + (108*mh2*
     power2(logmh2Q2)*power4(sa))/mH2 + 276*invdmhH*mH2*power2(logmh2Q2)*power4
     (sa) - (48*mH2*power2(logmh2Q2)*power4(sa))/mh2 + 288*invdmhH*mu2*power2(
     logmh2Q2)*power4(sa) - (144*mu2*power2(logmh2Q2)*power4(sa))/mh2 - (144*
     mu2*power2(logmh2Q2)*power4(sa))/mstau12 + 504*power2(logmH2Q2)*power4(sa)
     - (36*mh2*power2(logmH2Q2)*power4(sa))/mH2 + 180*invdmhH*mH2*power2(
     logmH2Q2)*power4(sa) + (576*mH2*power2(logmH2Q2)*power4(sa))/mh2 + 288*
     invdmhH*mstau12*power2(logmH2Q2)*power4(sa) + (288*mstau12*power2(logmH2Q2
     )*power4(sa))/mh2 - (768*invdmhH*mH2*mstau12*power2(logmH2Q2)*power4(sa))/
     mh2 + 288*invdmhH*mstau22*power2(logmH2Q2)*power4(sa) + (288*mstau22*
     power2(logmH2Q2)*power4(sa))/mh2 - (768*invdmhH*mH2*mstau22*power2(
     logmH2Q2)*power4(sa))/mh2 - 288*invdmhH*mu2*power2(logmH2Q2)*power4(sa) -
     (144*mu2*power2(logmH2Q2)*power4(sa))/mH2 + (384*invdmhH*mH2*mu2*power2(
     logmH2Q2)*power4(sa))/mh2 - (144*mu2*power2(logmH2Q2)*power4(sa))/mstau12
     - 384*mH2*mstau12*power2(invdmhH)*power2(logmH2Q2)*power4(sa) - 384*mH2*
     mstau22*power2(invdmhH)*power2(logmH2Q2)*power4(sa) + 384*mH2*mu2*power2(
     invdmhH)*power2(logmH2Q2)*power4(sa) + 288*invdmhH*mstau12*power2(
     logmstau12Q2)*power4(sa) + (288*mstau12*power2(logmstau12Q2)*power4(sa))/
     mh2 - (384*invdmhH*mH2*mstau12*power2(logmstau12Q2)*power4(sa))/mh2 - (144
     *mu2*power2(logmstau12Q2)*power4(sa))/mh2 - (144*mu2*power2(logmstau12Q2)*
     power4(sa))/mH2 - (288*mu2*power2(logmstau12Q2)*power4(sa))/mstau12 - (432
     *mstau22*mu2*power2(logmstau12Q2)*power4(sa))/(mh2*mstau12) - (432*mstau22
     *mu2*power2(logmstau12Q2)*power4(sa))/(mH2*mstau12) + 288*invdmhH*mstau22*
     power2(logmstau22Q2)*power4(sa) + (288*mstau22*power2(logmstau22Q2)*power4
     (sa))/mh2 - (384*invdmhH*mH2*mstau22*power2(logmstau22Q2)*power4(sa))/mh2
     - (144*mstau22*mu2*power2(logmstau22Q2)*power4(sa))/(mh2*mstau12) - (144*
     mstau22*mu2*power2(logmstau22Q2)*power4(sa))/(mH2*mstau12) - (528*mH2*
     mstau12*power4(sa))/power2(mh2) + (288*logmH2Q2*mH2*mstau12*power4(sa))/
     power2(mh2) + (192*logmstau12Q2*mH2*mstau12*power4(sa))/power2(mh2) - (528
     *mH2*mstau22*power4(sa))/power2(mh2) + (288*logmH2Q2*mH2*mstau22*power4(sa
     ))/power2(mh2) + (192*logmstau22Q2*mH2*mstau22*power4(sa))/power2(mh2) - (
     120*mH2*mu2*power4(sa))/power2(mh2) + (96*logmH2Q2*mH2*mu2*power4(sa))/
     power2(mh2) + (288*logmstau12Q2*mH2*mu2*power4(sa))/power2(mh2) - (192*
     logmH2Q2*logmstau12Q2*mH2*mu2*power4(sa))/power2(mh2) + (168*mstau12*mu2*
     power4(sa))/power2(mh2) - (96*logmstau12Q2*mstau12*mu2*power4(sa))/power2(
     mh2) - (168*invdmhH*mH2*mstau12*mu2*power4(sa))/power2(mh2) + (96*invdmhH*
     logmstau12Q2*mH2*mstau12*mu2*power4(sa))/power2(mh2) + (168*mstau22*mu2*
     power4(sa))/power2(mh2) - (96*logmstau22Q2*mstau22*mu2*power4(sa))/power2(
     mh2) - (168*invdmhH*mH2*mstau22*mu2*power4(sa))/power2(mh2) + (288*
     invdmstau*logmstau12Q2*mH2*mstau22*mu2*power4(sa))/power2(mh2) - (192*
     invdmstau*logmH2Q2*logmstau12Q2*mH2*mstau22*mu2*power4(sa))/power2(mh2) +
     (96*invdmhH*logmstau22Q2*mH2*mstau22*mu2*power4(sa))/power2(mh2) - (288*
     invdmstau*logmstau22Q2*mH2*mstau22*mu2*power4(sa))/power2(mh2) + (192*
     invdmstau*logmH2Q2*logmstau22Q2*mH2*mstau22*mu2*power4(sa))/power2(mh2) -
     (96*mH2*mstau12*power2(logmH2Q2)*power4(sa))/power2(mh2) - (96*mH2*mstau22
     *power2(logmH2Q2)*power4(sa))/power2(mh2) - (96*mH2*mstau12*power2(
     logmstau12Q2)*power4(sa))/power2(mh2) - (96*mH2*mstau22*power2(
     logmstau22Q2)*power4(sa))/power2(mh2) - (252*power2(mh2)*power4(sa))/
     power2(mH2) + (216*logmh2Q2*power2(mh2)*power4(sa))/power2(mH2) - (72*
     power2(logmh2Q2)*power2(mh2)*power4(sa))/power2(mH2) - (2472*invdmhH*
     power2(mH2)*power4(sa))/mh2 + (2064*invdmhH*logmH2Q2*power2(mH2)*power4(sa
     ))/mh2 + 21*power2(invdmhH)*power2(mH2)*power4(sa) + 72*logmh2Q2*power2(
     invdmhH)*power2(mH2)*power4(sa) - 36*logmH2Q2*power2(invdmhH)*power2(mH2)*
     power4(sa) - 24*logmh2Q2*logmH2Q2*power2(invdmhH)*power2(mH2)*power4(sa) +
     (1104*mstau12*power2(invdmhH)*power2(mH2)*power4(sa))/mh2 - (1056*logmH2Q2
     *mstau12*power2(invdmhH)*power2(mH2)*power4(sa))/mh2 + (1104*mstau22*
     power2(invdmhH)*power2(mH2)*power4(sa))/mh2 - (1056*logmH2Q2*mstau22*
     power2(invdmhH)*power2(mH2)*power4(sa))/mh2 - (1104*mu2*power2(invdmhH)*
     power2(mH2)*power4(sa))/mh2 + (1056*logmH2Q2*mu2*power2(invdmhH)*power2(
     mH2)*power4(sa))/mh2 - 12*power2(invdmhH)*power2(logmh2Q2)*power2(mH2)*
     power4(sa) - (624*invdmhH*power2(logmH2Q2)*power2(mH2)*power4(sa))/mh2 -
     12*power2(invdmhH)*power2(logmH2Q2)*power2(mH2)*power4(sa) + (480*mstau12*
     power2(invdmhH)*power2(logmH2Q2)*power2(mH2)*power4(sa))/mh2 + (480*
     mstau22*power2(invdmhH)*power2(logmH2Q2)*power2(mH2)*power4(sa))/mh2 - (
     480*mu2*power2(invdmhH)*power2(logmH2Q2)*power2(mH2)*power4(sa))/mh2 - (
     987*power2(mH2)*power4(sa))/power2(mh2) + (780*logmH2Q2*power2(mH2)*power4
     (sa))/power2(mh2) + (864*invdmhH*mstau12*power2(mH2)*power4(sa))/power2(
     mh2) - (576*invdmhH*logmH2Q2*mstau12*power2(mH2)*power4(sa))/power2(mh2) -
     (192*invdmhH*logmstau12Q2*mstau12*power2(mH2)*power4(sa))/power2(mh2) + (
     864*invdmhH*mstau22*power2(mH2)*power4(sa))/power2(mh2) - (576*invdmhH*
     logmH2Q2*mstau22*power2(mH2)*power4(sa))/power2(mh2) - (192*invdmhH*
     logmstau22Q2*mstau22*power2(mH2)*power4(sa))/power2(mh2) - (216*invdmhH*
     mu2*power2(mH2)*power4(sa))/power2(mh2) + (192*invdmhH*logmH2Q2*mu2*power2
     (mH2)*power4(sa))/power2(mh2) - (288*invdmhH*logmstau12Q2*mu2*power2(mH2)*
     power4(sa))/power2(mh2) + (192*invdmhH*logmH2Q2*logmstau12Q2*mu2*power2(
     mH2)*power4(sa))/power2(mh2) - (288*invdmhH*invdmstau*logmstau12Q2*mstau22
     *mu2*power2(mH2)*power4(sa))/power2(mh2) + (192*invdmhH*invdmstau*logmH2Q2
     *logmstau12Q2*mstau22*mu2*power2(mH2)*power4(sa))/power2(mh2) + (288*
     invdmhH*invdmstau*logmstau22Q2*mstau22*mu2*power2(mH2)*power4(sa))/power2(
     mh2) - (192*invdmhH*invdmstau*logmH2Q2*logmstau22Q2*mstau22*mu2*power2(mH2
     )*power4(sa))/power2(mh2) - (216*power2(logmH2Q2)*power2(mH2)*power4(sa))/
     power2(mh2) + (192*invdmhH*mstau12*power2(logmH2Q2)*power2(mH2)*power4(sa)
     )/power2(mh2) + (192*invdmhH*mstau22*power2(logmH2Q2)*power2(mH2)*power4(
     sa))/power2(mh2) - (96*invdmhH*mu2*power2(logmH2Q2)*power2(mH2)*power4(sa)
     )/power2(mh2) + (96*invdmhH*mstau12*power2(logmstau12Q2)*power2(mH2)*
     power4(sa))/power2(mh2) + (96*invdmhH*mstau22*power2(logmstau22Q2)*power2(
     mH2)*power4(sa))/power2(mh2) - 12*power2(Pi)*power4(sa) + (24*mh2*power2(
     Pi)*power4(sa))/mH2 + 124*invdmhH*mH2*power2(Pi)*power4(sa) + (148*mH2*
     power2(Pi)*power4(sa))/mh2 + 96*invdmhH*mstau12*power2(Pi)*power4(sa) + (
     96*mstau12*power2(Pi)*power4(sa))/mh2 - (192*invdmhH*mH2*mstau12*power2(Pi
     )*power4(sa))/mh2 + 96*invdmhH*mstau22*power2(Pi)*power4(sa) + (96*mstau22
     *power2(Pi)*power4(sa))/mh2 - (192*invdmhH*mH2*mstau22*power2(Pi)*power4(
     sa))/mh2 - (48*mu2*power2(Pi)*power4(sa))/mh2 - (48*mu2*power2(Pi)*power4(
     sa))/mH2 + (64*invdmhH*mH2*mu2*power2(Pi)*power4(sa))/mh2 - (96*mu2*power2
     (Pi)*power4(sa))/mstau12 - (96*mstau22*mu2*power2(Pi)*power4(sa))/(mh2*
     mstau12) - (96*mstau22*mu2*power2(Pi)*power4(sa))/(mH2*mstau12) - 64*mH2*
     mstau12*power2(invdmhH)*power2(Pi)*power4(sa) - 64*mH2*mstau22*power2(
     invdmhH)*power2(Pi)*power4(sa) + 64*mH2*mu2*power2(invdmhH)*power2(Pi)*
     power4(sa) - (32*mH2*mstau12*power2(Pi)*power4(sa))/power2(mh2) - (32*mH2*
     mstau22*power2(Pi)*power4(sa))/power2(mh2) - (18*power2(mh2)*power2(Pi)*
     power4(sa))/power2(mH2) - (168*invdmhH*power2(mH2)*power2(Pi)*power4(sa))/
     mh2 - 4*power2(invdmhH)*power2(mH2)*power2(Pi)*power4(sa) + (80*mstau12*
     power2(invdmhH)*power2(mH2)*power2(Pi)*power4(sa))/mh2 + (80*mstau22*
     power2(invdmhH)*power2(mH2)*power2(Pi)*power4(sa))/mh2 - (80*mu2*power2(
     invdmhH)*power2(mH2)*power2(Pi)*power4(sa))/mh2 - (58*power2(mH2)*power2(
     Pi)*power4(sa))/power2(mh2) + (48*invdmhH*mstau12*power2(mH2)*power2(Pi)*
     power4(sa))/power2(mh2) + (48*invdmhH*mstau22*power2(mH2)*power2(Pi)*
     power4(sa))/power2(mh2) - (16*invdmhH*mu2*power2(mH2)*power2(Pi)*power4(sa
     ))/power2(mh2) - (84*power2(invdmhH)*power3(mH2)*power4(sa))/mh2 - (72*
     logmh2Q2*power2(invdmhH)*power3(mH2)*power4(sa))/mh2 + (72*logmH2Q2*power2
     (invdmhH)*power3(mH2)*power4(sa))/mh2 + (24*logmh2Q2*logmH2Q2*power2(
     invdmhH)*power3(mH2)*power4(sa))/mh2 + (12*power2(invdmhH)*power2(logmh2Q2
     )*power3(mH2)*power4(sa))/mh2 + (12*power2(invdmhH)*power2(logmH2Q2)*
     power3(mH2)*power4(sa))/mh2 + (672*invdmhH*power3(mH2)*power4(sa))/power2(
     mh2) - (528*invdmhH*logmH2Q2*power3(mH2)*power4(sa))/power2(mh2) - (336*
     mstau12*power2(invdmhH)*power3(mH2)*power4(sa))/power2(mh2) + (288*
     logmH2Q2*mstau12*power2(invdmhH)*power3(mH2)*power4(sa))/power2(mh2) - (
     336*mstau22*power2(invdmhH)*power3(mH2)*power4(sa))/power2(mh2) + (288*
     logmH2Q2*mstau22*power2(invdmhH)*power3(mH2)*power4(sa))/power2(mh2) + (
     336*mu2*power2(invdmhH)*power3(mH2)*power4(sa))/power2(mh2) - (288*
     logmH2Q2*mu2*power2(invdmhH)*power3(mH2)*power4(sa))/power2(mh2) + (144*
     invdmhH*power2(logmH2Q2)*power3(mH2)*power4(sa))/power2(mh2) - (96*mstau12
     *power2(invdmhH)*power2(logmH2Q2)*power3(mH2)*power4(sa))/power2(mh2) - (
     96*mstau22*power2(invdmhH)*power2(logmH2Q2)*power3(mH2)*power4(sa))/power2
     (mh2) + (96*mu2*power2(invdmhH)*power2(logmH2Q2)*power3(mH2)*power4(sa))/
     power2(mh2) + (4*power2(invdmhH)*power2(Pi)*power3(mH2)*power4(sa))/mh2 +
     (40*invdmhH*power2(Pi)*power3(mH2)*power4(sa))/power2(mh2) - (16*mstau12*
     power2(invdmhH)*power2(Pi)*power3(mH2)*power4(sa))/power2(mh2) - (16*
     mstau22*power2(invdmhH)*power2(Pi)*power3(mH2)*power4(sa))/power2(mh2) + (
     16*mu2*power2(invdmhH)*power2(Pi)*power3(mH2)*power4(sa))/power2(mh2) + (
     63*power2(invdmhH)*power4(mH2)*power4(sa))/power2(mh2) - (36*logmH2Q2*
     power2(invdmhH)*power4(mH2)*power4(sa))/power2(mh2) - (12*logmA2Q2*(-2*mA2
     *power2(mh2)*power2(mH2)*(2*(3 - 2*logmsntau2Q2 + logmu2Q2)*mstau12*
     mstau22*mu2 + msntau2*(2*(3 - 2*logmstau22Q2 + logmu2Q2)*mstau12*mu2 + 2*(
     3 - 2*logmstau12Q2 + logmu2Q2)*mstau22*mu2 + mstau12*mstau22*(7 + invdmAC*
     (-3 + logmC2Q2)*mC2 + (-3 + logmC2Q2)*power2(invdmAC)*power2(mC2) + 8*
     invdmstau*mu2*(logmstau22Q2 - invdmstau2mu*logmstau22Q2*mu2 + (-
     invdmstau1mu + invdmstau2mu)*logmu2Q2*mu2 + logmstau12Q2*(-1 +
     invdmstau1mu*mu2))*power2(sb)))) + 2*mh2*mH2*msntau2*mstau12*mstau22*
     power2(mA2)*((33 - 4*logmh2Q2)*mH2*power2(sa) + mh2*(33 + 2*invdmstau1mu*
     mH2 + 2*invdmstau2mu*mH2 + 4*mH2*mu2*power2(invdmstau1mu) + 4*mH2*mu2*
     power2(invdmstau2mu) + 3*invdmAH*mH2*(-1 + power2(sa)) - logmH2Q2*(-4 +
     invdmAH*mH2)*(-1 + power2(sa)) - 33*power2(sa) - 3*invdmAh*mH2*power2(sa)
     + invdmAh*logmh2Q2*mH2*power2(sa) + 4*logmstau12Q2*mH2*power2(mu2)*power3(
     invdmstau1mu) - 4*logmu2Q2*mH2*power2(mu2)*power3(invdmstau1mu) + 4*
     logmstau22Q2*mH2*power2(mu2)*power3(invdmstau2mu) - 4*logmu2Q2*mH2*power2(
     mu2)*power3(invdmstau2mu))) + msntau2*mstau12*mstau22*(21*power2(mH2)*
     power2(sa) + power2(mh2)*(21 + (-3 + 2*logmH2Q2)*power2(invdmAH)*power2(
     mH2)*(-1 + power2(sa)) + (-21 + (3 - 2*logmh2Q2)*power2(invdmAh)*power2(
     mH2))*power2(sa)))*power3(mA2) + 2*msntau2*mstau12*mstau22*power2(mh2)*
     power2(mH2)*((-7 + 4*logmC2Q2)*mC2 + 5*((-2 + logmsntau2Q2)*msntau2 - 2*
     mstau12 + logmstau12Q2*mstau12 - 2*mstau22 + logmstau22Q2*mstau22 + 6*mu2
     - 3*logmu2Q2*mu2 - (-2 + logmH2Q2)*mH2*(-1 + power2(sa)) - 2*mh2*power2(sa
     ) + logmh2Q2*mh2*power2(sa)) + (-3 + logmC2Q2)*power2(invdmAC)*power3(mC2)
     ) + 2*mh2*mH2*msntau2*mstau12*mstau22*(-(logmH2Q2*mh2*power2(invdmAH)*(-1
     + power2(sa))) + logmh2Q2*mH2*power2(invdmAh)*power2(sa))*power4(mA2) + 3*
     msntau2*mstau12*mstau22*(power2(invdmAH)*power2(mh2)*(-1 + power2(sa)) -
     power2(invdmAh)*power2(mH2)*power2(sa))*power5(mA2)))/(mA2*msntau2*mstau12
     *mstau22*power2(mh2)*power2(mH2)) - (12*logmC2Q2*(4*mstau22*((-2 +
     logmstau22Q2)*mstau22 - (-2 + logmu2Q2)*mu2)*power2(mh2)*power2(mH2) + 2*
     mC2*power2(mh2)*power2(mH2)*(4*(-3 + 2*logmstau22Q2 - logmu2Q2)*mu2 +
     mstau22*(9 - 3*logmsntau2Q2 - 3*logmstau12Q2 - logmstau22Q2 + 3*
     invdmsntau2mu*mu2 + 3*invdmstau1mu*mu2 + invdmstau2mu*mu2 + 3*logmsntau2Q2
     *power2(invdmsntau2mu)*power2(mu2) - 3*logmu2Q2*power2(invdmsntau2mu)*
     power2(mu2) + 3*logmstau12Q2*power2(invdmstau1mu)*power2(mu2) - 3*logmu2Q2
     *power2(invdmstau1mu)*power2(mu2) + logmstau22Q2*power2(invdmstau2mu)*
     power2(mu2) - logmu2Q2*power2(invdmstau2mu)*power2(mu2) + 4*logmH2Q2*(-1 +
     power2(sa)) - 4*logmh2Q2*power2(sa))) + mh2*mH2*mstau22*power2(mC2)*(-((21
      + 8*logmh2Q2)*mH2*power2(sa)) + mh2*(-21 + 3*invdmAC*mH2 + 3*invdmCH*mH2
     + 8*invdmsntau2mu*mH2 + 8*invdmstau2mu*mH2 + 16*mH2*mu2*power2(
     invdmsntau2mu) + 16*mH2*mu2*power2(invdmstau2mu) + 2*logmH2Q2*(4 + invdmCH
     *mH2)*(-1 + power2(sa)) + 21*power2(sa) + 3*invdmCh*mH2*power2(sa) - 3*
     invdmCH*mH2*power2(sa) - 2*invdmCh*logmh2Q2*mH2*power2(sa) + 16*
     logmsntau2Q2*mH2*power2(mu2)*power3(invdmsntau2mu) - 16*logmu2Q2*mH2*
     power2(mu2)*power3(invdmsntau2mu) + 16*logmstau22Q2*mH2*power2(mu2)*power3
     (invdmstau2mu) - 16*logmu2Q2*mH2*power2(mu2)*power3(invdmstau2mu))) +
     mstau22*(21*power2(mH2)*power2(sa) + power2(mh2)*(21 + 2*logmH2Q2*power2(
     invdmCH)*power2(mH2)*(-1 + power2(sa)) - (21 + 2*logmh2Q2*power2(invdmCh)*
     power2(mH2))*power2(sa)))*power3(mC2) + mh2*mH2*mstau22*(-((3 + 2*logmH2Q2
     )*mh2*power2(invdmCH)*(-1 + power2(sa))) + (3 + 2*logmh2Q2)*mH2*power2(
     invdmCh)*power2(sa))*power4(mC2) + 3*mstau22*(power2(invdmCH)*power2(mh2)*
     (-1 + power2(sa)) - power2(invdmCh)*power2(mH2)*power2(sa))*power5(mC2)))/
     (mC2*mstau22*power2(mh2)*power2(mH2)) - 156*invdmsntau2mu*power4(
     invdmstau2mu)*power5(mstau22) - 156*invdmstau1mu*power4(invdmstau2mu)*
     power5(mstau22) + 72*invdmsntau2mu*logmstau22Q2*power4(invdmstau2mu)*
     power5(mstau22) + 72*invdmstau1mu*logmstau22Q2*power4(invdmstau2mu)*power5
     (mstau22) + 168*invdmsntau2mu*logmu2Q2*power4(invdmstau2mu)*power5(mstau22
     ) + 168*invdmstau1mu*logmu2Q2*power4(invdmstau2mu)*power5(mstau22) - 48*
     invdmsntau2mu*logmstau22Q2*logmu2Q2*power4(invdmstau2mu)*power5(mstau22) -
     48*invdmstau1mu*logmstau22Q2*logmu2Q2*power4(invdmstau2mu)*power5(mstau22)
     - 24*invdmsntau2mu*power2(logmstau22Q2)*power4(invdmstau2mu)*power5(
     mstau22) - 24*invdmstau1mu*power2(logmstau22Q2)*power4(invdmstau2mu)*
     power5(mstau22) - 24*invdmsntau2mu*power2(logmu2Q2)*power4(invdmstau2mu)*
     power5(mstau22) - 24*invdmstau1mu*power2(logmu2Q2)*power4(invdmstau2mu)*
     power5(mstau22) - 8*invdmsntau2mu*power2(Pi)*power4(invdmstau2mu)*power5(
     mstau22) - 8*invdmstau1mu*power2(Pi)*power4(invdmstau2mu)*power5(mstau22)
     + 72*invdmsntau2mu*power4(invdmstau2mu)*power5(mu2) + 72*invdmstau1mu*
     power4(invdmstau2mu)*power5(mu2) - 96*invdmsntau2mu*logmu2Q2*power4(
     invdmstau2mu)*power5(mu2) - 96*invdmstau1mu*logmu2Q2*power4(invdmstau2mu)*
     power5(mu2))/384. + (mu2*DeltaInv(mstau22,msntau2,mw2)*((1 + invdmstau2mu*
     mu2)*(mstau24 + mw2*(-2*mstau22 + mw2))*(-42 + 42*invdmstau2mu*mu2 + 84*
     invdmstau2mu*mw2 - 18*invdmstau2mu*logmw2Q2*mw2 - 6*logmstau22Q2*(-3 + 3*
     invdmstau2mu*mu2 + invdmstau2mu*logmw2Q2*mw2) + 6*logmsntau2Q2*(3 - 3*
     invdmstau2mu*mu2 - 9*invdmstau2mu*mw2 + 2*invdmstau2mu*logmw2Q2*mw2 +
     logmstau22Q2*(-1 + invdmstau2mu*(mu2 + mw2))) + 3*(-1 + invdmstau2mu*(mu2
     + 3*mw2))*power2(logmsntau2Q2) + 3*(-1 + invdmstau2mu*mu2)*power2(
     logmstau22Q2) + 3*invdmstau2mu*mw2*power2(logmw2Q2) - power2(Pi) +
     invdmstau2mu*mu2*power2(Pi) + 2*invdmstau2mu*mw2*power2(Pi)) -
     invdmsntau2mu*msntau2*(3*mstau24*(42 + 6*logmsntau2Q2*(-3 + logmstau22Q2)
     - 18*logmstau22Q2 + 3*power2(logmsntau2Q2) + 3*power2(logmstau22Q2) +
     power2(Pi)) - (2*mu2 - 3*mw2)*mw2*(42 + 6*logmsntau2Q2*(-3 + logmw2Q2) -
     18*logmw2Q2 + 3*power2(logmsntau2Q2) + 3*power2(logmw2Q2) + power2(Pi)) +
     mstau22*(-2*mu2*(42 + 6*logmsntau2Q2*(-3 + logmstau22Q2) - 18*logmstau22Q2
      + 3*power2(logmsntau2Q2) + 3*power2(logmstau22Q2) + power2(Pi)) + mw2*(
     168 - 54*logmw2Q2 + 6*logmsntau2Q2*(-6 + logmstau22Q2 + logmw2Q2) + 6*
     logmstau22Q2*(-9 + 2*logmw2Q2) + 6*power2(logmsntau2Q2) + 9*power2(
     logmstau22Q2) + 9*power2(logmw2Q2) + 4*power2(Pi)))) - msntau2*(-(mw2*(42
     + 12*logmsntau2Q2*(-3 + 2*logmstau22Q2 - logmw2Q2) + 54*logmw2Q2 - 6*
     logmstau22Q2*(9 + logmw2Q2) + 6*power2(logmsntau2Q2) + 9*power2(
     logmstau22Q2) - 9*power2(logmw2Q2) + power2(Pi))) + mstau22*(42 - 18*
     logmstau22Q2 + 168*invdmstau2mu*mw2 - 36*invdmstau2mu*logmw2Q2*mw2 - 12*
     invdmstau2mu*logmstau22Q2*logmw2Q2*mw2 + 168*mu2*mw2*power2(invdmstau2mu)
     - 36*logmw2Q2*mu2*mw2*power2(invdmstau2mu) - 12*logmstau22Q2*logmw2Q2*mu2*
     mw2*power2(invdmstau2mu) + 6*logmsntau2Q2*(-3 + logmstau22Q2 + 2*
     invdmstau2mu*logmstau22Q2*mw2 + 2*invdmstau2mu*(-9 + 2*logmw2Q2)*mw2 + 2*
     logmstau22Q2*mu2*(mu2 + mw2)*power2(invdmstau2mu) - 2*mu2*(3*mu2 + (9 - 2*
     logmw2Q2)*mw2)*power2(invdmstau2mu)) + 3*(1 + 6*invdmstau2mu*mw2 + 2*mu2*(
     mu2 + 3*mw2)*power2(invdmstau2mu))*power2(logmsntau2Q2) + 3*power2(
     logmstau22Q2) + 6*invdmstau2mu*mw2*power2(logmw2Q2) + 6*mu2*mw2*power2(
     invdmstau2mu)*power2(logmw2Q2) + 84*power2(invdmstau2mu)*power2(mu2) - 36*
     logmstau22Q2*power2(invdmstau2mu)*power2(mu2) + 6*power2(invdmstau2mu)*
     power2(logmstau22Q2)*power2(mu2) + power2(Pi) + 4*invdmstau2mu*mw2*power2(
     Pi) + 4*mu2*mw2*power2(invdmstau2mu)*power2(Pi) + 2*power2(invdmstau2mu)*
     power2(mu2)*power2(Pi)) + invdmstau2mu*(power2(mu2)*(42 + 6*logmsntau2Q2*(
     -3 + logmstau22Q2) - 18*logmstau22Q2 + 3*power2(logmsntau2Q2) + 3*power2(
     logmstau22Q2) + power2(Pi)) + power2(mw2)*(-6*(9 + logmstau22Q2)*logmw2Q2
     + 6*logmsntau2Q2*(-15 + logmstau22Q2 + 4*logmw2Q2) + 15*power2(
     logmsntau2Q2) + 9*power2(logmw2Q2) + 4*(42 + power2(Pi)))) + mu2*power2(
     invdmstau2mu)*(-(power2(mu2)*(42 + 6*logmsntau2Q2*(-3 + logmstau22Q2) - 18
     *logmstau22Q2 + 3*power2(logmsntau2Q2) + 3*power2(logmstau22Q2) + power2(
     Pi))) + 3*mu2*mw2*(42 + 8*logmsntau2Q2*(-3 + logmstau22Q2) + 6*logmw2Q2 -
     2*logmstau22Q2*(9 + logmw2Q2) + 4*power2(logmsntau2Q2) + 3*power2(
     logmstau22Q2) - power2(logmw2Q2) + power2(Pi)) + power2(mw2)*(-6*(9 +
     logmstau22Q2)*logmw2Q2 + 6*logmsntau2Q2*(-15 + logmstau22Q2 + 4*logmw2Q2)
     + 15*power2(logmsntau2Q2) + 9*power2(logmw2Q2) + 4*(42 + power2(Pi))))) +
     msntau2*power2(invdmsntau2mu)*(mstau24*mu2*(42 + 6*logmsntau2Q2*(-3 +
     logmstau22Q2) - 18*logmstau22Q2 + 3*power2(logmsntau2Q2) + 3*power2(
     logmstau22Q2) + power2(Pi)) + mstau24*mw2*(42 - 36*logmstau22Q2 + 6*
     logmsntau2Q2*(-3 + 2*logmstau22Q2 - logmw2Q2) + 18*logmw2Q2 + 3*power2(
     logmsntau2Q2) + 6*power2(logmstau22Q2) - 3*power2(logmw2Q2) + power2(Pi))
     + (mu2 - mw2)*power2(mw2)*(42 + 6*logmsntau2Q2*(-3 + logmw2Q2) - 18*
     logmw2Q2 + 3*power2(logmsntau2Q2) + 3*power2(logmw2Q2) + power2(Pi)) +
     mstau22*mw2*(mw2*(42 + 18*logmstau22Q2 - 6*logmsntau2Q2*(3 + logmstau22Q2
     - 2*logmw2Q2) - 36*logmw2Q2 + 3*power2(logmsntau2Q2) - 3*power2(
     logmstau22Q2) + 6*power2(logmw2Q2) + power2(Pi)) + mu2*(168 - 18*logmw2Q2
     + 18*logmsntau2Q2*(-6 + logmstau22Q2 + logmw2Q2) - 6*logmstau22Q2*(3 + 2*
     logmw2Q2) + 18*power2(logmsntau2Q2) + 3*power2(logmstau22Q2) + 3*power2(
     logmw2Q2) + 4*power2(Pi))) - (42 + 6*logmsntau2Q2*(-3 + logmstau22Q2) - 18
     *logmstau22Q2 + 3*power2(logmsntau2Q2) + 3*power2(logmstau22Q2) + power2(
     Pi))*power6(mstau2))))/(24.*msntau2) + (mu2*DeltaInv(mstau22,mstau12,mz2)*
     ((1 + invdmstau2mu*mu2)*(mstau24 + mz2*(-2*mstau22 + mz2))*(-42 + 42*
     invdmstau2mu*mu2 + 84*invdmstau2mu*mz2 - 18*invdmstau2mu*logmz2Q2*mz2 - 6*
     logmstau22Q2*(-3 + 3*invdmstau2mu*mu2 + invdmstau2mu*logmz2Q2*mz2) + 6*
     logmstau12Q2*(3 - 3*invdmstau2mu*mu2 - 9*invdmstau2mu*mz2 + 2*invdmstau2mu
     *logmz2Q2*mz2 + logmstau22Q2*(-1 + invdmstau2mu*(mu2 + mz2))) + 3*(-1 +
     invdmstau2mu*(mu2 + 3*mz2))*power2(logmstau12Q2) + 3*(-1 + invdmstau2mu*
     mu2)*power2(logmstau22Q2) + 3*invdmstau2mu*mz2*power2(logmz2Q2) - power2(
     Pi) + invdmstau2mu*mu2*power2(Pi) + 2*invdmstau2mu*mz2*power2(Pi)) -
     invdmstau1mu*mstau12*(3*mstau24*(42 + 6*logmstau12Q2*(-3 + logmstau22Q2) -
     18*logmstau22Q2 + 3*power2(logmstau12Q2) + 3*power2(logmstau22Q2) + power2
     (Pi)) - (2*mu2 - 3*mz2)*mz2*(42 + 6*logmstau12Q2*(-3 + logmz2Q2) - 18*
     logmz2Q2 + 3*power2(logmstau12Q2) + 3*power2(logmz2Q2) + power2(Pi)) +
     mstau22*(-2*mu2*(42 + 6*logmstau12Q2*(-3 + logmstau22Q2) - 18*logmstau22Q2
      + 3*power2(logmstau12Q2) + 3*power2(logmstau22Q2) + power2(Pi)) + mz2*(
     168 - 54*logmz2Q2 + 6*logmstau12Q2*(-6 + logmstau22Q2 + logmz2Q2) + 6*
     logmstau22Q2*(-9 + 2*logmz2Q2) + 6*power2(logmstau12Q2) + 9*power2(
     logmstau22Q2) + 9*power2(logmz2Q2) + 4*power2(Pi)))) - mstau12*(-(mz2*(42
     + 12*logmstau12Q2*(-3 + 2*logmstau22Q2 - logmz2Q2) + 54*logmz2Q2 - 6*
     logmstau22Q2*(9 + logmz2Q2) + 6*power2(logmstau12Q2) + 9*power2(
     logmstau22Q2) - 9*power2(logmz2Q2) + power2(Pi))) + mstau22*(42 - 18*
     logmstau22Q2 + 168*invdmstau2mu*mz2 - 36*invdmstau2mu*logmz2Q2*mz2 - 12*
     invdmstau2mu*logmstau22Q2*logmz2Q2*mz2 + 168*mu2*mz2*power2(invdmstau2mu)
     - 36*logmz2Q2*mu2*mz2*power2(invdmstau2mu) - 12*logmstau22Q2*logmz2Q2*mu2*
     mz2*power2(invdmstau2mu) + 6*logmstau12Q2*(-3 + logmstau22Q2 + 2*
     invdmstau2mu*logmstau22Q2*mz2 + 2*invdmstau2mu*(-9 + 2*logmz2Q2)*mz2 + 2*
     logmstau22Q2*mu2*(mu2 + mz2)*power2(invdmstau2mu) - 2*mu2*(3*mu2 + (9 - 2*
     logmz2Q2)*mz2)*power2(invdmstau2mu)) + 3*(1 + 6*invdmstau2mu*mz2 + 2*mu2*(
     mu2 + 3*mz2)*power2(invdmstau2mu))*power2(logmstau12Q2) + 3*power2(
     logmstau22Q2) + 6*invdmstau2mu*mz2*power2(logmz2Q2) + 6*mu2*mz2*power2(
     invdmstau2mu)*power2(logmz2Q2) + 84*power2(invdmstau2mu)*power2(mu2) - 36*
     logmstau22Q2*power2(invdmstau2mu)*power2(mu2) + 6*power2(invdmstau2mu)*
     power2(logmstau22Q2)*power2(mu2) + power2(Pi) + 4*invdmstau2mu*mz2*power2(
     Pi) + 4*mu2*mz2*power2(invdmstau2mu)*power2(Pi) + 2*power2(invdmstau2mu)*
     power2(mu2)*power2(Pi)) + invdmstau2mu*(power2(mu2)*(42 + 6*logmstau12Q2*(
     -3 + logmstau22Q2) - 18*logmstau22Q2 + 3*power2(logmstau12Q2) + 3*power2(
     logmstau22Q2) + power2(Pi)) + power2(mz2)*(-6*(9 + logmstau22Q2)*logmz2Q2
     + 6*logmstau12Q2*(-15 + logmstau22Q2 + 4*logmz2Q2) + 15*power2(
     logmstau12Q2) + 9*power2(logmz2Q2) + 4*(42 + power2(Pi)))) + mu2*power2(
     invdmstau2mu)*(-(power2(mu2)*(42 + 6*logmstau12Q2*(-3 + logmstau22Q2) - 18
     *logmstau22Q2 + 3*power2(logmstau12Q2) + 3*power2(logmstau22Q2) + power2(
     Pi))) + 3*mu2*mz2*(42 + 8*logmstau12Q2*(-3 + logmstau22Q2) + 6*logmz2Q2 -
     2*logmstau22Q2*(9 + logmz2Q2) + 4*power2(logmstau12Q2) + 3*power2(
     logmstau22Q2) - power2(logmz2Q2) + power2(Pi)) + power2(mz2)*(-6*(9 +
     logmstau22Q2)*logmz2Q2 + 6*logmstau12Q2*(-15 + logmstau22Q2 + 4*logmz2Q2)
     + 15*power2(logmstau12Q2) + 9*power2(logmz2Q2) + 4*(42 + power2(Pi))))) +
     mstau12*power2(invdmstau1mu)*(mstau24*mu2*(42 + 6*logmstau12Q2*(-3 +
     logmstau22Q2) - 18*logmstau22Q2 + 3*power2(logmstau12Q2) + 3*power2(
     logmstau22Q2) + power2(Pi)) + mstau24*mz2*(42 - 36*logmstau22Q2 + 6*
     logmstau12Q2*(-3 + 2*logmstau22Q2 - logmz2Q2) + 18*logmz2Q2 + 3*power2(
     logmstau12Q2) + 6*power2(logmstau22Q2) - 3*power2(logmz2Q2) + power2(Pi))
     + (mu2 - mz2)*power2(mz2)*(42 + 6*logmstau12Q2*(-3 + logmz2Q2) - 18*
     logmz2Q2 + 3*power2(logmstau12Q2) + 3*power2(logmz2Q2) + power2(Pi)) +
     mstau22*mz2*(mz2*(42 + 18*logmstau22Q2 - 6*logmstau12Q2*(3 + logmstau22Q2
     - 2*logmz2Q2) - 36*logmz2Q2 + 3*power2(logmstau12Q2) - 3*power2(
     logmstau22Q2) + 6*power2(logmz2Q2) + power2(Pi)) + mu2*(168 - 18*logmz2Q2
     + 18*logmstau12Q2*(-6 + logmstau22Q2 + logmz2Q2) - 6*logmstau22Q2*(3 + 2*
     logmz2Q2) + 18*power2(logmstau12Q2) + 3*power2(logmstau22Q2) + 3*power2(
     logmz2Q2) + 4*power2(Pi))) - (42 + 6*logmstau12Q2*(-3 + logmstau22Q2) - 18
     *logmstau22Q2 + 3*power2(logmstau12Q2) + 3*power2(logmstau22Q2) + power2(
     Pi))*power6(mstau2))))/(48.*mstau12) + Fin20(mstau22,mstau12,Q2)*((Al*mu*(
     (invdmstau*(invdmstau1mu*(mstau22 - mu2) + invdmstau2mu*(-mstau12 + mu2))*
     sb)/(cb*mA2) - (3*ca*sa*(power2(mh2)*(-1 + power2(sa)) + power2(mH2)*
     power2(sa) + mh2*(mH2 - 2*mH2*power2(sa))))/(power2(mh2)*power2(mH2))))/2.
      + (power2(Al)*((2*invdmstau*(invdmstau2mu*(mstau12 - mu2) + invdmstau1mu*
     (-mstau22 + mu2)))/mA2 + 1/power2(mA2) - (3*power2(mh2 - mh2*power2(sa) +
     mH2*power2(sa)))/(power2(mh2)*power2(mH2))))/4. + (invdmstau1mu - 3*
     invdmstau2mu + 4*invdmstau1mu*invdmstau2mu*(mstau22 + mu2) + invdmstau1mu*
     mstau24*power2(invdmstau) + (-1 + invdmstau2mu*mu2)*(-(invdmstau2mu*
     mstau24) + mstau22*(-1 + invdmstau2mu*mu2))*power2(invdmstau) + mu2*power2
     (invdmstau1mu) + (-mstau12 + mu2)*power2(invdmstau2mu) + (6*mu2*power2(mh2
      - mH2)*(-1 + power2(sa))*power2(sa))/(power2(mh2)*power2(mH2)) + (4*
     invdmstau*mu2*(invdmstau1mu*(mstau22 - mu2) + invdmstau2mu*(-mstau12 + mu2
     ))*power2(sb))/mA2 + invdmstau*(invdmstau1mu*mstau22 - (1 + invdmstau2mu*(
     mstau22 - mu2))*(-1 + invdmstau2mu*mu2) + (2*mstau24 - mstau22*mu2)*power2
     (invdmstau1mu) + (4*mu2*(invdmstau2mu*(mstau12 - mu2) + invdmstau1mu*(-
     mstau22 + mu2))*power2(sb))/mz2) + 2*(mstau22 - mu2)*mu2*power3(
     invdmstau1mu) + 2*(mstau12 - mu2)*mu2*power3(invdmstau2mu) + power2(
     invdmstau)*power2(invdmstau1mu)*(-(mstau24*mu2) + power6(mstau2)))/8.) + (
     DeltaInv(mstau22,mC2,msntau2)*power2(Al)*(msntau2*mstau24*(-6*(9 +
     logmC2Q2)*logmstau22Q2 + 6*logmsntau2Q2*(-15 + logmC2Q2 + 4*logmstau22Q2)
     + 15*power2(logmsntau2Q2) + 9*power2(logmstau22Q2) + 4*(42 + power2(Pi)))
     + power2(mC2)*(msntau2*(84 - 54*logmstau22Q2 - 168*invdmsntau2mu*mstau22 -
     168*invdmstau2mu*mstau22 + 54*invdmsntau2mu*logmstau22Q2*mstau22 + 84*
     invdmsntau2mu*mu2 + 42*mstau24*power2(invdmsntau2mu) - 36*logmstau22Q2*
     mstau24*power2(invdmsntau2mu) + 168*mstau22*mu2*power2(invdmsntau2mu) - 18
     *logmstau22Q2*mstau22*mu2*power2(invdmsntau2mu) - 168*mstau22*mu2*power2(
     invdmstau2mu) + 6*logmsntau2Q2*(-9 + 18*invdmstau2mu*mstau22 + 6*
     invdmsntau2mu*(mstau22 - mu2) - 3*(mstau24 + 6*mstau22*mu2)*power2(
     invdmsntau2mu) + logmstau22Q2*(-(invdmsntau2mu*mstau22) - 2*(1 +
     invdmstau2mu*mu2)*(-2 + invdmstau2mu*(mstau22 + 2*mu2)) + (2*mstau24 + 3*
     mstau22*mu2)*power2(invdmsntau2mu)) + 6*mu2*(3*mstau22 + 2*mu2)*power2(
     invdmstau2mu)) - 3*(2 + 2*invdmstau2mu*mstau22 + invdmsntau2mu*(3*mstau22
     - 2*mu2) + (mstau24 - mstau22*mu2)*power2(invdmsntau2mu) + (2*mstau22 -
     mu2)*mu2*power2(invdmstau2mu))*power2(logmC2Q2) + 3*(3 - 6*invdmstau2mu*
     mstau22 - 2*invdmsntau2mu*(mstau22 - mu2) + (mstau24 + 6*mstau22*mu2)*
     power2(invdmsntau2mu) - 2*mu2*(3*mstau22 + 2*mu2)*power2(invdmstau2mu))*
     power2(logmsntau2Q2) + 9*power2(logmstau22Q2) - 9*invdmsntau2mu*mstau22*
     power2(logmstau22Q2) + 6*mstau24*power2(invdmsntau2mu)*power2(logmstau22Q2
     ) + 3*mstau22*mu2*power2(invdmsntau2mu)*power2(logmstau22Q2) - 126*power2(
     invdmstau2mu)*power2(mu2) + 54*logmstau22Q2*power2(invdmstau2mu)*power2(
     mu2) - 9*power2(invdmstau2mu)*power2(logmstau22Q2)*power2(mu2) + 6*
     logmC2Q2*(6 + 9*invdmsntau2mu*mstau22 + 6*invdmstau2mu*mstau22 - 6*
     invdmsntau2mu*mu2 + 3*mstau24*power2(invdmsntau2mu) - 3*mstau22*mu2*power2
     (invdmsntau2mu) + logmstau22Q2*(-2*invdmsntau2mu*mstau22 + (1 +
     invdmstau2mu*mu2)*(-1 + 2*invdmstau2mu*mstau22 + invdmstau2mu*mu2) - 2*
     mstau22*mu2*power2(invdmsntau2mu)) + 6*mstau22*mu2*power2(invdmstau2mu) -
     logmsntau2Q2*(1 + 4*invdmstau2mu*mstau22 + invdmsntau2mu*(mstau22 - 2*mu2)
     + (mstau24 - 3*mstau22*mu2)*power2(invdmsntau2mu) + 4*mstau22*mu2*power2(
     invdmstau2mu)) - 3*power2(invdmstau2mu)*power2(mu2)) + 2*power2(Pi) - 4*
     invdmsntau2mu*mstau22*power2(Pi) - 4*invdmstau2mu*mstau22*power2(Pi) + 2*
     invdmsntau2mu*mu2*power2(Pi) + mstau24*power2(invdmsntau2mu)*power2(Pi) +
     4*mstau22*mu2*power2(invdmsntau2mu)*power2(Pi) - 4*mstau22*mu2*power2(
     invdmstau2mu)*power2(Pi) - 3*power2(invdmstau2mu)*power2(mu2)*power2(Pi))
     + mstau22*(84 - 18*logmstau22Q2 + 6*logmC2Q2*(-6 + logmsntau2Q2 +
     logmstau22Q2) + 6*power2(logmC2Q2) + 3*power2(logmstau22Q2) - 84*power2(
     invdmstau2mu)*power2(mu2) + 36*logmstau22Q2*power2(invdmstau2mu)*power2(
     mu2) - 6*power2(invdmstau2mu)*power2(logmstau22Q2)*power2(mu2) + power2(
     logmsntau2Q2)*(3 - 6*power2(invdmstau2mu)*power2(mu2)) - 6*logmsntau2Q2*(3
      + 2*(-3 + logmstau22Q2)*power2(invdmstau2mu)*power2(mu2)) + 2*power2(Pi)
     - 2*power2(invdmstau2mu)*power2(mu2)*power2(Pi)) + invdmstau2mu*mstau24*(1
      + invdmstau2mu*mu2)*(6*logmC2Q2*(-3 + 2*logmsntau2Q2 - logmstau22Q2) + 6*
     logmsntau2Q2*(-9 + logmstau22Q2) + 3*power2(logmC2Q2) + 9*power2(
     logmsntau2Q2) + 2*(42 + power2(Pi)))) + (-84 + 18*logmstau22Q2 - 126*
     invdmsntau2mu*msntau2 - 168*invdmstau2mu*msntau2 - 168*invdmstau2mu*
     mstau22 + 42*msntau2*mstau22*power2(invdmsntau2mu) + 18*logmstau22Q2*
     msntau2*mstau22*power2(invdmsntau2mu) + 42*msntau2*mu2*power2(
     invdmsntau2mu) - 168*msntau2*mu2*power2(invdmstau2mu) - 168*mstau22*mu2*
     power2(invdmstau2mu) - 6*logmsntau2Q2*(-6 + logmstau22Q2 - 9*invdmsntau2mu
     *msntau2 + invdmstau2mu*logmstau22Q2*(msntau2 + 2*mstau22) - 3*
     invdmstau2mu*(5*msntau2 + 6*mstau22) + logmstau22Q2*msntau2*mstau22*power2
     (invdmsntau2mu) + 3*msntau2*(mstau22 + mu2)*power2(invdmsntau2mu) +
     logmstau22Q2*(msntau2 + 2*mstau22 - mu2)*mu2*power2(invdmstau2mu) + 3*mu2*
     (-5*msntau2 - 6*mstau22 + mu2)*power2(invdmstau2mu)) + 6*logmC2Q2*(3 + 9*
     invdmsntau2mu*msntau2 + 9*invdmstau2mu*msntau2 + invdmstau2mu*logmstau22Q2
     *msntau2 + 6*invdmstau2mu*mstau22 + 2*invdmstau2mu*logmstau22Q2*mstau22 -
     6*msntau2*mstau22*power2(invdmsntau2mu) - 3*msntau2*mu2*power2(
     invdmsntau2mu) + 9*msntau2*mu2*power2(invdmstau2mu) + logmstau22Q2*msntau2
     *mu2*power2(invdmstau2mu) + 6*mstau22*mu2*power2(invdmstau2mu) + 2*
     logmstau22Q2*mstau22*mu2*power2(invdmstau2mu) + logmsntau2Q2*(-1 - 3*
     invdmsntau2mu*msntau2 - 4*invdmstau2mu*(msntau2 + mstau22) + msntau2*(2*
     mstau22 + mu2)*power2(invdmsntau2mu) - 4*(msntau2 + mstau22)*mu2*power2(
     invdmstau2mu))) + 3*(-1 - 3*invdmsntau2mu*msntau2 - invdmstau2mu*(3*
     msntau2 + 2*mstau22) + msntau2*(2*mstau22 + mu2)*power2(invdmsntau2mu) - (
     3*msntau2 + 2*mstau22)*mu2*power2(invdmstau2mu))*power2(logmC2Q2) + 3*(-2
     - 3*invdmsntau2mu*msntau2 - invdmstau2mu*(5*msntau2 + 6*mstau22) + msntau2
     *(mstau22 + mu2)*power2(invdmsntau2mu) + mu2*(-5*msntau2 - 6*mstau22 + mu2
     )*power2(invdmstau2mu))*power2(logmsntau2Q2) - 3*power2(logmstau22Q2) - 3*
     msntau2*mstau22*power2(invdmsntau2mu)*power2(logmstau22Q2) + 42*power2(
     invdmstau2mu)*power2(mu2) - 18*logmstau22Q2*power2(invdmstau2mu)*power2(
     mu2) + 3*power2(invdmstau2mu)*power2(logmstau22Q2)*power2(mu2) - 2*power2(
     Pi) - 3*invdmsntau2mu*msntau2*power2(Pi) - 4*invdmstau2mu*msntau2*power2(
     Pi) - 4*invdmstau2mu*mstau22*power2(Pi) + msntau2*mstau22*power2(
     invdmsntau2mu)*power2(Pi) + msntau2*mu2*power2(invdmsntau2mu)*power2(Pi) -
     4*msntau2*mu2*power2(invdmstau2mu)*power2(Pi) - 4*mstau22*mu2*power2(
     invdmstau2mu)*power2(Pi) + power2(invdmstau2mu)*power2(mu2)*power2(Pi))*
     power3(mC2) + (-(msntau2*power2(invdmsntau2mu)*(42 + 6*logmC2Q2*(-3 +
     logmsntau2Q2) - 18*logmsntau2Q2 + 3*power2(logmC2Q2) + 3*power2(
     logmsntau2Q2) + power2(Pi))) + invdmstau2mu*(6*logmC2Q2*(-3 + 2*
     logmsntau2Q2 - logmstau22Q2) + 6*logmsntau2Q2*(-9 + logmstau22Q2) + 3*
     power2(logmC2Q2) + 9*power2(logmsntau2Q2) + 2*(42 + power2(Pi))) + mu2*
     power2(invdmstau2mu)*(6*logmC2Q2*(-3 + 2*logmsntau2Q2 - logmstau22Q2) + 6*
     logmsntau2Q2*(-9 + logmstau22Q2) + 3*power2(logmC2Q2) + 9*power2(
     logmsntau2Q2) + 2*(42 + power2(Pi))))*power4(mC2) - (-6*(3 + logmC2Q2)*
     logmstau22Q2 + 6*logmsntau2Q2*(-9 + logmC2Q2 + 2*logmstau22Q2) + 9*power2(
     logmsntau2Q2) + 3*power2(logmstau22Q2) + 2*(42 + power2(Pi)))*power6(
     mstau2) + mC2*(-(invdmsntau2mu*msntau2*(3*mstau24 - 2*mstau22*mu2)*(42 + 6
     *logmsntau2Q2*(-3 + logmstau22Q2) - 18*logmstau22Q2 + 3*power2(
     logmsntau2Q2) + 3*power2(logmstau22Q2) + power2(Pi))) + mstau24*(84 + 6*
     logmC2Q2*(3 + logmsntau2Q2 - 2*logmstau22Q2) - 18*logmstau22Q2 - 3*power2(
     logmC2Q2) + 3*power2(logmstau22Q2) + 42*power2(invdmstau2mu)*power2(mu2) -
     18*logmstau22Q2*power2(invdmstau2mu)*power2(mu2) + 3*power2(invdmstau2mu)*
     power2(logmstau22Q2)*power2(mu2) + 3*power2(logmsntau2Q2)*(4 + power2(
     invdmstau2mu)*power2(mu2)) + 6*logmsntau2Q2*(-12 + 3*logmstau22Q2 - 3*
     power2(invdmstau2mu)*power2(mu2) + logmstau22Q2*power2(invdmstau2mu)*
     power2(mu2)) + 2*power2(Pi) + power2(invdmstau2mu)*power2(mu2)*power2(Pi))
     + msntau2*(invdmstau2mu*(-1 + invdmstau2mu*mu2)*power2(mu2)*(42 + 6*
     logmsntau2Q2*(-3 + logmstau22Q2) - 18*logmstau22Q2 + 3*power2(logmsntau2Q2
     ) + 3*power2(logmstau22Q2) + power2(Pi)) + mstau22*(252 + 18*logmC2Q2*(-3
     + 2*logmsntau2Q2 - logmstau22Q2) + 9*power2(logmC2Q2) - 84*power2(
     invdmstau2mu)*power2(mu2) + 36*logmstau22Q2*power2(invdmstau2mu)*power2(
     mu2) - 6*power2(invdmstau2mu)*power2(logmstau22Q2)*power2(mu2) + power2(
     logmsntau2Q2)*(27 - 6*power2(invdmstau2mu)*power2(mu2)) - 6*logmsntau2Q2*(
     27 - 3*logmstau22Q2 - 6*power2(invdmstau2mu)*power2(mu2) + 2*logmstau22Q2*
     power2(invdmstau2mu)*power2(mu2)) + 6*power2(Pi) - 2*power2(invdmstau2mu)*
     power2(mu2)*power2(Pi))) - msntau2*power2(invdmsntau2mu)*(42 + 6*
     logmsntau2Q2*(-3 + logmstau22Q2) - 18*logmstau22Q2 + 3*power2(logmsntau2Q2
     ) + 3*power2(logmstau22Q2) + power2(Pi))*(-(mstau24*mu2) + power6(mstau2))
     )))/(24.*mC2*msntau2) + (DeltaInv(mA2,mstau22,mstau12)*power2(Al)*(-2*
     mstau12*mstau24*(-6*(9 + logmA2Q2)*logmstau22Q2 + 6*logmstau12Q2*(-15 +
     logmA2Q2 + 4*logmstau22Q2) + 15*power2(logmstau12Q2) + 9*power2(
     logmstau22Q2) + 4*(42 + power2(Pi))) + power2(mA2)*(mstau12*(-42 - 168*
     invdmstau1mu*mstau22 - 168*invdmstau2mu*mstau22 + 36*invdmstau1mu*
     logmstau12Q2*mstau22 + 108*invdmstau2mu*logmstau12Q2*mstau22 + 84*
     invdmstau1mu*mu2 - 36*invdmstau1mu*logmstau12Q2*mu2 + 42*mstau24*power2(
     invdmstau1mu) - 18*logmstau12Q2*mstau24*power2(invdmstau1mu) + 168*mstau22
     *mu2*power2(invdmstau1mu) - 108*logmstau12Q2*mstau22*mu2*power2(
     invdmstau1mu) - 168*mstau22*mu2*power2(invdmstau2mu) + 108*logmstau12Q2*
     mstau22*mu2*power2(invdmstau2mu) - 3*(5 + 2*invdmstau2mu*mstau22 +
     invdmstau1mu*(3*mstau22 - 2*mu2) + (mstau24 - mstau22*mu2)*power2(
     invdmstau1mu) + (2*mstau22 - mu2)*mu2*power2(invdmstau2mu))*power2(
     logmA2Q2) - 6*invdmstau1mu*mstau22*power2(logmstau12Q2) - 18*invdmstau2mu*
     mstau22*power2(logmstau12Q2) + 6*invdmstau1mu*mu2*power2(logmstau12Q2) + 3
     *mstau24*power2(invdmstau1mu)*power2(logmstau12Q2) + 18*mstau22*mu2*power2
     (invdmstau1mu)*power2(logmstau12Q2) - 18*mstau22*mu2*power2(invdmstau2mu)*
     power2(logmstau12Q2) - 126*power2(invdmstau2mu)*power2(mu2) + 72*
     logmstau12Q2*power2(invdmstau2mu)*power2(mu2) - 12*power2(invdmstau2mu)*
     power2(logmstau12Q2)*power2(mu2) + 3*power2(logmstau22Q2)*(3 - 3*
     invdmstau1mu*mstau22 + (2*mstau24 + mstau22*mu2)*power2(invdmstau1mu) - 3*
     power2(invdmstau2mu)*power2(mu2)) + 6*logmA2Q2*(15 + 9*invdmstau1mu*
     mstau22 + 6*invdmstau2mu*mstau22 - 6*invdmstau1mu*mu2 + 3*mstau24*power2(
     invdmstau1mu) - 3*mstau22*mu2*power2(invdmstau1mu) + logmstau22Q2*(-2*
     invdmstau1mu*mstau22 + (1 + invdmstau2mu*mu2)*(-1 + 2*invdmstau2mu*mstau22
      + invdmstau2mu*mu2) - 2*mstau22*mu2*power2(invdmstau1mu)) + 6*mstau22*mu2
     *power2(invdmstau2mu) + logmstau12Q2*(-(invdmstau1mu*(mstau22 - 2*mu2)) -
     (mstau24 - 3*mstau22*mu2)*power2(invdmstau1mu) - 4*(1 + invdmstau2mu*
     mstau22 + mstau22*mu2*power2(invdmstau2mu))) - 3*power2(invdmstau2mu)*
     power2(mu2)) + 6*logmstau22Q2*(-9 + 9*invdmstau1mu*mstau22 - 3*(2*mstau24
     + mstau22*mu2)*power2(invdmstau1mu) + logmstau12Q2*(-(invdmstau1mu*mstau22
     ) - 2*(1 + invdmstau2mu*mu2)*(-2 + invdmstau2mu*mstau22 + 2*invdmstau2mu*
     mu2) + (2*mstau24 + 3*mstau22*mu2)*power2(invdmstau1mu)) + 9*power2(
     invdmstau2mu)*power2(mu2)) - power2(Pi) - 4*invdmstau1mu*mstau22*power2(Pi
     ) - 4*invdmstau2mu*mstau22*power2(Pi) + 2*invdmstau1mu*mu2*power2(Pi) +
     mstau24*power2(invdmstau1mu)*power2(Pi) + 4*mstau22*mu2*power2(
     invdmstau1mu)*power2(Pi) - 4*mstau22*mu2*power2(invdmstau2mu)*power2(Pi) -
     3*power2(invdmstau2mu)*power2(mu2)*power2(Pi)) - 2*mstau22*(-42 + 36*
     logmstau22Q2 + 6*logmA2Q2*(-6 + logmstau12Q2 + logmstau22Q2) + 6*power2(
     logmA2Q2) - 6*power2(logmstau22Q2) + 42*power2(invdmstau2mu)*power2(mu2) -
     18*logmstau22Q2*power2(invdmstau2mu)*power2(mu2) + 3*power2(invdmstau2mu)*
     power2(logmstau22Q2)*power2(mu2) + 3*power2(logmstau12Q2)*(-2 + power2(
     invdmstau2mu)*power2(mu2)) + 6*logmstau12Q2*(6 - 3*logmstau22Q2 - 3*power2
     (invdmstau2mu)*power2(mu2) + logmstau22Q2*power2(invdmstau2mu)*power2(mu2)
     ) - power2(Pi) + power2(invdmstau2mu)*power2(mu2)*power2(Pi)) +
     invdmstau2mu*mstau24*(1 + invdmstau2mu*mu2)*(6*logmA2Q2*(-3 + 2*
     logmstau12Q2 - logmstau22Q2) + 6*logmstau12Q2*(-9 + logmstau22Q2) + 3*
     power2(logmA2Q2) + 9*power2(logmstau12Q2) + 2*(42 + power2(Pi)))) + (42 +
     18*logmstau22Q2 - 126*invdmstau1mu*mstau12 - 168*invdmstau2mu*mstau12 -
     168*invdmstau2mu*mstau22 + 42*mstau12*mstau22*power2(invdmstau1mu) + 18*
     logmstau22Q2*mstau12*mstau22*power2(invdmstau1mu) + 42*mstau12*mu2*power2(
     invdmstau1mu) - 168*mstau12*mu2*power2(invdmstau2mu) - 168*mstau22*mu2*
     power2(invdmstau2mu) - 6*logmstau12Q2*(3 + logmstau22Q2 - 9*invdmstau1mu*
     mstau12 + invdmstau2mu*logmstau22Q2*(mstau12 + 2*mstau22) - 3*invdmstau2mu
     *(5*mstau12 + 6*mstau22) + logmstau22Q2*mstau12*mstau22*power2(
     invdmstau1mu) + 3*mstau12*(mstau22 + mu2)*power2(invdmstau1mu) +
     logmstau22Q2*(mstau12 + 2*mstau22 - mu2)*mu2*power2(invdmstau2mu) + 3*mu2*
     (-5*mstau12 - 6*mstau22 + mu2)*power2(invdmstau2mu)) + 6*logmA2Q2*(-6 + 9*
     invdmstau1mu*mstau12 + 9*invdmstau2mu*mstau12 + invdmstau2mu*logmstau22Q2*
     mstau12 + 6*invdmstau2mu*mstau22 + 2*invdmstau2mu*logmstau22Q2*mstau22 - 3
     *mstau12*(2*mstau22 + mu2)*power2(invdmstau1mu) + 9*mstau12*mu2*power2(
     invdmstau2mu) + logmstau22Q2*mstau12*mu2*power2(invdmstau2mu) + 6*mstau22*
     mu2*power2(invdmstau2mu) + 2*logmstau22Q2*mstau22*mu2*power2(invdmstau2mu)
     + logmstau12Q2*(2 - 3*invdmstau1mu*mstau12 - 4*invdmstau2mu*(mstau12 +
     mstau22) + mstau12*(2*mstau22 + mu2)*power2(invdmstau1mu) - 4*(mstau12 +
     mstau22)*mu2*power2(invdmstau2mu))) + 3*(2 - 3*invdmstau1mu*mstau12 -
     invdmstau2mu*(3*mstau12 + 2*mstau22) + mstau12*(2*mstau22 + mu2)*power2(
     invdmstau1mu) - (3*mstau12 + 2*mstau22)*mu2*power2(invdmstau2mu))*power2(
     logmA2Q2) + 3*(1 - 3*invdmstau1mu*mstau12 - invdmstau2mu*(5*mstau12 + 6*
     mstau22) + mstau12*(mstau22 + mu2)*power2(invdmstau1mu) + mu2*(-5*mstau12
     - 6*mstau22 + mu2)*power2(invdmstau2mu))*power2(logmstau12Q2) - 3*power2(
     logmstau22Q2) - 3*mstau12*mstau22*power2(invdmstau1mu)*power2(logmstau22Q2
     ) + 42*power2(invdmstau2mu)*power2(mu2) - 18*logmstau22Q2*power2(
     invdmstau2mu)*power2(mu2) + 3*power2(invdmstau2mu)*power2(logmstau22Q2)*
     power2(mu2) + power2(Pi) - 3*invdmstau1mu*mstau12*power2(Pi) - 4*
     invdmstau2mu*mstau12*power2(Pi) - 4*invdmstau2mu*mstau22*power2(Pi) +
     mstau12*mstau22*power2(invdmstau1mu)*power2(Pi) + mstau12*mu2*power2(
     invdmstau1mu)*power2(Pi) - 4*mstau12*mu2*power2(invdmstau2mu)*power2(Pi) -
     4*mstau22*mu2*power2(invdmstau2mu)*power2(Pi) + power2(invdmstau2mu)*
     power2(mu2)*power2(Pi))*power3(mA2) + (-(mstau12*power2(invdmstau1mu)*(42
     + 6*logmA2Q2*(-3 + logmstau12Q2) - 18*logmstau12Q2 + 3*power2(logmA2Q2) +
     3*power2(logmstau12Q2) + power2(Pi))) + invdmstau2mu*(6*logmA2Q2*(-3 + 2*
     logmstau12Q2 - logmstau22Q2) + 6*logmstau12Q2*(-9 + logmstau22Q2) + 3*
     power2(logmA2Q2) + 9*power2(logmstau12Q2) + 2*(42 + power2(Pi))) + mu2*
     power2(invdmstau2mu)*(6*logmA2Q2*(-3 + 2*logmstau12Q2 - logmstau22Q2) + 6*
     logmstau12Q2*(-9 + logmstau22Q2) + 3*power2(logmA2Q2) + 9*power2(
     logmstau12Q2) + 2*(42 + power2(Pi))))*power4(mA2) + 2*(-6*(3 + logmA2Q2)*
     logmstau22Q2 + 6*logmstau12Q2*(-9 + logmA2Q2 + 2*logmstau22Q2) + 9*power2(
     logmstau12Q2) + 3*power2(logmstau22Q2) + 2*(42 + power2(Pi)))*power6(
     mstau2) - mA2*(invdmstau1mu*mstau12*(3*mstau24 - 2*mstau22*mu2)*(42 + 6*
     logmstau12Q2*(-3 + logmstau22Q2) - 18*logmstau22Q2 + 3*power2(logmstau12Q2
     ) + 3*power2(logmstau22Q2) + power2(Pi)) + mstau24*(294 + 12*logmA2Q2*(3 +
     logmstau12Q2 - 2*logmstau22Q2) - 90*logmstau22Q2 - 6*power2(logmA2Q2) + 15
     *power2(logmstau22Q2) - 42*power2(invdmstau2mu)*power2(mu2) + 18*
     logmstau22Q2*power2(invdmstau2mu)*power2(mu2) - 3*power2(invdmstau2mu)*
     power2(logmstau22Q2)*power2(mu2) + power2(logmstau12Q2)*(33 - 3*power2(
     invdmstau2mu)*power2(mu2)) - 6*logmstau12Q2*(33 - 9*logmstau22Q2 - 3*
     power2(invdmstau2mu)*power2(mu2) + logmstau22Q2*power2(invdmstau2mu)*
     power2(mu2)) + 7*power2(Pi) - power2(invdmstau2mu)*power2(mu2)*power2(Pi))
     + mstau12*(-(invdmstau2mu*(-1 + invdmstau2mu*mu2)*power2(mu2)*(42 + 6*
     logmstau12Q2*(-3 + logmstau22Q2) - 18*logmstau22Q2 + 3*power2(logmstau12Q2
     ) + 3*power2(logmstau22Q2) + power2(Pi))) + mstau22*(630 + 36*logmA2Q2*(-3
      + 2*logmstau12Q2 - logmstau22Q2) - 54*logmstau22Q2 + 18*power2(logmA2Q2)
     + 9*power2(logmstau22Q2) + 84*power2(invdmstau2mu)*power2(mu2) - 36*
     logmstau22Q2*power2(invdmstau2mu)*power2(mu2) + 6*power2(invdmstau2mu)*
     power2(logmstau22Q2)*power2(mu2) + power2(logmstau12Q2)*(63 + 6*power2(
     invdmstau2mu)*power2(mu2)) + 6*logmstau12Q2*(-63 + 9*logmstau22Q2 - 6*
     power2(invdmstau2mu)*power2(mu2) + 2*logmstau22Q2*power2(invdmstau2mu)*
     power2(mu2)) + 15*power2(Pi) + 2*power2(invdmstau2mu)*power2(mu2)*power2(
     Pi))) + mstau12*power2(invdmstau1mu)*(42 + 6*logmstau12Q2*(-3 +
     logmstau22Q2) - 18*logmstau22Q2 + 3*power2(logmstau12Q2) + 3*power2(
     logmstau22Q2) + power2(Pi))*(-(mstau24*mu2) + power6(mstau2)))))/(48.*mA2*
     mstau12) + Fin20(mstau22,msntau2,Q2)*(-power2(Al)/(4.*power2(mC2)) + (
     invdmsntau2mu*(mstau24 + (1 + invdmsntau2mu*mu2 + 2*(mstau22 - mu2)*mu2*
     power2(invdmsntau2mu))*power2(msntau2) + 2*(mstau22 - mu2)*mu2*power2(
     invdmsntau2mu)*power2(mstau22) - msntau2*(mstau22 - 2*invdmsntau2mu*
     mstau24 + 3*invdmsntau2mu*mstau22*mu2 + 4*mu2*power2(invdmsntau2mu)*power2
     (mstau22) - 4*mstau22*power2(invdmsntau2mu)*power2(mu2)) + invdmsntau2mu*(
     -2*mstau22*mstau24 - mstau24*mu2 + 2*mu2*power2(mstau22) + power6(mstau2))
     ))/(8.*power2(msntau2 - mstau22))) + DeltaInv(mstau22,mstau12,mh2)*(-(Al*
     ca*mu*sa*(-(power2(mh2)*(invdmstau2mu*mstau24*(1 + invdmstau2mu*mu2)*(6*
     logmh2Q2*(-3 + 2*logmstau12Q2 - logmstau22Q2) + 6*logmstau12Q2*(-9 +
     logmstau22Q2) + 3*power2(logmh2Q2) + 9*power2(logmstau12Q2) + 2*(42 +
     power2(Pi))) + mstau12*(42 - 54*logmstau22Q2 - 168*invdmstau1mu*mstau22 -
     168*invdmstau2mu*mstau22 + 54*invdmstau1mu*logmstau22Q2*mstau22 + 84*
     invdmstau1mu*mu2 + 42*mstau24*power2(invdmstau1mu) - 36*logmstau22Q2*
     mstau24*power2(invdmstau1mu) + 168*mstau22*mu2*power2(invdmstau1mu) - 18*
     logmstau22Q2*mstau22*mu2*power2(invdmstau1mu) - 168*mstau22*mu2*power2(
     invdmstau2mu) + 9*power2(logmstau22Q2) - 9*invdmstau1mu*mstau22*power2(
     logmstau22Q2) + 6*mstau24*power2(invdmstau1mu)*power2(logmstau22Q2) + 3*
     mstau22*mu2*power2(invdmstau1mu)*power2(logmstau22Q2) - 126*power2(
     invdmstau2mu)*power2(mu2) + 54*logmstau22Q2*power2(invdmstau2mu)*power2(
     mu2) - 9*power2(invdmstau2mu)*power2(logmstau22Q2)*power2(mu2) + power2(Pi
     ) - 4*invdmstau1mu*mstau22*power2(Pi) - 4*invdmstau2mu*mstau22*power2(Pi)
     + 2*invdmstau1mu*mu2*power2(Pi) + mstau24*power2(invdmstau1mu)*power2(Pi)
     + 4*mstau22*mu2*power2(invdmstau1mu)*power2(Pi) - 4*mstau22*mu2*power2(
     invdmstau2mu)*power2(Pi) - 3*power2(invdmstau2mu)*power2(mu2)*power2(Pi) +
     6*logmstau12Q2*(-6 + 18*invdmstau2mu*mstau22 - invdmstau1mu*logmstau22Q2*
     mstau22 + 6*invdmstau1mu*(mstau22 - mu2) - 2*logmstau22Q2*(1 +
     invdmstau2mu*mu2)*(-2 + invdmstau2mu*(mstau22 + 2*mu2)) + logmstau22Q2*(2*
     mstau24 + 3*mstau22*mu2)*power2(invdmstau1mu) - 3*(mstau24 + 6*mstau22*mu2
     )*power2(invdmstau1mu) + 6*mu2*(3*mstau22 + 2*mu2)*power2(invdmstau2mu) -
     18*power2(sa)) + 6*logmh2Q2*(9 + 9*invdmstau1mu*mstau22 + 6*invdmstau2mu*
     mstau22 - 6*invdmstau1mu*mu2 + 3*mstau24*power2(invdmstau1mu) - 3*mstau22*
     mu2*power2(invdmstau1mu) + logmstau22Q2*(-2*invdmstau1mu*mstau22 + (1 +
     invdmstau2mu*mu2)*(-1 + 2*invdmstau2mu*mstau22 + invdmstau2mu*mu2) - 2*
     mstau22*mu2*power2(invdmstau1mu)) + 6*mstau22*mu2*power2(invdmstau2mu) - 3
     *power2(invdmstau2mu)*power2(mu2) - logmstau12Q2*(2 + 4*invdmstau2mu*
     mstau22 + invdmstau1mu*(mstau22 - 2*mu2) + (mstau24 - 3*mstau22*mu2)*
     power2(invdmstau1mu) + 4*mstau22*mu2*power2(invdmstau2mu) - 6*power2(sa))
     - 18*power2(sa)) + 252*power2(sa) + 6*power2(Pi)*power2(sa) + 3*power2(
     logmh2Q2)*(-3 - 2*invdmstau2mu*mstau22 + invdmstau1mu*(-3*mstau22 + 2*mu2)
     + (-mstau24 + mstau22*mu2)*power2(invdmstau1mu) + mu2*(-2*mstau22 + mu2)*
     power2(invdmstau2mu) + 6*power2(sa)) + 3*power2(logmstau12Q2)*(2 - 6*
     invdmstau2mu*mstau22 - 2*invdmstau1mu*(mstau22 - mu2) + (mstau24 + 6*
     mstau22*mu2)*power2(invdmstau1mu) - 2*mu2*(3*mstau22 + 2*mu2)*power2(
     invdmstau2mu) + 6*power2(sa))) - 2*mstau22*(-42 + 42*power2(invdmstau2mu)*
     power2(mu2) - power2(Pi) + power2(invdmstau2mu)*power2(mu2)*power2(Pi) +
     108*logmh2Q2*power2(sa) - 18*power2(logmh2Q2)*power2(sa) + 3*power2(
     logmstau12Q2)*(-1 + power2(invdmstau2mu)*power2(mu2) + 3*power2(sa)) + 3*
     power2(logmstau22Q2)*(-1 + power2(invdmstau2mu)*power2(mu2) + 3*power2(sa)
     ) - 18*logmstau22Q2*(-1 + power2(invdmstau2mu)*power2(mu2) + (3 + logmh2Q2
     )*power2(sa)) + 6*logmstau12Q2*(3 - 3*power2(invdmstau2mu)*power2(mu2) - 3
     *(3 + logmh2Q2)*power2(sa) + logmstau22Q2*(-1 + power2(invdmstau2mu)*
     power2(mu2) + 6*power2(sa)))))) + (42 + 126*invdmstau1mu*mstau12 + 168*
     invdmstau2mu*mstau12 - 54*invdmstau1mu*logmh2Q2*mstau12 - 54*invdmstau2mu*
     logmh2Q2*mstau12 + 168*invdmstau2mu*mstau22 - 36*invdmstau2mu*logmh2Q2*
     mstau22 - 42*mstau12*mstau22*power2(invdmstau1mu) + 36*logmh2Q2*mstau12*
     mstau22*power2(invdmstau1mu) - 42*mstau12*mu2*power2(invdmstau1mu) + 18*
     logmh2Q2*mstau12*mu2*power2(invdmstau1mu) + 168*mstau12*mu2*power2(
     invdmstau2mu) - 54*logmh2Q2*mstau12*mu2*power2(invdmstau2mu) + 168*mstau22
     *mu2*power2(invdmstau2mu) - 36*logmh2Q2*mstau22*mu2*power2(invdmstau2mu) -
     6*logmstau22Q2*(3 + invdmstau2mu*logmh2Q2*(mstau12 + 2*mstau22) + 3*
     mstau12*mstau22*power2(invdmstau1mu) + (logmh2Q2*(mstau12 + 2*mstau22) - 3
     *mu2)*mu2*power2(invdmstau2mu)) + 9*invdmstau1mu*mstau12*power2(logmh2Q2)
     + 9*invdmstau2mu*mstau12*power2(logmh2Q2) + 6*invdmstau2mu*mstau22*power2(
     logmh2Q2) - 6*mstau12*mstau22*power2(invdmstau1mu)*power2(logmh2Q2) - 3*
     mstau12*mu2*power2(invdmstau1mu)*power2(logmh2Q2) + 9*mstau12*mu2*power2(
     invdmstau2mu)*power2(logmh2Q2) + 6*mstau22*mu2*power2(invdmstau2mu)*power2
     (logmh2Q2) - 42*power2(invdmstau2mu)*power2(mu2) + 3*power2(logmstau22Q2)*
     (1 + mstau12*mstau22*power2(invdmstau1mu) - power2(invdmstau2mu)*power2(
     mu2)) + power2(Pi) + 3*invdmstau1mu*mstau12*power2(Pi) + 4*invdmstau2mu*
     mstau12*power2(Pi) + 4*invdmstau2mu*mstau22*power2(Pi) - mstau12*mstau22*
     power2(invdmstau1mu)*power2(Pi) - mstau12*mu2*power2(invdmstau1mu)*power2(
     Pi) + 4*mstau12*mu2*power2(invdmstau2mu)*power2(Pi) + 4*mstau22*mu2*power2
     (invdmstau2mu)*power2(Pi) - power2(invdmstau2mu)*power2(mu2)*power2(Pi) -
     3*power2(logmstau12Q2)*(-1 - 3*invdmstau1mu*mstau12 - invdmstau2mu*(5*
     mstau12 + 6*mstau22) + mstau12*(mstau22 + mu2)*power2(invdmstau1mu) + mu2*
     (-5*mstau12 - 6*mstau22 + mu2)*power2(invdmstau2mu) - 6*power2(sa)) + 252*
     power2(sa) - 108*logmh2Q2*power2(sa) + 18*power2(logmh2Q2)*power2(sa) + 6*
     power2(Pi)*power2(sa) + 6*logmstau12Q2*(-3 - 15*invdmstau2mu*mstau12 + 3*
     invdmstau1mu*(-3 + logmh2Q2)*mstau12 + 4*invdmstau2mu*logmh2Q2*mstau12 -
     18*invdmstau2mu*mstau22 + 4*invdmstau2mu*logmh2Q2*mstau22 + mstau12*(3*
     mstau22 - 2*logmh2Q2*mstau22 + 3*mu2 - logmh2Q2*mu2)*power2(invdmstau1mu)
     - 15*mstau12*mu2*power2(invdmstau2mu) + 4*logmh2Q2*mstau12*mu2*power2(
     invdmstau2mu) - 18*mstau22*mu2*power2(invdmstau2mu) + 4*logmh2Q2*mstau22*
     mu2*power2(invdmstau2mu) + logmstau22Q2*(1 + invdmstau2mu*(mstau12 + 2*
     mstau22) + mstau12*mstau22*power2(invdmstau1mu) + (mstau12 + 2*mstau22 -
     mu2)*mu2*power2(invdmstau2mu)) + 3*power2(invdmstau2mu)*power2(mu2) - 18*
     power2(sa) + 6*logmh2Q2*power2(sa)))*power3(mh2) - (-(mstau12*power2(
     invdmstau1mu)*(42 + 6*logmh2Q2*(-3 + logmstau12Q2) - 18*logmstau12Q2 + 3*
     power2(logmh2Q2) + 3*power2(logmstau12Q2) + power2(Pi))) + invdmstau2mu*(6
     *logmh2Q2*(-3 + 2*logmstau12Q2 - logmstau22Q2) + 6*logmstau12Q2*(-9 +
     logmstau22Q2) + 3*power2(logmh2Q2) + 9*power2(logmstau12Q2) + 2*(42 +
     power2(Pi))) + mu2*power2(invdmstau2mu)*(6*logmh2Q2*(-3 + 2*logmstau12Q2 -
     logmstau22Q2) + 6*logmstau12Q2*(-9 + logmstau22Q2) + 3*power2(logmh2Q2) +
     9*power2(logmstau12Q2) + 2*(42 + power2(Pi))))*power4(mh2) + 6*power2(sa)*
     (-(mstau12*mstau24*(-6*(9 + logmh2Q2)*logmstau22Q2 + 6*logmstau12Q2*(-15 +
     logmh2Q2 + 4*logmstau22Q2) + 15*power2(logmstau12Q2) + 9*power2(
     logmstau22Q2) + 4*(42 + power2(Pi)))) + (-6*(3 + logmh2Q2)*logmstau22Q2 +
     6*logmstau12Q2*(-9 + logmh2Q2 + 2*logmstau22Q2) + 9*power2(logmstau12Q2) +
     3*power2(logmstau22Q2) + 2*(42 + power2(Pi)))*power6(mstau2)) + mh2*(
     invdmstau1mu*mstau12*(3*mstau24 - 2*mstau22*mu2)*(42 + 6*logmstau12Q2*(-3
     + logmstau22Q2) - 18*logmstau22Q2 + 3*power2(logmstau12Q2) + 3*power2(
     logmstau22Q2) + power2(Pi)) - mstau24*(-42 + 42*power2(invdmstau2mu)*
     power2(mu2) - power2(Pi) + power2(invdmstau2mu)*power2(mu2)*power2(Pi) +
     756*power2(sa) + 108*logmh2Q2*power2(sa) - 18*power2(logmh2Q2)*power2(sa)
     + 18*power2(Pi)*power2(sa) + 3*power2(logmstau22Q2)*(-1 + power2(
     invdmstau2mu)*power2(mu2) + 12*power2(sa)) + 3*power2(logmstau12Q2)*(-1 +
     power2(invdmstau2mu)*power2(mu2) + 30*power2(sa)) - 18*logmstau22Q2*(-1 +
     power2(invdmstau2mu)*power2(mu2) + 4*(3 + logmh2Q2)*power2(sa)) + 6*
     logmstau12Q2*(3 - 3*power2(invdmstau2mu)*power2(mu2) - 90*power2(sa) + 6*
     logmh2Q2*power2(sa) + logmstau22Q2*(-1 + power2(invdmstau2mu)*power2(mu2)
     + 24*power2(sa)))) + mstau12*(-(invdmstau2mu*(-1 + invdmstau2mu*mu2)*
     power2(mu2)*(42 + 6*logmstau12Q2*(-3 + logmstau22Q2) - 18*logmstau22Q2 + 3
     *power2(logmstau12Q2) + 3*power2(logmstau22Q2) + power2(Pi))) + mstau22*(
     42 + 84*power2(invdmstau2mu)*power2(mu2) + power2(Pi) + 2*power2(
     invdmstau2mu)*power2(mu2)*power2(Pi) + 3*power2(logmstau12Q2)*(1 + 2*
     power2(invdmstau2mu)*power2(mu2) - 60*power2(sa)) + 3*power2(logmstau22Q2)
     *(1 + 2*power2(invdmstau2mu)*power2(mu2) - 6*power2(sa)) - 1764*power2(sa)
     + 324*logmh2Q2*power2(sa) - 54*power2(logmh2Q2)*power2(sa) - 42*power2(Pi)
     *power2(sa) - 18*logmstau22Q2*(1 + 2*power2(invdmstau2mu)*power2(mu2) - 6*
     (1 + logmh2Q2)*power2(sa)) + 6*logmstau12Q2*(-3 + logmstau22Q2 - 6*power2(
     invdmstau2mu)*power2(mu2) + 2*logmstau22Q2*power2(invdmstau2mu)*power2(mu2
     ) + 180*power2(sa) - 36*logmh2Q2*power2(sa) - 24*logmstau22Q2*power2(sa)))
     ) + mstau12*power2(invdmstau1mu)*(42 + 6*logmstau12Q2*(-3 + logmstau22Q2)
     - 18*logmstau22Q2 + 3*power2(logmstau12Q2) + 3*power2(logmstau22Q2) +
     power2(Pi))*(-(mstau24*mu2) + power6(mstau2)))))/(24.*mh2*mstau12) + (mu2*
     (-1 + power2(sa))*(-(power2(mh2)*(invdmstau2mu*mstau24*(1 + invdmstau2mu*
     mu2)*(6*logmh2Q2*(-3 + 2*logmstau12Q2 - logmstau22Q2) + 6*logmstau12Q2*(-9
      + logmstau22Q2) + 3*power2(logmh2Q2) + 9*power2(logmstau12Q2) + 2*(42 +
     power2(Pi))) + mstau12*(42 - 54*logmstau22Q2 - 168*invdmstau1mu*mstau22 -
     168*invdmstau2mu*mstau22 + 54*invdmstau1mu*logmstau22Q2*mstau22 + 84*
     invdmstau1mu*mu2 + 42*mstau24*power2(invdmstau1mu) - 36*logmstau22Q2*
     mstau24*power2(invdmstau1mu) + 168*mstau22*mu2*power2(invdmstau1mu) - 18*
     logmstau22Q2*mstau22*mu2*power2(invdmstau1mu) - 168*mstau22*mu2*power2(
     invdmstau2mu) + 9*power2(logmstau22Q2) - 9*invdmstau1mu*mstau22*power2(
     logmstau22Q2) + 6*mstau24*power2(invdmstau1mu)*power2(logmstau22Q2) + 3*
     mstau22*mu2*power2(invdmstau1mu)*power2(logmstau22Q2) - 126*power2(
     invdmstau2mu)*power2(mu2) + 54*logmstau22Q2*power2(invdmstau2mu)*power2(
     mu2) - 9*power2(invdmstau2mu)*power2(logmstau22Q2)*power2(mu2) + power2(Pi
     ) - 4*invdmstau1mu*mstau22*power2(Pi) - 4*invdmstau2mu*mstau22*power2(Pi)
     + 2*invdmstau1mu*mu2*power2(Pi) + mstau24*power2(invdmstau1mu)*power2(Pi)
     + 4*mstau22*mu2*power2(invdmstau1mu)*power2(Pi) - 4*mstau22*mu2*power2(
     invdmstau2mu)*power2(Pi) - 3*power2(invdmstau2mu)*power2(mu2)*power2(Pi) +
     6*logmstau12Q2*(-6 + 18*invdmstau2mu*mstau22 - invdmstau1mu*logmstau22Q2*
     mstau22 + 6*invdmstau1mu*(mstau22 - mu2) - 2*logmstau22Q2*(1 +
     invdmstau2mu*mu2)*(-2 + invdmstau2mu*(mstau22 + 2*mu2)) + logmstau22Q2*(2*
     mstau24 + 3*mstau22*mu2)*power2(invdmstau1mu) - 3*(mstau24 + 6*mstau22*mu2
     )*power2(invdmstau1mu) + 6*mu2*(3*mstau22 + 2*mu2)*power2(invdmstau2mu) -
     18*power2(sa)) + 6*logmh2Q2*(9 + 9*invdmstau1mu*mstau22 + 6*invdmstau2mu*
     mstau22 - 6*invdmstau1mu*mu2 + 3*mstau24*power2(invdmstau1mu) - 3*mstau22*
     mu2*power2(invdmstau1mu) + logmstau22Q2*(-2*invdmstau1mu*mstau22 + (1 +
     invdmstau2mu*mu2)*(-1 + 2*invdmstau2mu*mstau22 + invdmstau2mu*mu2) - 2*
     mstau22*mu2*power2(invdmstau1mu)) + 6*mstau22*mu2*power2(invdmstau2mu) - 3
     *power2(invdmstau2mu)*power2(mu2) - logmstau12Q2*(2 + 4*invdmstau2mu*
     mstau22 + invdmstau1mu*(mstau22 - 2*mu2) + (mstau24 - 3*mstau22*mu2)*
     power2(invdmstau1mu) + 4*mstau22*mu2*power2(invdmstau2mu) - 6*power2(sa))
     - 18*power2(sa)) + 252*power2(sa) + 6*power2(Pi)*power2(sa) + 3*power2(
     logmh2Q2)*(-3 - 2*invdmstau2mu*mstau22 + invdmstau1mu*(-3*mstau22 + 2*mu2)
     + (-mstau24 + mstau22*mu2)*power2(invdmstau1mu) + mu2*(-2*mstau22 + mu2)*
     power2(invdmstau2mu) + 6*power2(sa)) + 3*power2(logmstau12Q2)*(2 - 6*
     invdmstau2mu*mstau22 - 2*invdmstau1mu*(mstau22 - mu2) + (mstau24 + 6*
     mstau22*mu2)*power2(invdmstau1mu) - 2*mu2*(3*mstau22 + 2*mu2)*power2(
     invdmstau2mu) + 6*power2(sa))) - 2*mstau22*(-42 + 42*power2(invdmstau2mu)*
     power2(mu2) - power2(Pi) + power2(invdmstau2mu)*power2(mu2)*power2(Pi) +
     108*logmh2Q2*power2(sa) - 18*power2(logmh2Q2)*power2(sa) + 3*power2(
     logmstau12Q2)*(-1 + power2(invdmstau2mu)*power2(mu2) + 3*power2(sa)) + 3*
     power2(logmstau22Q2)*(-1 + power2(invdmstau2mu)*power2(mu2) + 3*power2(sa)
     ) - 18*logmstau22Q2*(-1 + power2(invdmstau2mu)*power2(mu2) + (3 + logmh2Q2
     )*power2(sa)) + 6*logmstau12Q2*(3 - 3*power2(invdmstau2mu)*power2(mu2) - 3
     *(3 + logmh2Q2)*power2(sa) + logmstau22Q2*(-1 + power2(invdmstau2mu)*
     power2(mu2) + 6*power2(sa)))))) + (42 + 126*invdmstau1mu*mstau12 + 168*
     invdmstau2mu*mstau12 - 54*invdmstau1mu*logmh2Q2*mstau12 - 54*invdmstau2mu*
     logmh2Q2*mstau12 + 168*invdmstau2mu*mstau22 - 36*invdmstau2mu*logmh2Q2*
     mstau22 - 42*mstau12*mstau22*power2(invdmstau1mu) + 36*logmh2Q2*mstau12*
     mstau22*power2(invdmstau1mu) - 42*mstau12*mu2*power2(invdmstau1mu) + 18*
     logmh2Q2*mstau12*mu2*power2(invdmstau1mu) + 168*mstau12*mu2*power2(
     invdmstau2mu) - 54*logmh2Q2*mstau12*mu2*power2(invdmstau2mu) + 168*mstau22
     *mu2*power2(invdmstau2mu) - 36*logmh2Q2*mstau22*mu2*power2(invdmstau2mu) -
     6*logmstau22Q2*(3 + invdmstau2mu*logmh2Q2*(mstau12 + 2*mstau22) + 3*
     mstau12*mstau22*power2(invdmstau1mu) + (logmh2Q2*(mstau12 + 2*mstau22) - 3
     *mu2)*mu2*power2(invdmstau2mu)) + 9*invdmstau1mu*mstau12*power2(logmh2Q2)
     + 9*invdmstau2mu*mstau12*power2(logmh2Q2) + 6*invdmstau2mu*mstau22*power2(
     logmh2Q2) - 6*mstau12*mstau22*power2(invdmstau1mu)*power2(logmh2Q2) - 3*
     mstau12*mu2*power2(invdmstau1mu)*power2(logmh2Q2) + 9*mstau12*mu2*power2(
     invdmstau2mu)*power2(logmh2Q2) + 6*mstau22*mu2*power2(invdmstau2mu)*power2
     (logmh2Q2) - 42*power2(invdmstau2mu)*power2(mu2) + 3*power2(logmstau22Q2)*
     (1 + mstau12*mstau22*power2(invdmstau1mu) - power2(invdmstau2mu)*power2(
     mu2)) + power2(Pi) + 3*invdmstau1mu*mstau12*power2(Pi) + 4*invdmstau2mu*
     mstau12*power2(Pi) + 4*invdmstau2mu*mstau22*power2(Pi) - mstau12*mstau22*
     power2(invdmstau1mu)*power2(Pi) - mstau12*mu2*power2(invdmstau1mu)*power2(
     Pi) + 4*mstau12*mu2*power2(invdmstau2mu)*power2(Pi) + 4*mstau22*mu2*power2
     (invdmstau2mu)*power2(Pi) - power2(invdmstau2mu)*power2(mu2)*power2(Pi) -
     3*power2(logmstau12Q2)*(-1 - 3*invdmstau1mu*mstau12 - invdmstau2mu*(5*
     mstau12 + 6*mstau22) + mstau12*(mstau22 + mu2)*power2(invdmstau1mu) + mu2*
     (-5*mstau12 - 6*mstau22 + mu2)*power2(invdmstau2mu) - 6*power2(sa)) + 252*
     power2(sa) - 108*logmh2Q2*power2(sa) + 18*power2(logmh2Q2)*power2(sa) + 6*
     power2(Pi)*power2(sa) + 6*logmstau12Q2*(-3 - 15*invdmstau2mu*mstau12 + 3*
     invdmstau1mu*(-3 + logmh2Q2)*mstau12 + 4*invdmstau2mu*logmh2Q2*mstau12 -
     18*invdmstau2mu*mstau22 + 4*invdmstau2mu*logmh2Q2*mstau22 + mstau12*(3*
     mstau22 - 2*logmh2Q2*mstau22 + 3*mu2 - logmh2Q2*mu2)*power2(invdmstau1mu)
     - 15*mstau12*mu2*power2(invdmstau2mu) + 4*logmh2Q2*mstau12*mu2*power2(
     invdmstau2mu) - 18*mstau22*mu2*power2(invdmstau2mu) + 4*logmh2Q2*mstau22*
     mu2*power2(invdmstau2mu) + logmstau22Q2*(1 + invdmstau2mu*(mstau12 + 2*
     mstau22) + mstau12*mstau22*power2(invdmstau1mu) + (mstau12 + 2*mstau22 -
     mu2)*mu2*power2(invdmstau2mu)) + 3*power2(invdmstau2mu)*power2(mu2) - 18*
     power2(sa) + 6*logmh2Q2*power2(sa)))*power3(mh2) - (-(mstau12*power2(
     invdmstau1mu)*(42 + 6*logmh2Q2*(-3 + logmstau12Q2) - 18*logmstau12Q2 + 3*
     power2(logmh2Q2) + 3*power2(logmstau12Q2) + power2(Pi))) + invdmstau2mu*(6
     *logmh2Q2*(-3 + 2*logmstau12Q2 - logmstau22Q2) + 6*logmstau12Q2*(-9 +
     logmstau22Q2) + 3*power2(logmh2Q2) + 9*power2(logmstau12Q2) + 2*(42 +
     power2(Pi))) + mu2*power2(invdmstau2mu)*(6*logmh2Q2*(-3 + 2*logmstau12Q2 -
     logmstau22Q2) + 6*logmstau12Q2*(-9 + logmstau22Q2) + 3*power2(logmh2Q2) +
     9*power2(logmstau12Q2) + 2*(42 + power2(Pi))))*power4(mh2) + 6*power2(sa)*
     (-(mstau12*mstau24*(-6*(9 + logmh2Q2)*logmstau22Q2 + 6*logmstau12Q2*(-15 +
     logmh2Q2 + 4*logmstau22Q2) + 15*power2(logmstau12Q2) + 9*power2(
     logmstau22Q2) + 4*(42 + power2(Pi)))) + (-6*(3 + logmh2Q2)*logmstau22Q2 +
     6*logmstau12Q2*(-9 + logmh2Q2 + 2*logmstau22Q2) + 9*power2(logmstau12Q2) +
     3*power2(logmstau22Q2) + 2*(42 + power2(Pi)))*power6(mstau2)) + mh2*(
     invdmstau1mu*mstau12*(3*mstau24 - 2*mstau22*mu2)*(42 + 6*logmstau12Q2*(-3
     + logmstau22Q2) - 18*logmstau22Q2 + 3*power2(logmstau12Q2) + 3*power2(
     logmstau22Q2) + power2(Pi)) - mstau24*(-42 + 42*power2(invdmstau2mu)*
     power2(mu2) - power2(Pi) + power2(invdmstau2mu)*power2(mu2)*power2(Pi) +
     756*power2(sa) + 108*logmh2Q2*power2(sa) - 18*power2(logmh2Q2)*power2(sa)
     + 18*power2(Pi)*power2(sa) + 3*power2(logmstau22Q2)*(-1 + power2(
     invdmstau2mu)*power2(mu2) + 12*power2(sa)) + 3*power2(logmstau12Q2)*(-1 +
     power2(invdmstau2mu)*power2(mu2) + 30*power2(sa)) - 18*logmstau22Q2*(-1 +
     power2(invdmstau2mu)*power2(mu2) + 4*(3 + logmh2Q2)*power2(sa)) + 6*
     logmstau12Q2*(3 - 3*power2(invdmstau2mu)*power2(mu2) - 90*power2(sa) + 6*
     logmh2Q2*power2(sa) + logmstau22Q2*(-1 + power2(invdmstau2mu)*power2(mu2)
     + 24*power2(sa)))) + mstau12*(-(invdmstau2mu*(-1 + invdmstau2mu*mu2)*
     power2(mu2)*(42 + 6*logmstau12Q2*(-3 + logmstau22Q2) - 18*logmstau22Q2 + 3
     *power2(logmstau12Q2) + 3*power2(logmstau22Q2) + power2(Pi))) + mstau22*(
     42 + 84*power2(invdmstau2mu)*power2(mu2) + power2(Pi) + 2*power2(
     invdmstau2mu)*power2(mu2)*power2(Pi) + 3*power2(logmstau12Q2)*(1 + 2*
     power2(invdmstau2mu)*power2(mu2) - 60*power2(sa)) + 3*power2(logmstau22Q2)
     *(1 + 2*power2(invdmstau2mu)*power2(mu2) - 6*power2(sa)) - 1764*power2(sa)
     + 324*logmh2Q2*power2(sa) - 54*power2(logmh2Q2)*power2(sa) - 42*power2(Pi)
     *power2(sa) - 18*logmstau22Q2*(1 + 2*power2(invdmstau2mu)*power2(mu2) - 6*
     (1 + logmh2Q2)*power2(sa)) + 6*logmstau12Q2*(-3 + logmstau22Q2 - 6*power2(
     invdmstau2mu)*power2(mu2) + 2*logmstau22Q2*power2(invdmstau2mu)*power2(mu2
     ) + 180*power2(sa) - 36*logmh2Q2*power2(sa) - 24*logmstau22Q2*power2(sa)))
     ) + mstau12*power2(invdmstau1mu)*(42 + 6*logmstau12Q2*(-3 + logmstau22Q2)
     - 18*logmstau22Q2 + 3*power2(logmstau12Q2) + 3*power2(logmstau22Q2) +
     power2(Pi))*(-(mstau24*mu2) + power6(mstau2)))))/(48.*mh2*mstau12) - (
     power2(Al)*power2(sa)*(-(power2(mh2)*(invdmstau2mu*mstau24*(1 +
     invdmstau2mu*mu2)*(6*logmh2Q2*(-3 + 2*logmstau12Q2 - logmstau22Q2) + 6*
     logmstau12Q2*(-9 + logmstau22Q2) + 3*power2(logmh2Q2) + 9*power2(
     logmstau12Q2) + 2*(42 + power2(Pi))) + mstau12*(42 - 54*logmstau22Q2 - 168
     *invdmstau1mu*mstau22 - 168*invdmstau2mu*mstau22 + 54*invdmstau1mu*
     logmstau22Q2*mstau22 + 84*invdmstau1mu*mu2 + 42*mstau24*power2(
     invdmstau1mu) - 36*logmstau22Q2*mstau24*power2(invdmstau1mu) + 168*mstau22
     *mu2*power2(invdmstau1mu) - 18*logmstau22Q2*mstau22*mu2*power2(
     invdmstau1mu) - 168*mstau22*mu2*power2(invdmstau2mu) + 9*power2(
     logmstau22Q2) - 9*invdmstau1mu*mstau22*power2(logmstau22Q2) + 6*mstau24*
     power2(invdmstau1mu)*power2(logmstau22Q2) + 3*mstau22*mu2*power2(
     invdmstau1mu)*power2(logmstau22Q2) - 126*power2(invdmstau2mu)*power2(mu2)
     + 54*logmstau22Q2*power2(invdmstau2mu)*power2(mu2) - 9*power2(invdmstau2mu
     )*power2(logmstau22Q2)*power2(mu2) + power2(Pi) - 4*invdmstau1mu*mstau22*
     power2(Pi) - 4*invdmstau2mu*mstau22*power2(Pi) + 2*invdmstau1mu*mu2*power2
     (Pi) + mstau24*power2(invdmstau1mu)*power2(Pi) + 4*mstau22*mu2*power2(
     invdmstau1mu)*power2(Pi) - 4*mstau22*mu2*power2(invdmstau2mu)*power2(Pi) -
     3*power2(invdmstau2mu)*power2(mu2)*power2(Pi) + 6*logmstau12Q2*(-6 + 18*
     invdmstau2mu*mstau22 - invdmstau1mu*logmstau22Q2*mstau22 + 6*invdmstau1mu*
     (mstau22 - mu2) - 2*logmstau22Q2*(1 + invdmstau2mu*mu2)*(-2 + invdmstau2mu
     *(mstau22 + 2*mu2)) + logmstau22Q2*(2*mstau24 + 3*mstau22*mu2)*power2(
     invdmstau1mu) - 3*(mstau24 + 6*mstau22*mu2)*power2(invdmstau1mu) + 6*mu2*(
     3*mstau22 + 2*mu2)*power2(invdmstau2mu) - 18*power2(sa)) + 6*logmh2Q2*(9 +
     9*invdmstau1mu*mstau22 + 6*invdmstau2mu*mstau22 - 6*invdmstau1mu*mu2 + 3*
     mstau24*power2(invdmstau1mu) - 3*mstau22*mu2*power2(invdmstau1mu) +
     logmstau22Q2*(-2*invdmstau1mu*mstau22 + (1 + invdmstau2mu*mu2)*(-1 + 2*
     invdmstau2mu*mstau22 + invdmstau2mu*mu2) - 2*mstau22*mu2*power2(
     invdmstau1mu)) + 6*mstau22*mu2*power2(invdmstau2mu) - 3*power2(
     invdmstau2mu)*power2(mu2) - logmstau12Q2*(2 + 4*invdmstau2mu*mstau22 +
     invdmstau1mu*(mstau22 - 2*mu2) + (mstau24 - 3*mstau22*mu2)*power2(
     invdmstau1mu) + 4*mstau22*mu2*power2(invdmstau2mu) - 6*power2(sa)) - 18*
     power2(sa)) + 252*power2(sa) + 6*power2(Pi)*power2(sa) + 3*power2(logmh2Q2
     )*(-3 - 2*invdmstau2mu*mstau22 + invdmstau1mu*(-3*mstau22 + 2*mu2) + (-
     mstau24 + mstau22*mu2)*power2(invdmstau1mu) + mu2*(-2*mstau22 + mu2)*
     power2(invdmstau2mu) + 6*power2(sa)) + 3*power2(logmstau12Q2)*(2 - 6*
     invdmstau2mu*mstau22 - 2*invdmstau1mu*(mstau22 - mu2) + (mstau24 + 6*
     mstau22*mu2)*power2(invdmstau1mu) - 2*mu2*(3*mstau22 + 2*mu2)*power2(
     invdmstau2mu) + 6*power2(sa))) - 2*mstau22*(-42 + 42*power2(invdmstau2mu)*
     power2(mu2) - power2(Pi) + power2(invdmstau2mu)*power2(mu2)*power2(Pi) +
     108*logmh2Q2*power2(sa) - 18*power2(logmh2Q2)*power2(sa) + 3*power2(
     logmstau12Q2)*(-1 + power2(invdmstau2mu)*power2(mu2) + 3*power2(sa)) + 3*
     power2(logmstau22Q2)*(-1 + power2(invdmstau2mu)*power2(mu2) + 3*power2(sa)
     ) - 18*logmstau22Q2*(-1 + power2(invdmstau2mu)*power2(mu2) + (3 + logmh2Q2
     )*power2(sa)) + 6*logmstau12Q2*(3 - 3*power2(invdmstau2mu)*power2(mu2) - 3
     *(3 + logmh2Q2)*power2(sa) + logmstau22Q2*(-1 + power2(invdmstau2mu)*
     power2(mu2) + 6*power2(sa)))))) + (42 + 126*invdmstau1mu*mstau12 + 168*
     invdmstau2mu*mstau12 - 54*invdmstau1mu*logmh2Q2*mstau12 - 54*invdmstau2mu*
     logmh2Q2*mstau12 + 168*invdmstau2mu*mstau22 - 36*invdmstau2mu*logmh2Q2*
     mstau22 - 42*mstau12*mstau22*power2(invdmstau1mu) + 36*logmh2Q2*mstau12*
     mstau22*power2(invdmstau1mu) - 42*mstau12*mu2*power2(invdmstau1mu) + 18*
     logmh2Q2*mstau12*mu2*power2(invdmstau1mu) + 168*mstau12*mu2*power2(
     invdmstau2mu) - 54*logmh2Q2*mstau12*mu2*power2(invdmstau2mu) + 168*mstau22
     *mu2*power2(invdmstau2mu) - 36*logmh2Q2*mstau22*mu2*power2(invdmstau2mu) -
     6*logmstau22Q2*(3 + invdmstau2mu*logmh2Q2*(mstau12 + 2*mstau22) + 3*
     mstau12*mstau22*power2(invdmstau1mu) + (logmh2Q2*(mstau12 + 2*mstau22) - 3
     *mu2)*mu2*power2(invdmstau2mu)) + 9*invdmstau1mu*mstau12*power2(logmh2Q2)
     + 9*invdmstau2mu*mstau12*power2(logmh2Q2) + 6*invdmstau2mu*mstau22*power2(
     logmh2Q2) - 6*mstau12*mstau22*power2(invdmstau1mu)*power2(logmh2Q2) - 3*
     mstau12*mu2*power2(invdmstau1mu)*power2(logmh2Q2) + 9*mstau12*mu2*power2(
     invdmstau2mu)*power2(logmh2Q2) + 6*mstau22*mu2*power2(invdmstau2mu)*power2
     (logmh2Q2) - 42*power2(invdmstau2mu)*power2(mu2) + 3*power2(logmstau22Q2)*
     (1 + mstau12*mstau22*power2(invdmstau1mu) - power2(invdmstau2mu)*power2(
     mu2)) + power2(Pi) + 3*invdmstau1mu*mstau12*power2(Pi) + 4*invdmstau2mu*
     mstau12*power2(Pi) + 4*invdmstau2mu*mstau22*power2(Pi) - mstau12*mstau22*
     power2(invdmstau1mu)*power2(Pi) - mstau12*mu2*power2(invdmstau1mu)*power2(
     Pi) + 4*mstau12*mu2*power2(invdmstau2mu)*power2(Pi) + 4*mstau22*mu2*power2
     (invdmstau2mu)*power2(Pi) - power2(invdmstau2mu)*power2(mu2)*power2(Pi) -
     3*power2(logmstau12Q2)*(-1 - 3*invdmstau1mu*mstau12 - invdmstau2mu*(5*
     mstau12 + 6*mstau22) + mstau12*(mstau22 + mu2)*power2(invdmstau1mu) + mu2*
     (-5*mstau12 - 6*mstau22 + mu2)*power2(invdmstau2mu) - 6*power2(sa)) + 252*
     power2(sa) - 108*logmh2Q2*power2(sa) + 18*power2(logmh2Q2)*power2(sa) + 6*
     power2(Pi)*power2(sa) + 6*logmstau12Q2*(-3 - 15*invdmstau2mu*mstau12 + 3*
     invdmstau1mu*(-3 + logmh2Q2)*mstau12 + 4*invdmstau2mu*logmh2Q2*mstau12 -
     18*invdmstau2mu*mstau22 + 4*invdmstau2mu*logmh2Q2*mstau22 + mstau12*(3*
     mstau22 - 2*logmh2Q2*mstau22 + 3*mu2 - logmh2Q2*mu2)*power2(invdmstau1mu)
     - 15*mstau12*mu2*power2(invdmstau2mu) + 4*logmh2Q2*mstau12*mu2*power2(
     invdmstau2mu) - 18*mstau22*mu2*power2(invdmstau2mu) + 4*logmh2Q2*mstau22*
     mu2*power2(invdmstau2mu) + logmstau22Q2*(1 + invdmstau2mu*(mstau12 + 2*
     mstau22) + mstau12*mstau22*power2(invdmstau1mu) + (mstau12 + 2*mstau22 -
     mu2)*mu2*power2(invdmstau2mu)) + 3*power2(invdmstau2mu)*power2(mu2) - 18*
     power2(sa) + 6*logmh2Q2*power2(sa)))*power3(mh2) - (-(mstau12*power2(
     invdmstau1mu)*(42 + 6*logmh2Q2*(-3 + logmstau12Q2) - 18*logmstau12Q2 + 3*
     power2(logmh2Q2) + 3*power2(logmstau12Q2) + power2(Pi))) + invdmstau2mu*(6
     *logmh2Q2*(-3 + 2*logmstau12Q2 - logmstau22Q2) + 6*logmstau12Q2*(-9 +
     logmstau22Q2) + 3*power2(logmh2Q2) + 9*power2(logmstau12Q2) + 2*(42 +
     power2(Pi))) + mu2*power2(invdmstau2mu)*(6*logmh2Q2*(-3 + 2*logmstau12Q2 -
     logmstau22Q2) + 6*logmstau12Q2*(-9 + logmstau22Q2) + 3*power2(logmh2Q2) +
     9*power2(logmstau12Q2) + 2*(42 + power2(Pi))))*power4(mh2) + 6*power2(sa)*
     (-(mstau12*mstau24*(-6*(9 + logmh2Q2)*logmstau22Q2 + 6*logmstau12Q2*(-15 +
     logmh2Q2 + 4*logmstau22Q2) + 15*power2(logmstau12Q2) + 9*power2(
     logmstau22Q2) + 4*(42 + power2(Pi)))) + (-6*(3 + logmh2Q2)*logmstau22Q2 +
     6*logmstau12Q2*(-9 + logmh2Q2 + 2*logmstau22Q2) + 9*power2(logmstau12Q2) +
     3*power2(logmstau22Q2) + 2*(42 + power2(Pi)))*power6(mstau2)) + mh2*(
     invdmstau1mu*mstau12*(3*mstau24 - 2*mstau22*mu2)*(42 + 6*logmstau12Q2*(-3
     + logmstau22Q2) - 18*logmstau22Q2 + 3*power2(logmstau12Q2) + 3*power2(
     logmstau22Q2) + power2(Pi)) - mstau24*(-42 + 42*power2(invdmstau2mu)*
     power2(mu2) - power2(Pi) + power2(invdmstau2mu)*power2(mu2)*power2(Pi) +
     756*power2(sa) + 108*logmh2Q2*power2(sa) - 18*power2(logmh2Q2)*power2(sa)
     + 18*power2(Pi)*power2(sa) + 3*power2(logmstau22Q2)*(-1 + power2(
     invdmstau2mu)*power2(mu2) + 12*power2(sa)) + 3*power2(logmstau12Q2)*(-1 +
     power2(invdmstau2mu)*power2(mu2) + 30*power2(sa)) - 18*logmstau22Q2*(-1 +
     power2(invdmstau2mu)*power2(mu2) + 4*(3 + logmh2Q2)*power2(sa)) + 6*
     logmstau12Q2*(3 - 3*power2(invdmstau2mu)*power2(mu2) - 90*power2(sa) + 6*
     logmh2Q2*power2(sa) + logmstau22Q2*(-1 + power2(invdmstau2mu)*power2(mu2)
     + 24*power2(sa)))) + mstau12*(-(invdmstau2mu*(-1 + invdmstau2mu*mu2)*
     power2(mu2)*(42 + 6*logmstau12Q2*(-3 + logmstau22Q2) - 18*logmstau22Q2 + 3
     *power2(logmstau12Q2) + 3*power2(logmstau22Q2) + power2(Pi))) + mstau22*(
     42 + 84*power2(invdmstau2mu)*power2(mu2) + power2(Pi) + 2*power2(
     invdmstau2mu)*power2(mu2)*power2(Pi) + 3*power2(logmstau12Q2)*(1 + 2*
     power2(invdmstau2mu)*power2(mu2) - 60*power2(sa)) + 3*power2(logmstau22Q2)
     *(1 + 2*power2(invdmstau2mu)*power2(mu2) - 6*power2(sa)) - 1764*power2(sa)
     + 324*logmh2Q2*power2(sa) - 54*power2(logmh2Q2)*power2(sa) - 42*power2(Pi)
     *power2(sa) - 18*logmstau22Q2*(1 + 2*power2(invdmstau2mu)*power2(mu2) - 6*
     (1 + logmh2Q2)*power2(sa)) + 6*logmstau12Q2*(-3 + logmstau22Q2 - 6*power2(
     invdmstau2mu)*power2(mu2) + 2*logmstau22Q2*power2(invdmstau2mu)*power2(mu2
     ) + 180*power2(sa) - 36*logmh2Q2*power2(sa) - 24*logmstau22Q2*power2(sa)))
     ) + mstau12*power2(invdmstau1mu)*(42 + 6*logmstau12Q2*(-3 + logmstau22Q2)
     - 18*logmstau22Q2 + 3*power2(logmstau12Q2) + 3*power2(logmstau22Q2) +
     power2(Pi))*(-(mstau24*mu2) + power6(mstau2)))))/(48.*mh2*mstau12)) +
     DeltaInv(mstau22,mH2,mstau12)*(-(power2(Al)*(-1 + power2(sa))*(power2(mH2)
     *(invdmstau2mu*mstau24*(1 + invdmstau2mu*mu2)*(6*logmH2Q2*(-3 + 2*
     logmstau12Q2 - logmstau22Q2) + 6*logmstau12Q2*(-9 + logmstau22Q2) + 3*
     power2(logmH2Q2) + 9*power2(logmstau12Q2) + 2*(42 + power2(Pi))) - 2*
     mstau22*(-42 - 36*logmstau22Q2 + 6*power2(logmstau22Q2) + 42*power2(
     invdmstau2mu)*power2(mu2) - 18*logmstau22Q2*power2(invdmstau2mu)*power2(
     mu2) + 3*power2(invdmstau2mu)*power2(logmstau22Q2)*power2(mu2) - power2(Pi
     ) + power2(invdmstau2mu)*power2(mu2)*power2(Pi) + 3*power2(logmstau12Q2)*(
     2 + power2(invdmstau2mu)*power2(mu2) - 3*power2(sa)) + 18*logmH2Q2*(-6 +
     logmstau12Q2 + logmstau22Q2)*(-1 + power2(sa)) + 18*power2(logmH2Q2)*(-1 +
     power2(sa)) + 54*logmstau22Q2*power2(sa) - 9*power2(logmstau22Q2)*power2(
     sa) + 6*logmstau12Q2*(-6 - 3*power2(invdmstau2mu)*power2(mu2) +
     logmstau22Q2*(5 + power2(invdmstau2mu)*power2(mu2) - 6*power2(sa)) + 9*
     power2(sa))) + mstau12*(294 - 54*logmstau22Q2 - 168*invdmstau1mu*mstau22 -
     168*invdmstau2mu*mstau22 + 54*invdmstau1mu*logmstau22Q2*mstau22 + 84*
     invdmstau1mu*mu2 + 42*mstau24*power2(invdmstau1mu) - 36*logmstau22Q2*
     mstau24*power2(invdmstau1mu) + 168*mstau22*mu2*power2(invdmstau1mu) - 18*
     logmstau22Q2*mstau22*mu2*power2(invdmstau1mu) - 168*mstau22*mu2*power2(
     invdmstau2mu) + 9*power2(logmstau22Q2) - 9*invdmstau1mu*mstau22*power2(
     logmstau22Q2) + 6*mstau24*power2(invdmstau1mu)*power2(logmstau22Q2) + 3*
     mstau22*mu2*power2(invdmstau1mu)*power2(logmstau22Q2) - 126*power2(
     invdmstau2mu)*power2(mu2) + 54*logmstau22Q2*power2(invdmstau2mu)*power2(
     mu2) - 9*power2(invdmstau2mu)*power2(logmstau22Q2)*power2(mu2) + 7*power2(
     Pi) - 4*invdmstau1mu*mstau22*power2(Pi) - 4*invdmstau2mu*mstau22*power2(Pi
     ) + 2*invdmstau1mu*mu2*power2(Pi) + mstau24*power2(invdmstau1mu)*power2(Pi
     ) + 4*mstau22*mu2*power2(invdmstau1mu)*power2(Pi) - 4*mstau22*mu2*power2(
     invdmstau2mu)*power2(Pi) - 3*power2(invdmstau2mu)*power2(mu2)*power2(Pi) +
     3*power2(logmH2Q2)*(3 - 2*invdmstau2mu*mstau22 + invdmstau1mu*(-3*mstau22
     + 2*mu2) + (-mstau24 + mstau22*mu2)*power2(invdmstau1mu) + mu2*(-2*mstau22
      + mu2)*power2(invdmstau2mu) - 6*power2(sa)) + 3*power2(logmstau12Q2)*(8 -
     6*invdmstau2mu*mstau22 - 2*invdmstau1mu*(mstau22 - mu2) + (mstau24 + 6*
     mstau22*mu2)*power2(invdmstau1mu) - 2*mu2*(3*mstau22 + 2*mu2)*power2(
     invdmstau2mu) - 6*power2(sa)) - 252*power2(sa) - 6*power2(Pi)*power2(sa) +
     6*logmstau12Q2*(6*invdmstau1mu*(mstau22 - mu2) - 3*(mstau24 + 6*mstau22*
     mu2)*power2(invdmstau1mu) + logmstau22Q2*(-(invdmstau1mu*mstau22) - 2*(1 +
     invdmstau2mu*mu2)*(-2 + invdmstau2mu*mstau22 + 2*invdmstau2mu*mu2) + (2*
     mstau24 + 3*mstau22*mu2)*power2(invdmstau1mu)) + 6*(-4 + 3*invdmstau2mu*
     mstau22 + mu2*(3*mstau22 + 2*mu2)*power2(invdmstau2mu) + 3*power2(sa))) +
     6*logmH2Q2*(-9 + 9*invdmstau1mu*mstau22 + 6*invdmstau2mu*mstau22 - 6*
     invdmstau1mu*mu2 + 3*mstau24*power2(invdmstau1mu) - 3*mstau22*mu2*power2(
     invdmstau1mu) + logmstau22Q2*(-2*invdmstau1mu*mstau22 + (1 + invdmstau2mu*
     mu2)*(-1 + 2*invdmstau2mu*mstau22 + invdmstau2mu*mu2) - 2*mstau22*mu2*
     power2(invdmstau1mu)) + 6*mstau22*mu2*power2(invdmstau2mu) - 3*power2(
     invdmstau2mu)*power2(mu2) + 18*power2(sa) - logmstau12Q2*(-4 + 4*
     invdmstau2mu*mstau22 + invdmstau1mu*(mstau22 - 2*mu2) + (mstau24 - 3*
     mstau22*mu2)*power2(invdmstau1mu) + 4*mstau22*mu2*power2(invdmstau2mu) + 6
     *power2(sa))))) + (-294 + 18*logmstau22Q2 - 126*invdmstau1mu*mstau12 - 168
     *invdmstau2mu*mstau12 - 168*invdmstau2mu*mstau22 + 42*mstau12*mstau22*
     power2(invdmstau1mu) + 18*logmstau22Q2*mstau12*mstau22*power2(invdmstau1mu
     ) + 42*mstau12*mu2*power2(invdmstau1mu) - 168*mstau12*mu2*power2(
     invdmstau2mu) - 168*mstau22*mu2*power2(invdmstau2mu) - 3*power2(
     logmstau22Q2) - 3*mstau12*mstau22*power2(invdmstau1mu)*power2(logmstau22Q2
     ) + 42*power2(invdmstau2mu)*power2(mu2) - 18*logmstau22Q2*power2(
     invdmstau2mu)*power2(mu2) + 3*power2(invdmstau2mu)*power2(logmstau22Q2)*
     power2(mu2) - 7*power2(Pi) - 3*invdmstau1mu*mstau12*power2(Pi) - 4*
     invdmstau2mu*mstau12*power2(Pi) - 4*invdmstau2mu*mstau22*power2(Pi) +
     mstau12*mstau22*power2(invdmstau1mu)*power2(Pi) + mstau12*mu2*power2(
     invdmstau1mu)*power2(Pi) - 4*mstau12*mu2*power2(invdmstau2mu)*power2(Pi) -
     4*mstau22*mu2*power2(invdmstau2mu)*power2(Pi) + power2(invdmstau2mu)*
     power2(mu2)*power2(Pi) + 3*power2(logmH2Q2)*(-3*invdmstau1mu*mstau12 -
     invdmstau2mu*(3*mstau12 + 2*mstau22) + mstau12*(2*mstau22 + mu2)*power2(
     invdmstau1mu) - (3*mstau12 + 2*mstau22)*mu2*power2(invdmstau2mu) + 6*(-1 +
     power2(sa))) + 252*power2(sa) + 6*power2(Pi)*power2(sa) + 3*power2(
     logmstau12Q2)*(-7 - 3*invdmstau1mu*mstau12 - invdmstau2mu*(5*mstau12 + 6*
     mstau22) + mstau12*(mstau22 + mu2)*power2(invdmstau1mu) + mu2*(-5*mstau12
     - 6*mstau22 + mu2)*power2(invdmstau2mu) + 6*power2(sa)) + 6*logmH2Q2*(18 +
     9*invdmstau1mu*mstau12 + 9*invdmstau2mu*mstau12 + invdmstau2mu*
     logmstau22Q2*mstau12 + 6*invdmstau2mu*mstau22 + 2*invdmstau2mu*
     logmstau22Q2*mstau22 - 3*mstau12*(2*mstau22 + mu2)*power2(invdmstau1mu) +
     9*mstau12*mu2*power2(invdmstau2mu) + logmstau22Q2*mstau12*mu2*power2(
     invdmstau2mu) + 6*mstau22*mu2*power2(invdmstau2mu) + 2*logmstau22Q2*
     mstau22*mu2*power2(invdmstau2mu) - 18*power2(sa) + logmstau12Q2*(-6 - 3*
     invdmstau1mu*mstau12 - 4*invdmstau2mu*(mstau12 + mstau22) + mstau12*(2*
     mstau22 + mu2)*power2(invdmstau1mu) - 4*(mstau12 + mstau22)*mu2*power2(
     invdmstau2mu) + 6*power2(sa))) - 6*logmstau12Q2*(logmstau22Q2*(1 +
     invdmstau2mu*(mstau12 + 2*mstau22) + mstau12*mstau22*power2(invdmstau1mu)
     + (mstau12 + 2*mstau22 - mu2)*mu2*power2(invdmstau2mu)) + 3*(-7 - 3*
     invdmstau1mu*mstau12 - invdmstau2mu*(5*mstau12 + 6*mstau22) + mstau12*(
     mstau22 + mu2)*power2(invdmstau1mu) + mu2*(-5*mstau12 - 6*mstau22 + mu2)*
     power2(invdmstau2mu) + 6*power2(sa))))*power3(mH2) + (-(mstau12*power2(
     invdmstau1mu)*(42 + 6*logmH2Q2*(-3 + logmstau12Q2) - 18*logmstau12Q2 + 3*
     power2(logmH2Q2) + 3*power2(logmstau12Q2) + power2(Pi))) + invdmstau2mu*(6
     *logmH2Q2*(-3 + 2*logmstau12Q2 - logmstau22Q2) + 6*logmstau12Q2*(-9 +
     logmstau22Q2) + 3*power2(logmH2Q2) + 9*power2(logmstau12Q2) + 2*(42 +
     power2(Pi))) + mu2*power2(invdmstau2mu)*(6*logmH2Q2*(-3 + 2*logmstau12Q2 -
     logmstau22Q2) + 6*logmstau12Q2*(-9 + logmstau22Q2) + 3*power2(logmH2Q2) +
     9*power2(logmstau12Q2) + 2*(42 + power2(Pi))))*power4(mH2) + 6*(-1 +
     power2(sa))*(-(mstau12*mstau24*(-6*(9 + logmH2Q2)*logmstau22Q2 + 6*
     logmstau12Q2*(-15 + logmH2Q2 + 4*logmstau22Q2) + 15*power2(logmstau12Q2) +
     9*power2(logmstau22Q2) + 4*(42 + power2(Pi)))) + (-6*(3 + logmH2Q2)*
     logmstau22Q2 + 6*logmstau12Q2*(-9 + logmH2Q2 + 2*logmstau22Q2) + 9*power2(
     logmstau12Q2) + 3*power2(logmstau22Q2) + 2*(42 + power2(Pi)))*power6(
     mstau2)) + mH2*(-(invdmstau1mu*mstau12*(3*mstau24 - 2*mstau22*mu2)*(42 + 6
     *logmstau12Q2*(-3 + logmstau22Q2) - 18*logmstau22Q2 + 3*power2(
     logmstau12Q2) + 3*power2(logmstau22Q2) + power2(Pi))) + invdmstau2mu*
     mstau12*(-1 + invdmstau2mu*mu2)*power2(mu2)*(42 + 6*logmstau12Q2*(-3 +
     logmstau22Q2) - 18*logmstau22Q2 + 3*power2(logmstau12Q2) + 3*power2(
     logmstau22Q2) + power2(Pi)) + mstau24*(714 - 198*logmstau22Q2 + 33*power2(
     logmstau22Q2) + 42*power2(invdmstau2mu)*power2(mu2) - 18*logmstau22Q2*
     power2(invdmstau2mu)*power2(mu2) + 3*power2(invdmstau2mu)*power2(
     logmstau22Q2)*power2(mu2) + 17*power2(Pi) + power2(invdmstau2mu)*power2(
     mu2)*power2(Pi) + 3*power2(logmstau12Q2)*(29 + power2(invdmstau2mu)*power2
     (mu2) - 30*power2(sa)) - 36*logmH2Q2*(3 + logmstau12Q2 - 2*logmstau22Q2)*(
     -1 + power2(sa)) + 18*power2(logmH2Q2)*(-1 + power2(sa)) - 756*power2(sa)
     + 216*logmstau22Q2*power2(sa) - 36*power2(logmstau22Q2)*power2(sa) - 18*
     power2(Pi)*power2(sa) + 6*logmstau12Q2*(-87 - 3*power2(invdmstau2mu)*
     power2(mu2) + logmstau22Q2*(23 + power2(invdmstau2mu)*power2(mu2) - 24*
     power2(sa)) + 90*power2(sa))) - mstau12*mstau22*(-1722 + 90*logmstau22Q2 -
     15*power2(logmstau22Q2) + 84*power2(invdmstau2mu)*power2(mu2) - 36*
     logmstau22Q2*power2(invdmstau2mu)*power2(mu2) + 6*power2(invdmstau2mu)*
     power2(logmstau22Q2)*power2(mu2) - 41*power2(Pi) + 2*power2(invdmstau2mu)*
     power2(mu2)*power2(Pi) + 108*logmH2Q2*(-3 + 2*logmstau12Q2 - logmstau22Q2)
     *(-1 + power2(sa)) + 54*power2(logmH2Q2)*(-1 + power2(sa)) + 1764*power2(
     sa) - 108*logmstau22Q2*power2(sa) + 18*power2(logmstau22Q2)*power2(sa) +
     42*power2(Pi)*power2(sa) + 3*power2(logmstau12Q2)*(-59 + 2*power2(
     invdmstau2mu)*power2(mu2) + 60*power2(sa)) + 6*logmstau12Q2*(177 - 23*
     logmstau22Q2 - 6*power2(invdmstau2mu)*power2(mu2) + 2*logmstau22Q2*power2(
     invdmstau2mu)*power2(mu2) - 180*power2(sa) + 24*logmstau22Q2*power2(sa)))
     - mstau12*power2(invdmstau1mu)*(42 + 6*logmstau12Q2*(-3 + logmstau22Q2) -
     18*logmstau22Q2 + 3*power2(logmstau12Q2) + 3*power2(logmstau22Q2) + power2
     (Pi))*(-(mstau24*mu2) + power6(mstau2)))))/(48.*mH2*mstau12) + (mu2*power2
     (sa)*(power2(mH2)*(invdmstau2mu*mstau24*(1 + invdmstau2mu*mu2)*(6*logmH2Q2
     *(-3 + 2*logmstau12Q2 - logmstau22Q2) + 6*logmstau12Q2*(-9 + logmstau22Q2)
     + 3*power2(logmH2Q2) + 9*power2(logmstau12Q2) + 2*(42 + power2(Pi))) - 2*
     mstau22*(-42 - 36*logmstau22Q2 + 6*power2(logmstau22Q2) + 42*power2(
     invdmstau2mu)*power2(mu2) - 18*logmstau22Q2*power2(invdmstau2mu)*power2(
     mu2) + 3*power2(invdmstau2mu)*power2(logmstau22Q2)*power2(mu2) - power2(Pi
     ) + power2(invdmstau2mu)*power2(mu2)*power2(Pi) + 3*power2(logmstau12Q2)*(
     2 + power2(invdmstau2mu)*power2(mu2) - 3*power2(sa)) + 18*logmH2Q2*(-6 +
     logmstau12Q2 + logmstau22Q2)*(-1 + power2(sa)) + 18*power2(logmH2Q2)*(-1 +
     power2(sa)) + 54*logmstau22Q2*power2(sa) - 9*power2(logmstau22Q2)*power2(
     sa) + 6*logmstau12Q2*(-6 - 3*power2(invdmstau2mu)*power2(mu2) +
     logmstau22Q2*(5 + power2(invdmstau2mu)*power2(mu2) - 6*power2(sa)) + 9*
     power2(sa))) + mstau12*(294 - 54*logmstau22Q2 - 168*invdmstau1mu*mstau22 -
     168*invdmstau2mu*mstau22 + 54*invdmstau1mu*logmstau22Q2*mstau22 + 84*
     invdmstau1mu*mu2 + 42*mstau24*power2(invdmstau1mu) - 36*logmstau22Q2*
     mstau24*power2(invdmstau1mu) + 168*mstau22*mu2*power2(invdmstau1mu) - 18*
     logmstau22Q2*mstau22*mu2*power2(invdmstau1mu) - 168*mstau22*mu2*power2(
     invdmstau2mu) + 9*power2(logmstau22Q2) - 9*invdmstau1mu*mstau22*power2(
     logmstau22Q2) + 6*mstau24*power2(invdmstau1mu)*power2(logmstau22Q2) + 3*
     mstau22*mu2*power2(invdmstau1mu)*power2(logmstau22Q2) - 126*power2(
     invdmstau2mu)*power2(mu2) + 54*logmstau22Q2*power2(invdmstau2mu)*power2(
     mu2) - 9*power2(invdmstau2mu)*power2(logmstau22Q2)*power2(mu2) + 7*power2(
     Pi) - 4*invdmstau1mu*mstau22*power2(Pi) - 4*invdmstau2mu*mstau22*power2(Pi
     ) + 2*invdmstau1mu*mu2*power2(Pi) + mstau24*power2(invdmstau1mu)*power2(Pi
     ) + 4*mstau22*mu2*power2(invdmstau1mu)*power2(Pi) - 4*mstau22*mu2*power2(
     invdmstau2mu)*power2(Pi) - 3*power2(invdmstau2mu)*power2(mu2)*power2(Pi) +
     3*power2(logmH2Q2)*(3 - 2*invdmstau2mu*mstau22 + invdmstau1mu*(-3*mstau22
     + 2*mu2) + (-mstau24 + mstau22*mu2)*power2(invdmstau1mu) + mu2*(-2*mstau22
      + mu2)*power2(invdmstau2mu) - 6*power2(sa)) + 3*power2(logmstau12Q2)*(8 -
     6*invdmstau2mu*mstau22 - 2*invdmstau1mu*(mstau22 - mu2) + (mstau24 + 6*
     mstau22*mu2)*power2(invdmstau1mu) - 2*mu2*(3*mstau22 + 2*mu2)*power2(
     invdmstau2mu) - 6*power2(sa)) - 252*power2(sa) - 6*power2(Pi)*power2(sa) +
     6*logmstau12Q2*(6*invdmstau1mu*(mstau22 - mu2) - 3*(mstau24 + 6*mstau22*
     mu2)*power2(invdmstau1mu) + logmstau22Q2*(-(invdmstau1mu*mstau22) - 2*(1 +
     invdmstau2mu*mu2)*(-2 + invdmstau2mu*mstau22 + 2*invdmstau2mu*mu2) + (2*
     mstau24 + 3*mstau22*mu2)*power2(invdmstau1mu)) + 6*(-4 + 3*invdmstau2mu*
     mstau22 + mu2*(3*mstau22 + 2*mu2)*power2(invdmstau2mu) + 3*power2(sa))) +
     6*logmH2Q2*(-9 + 9*invdmstau1mu*mstau22 + 6*invdmstau2mu*mstau22 - 6*
     invdmstau1mu*mu2 + 3*mstau24*power2(invdmstau1mu) - 3*mstau22*mu2*power2(
     invdmstau1mu) + logmstau22Q2*(-2*invdmstau1mu*mstau22 + (1 + invdmstau2mu*
     mu2)*(-1 + 2*invdmstau2mu*mstau22 + invdmstau2mu*mu2) - 2*mstau22*mu2*
     power2(invdmstau1mu)) + 6*mstau22*mu2*power2(invdmstau2mu) - 3*power2(
     invdmstau2mu)*power2(mu2) + 18*power2(sa) - logmstau12Q2*(-4 + 4*
     invdmstau2mu*mstau22 + invdmstau1mu*(mstau22 - 2*mu2) + (mstau24 - 3*
     mstau22*mu2)*power2(invdmstau1mu) + 4*mstau22*mu2*power2(invdmstau2mu) + 6
     *power2(sa))))) + (-294 + 18*logmstau22Q2 - 126*invdmstau1mu*mstau12 - 168
     *invdmstau2mu*mstau12 - 168*invdmstau2mu*mstau22 + 42*mstau12*mstau22*
     power2(invdmstau1mu) + 18*logmstau22Q2*mstau12*mstau22*power2(invdmstau1mu
     ) + 42*mstau12*mu2*power2(invdmstau1mu) - 168*mstau12*mu2*power2(
     invdmstau2mu) - 168*mstau22*mu2*power2(invdmstau2mu) - 3*power2(
     logmstau22Q2) - 3*mstau12*mstau22*power2(invdmstau1mu)*power2(logmstau22Q2
     ) + 42*power2(invdmstau2mu)*power2(mu2) - 18*logmstau22Q2*power2(
     invdmstau2mu)*power2(mu2) + 3*power2(invdmstau2mu)*power2(logmstau22Q2)*
     power2(mu2) - 7*power2(Pi) - 3*invdmstau1mu*mstau12*power2(Pi) - 4*
     invdmstau2mu*mstau12*power2(Pi) - 4*invdmstau2mu*mstau22*power2(Pi) +
     mstau12*mstau22*power2(invdmstau1mu)*power2(Pi) + mstau12*mu2*power2(
     invdmstau1mu)*power2(Pi) - 4*mstau12*mu2*power2(invdmstau2mu)*power2(Pi) -
     4*mstau22*mu2*power2(invdmstau2mu)*power2(Pi) + power2(invdmstau2mu)*
     power2(mu2)*power2(Pi) + 3*power2(logmH2Q2)*(-3*invdmstau1mu*mstau12 -
     invdmstau2mu*(3*mstau12 + 2*mstau22) + mstau12*(2*mstau22 + mu2)*power2(
     invdmstau1mu) - (3*mstau12 + 2*mstau22)*mu2*power2(invdmstau2mu) + 6*(-1 +
     power2(sa))) + 252*power2(sa) + 6*power2(Pi)*power2(sa) + 3*power2(
     logmstau12Q2)*(-7 - 3*invdmstau1mu*mstau12 - invdmstau2mu*(5*mstau12 + 6*
     mstau22) + mstau12*(mstau22 + mu2)*power2(invdmstau1mu) + mu2*(-5*mstau12
     - 6*mstau22 + mu2)*power2(invdmstau2mu) + 6*power2(sa)) + 6*logmH2Q2*(18 +
     9*invdmstau1mu*mstau12 + 9*invdmstau2mu*mstau12 + invdmstau2mu*
     logmstau22Q2*mstau12 + 6*invdmstau2mu*mstau22 + 2*invdmstau2mu*
     logmstau22Q2*mstau22 - 3*mstau12*(2*mstau22 + mu2)*power2(invdmstau1mu) +
     9*mstau12*mu2*power2(invdmstau2mu) + logmstau22Q2*mstau12*mu2*power2(
     invdmstau2mu) + 6*mstau22*mu2*power2(invdmstau2mu) + 2*logmstau22Q2*
     mstau22*mu2*power2(invdmstau2mu) - 18*power2(sa) + logmstau12Q2*(-6 - 3*
     invdmstau1mu*mstau12 - 4*invdmstau2mu*(mstau12 + mstau22) + mstau12*(2*
     mstau22 + mu2)*power2(invdmstau1mu) - 4*(mstau12 + mstau22)*mu2*power2(
     invdmstau2mu) + 6*power2(sa))) - 6*logmstau12Q2*(logmstau22Q2*(1 +
     invdmstau2mu*(mstau12 + 2*mstau22) + mstau12*mstau22*power2(invdmstau1mu)
     + (mstau12 + 2*mstau22 - mu2)*mu2*power2(invdmstau2mu)) + 3*(-7 - 3*
     invdmstau1mu*mstau12 - invdmstau2mu*(5*mstau12 + 6*mstau22) + mstau12*(
     mstau22 + mu2)*power2(invdmstau1mu) + mu2*(-5*mstau12 - 6*mstau22 + mu2)*
     power2(invdmstau2mu) + 6*power2(sa))))*power3(mH2) + (-(mstau12*power2(
     invdmstau1mu)*(42 + 6*logmH2Q2*(-3 + logmstau12Q2) - 18*logmstau12Q2 + 3*
     power2(logmH2Q2) + 3*power2(logmstau12Q2) + power2(Pi))) + invdmstau2mu*(6
     *logmH2Q2*(-3 + 2*logmstau12Q2 - logmstau22Q2) + 6*logmstau12Q2*(-9 +
     logmstau22Q2) + 3*power2(logmH2Q2) + 9*power2(logmstau12Q2) + 2*(42 +
     power2(Pi))) + mu2*power2(invdmstau2mu)*(6*logmH2Q2*(-3 + 2*logmstau12Q2 -
     logmstau22Q2) + 6*logmstau12Q2*(-9 + logmstau22Q2) + 3*power2(logmH2Q2) +
     9*power2(logmstau12Q2) + 2*(42 + power2(Pi))))*power4(mH2) + 6*(-1 +
     power2(sa))*(-(mstau12*mstau24*(-6*(9 + logmH2Q2)*logmstau22Q2 + 6*
     logmstau12Q2*(-15 + logmH2Q2 + 4*logmstau22Q2) + 15*power2(logmstau12Q2) +
     9*power2(logmstau22Q2) + 4*(42 + power2(Pi)))) + (-6*(3 + logmH2Q2)*
     logmstau22Q2 + 6*logmstau12Q2*(-9 + logmH2Q2 + 2*logmstau22Q2) + 9*power2(
     logmstau12Q2) + 3*power2(logmstau22Q2) + 2*(42 + power2(Pi)))*power6(
     mstau2)) + mH2*(-(invdmstau1mu*mstau12*(3*mstau24 - 2*mstau22*mu2)*(42 + 6
     *logmstau12Q2*(-3 + logmstau22Q2) - 18*logmstau22Q2 + 3*power2(
     logmstau12Q2) + 3*power2(logmstau22Q2) + power2(Pi))) + invdmstau2mu*
     mstau12*(-1 + invdmstau2mu*mu2)*power2(mu2)*(42 + 6*logmstau12Q2*(-3 +
     logmstau22Q2) - 18*logmstau22Q2 + 3*power2(logmstau12Q2) + 3*power2(
     logmstau22Q2) + power2(Pi)) + mstau24*(714 - 198*logmstau22Q2 + 33*power2(
     logmstau22Q2) + 42*power2(invdmstau2mu)*power2(mu2) - 18*logmstau22Q2*
     power2(invdmstau2mu)*power2(mu2) + 3*power2(invdmstau2mu)*power2(
     logmstau22Q2)*power2(mu2) + 17*power2(Pi) + power2(invdmstau2mu)*power2(
     mu2)*power2(Pi) + 3*power2(logmstau12Q2)*(29 + power2(invdmstau2mu)*power2
     (mu2) - 30*power2(sa)) - 36*logmH2Q2*(3 + logmstau12Q2 - 2*logmstau22Q2)*(
     -1 + power2(sa)) + 18*power2(logmH2Q2)*(-1 + power2(sa)) - 756*power2(sa)
     + 216*logmstau22Q2*power2(sa) - 36*power2(logmstau22Q2)*power2(sa) - 18*
     power2(Pi)*power2(sa) + 6*logmstau12Q2*(-87 - 3*power2(invdmstau2mu)*
     power2(mu2) + logmstau22Q2*(23 + power2(invdmstau2mu)*power2(mu2) - 24*
     power2(sa)) + 90*power2(sa))) - mstau12*mstau22*(-1722 + 90*logmstau22Q2 -
     15*power2(logmstau22Q2) + 84*power2(invdmstau2mu)*power2(mu2) - 36*
     logmstau22Q2*power2(invdmstau2mu)*power2(mu2) + 6*power2(invdmstau2mu)*
     power2(logmstau22Q2)*power2(mu2) - 41*power2(Pi) + 2*power2(invdmstau2mu)*
     power2(mu2)*power2(Pi) + 108*logmH2Q2*(-3 + 2*logmstau12Q2 - logmstau22Q2)
     *(-1 + power2(sa)) + 54*power2(logmH2Q2)*(-1 + power2(sa)) + 1764*power2(
     sa) - 108*logmstau22Q2*power2(sa) + 18*power2(logmstau22Q2)*power2(sa) +
     42*power2(Pi)*power2(sa) + 3*power2(logmstau12Q2)*(-59 + 2*power2(
     invdmstau2mu)*power2(mu2) + 60*power2(sa)) + 6*logmstau12Q2*(177 - 23*
     logmstau22Q2 - 6*power2(invdmstau2mu)*power2(mu2) + 2*logmstau22Q2*power2(
     invdmstau2mu)*power2(mu2) - 180*power2(sa) + 24*logmstau22Q2*power2(sa)))
     - mstau12*power2(invdmstau1mu)*(42 + 6*logmstau12Q2*(-3 + logmstau22Q2) -
     18*logmstau22Q2 + 3*power2(logmstau12Q2) + 3*power2(logmstau22Q2) + power2
     (Pi))*(-(mstau24*mu2) + power6(mstau2)))))/(48.*mH2*mstau12) + (Al*ca*mu*
     sa*(power2(mH2)*(-(invdmstau2mu*mstau24*(1 + invdmstau2mu*mu2)*(6*logmH2Q2
     *(-3 + 2*logmstau12Q2 - logmstau22Q2) + 6*logmstau12Q2*(-9 + logmstau22Q2)
     + 3*power2(logmH2Q2) + 9*power2(logmstau12Q2) + 2*(42 + power2(Pi)))) + 2*
     mstau22*(-42 - 36*logmstau22Q2 + 6*power2(logmstau22Q2) + 42*power2(
     invdmstau2mu)*power2(mu2) - 18*logmstau22Q2*power2(invdmstau2mu)*power2(
     mu2) + 3*power2(invdmstau2mu)*power2(logmstau22Q2)*power2(mu2) - power2(Pi
     ) + power2(invdmstau2mu)*power2(mu2)*power2(Pi) + 3*power2(logmstau12Q2)*(
     2 + power2(invdmstau2mu)*power2(mu2) - 3*power2(sa)) + 18*logmH2Q2*(-6 +
     logmstau12Q2 + logmstau22Q2)*(-1 + power2(sa)) + 18*power2(logmH2Q2)*(-1 +
     power2(sa)) + 54*logmstau22Q2*power2(sa) - 9*power2(logmstau22Q2)*power2(
     sa) + 6*logmstau12Q2*(-6 - 3*power2(invdmstau2mu)*power2(mu2) +
     logmstau22Q2*(5 + power2(invdmstau2mu)*power2(mu2) - 6*power2(sa)) + 9*
     power2(sa))) + mstau12*(-294 + 54*logmstau22Q2 + 168*invdmstau1mu*mstau22
     + 168*invdmstau2mu*mstau22 - 54*invdmstau1mu*logmstau22Q2*mstau22 - 84*
     invdmstau1mu*mu2 - 42*mstau24*power2(invdmstau1mu) + 36*logmstau22Q2*
     mstau24*power2(invdmstau1mu) - 168*mstau22*mu2*power2(invdmstau1mu) + 18*
     logmstau22Q2*mstau22*mu2*power2(invdmstau1mu) + 168*mstau22*mu2*power2(
     invdmstau2mu) - 9*power2(logmstau22Q2) + 9*invdmstau1mu*mstau22*power2(
     logmstau22Q2) - 6*mstau24*power2(invdmstau1mu)*power2(logmstau22Q2) - 3*
     mstau22*mu2*power2(invdmstau1mu)*power2(logmstau22Q2) + 126*power2(
     invdmstau2mu)*power2(mu2) - 54*logmstau22Q2*power2(invdmstau2mu)*power2(
     mu2) + 9*power2(invdmstau2mu)*power2(logmstau22Q2)*power2(mu2) - 7*power2(
     Pi) + 4*invdmstau1mu*mstau22*power2(Pi) + 4*invdmstau2mu*mstau22*power2(Pi
     ) - 2*invdmstau1mu*mu2*power2(Pi) - mstau24*power2(invdmstau1mu)*power2(Pi
     ) - 4*mstau22*mu2*power2(invdmstau1mu)*power2(Pi) + 4*mstau22*mu2*power2(
     invdmstau2mu)*power2(Pi) + 3*power2(invdmstau2mu)*power2(mu2)*power2(Pi) -
     3*power2(logmstau12Q2)*(8 - 6*invdmstau2mu*mstau22 - 2*invdmstau1mu*(
     mstau22 - mu2) + (mstau24 + 6*mstau22*mu2)*power2(invdmstau1mu) - 2*mu2*(3
     *mstau22 + 2*mu2)*power2(invdmstau2mu) - 6*power2(sa)) + 252*power2(sa) +
     6*power2(Pi)*power2(sa) + 3*power2(logmH2Q2)*(-3 + 2*invdmstau2mu*mstau22
     + invdmstau1mu*(3*mstau22 - 2*mu2) + (mstau24 - mstau22*mu2)*power2(
     invdmstau1mu) + (2*mstau22 - mu2)*mu2*power2(invdmstau2mu) + 6*power2(sa))
     - 6*logmstau12Q2*(6*invdmstau1mu*(mstau22 - mu2) - 3*(mstau24 + 6*mstau22*
     mu2)*power2(invdmstau1mu) + logmstau22Q2*(-(invdmstau1mu*mstau22) - 2*(1 +
     invdmstau2mu*mu2)*(-2 + invdmstau2mu*mstau22 + 2*invdmstau2mu*mu2) + (2*
     mstau24 + 3*mstau22*mu2)*power2(invdmstau1mu)) + 6*(-4 + 3*invdmstau2mu*
     mstau22 + mu2*(3*mstau22 + 2*mu2)*power2(invdmstau2mu) + 3*power2(sa))) -
     6*logmH2Q2*(-9 + 9*invdmstau1mu*mstau22 + 6*invdmstau2mu*mstau22 - 6*
     invdmstau1mu*mu2 + 3*mstau24*power2(invdmstau1mu) - 3*mstau22*mu2*power2(
     invdmstau1mu) + logmstau22Q2*(-2*invdmstau1mu*mstau22 + (1 + invdmstau2mu*
     mu2)*(-1 + 2*invdmstau2mu*mstau22 + invdmstau2mu*mu2) - 2*mstau22*mu2*
     power2(invdmstau1mu)) + 6*mstau22*mu2*power2(invdmstau2mu) - 3*power2(
     invdmstau2mu)*power2(mu2) + 18*power2(sa) - logmstau12Q2*(-4 + 4*
     invdmstau2mu*mstau22 + invdmstau1mu*(mstau22 - 2*mu2) + (mstau24 - 3*
     mstau22*mu2)*power2(invdmstau1mu) + 4*mstau22*mu2*power2(invdmstau2mu) + 6
     *power2(sa))))) + (294 - 18*logmstau22Q2 + 126*invdmstau1mu*mstau12 + 168*
     invdmstau2mu*mstau12 + 168*invdmstau2mu*mstau22 - 42*mstau12*mstau22*
     power2(invdmstau1mu) - 18*logmstau22Q2*mstau12*mstau22*power2(invdmstau1mu
     ) - 42*mstau12*mu2*power2(invdmstau1mu) + 168*mstau12*mu2*power2(
     invdmstau2mu) + 168*mstau22*mu2*power2(invdmstau2mu) + 3*power2(
     logmstau22Q2) + 3*mstau12*mstau22*power2(invdmstau1mu)*power2(logmstau22Q2
     ) - 42*power2(invdmstau2mu)*power2(mu2) + 18*logmstau22Q2*power2(
     invdmstau2mu)*power2(mu2) - 3*power2(invdmstau2mu)*power2(logmstau22Q2)*
     power2(mu2) + 7*power2(Pi) + 3*invdmstau1mu*mstau12*power2(Pi) + 4*
     invdmstau2mu*mstau12*power2(Pi) + 4*invdmstau2mu*mstau22*power2(Pi) -
     mstau12*mstau22*power2(invdmstau1mu)*power2(Pi) - mstau12*mu2*power2(
     invdmstau1mu)*power2(Pi) + 4*mstau12*mu2*power2(invdmstau2mu)*power2(Pi) +
     4*mstau22*mu2*power2(invdmstau2mu)*power2(Pi) - power2(invdmstau2mu)*
     power2(mu2)*power2(Pi) + 3*power2(logmH2Q2)*(6 + 3*invdmstau1mu*mstau12 +
     invdmstau2mu*(3*mstau12 + 2*mstau22) - mstau12*(2*mstau22 + mu2)*power2(
     invdmstau1mu) + (3*mstau12 + 2*mstau22)*mu2*power2(invdmstau2mu) - 6*
     power2(sa)) - 252*power2(sa) - 6*power2(Pi)*power2(sa) - 3*power2(
     logmstau12Q2)*(-7 - 3*invdmstau1mu*mstau12 - invdmstau2mu*(5*mstau12 + 6*
     mstau22) + mstau12*(mstau22 + mu2)*power2(invdmstau1mu) + mu2*(-5*mstau12
     - 6*mstau22 + mu2)*power2(invdmstau2mu) + 6*power2(sa)) - 6*logmH2Q2*(18 +
     9*invdmstau1mu*mstau12 + 9*invdmstau2mu*mstau12 + invdmstau2mu*
     logmstau22Q2*mstau12 + 6*invdmstau2mu*mstau22 + 2*invdmstau2mu*
     logmstau22Q2*mstau22 - 3*mstau12*(2*mstau22 + mu2)*power2(invdmstau1mu) +
     9*mstau12*mu2*power2(invdmstau2mu) + logmstau22Q2*mstau12*mu2*power2(
     invdmstau2mu) + 6*mstau22*mu2*power2(invdmstau2mu) + 2*logmstau22Q2*
     mstau22*mu2*power2(invdmstau2mu) - 18*power2(sa) + logmstau12Q2*(-6 - 3*
     invdmstau1mu*mstau12 - 4*invdmstau2mu*(mstau12 + mstau22) + mstau12*(2*
     mstau22 + mu2)*power2(invdmstau1mu) - 4*(mstau12 + mstau22)*mu2*power2(
     invdmstau2mu) + 6*power2(sa))) + 6*logmstau12Q2*(logmstau22Q2*(1 +
     invdmstau2mu*(mstau12 + 2*mstau22) + mstau12*mstau22*power2(invdmstau1mu)
     + (mstau12 + 2*mstau22 - mu2)*mu2*power2(invdmstau2mu)) + 3*(-7 - 3*
     invdmstau1mu*mstau12 - invdmstau2mu*(5*mstau12 + 6*mstau22) + mstau12*(
     mstau22 + mu2)*power2(invdmstau1mu) + mu2*(-5*mstau12 - 6*mstau22 + mu2)*
     power2(invdmstau2mu) + 6*power2(sa))))*power3(mH2) - (-(mstau12*power2(
     invdmstau1mu)*(42 + 6*logmH2Q2*(-3 + logmstau12Q2) - 18*logmstau12Q2 + 3*
     power2(logmH2Q2) + 3*power2(logmstau12Q2) + power2(Pi))) + invdmstau2mu*(6
     *logmH2Q2*(-3 + 2*logmstau12Q2 - logmstau22Q2) + 6*logmstau12Q2*(-9 +
     logmstau22Q2) + 3*power2(logmH2Q2) + 9*power2(logmstau12Q2) + 2*(42 +
     power2(Pi))) + mu2*power2(invdmstau2mu)*(6*logmH2Q2*(-3 + 2*logmstau12Q2 -
     logmstau22Q2) + 6*logmstau12Q2*(-9 + logmstau22Q2) + 3*power2(logmH2Q2) +
     9*power2(logmstau12Q2) + 2*(42 + power2(Pi))))*power4(mH2) - 6*(-1 +
     power2(sa))*(-(mstau12*mstau24*(-6*(9 + logmH2Q2)*logmstau22Q2 + 6*
     logmstau12Q2*(-15 + logmH2Q2 + 4*logmstau22Q2) + 15*power2(logmstau12Q2) +
     9*power2(logmstau22Q2) + 4*(42 + power2(Pi)))) + (-6*(3 + logmH2Q2)*
     logmstau22Q2 + 6*logmstau12Q2*(-9 + logmH2Q2 + 2*logmstau22Q2) + 9*power2(
     logmstau12Q2) + 3*power2(logmstau22Q2) + 2*(42 + power2(Pi)))*power6(
     mstau2)) + mH2*(invdmstau1mu*mstau12*(3*mstau24 - 2*mstau22*mu2)*(42 + 6*
     logmstau12Q2*(-3 + logmstau22Q2) - 18*logmstau22Q2 + 3*power2(logmstau12Q2
     ) + 3*power2(logmstau22Q2) + power2(Pi)) - mstau24*(714 - 198*logmstau22Q2
      + 33*power2(logmstau22Q2) + 42*power2(invdmstau2mu)*power2(mu2) - 18*
     logmstau22Q2*power2(invdmstau2mu)*power2(mu2) + 3*power2(invdmstau2mu)*
     power2(logmstau22Q2)*power2(mu2) + 17*power2(Pi) + power2(invdmstau2mu)*
     power2(mu2)*power2(Pi) + 3*power2(logmstau12Q2)*(29 + power2(invdmstau2mu)
     *power2(mu2) - 30*power2(sa)) - 36*logmH2Q2*(3 + logmstau12Q2 - 2*
     logmstau22Q2)*(-1 + power2(sa)) + 18*power2(logmH2Q2)*(-1 + power2(sa)) -
     756*power2(sa) + 216*logmstau22Q2*power2(sa) - 36*power2(logmstau22Q2)*
     power2(sa) - 18*power2(Pi)*power2(sa) + 6*logmstau12Q2*(-87 - 3*power2(
     invdmstau2mu)*power2(mu2) + logmstau22Q2*(23 + power2(invdmstau2mu)*power2
     (mu2) - 24*power2(sa)) + 90*power2(sa))) + mstau12*(-(invdmstau2mu*(-1 +
     invdmstau2mu*mu2)*power2(mu2)*(42 + 6*logmstau12Q2*(-3 + logmstau22Q2) -
     18*logmstau22Q2 + 3*power2(logmstau12Q2) + 3*power2(logmstau22Q2) + power2
     (Pi))) + mstau22*(-1722 + 90*logmstau22Q2 - 15*power2(logmstau22Q2) + 84*
     power2(invdmstau2mu)*power2(mu2) - 36*logmstau22Q2*power2(invdmstau2mu)*
     power2(mu2) + 6*power2(invdmstau2mu)*power2(logmstau22Q2)*power2(mu2) - 41
     *power2(Pi) + 2*power2(invdmstau2mu)*power2(mu2)*power2(Pi) + 108*logmH2Q2
     *(-3 + 2*logmstau12Q2 - logmstau22Q2)*(-1 + power2(sa)) + 54*power2(
     logmH2Q2)*(-1 + power2(sa)) + 1764*power2(sa) - 108*logmstau22Q2*power2(sa
     ) + 18*power2(logmstau22Q2)*power2(sa) + 42*power2(Pi)*power2(sa) + 3*
     power2(logmstau12Q2)*(-59 + 2*power2(invdmstau2mu)*power2(mu2) + 60*power2
     (sa)) + 6*logmstau12Q2*(177 - 23*logmstau22Q2 - 6*power2(invdmstau2mu)*
     power2(mu2) + 2*logmstau22Q2*power2(invdmstau2mu)*power2(mu2) - 180*power2
     (sa) + 24*logmstau22Q2*power2(sa)))) + mstau12*power2(invdmstau1mu)*(42 +
     6*logmstau12Q2*(-3 + logmstau22Q2) - 18*logmstau22Q2 + 3*power2(
     logmstau12Q2) + 3*power2(logmstau22Q2) + power2(Pi))*(-(mstau24*mu2) +
     power6(mstau2)))))/(24.*mH2*mstau12)) + Fin20(mstau22,mu2,Q2)*((invdmstau*
     power2(Al)*((invdmstau1mu*(mstau22 - mu2))/mA2 + (-1 + power2(sa))/mH2 -
     power2(sa)/mh2))/2. + (Al*invdmstau*mu*(ca*cb*mA2*(mh2 - mH2)*sa +
     invdmstau1mu*mh2*mH2*(-mstau22 + mu2)*sb + mA2*sb*(mh2 - mh2*power2(sa) +
     mH2*power2(sa))))/(2.*cb*mA2*mh2*mH2) + (6*invdmstau2mu - 8/mA2 + 16/mH2 +
     (8*invdmstau2mu*mu2)/mA2 - (8*invdmstau2mu*mu2)/mH2 + (8*ca*invdmstau*mu2*
     sa*sb)/(cb*mh2) - (8*ca*invdmstau*mu2*sa*sb)/(cb*mH2) + (5*mstau22)/power2
     (mA2) - (5*mu2)/power2(mA2) + (2*mstau22)/power2(mC2) - (2*mu2)/power2(mC2
     ) - (3*mstau22)/power2(mH2) + (3*mu2)/power2(mH2) + 2*invdmsntau2mu*(-1 +
     power2(invdmstau2mu)*(mstau24 - 2*mstau22*mu2 + power2(mu2))) + (16*power2
     (sa))/mh2 - (16*power2(sa))/mH2 - (8*invdmstau2mu*mu2*power2(sa))/mh2 + (8
     *invdmstau2mu*mu2*power2(sa))/mH2 - (3*mstau22*power2(sa))/power2(mh2) + (
     3*mu2*power2(sa))/power2(mh2) + (3*mstau22*power2(sa))/power2(mH2) - (3*
     mu2*power2(sa))/power2(mH2) - (8*invdmstau*invdmstau1mu*(mstau22 - mu2)*
     mu2*power2(sb))/mA2 - (2*invdmstau1mu*(mz2 + 4*invdmstau2mu*(mstau22 + mu2
     )*mz2 - (mstau24 + mu2*(-2*mstau22 + mu2))*mz2*power2(invdmstau2mu) + 4*
     invdmstau*mu2*(-mstau22 + mu2)*power2(sb)))/mz2 + 4*mu2*(-mstau22 + mu2)*
     power3(invdmsntau2mu) + 4*mu2*(-mstau22 + mu2)*power3(invdmstau1mu) + 2*
     power2(invdmsntau2mu)*(-mu2 + power2(invdmstau2mu)*(-4*mstau24*mu2 + 5*
     mstau22*power2(mu2) - 2*power3(mu2) + power6(mstau2))) + 2*power2(
     invdmstau1mu)*(-mu2 + power2(invdmstau2mu)*(-4*mstau24*mu2 + 5*mstau22*
     power2(mu2) - 2*power3(mu2) + power6(mstau2))))/16.);

   return result * power4(ytau) * twoLoop;
}

double delta_mtau_2loop_atau_at(const Parameters& pars)
{
   const Real ytau  = pars.ytau;
   const Real yt    = pars.yt;
   const Real xt    = pars.xt;
   const Real xtau  = pars.xtau;
   const Real mt    = pars.mt;
   const Real mt2   = power2(pars.mt);
   const Real mst1  = pars.mst1;
   const Real mst12 = power2(pars.mst1);
   const Real mst14 = power4(pars.mst1);
   const Real mst2  = pars.mst2;
   const Real mst22 = power2(pars.mst2);
   const Real mst24 = power4(pars.mst2);
   const Real msb1  = pars.msb1;
   const Real msb12 = power2(pars.msb1);
   const Real msb14 = power4(pars.msb1);
   const Real mstau12 = power2(pars.mstau1);
   const Real mstau22 = power2(pars.mstau2);
   const Real msntau2 = power2(pars.msntau);
   const Real mh2   = power2(pars.mh);
   const Real mH2   = power2(pars.mH);
   const Real mC2   = power2(pars.mC);
   const Real mA2   = power2(pars.mA);
   const Real mu    = pars.mu;
   const Real mu2   = power2(pars.mu);
   const Real tb    = pars.tb;
   const Real sb    = tb / std::sqrt(1 + power2(tb));
   const Real cb    = 1  / std::sqrt(1 + power2(tb));
   const Real Q2    = power2(pars.Q);
   const Real snt   = calc_sin_theta(mt, xt, mst12, mst22);
   const Real alpha = calc_alpha(mh2, mH2, tb);
   const Real sa    = std::sin(alpha);
   const Real ca    = std::cos(alpha);
   const Real At    = xt + mu/tb;
   const Real Al    = xtau + mu*tb;

   const Real invdmst       = 1/(-mst12 + mst22);
   const Real invdmsntaust1 = 1/(-msntau2 + mst12);
   const Real invdmsntaust2 = 1/(-msntau2 + mst22);
   const Real invdmstau     = 1/(mstau12 - mstau22);
   const Real invdmstau1mu  = 1/(-mstau12 + mu2);
   const Real invdmstau2mu  = 1/(-mstau22 + mu2);
   const Real invdmsntau2mu = 1/(-msntau2 + mu2);
   const Real invdmhH       = 1/(-mh2 + mH2);

   const Real logmstau12Q2  = std::log(mstau12/Q2);
   const Real logmstau22Q2  = std::log(mstau22/Q2);
   const Real logmH2Q2      = std::log(mH2/Q2);
   const Real logmA2Q2      = std::log(mA2/Q2);
   const Real logmh2Q2      = std::log(mh2/Q2);
   const Real logmu2Q2      = std::log(mu2/Q2);
   const Real logmC2Q2      = std::log(mC2/Q2);
   const Real logmsntau2Q2  = std::log(msntau2/Q2);
   const Real logmst12Q2    = std::log(mst12/Q2);
   const Real logmst22Q2    = std::log(mst22/Q2);
   const Real logmsb12Q2    = std::log(msb12/Q2);
   const Real logmt2Q2      = std::log(mt2/Q2);

   const double result =
   Fin3(msntau2,msb12,mt2,Q2)*((-3*mu2*(-1 + invdmsntau2mu*(-msb12 + mt2 + mu2)
     )*power2(invdmsntau2mu))/4. - (3*mu2*DeltaInv(msntau2,msb12,mt2)*(1 +
     invdmsntau2mu*(2*msb12 + mt2 - 2*mu2) + (msb14 + mu2*(-mt2 + mu2) - msb12*
     (mt2 + 2*mu2))*power2(invdmsntau2mu)))/4.) + Fin3(mA2,mst12,mst22,Q2)*((3*
     (mA2 - mst12 - mst22)*mu2*DeltaInv(mA2,mst12,mst22))/(4.*mA2) - (3*mu2)/(
     4.*power2(mA2))) - Al*At*invdmst*invdmstau*(invdmstau1mu - invdmstau2mu)*
     mu2*(mst12*(6*logmst12Q2*(-3 + 2*logmt2Q2 - logmu2Q2) + 6*logmt2Q2*(-9 +
     logmu2Q2) + 3*power2(logmst12Q2) + 9*power2(logmt2Q2) + 2*(42 + power2(Pi)
     )) - mst22*(6*logmst22Q2*(-3 + 2*logmt2Q2 - logmu2Q2) + 6*logmt2Q2*(-9 +
     logmu2Q2) + 3*power2(logmst22Q2) + 9*power2(logmt2Q2) + 2*(42 + power2(Pi)
     ))) + (mu2*DeltaInv(mA2,mst12,mst22)*(3*(logmst12Q2 - logmst22Q2)*(-6 + 2*
     logmA2Q2 + logmst12Q2 + logmst22Q2)*mst24 + power2(mA2)*(42 + 6*logmA2Q2*(
     -3 + logmst12Q2) - 18*logmst12Q2 + 3*power2(logmA2Q2) + 3*power2(
     logmst12Q2) + power2(Pi)) - mst12*mst22*(-6*(3 + logmA2Q2)*logmst22Q2 + 6*
     logmst12Q2*(-9 + logmA2Q2 + 2*logmst22Q2) + 9*power2(logmst12Q2) + 3*
     power2(logmst22Q2) + 2*(42 + power2(Pi))) - mA2*(mst12*(42 + 6*logmA2Q2*(-
     3 + logmst12Q2) - 18*logmst12Q2 + 3*power2(logmA2Q2) + 3*power2(logmst12Q2
     ) + power2(Pi)) + mst22*(42 - 36*logmst12Q2 + 6*logmA2Q2*(-3 + 2*
     logmst12Q2 - logmst22Q2) + 18*logmst22Q2 + 3*power2(logmA2Q2) + 6*power2(
     logmst12Q2) - 3*power2(logmst22Q2) + power2(Pi)))))/(8.*mA2) + (3*Fin3(mH2
     ,mt2,mt2,Q2)*(-2*mh2*mH2*(-3 + 4*invdmhH*mH2)*(mH2 - 4*mt2) + 3*(mH2 + 2*
     mt2)*power2(mh2) + 2*(-1 + invdmhH*mH2)*(mH2 - 4*mt2)*power2(mH2))*(-1 +
     power2(sa))*power2(sa))/(4.*power2(mh2)*power2(mH2)) + Fin3(mh2,mt2,mt2,Q2
     )*((3*(mh2*(3 - 20*invdmhH*mt2) + 2*mt2*(5 - 2*invdmhH*mH2 + 8*invdmhH*mt2
     ) + 6*invdmhH*power2(mh2))*(-1 + power2(sa))*power2(sa))/(4.*power2(mh2))
     + (12*(1 - invdmhH*mH2 + 4*invdmhH*mt2)*DeltaInv(mh2,mt2,mt2)*power2(mt2)*
     (-1 + power2(sa))*power2(sa))/mh2) + Fin3(mC2,mst22,msb12,Q2)*((3*(mC2 -
     msb12 - mst22)*mu2*DeltaInv(mC2,mst22,msb12)*(-1 + power2(snt)))/(4.*mC2)
     - (3*mu2*(-1 + power2(snt)))/(4.*power2(mC2))) - (3*mu2*(-1 +
     invdmsntau2mu*(-mst22 + mu2))*Fin20(mst22,mu2,Q2)*power2(invdmsntau2mu)*(-
     1 + power2(snt)))/4. + (3*invdmsntau2mu*mu2*Fin20(msntau2,mst22,Q2)*(
     invdmsntaust2 - 2*(mst22 - mu2)*power2(invdmsntau2mu) - mst22*power2(
     invdmsntaust2) + invdmsntau2mu*(-1 + invdmsntaust2*(2*mst22 - mu2) + (-
     mst24 + mst22*mu2)*power2(invdmsntaust2)))*(-1 + power2(snt)))/8. + (3*mu2
     *Fin20(mst22,msb12,Q2)*(-1 + power2(snt)))/(4.*power2(mC2)) + (mu2*
     DeltaInv(mC2,mst22,msb12)*(3*(logmsb12Q2 - logmst22Q2)*(-6 + 2*logmC2Q2 +
     logmsb12Q2 + logmst22Q2)*mst24 + power2(mC2)*(42 + 6*logmC2Q2*(-3 +
     logmsb12Q2) - 18*logmsb12Q2 + 3*power2(logmC2Q2) + 3*power2(logmsb12Q2) +
     power2(Pi)) - msb12*mst22*(-6*(3 + logmC2Q2)*logmst22Q2 + 6*logmsb12Q2*(-9
      + logmC2Q2 + 2*logmst22Q2) + 9*power2(logmsb12Q2) + 3*power2(logmst22Q2)
     + 2*(42 + power2(Pi))) - mC2*(msb12*(42 + 6*logmC2Q2*(-3 + logmsb12Q2) -
     18*logmsb12Q2 + 3*power2(logmC2Q2) + 3*power2(logmsb12Q2) + power2(Pi)) +
     mst22*(42 - 36*logmsb12Q2 + 6*logmC2Q2*(-3 + 2*logmsb12Q2 - logmst22Q2) +
     18*logmst22Q2 + 3*power2(logmC2Q2) + 6*power2(logmsb12Q2) - 3*power2(
     logmst22Q2) + power2(Pi))))*(-1 + power2(snt)))/(8.*mC2) + (3*mu2*(-1 +
     invdmsntau2mu*(-mst12 + mu2))*Fin20(mst12,mu2,Q2)*power2(invdmsntau2mu)*
     power2(snt))/4. - (3*invdmsntau2mu*mu2*Fin20(mst12,msntau2,Q2)*(
     invdmsntaust1 - 2*(mst12 - mu2)*power2(invdmsntau2mu) - mst12*power2(
     invdmsntaust1) + invdmsntau2mu*(-1 + invdmsntaust1*(2*mst12 - mu2) + (-
     mst14 + mst12*mu2)*power2(invdmsntaust1)))*power2(snt))/8. - (3*mu2*Fin20(
     mst12,msb12,Q2)*power2(snt))/(4.*power2(mC2)) + (mu2*DeltaInv(mst12,mC2,
     msb12)*(-3*(logmsb12Q2 - logmst12Q2)*(-6 + 2*logmC2Q2 + logmsb12Q2 +
     logmst12Q2)*mst14 - power2(mC2)*(42 + 6*logmC2Q2*(-3 + logmsb12Q2) - 18*
     logmsb12Q2 + 3*power2(logmC2Q2) + 3*power2(logmsb12Q2) + power2(Pi)) +
     msb12*mst12*(-6*(3 + logmC2Q2)*logmst12Q2 + 6*logmsb12Q2*(-9 + logmC2Q2 +
     2*logmst12Q2) + 9*power2(logmsb12Q2) + 3*power2(logmst12Q2) + 2*(42 +
     power2(Pi))) + mC2*(msb12*(42 + 6*logmC2Q2*(-3 + logmsb12Q2) - 18*
     logmsb12Q2 + 3*power2(logmC2Q2) + 3*power2(logmsb12Q2) + power2(Pi)) +
     mst12*(42 - 36*logmsb12Q2 + 6*logmC2Q2*(-3 + 2*logmsb12Q2 - logmst12Q2) +
     18*logmst12Q2 + 3*power2(logmC2Q2) + 6*power2(logmsb12Q2) - 3*power2(
     logmst12Q2) + power2(Pi))))*power2(snt))/(8.*mC2) + Fin3(mst12,mC2,msb12,
     Q2)*((3*(-mC2 + msb12 + mst12)*mu2*DeltaInv(mst12,mC2,msb12)*power2(snt))/
     (4.*mC2) + (3*mu2*power2(snt))/(4.*power2(mC2))) + Fin20(mst12,mst22,Q2)*(
     (9*power2(At)*power2(mh2 - mH2)*(-1 + power2(sa))*power2(sa)*power2(1 - 2*
     power2(snt)))/(4.*power2(mh2)*power2(mH2)) - (9*At*ca*mu*sa*(power2(mh2)*(
     -1 + power2(sa)) + power2(mH2)*power2(sa) + mh2*(mH2 - 2*mH2*power2(sa)))*
     power2(1 - 2*power2(snt)))/(2.*power2(mh2)*power2(mH2)) - (3*mu2*(-(power2
     (mh2)*power2(mH2)) + 3*power2(mA2)*power2(mh2 - mh2*power2(sa) + mH2*
     power2(sa))*power2(1 - 2*power2(snt))))/(4.*power2(mA2)*power2(mh2)*power2
     (mH2))) + Fin3(mst12,mH2,mst22,Q2)*((3*mu2*(-1 + power2(sa))*(3*power2(mh2
     )*(-1 + power2(sa)) + 2*mh2*mH2*(-3 + 4*invdmhH*mH2)*power2(sa) - 2*(-1 +
     invdmhH*mH2)*power2(mH2)*power2(sa))*power2(1 - 2*power2(snt)))/(4.*power2
     (mh2)*power2(mH2)) + (3*At*ca*mu*sa*(3*power2(mh2)*(-1 + power2(sa)) + mh2
     *mH2*(-3 + 4*invdmhH*mH2)*(-1 + 2*power2(sa)) - (-1 + invdmhH*mH2)*power2(
     mH2)*(-1 + 2*power2(sa)))*power2(1 - 2*power2(snt)))/(2.*power2(mh2)*
     power2(mH2)) + DeltaInv(mst12,mH2,mst22)*((-9*At*ca*(mH2 - mst12 - mst22)*
     mu*sa*(-1 + power2(sa))*power2(1 - 2*power2(snt)))/(2.*mH2) + (9*(mH2 -
     mst12 - mst22)*power2(At)*(-1 + power2(sa))*power2(sa)*power2(1 - 2*power2
     (snt)))/(4.*mH2) - (9*(mH2 - mst12 - mst22)*mu2*power2(-1 + power2(sa))*
     power2(1 - 2*power2(snt)))/(4.*mH2)) - (3*power2(At)*(2*mh2*mH2*(-3 + 4*
     invdmhH*mH2) + 3*power2(mh2) - 2*(-1 + invdmhH*mH2)*power2(mH2))*(-1 +
     power2(sa))*power2(sa - 2*sa*power2(snt)))/(4.*power2(mh2)*power2(mH2))) +
     DeltaInv(mst12,mH2,mst22)*((-3*At*ca*mu*sa*(3*(logmst12Q2 - logmst22Q2)*(-
     6 + 2*logmH2Q2 + logmst12Q2 + logmst22Q2)*mst24 + power2(mH2)*(42 + 6*
     logmH2Q2*(-3 + logmst12Q2) - 18*logmst12Q2 + 3*power2(logmH2Q2) + 3*power2
     (logmst12Q2) + power2(Pi)) - mst12*mst22*(-6*(3 + logmH2Q2)*logmst22Q2 + 6
     *logmst12Q2*(-9 + logmH2Q2 + 2*logmst22Q2) + 9*power2(logmst12Q2) + 3*
     power2(logmst22Q2) + 2*(42 + power2(Pi))) - mH2*(mst12*(42 + 6*logmH2Q2*(-
     3 + logmst12Q2) - 18*logmst12Q2 + 3*power2(logmH2Q2) + 3*power2(logmst12Q2
     ) + power2(Pi)) + mst22*(42 - 36*logmst12Q2 + 6*logmH2Q2*(-3 + 2*
     logmst12Q2 - logmst22Q2) + 18*logmst22Q2 + 3*power2(logmH2Q2) + 6*power2(
     logmst12Q2) - 3*power2(logmst22Q2) + power2(Pi))))*(-1 + power2(sa))*
     power2(1 - 2*power2(snt)))/(4.*mH2) - (3*mu2*(3*(logmst12Q2 - logmst22Q2)*
     (-6 + 2*logmH2Q2 + logmst12Q2 + logmst22Q2)*mst24 + power2(mH2)*(42 + 6*
     logmH2Q2*(-3 + logmst12Q2) - 18*logmst12Q2 + 3*power2(logmH2Q2) + 3*power2
     (logmst12Q2) + power2(Pi)) - mst12*mst22*(-6*(3 + logmH2Q2)*logmst22Q2 + 6
     *logmst12Q2*(-9 + logmH2Q2 + 2*logmst22Q2) + 9*power2(logmst12Q2) + 3*
     power2(logmst22Q2) + 2*(42 + power2(Pi))) - mH2*(mst12*(42 + 6*logmH2Q2*(-
     3 + logmst12Q2) - 18*logmst12Q2 + 3*power2(logmH2Q2) + 3*power2(logmst12Q2
     ) + power2(Pi)) + mst22*(42 - 36*logmst12Q2 + 6*logmH2Q2*(-3 + 2*
     logmst12Q2 - logmst22Q2) + 18*logmst22Q2 + 3*power2(logmH2Q2) + 6*power2(
     logmst12Q2) - 3*power2(logmst22Q2) + power2(Pi))))*power2(-1 + power2(sa))
     *power2(1 - 2*power2(snt)))/(8.*mH2) + (3*power2(At)*(3*(logmst12Q2 -
     logmst22Q2)*(-6 + 2*logmH2Q2 + logmst12Q2 + logmst22Q2)*mst24 + power2(mH2
     )*(42 + 6*logmH2Q2*(-3 + logmst12Q2) - 18*logmst12Q2 + 3*power2(logmH2Q2)
     + 3*power2(logmst12Q2) + power2(Pi)) - mst12*mst22*(-6*(3 + logmH2Q2)*
     logmst22Q2 + 6*logmst12Q2*(-9 + logmH2Q2 + 2*logmst22Q2) + 9*power2(
     logmst12Q2) + 3*power2(logmst22Q2) + 2*(42 + power2(Pi))) - mH2*(mst12*(42
      + 6*logmH2Q2*(-3 + logmst12Q2) - 18*logmst12Q2 + 3*power2(logmH2Q2) + 3*
     power2(logmst12Q2) + power2(Pi)) + mst22*(42 - 36*logmst12Q2 + 6*logmH2Q2*
     (-3 + 2*logmst12Q2 - logmst22Q2) + 18*logmst22Q2 + 3*power2(logmH2Q2) + 6*
     power2(logmst12Q2) - 3*power2(logmst22Q2) + power2(Pi))))*(-1 + power2(sa)
     )*power2(sa - 2*sa*power2(snt)))/(8.*mH2)) + Fin3(msb12,mt2,mu2,Q2)*((3*
     mu2*(-msb14 + (mt2 - mu2)*mu2 + msb12*(mt2 + 2*mu2))*DeltaInv(msb12,mt2,
     mu2)*power2(invdmsntau2mu))/4. - (3*(msb12 - mt2 - mu2)*mu2*power3(
     invdmsntau2mu))/4.) + Fin3(mst12,mt2,mu2,Q2)*(-6*Al*At*invdmst*invdmstau*(
     invdmstau1mu - invdmstau2mu)*(1 + invdmstau1mu*mt2 + invdmstau2mu*mt2)*mu2
      + (6*At*invdmst*invdmstau*(invdmstau1mu - invdmstau2mu)*(1 + invdmstau1mu
     *mt2 + invdmstau2mu*mt2)*mu*mu2*sb)/cb + DeltaInv(mst12,mt2,mu2)*(6*Al*At*
     invdmst*invdmstau*(invdmstau1mu - invdmstau2mu)*mu2*(mst14 - 3*mst12*mt2 -
     2*mst12*mu2 - mt2*mu2 + power2(mu2)) - (6*At*invdmst*invdmstau*(
     invdmstau1mu - invdmstau2mu)*mu*mu2*sb*(mst14 - 3*mst12*mt2 - 2*mst12*mu2
     - mt2*mu2 + power2(mu2)))/cb + (3*mu2*((-mst14 + (mt2 - mu2)*mu2 + mst12*(
     mt2 + 2*mu2))*power2(invdmstau1mu) + 8*invdmst*invdmstau*invdmstau1mu*mu2*
     (mst14 - 3*mst12*mt2 - 2*mst12*mu2 - mt2*mu2 + power2(mu2)) + invdmstau2mu
     *(invdmstau2mu*(-mst14 + (mt2 - mu2)*mu2 + mst12*(mt2 + 2*mu2)) - 8*
     invdmst*invdmstau*mu2*(mst14 - 3*mst12*mt2 - 2*mst12*mu2 - mt2*mu2 +
     power2(mu2)))))/4.) + (3*mu2*(-8*invdmst*invdmstau*invdmstau1mu*mu2 - 8*
     invdmst*invdmstau*mt2*mu2*power2(invdmstau1mu) + invdmstau2mu*(8*invdmst*
     invdmstau*mu2 + 8*invdmst*invdmstau*invdmstau2mu*mt2*mu2 + (-mst12 + mt2 +
     mu2)*power2(invdmstau2mu)) + (-mst12 + mt2 + mu2)*power3(invdmstau1mu)))/
     4.) + Fin3(mst22,mt2,mu2,Q2)*(6*Al*At*invdmst*invdmstau*(invdmstau1mu -
     invdmstau2mu)*(1 + invdmstau1mu*mt2 + invdmstau2mu*mt2)*mu2 - (6*At*
     invdmst*invdmstau*(invdmstau1mu - invdmstau2mu)*(1 + invdmstau1mu*mt2 +
     invdmstau2mu*mt2)*mu*mu2*sb)/cb + DeltaInv(mst22,mt2,mu2)*(-6*Al*At*
     invdmst*invdmstau*(invdmstau1mu - invdmstau2mu)*mu2*(mst24 - 3*mst22*mt2 -
     2*mst22*mu2 - mt2*mu2 + power2(mu2)) + (6*At*invdmst*invdmstau*(
     invdmstau1mu - invdmstau2mu)*mu*mu2*sb*(mst24 - 3*mst22*mt2 - 2*mst22*mu2
     - mt2*mu2 + power2(mu2)))/cb - (3*mu2*((mst24 + mu2*(-mt2 + mu2) - mst22*(
     mt2 + 2*mu2))*power2(invdmstau1mu) + 8*invdmst*invdmstau*invdmstau1mu*mu2*
     (mst24 - 3*mst22*mt2 - 2*mst22*mu2 - mt2*mu2 + power2(mu2)) + invdmstau2mu
     *(invdmstau2mu*(mst24 + mu2*(-mt2 + mu2) - mst22*(mt2 + 2*mu2)) - 8*
     invdmst*invdmstau*mu2*(mst24 - 3*mst22*mt2 - 2*mst22*mu2 - mt2*mu2 +
     power2(mu2)))))/4.) + (3*mu2*(8*invdmst*invdmstau*invdmstau1mu*mu2 + 8*
     invdmst*invdmstau*mt2*mu2*power2(invdmstau1mu) + invdmstau2mu*(-8*invdmst*
     invdmstau*mu2 - 8*invdmst*invdmstau*invdmstau2mu*mt2*mu2 + (-mst22 + mt2 +
     mu2)*power2(invdmstau2mu)) + (-mst22 + mt2 + mu2)*power3(invdmstau1mu)))/
     4.) + (2*(1 - invdmhH*mH2 + 4*invdmhH*mt2)*DeltaInv(mh2,mt2,mt2)*(42 - 36*
     logmt2Q2 + 12*power2(logmt2Q2) + power2(Pi))*(-1 + power2(sa))*power2(sa)*
     power3(mt2))/mh2 + Fin3(mst12,mstau12,mt2,Q2)*(6*Al*At*invdmst*invdmstau*
     mt2*mu2*power2(invdmstau1mu) + (3*mu2*(1 + invdmstau1mu*(mst12 - mt2 - mu2
     ) + 8*invdmst*invdmstau*mt2*mu2)*power2(invdmstau1mu))/4. - (3*mu2*
     DeltaInv(mst12,mstau12,mt2)*(1 + invdmstau1mu*(2*mst12 + mt2 - 2*mu2) + (
     mst14 + mu2*(-mt2 + mu2) - mst12*(mt2 + 2*mu2))*power2(invdmstau1mu)))/4.
     - (6*At*invdmst*invdmstau*mt2*sb*power2(invdmstau1mu)*power3(mu))/cb) +
     Fin3(mst22,mstau12,mt2,Q2)*(-6*Al*At*invdmst*invdmstau*mt2*mu2*power2(
     invdmstau1mu) - (3*mu2*(-1 + 8*invdmst*invdmstau*mt2*mu2 + invdmstau1mu*(-
     mst22 + mt2 + mu2))*power2(invdmstau1mu))/4. - (3*mu2*DeltaInv(mst22,
     mstau12,mt2)*(1 + invdmstau1mu*(2*mst22 + mt2 - 2*mu2) + (mst24 + mu2*(-
     mt2 + mu2) - mst22*(mt2 + 2*mu2))*power2(invdmstau1mu)))/4. + (6*At*
     invdmst*invdmstau*mt2*sb*power2(invdmstau1mu)*power3(mu))/cb) + Fin3(
     mstau22,mst22,mt2,Q2)*(6*Al*At*invdmst*invdmstau*mt2*mu2*power2(
     invdmstau2mu) + (3*mu2*(1 + invdmstau2mu*(mst22 - mt2 - mu2) + 8*invdmst*
     invdmstau*mt2*mu2)*power2(invdmstau2mu))/4. - (3*mu2*DeltaInv(mstau22,
     mst22,mt2)*(1 + invdmstau2mu*(2*mst22 + mt2 - 2*mu2) + (mst24 + mu2*(-mt2
     + mu2) - mst22*(mt2 + 2*mu2))*power2(invdmstau2mu)))/4. - (6*At*invdmst*
     invdmstau*mt2*sb*power2(invdmstau2mu)*power3(mu))/cb) + Fin3(mstau22,mst12
     ,mt2,Q2)*(-6*Al*At*invdmst*invdmstau*mt2*mu2*power2(invdmstau2mu) - (3*mu2
     *(-1 + 8*invdmst*invdmstau*mt2*mu2 + invdmstau2mu*(-mst12 + mt2 + mu2))*
     power2(invdmstau2mu))/4. - (3*mu2*DeltaInv(mstau22,mst12,mt2)*(1 +
     invdmstau2mu*(2*mst12 + mt2 - 2*mu2) + (mst14 + mu2*(-mt2 + mu2) - mst12*(
     mt2 + 2*mu2))*power2(invdmstau2mu)))/4. + (6*At*invdmst*invdmstau*mt2*sb*
     power2(invdmstau2mu)*power3(mu))/cb) + DeltaInv(mst12,mst22,mh2)*((3*
     power2(At)*(3*(logmst12Q2 - logmst22Q2)*(-6 + 2*logmh2Q2 + logmst12Q2 +
     logmst22Q2)*mst24 + power2(mh2)*(42 + 6*logmh2Q2*(-3 + logmst12Q2) - 18*
     logmst12Q2 + 3*power2(logmh2Q2) + 3*power2(logmst12Q2) + power2(Pi)) -
     mst12*mst22*(-6*(3 + logmh2Q2)*logmst22Q2 + 6*logmst12Q2*(-9 + logmh2Q2 +
     2*logmst22Q2) + 9*power2(logmst12Q2) + 3*power2(logmst22Q2) + 2*(42 +
     power2(Pi))) - mh2*(mst12*(42 + 6*logmh2Q2*(-3 + logmst12Q2) - 18*
     logmst12Q2 + 3*power2(logmh2Q2) + 3*power2(logmst12Q2) + power2(Pi)) +
     mst22*(42 - 36*logmst12Q2 + 6*logmh2Q2*(-3 + 2*logmst12Q2 - logmst22Q2) +
     18*logmst22Q2 + 3*power2(logmh2Q2) + 6*power2(logmst12Q2) - 3*power2(
     logmst22Q2) + power2(Pi))))*(-1 + power2(sa))*power2(sa - 2*sa*power2(snt)
     ))/(8.*mh2) - (3*At*ca*mu*(3*(logmst12Q2 - logmst22Q2)*(-6 + 2*logmh2Q2 +
     logmst12Q2 + logmst22Q2)*mst24 + power2(mh2)*(42 + 6*logmh2Q2*(-3 +
     logmst12Q2) - 18*logmst12Q2 + 3*power2(logmh2Q2) + 3*power2(logmst12Q2) +
     power2(Pi)) - mst12*mst22*(-6*(3 + logmh2Q2)*logmst22Q2 + 6*logmst12Q2*(-9
      + logmh2Q2 + 2*logmst22Q2) + 9*power2(logmst12Q2) + 3*power2(logmst22Q2)
     + 2*(42 + power2(Pi))) - mh2*(mst12*(42 + 6*logmh2Q2*(-3 + logmst12Q2) -
     18*logmst12Q2 + 3*power2(logmh2Q2) + 3*power2(logmst12Q2) + power2(Pi)) +
     mst22*(42 - 36*logmst12Q2 + 6*logmh2Q2*(-3 + 2*logmst12Q2 - logmst22Q2) +
     18*logmst22Q2 + 3*power2(logmh2Q2) + 6*power2(logmst12Q2) - 3*power2(
     logmst22Q2) + power2(Pi))))*power2(1 - 2*power2(snt))*power3(sa))/(4.*mh2)
     - (3*mu2*(3*(logmst12Q2 - logmst22Q2)*(-6 + 2*logmh2Q2 + logmst12Q2 +
     logmst22Q2)*mst24 + power2(mh2)*(42 + 6*logmh2Q2*(-3 + logmst12Q2) - 18*
     logmst12Q2 + 3*power2(logmh2Q2) + 3*power2(logmst12Q2) + power2(Pi)) -
     mst12*mst22*(-6*(3 + logmh2Q2)*logmst22Q2 + 6*logmst12Q2*(-9 + logmh2Q2 +
     2*logmst22Q2) + 9*power2(logmst12Q2) + 3*power2(logmst22Q2) + 2*(42 +
     power2(Pi))) - mh2*(mst12*(42 + 6*logmh2Q2*(-3 + logmst12Q2) - 18*
     logmst12Q2 + 3*power2(logmh2Q2) + 3*power2(logmst12Q2) + power2(Pi)) +
     mst22*(42 - 36*logmst12Q2 + 6*logmh2Q2*(-3 + 2*logmst12Q2 - logmst22Q2) +
     18*logmst22Q2 + 3*power2(logmh2Q2) + 6*power2(logmst12Q2) - 3*power2(
     logmst22Q2) + power2(Pi))))*power2(1 - 2*power2(snt))*power4(sa))/(8.*mh2)
     ) + Fin3(mst12,mst22,mh2,Q2)*((9*(-1 + 2*invdmhH*mh2)*power2(At)*(-1 +
     power2(sa))*power2(sa)*power2(1 - 2*power2(snt)))/(4.*power2(mh2)) - (9*At
     *ca*mu*sa*(-power2(sa) + invdmhH*mh2*(-1 + 2*power2(sa)))*power2(1 - 2*
     power2(snt)))/(2.*power2(mh2)) - (9*mu2*(2*invdmhH*mh2*(-1 + power2(sa)) -
     power2(sa))*power2(sa - 2*sa*power2(snt)))/(4.*power2(mh2)) + DeltaInv(
     mst12,mst22,mh2)*((9*(mh2 - mst12 - mst22)*power2(At)*(-1 + power2(sa))*
     power2(sa)*power2(1 - 2*power2(snt)))/(4.*mh2) - (9*At*ca*(mh2 - mst12 -
     mst22)*mu*power2(1 - 2*power2(snt))*power3(sa))/(2.*mh2) - (9*(mh2 - mst12
      - mst22)*mu2*power2(1 - 2*power2(snt))*power4(sa))/(4.*mh2))) + Fin3(mH2,
     mst22,mst22,Q2)*((-3*(mH2 + mh2*(-3 + 4*invdmhH*mH2) - invdmhH*power2(mH2)
     )*(-1 + power2(sa))*power2(sa)*(mt2 + mu2*(-1 + power2(snt))*power2(snt)))
     /(mH2*power2(mh2)) + DeltaInv(mH2,mst22,mst22)*((9*mst22*(-1 + power2(sa))
     *(mt2*power2(sa) + mu2*(-1 + power2(sa))*(-1 + power2(snt))*power2(snt)))/
     mH2 - (18*At*ca*mst22*mu*sa*(-1 + power2(sa))*(invdmst*mt2 + power2(snt) -
     power4(snt)))/mH2 + (9*mst22*power2(At)*(-1 + power2(sa))*power2(sa)*(2*
     invdmst*mt2 + power2(snt) - power4(snt)))/mH2) + (3*At*ca*mu*sa*(mH2 + mh2
     *(-3 + 4*invdmhH*mH2) - invdmhH*power2(mH2))*(-1 + 2*power2(sa))*(invdmst*
     mt2 + power2(snt) - power4(snt)))/(mH2*power2(mh2)) + (3*(mh2*(3 - 4*
     invdmhH*mH2) + mH2*(-1 + invdmhH*mH2))*power2(At)*(-1 + power2(sa))*power2
     (sa)*(2*invdmst*mt2 + power2(snt) - power4(snt)))/(mH2*power2(mh2))) +
     Fin3(mst22,mst22,mh2,Q2)*((3*(-1 + invdmhH*(5*mh2 + mH2 - 4*mst22))*(-1 +
     power2(sa))*power2(sa)*(mt2 + mu2*(-1 + power2(snt))*power2(snt)))/(2.*
     power2(mh2)) + DeltaInv(mst22,mst22,mh2)*((3*power2(sa)*(-8*invdmhH*mst24*
     (-1 + power2(sa))*(mt2 + mu2*(-1 + power2(snt))*power2(snt)) + mst22*((1 +
     2*invdmhH*mH2)*mt2*(-1 + power2(sa)) + mu2*(2 + 2*invdmhH*mH2*(-1 + power2
     (sa)) + power2(sa))*(-1 + power2(snt))*power2(snt))))/mh2 - (6*At*ca*mu*sa
     *(4*invdmhH*mst24*(1 - 2*power2(sa)) + mst22*(1 + power2(sa) + invdmhH*mH2
     *(-1 + 2*power2(sa))))*(invdmst*mt2 + power2(snt) - power4(snt)))/mh2 + (3
     *(mst22 + 2*invdmhH*mH2*mst22 - 8*invdmhH*mst24)*power2(At)*(-1 + power2(
     sa))*power2(sa)*(2*invdmst*mt2 + power2(snt) - power4(snt)))/mh2) - (3*At*
     ca*(-1 + invdmhH*(5*mh2 + mH2 - 4*mst22))*mu*sa*(-1 + 2*power2(sa))*(
     invdmst*mt2 + power2(snt) - power4(snt)))/(2.*power2(mh2)) + (3*(-1 +
     invdmhH*(5*mh2 + mH2 - 4*mst22))*power2(At)*(-1 + power2(sa))*power2(sa)*(
     2*invdmst*mt2 + power2(snt) - power4(snt)))/(2.*power2(mh2))) + DeltaInv(
     mH2,mst22,mst22)*((9*mst24*(42 + 8*logmH2Q2*(-3 + logmst22Q2) - 12*
     logmst22Q2 + 4*power2(logmH2Q2) + power2(Pi))*(-1 + power2(sa))*(mt2*
     power2(sa) + mu2*(-1 + power2(sa))*(-1 + power2(snt))*power2(snt)))/(2.*
     mH2) - (9*At*ca*mst24*mu*sa*(42 + 8*logmH2Q2*(-3 + logmst22Q2) - 12*
     logmst22Q2 + 4*power2(logmH2Q2) + power2(Pi))*(-1 + power2(sa))*(invdmst*
     mt2 + power2(snt) - power4(snt)))/mH2 + (9*mst24*power2(At)*(42 + 8*
     logmH2Q2*(-3 + logmst22Q2) - 12*logmst22Q2 + 4*power2(logmH2Q2) + power2(
     Pi))*(-1 + power2(sa))*power2(sa)*(2*invdmst*mt2 + power2(snt) - power4(
     snt)))/(2.*mH2)) + DeltaInv(mst12,mst12,mH2)*((9*mst14*(42 + 8*logmH2Q2*(-
     3 + logmst12Q2) - 12*logmst12Q2 + 4*power2(logmH2Q2) + power2(Pi))*(-1 +
     power2(sa))*(mt2*power2(sa) + mu2*(-1 + power2(sa))*(-1 + power2(snt))*
     power2(snt)))/(2.*mH2) + (9*At*ca*mst14*mu*sa*(42 + 8*logmH2Q2*(-3 +
     logmst12Q2) - 12*logmst12Q2 + 4*power2(logmH2Q2) + power2(Pi))*(-1 +
     power2(sa))*(invdmst*mt2 - power2(snt) + power4(snt)))/mH2 - (9*mst14*
     power2(At)*(42 + 8*logmH2Q2*(-3 + logmst12Q2) - 12*logmst12Q2 + 4*power2(
     logmH2Q2) + power2(Pi))*(-1 + power2(sa))*power2(sa)*(2*invdmst*mt2 -
     power2(snt) + power4(snt)))/(2.*mH2)) + Fin3(mst12,mst12,mH2,Q2)*((-3*(mH2
      + mh2*(-3 + 4*invdmhH*mH2) - invdmhH*power2(mH2))*(-1 + power2(sa))*
     power2(sa)*(mt2 + mu2*(-1 + power2(snt))*power2(snt)))/(mH2*power2(mh2)) +
     (3*At*ca*(mh2*(3 - 4*invdmhH*mH2) + mH2*(-1 + invdmhH*mH2))*mu*sa*(-1 + 2*
     power2(sa))*(invdmst*mt2 - power2(snt) + power4(snt)))/(mH2*power2(mh2)) +
     (3*power2(At)*(mH2 + mh2*(-3 + 4*invdmhH*mH2) - invdmhH*power2(mH2))*(-1 +
     power2(sa))*power2(sa)*(2*invdmst*mt2 - power2(snt) + power4(snt)))/(mH2*
     power2(mh2)) + DeltaInv(mst12,mst12,mH2)*((9*mst12*(-1 + power2(sa))*(mt2*
     power2(sa) + mu2*(-1 + power2(sa))*(-1 + power2(snt))*power2(snt)))/mH2 +
     (18*At*ca*mst12*mu*sa*(-1 + power2(sa))*(invdmst*mt2 - power2(snt) +
     power4(snt)))/mH2 - (9*mst12*power2(At)*(-1 + power2(sa))*power2(sa)*(2*
     invdmst*mt2 - power2(snt) + power4(snt)))/mH2)) + Fin3(mst12,mst12,mh2,Q2)
     *((3*(-1 + invdmhH*(5*mh2 + mH2 - 4*mst12))*(-1 + power2(sa))*power2(sa)*(
     mt2 + mu2*(-1 + power2(snt))*power2(snt)))/(2.*power2(mh2)) + (3*At*ca*(-1
      + invdmhH*(5*mh2 + mH2 - 4*mst12))*mu*sa*(-1 + 2*power2(sa))*(invdmst*mt2
      - power2(snt) + power4(snt)))/(2.*power2(mh2)) - (3*(-1 + invdmhH*(5*mh2
     + mH2 - 4*mst12))*power2(At)*(-1 + power2(sa))*power2(sa)*(2*invdmst*mt2 -
     power2(snt) + power4(snt)))/(2.*power2(mh2)) + DeltaInv(mst12,mst12,mh2)*(
     (3*power2(sa)*(-8*invdmhH*mst14*(-1 + power2(sa))*(mt2 + mu2*(-1 + power2(
     snt))*power2(snt)) + mst12*((1 + 2*invdmhH*mH2)*mt2*(-1 + power2(sa)) +
     mu2*(2 + 2*invdmhH*mH2*(-1 + power2(sa)) + power2(sa))*(-1 + power2(snt))*
     power2(snt))))/mh2 + (6*At*ca*mu*sa*(4*invdmhH*mst14*(1 - 2*power2(sa)) +
     mst12*(1 + power2(sa) + invdmhH*mH2*(-1 + 2*power2(sa))))*(invdmst*mt2 -
     power2(snt) + power4(snt)))/mh2 - (3*(mst12 + 2*invdmhH*mH2*mst12 - 8*
     invdmhH*mst14)*power2(At)*(-1 + power2(sa))*power2(sa)*(2*invdmst*mt2 -
     power2(snt) + power4(snt)))/mh2)) + (power2(At)*(-1 + power2(sa))*power2(
     sa)*(-(mh2*mH2*(-8*invdmhH*power2(mH2)*(12*invdmst*logmst22Q2*mst22 - 12*
     logmst12Q2*(-1 + invdmst*mst22) - 12*logmH2Q2*(logmst12Q2 - invdmst*
     logmst12Q2*mst22 + invdmst*logmst22Q2*mst22) + 6*power2(logmH2Q2) + power2
     (Pi)) + 24*(2*invdmst*mt2*(mst12*(42 - 36*logmst12Q2 + 12*power2(
     logmst12Q2) + power2(Pi)) - mst22*(42 - 36*logmst22Q2 + 12*power2(
     logmst22Q2) + power2(Pi))) + (mst12*(42 - 36*logmst12Q2 + 12*power2(
     logmst12Q2) + power2(Pi)) + mst22*(42 - 36*logmst22Q2 + 12*power2(
     logmst22Q2) + power2(Pi)))*(-1 + power2(snt))*power2(snt)) + 2*power2(
     invdmhH)*(69 - 66*logmH2Q2 + 30*power2(logmH2Q2) + 5*power2(Pi))*power3(
     mH2) + mH2*(135 - 84*logmst12Q2 - 21*invdmhH*mst12 + 12*invdmhH*logmst12Q2
     *mst12 - 21*invdmhH*mst22 - 24*invdmst*logmst12Q2*mst22 + 12*invdmhH*
     logmst22Q2*mst22 + 24*invdmst*logmst22Q2*mst22 + 12*logmH2Q2*(7 - 6*
     logmst12Q2 + 6*invdmst*logmst12Q2*mst22 - 6*invdmst*logmst22Q2*mst22) -
     168*invdmhH*invdmst*mst12*mt2 + 192*invdmhH*invdmst*logmst12Q2*mst12*mt2 +
     168*invdmhH*invdmst*mst22*mt2 - 192*invdmhH*invdmst*logmst22Q2*mst22*mt2 +
     18*power2(logmst12Q2) - 96*invdmhH*invdmst*mst12*mt2*power2(logmst12Q2) +
     96*invdmhH*invdmst*mst22*mt2*power2(logmst22Q2) + 6*power2(Pi) - 8*invdmhH
     *invdmst*mst12*mt2*power2(Pi) + 8*invdmhH*invdmst*mst22*mt2*power2(Pi) -
     1008*power2(snt) + 432*logmst12Q2*power2(snt) + 168*invdmhH*mst12*power2(
     snt) - 144*invdmhH*logmst12Q2*mst12*power2(snt) + 168*invdmhH*mst22*power2
     (snt) - 144*invdmhH*logmst22Q2*mst22*power2(snt) - 72*power2(logmst12Q2)*
     power2(snt) + 48*invdmhH*mst12*power2(logmst12Q2)*power2(snt) + 48*invdmhH
     *mst22*power2(logmst22Q2)*power2(snt) - 24*power2(Pi)*power2(snt) + 4*
     invdmhH*mst12*power2(Pi)*power2(snt) + 4*invdmhH*mst22*power2(Pi)*power2(
     snt) + 18*power2(logmh2Q2)*power2(1 - 2*power2(snt)) + 1008*power4(snt) -
     432*logmst12Q2*power4(snt) - 168*invdmhH*mst12*power4(snt) + 144*invdmhH*
     logmst12Q2*mst12*power4(snt) - 168*invdmhH*mst22*power4(snt) + 144*invdmhH
     *logmst22Q2*mst22*power4(snt) + 72*power2(logmst12Q2)*power4(snt) - 48*
     invdmhH*mst12*power2(logmst12Q2)*power4(snt) - 48*invdmhH*mst22*power2(
     logmst22Q2)*power4(snt) + 24*power2(Pi)*power4(snt) - 4*invdmhH*mst12*
     power2(Pi)*power4(snt) - 4*invdmhH*mst22*power2(Pi)*power4(snt) + 36*
     logmh2Q2*(logmst12Q2*power2(1 - 2*power2(snt)) - 2*(1 - 6*power2(snt) + 6*
     power4(snt)))))) + power2(mh2)*(3*invdmhH*(27 + 20*logmH2Q2 - 20*
     logmst12Q2 + 20*invdmst*logmst12Q2*mst22 - 20*invdmst*logmst22Q2*mst22 +
     24*logmh2Q2*(-2 + logmst12Q2 - invdmst*logmst12Q2*mst22 + invdmst*
     logmst22Q2*mst22) + 12*power2(logmh2Q2) - 12*power2(logmH2Q2))*power2(mH2)
     - 12*(2*invdmst*mt2*(mst12*(42 + 12*logmH2Q2*(-3 + logmst12Q2) + 6*power2(
     logmH2Q2) - 6*power2(logmst12Q2) + power2(Pi)) - mst22*(42 + 12*logmH2Q2*(
     -3 + logmst22Q2) + 6*power2(logmH2Q2) - 6*power2(logmst22Q2) + power2(Pi))
     ) + (mst12*(42 + 12*logmH2Q2*(-3 + logmst12Q2) + 6*power2(logmH2Q2) - 6*
     power2(logmst12Q2) + power2(Pi)) + mst22*(42 + 12*logmH2Q2*(-3 +
     logmst22Q2) + 6*power2(logmH2Q2) - 6*power2(logmst22Q2) + power2(Pi)))*(-1
      + power2(snt))*power2(snt)) + 8*power2(invdmhH)*(12 - 12*logmH2Q2 + 6*
     power2(logmH2Q2) + power2(Pi))*power3(mH2) - 6*mH2*(36 - 6*invdmst*
     logmst22Q2*mst22 + power2(Pi) - 168*power2(snt) - 4*power2(Pi)*power2(snt)
     + 3*power2(logmH2Q2)*power2(1 - 2*power2(snt)) + 3*power2(logmst12Q2)*
     power2(1 - 2*power2(snt)) + 6*logmst12Q2*(-4 + invdmst*mst22 + 12*power2(
     snt) - 12*power4(snt)) + 168*power4(snt) + 4*power2(Pi)*power4(snt) + 6*
     logmH2Q2*(logmst12Q2*power2(1 - 2*power2(snt)) - 2*(1 - 6*power2(snt) + 6*
     power4(snt))))) + power2(mH2)*(21*mst22 - 12*logmst22Q2*mst22 + 1344*
     invdmhH*invdmst*mst14*mt2 - 1152*invdmhH*invdmst*logmst12Q2*mst14*mt2 +
     840*invdmst*mst22*mt2 - 864*invdmst*logmh2Q2*mst22*mt2 + 192*invdmst*
     logmst22Q2*mst22*mt2 + 288*invdmst*logmh2Q2*logmst22Q2*mst22*mt2 - 1344*
     invdmhH*invdmst*mst24*mt2 + 1152*invdmhH*invdmst*logmst22Q2*mst24*mt2 +
     144*invdmst*mst22*mt2*power2(logmh2Q2) + 384*invdmhH*invdmst*mst14*mt2*
     power2(logmst12Q2) - 240*invdmst*mst22*mt2*power2(logmst22Q2) - 384*
     invdmhH*invdmst*mst24*mt2*power2(logmst22Q2) + 32*invdmhH*invdmst*mst14*
     mt2*power2(Pi) + 16*invdmst*mst22*mt2*power2(Pi) - 32*invdmhH*invdmst*
     mst24*mt2*power2(Pi) - invdmhH*power2(mH2)*(27 + 36*invdmst*logmst22Q2*
     mst22 - 36*logmst12Q2*(-1 + invdmst*mst22) - 24*logmH2Q2*(1 + logmst12Q2 -
     invdmst*logmst12Q2*mst22 + invdmst*logmst22Q2*mst22) + 12*power2(logmH2Q2)
     + 2*power2(Pi)) - 672*invdmhH*mst14*power2(snt) + 576*invdmhH*logmst12Q2*
     mst14*power2(snt) + 336*mst22*power2(snt) - 432*logmh2Q2*mst22*power2(snt)
     + 144*logmst22Q2*mst22*power2(snt) + 144*logmh2Q2*logmst22Q2*mst22*power2(
     snt) - 672*invdmhH*mst24*power2(snt) + 576*invdmhH*logmst22Q2*mst24*power2
     (snt) + 72*mst22*power2(logmh2Q2)*power2(snt) - 192*invdmhH*mst14*power2(
     logmst12Q2)*power2(snt) - 120*mst22*power2(logmst22Q2)*power2(snt) - 192*
     invdmhH*mst24*power2(logmst22Q2)*power2(snt) - 16*invdmhH*mst14*power2(Pi)
     *power2(snt) + 8*mst22*power2(Pi)*power2(snt) - 16*invdmhH*mst24*power2(Pi
     )*power2(snt) + 2*power2(invdmhH)*(21 - 18*logmH2Q2 + 6*power2(logmH2Q2) +
     power2(Pi))*power3(mH2) + 672*invdmhH*mst14*power4(snt) - 576*invdmhH*
     logmst12Q2*mst14*power4(snt) - 336*mst22*power4(snt) + 432*logmh2Q2*mst22*
     power4(snt) - 144*logmst22Q2*mst22*power4(snt) - 144*logmh2Q2*logmst22Q2*
     mst22*power4(snt) + 672*invdmhH*mst24*power4(snt) - 576*invdmhH*logmst22Q2
     *mst24*power4(snt) - 72*mst22*power2(logmh2Q2)*power4(snt) + 192*invdmhH*
     mst14*power2(logmst12Q2)*power4(snt) + 120*mst22*power2(logmst22Q2)*power4
     (snt) + 192*invdmhH*mst24*power2(logmst22Q2)*power4(snt) + 16*invdmhH*
     mst14*power2(Pi)*power4(snt) - 8*mst22*power2(Pi)*power4(snt) + 16*invdmhH
     *mst24*power2(Pi)*power4(snt) + mst12*(21 - 840*invdmst*mt2 + 864*invdmst*
     logmh2Q2*mt2 - 144*invdmst*mt2*power2(logmh2Q2) - 16*invdmst*mt2*power2(Pi
     ) + 336*power2(snt) - 432*logmh2Q2*power2(snt) + 72*power2(logmh2Q2)*
     power2(snt) + 8*power2(Pi)*power2(snt) - 336*power4(snt) + 432*logmh2Q2*
     power4(snt) - 72*power2(logmh2Q2)*power4(snt) - 8*power2(Pi)*power4(snt) +
     120*power2(logmst12Q2)*(2*invdmst*mt2 - power2(snt) + power4(snt)) - 12*
     logmst12Q2*(1 + 8*invdmst*(2 + 3*logmh2Q2)*mt2 - 12*(1 + logmh2Q2)*power2(
     snt) + 12*(1 + logmh2Q2)*power4(snt))) + mH2*(-15 - 21*invdmhH*mst12 - 21*
     invdmhH*mst22 + 12*invdmhH*logmst22Q2*mst22 + 36*invdmst*logmst22Q2*mst22
     + 12*logmH2Q2*(1 - 2*invdmst*logmst22Q2*mst22 + 2*logmst12Q2*(-1 + invdmst
     *mst22)) - 168*invdmhH*invdmst*mst12*mt2 + 168*invdmhH*invdmst*mst22*mt2 -
     192*invdmhH*invdmst*logmst22Q2*mst22*mt2 + 96*invdmhH*invdmst*mst22*mt2*
     power2(logmst22Q2) - 8*invdmhH*invdmst*mst12*mt2*power2(Pi) + 8*invdmhH*
     invdmst*mst22*mt2*power2(Pi) + 168*invdmhH*mst12*power2(snt) + 168*invdmhH
     *mst22*power2(snt) - 144*invdmhH*logmst22Q2*mst22*power2(snt) + 48*invdmhH
     *mst22*power2(logmst22Q2)*power2(snt) + 4*invdmhH*mst12*power2(Pi)*power2(
     snt) + 4*invdmhH*mst22*power2(Pi)*power2(snt) - 168*invdmhH*mst12*power4(
     snt) - 168*invdmhH*mst22*power4(snt) + 144*invdmhH*logmst22Q2*mst22*power4
     (snt) - 48*invdmhH*mst22*power2(logmst22Q2)*power4(snt) - 4*invdmhH*mst12*
     power2(Pi)*power4(snt) - 4*invdmhH*mst22*power2(Pi)*power4(snt) - 48*
     invdmhH*mst12*power2(logmst12Q2)*(2*invdmst*mt2 - power2(snt) + power4(snt
     )) + 12*logmst12Q2*(3 - 3*invdmst*mst22 + invdmhH*mst12*(1 + 16*invdmst*
     mt2 - 12*power2(snt) + 12*power4(snt)))))))/(16.*power2(mh2)*power2(mH2))
     + (At*mu*(-4*mH2*sb*(-3*mA2*mH2*(-1 + logmst12Q2 - invdmst*logmst12Q2*
     mst22 + invdmst*logmst22Q2*mst22)*power2(mh2) + 3*mC2*mH2*(-1 + logmst12Q2
      - invdmst*logmst12Q2*mst22 + invdmst*logmst22Q2*mst22)*power2(mh2) + mA2*
     mC2*(-3*(-3 + 2*logmH2Q2)*(-1 + invdmhH*mH2)*(-1 + logmst12Q2 - invdmst*
     logmst12Q2*mst22 + invdmst*logmst22Q2*mst22)*power2(mH2)*(-1 + power2(sa))
     *power2(sa) - 3*mh2*mH2*(1 - invdmst*logmst22Q2*mst22 + logmst12Q2*(-1 +
     invdmst*mst22))*power2(sa)*(-5 - 8*invdmhH*mH2*(-1 + power2(sa)) + 2*
     logmH2Q2*(-3 + 4*invdmhH*mH2)*(-1 + power2(sa)) + 2*power2(sa)) + power2(
     mh2)*(-12*invdmst*invdmstau*(invdmstau1mu - invdmstau2mu)*mH2*mst12*mu2*
     power2(logmst12Q2) + 3*(-1 + power2(sa))*(-3 + (3 + invdmhH*(-5 + 6*
     logmh2Q2)*mH2)*power2(sa)) + invdmst*(12*invdmstau*(invdmstau1mu -
     invdmstau2mu)*mH2*mst22*mu2*power2(logmst22Q2) - 4*invdmstau*(invdmstau1mu
      - invdmstau2mu)*mH2*(mst12 - mst22)*mu2*(6*logmt2Q2*(-9 + logmu2Q2) + 9*
     power2(logmt2Q2) + 2*(42 + power2(Pi))) + 3*logmst22Q2*mst22*(8*invdmstau*
     (invdmstau1mu - invdmstau2mu)*(-3 + 2*logmt2Q2 - logmu2Q2)*mH2*mu2 - (-1 +
     power2(sa))*(-3 + (3 + invdmhH*(-5 + 6*logmh2Q2)*mH2)*power2(sa)))) + 3*
     logmst12Q2*(-((-1 + power2(sa))*(-3 + (3 + invdmhH*(-5 + 6*logmh2Q2)*mH2)*
     power2(sa))) + invdmst*(-8*invdmstau*(invdmstau1mu - invdmstau2mu)*(-3 + 2
     *logmt2Q2 - logmu2Q2)*mH2*mst12*mu2 + mst22*(-1 + power2(sa))*(-3 + (3 +
     invdmhH*(-5 + 6*logmh2Q2)*mH2)*power2(sa))))))) + ca*cb*mA2*mC2*sa*(power2
     (mh2)*(-3*invdmhH*(7 - 24*logmh2Q2 + 20*logmH2Q2 + 12*power2(logmh2Q2) -
     12*power2(logmH2Q2))*power2(mH2)*(-1 + 2*power2(sa)) + 24*(-1 + power2(sa)
     )*(invdmst*mt2*(mst12*(42 + 12*logmH2Q2*(-3 + logmst12Q2) + 6*power2(
     logmH2Q2) - 6*power2(logmst12Q2) + power2(Pi)) - mst22*(42 + 12*logmH2Q2*(
     -3 + logmst22Q2) + 6*power2(logmH2Q2) - 6*power2(logmst22Q2) + power2(Pi))
     ) + (mst12*(42 + 12*logmH2Q2*(-3 + logmst12Q2) + 6*power2(logmH2Q2) - 6*
     power2(logmst12Q2) + power2(Pi)) + mst22*(42 + 12*logmH2Q2*(-3 +
     logmst22Q2) + 6*power2(logmH2Q2) - 6*power2(logmst22Q2) + power2(Pi)))*(-1
      + power2(snt))*power2(snt)) - 8*power2(invdmhH)*(12 - 12*logmH2Q2 + 6*
     power2(logmH2Q2) + power2(Pi))*(-1 + 2*power2(sa))*power3(mH2) + 12*mH2*(-
     1 + power2(sa))*(30 + power2(Pi) - 168*power2(snt) - 4*power2(Pi)*power2(
     snt) - 18*logmst12Q2*power2(1 - 2*power2(snt)) + 3*power2(logmH2Q2)*power2
     (1 - 2*power2(snt)) + 3*power2(logmst12Q2)*power2(1 - 2*power2(snt)) + 168
     *power4(snt) + 4*power2(Pi)*power4(snt) + 6*logmH2Q2*(logmst12Q2*power2(1
     - 2*power2(snt)) - 2*(1 - 6*power2(snt) + 6*power4(snt))))) + power2(mH2)*
     (21*mst22 - 12*logmst22Q2*mst22 + 672*invdmhH*invdmst*mst14*mt2 - 576*
     invdmhH*invdmst*logmst12Q2*mst14*mt2 - 84*invdmst*mst22*mt2 + 96*invdmst*
     logmst22Q2*mst22*mt2 - 672*invdmhH*invdmst*mst24*mt2 + 576*invdmhH*invdmst
     *logmst22Q2*mst24*mt2 + 192*invdmhH*invdmst*mst14*mt2*power2(logmst12Q2) -
     48*invdmst*mst22*mt2*power2(logmst22Q2) - 192*invdmhH*invdmst*mst24*mt2*
     power2(logmst22Q2) + 16*invdmhH*invdmst*mst14*mt2*power2(Pi) - 4*invdmst*
     mst22*mt2*power2(Pi) - 16*invdmhH*invdmst*mst24*mt2*power2(Pi) - 42*mst22*
     power2(sa) + 24*logmst22Q2*mst22*power2(sa) - 1344*invdmhH*invdmst*mst14*
     mt2*power2(sa) + 1152*invdmhH*invdmst*logmst12Q2*mst14*mt2*power2(sa) -
     840*invdmst*mst22*mt2*power2(sa) + 864*invdmst*logmh2Q2*mst22*mt2*power2(
     sa) - 192*invdmst*logmst22Q2*mst22*mt2*power2(sa) - 288*invdmst*logmh2Q2*
     logmst22Q2*mst22*mt2*power2(sa) + 1344*invdmhH*invdmst*mst24*mt2*power2(sa
     ) - 1152*invdmhH*invdmst*logmst22Q2*mst24*mt2*power2(sa) - 144*invdmst*
     mst22*mt2*power2(logmh2Q2)*power2(sa) - 384*invdmhH*invdmst*mst14*mt2*
     power2(logmst12Q2)*power2(sa) + 240*invdmst*mst22*mt2*power2(logmst22Q2)*
     power2(sa) + 384*invdmhH*invdmst*mst24*mt2*power2(logmst22Q2)*power2(sa) -
     32*invdmhH*invdmst*mst14*mt2*power2(Pi)*power2(sa) - 16*invdmst*mst22*mt2*
     power2(Pi)*power2(sa) + 32*invdmhH*invdmst*mst24*mt2*power2(Pi)*power2(sa)
     + invdmhH*power2(mH2)*(63 - 48*logmH2Q2 + 12*power2(logmH2Q2) + 2*power2(
     Pi))*(-1 + 2*power2(sa)) - 672*invdmhH*mst14*power2(snt) + 576*invdmhH*
     logmst12Q2*mst14*power2(snt) - 168*mst22*power2(snt) + 144*logmst22Q2*
     mst22*power2(snt) - 672*invdmhH*mst24*power2(snt) + 576*invdmhH*logmst22Q2
     *mst24*power2(snt) - 192*invdmhH*mst14*power2(logmst12Q2)*power2(snt) - 48
     *mst22*power2(logmst22Q2)*power2(snt) - 192*invdmhH*mst24*power2(
     logmst22Q2)*power2(snt) - 16*invdmhH*mst14*power2(Pi)*power2(snt) - 4*
     mst22*power2(Pi)*power2(snt) - 16*invdmhH*mst24*power2(Pi)*power2(snt) +
     1344*invdmhH*mst14*power2(sa)*power2(snt) - 1152*invdmhH*logmst12Q2*mst14*
     power2(sa)*power2(snt) - 672*mst22*power2(sa)*power2(snt) + 864*logmh2Q2*
     mst22*power2(sa)*power2(snt) - 288*logmst22Q2*mst22*power2(sa)*power2(snt)
     - 288*logmh2Q2*logmst22Q2*mst22*power2(sa)*power2(snt) + 1344*invdmhH*
     mst24*power2(sa)*power2(snt) - 1152*invdmhH*logmst22Q2*mst24*power2(sa)*
     power2(snt) - 144*mst22*power2(logmh2Q2)*power2(sa)*power2(snt) + 384*
     invdmhH*mst14*power2(logmst12Q2)*power2(sa)*power2(snt) + 240*mst22*power2
     (logmst22Q2)*power2(sa)*power2(snt) + 384*invdmhH*mst24*power2(logmst22Q2)
     *power2(sa)*power2(snt) + 32*invdmhH*mst14*power2(Pi)*power2(sa)*power2(
     snt) - 16*mst22*power2(Pi)*power2(sa)*power2(snt) + 32*invdmhH*mst24*
     power2(Pi)*power2(sa)*power2(snt) - 2*power2(invdmhH)*(21 - 18*logmH2Q2 +
     6*power2(logmH2Q2) + power2(Pi))*(-1 + 2*power2(sa))*power3(mH2) + 672*
     invdmhH*mst14*power4(snt) - 576*invdmhH*logmst12Q2*mst14*power4(snt) + 168
     *mst22*power4(snt) - 144*logmst22Q2*mst22*power4(snt) + 672*invdmhH*mst24*
     power4(snt) - 576*invdmhH*logmst22Q2*mst24*power4(snt) + 192*invdmhH*mst14
     *power2(logmst12Q2)*power4(snt) + 48*mst22*power2(logmst22Q2)*power4(snt)
     + 192*invdmhH*mst24*power2(logmst22Q2)*power4(snt) + 16*invdmhH*mst14*
     power2(Pi)*power4(snt) + 4*mst22*power2(Pi)*power4(snt) + 16*invdmhH*mst24
     *power2(Pi)*power4(snt) - 1344*invdmhH*mst14*power2(sa)*power4(snt) + 1152
     *invdmhH*logmst12Q2*mst14*power2(sa)*power4(snt) + 672*mst22*power2(sa)*
     power4(snt) - 864*logmh2Q2*mst22*power2(sa)*power4(snt) + 288*logmst22Q2*
     mst22*power2(sa)*power4(snt) + 288*logmh2Q2*logmst22Q2*mst22*power2(sa)*
     power4(snt) - 1344*invdmhH*mst24*power2(sa)*power4(snt) + 1152*invdmhH*
     logmst22Q2*mst24*power2(sa)*power4(snt) + 144*mst22*power2(logmh2Q2)*
     power2(sa)*power4(snt) - 384*invdmhH*mst14*power2(logmst12Q2)*power2(sa)*
     power4(snt) - 240*mst22*power2(logmst22Q2)*power2(sa)*power4(snt) - 384*
     invdmhH*mst24*power2(logmst22Q2)*power2(sa)*power4(snt) - 32*invdmhH*mst14
     *power2(Pi)*power2(sa)*power4(snt) + 16*mst22*power2(Pi)*power2(sa)*power4
     (snt) - 32*invdmhH*mst24*power2(Pi)*power2(sa)*power4(snt) + mH2*(-1 + 2*
     power2(sa))*(-21 + 12*logmH2Q2 - invdmhH*mst22*(-21 + 4*invdmst*mt2*(21 +
     power2(Pi)) + 168*power2(snt) + 4*power2(Pi)*power2(snt) - 12*logmst22Q2*(
     -1 + 8*invdmst*mt2 + 12*power2(snt) - 12*power4(snt)) + 48*power2(
     logmst22Q2)*(invdmst*mt2 + power2(snt) - power4(snt)) - 168*power4(snt) -
     4*power2(Pi)*power4(snt)) + invdmhH*mst12*(21 + 4*invdmst*mt2*(21 + power2
     (Pi)) - 168*power2(snt) - 4*power2(Pi)*power2(snt) + 168*power4(snt) + 4*
     power2(Pi)*power4(snt) + 48*power2(logmst12Q2)*(invdmst*mt2 - power2(snt)
     + power4(snt)) - 12*logmst12Q2*(1 + 8*invdmst*mt2 - 12*power2(snt) + 12*
     power4(snt)))) + mst12*(21 + 84*invdmst*mt2 + 4*invdmst*mt2*power2(Pi) -
     42*power2(sa) + 840*invdmst*mt2*power2(sa) - 864*invdmst*logmh2Q2*mt2*
     power2(sa) + 144*invdmst*mt2*power2(logmh2Q2)*power2(sa) + 16*invdmst*mt2*
     power2(Pi)*power2(sa) - 168*power2(snt) - 4*power2(Pi)*power2(snt) - 672*
     power2(sa)*power2(snt) + 864*logmh2Q2*power2(sa)*power2(snt) - 144*power2(
     logmh2Q2)*power2(sa)*power2(snt) - 16*power2(Pi)*power2(sa)*power2(snt) +
     168*power4(snt) + 4*power2(Pi)*power4(snt) + 672*power2(sa)*power4(snt) -
     864*logmh2Q2*power2(sa)*power4(snt) + 144*power2(logmh2Q2)*power2(sa)*
     power4(snt) + 16*power2(Pi)*power2(sa)*power4(snt) - 48*power2(logmst12Q2)
     *(-1 + 5*power2(sa))*(invdmst*mt2 - power2(snt) + power4(snt)) + 12*
     logmst12Q2*(-1 + 8*invdmst*mt2*(-1 + (2 + 3*logmh2Q2)*power2(sa)) + 12*
     power2(snt) - 12*power4(snt) + power2(sa)*(2 - 24*(1 + logmh2Q2)*power2(
     snt) + 24*(1 + logmh2Q2)*power4(snt))))) + mh2*mH2*(-8*invdmhH*power2(mH2)
     *(12 - 12*logmH2Q2 + 6*power2(logmH2Q2) + power2(Pi))*(-1 + 2*power2(sa))
     + 24*(-1 + 2*power2(sa))*(invdmst*mt2*(mst12*(42 - 36*logmst12Q2 + 12*
     power2(logmst12Q2) + power2(Pi)) - mst22*(42 - 36*logmst22Q2 + 12*power2(
     logmst22Q2) + power2(Pi))) + (mst12*(42 - 36*logmst12Q2 + 12*power2(
     logmst12Q2) + power2(Pi)) + mst22*(42 - 36*logmst22Q2 + 12*power2(
     logmst22Q2) + power2(Pi)))*(-1 + power2(snt))*power2(snt)) + 2*power2(
     invdmhH)*(69 - 66*logmH2Q2 + 30*power2(logmH2Q2) + 5*power2(Pi))*(-1 + 2*
     power2(sa))*power3(mH2) + mH2*(21 + 12*logmH2Q2*(-1 + 2*power2(sa)) + 6*
     power2(sa)*(53 + 2*power2(Pi) - 336*power2(snt) - 8*power2(Pi)*power2(snt)
     - 36*logmst12Q2*power2(1 - 2*power2(snt)) + 6*power2(logmh2Q2)*power2(1 -
     2*power2(snt)) + 6*power2(logmst12Q2)*power2(1 - 2*power2(snt)) + 336*
     power4(snt) + 8*power2(Pi)*power4(snt) + 12*logmh2Q2*(logmst12Q2*power2(1
     - 2*power2(snt)) - 2*(1 - 6*power2(snt) + 6*power4(snt)))) - invdmhH*(-1 +
     2*power2(sa))*(-(mst22*(-21 + 4*invdmst*mt2*(21 + power2(Pi)) + 168*power2
     (snt) + 4*power2(Pi)*power2(snt) - 12*logmst22Q2*(-1 + 8*invdmst*mt2 + 12*
     power2(snt) - 12*power4(snt)) + 48*power2(logmst22Q2)*(invdmst*mt2 +
     power2(snt) - power4(snt)) - 168*power4(snt) - 4*power2(Pi)*power4(snt)))
     + mst12*(21 + 4*invdmst*mt2*(21 + power2(Pi)) - 168*power2(snt) - 4*power2
     (Pi)*power2(snt) + 168*power4(snt) + 4*power2(Pi)*power4(snt) + 48*power2(
     logmst12Q2)*(invdmst*mt2 - power2(snt) + power4(snt)) - 12*logmst12Q2*(1 +
     8*invdmst*mt2 - 12*power2(snt) + 12*power4(snt)))))))))/(16.*cb*mA2*mC2*
     power2(mh2)*power2(mH2)) + (mu2*DeltaInv(msb12,mt2,mu2)*power2(
     invdmsntau2mu)*(msb14*mt2*(42 + 6*logmsb12Q2*(-3 + logmt2Q2) - 18*logmt2Q2
      + 3*power2(logmsb12Q2) + 3*power2(logmt2Q2) + power2(Pi)) + msb14*mu2*(42
      + 12*logmsb12Q2*(-3 + logmt2Q2) + 18*logmu2Q2 - 6*logmt2Q2*(3 + logmu2Q2)
     + 6*power2(logmsb12Q2) + 3*power2(logmt2Q2) - 3*power2(logmu2Q2) + power2(
     Pi)) + (mt2 - mu2)*power2(mu2)*(42 + 6*logmt2Q2*(-3 + logmu2Q2) - 18*
     logmu2Q2 + 3*power2(logmt2Q2) + 3*power2(logmu2Q2) + power2(Pi)) + msb12*
     mu2*(mu2*(42 - 6*logmsb12Q2*(-3 + logmt2Q2) - 36*logmu2Q2 + 6*logmt2Q2*(-3
      + 2*logmu2Q2) - 3*power2(logmsb12Q2) + 3*power2(logmt2Q2) + 6*power2(
     logmu2Q2) + power2(Pi)) + mt2*(168 + 6*logmsb12Q2*(-3 + 3*logmt2Q2 - 2*
     logmu2Q2) + 18*logmt2Q2*(-6 + logmu2Q2) - 18*logmu2Q2 + 3*power2(
     logmsb12Q2) + 18*power2(logmt2Q2) + 3*power2(logmu2Q2) + 4*power2(Pi))) -
     (42 + 6*logmsb12Q2*(-3 + logmt2Q2) - 18*logmt2Q2 + 3*power2(logmsb12Q2) +
     3*power2(logmt2Q2) + power2(Pi))*power6(msb1)))/8. + (mu2*DeltaInv(msntau2
     ,msb12,mt2)*(-((msntau2 - mt2 + 2*mu2)*(42 + 6*logmsntau2Q2*(-3 + logmt2Q2
     ) - 18*logmt2Q2 + 3*power2(logmsntau2Q2) + 3*power2(logmt2Q2) + power2(Pi)
     )) + msb12*(42 - 6*logmsb12Q2*(-3 + logmt2Q2) + 12*logmsntau2Q2*(-3 +
     logmt2Q2) - 18*logmt2Q2 - 3*power2(logmsb12Q2) + 6*power2(logmsntau2Q2) +
     3*power2(logmt2Q2) + power2(Pi)) - invdmsntau2mu*(msb14*(42 + 12*
     logmsb12Q2*(-3 + logmt2Q2) - 6*logmsntau2Q2*(-3 + logmt2Q2) - 18*logmt2Q2
     + 6*power2(logmsb12Q2) - 3*power2(logmsntau2Q2) + 3*power2(logmt2Q2) +
     power2(Pi)) + (2*mt2 - 3*mu2)*mu2*(42 + 6*logmsntau2Q2*(-3 + logmt2Q2) -
     18*logmt2Q2 + 3*power2(logmsntau2Q2) + 3*power2(logmt2Q2) + power2(Pi)) +
     msb12*(2*mu2*(42 - 6*logmsb12Q2*(-3 + logmt2Q2) + 12*logmsntau2Q2*(-3 +
     logmt2Q2) - 18*logmt2Q2 - 3*power2(logmsb12Q2) + 6*power2(logmsntau2Q2) +
     3*power2(logmt2Q2) + power2(Pi)) + mt2*(-6*logmsb12Q2*(3 + 2*logmsntau2Q2
     - 3*logmt2Q2) + 18*logmsntau2Q2*(-1 + logmt2Q2) + 3*power2(logmsb12Q2) + 3
     *power2(logmsntau2Q2) + 2*(84 - 54*logmt2Q2 + 9*power2(logmt2Q2) + 2*
     power2(Pi))))) + power2(invdmsntau2mu)*(msb14*mt2*(42 + 6*logmsb12Q2*(-3 +
     logmt2Q2) - 18*logmt2Q2 + 3*power2(logmsb12Q2) + 3*power2(logmt2Q2) +
     power2(Pi)) + msb14*mu2*(42 + 12*logmsb12Q2*(-3 + logmt2Q2) - 6*
     logmsntau2Q2*(-3 + logmt2Q2) - 18*logmt2Q2 + 6*power2(logmsb12Q2) - 3*
     power2(logmsntau2Q2) + 3*power2(logmt2Q2) + power2(Pi)) + (mt2 - mu2)*
     power2(mu2)*(42 + 6*logmsntau2Q2*(-3 + logmt2Q2) - 18*logmt2Q2 + 3*power2(
     logmsntau2Q2) + 3*power2(logmt2Q2) + power2(Pi)) + msb12*mu2*(mu2*(42 - 6*
     logmsb12Q2*(-3 + logmt2Q2) + 12*logmsntau2Q2*(-3 + logmt2Q2) - 18*logmt2Q2
      - 3*power2(logmsb12Q2) + 6*power2(logmsntau2Q2) + 3*power2(logmt2Q2) +
     power2(Pi)) + mt2*(-6*logmsb12Q2*(3 + 2*logmsntau2Q2 - 3*logmt2Q2) + 18*
     logmsntau2Q2*(-1 + logmt2Q2) + 3*power2(logmsb12Q2) + 3*power2(
     logmsntau2Q2) + 2*(84 - 54*logmt2Q2 + 9*power2(logmt2Q2) + 2*power2(Pi))))
     - (42 + 6*logmsb12Q2*(-3 + logmt2Q2) - 18*logmt2Q2 + 3*power2(logmsb12Q2)
     + 3*power2(logmt2Q2) + power2(Pi))*power6(msb1))))/8. + (mu2*DeltaInv(
     mst12,mstau12,mt2)*(-((mstau12 - mt2 + 2*mu2)*(42 + 6*logmstau12Q2*(-3 +
     logmt2Q2) - 18*logmt2Q2 + 3*power2(logmstau12Q2) + 3*power2(logmt2Q2) +
     power2(Pi))) + mst12*(42 - 6*logmst12Q2*(-3 + logmt2Q2) + 12*logmstau12Q2*
     (-3 + logmt2Q2) - 18*logmt2Q2 - 3*power2(logmst12Q2) + 6*power2(
     logmstau12Q2) + 3*power2(logmt2Q2) + power2(Pi)) - invdmstau1mu*(mst14*(42
      + 12*logmst12Q2*(-3 + logmt2Q2) - 6*logmstau12Q2*(-3 + logmt2Q2) - 18*
     logmt2Q2 + 6*power2(logmst12Q2) - 3*power2(logmstau12Q2) + 3*power2(
     logmt2Q2) + power2(Pi)) + (2*mt2 - 3*mu2)*mu2*(42 + 6*logmstau12Q2*(-3 +
     logmt2Q2) - 18*logmt2Q2 + 3*power2(logmstau12Q2) + 3*power2(logmt2Q2) +
     power2(Pi)) + mst12*(2*mu2*(42 - 6*logmst12Q2*(-3 + logmt2Q2) + 12*
     logmstau12Q2*(-3 + logmt2Q2) - 18*logmt2Q2 - 3*power2(logmst12Q2) + 6*
     power2(logmstau12Q2) + 3*power2(logmt2Q2) + power2(Pi)) + mt2*(-6*
     logmst12Q2*(3 + 2*logmstau12Q2 - 3*logmt2Q2) + 18*logmstau12Q2*(-1 +
     logmt2Q2) + 3*power2(logmst12Q2) + 3*power2(logmstau12Q2) + 2*(84 - 54*
     logmt2Q2 + 9*power2(logmt2Q2) + 2*power2(Pi))))) + power2(invdmstau1mu)*(
     mst14*mt2*(42 + 6*logmst12Q2*(-3 + logmt2Q2) - 18*logmt2Q2 + 3*power2(
     logmst12Q2) + 3*power2(logmt2Q2) + power2(Pi)) + mst14*mu2*(42 + 12*
     logmst12Q2*(-3 + logmt2Q2) - 6*logmstau12Q2*(-3 + logmt2Q2) - 18*logmt2Q2
     + 6*power2(logmst12Q2) - 3*power2(logmstau12Q2) + 3*power2(logmt2Q2) +
     power2(Pi)) + (mt2 - mu2)*power2(mu2)*(42 + 6*logmstau12Q2*(-3 + logmt2Q2)
     - 18*logmt2Q2 + 3*power2(logmstau12Q2) + 3*power2(logmt2Q2) + power2(Pi))
     + mst12*mu2*(mu2*(42 - 6*logmst12Q2*(-3 + logmt2Q2) + 12*logmstau12Q2*(-3
     + logmt2Q2) - 18*logmt2Q2 - 3*power2(logmst12Q2) + 6*power2(logmstau12Q2)
     + 3*power2(logmt2Q2) + power2(Pi)) + mt2*(-6*logmst12Q2*(3 + 2*
     logmstau12Q2 - 3*logmt2Q2) + 18*logmstau12Q2*(-1 + logmt2Q2) + 3*power2(
     logmst12Q2) + 3*power2(logmstau12Q2) + 2*(84 - 54*logmt2Q2 + 9*power2(
     logmt2Q2) + 2*power2(Pi)))) - (42 + 6*logmst12Q2*(-3 + logmt2Q2) - 18*
     logmt2Q2 + 3*power2(logmst12Q2) + 3*power2(logmt2Q2) + power2(Pi))*power6(
     mst1))))/8. + (mu2*DeltaInv(mstau22,mst12,mt2)*(-((mstau22 - mt2 + 2*mu2)*
     (42 + 6*logmstau22Q2*(-3 + logmt2Q2) - 18*logmt2Q2 + 3*power2(logmstau22Q2
     ) + 3*power2(logmt2Q2) + power2(Pi))) + mst12*(42 - 6*logmst12Q2*(-3 +
     logmt2Q2) + 12*logmstau22Q2*(-3 + logmt2Q2) - 18*logmt2Q2 - 3*power2(
     logmst12Q2) + 6*power2(logmstau22Q2) + 3*power2(logmt2Q2) + power2(Pi)) -
     invdmstau2mu*(mst14*(42 + 12*logmst12Q2*(-3 + logmt2Q2) - 6*logmstau22Q2*(
     -3 + logmt2Q2) - 18*logmt2Q2 + 6*power2(logmst12Q2) - 3*power2(
     logmstau22Q2) + 3*power2(logmt2Q2) + power2(Pi)) + (2*mt2 - 3*mu2)*mu2*(42
      + 6*logmstau22Q2*(-3 + logmt2Q2) - 18*logmt2Q2 + 3*power2(logmstau22Q2) +
     3*power2(logmt2Q2) + power2(Pi)) + mst12*(2*mu2*(42 - 6*logmst12Q2*(-3 +
     logmt2Q2) + 12*logmstau22Q2*(-3 + logmt2Q2) - 18*logmt2Q2 - 3*power2(
     logmst12Q2) + 6*power2(logmstau22Q2) + 3*power2(logmt2Q2) + power2(Pi)) +
     mt2*(-6*logmst12Q2*(3 + 2*logmstau22Q2 - 3*logmt2Q2) + 18*logmstau22Q2*(-1
      + logmt2Q2) + 3*power2(logmst12Q2) + 3*power2(logmstau22Q2) + 2*(84 - 54*
     logmt2Q2 + 9*power2(logmt2Q2) + 2*power2(Pi))))) + power2(invdmstau2mu)*(
     mst14*mt2*(42 + 6*logmst12Q2*(-3 + logmt2Q2) - 18*logmt2Q2 + 3*power2(
     logmst12Q2) + 3*power2(logmt2Q2) + power2(Pi)) + mst14*mu2*(42 + 12*
     logmst12Q2*(-3 + logmt2Q2) - 6*logmstau22Q2*(-3 + logmt2Q2) - 18*logmt2Q2
     + 6*power2(logmst12Q2) - 3*power2(logmstau22Q2) + 3*power2(logmt2Q2) +
     power2(Pi)) + (mt2 - mu2)*power2(mu2)*(42 + 6*logmstau22Q2*(-3 + logmt2Q2)
     - 18*logmt2Q2 + 3*power2(logmstau22Q2) + 3*power2(logmt2Q2) + power2(Pi))
     + mst12*mu2*(mu2*(42 - 6*logmst12Q2*(-3 + logmt2Q2) + 12*logmstau22Q2*(-3
     + logmt2Q2) - 18*logmt2Q2 - 3*power2(logmst12Q2) + 6*power2(logmstau22Q2)
     + 3*power2(logmt2Q2) + power2(Pi)) + mt2*(-6*logmst12Q2*(3 + 2*
     logmstau22Q2 - 3*logmt2Q2) + 18*logmstau22Q2*(-1 + logmt2Q2) + 3*power2(
     logmst12Q2) + 3*power2(logmstau22Q2) + 2*(84 - 54*logmt2Q2 + 9*power2(
     logmt2Q2) + 2*power2(Pi)))) - (42 + 6*logmst12Q2*(-3 + logmt2Q2) - 18*
     logmt2Q2 + 3*power2(logmst12Q2) + 3*power2(logmt2Q2) + power2(Pi))*power6(
     mst1))))/8. + DeltaInv(mst12,mst12,mh2)*((power2(At)*(-1 + power2(sa))*
     power2(sa)*(2*invdmst*mt2 - power2(snt) + power4(snt))*(mst14*(-72*
     logmh2Q2*(-3 + logmst12Q2) + 36*logmst12Q2 - 36*power2(logmh2Q2) + 24*
     power2(logmst12Q2) - 7*(42 + power2(Pi))) + 2*invdmhH*(42 - 36*logmst12Q2
     + 12*power2(logmst12Q2) + power2(Pi))*(-(mH2*mst14) + 4*power6(mst1))))/(
     2.*mh2) - (At*ca*mu*sa*(invdmst*mt2 - power2(snt) + power4(snt))*(mst14*(
     12*power2(logmst12Q2)*(-1 + 2*power2(sa)) - power2(Pi)*(1 + 7*power2(sa))
     + 36*logmst12Q2*(1 + (1 - 2*logmh2Q2)*power2(sa)) - 6*(7 + (49 - 36*
     logmh2Q2 + 6*power2(logmh2Q2))*power2(sa))) + invdmhH*(42 - 36*logmst12Q2
     + 12*power2(logmst12Q2) + power2(Pi))*(-1 + 2*power2(sa))*(-(mH2*mst14) +
     4*power6(mst1))))/mh2 + (power2(sa)*(mst14*(mt2*(72*logmh2Q2*(-3 +
     logmst12Q2) - 36*logmst12Q2 + 36*power2(logmh2Q2) - 24*power2(logmst12Q2)
     + 7*(42 + power2(Pi)))*(-1 + power2(sa)) + mu2*(84 - 24*power2(logmst12Q2)
     *(-1 + power2(sa)) + 6*(49 - 36*logmh2Q2 + 6*power2(logmh2Q2))*power2(sa)
     + power2(Pi)*(2 + 7*power2(sa)) + 36*logmst12Q2*(-2 + (-1 + 2*logmh2Q2)*
     power2(sa)))*(-1 + power2(snt))*power2(snt)) - 2*invdmhH*(42 - 36*
     logmst12Q2 + 12*power2(logmst12Q2) + power2(Pi))*(-1 + power2(sa))*(mt2 +
     mu2*(-1 + power2(snt))*power2(snt))*(-(mH2*mst14) + 4*power6(mst1))))/(2.*
     mh2)) + DeltaInv(mst12,mt2,mu2)*(Al*At*invdmst*invdmstau*(invdmstau1mu -
     invdmstau2mu)*mu2*(mu2*(3*mst12*mu2*(4*logmst12Q2*logmt2Q2 - 2*(-6 +
     logmu2Q2)*logmu2Q2 - 2*logmst12Q2*(3 + logmu2Q2) - 2*logmt2Q2*(3 +
     logmu2Q2) + power2(logmst12Q2) + power2(logmt2Q2)) - (mt2 - mu2)*mu2*(42 +
     6*logmt2Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*power2(logmt2Q2) + 3*power2(
     logmu2Q2) + power2(Pi)) - mst12*mt2*(294 + 36*logmt2Q2*(-5 + logmu2Q2) -
     54*logmu2Q2 + 6*logmst12Q2*(4*logmt2Q2 - 3*(1 + logmu2Q2)) + 3*power2(
     logmst12Q2) + 30*power2(logmt2Q2) + 9*power2(logmu2Q2) + 7*power2(Pi))) -
     mst14*(3*mu2*(42 + 4*logmst12Q2*(-3 + 2*logmt2Q2 - logmu2Q2) + 2*logmt2Q2*
     (-15 + logmu2Q2) + 6*logmu2Q2 + 2*power2(logmst12Q2) + 5*power2(logmt2Q2)
     - power2(logmu2Q2) + power2(Pi)) + mt2*(6*logmst12Q2*(-9 + 4*logmt2Q2 -
     logmu2Q2) + 6*logmt2Q2*(-15 + logmu2Q2) + 9*power2(logmst12Q2) + 15*power2
     (logmt2Q2) + 4*(42 + power2(Pi)))) + (6*logmst12Q2*(-3 + 2*logmt2Q2 -
     logmu2Q2) + 6*logmt2Q2*(-9 + logmu2Q2) + 3*power2(logmst12Q2) + 9*power2(
     logmt2Q2) + 2*(42 + power2(Pi)))*power6(mst1)) - (At*invdmst*invdmstau*(
     invdmstau1mu - invdmstau2mu)*mu*mu2*sb*(mu2*(3*mst12*mu2*(4*logmst12Q2*
     logmt2Q2 - 2*(-6 + logmu2Q2)*logmu2Q2 - 2*logmst12Q2*(3 + logmu2Q2) - 2*
     logmt2Q2*(3 + logmu2Q2) + power2(logmst12Q2) + power2(logmt2Q2)) - (mt2 -
     mu2)*mu2*(42 + 6*logmt2Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*power2(
     logmt2Q2) + 3*power2(logmu2Q2) + power2(Pi)) - mst12*mt2*(294 + 36*
     logmt2Q2*(-5 + logmu2Q2) - 54*logmu2Q2 + 6*logmst12Q2*(4*logmt2Q2 - 3*(1 +
     logmu2Q2)) + 3*power2(logmst12Q2) + 30*power2(logmt2Q2) + 9*power2(
     logmu2Q2) + 7*power2(Pi))) - mst14*(3*mu2*(42 + 4*logmst12Q2*(-3 + 2*
     logmt2Q2 - logmu2Q2) + 2*logmt2Q2*(-15 + logmu2Q2) + 6*logmu2Q2 + 2*power2
     (logmst12Q2) + 5*power2(logmt2Q2) - power2(logmu2Q2) + power2(Pi)) + mt2*(
     6*logmst12Q2*(-9 + 4*logmt2Q2 - logmu2Q2) + 6*logmt2Q2*(-15 + logmu2Q2) +
     9*power2(logmst12Q2) + 15*power2(logmt2Q2) + 4*(42 + power2(Pi)))) + (6*
     logmst12Q2*(-3 + 2*logmt2Q2 - logmu2Q2) + 6*logmt2Q2*(-9 + logmu2Q2) + 3*
     power2(logmst12Q2) + 9*power2(logmt2Q2) + 2*(42 + power2(Pi)))*power6(mst1
     )))/cb + (mu2*(power2(invdmstau1mu)*(mst14*mt2*(42 + 6*logmst12Q2*(-3 +
     logmt2Q2) - 18*logmt2Q2 + 3*power2(logmst12Q2) + 3*power2(logmt2Q2) +
     power2(Pi)) + mst14*mu2*(42 + 12*logmst12Q2*(-3 + logmt2Q2) + 18*logmu2Q2
     - 6*logmt2Q2*(3 + logmu2Q2) + 6*power2(logmst12Q2) + 3*power2(logmt2Q2) -
     3*power2(logmu2Q2) + power2(Pi)) + (mt2 - mu2)*power2(mu2)*(42 + 6*
     logmt2Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*power2(logmt2Q2) + 3*power2(
     logmu2Q2) + power2(Pi)) + mst12*mu2*(mu2*(42 - 6*logmst12Q2*(-3 + logmt2Q2
     ) - 36*logmu2Q2 + 6*logmt2Q2*(-3 + 2*logmu2Q2) - 3*power2(logmst12Q2) + 3*
     power2(logmt2Q2) + 6*power2(logmu2Q2) + power2(Pi)) + mt2*(168 + 6*
     logmst12Q2*(-3 + 3*logmt2Q2 - 2*logmu2Q2) + 18*logmt2Q2*(-6 + logmu2Q2) -
     18*logmu2Q2 + 3*power2(logmst12Q2) + 18*power2(logmt2Q2) + 3*power2(
     logmu2Q2) + 4*power2(Pi))) - (42 + 6*logmst12Q2*(-3 + logmt2Q2) - 18*
     logmt2Q2 + 3*power2(logmst12Q2) + 3*power2(logmt2Q2) + power2(Pi))*power6(
     mst1)) + 8*invdmst*invdmstau*invdmstau1mu*mu2*(mu2*(3*mst12*mu2*(4*
     logmst12Q2*logmt2Q2 - 2*(-6 + logmu2Q2)*logmu2Q2 - 2*logmst12Q2*(3 +
     logmu2Q2) - 2*logmt2Q2*(3 + logmu2Q2) + power2(logmst12Q2) + power2(
     logmt2Q2)) - (mt2 - mu2)*mu2*(42 + 6*logmt2Q2*(-3 + logmu2Q2) - 18*
     logmu2Q2 + 3*power2(logmt2Q2) + 3*power2(logmu2Q2) + power2(Pi)) - mst12*
     mt2*(294 + 36*logmt2Q2*(-5 + logmu2Q2) - 54*logmu2Q2 + 6*logmst12Q2*(4*
     logmt2Q2 - 3*(1 + logmu2Q2)) + 3*power2(logmst12Q2) + 30*power2(logmt2Q2)
     + 9*power2(logmu2Q2) + 7*power2(Pi))) - mst14*(3*mu2*(42 + 4*logmst12Q2*(-
     3 + 2*logmt2Q2 - logmu2Q2) + 2*logmt2Q2*(-15 + logmu2Q2) + 6*logmu2Q2 + 2*
     power2(logmst12Q2) + 5*power2(logmt2Q2) - power2(logmu2Q2) + power2(Pi)) +
     mt2*(6*logmst12Q2*(-9 + 4*logmt2Q2 - logmu2Q2) + 6*logmt2Q2*(-15 +
     logmu2Q2) + 9*power2(logmst12Q2) + 15*power2(logmt2Q2) + 4*(42 + power2(Pi
     )))) + (6*logmst12Q2*(-3 + 2*logmt2Q2 - logmu2Q2) + 6*logmt2Q2*(-9 +
     logmu2Q2) + 3*power2(logmst12Q2) + 9*power2(logmt2Q2) + 2*(42 + power2(Pi)
     ))*power6(mst1)) + invdmstau2mu*(invdmstau2mu*(mst14*mt2*(42 + 6*
     logmst12Q2*(-3 + logmt2Q2) - 18*logmt2Q2 + 3*power2(logmst12Q2) + 3*power2
     (logmt2Q2) + power2(Pi)) + mst14*mu2*(42 + 12*logmst12Q2*(-3 + logmt2Q2) +
     18*logmu2Q2 - 6*logmt2Q2*(3 + logmu2Q2) + 6*power2(logmst12Q2) + 3*power2(
     logmt2Q2) - 3*power2(logmu2Q2) + power2(Pi)) + (mt2 - mu2)*power2(mu2)*(42
      + 6*logmt2Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*power2(logmt2Q2) + 3*
     power2(logmu2Q2) + power2(Pi)) + mst12*mu2*(mu2*(42 - 6*logmst12Q2*(-3 +
     logmt2Q2) - 36*logmu2Q2 + 6*logmt2Q2*(-3 + 2*logmu2Q2) - 3*power2(
     logmst12Q2) + 3*power2(logmt2Q2) + 6*power2(logmu2Q2) + power2(Pi)) + mt2*
     (168 + 6*logmst12Q2*(-3 + 3*logmt2Q2 - 2*logmu2Q2) + 18*logmt2Q2*(-6 +
     logmu2Q2) - 18*logmu2Q2 + 3*power2(logmst12Q2) + 18*power2(logmt2Q2) + 3*
     power2(logmu2Q2) + 4*power2(Pi))) - (42 + 6*logmst12Q2*(-3 + logmt2Q2) -
     18*logmt2Q2 + 3*power2(logmst12Q2) + 3*power2(logmt2Q2) + power2(Pi))*
     power6(mst1)) - 8*invdmst*invdmstau*mu2*(mu2*(3*mst12*mu2*(4*logmst12Q2*
     logmt2Q2 - 2*(-6 + logmu2Q2)*logmu2Q2 - 2*logmst12Q2*(3 + logmu2Q2) - 2*
     logmt2Q2*(3 + logmu2Q2) + power2(logmst12Q2) + power2(logmt2Q2)) - (mt2 -
     mu2)*mu2*(42 + 6*logmt2Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*power2(
     logmt2Q2) + 3*power2(logmu2Q2) + power2(Pi)) - mst12*mt2*(294 + 36*
     logmt2Q2*(-5 + logmu2Q2) - 54*logmu2Q2 + 6*logmst12Q2*(4*logmt2Q2 - 3*(1 +
     logmu2Q2)) + 3*power2(logmst12Q2) + 30*power2(logmt2Q2) + 9*power2(
     logmu2Q2) + 7*power2(Pi))) - mst14*(3*mu2*(42 + 4*logmst12Q2*(-3 + 2*
     logmt2Q2 - logmu2Q2) + 2*logmt2Q2*(-15 + logmu2Q2) + 6*logmu2Q2 + 2*power2
     (logmst12Q2) + 5*power2(logmt2Q2) - power2(logmu2Q2) + power2(Pi)) + mt2*(
     6*logmst12Q2*(-9 + 4*logmt2Q2 - logmu2Q2) + 6*logmt2Q2*(-15 + logmu2Q2) +
     9*power2(logmst12Q2) + 15*power2(logmt2Q2) + 4*(42 + power2(Pi)))) + (6*
     logmst12Q2*(-3 + 2*logmt2Q2 - logmu2Q2) + 6*logmt2Q2*(-9 + logmu2Q2) + 3*
     power2(logmst12Q2) + 9*power2(logmt2Q2) + 2*(42 + power2(Pi)))*power6(mst1
     )))))/8.) + (mu2*DeltaInv(mst22,mstau12,mt2)*(-((mstau12 - mt2 + 2*mu2)*(
     42 + 6*logmstau12Q2*(-3 + logmt2Q2) - 18*logmt2Q2 + 3*power2(logmstau12Q2)
     + 3*power2(logmt2Q2) + power2(Pi))) + mst22*(42 - 6*logmst22Q2*(-3 +
     logmt2Q2) + 12*logmstau12Q2*(-3 + logmt2Q2) - 18*logmt2Q2 - 3*power2(
     logmst22Q2) + 6*power2(logmstau12Q2) + 3*power2(logmt2Q2) + power2(Pi)) -
     invdmstau1mu*(mst24*(42 + 12*logmst22Q2*(-3 + logmt2Q2) - 6*logmstau12Q2*(
     -3 + logmt2Q2) - 18*logmt2Q2 + 6*power2(logmst22Q2) - 3*power2(
     logmstau12Q2) + 3*power2(logmt2Q2) + power2(Pi)) + (2*mt2 - 3*mu2)*mu2*(42
      + 6*logmstau12Q2*(-3 + logmt2Q2) - 18*logmt2Q2 + 3*power2(logmstau12Q2) +
     3*power2(logmt2Q2) + power2(Pi)) + mst22*(2*mu2*(42 - 6*logmst22Q2*(-3 +
     logmt2Q2) + 12*logmstau12Q2*(-3 + logmt2Q2) - 18*logmt2Q2 - 3*power2(
     logmst22Q2) + 6*power2(logmstau12Q2) + 3*power2(logmt2Q2) + power2(Pi)) +
     mt2*(-6*logmst22Q2*(3 + 2*logmstau12Q2 - 3*logmt2Q2) + 18*logmstau12Q2*(-1
      + logmt2Q2) + 3*power2(logmst22Q2) + 3*power2(logmstau12Q2) + 2*(84 - 54*
     logmt2Q2 + 9*power2(logmt2Q2) + 2*power2(Pi))))) + power2(invdmstau1mu)*(
     mst24*mt2*(42 + 6*logmst22Q2*(-3 + logmt2Q2) - 18*logmt2Q2 + 3*power2(
     logmst22Q2) + 3*power2(logmt2Q2) + power2(Pi)) + mst24*mu2*(42 + 12*
     logmst22Q2*(-3 + logmt2Q2) - 6*logmstau12Q2*(-3 + logmt2Q2) - 18*logmt2Q2
     + 6*power2(logmst22Q2) - 3*power2(logmstau12Q2) + 3*power2(logmt2Q2) +
     power2(Pi)) + (mt2 - mu2)*power2(mu2)*(42 + 6*logmstau12Q2*(-3 + logmt2Q2)
     - 18*logmt2Q2 + 3*power2(logmstau12Q2) + 3*power2(logmt2Q2) + power2(Pi))
     + mst22*mu2*(mu2*(42 - 6*logmst22Q2*(-3 + logmt2Q2) + 12*logmstau12Q2*(-3
     + logmt2Q2) - 18*logmt2Q2 - 3*power2(logmst22Q2) + 6*power2(logmstau12Q2)
     + 3*power2(logmt2Q2) + power2(Pi)) + mt2*(-6*logmst22Q2*(3 + 2*
     logmstau12Q2 - 3*logmt2Q2) + 18*logmstau12Q2*(-1 + logmt2Q2) + 3*power2(
     logmst22Q2) + 3*power2(logmstau12Q2) + 2*(84 - 54*logmt2Q2 + 9*power2(
     logmt2Q2) + 2*power2(Pi)))) - (42 + 6*logmst22Q2*(-3 + logmt2Q2) - 18*
     logmt2Q2 + 3*power2(logmst22Q2) + 3*power2(logmt2Q2) + power2(Pi))*power6(
     mst2))))/8. + (mu2*DeltaInv(mstau22,mst22,mt2)*(-((mstau22 - mt2 + 2*mu2)*
     (42 + 6*logmstau22Q2*(-3 + logmt2Q2) - 18*logmt2Q2 + 3*power2(logmstau22Q2
     ) + 3*power2(logmt2Q2) + power2(Pi))) + mst22*(42 - 6*logmst22Q2*(-3 +
     logmt2Q2) + 12*logmstau22Q2*(-3 + logmt2Q2) - 18*logmt2Q2 - 3*power2(
     logmst22Q2) + 6*power2(logmstau22Q2) + 3*power2(logmt2Q2) + power2(Pi)) -
     invdmstau2mu*(mst24*(42 + 12*logmst22Q2*(-3 + logmt2Q2) - 6*logmstau22Q2*(
     -3 + logmt2Q2) - 18*logmt2Q2 + 6*power2(logmst22Q2) - 3*power2(
     logmstau22Q2) + 3*power2(logmt2Q2) + power2(Pi)) + (2*mt2 - 3*mu2)*mu2*(42
      + 6*logmstau22Q2*(-3 + logmt2Q2) - 18*logmt2Q2 + 3*power2(logmstau22Q2) +
     3*power2(logmt2Q2) + power2(Pi)) + mst22*(2*mu2*(42 - 6*logmst22Q2*(-3 +
     logmt2Q2) + 12*logmstau22Q2*(-3 + logmt2Q2) - 18*logmt2Q2 - 3*power2(
     logmst22Q2) + 6*power2(logmstau22Q2) + 3*power2(logmt2Q2) + power2(Pi)) +
     mt2*(-6*logmst22Q2*(3 + 2*logmstau22Q2 - 3*logmt2Q2) + 18*logmstau22Q2*(-1
      + logmt2Q2) + 3*power2(logmst22Q2) + 3*power2(logmstau22Q2) + 2*(84 - 54*
     logmt2Q2 + 9*power2(logmt2Q2) + 2*power2(Pi))))) + power2(invdmstau2mu)*(
     mst24*mt2*(42 + 6*logmst22Q2*(-3 + logmt2Q2) - 18*logmt2Q2 + 3*power2(
     logmst22Q2) + 3*power2(logmt2Q2) + power2(Pi)) + mst24*mu2*(42 + 12*
     logmst22Q2*(-3 + logmt2Q2) - 6*logmstau22Q2*(-3 + logmt2Q2) - 18*logmt2Q2
     + 6*power2(logmst22Q2) - 3*power2(logmstau22Q2) + 3*power2(logmt2Q2) +
     power2(Pi)) + (mt2 - mu2)*power2(mu2)*(42 + 6*logmstau22Q2*(-3 + logmt2Q2)
     - 18*logmt2Q2 + 3*power2(logmstau22Q2) + 3*power2(logmt2Q2) + power2(Pi))
     + mst22*mu2*(mu2*(42 - 6*logmst22Q2*(-3 + logmt2Q2) + 12*logmstau22Q2*(-3
     + logmt2Q2) - 18*logmt2Q2 - 3*power2(logmst22Q2) + 6*power2(logmstau22Q2)
     + 3*power2(logmt2Q2) + power2(Pi)) + mt2*(-6*logmst22Q2*(3 + 2*
     logmstau22Q2 - 3*logmt2Q2) + 18*logmstau22Q2*(-1 + logmt2Q2) + 3*power2(
     logmst22Q2) + 3*power2(logmstau22Q2) + 2*(84 - 54*logmt2Q2 + 9*power2(
     logmt2Q2) + 2*power2(Pi)))) - (42 + 6*logmst22Q2*(-3 + logmt2Q2) - 18*
     logmt2Q2 + 3*power2(logmst22Q2) + 3*power2(logmt2Q2) + power2(Pi))*power6(
     mst2))))/8. + ((8*mu2*(36 + 6*logmC2Q2*(-2 + logmsb12Q2) - 18*logmsb12Q2 -
     6*logmst12Q2 + 6*invdmst*logmst12Q2*mst22 - 6*invdmst*logmst22Q2*mst22 + 3
     *power2(logmC2Q2) + 3*power2(logmsb12Q2) + power2(Pi)))/mC2 - (8*mu2*(36 +
     6*logmA2Q2*(-2 + logmst12Q2) - 6*invdmst*logmst22Q2*mst22 + 6*logmst12Q2*(
     -4 + invdmst*mst22) + 3*power2(logmA2Q2) + 3*power2(logmst12Q2) + power2(
     Pi)))/mA2 + 4*invdmhH*mH2*(27 - 72*logmh2Q2 + 96*invdmhH*mst12 + 96*
     invdmhH*mst22 - 192*invdmhH*mt2 - 48*logmH2Q2*(-1 + 2*invdmhH*(mst12 +
     mst22 - 2*mt2 - mu2)) - 96*invdmhH*mu2 + 36*power2(logmh2Q2) + 24*(-1 + 2*
     invdmhH*(mst12 + mst22 - 2*mt2 - mu2))*power2(logmH2Q2) + 2*power2(Pi) + 8
     *invdmhH*mst12*power2(Pi) + 8*invdmhH*mst22*power2(Pi) - 16*invdmhH*mt2*
     power2(Pi) - 8*invdmhH*mu2*power2(Pi))*(-1 + power2(sa))*power2(sa) +
     invdmsntau2mu*mu2*(-240 + 144*logmt2Q2 - 63*invdmsntaust2*mst22 + 36*
     invdmsntaust2*logmst22Q2*mst22 - 42*mst24*power2(invdmsntaust2) - 12*mst24
     *power2(invdmsntaust2)*power2(logmst22Q2) - 24*power2(logmt2Q2) - 8*power2
     (Pi) - 4*mst24*power2(invdmsntaust2)*power2(Pi) - 63*invdmsntaust1*mst12*
     power2(snt) + 36*invdmsntaust1*logmst12Q2*mst12*power2(snt) + 63*
     invdmsntaust2*mst22*power2(snt) - 36*invdmsntaust2*logmst22Q2*mst22*power2
     (snt) - 42*mst14*power2(invdmsntaust1)*power2(snt) + 42*mst24*power2(
     invdmsntaust2)*power2(snt) - 12*mst14*power2(invdmsntaust1)*power2(
     logmst12Q2)*power2(snt) + 12*mst24*power2(invdmsntaust2)*power2(logmst22Q2
     )*power2(snt) - 4*mst14*power2(invdmsntaust1)*power2(Pi)*power2(snt) + 4*
     mst24*power2(invdmsntaust2)*power2(Pi)*power2(snt) + 12*power2(
     logmsntau2Q2)*(-2 + mst24*power2(invdmsntaust2)*(-1 + power2(snt)) - mst14
     *power2(invdmsntaust1)*power2(snt)) - 24*logmsntau2Q2*(-4 + 2*logmt2Q2 - (
     -3 + logmst22Q2)*mst24*power2(invdmsntaust2)*(-1 + power2(snt)) - 3*mst14*
     power2(invdmsntaust1)*power2(snt) + logmst12Q2*mst14*power2(invdmsntaust1)
     *power2(snt))) + (48*(-1 + power2(sa))*(power2(mt2)*(42 - 36*logmt2Q2 + 12
     *power2(logmt2Q2) + power2(Pi))*power2(sa) + mst12*(42 + 12*logmH2Q2*(-3 +
     logmst12Q2) + 6*power2(logmH2Q2) - 6*power2(logmst12Q2) + power2(Pi))*(mt2
     *power2(sa) + mu2*(-1 + power2(sa))*(-1 + power2(snt))*power2(snt)) +
     mst22*(42 + 12*logmH2Q2*(-3 + logmst22Q2) + 6*power2(logmH2Q2) - 6*power2(
     logmst22Q2) + power2(Pi))*(mt2*power2(sa) + mu2*(-1 + power2(sa))*(-1 +
     power2(snt))*power2(snt))))/power2(mH2) + 24*(logmsntau2Q2 - logmu2Q2)*
     power2(mu2)*((-4 + 2*logmsb12Q2 + logmsntau2Q2 + logmu2Q2)*msb12 + 4*mt2 -
     logmsntau2Q2*mt2 - 2*logmt2Q2*mt2 - logmu2Q2*mt2 + 4*mu2 - 2*logmsntau2Q2*
     mu2 - 2*logmu2Q2*mu2 - (-4 + logmsntau2Q2 + 2*logmst22Q2 + logmu2Q2)*mst22
     *(-1 + power2(snt)) - 4*mst12*power2(snt) + logmsntau2Q2*mst12*power2(snt)
     + 2*logmst12Q2*mst12*power2(snt) + logmu2Q2*mst12*power2(snt))*power3(
     invdmsntau2mu) + 4*(4*invdmstau2mu*mu2*(-36 + logmstau22Q2*(15 - 6*
     logmt2Q2) + 336*invdmst*invdmstau*mst12*mu2 - 72*invdmst*invdmstau*
     logmst12Q2*mst12*mu2 - 24*invdmst*invdmstau*logmst12Q2*logmu2Q2*mst12*mu2
     - 336*invdmst*invdmstau*mst22*mu2 + 72*invdmst*invdmstau*logmst22Q2*mst22*
     mu2 + 24*invdmst*invdmstau*logmst22Q2*logmu2Q2*mst22*mu2 + 6*logmt2Q2*(3 +
     4*invdmst*invdmstau*((-9 + 2*logmst12Q2 + logmu2Q2)*mst12 - (-9 + 2*
     logmst22Q2 + logmu2Q2)*mst22)*mu2) + 12*invdmst*invdmstau*mst12*mu2*power2
     (logmst12Q2) - 12*invdmst*invdmstau*mst22*mu2*power2(logmst22Q2) - 3*
     power2(logmstau22Q2) + (-3 + 36*invdmst*invdmstau*(mst12 - mst22)*mu2)*
     power2(logmt2Q2) - power2(Pi) + 8*invdmst*invdmstau*mst12*mu2*power2(Pi) -
     8*invdmst*invdmstau*mst22*mu2*power2(Pi)) - 4*invdmstau1mu*mu2*(36 + 3*
     logmstau12Q2*(-5 + 2*logmt2Q2) + 336*invdmst*invdmstau*mst12*mu2 - 72*
     invdmst*invdmstau*logmst12Q2*mst12*mu2 - 24*invdmst*invdmstau*logmst12Q2*
     logmu2Q2*mst12*mu2 - 336*invdmst*invdmstau*mst22*mu2 + 72*invdmst*
     invdmstau*logmst22Q2*mst22*mu2 + 24*invdmst*invdmstau*logmst22Q2*logmu2Q2*
     mst22*mu2 + 6*logmt2Q2*(-3 + 4*invdmst*invdmstau*((-9 + 2*logmst12Q2 +
     logmu2Q2)*mst12 - (-9 + 2*logmst22Q2 + logmu2Q2)*mst22)*mu2) + 12*invdmst*
     invdmstau*mst12*mu2*power2(logmst12Q2) - 12*invdmst*invdmstau*mst22*mu2*
     power2(logmst22Q2) + 3*power2(logmstau12Q2) + (3 + 36*invdmst*invdmstau*(
     mst12 - mst22)*mu2)*power2(logmt2Q2) + power2(Pi) + 8*invdmst*invdmstau*
     mst12*mu2*power2(Pi) - 8*invdmst*invdmstau*mst22*mu2*power2(Pi)) + 2*mu2*
     power2(invdmstau1mu)*(mst12*(48 + 12*logmstau12Q2 - 6*logmst12Q2*(2 +
     logmstau12Q2 - 2*logmt2Q2) - 36*logmt2Q2 + 3*power2(logmst12Q2) - 3*power2
     (logmstau12Q2) + 6*power2(logmt2Q2) + power2(Pi)) + mst22*(48 + 12*
     logmstau12Q2 - 6*logmst22Q2*(2 + logmstau12Q2 - 2*logmt2Q2) - 36*logmt2Q2
     + 3*power2(logmst22Q2) - 3*power2(logmstau12Q2) + 6*power2(logmt2Q2) +
     power2(Pi)) + 2*(mt2*(36 + 6*logmstau12Q2*(-2 + logmt2Q2) - 24*logmt2Q2 +
     3*power2(logmstau12Q2) + 3*power2(logmt2Q2) + power2(Pi)) + mu2*(72 + 3*
     logmstau12Q2*(-7 + 2*logmt2Q2) + 6*logmt2Q2*(-6 + logmu2Q2) - 9*logmu2Q2 +
     6*power2(logmstau12Q2) + 6*power2(logmt2Q2) + 2*power2(Pi)))) + 2*mu2*
     power2(invdmstau2mu)*(mst12*(48 + 12*logmstau22Q2 - 6*logmst12Q2*(2 +
     logmstau22Q2 - 2*logmt2Q2) - 36*logmt2Q2 + 3*power2(logmst12Q2) - 3*power2
     (logmstau22Q2) + 6*power2(logmt2Q2) + power2(Pi)) + mst22*(48 + 12*
     logmstau22Q2 - 6*logmst22Q2*(2 + logmstau22Q2 - 2*logmt2Q2) - 36*logmt2Q2
     + 3*power2(logmst22Q2) - 3*power2(logmstau22Q2) + 6*power2(logmt2Q2) +
     power2(Pi)) + 2*(mt2*(36 + 6*logmstau22Q2*(-2 + logmt2Q2) - 24*logmt2Q2 +
     3*power2(logmstau22Q2) + 3*power2(logmt2Q2) + power2(Pi)) + mu2*(72 + 3*
     logmstau22Q2*(-7 + 2*logmt2Q2) + 6*logmt2Q2*(-6 + logmu2Q2) - 9*logmu2Q2 +
     6*power2(logmstau22Q2) + 6*power2(logmt2Q2) + 2*power2(Pi)))) - 3*(-48 +
     42*invdmhH*mst12 - 24*invdmhH*logmst12Q2*mst12 + 42*invdmhH*mst22 - 24*
     invdmhH*logmst22Q2*mst22 - 98*invdmhH*mt2 + 56*invdmhH*logmt2Q2*mt2 - 4*
     logmH2Q2*(-6 + 5*invdmhH*(mst12 + mst22 - 2*mt2 - mu2)) + 27*invdmhH*mu2 -
     20*invdmhH*logmst12Q2*mu2 + 20*invdmhH*invdmst*logmst12Q2*mst22*mu2 - 20*
     invdmhH*invdmst*logmst22Q2*mst22*mu2 + 24*invdmhH*logmh2Q2*(-2 +
     logmst12Q2 - invdmst*logmst12Q2*mst22 + invdmst*logmst22Q2*mst22)*mu2 + 6*
     (1 + 2*invdmhH*mu2)*power2(logmh2Q2) + 6*(-1 + 2*invdmhH*(mst12 + mst22 -
     2*mt2 - mu2))*power2(logmH2Q2) + 12*invdmhH*mst12*power2(logmst12Q2) + 12*
     invdmhH*mst22*power2(logmst22Q2) - 24*invdmhH*mt2*power2(logmt2Q2) + 4*
     invdmhH*mst12*power2(Pi) + 4*invdmhH*mst22*power2(Pi) - 8*invdmhH*mt2*
     power2(Pi))*(-1 + power2(sa))*power2(sa) + 6*(logmstau12Q2 - logmu2Q2)*((-
     4 + 2*logmst12Q2 + logmstau12Q2 + logmu2Q2)*mst12 + (-4 + 2*logmst22Q2 +
     logmstau12Q2 + logmu2Q2)*mst22 - 2*((-4 + logmstau12Q2 + 2*logmt2Q2 +
     logmu2Q2)*mt2 + (-2 + logmstau12Q2 + logmu2Q2)*mu2))*power2(mu2)*power3(
     invdmstau1mu) + 6*(logmstau22Q2 - logmu2Q2)*((-4 + 2*logmst12Q2 +
     logmstau22Q2 + logmu2Q2)*mst12 + (-4 + 2*logmst22Q2 + logmstau22Q2 +
     logmu2Q2)*mst22 - 2*((-4 + logmstau22Q2 + 2*logmt2Q2 + logmu2Q2)*mt2 + (-2
      + logmstau22Q2 + logmu2Q2)*mu2))*power2(mu2)*power3(invdmstau2mu)) + (24*
     (-1 + power2(sa))*(2*mt2*(42 + 6*logmH2Q2*(-3 + logmt2Q2) - 18*logmt2Q2 +
     3*power2(logmH2Q2) + 3*power2(logmt2Q2) + power2(Pi))*power2(sa) + mu2*(-1
      + power2(sa))*(36 - 6*invdmst*logmst22Q2*mst22 + power2(Pi) - 168*power2(
     snt) - 4*power2(Pi)*power2(snt) + 3*power2(logmH2Q2)*power2(1 - 2*power2(
     snt)) + 3*power2(logmst12Q2)*power2(1 - 2*power2(snt)) + 6*logmst12Q2*(-4
     + invdmst*mst22 + 12*power2(snt) - 12*power4(snt)) + 168*power4(snt) + 4*
     power2(Pi)*power4(snt) + 6*logmH2Q2*(logmst12Q2*power2(1 - 2*power2(snt))
     - 2*(1 - 6*power2(snt) + 6*power4(snt))))))/mH2 + (4*power2(sa)*(-420*
     mst12*mt2 + 432*logmh2Q2*mst12*mt2 - 96*logmst12Q2*mst12*mt2 - 144*
     logmh2Q2*logmst12Q2*mst12*mt2 + 672*invdmhH*mst14*mt2 - 576*invdmhH*
     logmst12Q2*mst14*mt2 - 420*mst22*mt2 + 432*logmh2Q2*mst22*mt2 - 96*
     logmst22Q2*mst22*mt2 - 144*logmh2Q2*logmst22Q2*mst22*mt2 + 672*invdmhH*
     mst24*mt2 - 576*invdmhH*logmst22Q2*mst24*mt2 + 21*mst12*mu2 - 12*
     logmst12Q2*mst12*mu2 + 21*mst22*mu2 - 12*logmst22Q2*mst22*mu2 - 72*mst12*
     mt2*power2(logmh2Q2) - 72*mst22*mt2*power2(logmh2Q2) + 120*mst12*mt2*
     power2(logmst12Q2) + 192*invdmhH*mst14*mt2*power2(logmst12Q2) + 120*mst22*
     mt2*power2(logmst22Q2) + 192*invdmhH*mst24*mt2*power2(logmst22Q2) - 774*
     power2(mt2) + 744*logmt2Q2*power2(mt2) - 288*power2(logmt2Q2)*power2(mt2)
     - 8*mst12*mt2*power2(Pi) + 16*invdmhH*mst14*mt2*power2(Pi) - 8*mst22*mt2*
     power2(Pi) + 16*invdmhH*mst24*mt2*power2(Pi) - 24*power2(mt2)*power2(Pi) -
     6*logmH2Q2*mH2*(-1 + invdmhH*mH2)*(mH2 + 6*invdmhH*mH2*(mst12 + mst22 - 2*
     mt2 - mu2) - 2*(3*mst12 - 6*mt2 + mu2 - 2*logmst12Q2*mu2 + mst22*(3 + 2*
     invdmst*(logmst12Q2 - logmst22Q2)*mu2)))*(-1 + power2(sa)) + 6*mH2*(-1 +
     invdmhH*mH2)*(mH2 - 2*(mst12 + mst22 - 2*mt2) + 2*invdmhH*mH2*(mst12 +
     mst22 - 2*mt2 - mu2))*power2(logmH2Q2)*(-1 + power2(sa)) - power2(mH2)*(
     power2(Pi) + invdmhH*(-174*mt2 + 24*logmt2Q2*mt2 - 27*mu2 - 36*logmst12Q2*
     mu2 - 24*mt2*power2(logmt2Q2) - 12*mt2*power2(Pi) - 2*mu2*power2(Pi) + 6*
     mst12*(18 - 4*logmst12Q2 + 2*power2(logmst12Q2) + power2(Pi)) + 6*mst22*(
     18 - 4*logmst22Q2 + 6*invdmst*logmst12Q2*mu2 - 6*invdmst*logmst22Q2*mu2 +
     2*power2(logmst22Q2) + power2(Pi))))*(-1 + power2(sa)) + 420*mst12*mt2*
     power2(sa) - 432*logmh2Q2*mst12*mt2*power2(sa) + 96*logmst12Q2*mst12*mt2*
     power2(sa) + 144*logmh2Q2*logmst12Q2*mst12*mt2*power2(sa) - 672*invdmhH*
     mst14*mt2*power2(sa) + 576*invdmhH*logmst12Q2*mst14*mt2*power2(sa) + 420*
     mst22*mt2*power2(sa) - 432*logmh2Q2*mst22*mt2*power2(sa) + 96*logmst22Q2*
     mst22*mt2*power2(sa) + 144*logmh2Q2*logmst22Q2*mst22*mt2*power2(sa) - 672*
     invdmhH*mst24*mt2*power2(sa) + 576*invdmhH*logmst22Q2*mst24*mt2*power2(sa)
     - 21*mst12*mu2*power2(sa) + 12*logmst12Q2*mst12*mu2*power2(sa) - 21*mst22*
     mu2*power2(sa) + 12*logmst22Q2*mst22*mu2*power2(sa) + 72*mst12*mt2*power2(
     logmh2Q2)*power2(sa) + 72*mst22*mt2*power2(logmh2Q2)*power2(sa) - 120*
     mst12*mt2*power2(logmst12Q2)*power2(sa) - 192*invdmhH*mst14*mt2*power2(
     logmst12Q2)*power2(sa) - 120*mst22*mt2*power2(logmst22Q2)*power2(sa) - 192
     *invdmhH*mst24*mt2*power2(logmst22Q2)*power2(sa) + 774*power2(mt2)*power2(
     sa) - 744*logmt2Q2*power2(mt2)*power2(sa) + 288*power2(logmt2Q2)*power2(
     mt2)*power2(sa) + 8*mst12*mt2*power2(Pi)*power2(sa) - 16*invdmhH*mst14*mt2
     *power2(Pi)*power2(sa) + 8*mst22*mt2*power2(Pi)*power2(sa) - 16*invdmhH*
     mst24*mt2*power2(Pi)*power2(sa) + 24*power2(mt2)*power2(Pi)*power2(sa) -
     168*mst12*mu2*power2(snt) + 144*logmst12Q2*mst12*mu2*power2(snt) - 672*
     invdmhH*mst14*mu2*power2(snt) + 576*invdmhH*logmst12Q2*mst14*mu2*power2(
     snt) - 168*mst22*mu2*power2(snt) + 144*logmst22Q2*mst22*mu2*power2(snt) -
     672*invdmhH*mst24*mu2*power2(snt) + 576*invdmhH*logmst22Q2*mst24*mu2*
     power2(snt) - 48*mst12*mu2*power2(logmst12Q2)*power2(snt) - 192*invdmhH*
     mst14*mu2*power2(logmst12Q2)*power2(snt) - 48*mst22*mu2*power2(logmst22Q2)
     *power2(snt) - 192*invdmhH*mst24*mu2*power2(logmst22Q2)*power2(snt) - 4*
     mst12*mu2*power2(Pi)*power2(snt) - 16*invdmhH*mst14*mu2*power2(Pi)*power2(
     snt) - 4*mst22*mu2*power2(Pi)*power2(snt) - 16*invdmhH*mst24*mu2*power2(Pi
     )*power2(snt) - 336*mst12*mu2*power2(sa)*power2(snt) + 432*logmh2Q2*mst12*
     mu2*power2(sa)*power2(snt) - 144*logmst12Q2*mst12*mu2*power2(sa)*power2(
     snt) - 144*logmh2Q2*logmst12Q2*mst12*mu2*power2(sa)*power2(snt) + 672*
     invdmhH*mst14*mu2*power2(sa)*power2(snt) - 576*invdmhH*logmst12Q2*mst14*
     mu2*power2(sa)*power2(snt) - 336*mst22*mu2*power2(sa)*power2(snt) + 432*
     logmh2Q2*mst22*mu2*power2(sa)*power2(snt) - 144*logmst22Q2*mst22*mu2*
     power2(sa)*power2(snt) - 144*logmh2Q2*logmst22Q2*mst22*mu2*power2(sa)*
     power2(snt) + 672*invdmhH*mst24*mu2*power2(sa)*power2(snt) - 576*invdmhH*
     logmst22Q2*mst24*mu2*power2(sa)*power2(snt) - 72*mst12*mu2*power2(logmh2Q2
     )*power2(sa)*power2(snt) - 72*mst22*mu2*power2(logmh2Q2)*power2(sa)*power2
     (snt) + 120*mst12*mu2*power2(logmst12Q2)*power2(sa)*power2(snt) + 192*
     invdmhH*mst14*mu2*power2(logmst12Q2)*power2(sa)*power2(snt) + 120*mst22*
     mu2*power2(logmst22Q2)*power2(sa)*power2(snt) + 192*invdmhH*mst24*mu2*
     power2(logmst22Q2)*power2(sa)*power2(snt) - 8*mst12*mu2*power2(Pi)*power2(
     sa)*power2(snt) + 16*invdmhH*mst14*mu2*power2(Pi)*power2(sa)*power2(snt) -
     8*mst22*mu2*power2(Pi)*power2(sa)*power2(snt) + 16*invdmhH*mst24*mu2*
     power2(Pi)*power2(sa)*power2(snt) + invdmhH*(power2(Pi) + 2*invdmhH*(mst12
      + mst22 - 2*mt2 - mu2)*(21 + power2(Pi)))*(-1 + power2(sa))*power3(mH2) -
     1344*invdmhH*power3(mt2) + 1152*invdmhH*logmt2Q2*power3(mt2) - 384*invdmhH
     *power2(logmt2Q2)*power3(mt2) - 32*invdmhH*power2(Pi)*power3(mt2) + 1344*
     invdmhH*power2(sa)*power3(mt2) - 1152*invdmhH*logmt2Q2*power2(sa)*power3(
     mt2) + 384*invdmhH*power2(logmt2Q2)*power2(sa)*power3(mt2) + 32*invdmhH*
     power2(Pi)*power2(sa)*power3(mt2) + 168*mst12*mu2*power4(snt) - 144*
     logmst12Q2*mst12*mu2*power4(snt) + 672*invdmhH*mst14*mu2*power4(snt) - 576
     *invdmhH*logmst12Q2*mst14*mu2*power4(snt) + 168*mst22*mu2*power4(snt) -
     144*logmst22Q2*mst22*mu2*power4(snt) + 672*invdmhH*mst24*mu2*power4(snt) -
     576*invdmhH*logmst22Q2*mst24*mu2*power4(snt) + 48*mst12*mu2*power2(
     logmst12Q2)*power4(snt) + 192*invdmhH*mst14*mu2*power2(logmst12Q2)*power4(
     snt) + 48*mst22*mu2*power2(logmst22Q2)*power4(snt) + 192*invdmhH*mst24*mu2
     *power2(logmst22Q2)*power4(snt) + 4*mst12*mu2*power2(Pi)*power4(snt) + 16*
     invdmhH*mst14*mu2*power2(Pi)*power4(snt) + 4*mst22*mu2*power2(Pi)*power4(
     snt) + 16*invdmhH*mst24*mu2*power2(Pi)*power4(snt) + 336*mst12*mu2*power2(
     sa)*power4(snt) - 432*logmh2Q2*mst12*mu2*power2(sa)*power4(snt) + 144*
     logmst12Q2*mst12*mu2*power2(sa)*power4(snt) + 144*logmh2Q2*logmst12Q2*
     mst12*mu2*power2(sa)*power4(snt) - 672*invdmhH*mst14*mu2*power2(sa)*power4
     (snt) + 576*invdmhH*logmst12Q2*mst14*mu2*power2(sa)*power4(snt) + 336*
     mst22*mu2*power2(sa)*power4(snt) - 432*logmh2Q2*mst22*mu2*power2(sa)*
     power4(snt) + 144*logmst22Q2*mst22*mu2*power2(sa)*power4(snt) + 144*
     logmh2Q2*logmst22Q2*mst22*mu2*power2(sa)*power4(snt) - 672*invdmhH*mst24*
     mu2*power2(sa)*power4(snt) + 576*invdmhH*logmst22Q2*mst24*mu2*power2(sa)*
     power4(snt) + 72*mst12*mu2*power2(logmh2Q2)*power2(sa)*power4(snt) + 72*
     mst22*mu2*power2(logmh2Q2)*power2(sa)*power4(snt) - 120*mst12*mu2*power2(
     logmst12Q2)*power2(sa)*power4(snt) - 192*invdmhH*mst14*mu2*power2(
     logmst12Q2)*power2(sa)*power4(snt) - 120*mst22*mu2*power2(logmst22Q2)*
     power2(sa)*power4(snt) - 192*invdmhH*mst24*mu2*power2(logmst22Q2)*power2(
     sa)*power4(snt) + 8*mst12*mu2*power2(Pi)*power2(sa)*power4(snt) - 16*
     invdmhH*mst14*mu2*power2(Pi)*power2(sa)*power4(snt) + 8*mst22*mu2*power2(
     Pi)*power2(sa)*power4(snt) - 16*invdmhH*mst24*mu2*power2(Pi)*power2(sa)*
     power4(snt) + mH2*(-1 + power2(sa))*(-90*mt2 + 24*logmt2Q2*mt2 + 15*mu2 -
     36*logmst12Q2*mu2 - 24*mt2*power2(logmt2Q2) - 270*invdmhH*power2(mt2) +
     312*invdmhH*logmt2Q2*power2(mt2) - 144*invdmhH*power2(logmt2Q2)*power2(mt2
     ) - 8*mt2*power2(Pi) - 12*invdmhH*power2(mt2)*power2(Pi) + mst12*(66 + 4*
     power2(Pi) + 4*invdmhH*mt2*(21 + power2(Pi)) + 12*power2(logmst12Q2)*(1 +
     4*invdmhH*(mt2 + mu2*(-1 + power2(snt))*power2(snt))) + invdmhH*mu2*(21 -
     4*(42 + power2(Pi))*power2(snt) + 4*(42 + power2(Pi))*power4(snt)) - 12*
     logmst12Q2*(2 + invdmhH*(8*mt2 + mu2 - 12*mu2*power2(snt) + 12*mu2*power4(
     snt)))) + mst22*(66 + 36*invdmst*logmst12Q2*mu2 + 4*power2(Pi) + 4*invdmhH
     *mt2*(21 + power2(Pi)) + 12*power2(logmst22Q2)*(1 + 4*invdmhH*(mt2 + mu2*(
     -1 + power2(snt))*power2(snt))) + invdmhH*mu2*(21 - 4*(42 + power2(Pi))*
     power2(snt) + 4*(42 + power2(Pi))*power4(snt)) - 12*logmst22Q2*(2 + 3*
     invdmst*mu2 + invdmhH*(8*mt2 + mu2 - 12*mu2*power2(snt) + 12*mu2*power4(
     snt)))))))/power2(mh2) - (4*power2(sa)*(-(power2(mH2)*(27 - 24*logmH2Q2*(1
      + 4*invdmhH*(2*mst12 - 4*mt2 - logmst12Q2*mu2 + mst22*(2 + invdmst*(
     logmst12Q2 - logmst22Q2)*mu2))) + 12*(1 + 4*invdmhH*(2*mst12 + 2*mst22 - 4
     *mt2 - mu2))*power2(logmH2Q2) + 2*power2(Pi) + 8*invdmhH*(-72*mt2 + 24*
     logmt2Q2*mt2 - 12*logmst12Q2*mu2 - 12*mt2*power2(logmt2Q2) - 6*mt2*power2(
     Pi) - mu2*power2(Pi) + 3*mst12*(12 - 4*logmst12Q2 + 2*power2(logmst12Q2) +
     power2(Pi)) + 3*mst22*(12 + 4*invdmst*logmst12Q2*mu2 - 4*logmst22Q2*(1 +
     invdmst*mu2) + 2*power2(logmst22Q2) + power2(Pi))))*(-1 + power2(sa))) -
     24*(-1 + power2(sa))*(-2*power2(mt2)*(42 - 36*logmt2Q2 + 12*power2(
     logmt2Q2) + power2(Pi)) + mst12*(42 - 36*logmst12Q2 + 12*power2(logmst12Q2
     ) + power2(Pi))*(mt2 + mu2*(-1 + power2(snt))*power2(snt)) + mst22*(42 -
     36*logmst22Q2 + 12*power2(logmst22Q2) + power2(Pi))*(mt2 + mu2*(-1 +
     power2(snt))*power2(snt))) + invdmhH*(-6*logmH2Q2*(5 + 22*invdmhH*(mst12 +
     mst22 - 2*mt2 - mu2)) + 6*(3 + 10*invdmhH*(mst12 + mst22 - 2*mt2 - mu2))*
     power2(logmH2Q2) + 3*(9 + power2(Pi)) + 2*invdmhH*(mst12 + mst22 - 2*mt2 -
     mu2)*(69 + 5*power2(Pi)))*(-1 + power2(sa))*power3(mH2) + mH2*(-3*(2*
     invdmhH*power2(mt2)*(45 - 52*logmt2Q2 + 24*power2(logmt2Q2) + 2*power2(Pi)
     )*(-1 + power2(sa)) + 2*mt2*(133 - 20*logmH2Q2 + 12*logmh2Q2*(-3 +
     logmt2Q2) - 64*logmt2Q2 + 6*power2(logmh2Q2) + 12*power2(logmH2Q2) + 18*
     power2(logmt2Q2) + 6*power2(Pi))*(-1 + power2(sa)) + mu2*(27 - 4*logmH2Q2*
     (-7 + 6*logmst12Q2)*(-1 + power2(sa)) + 45*power2(sa) - 24*logmh2Q2*power2
     (sa) + 6*power2(logmh2Q2)*power2(sa) + 2*power2(Pi)*power2(sa) - 336*
     power2(sa)*power2(snt) + 144*logmh2Q2*power2(sa)*power2(snt) - 24*power2(
     logmh2Q2)*power2(sa)*power2(snt) - 8*power2(Pi)*power2(sa)*power2(snt) + 6
     *power2(logmst12Q2)*power2(sa)*power2(1 - 2*power2(snt)) + 4*logmst12Q2*(-
     5 + power2(sa)*(-7 + 36*power2(snt) + 3*logmh2Q2*power2(1 - 2*power2(snt))
     - 36*power4(snt))) + 336*power2(sa)*power4(snt) - 144*logmh2Q2*power2(sa)*
     power4(snt) + 24*power2(logmh2Q2)*power2(sa)*power4(snt) + 8*power2(Pi)*
     power2(sa)*power4(snt))) + mst12*(-1 + power2(sa))*(126 - 60*logmH2Q2 + 84
     *invdmhH*mt2 + 21*invdmhH*mu2 + 36*power2(logmH2Q2) + 12*power2(Pi) + 4*
     invdmhH*mt2*power2(Pi) - 168*invdmhH*mu2*power2(snt) - 4*invdmhH*mu2*
     power2(Pi)*power2(snt) + 12*power2(logmst12Q2)*(3 + 4*invdmhH*(mt2 + mu2*(
     -1 + power2(snt))*power2(snt))) + 168*invdmhH*mu2*power4(snt) + 4*invdmhH*
     mu2*power2(Pi)*power4(snt) - 12*logmst12Q2*(6 + invdmhH*(8*mt2 + mu2 - 12*
     mu2*power2(snt) + 12*mu2*power4(snt)))) + mst22*(-126 - 84*invdmhH*mt2 -
     21*invdmhH*mu2 - 60*invdmst*logmst12Q2*mu2 - 12*power2(Pi) - 4*invdmhH*mt2
     *power2(Pi) - 12*logmH2Q2*(5 + 6*invdmst*(logmst12Q2 - logmst22Q2)*mu2)*(-
     1 + power2(sa)) + 36*power2(logmH2Q2)*(-1 + power2(sa)) + 126*power2(sa) +
     84*invdmhH*mt2*power2(sa) + 21*invdmhH*mu2*power2(sa) + 24*invdmst*
     logmst12Q2*mu2*power2(sa) + 12*power2(Pi)*power2(sa) + 4*invdmhH*mt2*
     power2(Pi)*power2(sa) + 168*invdmhH*mu2*power2(snt) + 4*invdmhH*mu2*power2
     (Pi)*power2(snt) - 168*invdmhH*mu2*power2(sa)*power2(snt) - 4*invdmhH*mu2*
     power2(Pi)*power2(sa)*power2(snt) + 12*power2(logmst22Q2)*(-1 + power2(sa)
     )*(3 + 4*invdmhH*(mt2 + mu2*(-1 + power2(snt))*power2(snt))) - 168*invdmhH
     *mu2*power4(snt) - 4*invdmhH*mu2*power2(Pi)*power4(snt) + 168*invdmhH*mu2*
     power2(sa)*power4(snt) + 4*invdmhH*mu2*power2(Pi)*power2(sa)*power4(snt) -
     12*logmst22Q2*(-6 - 5*invdmst*mu2 + 6*power2(sa) + 2*invdmst*mu2*power2(sa
     ) + invdmhH*(-1 + power2(sa))*(8*mt2 + mu2 - 12*mu2*power2(snt) + 12*mu2*
     power4(snt)))))))/(mh2*mH2) + mu2*power2(invdmsntau2mu)*(-225*mst22 + 96*
     logmsntau2Q2*mst22 + 156*logmst22Q2*mst22 - 48*logmsntau2Q2*logmst22Q2*
     mst22 + 288*mt2 - 96*logmsntau2Q2*mt2 - 192*logmt2Q2*mt2 + 48*logmsntau2Q2
     *logmt2Q2*mt2 + 480*mu2 - 192*logmsntau2Q2*mu2 - 288*logmt2Q2*mu2 + 48*
     logmsntau2Q2*logmt2Q2*mu2 + 48*logmt2Q2*logmu2Q2*mu2 - 24*mst22*power2(
     logmsntau2Q2) + 24*mt2*power2(logmsntau2Q2) + 72*mu2*power2(logmsntau2Q2)
     - 24*mst22*power2(logmst22Q2) + 24*mt2*power2(logmt2Q2) + 48*mu2*power2(
     logmt2Q2) - 24*mu2*power2(logmu2Q2) - 8*mst22*power2(Pi) + 8*mt2*power2(Pi
     ) + 16*mu2*power2(Pi) + 8*msb12*(48 + 12*logmsntau2Q2 - 6*logmsb12Q2*(2 +
     logmsntau2Q2 - 2*logmt2Q2) - 36*logmt2Q2 + 3*power2(logmsb12Q2) - 3*power2
     (logmsntau2Q2) + 6*power2(logmt2Q2) + power2(Pi)) - invdmsntaust2*(9*(7 -
     4*logmst22Q2)*mst22*mu2 + mst24*(-21 + 24*logmsntau2Q2*(-3 + logmst22Q2) +
     36*logmst22Q2 + 12*power2(logmsntau2Q2) + 12*power2(logmst22Q2) + 4*power2
     (Pi)))*(-1 + power2(snt)) - 225*mst12*power2(snt) + 96*logmsntau2Q2*mst12*
     power2(snt) + 156*logmst12Q2*mst12*power2(snt) - 48*logmsntau2Q2*
     logmst12Q2*mst12*power2(snt) - 21*invdmsntaust1*mst14*power2(snt) - 72*
     invdmsntaust1*logmsntau2Q2*mst14*power2(snt) + 36*invdmsntaust1*logmst12Q2
     *mst14*power2(snt) + 24*invdmsntaust1*logmsntau2Q2*logmst12Q2*mst14*power2
     (snt) + 225*mst22*power2(snt) - 96*logmsntau2Q2*mst22*power2(snt) - 156*
     logmst22Q2*mst22*power2(snt) + 48*logmsntau2Q2*logmst22Q2*mst22*power2(snt
     ) + 63*invdmsntaust1*mst12*mu2*power2(snt) - 36*invdmsntaust1*logmst12Q2*
     mst12*mu2*power2(snt) + 42*mst14*mu2*power2(invdmsntaust1)*power2(snt) -
     72*logmsntau2Q2*mst14*mu2*power2(invdmsntaust1)*power2(snt) + 24*
     logmsntau2Q2*logmst12Q2*mst14*mu2*power2(invdmsntaust1)*power2(snt) - 24*
     mst12*power2(logmsntau2Q2)*power2(snt) + 12*invdmsntaust1*mst14*power2(
     logmsntau2Q2)*power2(snt) + 24*mst22*power2(logmsntau2Q2)*power2(snt) + 12
     *mst14*mu2*power2(invdmsntaust1)*power2(logmsntau2Q2)*power2(snt) - 24*
     mst12*power2(logmst12Q2)*power2(snt) + 12*invdmsntaust1*mst14*power2(
     logmst12Q2)*power2(snt) + 12*mst14*mu2*power2(invdmsntaust1)*power2(
     logmst12Q2)*power2(snt) + 24*mst22*power2(logmst22Q2)*power2(snt) - 8*
     mst12*power2(Pi)*power2(snt) + 4*invdmsntaust1*mst14*power2(Pi)*power2(snt
     ) + 8*mst22*power2(Pi)*power2(snt) + 4*mst14*mu2*power2(invdmsntaust1)*
     power2(Pi)*power2(snt) - 42*power2(invdmsntaust1)*power2(snt)*power6(mst1)
     + 72*logmsntau2Q2*power2(invdmsntaust1)*power2(snt)*power6(mst1) - 24*
     logmsntau2Q2*logmst12Q2*power2(invdmsntaust1)*power2(snt)*power6(mst1) -
     12*power2(invdmsntaust1)*power2(logmsntau2Q2)*power2(snt)*power6(mst1) -
     12*power2(invdmsntaust1)*power2(logmst12Q2)*power2(snt)*power6(mst1) - 4*
     power2(invdmsntaust1)*power2(Pi)*power2(snt)*power6(mst1) + 2*power2(
     invdmsntaust2)*(21 + 12*logmsntau2Q2*(-3 + logmst22Q2) + 6*power2(
     logmsntau2Q2) + 6*power2(logmst22Q2) + 2*power2(Pi))*(-1 + power2(snt))*(-
     (mst24*mu2) + power6(mst2))))/64. + DeltaInv(mst22,mst22,mh2)*(-(power2(At
     )*(-1 + power2(sa))*power2(sa)*(2*invdmst*mt2 + power2(snt) - power4(snt))
     *(mst24*(-72*logmh2Q2*(-3 + logmst22Q2) + 36*logmst22Q2 - 36*power2(
     logmh2Q2) + 24*power2(logmst22Q2) - 7*(42 + power2(Pi))) + 2*invdmhH*(42 -
     36*logmst22Q2 + 12*power2(logmst22Q2) + power2(Pi))*(-(mH2*mst24) + 4*
     power6(mst2))))/(2.*mh2) + (At*ca*mu*sa*(invdmst*mt2 + power2(snt) -
     power4(snt))*(mst24*(12*power2(logmst22Q2)*(-1 + 2*power2(sa)) - power2(Pi
     )*(1 + 7*power2(sa)) + 36*logmst22Q2*(1 + (1 - 2*logmh2Q2)*power2(sa)) - 6
     *(7 + (49 - 36*logmh2Q2 + 6*power2(logmh2Q2))*power2(sa))) + invdmhH*(42 -
     36*logmst22Q2 + 12*power2(logmst22Q2) + power2(Pi))*(-1 + 2*power2(sa))*(-
     (mH2*mst24) + 4*power6(mst2))))/mh2 + (power2(sa)*(mst24*(mt2*(72*logmh2Q2
     *(-3 + logmst22Q2) - 36*logmst22Q2 + 36*power2(logmh2Q2) - 24*power2(
     logmst22Q2) + 7*(42 + power2(Pi)))*(-1 + power2(sa)) + mu2*(84 - 24*power2
     (logmst22Q2)*(-1 + power2(sa)) + 6*(49 - 36*logmh2Q2 + 6*power2(logmh2Q2))
     *power2(sa) + power2(Pi)*(2 + 7*power2(sa)) + 36*logmst22Q2*(-2 + (-1 + 2*
     logmh2Q2)*power2(sa)))*(-1 + power2(snt))*power2(snt)) - 2*invdmhH*(42 -
     36*logmst22Q2 + 12*power2(logmst22Q2) + power2(Pi))*(-1 + power2(sa))*(mt2
      + mu2*(-1 + power2(snt))*power2(snt))*(-(mH2*mst24) + 4*power6(mst2))))/(
     2.*mh2)) + DeltaInv(mst22,mt2,mu2)*(-(Al*At*invdmst*invdmstau*(
     invdmstau1mu - invdmstau2mu)*mu2*(mu2*(3*mst22*mu2*(4*logmst22Q2*logmt2Q2
     - 2*(-6 + logmu2Q2)*logmu2Q2 - 2*logmst22Q2*(3 + logmu2Q2) - 2*logmt2Q2*(3
      + logmu2Q2) + power2(logmst22Q2) + power2(logmt2Q2)) - (mt2 - mu2)*mu2*(
     42 + 6*logmt2Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*power2(logmt2Q2) + 3*
     power2(logmu2Q2) + power2(Pi)) - mst22*mt2*(294 + 36*logmt2Q2*(-5 +
     logmu2Q2) - 54*logmu2Q2 + 6*logmst22Q2*(4*logmt2Q2 - 3*(1 + logmu2Q2)) + 3
     *power2(logmst22Q2) + 30*power2(logmt2Q2) + 9*power2(logmu2Q2) + 7*power2(
     Pi))) - mst24*(3*mu2*(42 + 4*logmst22Q2*(-3 + 2*logmt2Q2 - logmu2Q2) + 2*
     logmt2Q2*(-15 + logmu2Q2) + 6*logmu2Q2 + 2*power2(logmst22Q2) + 5*power2(
     logmt2Q2) - power2(logmu2Q2) + power2(Pi)) + mt2*(6*logmst22Q2*(-9 + 4*
     logmt2Q2 - logmu2Q2) + 6*logmt2Q2*(-15 + logmu2Q2) + 9*power2(logmst22Q2)
     + 15*power2(logmt2Q2) + 4*(42 + power2(Pi)))) + (6*logmst22Q2*(-3 + 2*
     logmt2Q2 - logmu2Q2) + 6*logmt2Q2*(-9 + logmu2Q2) + 3*power2(logmst22Q2) +
     9*power2(logmt2Q2) + 2*(42 + power2(Pi)))*power6(mst2))) + (At*invdmst*
     invdmstau*(invdmstau1mu - invdmstau2mu)*mu*mu2*sb*(mu2*(3*mst22*mu2*(4*
     logmst22Q2*logmt2Q2 - 2*(-6 + logmu2Q2)*logmu2Q2 - 2*logmst22Q2*(3 +
     logmu2Q2) - 2*logmt2Q2*(3 + logmu2Q2) + power2(logmst22Q2) + power2(
     logmt2Q2)) - (mt2 - mu2)*mu2*(42 + 6*logmt2Q2*(-3 + logmu2Q2) - 18*
     logmu2Q2 + 3*power2(logmt2Q2) + 3*power2(logmu2Q2) + power2(Pi)) - mst22*
     mt2*(294 + 36*logmt2Q2*(-5 + logmu2Q2) - 54*logmu2Q2 + 6*logmst22Q2*(4*
     logmt2Q2 - 3*(1 + logmu2Q2)) + 3*power2(logmst22Q2) + 30*power2(logmt2Q2)
     + 9*power2(logmu2Q2) + 7*power2(Pi))) - mst24*(3*mu2*(42 + 4*logmst22Q2*(-
     3 + 2*logmt2Q2 - logmu2Q2) + 2*logmt2Q2*(-15 + logmu2Q2) + 6*logmu2Q2 + 2*
     power2(logmst22Q2) + 5*power2(logmt2Q2) - power2(logmu2Q2) + power2(Pi)) +
     mt2*(6*logmst22Q2*(-9 + 4*logmt2Q2 - logmu2Q2) + 6*logmt2Q2*(-15 +
     logmu2Q2) + 9*power2(logmst22Q2) + 15*power2(logmt2Q2) + 4*(42 + power2(Pi
     )))) + (6*logmst22Q2*(-3 + 2*logmt2Q2 - logmu2Q2) + 6*logmt2Q2*(-9 +
     logmu2Q2) + 3*power2(logmst22Q2) + 9*power2(logmt2Q2) + 2*(42 + power2(Pi)
     ))*power6(mst2)))/cb + (mu2*(power2(invdmstau1mu)*(mst24*mt2*(42 + 6*
     logmst22Q2*(-3 + logmt2Q2) - 18*logmt2Q2 + 3*power2(logmst22Q2) + 3*power2
     (logmt2Q2) + power2(Pi)) + mst24*mu2*(42 + 12*logmst22Q2*(-3 + logmt2Q2) +
     18*logmu2Q2 - 6*logmt2Q2*(3 + logmu2Q2) + 6*power2(logmst22Q2) + 3*power2(
     logmt2Q2) - 3*power2(logmu2Q2) + power2(Pi)) + (mt2 - mu2)*power2(mu2)*(42
      + 6*logmt2Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*power2(logmt2Q2) + 3*
     power2(logmu2Q2) + power2(Pi)) + mst22*mu2*(mu2*(42 - 6*logmst22Q2*(-3 +
     logmt2Q2) - 36*logmu2Q2 + 6*logmt2Q2*(-3 + 2*logmu2Q2) - 3*power2(
     logmst22Q2) + 3*power2(logmt2Q2) + 6*power2(logmu2Q2) + power2(Pi)) + mt2*
     (168 + 6*logmst22Q2*(-3 + 3*logmt2Q2 - 2*logmu2Q2) + 18*logmt2Q2*(-6 +
     logmu2Q2) - 18*logmu2Q2 + 3*power2(logmst22Q2) + 18*power2(logmt2Q2) + 3*
     power2(logmu2Q2) + 4*power2(Pi))) - (42 + 6*logmst22Q2*(-3 + logmt2Q2) -
     18*logmt2Q2 + 3*power2(logmst22Q2) + 3*power2(logmt2Q2) + power2(Pi))*
     power6(mst2)) - 8*invdmst*invdmstau*invdmstau1mu*mu2*(mu2*(3*mst22*mu2*(4*
     logmst22Q2*logmt2Q2 - 2*(-6 + logmu2Q2)*logmu2Q2 - 2*logmst22Q2*(3 +
     logmu2Q2) - 2*logmt2Q2*(3 + logmu2Q2) + power2(logmst22Q2) + power2(
     logmt2Q2)) - (mt2 - mu2)*mu2*(42 + 6*logmt2Q2*(-3 + logmu2Q2) - 18*
     logmu2Q2 + 3*power2(logmt2Q2) + 3*power2(logmu2Q2) + power2(Pi)) - mst22*
     mt2*(294 + 36*logmt2Q2*(-5 + logmu2Q2) - 54*logmu2Q2 + 6*logmst22Q2*(4*
     logmt2Q2 - 3*(1 + logmu2Q2)) + 3*power2(logmst22Q2) + 30*power2(logmt2Q2)
     + 9*power2(logmu2Q2) + 7*power2(Pi))) - mst24*(3*mu2*(42 + 4*logmst22Q2*(-
     3 + 2*logmt2Q2 - logmu2Q2) + 2*logmt2Q2*(-15 + logmu2Q2) + 6*logmu2Q2 + 2*
     power2(logmst22Q2) + 5*power2(logmt2Q2) - power2(logmu2Q2) + power2(Pi)) +
     mt2*(6*logmst22Q2*(-9 + 4*logmt2Q2 - logmu2Q2) + 6*logmt2Q2*(-15 +
     logmu2Q2) + 9*power2(logmst22Q2) + 15*power2(logmt2Q2) + 4*(42 + power2(Pi
     )))) + (6*logmst22Q2*(-3 + 2*logmt2Q2 - logmu2Q2) + 6*logmt2Q2*(-9 +
     logmu2Q2) + 3*power2(logmst22Q2) + 9*power2(logmt2Q2) + 2*(42 + power2(Pi)
     ))*power6(mst2)) + invdmstau2mu*(invdmstau2mu*(mst24*mt2*(42 + 6*
     logmst22Q2*(-3 + logmt2Q2) - 18*logmt2Q2 + 3*power2(logmst22Q2) + 3*power2
     (logmt2Q2) + power2(Pi)) + mst24*mu2*(42 + 12*logmst22Q2*(-3 + logmt2Q2) +
     18*logmu2Q2 - 6*logmt2Q2*(3 + logmu2Q2) + 6*power2(logmst22Q2) + 3*power2(
     logmt2Q2) - 3*power2(logmu2Q2) + power2(Pi)) + (mt2 - mu2)*power2(mu2)*(42
      + 6*logmt2Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*power2(logmt2Q2) + 3*
     power2(logmu2Q2) + power2(Pi)) + mst22*mu2*(mu2*(42 - 6*logmst22Q2*(-3 +
     logmt2Q2) - 36*logmu2Q2 + 6*logmt2Q2*(-3 + 2*logmu2Q2) - 3*power2(
     logmst22Q2) + 3*power2(logmt2Q2) + 6*power2(logmu2Q2) + power2(Pi)) + mt2*
     (168 + 6*logmst22Q2*(-3 + 3*logmt2Q2 - 2*logmu2Q2) + 18*logmt2Q2*(-6 +
     logmu2Q2) - 18*logmu2Q2 + 3*power2(logmst22Q2) + 18*power2(logmt2Q2) + 3*
     power2(logmu2Q2) + 4*power2(Pi))) - (42 + 6*logmst22Q2*(-3 + logmt2Q2) -
     18*logmt2Q2 + 3*power2(logmst22Q2) + 3*power2(logmt2Q2) + power2(Pi))*
     power6(mst2)) + 8*invdmst*invdmstau*mu2*(mu2*(3*mst22*mu2*(4*logmst22Q2*
     logmt2Q2 - 2*(-6 + logmu2Q2)*logmu2Q2 - 2*logmst22Q2*(3 + logmu2Q2) - 2*
     logmt2Q2*(3 + logmu2Q2) + power2(logmst22Q2) + power2(logmt2Q2)) - (mt2 -
     mu2)*mu2*(42 + 6*logmt2Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*power2(
     logmt2Q2) + 3*power2(logmu2Q2) + power2(Pi)) - mst22*mt2*(294 + 36*
     logmt2Q2*(-5 + logmu2Q2) - 54*logmu2Q2 + 6*logmst22Q2*(4*logmt2Q2 - 3*(1 +
     logmu2Q2)) + 3*power2(logmst22Q2) + 30*power2(logmt2Q2) + 9*power2(
     logmu2Q2) + 7*power2(Pi))) - mst24*(3*mu2*(42 + 4*logmst22Q2*(-3 + 2*
     logmt2Q2 - logmu2Q2) + 2*logmt2Q2*(-15 + logmu2Q2) + 6*logmu2Q2 + 2*power2
     (logmst22Q2) + 5*power2(logmt2Q2) - power2(logmu2Q2) + power2(Pi)) + mt2*(
     6*logmst22Q2*(-9 + 4*logmt2Q2 - logmu2Q2) + 6*logmt2Q2*(-15 + logmu2Q2) +
     9*power2(logmst22Q2) + 15*power2(logmt2Q2) + 4*(42 + power2(Pi)))) + (6*
     logmst22Q2*(-3 + 2*logmt2Q2 - logmu2Q2) + 6*logmt2Q2*(-9 + logmu2Q2) + 3*
     power2(logmst22Q2) + 9*power2(logmt2Q2) + 2*(42 + power2(Pi)))*power6(mst2
     )))))/8.);

   return result * power2(ytau) * power2(yt) * twoLoop;
}

double delta_mtau_2loop_atau_ab(const Parameters& pars)
{
   const Real ytau  = pars.ytau;
   const Real yb    = pars.yb;
   const Real xt    = pars.xt;
   const Real xb    = pars.xb;
   const Real mt    = pars.mt;
   const Real mt2   = power2(pars.mt);
   const Real mst1  = pars.mst1;
   const Real mst12 = power2(pars.mst1);
   const Real mst14 = power4(pars.mst1);
   const Real mst2  = pars.mst2;
   const Real mst22 = power2(pars.mst2);
   const Real mst24 = power4(pars.mst2);
   const Real msb1  = pars.msb1;
   const Real msb12 = power2(pars.msb1);
   const Real msb14 = power4(pars.msb1);
   const Real msb2  = pars.msb2;
   const Real msb22 = power2(pars.msb2);
   const Real msb24 = power4(pars.msb2);
   const Real mstau12 = power2(pars.mstau1);
   const Real mstau22 = power2(pars.mstau2);
   const Real msntau2 = power2(pars.msntau);
   const Real msntau4 = power4(pars.msntau);
   const Real mh2   = power2(pars.mh);
   const Real mH2   = power2(pars.mH);
   const Real mC2   = power2(pars.mC);
   const Real mA2   = power2(pars.mA);
   const Real mu    = pars.mu;
   const Real mu2   = power2(pars.mu);
   const Real tb    = pars.tb;
   const Real sb    = tb / std::sqrt(1 + power2(tb));
   const Real cb    = 1  / std::sqrt(1 + power2(tb));
   const Real Q2    = power2(pars.Q);
   const Real snt   = calc_sin_theta(mt, xt, mst12, mst22);
   const Real alpha = calc_alpha(mh2, mH2, tb);
   const Real sa    = std::sin(alpha);
   const Real ca    = std::cos(alpha);
   const Real At    = xt + mu/tb;
   const Real Ab    = xb + mu*tb;

   const Real invdmst       = 1/(-mst12 + mst22);
   const Real invdmsb       = 1/(msb12 - msb22);
   const Real invdmsb1stau1 = 1/(msb12 - mstau12);
   const Real invdmsb1stau2 = 1/(msb12 - mstau22);
   const Real invdmsb2stau2 = 1/(-msb22 + mstau22);
   const Real invdmsb2stau1 = 1/(-msb22 + mstau12);
   const Real invdmsntaust1 = 1/(-msntau2 + mst12);
   const Real invdmsntaust2 = 1/(-msntau2 + mst22);
   const Real invdmstau1mu  = 1/(-mstau12 + mu2);
   const Real invdmstau2mu  = 1/(-mstau22 + mu2);
   const Real invdmsntau2mu = 1/(-msntau2 + mu2);
   const Real invdmhH       = 1/(-mh2 + mH2);

   const Real logmstau12Q2  = std::log(mstau12/Q2);
   const Real logmstau22Q2  = std::log(mstau22/Q2);
   const Real logmH2Q2      = std::log(mH2/Q2);
   const Real logmA2Q2      = std::log(mA2/Q2);
   const Real logmh2Q2      = std::log(mh2/Q2);
   const Real logmu2Q2      = std::log(mu2/Q2);
   const Real logmC2Q2      = std::log(mC2/Q2);
   const Real logmsntau2Q2  = std::log(msntau2/Q2);
   const Real logmst12Q2    = std::log(mst12/Q2);
   const Real logmst22Q2    = std::log(mst22/Q2);
   const Real logmsb12Q2    = std::log(msb12/Q2);
   const Real logmsb22Q2    = std::log(msb22/Q2);
   const Real logmt2Q2      = std::log(mt2/Q2);

   const double result =
   Fin3(msb22,mt2,mu2,Q2)*((-3*(msb22 - mt2 - mu2)*(-1 + 2*invdmsntau2mu*mu2)*
     power2(invdmsntau2mu))/8. + (3*mu2*(-msb24 + (mt2 - mu2)*mu2 + msb22*(mt2
     + 2*mu2))*DeltaInv(msb22,mt2,mu2)*power2(invdmsntau2mu))/4.) + Fin3(mA2,
     msb22,msb12,Q2)*((3*(mA2 - msb12 - msb22)*DeltaInv(mA2,msb22,msb12)*power2
     (Ab))/(4.*mA2) - (3*power2(Ab))/(4.*power2(mA2))) - (3*(mC2 + mt2)*Fin20(
     mC2,mt2,Q2))/(4.*power2(mC2)) + Fin3(msb22,msntau2,mt2,Q2)*((3*
     invdmsntau2mu*(1 + invdmsntau2mu*(msb22 - mt2 - mu2))*(-1 + 2*
     invdmsntau2mu*mu2))/8. + (3*DeltaInv(msb22,msntau2,mt2)*(-msntau2 + mt2 -
     2*mu2 + msb22*(-1 + invdmsntau2mu*mu2)*(-2 + invdmsntau2mu*(mt2 + 2*mu2))
     - mu2*power2(invdmsntau2mu)*(msb24 - mt2*mu2 + power2(mu2)) +
     invdmsntau2mu*(msb24 - 2*mt2*mu2 + 3*power2(mu2))))/4.) + (DeltaInv(mA2,
     msb22,msb12)*power2(Ab)*(3*(logmsb12Q2 - logmsb22Q2)*(-6 + 2*logmA2Q2 +
     logmsb12Q2 + logmsb22Q2)*msb24 + power2(mA2)*(42 + 6*logmA2Q2*(-3 +
     logmsb12Q2) - 18*logmsb12Q2 + 3*power2(logmA2Q2) + 3*power2(logmsb12Q2) +
     power2(Pi)) - msb12*msb22*(-6*(3 + logmA2Q2)*logmsb22Q2 + 6*logmsb12Q2*(-9
      + logmA2Q2 + 2*logmsb22Q2) + 9*power2(logmsb12Q2) + 3*power2(logmsb22Q2)
     + 2*(42 + power2(Pi))) - mA2*(msb12*(42 + 6*logmA2Q2*(-3 + logmsb12Q2) -
     18*logmsb12Q2 + 3*power2(logmA2Q2) + 3*power2(logmsb12Q2) + power2(Pi)) +
     msb22*(42 - 36*logmsb12Q2 + 6*logmA2Q2*(-3 + 2*logmsb12Q2 - logmsb22Q2) +
     18*logmsb22Q2 + 3*power2(logmA2Q2) + 6*power2(logmsb12Q2) - 3*power2(
     logmsb22Q2) + power2(Pi)))))/(8.*mA2) + (power2(Ab)*((-2*(36 + 6*logmA2Q2*
     (-2 + logmsb12Q2) + 6*invdmsb*logmsb22Q2*msb22 - 6*logmsb12Q2*(4 + invdmsb
     *msb22) + 3*power2(logmA2Q2) + 3*power2(logmsb12Q2) + power2(Pi)))/mA2 + (
     2*(36 + 6*logmC2Q2*(-2 + logmsb22Q2) - 18*logmsb22Q2 + 6*invdmsb*
     logmsb22Q2*msb22 - 6*logmsb12Q2*(1 + invdmsb*msb22) + 3*power2(logmC2Q2) +
     3*power2(logmsb22Q2) + power2(Pi)))/mC2 - ((-1 + invdmhH*mH2)*(3*((-7 + 4*
     logmsb12Q2)*msb12 + (-7 + 4*logmsb22Q2)*msb22) + 3*mH2*(5 + 12*invdmsb*
     logmsb22Q2*msb22 - 12*logmsb12Q2*(1 + invdmsb*msb22) + logmH2Q2*(-4 - 8*
     invdmsb*logmsb22Q2*msb22 + 8*logmsb12Q2*(1 + invdmsb*msb22))) + 2*invdmhH*
     power2(mH2)*(21 - 18*logmH2Q2 + 6*power2(logmH2Q2) + power2(Pi)))*(-1 +
     power2(sa))*power2(sa))/power2(mh2) + (power2(sa)*(81 + 21*invdmhH*msb12 +
     21*invdmhH*msb22 - 12*invdmhH*logmsb22Q2*msb22 + 60*invdmsb*logmsb22Q2*
     msb22 - 96*invdmhH*invdmsb*logmsb22Q2*mH2*msb22 - 138*power2(invdmhH)*
     power2(mH2) + 8*invdmhH*mH2*power2(Pi) - 10*power2(invdmhH)*power2(mH2)*
     power2(Pi) + 12*invdmhH*mH2*(-4 + 5*invdmhH*mH2)*power2(logmH2Q2)*(-1 +
     power2(sa)) + 12*logmH2Q2*(7 + 6*invdmsb*logmsb22Q2*msb22 - 8*invdmhH*
     invdmsb*logmsb22Q2*mH2*msb22 + 2*logmsb12Q2*(-3 + 4*invdmhH*mH2)*(1 +
     invdmsb*msb22) - 11*power2(invdmhH)*power2(mH2))*(-1 + power2(sa)) + 135*
     power2(sa) - 72*logmh2Q2*power2(sa) - 21*invdmhH*msb12*power2(sa) - 21*
     invdmhH*msb22*power2(sa) + 12*invdmhH*logmsb22Q2*msb22*power2(sa) - 24*
     invdmsb*logmsb22Q2*msb22*power2(sa) + 96*invdmhH*invdmsb*logmsb22Q2*mH2*
     msb22*power2(sa) + 18*power2(logmh2Q2)*power2(sa) + 18*power2(logmsb12Q2)*
     power2(sa) + 138*power2(invdmhH)*power2(mH2)*power2(sa) + 6*power2(Pi)*
     power2(sa) - 8*invdmhH*mH2*power2(Pi)*power2(sa) + 10*power2(invdmhH)*
     power2(mH2)*power2(Pi)*power2(sa) - 12*logmsb12Q2*(5 + invdmsb*msb22*(5 -
     2*power2(sa)) + invdmhH*(-msb12 + 8*mH2*(1 + invdmsb*msb22))*(-1 + power2(
     sa)) + 7*power2(sa) - 3*logmh2Q2*power2(sa))))/mh2 - ((-1 + power2(sa))*(
     216 + 36*invdmsb*logmsb22Q2*msb22 + 6*power2(Pi) - 18*power2(logmsb12Q2)*(
     -1 + power2(sa)) - 216*power2(sa) + 81*invdmhH*mH2*power2(sa) - 144*
     invdmhH*logmh2Q2*mH2*power2(sa) - 36*invdmsb*logmsb22Q2*msb22*power2(sa) +
     60*invdmhH*invdmsb*logmsb22Q2*mH2*msb22*power2(sa) - 72*invdmhH*invdmsb*
     logmh2Q2*logmsb22Q2*mH2*msb22*power2(sa) + 36*invdmhH*mH2*power2(logmh2Q2)
     *power2(sa) + 96*power2(invdmhH)*power2(mH2)*power2(sa) - 6*power2(Pi)*
     power2(sa) + 8*power2(invdmhH)*power2(mH2)*power2(Pi)*power2(sa) + 6*
     power2(logmH2Q2)*(3 + (-3 - 6*invdmhH*mH2 + 8*power2(invdmhH)*power2(mH2))
     *power2(sa)) - 12*logmH2Q2*(6 + 3*logmsb12Q2*(-1 + power2(sa)) + (-6 - 5*
     invdmhH*mH2 + 8*power2(invdmhH)*power2(mH2))*power2(sa)) + 12*logmsb12Q2*(
     -12 + (12 + invdmhH*(-5 + 6*logmh2Q2)*mH2)*power2(sa) + invdmsb*msb22*(-3
     + (3 + invdmhH*(-5 + 6*logmh2Q2)*mH2)*power2(sa)))))/mH2))/16. + (Ab*mu*(
     12*(-1 + logmsb12Q2 + invdmsb*logmsb12Q2*msb22 - invdmsb*logmsb22Q2*msb22)
     *sb*(mA2*mH2*power2(mh2) - mC2*mH2*power2(mh2) + mA2*mC2*(mh2*mH2*(5 + 8*
     invdmhH*mH2*(-1 + power2(sa)) - 2*logmH2Q2*(-3 + 4*invdmhH*mH2)*(-1 +
     power2(sa)) - 2*power2(sa))*power2(sa) + (-3 + 2*logmH2Q2)*(-1 + invdmhH*
     mH2)*power2(mH2)*(-1 + power2(sa))*power2(sa) + power2(mh2)*(-1 + power2(
     sa))*(-3 + (3 + invdmhH*(-5 + 6*logmh2Q2)*mH2)*power2(sa)))) + ca*cb*mA2*
     mC2*sa*(-(mH2*(-1 + invdmhH*mH2)*(3*(-7 + 4*logmH2Q2)*mH2 + 3*(-7 + 4*
     logmsb12Q2)*msb12 + 3*(-7 + 4*logmsb22Q2)*msb22 + 2*invdmhH*power2(mH2)*(
     21 - 18*logmH2Q2 + 6*power2(logmH2Q2) + power2(Pi)))*(-1 + 2*power2(sa)))
     + mh2*mH2*(21 + 6*(53 + 12*logmh2Q2*(-2 + logmsb12Q2) - 36*logmsb12Q2 + 6*
     power2(logmh2Q2) + 6*power2(logmsb12Q2) + 2*power2(Pi))*power2(sa) + 12*
     invdmhH*mH2*(-4 + 5*invdmhH*mH2)*power2(logmH2Q2)*(-1 + 2*power2(sa)) - 12
     *logmH2Q2*(-1 - 8*invdmhH*mH2 + 11*power2(invdmhH)*power2(mH2))*(-1 + 2*
     power2(sa)) + 2*power2(invdmhH)*power2(mH2)*(69 + 5*power2(Pi))*(-1 + 2*
     power2(sa)) - invdmhH*(-3*(-7*msb12 + 4*logmsb12Q2*msb12 - 7*msb22 + 4*
     logmsb22Q2*msb22) + 8*mH2*(12 + power2(Pi)))*(-1 + 2*power2(sa))) + power2
     (mh2)*(-360 + 21*invdmhH*mH2 - 72*invdmhH*logmh2Q2*mH2 + 36*invdmhH*mH2*
     power2(logmh2Q2) + 96*power2(invdmhH)*power2(mH2) - 12*power2(Pi) + 8*
     power2(invdmhH)*power2(mH2)*power2(Pi) - 216*logmsb12Q2*(-1 + power2(sa))
     + 36*power2(logmsb12Q2)*(-1 + power2(sa)) + 360*power2(sa) - 42*invdmhH*
     mH2*power2(sa) + 144*invdmhH*logmh2Q2*mH2*power2(sa) - 72*invdmhH*mH2*
     power2(logmh2Q2)*power2(sa) - 192*power2(invdmhH)*power2(mH2)*power2(sa) +
     12*power2(Pi)*power2(sa) - 16*power2(invdmhH)*power2(mH2)*power2(Pi)*
     power2(sa) - 12*power2(logmH2Q2)*(3 + invdmhH*mH2*(3 - 6*power2(sa)) - 3*
     power2(sa) + 4*power2(invdmhH)*power2(mH2)*(-1 + 2*power2(sa))) + 12*
     logmH2Q2*(5*invdmhH*mH2*(1 - 2*power2(sa)) - 12*(-1 + power2(sa)) + 6*
     logmsb12Q2*(-1 + power2(sa)) + 8*power2(invdmhH)*power2(mH2)*(-1 + 2*
     power2(sa)))))))/(16.*cb*mA2*mC2*mH2*power2(mh2)) - (3*(mst12 - 2*
     invdmsntau2mu*mst12*mu2 + mu2*(-3 + 2*invdmsntau2mu*mu2))*Fin20(mst12,mu2,
     Q2)*power2(invdmsntau2mu)*(-1 + power2(snt)))/8. + (3*(mst22 - 2*
     invdmsntau2mu*mst22*mu2 + mu2*(-3 + 2*invdmsntau2mu*mu2))*Fin20(mst22,mu2,
     Q2)*power2(invdmsntau2mu)*power2(snt))/8. + Fin20(msb22,mst12,Q2)*((3*Ab*
     At*invdmst*mt2)/(2.*power2(mC2)) + (3*power2(Ab)*(-1 + power2(snt)))/(4.*
     power2(mC2)) - (3*mt2*power2(snt))/(4.*power2(mC2))) + Fin20(msb22,mst22,
     Q2)*((-3*Ab*At*invdmst*mt2)/(2.*power2(mC2)) + (3*mt2*(-1 + power2(snt)))/
     (4.*power2(mC2)) - (3*power2(Ab)*power2(snt))/(4.*power2(mC2))) + DeltaInv
     (msb22,mst12,mC2)*((Ab*At*invdmst*mt2*(3*(logmsb22Q2 - logmst12Q2)*(-6 + 2
     *logmC2Q2 + logmsb22Q2 + logmst12Q2)*mst14 + power2(mC2)*(42 + 6*logmC2Q2*
     (-3 + logmsb22Q2) - 18*logmsb22Q2 + 3*power2(logmC2Q2) + 3*power2(
     logmsb22Q2) + power2(Pi)) - msb22*mst12*(-6*(3 + logmC2Q2)*logmst12Q2 + 6*
     logmsb22Q2*(-9 + logmC2Q2 + 2*logmst12Q2) + 9*power2(logmsb22Q2) + 3*
     power2(logmst12Q2) + 2*(42 + power2(Pi))) - mC2*(msb22*(42 + 6*logmC2Q2*(-
     3 + logmsb22Q2) - 18*logmsb22Q2 + 3*power2(logmC2Q2) + 3*power2(logmsb22Q2
     ) + power2(Pi)) + mst12*(42 - 36*logmsb22Q2 + 6*logmC2Q2*(-3 + 2*
     logmsb22Q2 - logmst12Q2) + 18*logmst12Q2 + 3*power2(logmC2Q2) + 6*power2(
     logmsb22Q2) - 3*power2(logmst12Q2) + power2(Pi)))))/(4.*mC2) + (power2(Ab)
     *(3*(logmsb22Q2 - logmst12Q2)*(-6 + 2*logmC2Q2 + logmsb22Q2 + logmst12Q2)*
     mst14 + power2(mC2)*(42 + 6*logmC2Q2*(-3 + logmsb22Q2) - 18*logmsb22Q2 + 3
     *power2(logmC2Q2) + 3*power2(logmsb22Q2) + power2(Pi)) - msb22*mst12*(-6*(
     3 + logmC2Q2)*logmst12Q2 + 6*logmsb22Q2*(-9 + logmC2Q2 + 2*logmst12Q2) + 9
     *power2(logmsb22Q2) + 3*power2(logmst12Q2) + 2*(42 + power2(Pi))) - mC2*(
     msb22*(42 + 6*logmC2Q2*(-3 + logmsb22Q2) - 18*logmsb22Q2 + 3*power2(
     logmC2Q2) + 3*power2(logmsb22Q2) + power2(Pi)) + mst12*(42 - 36*logmsb22Q2
      + 6*logmC2Q2*(-3 + 2*logmsb22Q2 - logmst12Q2) + 18*logmst12Q2 + 3*power2(
     logmC2Q2) + 6*power2(logmsb22Q2) - 3*power2(logmst12Q2) + power2(Pi))))*(-
     1 + power2(snt)))/(8.*mC2) + (mt2*(-3*(logmsb22Q2 - logmst12Q2)*(-6 + 2*
     logmC2Q2 + logmsb22Q2 + logmst12Q2)*mst14 - power2(mC2)*(42 + 6*logmC2Q2*(
     -3 + logmsb22Q2) - 18*logmsb22Q2 + 3*power2(logmC2Q2) + 3*power2(
     logmsb22Q2) + power2(Pi)) + msb22*mst12*(-6*(3 + logmC2Q2)*logmst12Q2 + 6*
     logmsb22Q2*(-9 + logmC2Q2 + 2*logmst12Q2) + 9*power2(logmsb22Q2) + 3*
     power2(logmst12Q2) + 2*(42 + power2(Pi))) + mC2*(msb22*(42 + 6*logmC2Q2*(-
     3 + logmsb22Q2) - 18*logmsb22Q2 + 3*power2(logmC2Q2) + 3*power2(logmsb22Q2
     ) + power2(Pi)) + mst12*(42 - 36*logmsb22Q2 + 6*logmC2Q2*(-3 + 2*
     logmsb22Q2 - logmst12Q2) + 18*logmst12Q2 + 3*power2(logmC2Q2) + 6*power2(
     logmsb22Q2) - 3*power2(logmst12Q2) + power2(Pi))))*power2(snt))/(8.*mC2))
     + DeltaInv(msb22,mC2,mst22)*((Ab*At*invdmst*mt2*(-3*(logmsb22Q2 -
     logmst22Q2)*(-6 + 2*logmC2Q2 + logmsb22Q2 + logmst22Q2)*mst24 - power2(mC2
     )*(42 + 6*logmC2Q2*(-3 + logmsb22Q2) - 18*logmsb22Q2 + 3*power2(logmC2Q2)
     + 3*power2(logmsb22Q2) + power2(Pi)) + msb22*mst22*(-6*(3 + logmC2Q2)*
     logmst22Q2 + 6*logmsb22Q2*(-9 + logmC2Q2 + 2*logmst22Q2) + 9*power2(
     logmsb22Q2) + 3*power2(logmst22Q2) + 2*(42 + power2(Pi))) + mC2*(msb22*(42
      + 6*logmC2Q2*(-3 + logmsb22Q2) - 18*logmsb22Q2 + 3*power2(logmC2Q2) + 3*
     power2(logmsb22Q2) + power2(Pi)) + mst22*(42 - 36*logmsb22Q2 + 6*logmC2Q2*
     (-3 + 2*logmsb22Q2 - logmst22Q2) + 18*logmst22Q2 + 3*power2(logmC2Q2) + 6*
     power2(logmsb22Q2) - 3*power2(logmst22Q2) + power2(Pi)))))/(4.*mC2) + (mt2
     *(3*(logmsb22Q2 - logmst22Q2)*(-6 + 2*logmC2Q2 + logmsb22Q2 + logmst22Q2)*
     mst24 + power2(mC2)*(42 + 6*logmC2Q2*(-3 + logmsb22Q2) - 18*logmsb22Q2 + 3
     *power2(logmC2Q2) + 3*power2(logmsb22Q2) + power2(Pi)) - msb22*mst22*(-6*(
     3 + logmC2Q2)*logmst22Q2 + 6*logmsb22Q2*(-9 + logmC2Q2 + 2*logmst22Q2) + 9
     *power2(logmsb22Q2) + 3*power2(logmst22Q2) + 2*(42 + power2(Pi))) - mC2*(
     msb22*(42 + 6*logmC2Q2*(-3 + logmsb22Q2) - 18*logmsb22Q2 + 3*power2(
     logmC2Q2) + 3*power2(logmsb22Q2) + power2(Pi)) + mst22*(42 - 36*logmsb22Q2
      + 6*logmC2Q2*(-3 + 2*logmsb22Q2 - logmst22Q2) + 18*logmst22Q2 + 3*power2(
     logmC2Q2) + 6*power2(logmsb22Q2) - 3*power2(logmst22Q2) + power2(Pi))))*(-
     1 + power2(snt)))/(8.*mC2) + (power2(Ab)*(-3*(logmsb22Q2 - logmst22Q2)*(-6
      + 2*logmC2Q2 + logmsb22Q2 + logmst22Q2)*mst24 - power2(mC2)*(42 + 6*
     logmC2Q2*(-3 + logmsb22Q2) - 18*logmsb22Q2 + 3*power2(logmC2Q2) + 3*power2
     (logmsb22Q2) + power2(Pi)) + msb22*mst22*(-6*(3 + logmC2Q2)*logmst22Q2 + 6
     *logmsb22Q2*(-9 + logmC2Q2 + 2*logmst22Q2) + 9*power2(logmsb22Q2) + 3*
     power2(logmst22Q2) + 2*(42 + power2(Pi))) + mC2*(msb22*(42 + 6*logmC2Q2*(-
     3 + logmsb22Q2) - 18*logmsb22Q2 + 3*power2(logmC2Q2) + 3*power2(logmsb22Q2
     ) + power2(Pi)) + mst22*(42 - 36*logmsb22Q2 + 6*logmC2Q2*(-3 + 2*
     logmsb22Q2 - logmst22Q2) + 18*logmst22Q2 + 3*power2(logmC2Q2) + 6*power2(
     logmsb22Q2) - 3*power2(logmst22Q2) + power2(Pi))))*power2(snt))/(8.*mC2))
     + Fin3(msb22,mst12,mC2,Q2)*((-3*Ab*At*invdmst*mt2)/(2.*power2(mC2)) - (3*
     power2(Ab)*(-1 + power2(snt)))/(4.*power2(mC2)) + (3*mt2*power2(snt))/(4.*
     power2(mC2)) + DeltaInv(msb22,mst12,mC2)*((3*Ab*At*invdmst*(mC2 - msb22 -
     mst12)*mt2)/(2.*mC2) + (3*(mC2 - msb22 - mst12)*power2(Ab)*(-1 + power2(
     snt)))/(4.*mC2) + (3*(-mC2 + msb22 + mst12)*mt2*power2(snt))/(4.*mC2))) +
     Fin3(msb22,mC2,mst22,Q2)*((3*Ab*At*invdmst*mt2)/(2.*power2(mC2)) - (3*mt2*
     (-1 + power2(snt)))/(4.*power2(mC2)) + (3*power2(Ab)*power2(snt))/(4.*
     power2(mC2)) + DeltaInv(msb22,mC2,mst22)*((3*Ab*At*invdmst*(-mC2 + msb22 +
     mst22)*mt2)/(2.*mC2) + (3*(mC2 - msb22 - mst22)*mt2*(-1 + power2(snt)))/(
     4.*mC2) + (3*(-mC2 + msb22 + mst22)*power2(Ab)*power2(snt))/(4.*mC2))) +
     DeltaInv(msb22,mH2,msb12)*((-3*Ab*ca*mu*sa*(3*(logmsb12Q2 - logmsb22Q2)*(-
     6 + 2*logmH2Q2 + logmsb12Q2 + logmsb22Q2)*msb24 + power2(mH2)*(42 + 6*
     logmH2Q2*(-3 + logmsb12Q2) - 18*logmsb12Q2 + 3*power2(logmH2Q2) + 3*power2
     (logmsb12Q2) + power2(Pi)) - msb12*msb22*(-6*(3 + logmH2Q2)*logmsb22Q2 + 6
     *logmsb12Q2*(-9 + logmH2Q2 + 2*logmsb22Q2) + 9*power2(logmsb12Q2) + 3*
     power2(logmsb22Q2) + 2*(42 + power2(Pi))) - mH2*(msb12*(42 + 6*logmH2Q2*(-
     3 + logmsb12Q2) - 18*logmsb12Q2 + 3*power2(logmH2Q2) + 3*power2(logmsb12Q2
     ) + power2(Pi)) + msb22*(42 - 36*logmsb12Q2 + 6*logmH2Q2*(-3 + 2*
     logmsb12Q2 - logmsb22Q2) + 18*logmsb22Q2 + 3*power2(logmH2Q2) + 6*power2(
     logmsb12Q2) - 3*power2(logmsb22Q2) + power2(Pi))))*(-1 + power2(sa)))/(4.*
     mH2) + (3*mu2*(3*(logmsb12Q2 - logmsb22Q2)*(-6 + 2*logmH2Q2 + logmsb12Q2 +
     logmsb22Q2)*msb24 + power2(mH2)*(42 + 6*logmH2Q2*(-3 + logmsb12Q2) - 18*
     logmsb12Q2 + 3*power2(logmH2Q2) + 3*power2(logmsb12Q2) + power2(Pi)) -
     msb12*msb22*(-6*(3 + logmH2Q2)*logmsb22Q2 + 6*logmsb12Q2*(-9 + logmH2Q2 +
     2*logmsb22Q2) + 9*power2(logmsb12Q2) + 3*power2(logmsb22Q2) + 2*(42 +
     power2(Pi))) - mH2*(msb12*(42 + 6*logmH2Q2*(-3 + logmsb12Q2) - 18*
     logmsb12Q2 + 3*power2(logmH2Q2) + 3*power2(logmsb12Q2) + power2(Pi)) +
     msb22*(42 - 36*logmsb12Q2 + 6*logmH2Q2*(-3 + 2*logmsb12Q2 - logmsb22Q2) +
     18*logmsb22Q2 + 3*power2(logmH2Q2) + 6*power2(logmsb12Q2) - 3*power2(
     logmsb22Q2) + power2(Pi))))*(-1 + power2(sa))*power2(sa))/(8.*mH2) - (3*
     power2(Ab)*(3*(logmsb12Q2 - logmsb22Q2)*(-6 + 2*logmH2Q2 + logmsb12Q2 +
     logmsb22Q2)*msb24 + power2(mH2)*(42 + 6*logmH2Q2*(-3 + logmsb12Q2) - 18*
     logmsb12Q2 + 3*power2(logmH2Q2) + 3*power2(logmsb12Q2) + power2(Pi)) -
     msb12*msb22*(-6*(3 + logmH2Q2)*logmsb22Q2 + 6*logmsb12Q2*(-9 + logmH2Q2 +
     2*logmsb22Q2) + 9*power2(logmsb12Q2) + 3*power2(logmsb22Q2) + 2*(42 +
     power2(Pi))) - mH2*(msb12*(42 + 6*logmH2Q2*(-3 + logmsb12Q2) - 18*
     logmsb12Q2 + 3*power2(logmH2Q2) + 3*power2(logmsb12Q2) + power2(Pi)) +
     msb22*(42 - 36*logmsb12Q2 + 6*logmH2Q2*(-3 + 2*logmsb12Q2 - logmsb22Q2) +
     18*logmsb22Q2 + 3*power2(logmH2Q2) + 6*power2(logmsb12Q2) - 3*power2(
     logmsb22Q2) + power2(Pi))))*power2(-1 + power2(sa)))/(8.*mH2)) + Fin3(
     msb22,mH2,msb12,Q2)*((-3*mu2*(2*mh2*mH2*(-3 + 4*invdmhH*mH2) + 3*power2(
     mh2) - 2*(-1 + invdmhH*mH2)*power2(mH2))*(-1 + power2(sa))*power2(sa))/(4.
     *power2(mh2)*power2(mH2)) + (3*power2(Ab)*(-1 + power2(sa))*(3*power2(mh2)
     *(-1 + power2(sa)) + 2*mh2*mH2*(-3 + 4*invdmhH*mH2)*power2(sa) - 2*(-1 +
     invdmhH*mH2)*power2(mH2)*power2(sa)))/(4.*power2(mh2)*power2(mH2)) + (3*Ab
     *ca*mu*sa*(3*power2(mh2)*(-1 + power2(sa)) + mh2*mH2*(-3 + 4*invdmhH*mH2)*
     (-1 + 2*power2(sa)) - (-1 + invdmhH*mH2)*power2(mH2)*(-1 + 2*power2(sa))))
     /(2.*power2(mh2)*power2(mH2)) + DeltaInv(msb22,mH2,msb12)*((-9*Ab*ca*(mH2
     - msb12 - msb22)*mu*sa*(-1 + power2(sa)))/(2.*mH2) + (9*(mH2 - msb12 -
     msb22)*mu2*(-1 + power2(sa))*power2(sa))/(4.*mH2) - (9*(mH2 - msb12 -
     msb22)*power2(Ab)*power2(-1 + power2(sa)))/(4.*mH2))) + Fin20(msb22,msb12,
     Q2)*((9*mu2*power2(mh2 - mH2)*(-1 + power2(sa))*power2(sa))/(4.*power2(mh2
     )*power2(mH2)) - (9*Ab*ca*mu*sa*(power2(mh2)*(-1 + power2(sa)) + power2(
     mH2)*power2(sa) + mh2*(mH2 - 2*mH2*power2(sa))))/(2.*power2(mh2)*power2(
     mH2)) - (3*power2(Ab)*(-(power2(mh2)*power2(mH2)) + 3*power2(mA2)*power2(
     mh2 - mh2*power2(sa) + mH2*power2(sa))))/(4.*power2(mA2)*power2(mh2)*
     power2(mH2))) + (Fin20(msb12,mu2,Q2)*(3*(msb12 - 3*mu2)*power2(
     invdmstau1mu) + 3*(msb12 - 2*invdmstau2mu*msb12*mu2 + mu2*(-3 + 2*
     invdmstau2mu*mu2))*power2(invdmstau2mu) + 6*mu2*(-msb12 + mu2)*power3(
     invdmstau1mu)))/8. + (Fin20(msb22,mu2,Q2)*(3*(msb22 - 3*mu2)*power2(
     invdmstau1mu) + 3*(msb22 - 2*invdmstau2mu*msb22*mu2 + mu2*(-3 + 2*
     invdmstau2mu*mu2))*power2(invdmstau2mu) + 6*mu2*(-msb22 + mu2)*power3(
     invdmstau1mu)))/8. + DeltaInv(msb22,mh2,msb12)*((3*mu2*(3*(logmsb12Q2 -
     logmsb22Q2)*(-6 + 2*logmh2Q2 + logmsb12Q2 + logmsb22Q2)*msb24 + power2(mh2
     )*(42 + 6*logmh2Q2*(-3 + logmsb12Q2) - 18*logmsb12Q2 + 3*power2(logmh2Q2)
     + 3*power2(logmsb12Q2) + power2(Pi)) - msb12*msb22*(-6*(3 + logmh2Q2)*
     logmsb22Q2 + 6*logmsb12Q2*(-9 + logmh2Q2 + 2*logmsb22Q2) + 9*power2(
     logmsb12Q2) + 3*power2(logmsb22Q2) + 2*(42 + power2(Pi))) - mh2*(msb12*(42
      + 6*logmh2Q2*(-3 + logmsb12Q2) - 18*logmsb12Q2 + 3*power2(logmh2Q2) + 3*
     power2(logmsb12Q2) + power2(Pi)) + msb22*(42 - 36*logmsb12Q2 + 6*logmh2Q2*
     (-3 + 2*logmsb12Q2 - logmsb22Q2) + 18*logmsb22Q2 + 3*power2(logmh2Q2) + 6*
     power2(logmsb12Q2) - 3*power2(logmsb22Q2) + power2(Pi))))*(-1 + power2(sa)
     )*power2(sa))/(8.*mh2) - (3*Ab*ca*mu*(3*(logmsb12Q2 - logmsb22Q2)*(-6 + 2*
     logmh2Q2 + logmsb12Q2 + logmsb22Q2)*msb24 + power2(mh2)*(42 + 6*logmh2Q2*(
     -3 + logmsb12Q2) - 18*logmsb12Q2 + 3*power2(logmh2Q2) + 3*power2(
     logmsb12Q2) + power2(Pi)) - msb12*msb22*(-6*(3 + logmh2Q2)*logmsb22Q2 + 6*
     logmsb12Q2*(-9 + logmh2Q2 + 2*logmsb22Q2) + 9*power2(logmsb12Q2) + 3*
     power2(logmsb22Q2) + 2*(42 + power2(Pi))) - mh2*(msb12*(42 + 6*logmh2Q2*(-
     3 + logmsb12Q2) - 18*logmsb12Q2 + 3*power2(logmh2Q2) + 3*power2(logmsb12Q2
     ) + power2(Pi)) + msb22*(42 - 36*logmsb12Q2 + 6*logmh2Q2*(-3 + 2*
     logmsb12Q2 - logmsb22Q2) + 18*logmsb22Q2 + 3*power2(logmh2Q2) + 6*power2(
     logmsb12Q2) - 3*power2(logmsb22Q2) + power2(Pi))))*power3(sa))/(4.*mh2) -
     (3*power2(Ab)*(3*(logmsb12Q2 - logmsb22Q2)*(-6 + 2*logmh2Q2 + logmsb12Q2 +
     logmsb22Q2)*msb24 + power2(mh2)*(42 + 6*logmh2Q2*(-3 + logmsb12Q2) - 18*
     logmsb12Q2 + 3*power2(logmh2Q2) + 3*power2(logmsb12Q2) + power2(Pi)) -
     msb12*msb22*(-6*(3 + logmh2Q2)*logmsb22Q2 + 6*logmsb12Q2*(-9 + logmh2Q2 +
     2*logmsb22Q2) + 9*power2(logmsb12Q2) + 3*power2(logmsb22Q2) + 2*(42 +
     power2(Pi))) - mh2*(msb12*(42 + 6*logmh2Q2*(-3 + logmsb12Q2) - 18*
     logmsb12Q2 + 3*power2(logmh2Q2) + 3*power2(logmsb12Q2) + power2(Pi)) +
     msb22*(42 - 36*logmsb12Q2 + 6*logmh2Q2*(-3 + 2*logmsb12Q2 - logmsb22Q2) +
     18*logmsb22Q2 + 3*power2(logmh2Q2) + 6*power2(logmsb12Q2) - 3*power2(
     logmsb22Q2) + power2(Pi))))*power4(sa))/(8.*mh2)) + Fin3(msb22,mh2,msb12,
     Q2)*((9*(-1 + 2*invdmhH*mh2)*mu2*(-1 + power2(sa))*power2(sa))/(4.*power2(
     mh2)) + (9*power2(Ab)*power2(sa)*(-2*invdmhH*mh2*(-1 + power2(sa)) +
     power2(sa)))/(4.*power2(mh2)) - (9*Ab*ca*mu*(invdmhH*mh2*sa*(-1 + 2*power2
     (sa)) - power3(sa)))/(2.*power2(mh2)) + DeltaInv(msb22,mh2,msb12)*((9*(mh2
      - msb12 - msb22)*mu2*(-1 + power2(sa))*power2(sa))/(4.*mh2) + (9*Ab*ca*(-
     mh2 + msb12 + msb22)*mu*power3(sa))/(2.*mh2) + (9*(-mh2 + msb12 + msb22)*
     power2(Ab)*power4(sa))/(4.*mh2))) + ((-12*power2(mt2)*(14 - 12*logmt2Q2 +
     4*power2(logmt2Q2) + power2(Pi)))/power2(mC2) - (4*(-1 + invdmhH*mH2)*(-3*
     ((-7 + 4*logmsb12Q2)*msb12 + (-7 + 4*logmsb22Q2)*msb22)*mu2 + power2(mH2)*
     (-84 - 6*logmH2Q2*(-11 + 6*invdmhH*(msb12 + msb22 - mu2)) + 6*(-3 + 2*
     invdmhH*(msb12 + msb22 - mu2))*power2(logmH2Q2) - 5*power2(Pi) + 2*invdmhH
     *(msb12 + msb22 - mu2)*(21 + power2(Pi))) - mH2*(3*(5 - 12*logmsb12Q2 +
     logmH2Q2*(-4 + 8*logmsb12Q2))*mu2 + 2*msb12*(33 - 18*logmH2Q2 - 12*
     logmsb12Q2 + 6*power2(logmH2Q2) + 6*power2(logmsb12Q2) + 2*power2(Pi)) + 2
     *msb22*(33 - 18*invdmsb*logmsb12Q2*mu2 + 6*logmsb22Q2*(-2 + 3*invdmsb*mu2)
     + 6*logmH2Q2*(-3 + 2*invdmsb*(logmsb12Q2 - logmsb22Q2)*mu2) + 6*power2(
     logmH2Q2) + 6*power2(logmsb22Q2) + 2*power2(Pi))))*(-1 + power2(sa))*
     power2(sa))/power2(mh2) + (4*(-3*mu2*(45 + 12*logmh2Q2*(-2 + logmsb12Q2) -
     28*logmsb12Q2 - 4*logmH2Q2*(-7 + 6*logmsb12Q2) + 6*power2(logmh2Q2) + 6*
     power2(logmsb12Q2) + 2*power2(Pi)) + 3*msb12*(42 - 20*logmH2Q2 + 7*invdmhH
     *mu2 - 4*logmsb12Q2*(6 + invdmhH*mu2) + 12*power2(logmH2Q2) + 12*power2(
     logmsb12Q2) + 4*power2(Pi)) + 3*msb22*(42 + 7*invdmhH*mu2 - 8*invdmsb*
     logmsb12Q2*mu2 - 4*logmsb22Q2*(6 + invdmhH*mu2 - 2*invdmsb*mu2) + 4*
     logmH2Q2*(-5 + 6*invdmsb*(logmsb12Q2 - logmsb22Q2)*mu2) + 12*power2(
     logmH2Q2) + 12*power2(logmsb22Q2) + 4*power2(Pi)) + invdmhH*power2(mH2)*(-
     6*logmH2Q2*(-43 + 22*invdmhH*(msb12 + msb22 - mu2)) + (-78 + 60*invdmhH*(
     msb12 + msb22 - mu2))*power2(logmH2Q2) + 2*invdmhH*(msb12 + msb22 - mu2)*(
     69 + 5*power2(Pi)) - 3*(103 + 7*power2(Pi))) + mH2*(225 + 96*logmH2Q2*(-2
     + invdmhH*(2*msb12 - logmsb12Q2*mu2 + msb22*(2 + invdmsb*(-logmsb12Q2 +
     logmsb22Q2)*mu2))) - 12*(-5 + invdmhH*(8*msb12 + 8*msb22 - 4*mu2))*power2(
     logmH2Q2) + 16*power2(Pi) - 8*invdmhH*(-(mu2*(12*logmsb12Q2 + power2(Pi)))
     + 3*msb12*(12 - 4*logmsb12Q2 + 2*power2(logmsb12Q2) + power2(Pi)) + 3*
     msb22*(12 - 4*logmsb22Q2 - 4*invdmsb*logmsb12Q2*mu2 + 4*invdmsb*logmsb22Q2
     *mu2 + 2*power2(logmsb22Q2) + power2(Pi)))))*(-1 + power2(sa))*power2(sa))
     /mh2 - (24*(2*(-1 + logmsb12Q2)*msb12 - 2*mt2 + 6*logmsb22Q2*mt2 - 2*
     logmC2Q2*logmsb22Q2*mt2 - 4*logmt2Q2*mt2 + 2*logmC2Q2*logmt2Q2*mt2 - mt2*
     power2(logmsb22Q2) + mt2*power2(logmt2Q2) + 2*(-1 + logmst12Q2)*mst12*(-1
     + power2(snt)) + 2*mst22*power2(snt) - 2*logmst22Q2*mst22*power2(snt)))/
     mC2 - (24*mu2*(36 + 6*logmH2Q2*(-2 + logmsb12Q2) + 6*invdmsb*logmsb22Q2*
     msb22 - 6*logmsb12Q2*(4 + invdmsb*msb22) + 3*power2(logmH2Q2) + 3*power2(
     logmsb12Q2) + power2(Pi))*(-1 + power2(sa))*power2(sa) + 4*invdmhH*power2(
     mH2)*(-225 + 144*logmh2Q2 + 96*invdmhH*msb12 + 96*invdmhH*msb22 - 48*
     logmH2Q2*(-1 + 2*invdmhH*(msb12 + msb22 - mu2)) - 96*invdmhH*mu2 - 36*
     power2(logmh2Q2) + 24*(-1 + 2*invdmhH*(msb12 + msb22 - mu2))*power2(
     logmH2Q2) - 16*power2(Pi) + 8*invdmhH*msb12*power2(Pi) + 8*invdmhH*msb22*
     power2(Pi) - 8*invdmhH*mu2*power2(Pi))*(-1 + power2(sa))*power2(sa) + mH2*
     (192 - 48*logmA2Q2 - 96*logmC2Q2 + 144*logmH2Q2 + 48*logmsntau2Q2 - 96*
     logmstau12Q2 - 96*logmstau22Q2 + 144*logmt2Q2 - 48*logmsntau2Q2*logmt2Q2 -
     144*invdmstau1mu*msb12 - 144*invdmstau2mu*msb12 + 96*invdmstau1mu*
     logmsb12Q2*msb12 + 96*invdmstau2mu*logmsb12Q2*msb12 + 48*invdmstau1mu*
     logmstau12Q2*msb12 - 24*invdmstau1mu*logmsb12Q2*logmstau12Q2*msb12 + 48*
     invdmstau2mu*logmstau22Q2*msb12 - 24*invdmstau2mu*logmsb12Q2*logmstau22Q2*
     msb12 + 63*invdmsb1stau1*invdmstau1mu*msb14 + 63*invdmsb1stau2*
     invdmstau2mu*msb14 - 36*invdmsb1stau1*invdmstau1mu*logmsb12Q2*msb14 - 36*
     invdmsb1stau2*invdmstau2mu*logmsb12Q2*msb14 + 192*invdmsntau2mu*msb22 -
     144*invdmstau1mu*msb22 - 144*invdmstau2mu*msb22 - 48*invdmsntau2mu*
     logmsb22Q2*msb22 + 96*invdmstau1mu*logmsb22Q2*msb22 + 96*invdmstau2mu*
     logmsb22Q2*msb22 + 48*invdmsntau2mu*logmsntau2Q2*msb22 - 24*invdmsntau2mu*
     logmsb22Q2*logmsntau2Q2*msb22 + 48*invdmstau1mu*logmstau12Q2*msb22 - 24*
     invdmstau1mu*logmsb22Q2*logmstau12Q2*msb22 + 48*invdmstau2mu*logmstau22Q2*
     msb22 - 24*invdmstau2mu*logmsb22Q2*logmstau22Q2*msb22 - 144*invdmsntau2mu*
     logmt2Q2*msb22 + 48*invdmsntau2mu*logmsb22Q2*logmt2Q2*msb22 - 63*
     invdmsb2stau1*invdmstau1mu*msb24 - 63*invdmsb2stau2*invdmstau2mu*msb24 +
     36*invdmsb2stau1*invdmstau1mu*logmsb22Q2*msb24 + 36*invdmsb2stau2*
     invdmstau2mu*logmsb22Q2*msb24 - 144*invdmsntau2mu*mst12 + 48*invdmsntau2mu
     *logmsntau2Q2*mst12 + 96*invdmsntau2mu*logmst12Q2*mst12 - 24*invdmsntau2mu
     *logmsntau2Q2*logmst12Q2*mst12 + 63*invdmsntau2mu*invdmsntaust1*mst14 - 36
     *invdmsntau2mu*invdmsntaust1*logmst12Q2*mst14 + 144*invdmsntau2mu*mt2 - 48
     *invdmsntau2mu*logmsntau2Q2*mt2 - 96*invdmsntau2mu*logmt2Q2*mt2 + 24*
     invdmsntau2mu*logmsntau2Q2*logmt2Q2*mt2 + 480*invdmsntau2mu*mu2 - 192*
     invdmstau1mu*mu2 - 192*invdmstau2mu*mu2 - 192*invdmsntau2mu*logmsntau2Q2*
     mu2 + 96*invdmstau1mu*logmstau12Q2*mu2 + 96*invdmstau2mu*logmstau22Q2*mu2
     - 288*invdmsntau2mu*logmt2Q2*mu2 + 96*invdmsntau2mu*logmsntau2Q2*logmt2Q2*
     mu2 - 63*mst14*power2(invdmsntau2mu) + 36*logmst12Q2*mst14*power2(
     invdmsntau2mu) - 384*msb22*mu2*power2(invdmsntau2mu) + 96*logmsb22Q2*msb22
     *mu2*power2(invdmsntau2mu) - 144*logmsntau2Q2*msb22*mu2*power2(
     invdmsntau2mu) + 72*logmsb22Q2*logmsntau2Q2*msb22*mu2*power2(invdmsntau2mu
     ) + 288*logmt2Q2*msb22*mu2*power2(invdmsntau2mu) - 96*logmsb22Q2*logmt2Q2*
     msb22*mu2*power2(invdmsntau2mu) + 48*logmu2Q2*msb22*mu2*power2(
     invdmsntau2mu) - 24*logmsb22Q2*logmu2Q2*msb22*mu2*power2(invdmsntau2mu) +
     288*mst12*mu2*power2(invdmsntau2mu) - 144*logmsntau2Q2*mst12*mu2*power2(
     invdmsntau2mu) - 192*logmst12Q2*mst12*mu2*power2(invdmsntau2mu) + 72*
     logmsntau2Q2*logmst12Q2*mst12*mu2*power2(invdmsntau2mu) + 48*logmu2Q2*
     mst12*mu2*power2(invdmsntau2mu) - 24*logmst12Q2*logmu2Q2*mst12*mu2*power2(
     invdmsntau2mu) - 63*invdmsntaust1*mst14*mu2*power2(invdmsntau2mu) + 36*
     invdmsntaust1*logmst12Q2*mst14*mu2*power2(invdmsntau2mu) - 288*mt2*mu2*
     power2(invdmsntau2mu) + 144*logmsntau2Q2*mt2*mu2*power2(invdmsntau2mu) +
     192*logmt2Q2*mt2*mu2*power2(invdmsntau2mu) - 72*logmsntau2Q2*logmt2Q2*mt2*
     mu2*power2(invdmsntau2mu) - 48*logmu2Q2*mt2*mu2*power2(invdmsntau2mu) + 24
     *logmt2Q2*logmu2Q2*mt2*mu2*power2(invdmsntau2mu) - 63*msb14*power2(
     invdmstau1mu) + 36*logmsb12Q2*msb14*power2(invdmstau1mu) - 63*msb24*power2
     (invdmstau1mu) + 36*logmsb22Q2*msb24*power2(invdmstau1mu) + 288*msb12*mu2*
     power2(invdmstau1mu) - 192*logmsb12Q2*msb12*mu2*power2(invdmstau1mu) - 144
     *logmstau12Q2*msb12*mu2*power2(invdmstau1mu) + 72*logmsb12Q2*logmstau12Q2*
     msb12*mu2*power2(invdmstau1mu) + 48*logmu2Q2*msb12*mu2*power2(invdmstau1mu
     ) - 24*logmsb12Q2*logmu2Q2*msb12*mu2*power2(invdmstau1mu) - 63*
     invdmsb1stau1*msb14*mu2*power2(invdmstau1mu) + 36*invdmsb1stau1*logmsb12Q2
     *msb14*mu2*power2(invdmstau1mu) + 288*msb22*mu2*power2(invdmstau1mu) - 192
     *logmsb22Q2*msb22*mu2*power2(invdmstau1mu) - 144*logmstau12Q2*msb22*mu2*
     power2(invdmstau1mu) + 72*logmsb22Q2*logmstau12Q2*msb22*mu2*power2(
     invdmstau1mu) + 48*logmu2Q2*msb22*mu2*power2(invdmstau1mu) - 24*logmsb22Q2
     *logmu2Q2*msb22*mu2*power2(invdmstau1mu) + 63*invdmsb2stau1*msb24*mu2*
     power2(invdmstau1mu) - 36*invdmsb2stau1*logmsb22Q2*msb24*mu2*power2(
     invdmstau1mu) - 63*msb14*power2(invdmstau2mu) + 36*logmsb12Q2*msb14*power2
     (invdmstau2mu) - 63*msb24*power2(invdmstau2mu) + 36*logmsb22Q2*msb24*
     power2(invdmstau2mu) + 288*msb12*mu2*power2(invdmstau2mu) - 192*logmsb12Q2
     *msb12*mu2*power2(invdmstau2mu) - 144*logmstau22Q2*msb12*mu2*power2(
     invdmstau2mu) + 72*logmsb12Q2*logmstau22Q2*msb12*mu2*power2(invdmstau2mu)
     + 48*logmu2Q2*msb12*mu2*power2(invdmstau2mu) - 24*logmsb12Q2*logmu2Q2*
     msb12*mu2*power2(invdmstau2mu) - 63*invdmsb1stau2*msb14*mu2*power2(
     invdmstau2mu) + 36*invdmsb1stau2*logmsb12Q2*msb14*mu2*power2(invdmstau2mu)
     + 288*msb22*mu2*power2(invdmstau2mu) - 192*logmsb22Q2*msb22*mu2*power2(
     invdmstau2mu) - 144*logmstau22Q2*msb22*mu2*power2(invdmstau2mu) + 72*
     logmsb22Q2*logmstau22Q2*msb22*mu2*power2(invdmstau2mu) + 48*logmu2Q2*msb22
     *mu2*power2(invdmstau2mu) - 24*logmsb22Q2*logmu2Q2*msb22*mu2*power2(
     invdmstau2mu) + 63*invdmsb2stau2*msb24*mu2*power2(invdmstau2mu) - 36*
     invdmsb2stau2*logmsb22Q2*msb24*mu2*power2(invdmstau2mu) + 24*power2(
     logmA2Q2) + 24*power2(logmC2Q2) - 72*power2(logmH2Q2) - 12*invdmstau1mu*
     msb12*power2(logmsb12Q2) - 12*invdmstau2mu*msb12*power2(logmsb12Q2) + 24*
     msb12*mu2*power2(invdmstau1mu)*power2(logmsb12Q2) + 24*msb12*mu2*power2(
     invdmstau2mu)*power2(logmsb12Q2) + 12*invdmsntau2mu*msb22*power2(
     logmsb22Q2) - 12*invdmstau1mu*msb22*power2(logmsb22Q2) - 12*invdmstau2mu*
     msb22*power2(logmsb22Q2) - 24*msb22*mu2*power2(invdmsntau2mu)*power2(
     logmsb22Q2) + 24*msb22*mu2*power2(invdmstau1mu)*power2(logmsb22Q2) + 24*
     msb22*mu2*power2(invdmstau2mu)*power2(logmsb22Q2) - 12*invdmsntau2mu*msb22
     *power2(logmsntau2Q2) - 12*invdmsntau2mu*mst12*power2(logmsntau2Q2) + 12*
     invdmsntau2mu*mt2*power2(logmsntau2Q2) + 48*invdmsntau2mu*mu2*power2(
     logmsntau2Q2) + 36*msb22*mu2*power2(invdmsntau2mu)*power2(logmsntau2Q2) +
     36*mst12*mu2*power2(invdmsntau2mu)*power2(logmsntau2Q2) - 36*mt2*mu2*
     power2(invdmsntau2mu)*power2(logmsntau2Q2) - 12*invdmsntau2mu*mst12*power2
     (logmst12Q2) + 24*mst12*mu2*power2(invdmsntau2mu)*power2(logmst12Q2) + 24*
     power2(logmstau12Q2) - 12*invdmstau1mu*msb12*power2(logmstau12Q2) - 12*
     invdmstau1mu*msb22*power2(logmstau12Q2) + 36*msb12*mu2*power2(invdmstau1mu
     )*power2(logmstau12Q2) + 36*msb22*mu2*power2(invdmstau1mu)*power2(
     logmstau12Q2) + 24*power2(logmstau22Q2) - 12*invdmstau2mu*msb12*power2(
     logmstau22Q2) - 12*invdmstau2mu*msb22*power2(logmstau22Q2) + 36*msb12*mu2*
     power2(invdmstau2mu)*power2(logmstau22Q2) + 36*msb22*mu2*power2(
     invdmstau2mu)*power2(logmstau22Q2) - 24*power2(logmt2Q2) + 24*
     invdmsntau2mu*msb22*power2(logmt2Q2) + 12*invdmsntau2mu*mt2*power2(
     logmt2Q2) + 48*invdmsntau2mu*mu2*power2(logmt2Q2) - 48*msb22*mu2*power2(
     invdmsntau2mu)*power2(logmt2Q2) - 24*mt2*mu2*power2(invdmsntau2mu)*power2(
     logmt2Q2) - 12*msb22*mu2*power2(invdmsntau2mu)*power2(logmu2Q2) - 12*mst12
     *mu2*power2(invdmsntau2mu)*power2(logmu2Q2) + 12*mt2*mu2*power2(
     invdmsntau2mu)*power2(logmu2Q2) - 12*msb12*mu2*power2(invdmstau1mu)*power2
     (logmu2Q2) - 12*msb22*mu2*power2(invdmstau1mu)*power2(logmu2Q2) - 12*msb12
     *mu2*power2(invdmstau2mu)*power2(logmu2Q2) - 12*msb22*mu2*power2(
     invdmstau2mu)*power2(logmu2Q2) - 480*power2(invdmsntau2mu)*power2(mu2) +
     240*logmsntau2Q2*power2(invdmsntau2mu)*power2(mu2) + 288*logmt2Q2*power2(
     invdmsntau2mu)*power2(mu2) - 48*logmsntau2Q2*logmt2Q2*power2(invdmsntau2mu
     )*power2(mu2) - 48*logmu2Q2*power2(invdmsntau2mu)*power2(mu2) - 48*
     logmt2Q2*logmu2Q2*power2(invdmsntau2mu)*power2(mu2) + 192*power2(
     invdmstau1mu)*power2(mu2) + 96*logmstau12Q2*power2(invdmstau1mu)*power2(
     mu2) - 192*logmu2Q2*power2(invdmstau1mu)*power2(mu2) + 192*power2(
     invdmstau2mu)*power2(mu2) + 96*logmstau22Q2*power2(invdmstau2mu)*power2(
     mu2) - 192*logmu2Q2*power2(invdmstau2mu)*power2(mu2) - 96*power2(
     invdmsntau2mu)*power2(logmsntau2Q2)*power2(mu2) - 72*power2(invdmstau1mu)*
     power2(logmstau12Q2)*power2(mu2) - 72*power2(invdmstau2mu)*power2(
     logmstau22Q2)*power2(mu2) - 48*power2(invdmsntau2mu)*power2(logmt2Q2)*
     power2(mu2) + 48*power2(invdmsntau2mu)*power2(logmu2Q2)*power2(mu2) + 72*
     power2(invdmstau1mu)*power2(logmu2Q2)*power2(mu2) + 72*power2(invdmstau2mu
     )*power2(logmu2Q2)*power2(mu2) - 8*power2(Pi) - 4*invdmstau1mu*msb12*
     power2(Pi) - 4*invdmstau2mu*msb12*power2(Pi) + 4*invdmsntau2mu*msb22*
     power2(Pi) - 4*invdmstau1mu*msb22*power2(Pi) - 4*invdmstau2mu*msb22*power2
     (Pi) - 4*invdmsntau2mu*mst12*power2(Pi) + 4*invdmsntau2mu*mt2*power2(Pi) +
     16*invdmsntau2mu*mu2*power2(Pi) - 8*msb22*mu2*power2(invdmsntau2mu)*power2
     (Pi) + 8*mst12*mu2*power2(invdmsntau2mu)*power2(Pi) - 8*mt2*mu2*power2(
     invdmsntau2mu)*power2(Pi) + 8*msb12*mu2*power2(invdmstau1mu)*power2(Pi) +
     8*msb22*mu2*power2(invdmstau1mu)*power2(Pi) + 8*msb12*mu2*power2(
     invdmstau2mu)*power2(Pi) + 8*msb22*mu2*power2(invdmstau2mu)*power2(Pi) -
     16*power2(invdmsntau2mu)*power2(mu2)*power2(Pi) - 576*power2(sa) + 576*
     logmh2Q2*power2(sa) - 288*logmH2Q2*power2(sa) + 504*invdmhH*msb12*power2(
     sa) - 240*invdmhH*logmH2Q2*msb12*power2(sa) - 288*invdmhH*logmsb12Q2*msb12
     *power2(sa) + 504*invdmhH*msb22*power2(sa) - 240*invdmhH*logmH2Q2*msb22*
     power2(sa) - 288*invdmhH*logmsb22Q2*msb22*power2(sa) + 324*invdmhH*mu2*
     power2(sa) - 576*invdmhH*logmh2Q2*mu2*power2(sa) + 240*invdmhH*logmH2Q2*
     mu2*power2(sa) - 240*invdmhH*logmsb12Q2*mu2*power2(sa) + 288*invdmhH*
     logmh2Q2*logmsb12Q2*mu2*power2(sa) - 240*invdmhH*invdmsb*logmsb12Q2*msb22*
     mu2*power2(sa) + 288*invdmhH*invdmsb*logmh2Q2*logmsb12Q2*msb22*mu2*power2(
     sa) + 240*invdmhH*invdmsb*logmsb22Q2*msb22*mu2*power2(sa) - 288*invdmhH*
     invdmsb*logmh2Q2*logmsb22Q2*msb22*mu2*power2(sa) - 144*power2(logmh2Q2)*
     power2(sa) + 144*invdmhH*mu2*power2(logmh2Q2)*power2(sa) + 144*power2(
     logmH2Q2)*power2(sa) + 144*invdmhH*msb12*power2(logmH2Q2)*power2(sa) + 144
     *invdmhH*msb22*power2(logmH2Q2)*power2(sa) - 144*invdmhH*mu2*power2(
     logmH2Q2)*power2(sa) + 144*invdmhH*msb12*power2(logmsb12Q2)*power2(sa) +
     144*invdmhH*msb22*power2(logmsb22Q2)*power2(sa) + 48*invdmhH*msb12*power2(
     Pi)*power2(sa) + 48*invdmhH*msb22*power2(Pi)*power2(sa) + 144*
     invdmsntau2mu*mst12*power2(snt) - 48*invdmsntau2mu*logmsntau2Q2*mst12*
     power2(snt) - 96*invdmsntau2mu*logmst12Q2*mst12*power2(snt) + 24*
     invdmsntau2mu*logmsntau2Q2*logmst12Q2*mst12*power2(snt) - 63*invdmsntau2mu
     *invdmsntaust1*mst14*power2(snt) + 36*invdmsntau2mu*invdmsntaust1*
     logmst12Q2*mst14*power2(snt) - 144*invdmsntau2mu*mst22*power2(snt) + 48*
     invdmsntau2mu*logmsntau2Q2*mst22*power2(snt) + 96*invdmsntau2mu*logmst22Q2
     *mst22*power2(snt) - 24*invdmsntau2mu*logmsntau2Q2*logmst22Q2*mst22*power2
     (snt) + 63*invdmsntau2mu*invdmsntaust2*mst24*power2(snt) - 36*
     invdmsntau2mu*invdmsntaust2*logmst22Q2*mst24*power2(snt) + 63*mst14*power2
     (invdmsntau2mu)*power2(snt) - 36*logmst12Q2*mst14*power2(invdmsntau2mu)*
     power2(snt) - 63*mst24*power2(invdmsntau2mu)*power2(snt) + 36*logmst22Q2*
     mst24*power2(invdmsntau2mu)*power2(snt) - 288*mst12*mu2*power2(
     invdmsntau2mu)*power2(snt) + 144*logmsntau2Q2*mst12*mu2*power2(
     invdmsntau2mu)*power2(snt) + 192*logmst12Q2*mst12*mu2*power2(invdmsntau2mu
     )*power2(snt) - 72*logmsntau2Q2*logmst12Q2*mst12*mu2*power2(invdmsntau2mu)
     *power2(snt) - 48*logmu2Q2*mst12*mu2*power2(invdmsntau2mu)*power2(snt) +
     24*logmst12Q2*logmu2Q2*mst12*mu2*power2(invdmsntau2mu)*power2(snt) + 63*
     invdmsntaust1*mst14*mu2*power2(invdmsntau2mu)*power2(snt) - 36*
     invdmsntaust1*logmst12Q2*mst14*mu2*power2(invdmsntau2mu)*power2(snt) + 288
     *mst22*mu2*power2(invdmsntau2mu)*power2(snt) - 144*logmsntau2Q2*mst22*mu2*
     power2(invdmsntau2mu)*power2(snt) - 192*logmst22Q2*mst22*mu2*power2(
     invdmsntau2mu)*power2(snt) + 72*logmsntau2Q2*logmst22Q2*mst22*mu2*power2(
     invdmsntau2mu)*power2(snt) + 48*logmu2Q2*mst22*mu2*power2(invdmsntau2mu)*
     power2(snt) - 24*logmst22Q2*logmu2Q2*mst22*mu2*power2(invdmsntau2mu)*
     power2(snt) - 63*invdmsntaust2*mst24*mu2*power2(invdmsntau2mu)*power2(snt)
     + 36*invdmsntaust2*logmst22Q2*mst24*mu2*power2(invdmsntau2mu)*power2(snt)
     + 12*invdmsntau2mu*mst12*power2(logmsntau2Q2)*power2(snt) - 12*
     invdmsntau2mu*mst22*power2(logmsntau2Q2)*power2(snt) - 36*mst12*mu2*power2
     (invdmsntau2mu)*power2(logmsntau2Q2)*power2(snt) + 36*mst22*mu2*power2(
     invdmsntau2mu)*power2(logmsntau2Q2)*power2(snt) + 12*invdmsntau2mu*mst12*
     power2(logmst12Q2)*power2(snt) - 24*mst12*mu2*power2(invdmsntau2mu)*power2
     (logmst12Q2)*power2(snt) - 12*invdmsntau2mu*mst22*power2(logmst22Q2)*
     power2(snt) + 24*mst22*mu2*power2(invdmsntau2mu)*power2(logmst22Q2)*power2
     (snt) + 12*mst12*mu2*power2(invdmsntau2mu)*power2(logmu2Q2)*power2(snt) -
     12*mst22*mu2*power2(invdmsntau2mu)*power2(logmu2Q2)*power2(snt) + 4*
     invdmsntau2mu*mst12*power2(Pi)*power2(snt) - 4*invdmsntau2mu*mst22*power2(
     Pi)*power2(snt) - 8*mst12*mu2*power2(invdmsntau2mu)*power2(Pi)*power2(snt)
     + 8*mst22*mu2*power2(invdmsntau2mu)*power2(Pi)*power2(snt) + 96*
     logmsntau2Q2*msb22*power2(mu2)*power3(invdmsntau2mu) - 48*logmsb22Q2*
     logmsntau2Q2*msb22*power2(mu2)*power3(invdmsntau2mu) - 96*logmu2Q2*msb22*
     power2(mu2)*power3(invdmsntau2mu) + 48*logmsb22Q2*logmu2Q2*msb22*power2(
     mu2)*power3(invdmsntau2mu) + 96*logmsntau2Q2*mst12*power2(mu2)*power3(
     invdmsntau2mu) - 48*logmsntau2Q2*logmst12Q2*mst12*power2(mu2)*power3(
     invdmsntau2mu) - 96*logmu2Q2*mst12*power2(mu2)*power3(invdmsntau2mu) + 48*
     logmst12Q2*logmu2Q2*mst12*power2(mu2)*power3(invdmsntau2mu) - 96*
     logmsntau2Q2*mt2*power2(mu2)*power3(invdmsntau2mu) + 48*logmsntau2Q2*
     logmt2Q2*mt2*power2(mu2)*power3(invdmsntau2mu) + 96*logmu2Q2*mt2*power2(
     mu2)*power3(invdmsntau2mu) - 48*logmt2Q2*logmu2Q2*mt2*power2(mu2)*power3(
     invdmsntau2mu) - 24*msb22*power2(logmsntau2Q2)*power2(mu2)*power3(
     invdmsntau2mu) - 24*mst12*power2(logmsntau2Q2)*power2(mu2)*power3(
     invdmsntau2mu) + 24*mt2*power2(logmsntau2Q2)*power2(mu2)*power3(
     invdmsntau2mu) + 24*msb22*power2(logmu2Q2)*power2(mu2)*power3(
     invdmsntau2mu) + 24*mst12*power2(logmu2Q2)*power2(mu2)*power3(
     invdmsntau2mu) - 24*mt2*power2(logmu2Q2)*power2(mu2)*power3(invdmsntau2mu)
     - 96*logmsntau2Q2*mst12*power2(mu2)*power2(snt)*power3(invdmsntau2mu) + 48
     *logmsntau2Q2*logmst12Q2*mst12*power2(mu2)*power2(snt)*power3(
     invdmsntau2mu) + 96*logmu2Q2*mst12*power2(mu2)*power2(snt)*power3(
     invdmsntau2mu) - 48*logmst12Q2*logmu2Q2*mst12*power2(mu2)*power2(snt)*
     power3(invdmsntau2mu) + 96*logmsntau2Q2*mst22*power2(mu2)*power2(snt)*
     power3(invdmsntau2mu) - 48*logmsntau2Q2*logmst22Q2*mst22*power2(mu2)*
     power2(snt)*power3(invdmsntau2mu) - 96*logmu2Q2*mst22*power2(mu2)*power2(
     snt)*power3(invdmsntau2mu) + 48*logmst22Q2*logmu2Q2*mst22*power2(mu2)*
     power2(snt)*power3(invdmsntau2mu) + 24*mst12*power2(logmsntau2Q2)*power2(
     mu2)*power2(snt)*power3(invdmsntau2mu) - 24*mst22*power2(logmsntau2Q2)*
     power2(mu2)*power2(snt)*power3(invdmsntau2mu) - 24*mst12*power2(logmu2Q2)*
     power2(mu2)*power2(snt)*power3(invdmsntau2mu) + 24*mst22*power2(logmu2Q2)*
     power2(mu2)*power2(snt)*power3(invdmsntau2mu) + 96*logmstau12Q2*msb12*
     power2(mu2)*power3(invdmstau1mu) - 48*logmsb12Q2*logmstau12Q2*msb12*power2
     (mu2)*power3(invdmstau1mu) - 96*logmu2Q2*msb12*power2(mu2)*power3(
     invdmstau1mu) + 48*logmsb12Q2*logmu2Q2*msb12*power2(mu2)*power3(
     invdmstau1mu) + 96*logmstau12Q2*msb22*power2(mu2)*power3(invdmstau1mu) -
     48*logmsb22Q2*logmstau12Q2*msb22*power2(mu2)*power3(invdmstau1mu) - 96*
     logmu2Q2*msb22*power2(mu2)*power3(invdmstau1mu) + 48*logmsb22Q2*logmu2Q2*
     msb22*power2(mu2)*power3(invdmstau1mu) - 24*msb12*power2(logmstau12Q2)*
     power2(mu2)*power3(invdmstau1mu) - 24*msb22*power2(logmstau12Q2)*power2(
     mu2)*power3(invdmstau1mu) + 24*msb12*power2(logmu2Q2)*power2(mu2)*power3(
     invdmstau1mu) + 24*msb22*power2(logmu2Q2)*power2(mu2)*power3(invdmstau1mu)
     + 96*logmstau22Q2*msb12*power2(mu2)*power3(invdmstau2mu) - 48*logmsb12Q2*
     logmstau22Q2*msb12*power2(mu2)*power3(invdmstau2mu) - 96*logmu2Q2*msb12*
     power2(mu2)*power3(invdmstau2mu) + 48*logmsb12Q2*logmu2Q2*msb12*power2(mu2
     )*power3(invdmstau2mu) + 96*logmstau22Q2*msb22*power2(mu2)*power3(
     invdmstau2mu) - 48*logmsb22Q2*logmstau22Q2*msb22*power2(mu2)*power3(
     invdmstau2mu) - 96*logmu2Q2*msb22*power2(mu2)*power3(invdmstau2mu) + 48*
     logmsb22Q2*logmu2Q2*msb22*power2(mu2)*power3(invdmstau2mu) - 24*msb12*
     power2(logmstau22Q2)*power2(mu2)*power3(invdmstau2mu) - 24*msb22*power2(
     logmstau22Q2)*power2(mu2)*power3(invdmstau2mu) + 24*msb12*power2(logmu2Q2)
     *power2(mu2)*power3(invdmstau2mu) + 24*msb22*power2(logmu2Q2)*power2(mu2)*
     power3(invdmstau2mu) + 42*invdmstau1mu*power2(invdmsb1stau1)*power3(msb12)
     - 72*invdmstau1mu*logmstau12Q2*power2(invdmsb1stau1)*power3(msb12) + 24*
     invdmstau1mu*logmsb12Q2*logmstau12Q2*power2(invdmsb1stau1)*power3(msb12) +
     42*invdmstau2mu*power2(invdmsb1stau2)*power3(msb12) - 72*invdmstau2mu*
     logmstau22Q2*power2(invdmsb1stau2)*power3(msb12) + 24*invdmstau2mu*
     logmsb12Q2*logmstau22Q2*power2(invdmsb1stau2)*power3(msb12) + 21*
     invdmsb1stau1*power2(invdmstau1mu)*power3(msb12) - 36*invdmsb1stau1*
     logmsb12Q2*power2(invdmstau1mu)*power3(msb12) + 72*invdmsb1stau1*
     logmstau12Q2*power2(invdmstau1mu)*power3(msb12) - 24*invdmsb1stau1*
     logmsb12Q2*logmstau12Q2*power2(invdmstau1mu)*power3(msb12) - 42*mu2*power2
     (invdmsb1stau1)*power2(invdmstau1mu)*power3(msb12) + 72*logmstau12Q2*mu2*
     power2(invdmsb1stau1)*power2(invdmstau1mu)*power3(msb12) - 24*logmsb12Q2*
     logmstau12Q2*mu2*power2(invdmsb1stau1)*power2(invdmstau1mu)*power3(msb12)
     + 21*invdmsb1stau2*power2(invdmstau2mu)*power3(msb12) - 36*invdmsb1stau2*
     logmsb12Q2*power2(invdmstau2mu)*power3(msb12) + 72*invdmsb1stau2*
     logmstau22Q2*power2(invdmstau2mu)*power3(msb12) - 24*invdmsb1stau2*
     logmsb12Q2*logmstau22Q2*power2(invdmstau2mu)*power3(msb12) - 42*mu2*power2
     (invdmsb1stau2)*power2(invdmstau2mu)*power3(msb12) + 72*logmstau22Q2*mu2*
     power2(invdmsb1stau2)*power2(invdmstau2mu)*power3(msb12) - 24*logmsb12Q2*
     logmstau22Q2*mu2*power2(invdmsb1stau2)*power2(invdmstau2mu)*power3(msb12)
     + 12*invdmstau1mu*power2(invdmsb1stau1)*power2(logmsb12Q2)*power3(msb12) +
     12*invdmstau2mu*power2(invdmsb1stau2)*power2(logmsb12Q2)*power3(msb12) -
     12*invdmsb1stau1*power2(invdmstau1mu)*power2(logmsb12Q2)*power3(msb12) -
     12*mu2*power2(invdmsb1stau1)*power2(invdmstau1mu)*power2(logmsb12Q2)*
     power3(msb12) - 12*invdmsb1stau2*power2(invdmstau2mu)*power2(logmsb12Q2)*
     power3(msb12) - 12*mu2*power2(invdmsb1stau2)*power2(invdmstau2mu)*power2(
     logmsb12Q2)*power3(msb12) + 12*invdmstau1mu*power2(invdmsb1stau1)*power2(
     logmstau12Q2)*power3(msb12) - 12*invdmsb1stau1*power2(invdmstau1mu)*power2
     (logmstau12Q2)*power3(msb12) - 12*mu2*power2(invdmsb1stau1)*power2(
     invdmstau1mu)*power2(logmstau12Q2)*power3(msb12) + 12*invdmstau2mu*power2(
     invdmsb1stau2)*power2(logmstau22Q2)*power3(msb12) - 12*invdmsb1stau2*
     power2(invdmstau2mu)*power2(logmstau22Q2)*power3(msb12) - 12*mu2*power2(
     invdmsb1stau2)*power2(invdmstau2mu)*power2(logmstau22Q2)*power3(msb12) + 4
     *invdmstau1mu*power2(invdmsb1stau1)*power2(Pi)*power3(msb12) + 4*
     invdmstau2mu*power2(invdmsb1stau2)*power2(Pi)*power3(msb12) - 4*
     invdmsb1stau1*power2(invdmstau1mu)*power2(Pi)*power3(msb12) - 4*mu2*power2
     (invdmsb1stau1)*power2(invdmstau1mu)*power2(Pi)*power3(msb12) - 4*
     invdmsb1stau2*power2(invdmstau2mu)*power2(Pi)*power3(msb12) - 4*mu2*power2
     (invdmsb1stau2)*power2(invdmstau2mu)*power2(Pi)*power3(msb12) + 42*
     invdmstau1mu*power2(invdmsb2stau1)*power3(msb22) - 72*invdmstau1mu*
     logmstau12Q2*power2(invdmsb2stau1)*power3(msb22) + 24*invdmstau1mu*
     logmsb22Q2*logmstau12Q2*power2(invdmsb2stau1)*power3(msb22) + 42*
     invdmstau2mu*power2(invdmsb2stau2)*power3(msb22) - 72*invdmstau2mu*
     logmstau22Q2*power2(invdmsb2stau2)*power3(msb22) + 24*invdmstau2mu*
     logmsb22Q2*logmstau22Q2*power2(invdmsb2stau2)*power3(msb22) - 21*
     invdmsb2stau1*power2(invdmstau1mu)*power3(msb22) + 36*invdmsb2stau1*
     logmsb22Q2*power2(invdmstau1mu)*power3(msb22) - 72*invdmsb2stau1*
     logmstau12Q2*power2(invdmstau1mu)*power3(msb22) + 24*invdmsb2stau1*
     logmsb22Q2*logmstau12Q2*power2(invdmstau1mu)*power3(msb22) - 42*mu2*power2
     (invdmsb2stau1)*power2(invdmstau1mu)*power3(msb22) + 72*logmstau12Q2*mu2*
     power2(invdmsb2stau1)*power2(invdmstau1mu)*power3(msb22) - 24*logmsb22Q2*
     logmstau12Q2*mu2*power2(invdmsb2stau1)*power2(invdmstau1mu)*power3(msb22)
     - 21*invdmsb2stau2*power2(invdmstau2mu)*power3(msb22) + 36*invdmsb2stau2*
     logmsb22Q2*power2(invdmstau2mu)*power3(msb22) - 72*invdmsb2stau2*
     logmstau22Q2*power2(invdmstau2mu)*power3(msb22) + 24*invdmsb2stau2*
     logmsb22Q2*logmstau22Q2*power2(invdmstau2mu)*power3(msb22) - 42*mu2*power2
     (invdmsb2stau2)*power2(invdmstau2mu)*power3(msb22) + 72*logmstau22Q2*mu2*
     power2(invdmsb2stau2)*power2(invdmstau2mu)*power3(msb22) - 24*logmsb22Q2*
     logmstau22Q2*mu2*power2(invdmsb2stau2)*power2(invdmstau2mu)*power3(msb22)
     + 12*invdmstau1mu*power2(invdmsb2stau1)*power2(logmsb22Q2)*power3(msb22) +
     12*invdmstau2mu*power2(invdmsb2stau2)*power2(logmsb22Q2)*power3(msb22) +
     12*invdmsb2stau1*power2(invdmstau1mu)*power2(logmsb22Q2)*power3(msb22) -
     12*mu2*power2(invdmsb2stau1)*power2(invdmstau1mu)*power2(logmsb22Q2)*
     power3(msb22) + 12*invdmsb2stau2*power2(invdmstau2mu)*power2(logmsb22Q2)*
     power3(msb22) - 12*mu2*power2(invdmsb2stau2)*power2(invdmstau2mu)*power2(
     logmsb22Q2)*power3(msb22) + 12*invdmstau1mu*power2(invdmsb2stau1)*power2(
     logmstau12Q2)*power3(msb22) + 12*invdmsb2stau1*power2(invdmstau1mu)*power2
     (logmstau12Q2)*power3(msb22) - 12*mu2*power2(invdmsb2stau1)*power2(
     invdmstau1mu)*power2(logmstau12Q2)*power3(msb22) + 12*invdmstau2mu*power2(
     invdmsb2stau2)*power2(logmstau22Q2)*power3(msb22) + 12*invdmsb2stau2*
     power2(invdmstau2mu)*power2(logmstau22Q2)*power3(msb22) - 12*mu2*power2(
     invdmsb2stau2)*power2(invdmstau2mu)*power2(logmstau22Q2)*power3(msb22) + 4
     *invdmstau1mu*power2(invdmsb2stau1)*power2(Pi)*power3(msb22) + 4*
     invdmstau2mu*power2(invdmsb2stau2)*power2(Pi)*power3(msb22) + 4*
     invdmsb2stau1*power2(invdmstau1mu)*power2(Pi)*power3(msb22) - 4*mu2*power2
     (invdmsb2stau1)*power2(invdmstau1mu)*power2(Pi)*power3(msb22) + 4*
     invdmsb2stau2*power2(invdmstau2mu)*power2(Pi)*power3(msb22) - 4*mu2*power2
     (invdmsb2stau2)*power2(invdmstau2mu)*power2(Pi)*power3(msb22) + 21*
     invdmsntaust1*power2(invdmsntau2mu)*power3(mst12) + 72*invdmsntaust1*
     logmsntau2Q2*power2(invdmsntau2mu)*power3(mst12) - 36*invdmsntaust1*
     logmst12Q2*power2(invdmsntau2mu)*power3(mst12) - 24*invdmsntaust1*
     logmsntau2Q2*logmst12Q2*power2(invdmsntau2mu)*power3(mst12) + 42*
     invdmsntau2mu*power2(invdmsntaust1)*power3(mst12) - 72*invdmsntau2mu*
     logmsntau2Q2*power2(invdmsntaust1)*power3(mst12) + 24*invdmsntau2mu*
     logmsntau2Q2*logmst12Q2*power2(invdmsntaust1)*power3(mst12) - 42*mu2*
     power2(invdmsntau2mu)*power2(invdmsntaust1)*power3(mst12) + 72*
     logmsntau2Q2*mu2*power2(invdmsntau2mu)*power2(invdmsntaust1)*power3(mst12)
     - 24*logmsntau2Q2*logmst12Q2*mu2*power2(invdmsntau2mu)*power2(
     invdmsntaust1)*power3(mst12) - 12*invdmsntaust1*power2(invdmsntau2mu)*
     power2(logmsntau2Q2)*power3(mst12) + 12*invdmsntau2mu*power2(invdmsntaust1
     )*power2(logmsntau2Q2)*power3(mst12) - 12*mu2*power2(invdmsntau2mu)*power2
     (invdmsntaust1)*power2(logmsntau2Q2)*power3(mst12) - 12*invdmsntaust1*
     power2(invdmsntau2mu)*power2(logmst12Q2)*power3(mst12) + 12*invdmsntau2mu*
     power2(invdmsntaust1)*power2(logmst12Q2)*power3(mst12) - 12*mu2*power2(
     invdmsntau2mu)*power2(invdmsntaust1)*power2(logmst12Q2)*power3(mst12) - 4*
     invdmsntaust1*power2(invdmsntau2mu)*power2(Pi)*power3(mst12) + 4*
     invdmsntau2mu*power2(invdmsntaust1)*power2(Pi)*power3(mst12) - 4*mu2*
     power2(invdmsntau2mu)*power2(invdmsntaust1)*power2(Pi)*power3(mst12) - 21*
     invdmsntaust1*power2(invdmsntau2mu)*power2(snt)*power3(mst12) - 72*
     invdmsntaust1*logmsntau2Q2*power2(invdmsntau2mu)*power2(snt)*power3(mst12)
     + 36*invdmsntaust1*logmst12Q2*power2(invdmsntau2mu)*power2(snt)*power3(
     mst12) + 24*invdmsntaust1*logmsntau2Q2*logmst12Q2*power2(invdmsntau2mu)*
     power2(snt)*power3(mst12) - 42*invdmsntau2mu*power2(invdmsntaust1)*power2(
     snt)*power3(mst12) + 72*invdmsntau2mu*logmsntau2Q2*power2(invdmsntaust1)*
     power2(snt)*power3(mst12) - 24*invdmsntau2mu*logmsntau2Q2*logmst12Q2*
     power2(invdmsntaust1)*power2(snt)*power3(mst12) + 42*mu2*power2(
     invdmsntau2mu)*power2(invdmsntaust1)*power2(snt)*power3(mst12) - 72*
     logmsntau2Q2*mu2*power2(invdmsntau2mu)*power2(invdmsntaust1)*power2(snt)*
     power3(mst12) + 24*logmsntau2Q2*logmst12Q2*mu2*power2(invdmsntau2mu)*
     power2(invdmsntaust1)*power2(snt)*power3(mst12) + 12*invdmsntaust1*power2(
     invdmsntau2mu)*power2(logmsntau2Q2)*power2(snt)*power3(mst12) - 12*
     invdmsntau2mu*power2(invdmsntaust1)*power2(logmsntau2Q2)*power2(snt)*
     power3(mst12) + 12*mu2*power2(invdmsntau2mu)*power2(invdmsntaust1)*power2(
     logmsntau2Q2)*power2(snt)*power3(mst12) + 12*invdmsntaust1*power2(
     invdmsntau2mu)*power2(logmst12Q2)*power2(snt)*power3(mst12) - 12*
     invdmsntau2mu*power2(invdmsntaust1)*power2(logmst12Q2)*power2(snt)*power3(
     mst12) + 12*mu2*power2(invdmsntau2mu)*power2(invdmsntaust1)*power2(
     logmst12Q2)*power2(snt)*power3(mst12) + 4*invdmsntaust1*power2(
     invdmsntau2mu)*power2(Pi)*power2(snt)*power3(mst12) - 4*invdmsntau2mu*
     power2(invdmsntaust1)*power2(Pi)*power2(snt)*power3(mst12) + 4*mu2*power2(
     invdmsntau2mu)*power2(invdmsntaust1)*power2(Pi)*power2(snt)*power3(mst12)
     + 21*invdmsntaust2*power2(invdmsntau2mu)*power2(snt)*power3(mst22) + 72*
     invdmsntaust2*logmsntau2Q2*power2(invdmsntau2mu)*power2(snt)*power3(mst22)
     - 36*invdmsntaust2*logmst22Q2*power2(invdmsntau2mu)*power2(snt)*power3(
     mst22) - 24*invdmsntaust2*logmsntau2Q2*logmst22Q2*power2(invdmsntau2mu)*
     power2(snt)*power3(mst22) + 42*invdmsntau2mu*power2(invdmsntaust2)*power2(
     snt)*power3(mst22) - 72*invdmsntau2mu*logmsntau2Q2*power2(invdmsntaust2)*
     power2(snt)*power3(mst22) + 24*invdmsntau2mu*logmsntau2Q2*logmst22Q2*
     power2(invdmsntaust2)*power2(snt)*power3(mst22) - 42*mu2*power2(
     invdmsntau2mu)*power2(invdmsntaust2)*power2(snt)*power3(mst22) + 72*
     logmsntau2Q2*mu2*power2(invdmsntau2mu)*power2(invdmsntaust2)*power2(snt)*
     power3(mst22) - 24*logmsntau2Q2*logmst22Q2*mu2*power2(invdmsntau2mu)*
     power2(invdmsntaust2)*power2(snt)*power3(mst22) - 12*invdmsntaust2*power2(
     invdmsntau2mu)*power2(logmsntau2Q2)*power2(snt)*power3(mst22) + 12*
     invdmsntau2mu*power2(invdmsntaust2)*power2(logmsntau2Q2)*power2(snt)*
     power3(mst22) - 12*mu2*power2(invdmsntau2mu)*power2(invdmsntaust2)*power2(
     logmsntau2Q2)*power2(snt)*power3(mst22) - 12*invdmsntaust2*power2(
     invdmsntau2mu)*power2(logmst22Q2)*power2(snt)*power3(mst22) + 12*
     invdmsntau2mu*power2(invdmsntaust2)*power2(logmst22Q2)*power2(snt)*power3(
     mst22) - 12*mu2*power2(invdmsntau2mu)*power2(invdmsntaust2)*power2(
     logmst22Q2)*power2(snt)*power3(mst22) - 4*invdmsntaust2*power2(
     invdmsntau2mu)*power2(Pi)*power2(snt)*power3(mst22) + 4*invdmsntau2mu*
     power2(invdmsntaust2)*power2(Pi)*power2(snt)*power3(mst22) - 4*mu2*power2(
     invdmsntau2mu)*power2(invdmsntaust2)*power2(Pi)*power2(snt)*power3(mst22)
     - 96*logmsntau2Q2*power3(invdmsntau2mu)*power3(mu2) + 96*logmu2Q2*power3(
     invdmsntau2mu)*power3(mu2) + 48*power2(logmsntau2Q2)*power3(invdmsntau2mu)
     *power3(mu2) - 48*power2(logmu2Q2)*power3(invdmsntau2mu)*power3(mu2) - 96*
     logmstau12Q2*power3(invdmstau1mu)*power3(mu2) + 96*logmu2Q2*power3(
     invdmstau1mu)*power3(mu2) + 48*power2(logmstau12Q2)*power3(invdmstau1mu)*
     power3(mu2) - 48*power2(logmu2Q2)*power3(invdmstau1mu)*power3(mu2) - 96*
     logmstau22Q2*power3(invdmstau2mu)*power3(mu2) + 96*logmu2Q2*power3(
     invdmstau2mu)*power3(mu2) + 48*power2(logmstau22Q2)*power3(invdmstau2mu)*
     power3(mu2) - 48*power2(logmu2Q2)*power3(invdmstau2mu)*power3(mu2) + 42*
     power2(invdmsb1stau1)*power2(invdmstau1mu)*power4(msb12) - 72*logmstau12Q2
     *power2(invdmsb1stau1)*power2(invdmstau1mu)*power4(msb12) + 24*logmsb12Q2*
     logmstau12Q2*power2(invdmsb1stau1)*power2(invdmstau1mu)*power4(msb12) + 42
     *power2(invdmsb1stau2)*power2(invdmstau2mu)*power4(msb12) - 72*
     logmstau22Q2*power2(invdmsb1stau2)*power2(invdmstau2mu)*power4(msb12) + 24
     *logmsb12Q2*logmstau22Q2*power2(invdmsb1stau2)*power2(invdmstau2mu)*power4
     (msb12) + 12*power2(invdmsb1stau1)*power2(invdmstau1mu)*power2(logmsb12Q2)
     *power4(msb12) + 12*power2(invdmsb1stau2)*power2(invdmstau2mu)*power2(
     logmsb12Q2)*power4(msb12) + 12*power2(invdmsb1stau1)*power2(invdmstau1mu)*
     power2(logmstau12Q2)*power4(msb12) + 12*power2(invdmsb1stau2)*power2(
     invdmstau2mu)*power2(logmstau22Q2)*power4(msb12) + 4*power2(invdmsb1stau1)
     *power2(invdmstau1mu)*power2(Pi)*power4(msb12) + 4*power2(invdmsb1stau2)*
     power2(invdmstau2mu)*power2(Pi)*power4(msb12) + 42*power2(invdmsb2stau1)*
     power2(invdmstau1mu)*power4(msb22) - 72*logmstau12Q2*power2(invdmsb2stau1)
     *power2(invdmstau1mu)*power4(msb22) + 24*logmsb22Q2*logmstau12Q2*power2(
     invdmsb2stau1)*power2(invdmstau1mu)*power4(msb22) + 42*power2(
     invdmsb2stau2)*power2(invdmstau2mu)*power4(msb22) - 72*logmstau22Q2*power2
     (invdmsb2stau2)*power2(invdmstau2mu)*power4(msb22) + 24*logmsb22Q2*
     logmstau22Q2*power2(invdmsb2stau2)*power2(invdmstau2mu)*power4(msb22) + 12
     *power2(invdmsb2stau1)*power2(invdmstau1mu)*power2(logmsb22Q2)*power4(
     msb22) + 12*power2(invdmsb2stau2)*power2(invdmstau2mu)*power2(logmsb22Q2)*
     power4(msb22) + 12*power2(invdmsb2stau1)*power2(invdmstau1mu)*power2(
     logmstau12Q2)*power4(msb22) + 12*power2(invdmsb2stau2)*power2(invdmstau2mu
     )*power2(logmstau22Q2)*power4(msb22) + 4*power2(invdmsb2stau1)*power2(
     invdmstau1mu)*power2(Pi)*power4(msb22) + 4*power2(invdmsb2stau2)*power2(
     invdmstau2mu)*power2(Pi)*power4(msb22) + 42*power2(invdmsntau2mu)*power2(
     invdmsntaust1)*power4(mst12) - 72*logmsntau2Q2*power2(invdmsntau2mu)*
     power2(invdmsntaust1)*power4(mst12) + 24*logmsntau2Q2*logmst12Q2*power2(
     invdmsntau2mu)*power2(invdmsntaust1)*power4(mst12) + 12*power2(
     invdmsntau2mu)*power2(invdmsntaust1)*power2(logmsntau2Q2)*power4(mst12) +
     12*power2(invdmsntau2mu)*power2(invdmsntaust1)*power2(logmst12Q2)*power4(
     mst12) + 4*power2(invdmsntau2mu)*power2(invdmsntaust1)*power2(Pi)*power4(
     mst12) - 42*power2(invdmsntau2mu)*power2(invdmsntaust1)*power2(snt)*power4
     (mst12) + 72*logmsntau2Q2*power2(invdmsntau2mu)*power2(invdmsntaust1)*
     power2(snt)*power4(mst12) - 24*logmsntau2Q2*logmst12Q2*power2(
     invdmsntau2mu)*power2(invdmsntaust1)*power2(snt)*power4(mst12) - 12*power2
     (invdmsntau2mu)*power2(invdmsntaust1)*power2(logmsntau2Q2)*power2(snt)*
     power4(mst12) - 12*power2(invdmsntau2mu)*power2(invdmsntaust1)*power2(
     logmst12Q2)*power2(snt)*power4(mst12) - 4*power2(invdmsntau2mu)*power2(
     invdmsntaust1)*power2(Pi)*power2(snt)*power4(mst12) + 42*power2(
     invdmsntau2mu)*power2(invdmsntaust2)*power2(snt)*power4(mst22) - 72*
     logmsntau2Q2*power2(invdmsntau2mu)*power2(invdmsntaust2)*power2(snt)*
     power4(mst22) + 24*logmsntau2Q2*logmst22Q2*power2(invdmsntau2mu)*power2(
     invdmsntaust2)*power2(snt)*power4(mst22) + 12*power2(invdmsntau2mu)*power2
     (invdmsntaust2)*power2(logmsntau2Q2)*power2(snt)*power4(mst22) + 12*power2
     (invdmsntau2mu)*power2(invdmsntaust2)*power2(logmst22Q2)*power2(snt)*
     power4(mst22) + 4*power2(invdmsntau2mu)*power2(invdmsntaust2)*power2(Pi)*
     power2(snt)*power4(mst22) + 576*power4(sa) - 432*logmh2Q2*power4(sa) + 144
     *logmH2Q2*power4(sa) - 504*invdmhH*msb12*power4(sa) + 240*invdmhH*logmH2Q2
     *msb12*power4(sa) + 288*invdmhH*logmsb12Q2*msb12*power4(sa) - 504*invdmhH*
     msb22*power4(sa) + 240*invdmhH*logmH2Q2*msb22*power4(sa) + 288*invdmhH*
     logmsb22Q2*msb22*power4(sa) - 324*invdmhH*mu2*power4(sa) + 576*invdmhH*
     logmh2Q2*mu2*power4(sa) - 240*invdmhH*logmH2Q2*mu2*power4(sa) + 240*
     invdmhH*logmsb12Q2*mu2*power4(sa) - 288*invdmhH*logmh2Q2*logmsb12Q2*mu2*
     power4(sa) + 240*invdmhH*invdmsb*logmsb12Q2*msb22*mu2*power4(sa) - 288*
     invdmhH*invdmsb*logmh2Q2*logmsb12Q2*msb22*mu2*power4(sa) - 240*invdmhH*
     invdmsb*logmsb22Q2*msb22*mu2*power4(sa) + 288*invdmhH*invdmsb*logmh2Q2*
     logmsb22Q2*msb22*mu2*power4(sa) + 72*power2(logmh2Q2)*power4(sa) - 144*
     invdmhH*mu2*power2(logmh2Q2)*power4(sa) - 72*power2(logmH2Q2)*power4(sa) -
     144*invdmhH*msb12*power2(logmH2Q2)*power4(sa) - 144*invdmhH*msb22*power2(
     logmH2Q2)*power4(sa) + 144*invdmhH*mu2*power2(logmH2Q2)*power4(sa) - 144*
     invdmhH*msb12*power2(logmsb12Q2)*power4(sa) - 144*invdmhH*msb22*power2(
     logmsb22Q2)*power4(sa) - 48*invdmhH*msb12*power2(Pi)*power4(sa) - 48*
     invdmhH*msb22*power2(Pi)*power4(sa)))/mH2)/64. + (3*invdmstau1mu*Fin20(
     mstau12,msb12,Q2)*(1 + invdmstau1mu*mu2 + invdmsb1stau1*(-2*invdmstau1mu*
     msb14 + msb12*(-1 + invdmstau1mu*mu2)) + 2*(msb12 - mu2)*mu2*power2(
     invdmstau1mu) + power2(invdmsb1stau1)*(msb14 + invdmstau1mu*(-(msb14*mu2)
     + power6(msb1)))))/8. + (3*invdmstau2mu*Fin20(mstau22,msb12,Q2)*(1 +
     invdmstau2mu*mu2 + invdmsb1stau2*(-2*invdmstau2mu*msb14 + msb12*(-1 +
     invdmstau2mu*mu2)) + 2*(msb12 - mu2)*mu2*power2(invdmstau2mu) + power2(
     invdmsb1stau2)*(msb14 + invdmstau2mu*(-(msb14*mu2) + power6(msb1)))))/8. +
     (mu2*DeltaInv(msb22,mt2,mu2)*power2(invdmsntau2mu)*(msb24*mt2*(42 + 6*
     logmsb22Q2*(-3 + logmt2Q2) - 18*logmt2Q2 + 3*power2(logmsb22Q2) + 3*power2
     (logmt2Q2) + power2(Pi)) + msb24*mu2*(42 + 12*logmsb22Q2*(-3 + logmt2Q2) +
     18*logmu2Q2 - 6*logmt2Q2*(3 + logmu2Q2) + 6*power2(logmsb22Q2) + 3*power2(
     logmt2Q2) - 3*power2(logmu2Q2) + power2(Pi)) + (mt2 - mu2)*power2(mu2)*(42
      + 6*logmt2Q2*(-3 + logmu2Q2) - 18*logmu2Q2 + 3*power2(logmt2Q2) + 3*
     power2(logmu2Q2) + power2(Pi)) + msb22*mu2*(mu2*(42 - 6*logmsb22Q2*(-3 +
     logmt2Q2) - 36*logmu2Q2 + 6*logmt2Q2*(-3 + 2*logmu2Q2) - 3*power2(
     logmsb22Q2) + 3*power2(logmt2Q2) + 6*power2(logmu2Q2) + power2(Pi)) + mt2*
     (168 + 6*logmsb22Q2*(-3 + 3*logmt2Q2 - 2*logmu2Q2) + 18*logmt2Q2*(-6 +
     logmu2Q2) - 18*logmu2Q2 + 3*power2(logmsb22Q2) + 18*power2(logmt2Q2) + 3*
     power2(logmu2Q2) + 4*power2(Pi))) - (42 + 6*logmsb22Q2*(-3 + logmt2Q2) -
     18*logmt2Q2 + 3*power2(logmsb22Q2) + 3*power2(logmt2Q2) + power2(Pi))*
     power6(msb2)))/8. + (DeltaInv(msb22,msntau2,mt2)*(42*msb22*msntau2 + 18*
     logmsb22Q2*msb22*msntau2 - 36*logmsntau2Q2*msb22*msntau2 - 18*logmt2Q2*
     msb22*msntau2 - 6*logmsb22Q2*logmt2Q2*msb22*msntau2 + 12*logmsntau2Q2*
     logmt2Q2*msb22*msntau2 - 42*msntau4 + 18*logmsntau2Q2*msntau4 + 18*
     logmt2Q2*msntau4 - 6*logmsntau2Q2*logmt2Q2*msntau4 + 168*msb22*mt2 - 18*
     logmsb22Q2*msb22*mt2 - 18*logmsntau2Q2*msb22*mt2 - 12*logmsb22Q2*
     logmsntau2Q2*msb22*mt2 - 108*logmt2Q2*msb22*mt2 + 18*logmsb22Q2*logmt2Q2*
     msb22*mt2 + 18*logmsntau2Q2*logmt2Q2*msb22*mt2 + 42*msntau2*mt2 - 18*
     logmsntau2Q2*msntau2*mt2 - 18*logmt2Q2*msntau2*mt2 + 6*logmsntau2Q2*
     logmt2Q2*msntau2*mt2 + 84*msb22*mu2 + 36*logmsb22Q2*msb22*mu2 - 72*
     logmsntau2Q2*msb22*mu2 - 36*logmt2Q2*msb22*mu2 - 12*logmsb22Q2*logmt2Q2*
     msb22*mu2 + 24*logmsntau2Q2*logmt2Q2*msb22*mu2 - 84*msntau2*mu2 + 36*
     logmsntau2Q2*msntau2*mu2 + 36*logmt2Q2*msntau2*mu2 - 12*logmsntau2Q2*
     logmt2Q2*msntau2*mu2 + 84*mt2*mu2 - 36*logmsntau2Q2*mt2*mu2 - 36*logmt2Q2*
     mt2*mu2 + 12*logmsntau2Q2*logmt2Q2*mt2*mu2 - 3*msb22*msntau2*power2(
     logmsb22Q2) + 3*msb22*mt2*power2(logmsb22Q2) - 6*msb22*mu2*power2(
     logmsb22Q2) + 6*msb22*msntau2*power2(logmsntau2Q2) - 3*msntau4*power2(
     logmsntau2Q2) + 3*msb22*mt2*power2(logmsntau2Q2) + 3*msntau2*mt2*power2(
     logmsntau2Q2) + 12*msb22*mu2*power2(logmsntau2Q2) - 6*msntau2*mu2*power2(
     logmsntau2Q2) + 6*mt2*mu2*power2(logmsntau2Q2) + 3*msb22*msntau2*power2(
     logmt2Q2) - 3*msntau4*power2(logmt2Q2) + 18*msb22*mt2*power2(logmt2Q2) + 3
     *msntau2*mt2*power2(logmt2Q2) + 6*msb22*mu2*power2(logmt2Q2) - 6*msntau2*
     mu2*power2(logmt2Q2) + 6*mt2*mu2*power2(logmt2Q2) - 126*power2(mu2) + 54*
     logmsntau2Q2*power2(mu2) + 54*logmt2Q2*power2(mu2) - 18*logmsntau2Q2*
     logmt2Q2*power2(mu2) - 9*power2(logmsntau2Q2)*power2(mu2) - 9*power2(
     logmt2Q2)*power2(mu2) + msb22*msntau2*power2(Pi) - msntau4*power2(Pi) + 4*
     msb22*mt2*power2(Pi) + msntau2*mt2*power2(Pi) + 2*msb22*mu2*power2(Pi) - 2
     *msntau2*mu2*power2(Pi) + 2*mt2*mu2*power2(Pi) - 3*power2(mu2)*power2(Pi)
     + msb24*(42 + 12*logmsb22Q2*(-3 + logmt2Q2) - 6*logmsntau2Q2*(-3 +
     logmt2Q2) - 18*logmt2Q2 + 6*power2(logmsb22Q2) - 3*power2(logmsntau2Q2) +
     3*power2(logmt2Q2) + power2(Pi)) + mu2*power2(invdmsntau2mu)*(msb24*mt2*(
     42 + 6*logmsb22Q2*(-3 + logmt2Q2) - 18*logmt2Q2 + 3*power2(logmsb22Q2) + 3
     *power2(logmt2Q2) + power2(Pi)) + msb24*mu2*(42 + 12*logmsb22Q2*(-3 +
     logmt2Q2) - 6*logmsntau2Q2*(-3 + logmt2Q2) - 18*logmt2Q2 + 6*power2(
     logmsb22Q2) - 3*power2(logmsntau2Q2) + 3*power2(logmt2Q2) + power2(Pi)) +
     (mt2 - mu2)*power2(mu2)*(42 + 6*logmsntau2Q2*(-3 + logmt2Q2) - 18*logmt2Q2
      + 3*power2(logmsntau2Q2) + 3*power2(logmt2Q2) + power2(Pi)) + msb22*mu2*(
     mu2*(42 - 6*logmsb22Q2*(-3 + logmt2Q2) + 12*logmsntau2Q2*(-3 + logmt2Q2) -
     18*logmt2Q2 - 3*power2(logmsb22Q2) + 6*power2(logmsntau2Q2) + 3*power2(
     logmt2Q2) + power2(Pi)) + mt2*(-6*logmsb22Q2*(3 + 2*logmsntau2Q2 - 3*
     logmt2Q2) + 18*logmsntau2Q2*(-1 + logmt2Q2) + 3*power2(logmsb22Q2) + 3*
     power2(logmsntau2Q2) + 2*(84 - 54*logmt2Q2 + 9*power2(logmt2Q2) + 2*power2
     (Pi)))) - (42 + 6*logmsb22Q2*(-3 + logmt2Q2) - 18*logmt2Q2 + 3*power2(
     logmsb22Q2) + 3*power2(logmt2Q2) + power2(Pi))*power6(msb2)) +
     invdmsntau2mu*(-(msb24*(mt2*(42 + 6*logmsb22Q2*(-3 + logmt2Q2) - 18*
     logmt2Q2 + 3*power2(logmsb22Q2) + 3*power2(logmt2Q2) + power2(Pi)) + 2*mu2
     *(42 + 12*logmsb22Q2*(-3 + logmt2Q2) - 6*logmsntau2Q2*(-3 + logmt2Q2) - 18
     *logmt2Q2 + 6*power2(logmsb22Q2) - 3*power2(logmsntau2Q2) + 3*power2(
     logmt2Q2) + power2(Pi)))) + mu2*(-((3*mt2 - 4*mu2)*mu2*(42 + 6*
     logmsntau2Q2*(-3 + logmt2Q2) - 18*logmt2Q2 + 3*power2(logmsntau2Q2) + 3*
     power2(logmt2Q2) + power2(Pi))) - msb22*(3*mu2*(42 - 6*logmsb22Q2*(-3 +
     logmt2Q2) + 12*logmsntau2Q2*(-3 + logmt2Q2) - 18*logmt2Q2 - 3*power2(
     logmsb22Q2) + 6*power2(logmsntau2Q2) + 3*power2(logmt2Q2) + power2(Pi)) +
     mt2*(-12*logmsb22Q2*(3 + 2*logmsntau2Q2 - 3*logmt2Q2) + 36*logmsntau2Q2*(-
     1 + logmt2Q2) + 6*power2(logmsb22Q2) + 6*power2(logmsntau2Q2) + 4*(84 - 54
     *logmt2Q2 + 9*power2(logmt2Q2) + 2*power2(Pi))))) + (42 + 6*logmsb22Q2*(-3
      + logmt2Q2) - 18*logmt2Q2 + 3*power2(logmsb22Q2) + 3*power2(logmt2Q2) +
     power2(Pi))*power6(msb2))))/8. + (3*invdmstau1mu*Fin20(msb22,mstau12,Q2)*(
     1 + invdmstau1mu*mu2 + invdmsb2stau1*(msb22 + 2*invdmstau1mu*msb24 -
     invdmstau1mu*msb22*mu2) + 2*(msb22 - mu2)*mu2*power2(invdmstau1mu) +
     power2(invdmsb2stau1)*(msb24 + invdmstau1mu*(-(msb24*mu2) + power6(msb2)))
     ))/8. + (3*invdmstau2mu*Fin20(mstau22,msb22,Q2)*(1 + invdmstau2mu*mu2 +
     invdmsb2stau2*(msb22 + 2*invdmstau2mu*msb24 - invdmstau2mu*msb22*mu2) + 2*
     (msb22 - mu2)*mu2*power2(invdmstau2mu) + power2(invdmsb2stau2)*(msb24 +
     invdmstau2mu*(-(msb24*mu2) + power6(msb2)))))/8. - (3*invdmsntau2mu*Fin20(
     mst12,msntau2,Q2)*(-1 + power2(snt))*(1 - invdmsntaust1*mst12 + 2*(mst12 -
     mu2)*mu2*power2(invdmsntau2mu) + mst14*power2(invdmsntaust1) +
     invdmsntau2mu*(mu2 + invdmsntaust1*(-2*mst14 + mst12*mu2) + power2(
     invdmsntaust1)*(-(mst14*mu2) + power6(mst1)))))/8. + (3*invdmsntau2mu*
     Fin20(msntau2,mst22,Q2)*power2(snt)*(1 - invdmsntaust2*mst22 + 2*(mst22 -
     mu2)*mu2*power2(invdmsntau2mu) + mst24*power2(invdmsntaust2) +
     invdmsntau2mu*(mu2 + invdmsntaust2*(-2*mst24 + mst22*mu2) + power2(
     invdmsntaust2)*(-(mst24*mu2) + power6(mst2)))))/8.;

   return result * power2(ytau) * power2(yb) * twoLoop;
}

std::ostream& operator<<(std::ostream& out, const Parameters& pars)
{
   out <<
      "Delta m_tau 2L parameters:\n"
      "yt     = " << pars.yt     << '\n' <<
      "yb     = " << pars.yb     << '\n' <<
      "ytau   = " << pars.ytau   << '\n' <<
      "mt     = " << pars.mt     << '\n' <<
      "mb     = " << pars.mb     << '\n' <<
      "mtau   = " << pars.mtau   << '\n' <<
      "mst1   = " << pars.mst1   << '\n' <<
      "mst2   = " << pars.mst2   << '\n' <<
      "msb1   = " << pars.msb1   << '\n' <<
      "msb2   = " << pars.msb2   << '\n' <<
      "mstau1 = " << pars.mstau1 << '\n' <<
      "mstau2 = " << pars.mstau2 << '\n' <<
      "msntau = " << pars.msntau << '\n' <<
      "xt     = " << pars.xt     << '\n' <<
      "xb     = " << pars.xb     << '\n' <<
      "xtau   = " << pars.xtau   << '\n' <<
      "mw     = " << pars.mw     << '\n' <<
      "mz     = " << pars.mz     << '\n' <<
      "mh     = " << pars.mh     << '\n' <<
      "mH     = " << pars.mH     << '\n' <<
      "mC     = " << pars.mC     << '\n' <<
      "mA     = " << pars.mA     << '\n' <<
      "mu     = " << pars.mu     << '\n' <<
      "tb     = " << pars.tb     << '\n' <<
      "Q      = " << pars.Q      << '\n';

   return out;
}

} // namespace mssm_twoloop_mtau
} // namespace flexiblesusy
