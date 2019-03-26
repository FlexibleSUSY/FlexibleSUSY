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

/** \file numerics.cpp
   - Project:     SOFTSUSY
   - Author:      Ben Allanach, Alexander Voigt
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

*/

/// Comment if you want default softsusy behaviour
// #define USE_LOOPTOOLS

#include "numerics.h"
#include "numerics2.hpp"
#ifdef USE_LOOPTOOLS
#include <clooptools.h>
#endif
#include <algorithm>
#include <cmath>

namespace softsusy {

namespace {

constexpr double EPSTOL = 1.0e-11; ///< underflow accuracy
constexpr double TOL = 1e-4;

constexpr double dabs(double a) noexcept { return a >= 0. ? a : -a; }
constexpr double sqr(double a) noexcept { return a*a; }
constexpr double pow3(double a) noexcept { return a*a*a; }
constexpr double pow6(double a) noexcept { return a*a*a*a*a*a; }

constexpr bool is_zero(double m, double tol) noexcept
{
   const double am = dabs(m);
   const double mtol = tol * am;

   if (mtol == 0.0 && am != 0.0 && tol != 0.0)
      return am <= tol;

   return am <= mtol;
}

constexpr bool is_close(double m1, double m2, double tol) noexcept
{
   const double mmax = std::max(dabs(m1), dabs(m2));
   const double mmin = std::min(dabs(m1), dabs(m2));
   const double max_tol = tol * mmax;

   if (max_tol == 0.0 && mmax != 0.0 && tol != 0.0)
      return mmax - mmin <= tol;

   return mmax - mmin <= max_tol;
}

/// returns a/b if a/b is finite, otherwise returns numeric_limits::max()
template <typename T>
constexpr T divide_finite(T a, T b) noexcept {
   T result = a / b;
   if (!std::isfinite(result))
      result = std::numeric_limits<T>::max();
   return result;
}

// can be made constexpr in C++20
double fB(const std::complex<double>& a) noexcept
{
   using flexiblesusy::fast_log;
   const double x = a.real();

   if (fabs(x) < EPSTOL)
      return -1. - x + sqr(x) * 0.5;

   if (is_close(x, 1., EPSTOL))
      return -1.;

   return std::real(fast_log(1. - a) - 1. - a * fast_log(1.0 - 1.0 / a));
}

} // anonymous namespace

double a0(double m, double q) noexcept {
   using std::fabs;
   using std::log;
   constexpr double TOL = 1e-4;
   if (fabs(m) < TOL) return 0.;
   return sqr(m) * (1.0 - 2. * log(fabs(m / q)));
}

double ffn(double p, double m1, double m2, double q) noexcept {
   return a0(m1, q) - 2.0 * a0(m2, q) -
      (2.0 * sqr(p) + 2.0 * sqr(m1) - sqr(m2)) *
      b0(p, m1, m2, q);
}

double gfn(double p, double m1, double m2, double q) noexcept {
   return (sqr(p) - sqr(m1) - sqr(m2)) * b0(p, m1, m2, q) - a0(m1, q)
      - a0(m2, q);
}

double hfn(double p, double m1, double m2, double q) noexcept {
   return 4.0 * b22(p, m1, m2, q) + gfn(p, m1, m2, q);
}

double b22bar(double p, double m1, double m2, double q) noexcept {
   return b22(p, m1, m2, q) - 0.25 * a0(m1, q) - 0.25 * a0(m2, q);
}

/**
 * Returns Re(B0(p,m1,m2,q)), from hep-ph/9606211
 */
double b0(double p, double m1, double m2, double q) noexcept
{
   using std::fabs;
   using std::log;
   using std::sqrt;

#ifdef USE_LOOPTOOLS
   setmudim(q*q);
   double b0l = B0(p*p, m1*m1, m2*m2).real();
   //  return B0(p*p, m1*m1, m2*m2).real();
#endif

   // protect against infrared divergence
   if (is_zero(p, EPSTOL) && is_zero(m1, EPSTOL) && is_zero(m2, EPSTOL))
      return 0.0;

   const double mMin = std::min(fabs(m1), fabs(m2));
   const double mMax = std::max(fabs(m1), fabs(m2));

   const double pSq = sqr(p), mMinSq = sqr(mMin), mMaxSq = sqr(mMax);
   /// Try to increase the accuracy of s
   const double dmSq = mMaxSq - mMinSq;
   const double s = pSq + dmSq;

   const double pTest = divide_finite(pSq, mMaxSq);
   /// Decides level at which one switches to p=0 limit of calculations
   const double pTolerance = 1.0e-10;

   /// p is not 0
   if (pTest > pTolerance) {
      const std::complex<double> ieps(0.0, EPSTOL * mMaxSq);
      const std::complex<double> x = s + sqrt(sqr(s) - 4. * pSq * (mMaxSq - ieps));
      const std::complex<double> xPlus = x / (2. * pSq);
      const std::complex<double> xMinus = 2. * (mMaxSq - ieps) / x;

      return -2.0 * log(p / q) - fB(xPlus) - fB(xMinus);
   }

   if (is_close(m1, m2, EPSTOL)) {
      return - log(sqr(m1 / q));
   }

   const double Mmax2 = mMaxSq, Mmin2 = mMinSq;

   if (Mmin2 < 1.e-30) {
      return 1.0 - log(Mmax2 / sqr(q));
   }

   return 1.0 - log(Mmax2 / sqr(q)) + Mmin2 * log(Mmax2 / Mmin2)
      / (Mmin2 - Mmax2);
}

/// Note that b1 is NOT symmetric in m1 <-> m2!!!
double b1(double p, double m1, double m2, double q) noexcept
{
   using std::fabs;
   using std::log;

#ifdef USE_LOOPTOOLS
   setmudim(q*q);
   double b1l = -B1(p*p, m1*m1, m2*m2).real();
   //    return b1l;
#endif

   // protect against infrared divergence
   if (is_zero(p, EPSTOL) && is_zero(m1, EPSTOL) && is_zero(m2, EPSTOL))
      return 0.0;

   const double p2 = sqr(p), m12 = sqr(m1), m22 = sqr(m2), q2 = sqr(q);
   const double pTest = divide_finite(p2, std::max(m12, m22));

   /// Decides level at which one switches to p=0 limit of calculations
   const double pTolerance = 1.0e-4;

   if (pTest > pTolerance) {
      return (a0(m2, q) - a0(m1, q) + (p2 + m12 - m22)
              * b0(p, m1, m2, q)) / (2.0 * p2);
   }

   if (fabs(m1) > 1.0e-15 && fabs(m2) > 1.0e-15) {
      const double m14 = sqr(m12), m24 = sqr(m22);
      const double m16 = m12*m14 , m26 = m22*m24;
      const double m18 = sqr(m14), m28 = sqr(m24);
      const double p4 = sqr(p2);

      if (fabs(m12 - m22) < pTolerance * std::max(m12, m22)) {
         return 0.08333333333333333*p2/m22
            + 0.008333333333333333*p4/m24
            + sqr(m12 - m22)*(0.041666666666666664/m24 +
                              0.016666666666666666*p2/m26 +
                              0.005357142857142856*p4/m28)
            + (m12 - m22)*(-0.16666666666666666/m22 -
                           0.03333333333333333*p2/m24 -
                           0.007142857142857142*p4/m26)
            - 0.5*log(m22/q2);
      }

      const double l12 = log(m12/m22);

      return (3*m14 - 4*m12*m22 + m24 - 2*m14*l12)/(4.*sqr(m12 - m22))
         + (p2*(4*pow3(m12 - m22)*
                (2*m14 + 5*m12*m22 - m24) +
                (3*m18 + 44*m16*m22 - 36*m14*m24 - 12*m12*m26 + m28)*p2
                - 12*m14*m22*(2*sqr(m12 - m22) + (2*m12 + 3*m22)*p2)*l12))/
         (24.*pow6(m12 - m22)) - 0.5*log(m22/q2);
   }

   return (m12 > m22)
      ? -0.5*log(m12/q2) + 0.75
      : -0.5*log(m22/q2) + 0.25;
}

double b22(double p,  double m1, double m2, double q) noexcept
{
   using std::fabs;
   using std::log;

#ifdef USE_LOOPTOOLS
   setmudim(q*q);
   double b22l = B00(p*p, m1*m1, m2*m2).real();
#endif

   // protect against infrared divergence
   if (is_zero(p, EPSTOL) && is_zero(m1, EPSTOL) && is_zero(m2, EPSTOL))
      return 0.0;

   /// Decides level at which one switches to p=0 limit of calculations
   const double p2 = sqr(p), m12 = sqr(m1), m22 = sqr(m2);
   const double pTolerance = 1.0e-10;

   if (p2 < pTolerance * std::max(m12, m22) ) {
      // m1 == m2 with good accuracy
      if (is_close(m1, m2, EPSTOL)) {
         return -m12 * log(sqr(m1 / q)) * 0.5 + m12 * 0.5;
      }
      // p == 0 limit
      if (fabs(m1) > EPSTOL && fabs(m2) > EPSTOL) {
         return 0.375 * (m12 + m22) - 0.25 *
            (sqr(m22) * log(sqr(m2 / q)) - sqr(m12) *
             log(sqr(m1 / q))) / (m22 - m12);
      }
      return (fabs(m1) < EPSTOL)
         ? 0.375 * m22 - 0.25 * m22 * log(sqr(m2 / q))
         : 0.375 * m12 - 0.25 * m12 * log(sqr(m1 / q));
   }

   const double b0Save = b0(p, m1, m2, q);
   const double a01 = a0(m1, q);
   const double a02 = a0(m2, q);

   return 1.0 / 6.0 *
      (0.5 * (a01 + a02) + (m12 + m22 - 0.5 * p2)
       * b0Save + (m22 - m12) / (2.0 * p2) *
       (a02 - a01 - (m22 - m12) * b0Save) +
       m12 + m22 - p2 / 3.0);
}

double d0(double m1, double m2, double m3, double m4) noexcept
{
   using std::log;

   const double m1sq = sqr(m1);
   const double m2sq = sqr(m2);

   if (is_close(m1, m2, EPSTOL)) {
      const double m3sq = sqr(m3);
      const double m4sq = sqr(m4);

      if (is_zero(m2, EPSTOL)) {
         // d0 is undefined for m1 == m2 == 0
         return 0.;
      } else if (is_zero(m3, EPSTOL)) {
         return (-m2sq + m4sq - m2sq * log(m4sq/m2sq))/
            sqr(m2 * m2sq - m2 * m4sq);
      } else if (is_zero(m4, EPSTOL)) {
         return (-m2sq + m3sq - m2sq * log(m3sq/m2sq))/
            sqr(m2 * m2sq - m2 * m3sq);
      } else if (is_close(m2, m3, EPSTOL) && is_close(m2, m4, EPSTOL)) {
         return 1.0 / (6.0 * sqr(m2sq));
      } else if (is_close(m2, m3, EPSTOL)) {
         return (sqr(m2sq) - sqr(m4sq) + 2.0 * m4sq * m2sq * log(m4sq / m2sq)) /
            (2.0 * m2sq * sqr(m2sq - m4sq) * (m2sq - m4sq));
      } else if (is_close(m2, m4, EPSTOL)) {
         return (sqr(m2sq) - sqr(m3sq) + 2.0 * m3sq * m2sq * log(m3sq / m2sq)) /
            (2.0 * m2sq * sqr(m2sq - m3sq) * (m2sq - m3sq));
      } else if (is_close(m3, m4, EPSTOL)) {
         return -1.0 / sqr(m2sq - m3sq) *
            ((m2sq + m3sq) / (m2sq - m3sq) * log(m3sq / m2sq) + 2.0);
      }

      return
         (m4sq / sqr(m2sq - m4sq) * log(m4sq / m2sq) +
          m4sq / (m2sq * (m2sq - m4sq)) -
          m3sq / sqr(m2sq - m3sq) * log(m3sq / m2sq) -
          m3sq / (m2sq * (m2sq - m3sq))) / (m3sq - m4sq);
   }
   return (c0(m1, m3, m4) - c0(m2, m3, m4)) / (m1sq - m2sq);
}

double d27(double m1, double m2, double m3, double m4) noexcept
{
   if (is_close(m1, m2, EPSTOL))
      m1 += TOL * 0.01;

   const double m12 = sqr(m1), m22 = sqr(m2);

   return (m12 * c0(m1, m3, m4) - m22 * c0(m2, m3, m4))
      / (4.0 * (m12 - m22));
}

double c0(double m1, double m2, double m3) noexcept
{
   using std::log;

#ifdef USE_LOOPTOOLS
   double q = 100.;
   setmudim(q*q);
   double psq = 0.;
   double c0l = C0(psq, psq, psq, m1*m1, m2*m2, m3*m3).real();
#endif

   const double m12 = sqr(m1), m22 = sqr(m2), m32 = sqr(m3);
   const bool m1_is_zero = is_zero(m1, EPSTOL);
   const bool m2_is_zero = is_zero(m2, EPSTOL);
   const bool m3_is_zero = is_zero(m3, EPSTOL);

   if ((m1_is_zero && m2_is_zero && m3_is_zero) ||
       (m2_is_zero && m3_is_zero) ||
       (m1_is_zero && m3_is_zero) ||
       (m1_is_zero && m2_is_zero)) {
      return 0.;
   }

   if (m1_is_zero) {
      if (is_close(m2,m3,EPSTOL))
         return -1./m22;
      return log(m32/m22)/(m22 - m32);
   }

   if (m2_is_zero) {
      if (is_close(m1,m3,EPSTOL))
         return -1./m12;
      return log(m32/m12)/(m12 - m32);
   }

   if (m3_is_zero) {
      if (is_close(m1,m2,EPSTOL))
         return -1./m12;
      return log(m22/m12)/(m12 - m22);
   }

   if (is_close(m2, m3, EPSTOL)) {
      if (is_close(m1, m2, EPSTOL))
         return - 0.5 / m22;
      return m12 / sqr(m12-m22) * log(m22/m12) + 1.0 / (m12 - m22);
   }

   if (is_close(m1, m2, EPSTOL)) {
      return - (1.0 + m32 / (m22-m32) * log(m32/m22)) / (m22-m32);
   }

   if (is_close(m1, m3, EPSTOL)) {
      return - (1.0 + m22 / (m32-m22) * log(m22/m32)) / (m32-m22);
   }

   return (m22 / (m12 - m22) * log(m22 / m12) -
           m32 / (m12 - m32) * log(m32 / m12)) / (m22 - m32);
}

/**
 * Derivative of B0(p^2, m1^2, m2^2, Q^2) w.r.t. p^2.
 *
 * @note Implemented only in the p^2 = 0 limit.
 *
 * @param p2 squared momentum
 * @param m2a squared mass
 * @param m2b squared mass
 *
 * @return derivative of B0 w.r.t. p^2
 */
double d1_b0(double /* p2 */, double m2a, double m2b) noexcept
{
   using std::abs;

   const double m4a = m2a * m2a;
   const double m4b = m2b * m2b;

   if ((std::abs(m2a) < 0.0001) != (std::abs(m2b) < 0.0001)) {
      return (m4a - m4b) / (2. * pow3(m2a - m2b));
   } else if (std::abs(m2a) < 0.0001 && std::abs(m2b) < 0.0001) {
      return 0.;
   } else if (std::abs(m2b - m2a) < 0.001) {
      return 1./(6. * m2a) + (m2a - m2b)/(12.* m4a);
   }

   return (m4a - m4b + 2. * m2a * m2b * std::log(m2b/m2a))
      /(2. * pow3(m2a - m2b));
}

} // namespace softsusy
