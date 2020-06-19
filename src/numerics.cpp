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

double sign(double x) noexcept
{
   return x >= 0.0 ? 1.0 : -1.0;
}

// can be made constexpr in C++20
double fB(const std::complex<double>& a) noexcept
{
   using flexiblesusy::fast_log;

   const double re = std::real(a);
   const double im = std::imag(a);

   if ((std::abs(re) == 0.0 || std::abs(re) == 1.0) && im == 0.0) {
      return -1.0;
   }

   return std::real(-1.0 + fast_log(1.0 - a) - a*fast_log(1.0 - 1.0/a));
}

} // anonymous namespace

double a0(double m, double q) noexcept
{
   constexpr double TOL = 1e-4;

   if (std::abs(m) < TOL) {
      return 0.0;
   }

   return sqr(m) * (1.0 - 2. * std::log(std::abs(m / q)));
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
#ifdef USE_LOOPTOOLS
   setmudim(q*q);
   double b0l = B0(p*p, m1*m1, m2*m2).real();
   //  return B0(p*p, m1*m1, m2*m2).real();
#endif

   m1 = std::abs(m1);
   m2 = std::abs(m2);
   p = std::abs(p);
   q = std::abs(q);

   if (m1 > m2) {
      std::swap(m1, m2);
   }

   // protect against infrared divergence
   if (is_zero(p, EPSTOL) && is_zero(m1, EPSTOL) && is_zero(m2, EPSTOL)) {
      return 0.0;
   }

   // p is not 0
   if (p > 1.0e-5 * m2) {
      const double m12 = sqr(m1), m22 = sqr(m2), p2 = sqr(p);
      const double s = p2 - m22 + m12;
      const std::complex<double> imin(m12, -EPSTOL);
      const std::complex<double> x = std::sqrt(sqr(s) - 4.0 * p2 * imin);
      const std::complex<double> xp  = (s + sign(s)*x) / (2*p2);
      const std::complex<double> xm = imin / (xp*p2);

      return -2.0*std::log(p/q) - fB(xp) - fB(xm);
   }

   if (is_close(m1, m2, EPSTOL)) {
      return -2.0*std::log(m1/q);
   }

   if (m1 < 1.0e-15) {
      return 1.0 - 2.0*std::log(m2/q);
   }

   const double m12 = sqr(m1), m22 = sqr(m2);

   return 1.0 - 2.0 * std::log(m2/q)
        + 2.0 * m12 * std::log(m2/m1) / (m12 - m22);
}

/// Note that b1 is NOT symmetric in m1 <-> m2!!!
double b1(double p, double m1, double m2, double q) noexcept
{
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

   if (std::abs(m1) > 1.0e-15 && std::abs(m2) > 1.0e-15) {
      const double m14 = sqr(m12), m24 = sqr(m22);
      const double m16 = m12*m14 , m26 = m22*m24;
      const double m18 = sqr(m14), m28 = sqr(m24);
      const double p4 = sqr(p2);
      const double diff = m12 - m22;

      if (std::abs(diff) < pTolerance * std::max(m12, m22)) {
         return
            - 0.5*std::log(m22/q2)
            + 1.0/12.0*p2/m22 + 1.0/120.0*p4/m24
            + diff*(-1.0/6.0/m22 - 1.0/30.0*p2/m24 - 1.0/140.0*p4/m26)
            + sqr(diff)*(1.0/24.0/m24 + 1.0/60.0*p2/m26 + 3.0/560.0*p4/m28);
      }

      const double l12 = std::log(m12/m22);

      return (3*m14 - 4*m12*m22 + m24 - 2*m14*l12)/(4.*sqr(diff))
         + (p2*(4*pow3(diff)*
                (2*m14 + 5*m12*m22 - m24) +
                (3*m18 + 44*m16*m22 - 36*m14*m24 - 12*m12*m26 + m28)*p2
                - 12*m14*m22*(2*sqr(diff) + (2*m12 + 3*m22)*p2)*l12))/
         (24.*pow6(diff)) - 0.5*std::log(m22/q2);
   }

   return (m12 > m22)
      ? -0.5*std::log(m12/q2) + 0.75
      : -0.5*std::log(m22/q2) + 0.25;
}

double b22(double p,  double m1, double m2, double q) noexcept
{
#ifdef USE_LOOPTOOLS
   setmudim(q*q);
   double b22l = B00(p*p, m1*m1, m2*m2).real();
#endif

   p = std::abs(p);
   m1 = std::abs(m1);
   m2 = std::abs(m2);
   q = std::abs(q);

   // protect against infrared divergence
   if (is_zero(p, EPSTOL) && is_zero(m1, EPSTOL) && is_zero(m2, EPSTOL))
      return 0.0;

   /// Decides level at which one switches to p=0 limit of calculations
   const double p2 = sqr(p), m12 = sqr(m1), m22 = sqr(m2);
   const double pTolerance = 1.0e-10;

   if (p2 < pTolerance * std::max(m12, m22)) {
      // m1 == m2 with good accuracy
      if (is_close(m1, m2, EPSTOL)) {
         return -m12 * std::log(m1/q) + m12 * 0.5;
      }
      // p == 0 limit
      if (m1 > EPSTOL && m2 > EPSTOL) {
         return 0.375 * (m12 + m22) - 0.5 *
            (sqr(m22) * std::log(m2/q) -
             sqr(m12) * std::log(m1/q)) / (m22 - m12);
      }
      return (m1 < EPSTOL)
         ? 0.375 * m22 - 0.5 * m22 * std::log(m2/q)
         : 0.375 * m12 - 0.5 * m12 * std::log(m1/q);
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
   const double m1sq = sqr(m1);
   const double m2sq = sqr(m2);

   if (is_close(m1, m2, EPSTOL)) {
      const double m3sq = sqr(m3);
      const double m4sq = sqr(m4);

      if (is_zero(m2, EPSTOL)) {
         // d0 is undefined for m1 == m2 == 0
         return 0.;
      } else if (is_zero(m3, EPSTOL)) {
         return (-m2sq + m4sq - m2sq * std::log(m4sq/m2sq))/
            sqr(m2 * m2sq - m2 * m4sq);
      } else if (is_zero(m4, EPSTOL)) {
         return (-m2sq + m3sq - m2sq * std::log(m3sq/m2sq))/
            sqr(m2 * m2sq - m2 * m3sq);
      } else if (is_close(m2, m3, EPSTOL) && is_close(m2, m4, EPSTOL)) {
         return 1.0 / (6.0 * sqr(m2sq));
      } else if (is_close(m2, m3, EPSTOL)) {
         return (sqr(m2sq) - sqr(m4sq) + 2.0 * m4sq * m2sq * std::log(m4sq / m2sq)) /
            (2.0 * m2sq * sqr(m2sq - m4sq) * (m2sq - m4sq));
      } else if (is_close(m2, m4, EPSTOL)) {
         return (sqr(m2sq) - sqr(m3sq) + 2.0 * m3sq * m2sq * std::log(m3sq / m2sq)) /
            (2.0 * m2sq * sqr(m2sq - m3sq) * (m2sq - m3sq));
      } else if (is_close(m3, m4, EPSTOL)) {
         return -1.0 / sqr(m2sq - m3sq) *
            ((m2sq + m3sq) / (m2sq - m3sq) * std::log(m3sq / m2sq) + 2.0);
      }

      return
         (m4sq / sqr(m2sq - m4sq) * std::log(m4sq / m2sq) +
          m4sq / (m2sq * (m2sq - m4sq)) -
          m3sq / sqr(m2sq - m3sq) * std::log(m3sq / m2sq) -
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
      return std::log(m32/m22)/(m22 - m32);
   }

   if (m2_is_zero) {
      if (is_close(m1,m3,EPSTOL))
         return -1./m12;
      return std::log(m32/m12)/(m12 - m32);
   }

   if (m3_is_zero) {
      if (is_close(m1,m2,EPSTOL))
         return -1./m12;
      return std::log(m22/m12)/(m12 - m22);
   }

   if (is_close(m2, m3, EPSTOL)) {
      if (is_close(m1, m2, EPSTOL))
         return - 0.5 / m22;
      return m12 / sqr(m12-m22) * std::log(m22/m12) + 1.0 / (m12 - m22);
   }

   if (is_close(m1, m2, EPSTOL)) {
      return - (1.0 + m32 / (m22-m32) * std::log(m32/m22)) / (m22-m32);
   }

   if (is_close(m1, m3, EPSTOL)) {
      return - (1.0 + m22 / (m32-m22) * std::log(m22/m32)) / (m32-m22);
   }

   return (m22 / (m12 - m22) * std::log(m22 / m12) -
           m32 / (m12 - m32) * std::log(m32 / m12)) / (m22 - m32);
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

double c00(double m1, double m2, double m3, double q) noexcept
{
  // taken from Package X
  const double m12 = sqr(m1), m22 = sqr(m2), m32 = sqr(m3), q2 = sqr(q);

  double ans = 0.;

  if (is_close(m1, 0., EPSTOL) && is_close(m2, 0., EPSTOL)
      && is_close(m3, 0., EPSTOL)) {
      // IR singularity
     ans = 0.;
  } else if (is_close(m2, 0., EPSTOL) && is_close(m3, 0., EPSTOL)) {
     ans = 3./8. + 1./4*std::log(q2/m12);
  } else if (is_close(m1, 0., EPSTOL) && is_close(m3, 0., EPSTOL)) {
     ans = 3./8. + 1./4*std::log(q2/m22);
  } else if (is_close(m1, 0., EPSTOL) && is_close(m2, 0., EPSTOL)) {
     ans = 3./8. + 1./4*std::log(q2/m32);
  } else if (is_close(m1, 0., EPSTOL)) {
     if (is_close(m2, m3, EPSTOL)) {
        ans = 1./8 + 1./4*std::log(q2/m22);
     } else {
        ans = 3./8 - m22*std::log(m22/m32)/(4*(m22-m32)) + 1./4*std::log(q2/m32);
     }
  } else if (is_close(m2, 0., EPSTOL)) {
     if (is_close(m1, m3, EPSTOL)) {
        ans = 1./8 + 1./4*std::log(q2/m12);
     } else {
        ans = 3./8 - m12*std::log(m12/m32)/(4*(m12-m32)) + 1./4*std::log(q2/m32);
     }
  } else if (is_close(m3, 0., EPSTOL)) {
     if (is_close(m1, m2, EPSTOL)) {
        ans = 1./8 + 1./4*std::log(q2/m12);
     } else {
        ans = 3./8 - m22*std::log(m12/m22)/(4*(m12-m22)) + 1./4*std::log(q2/m12);
     }
  } else if (is_close(m2, m3, EPSTOL)) {
    if (is_close(m1, m2, EPSTOL)) {
        ans = 1./4*std::log(q2/m12);
    } else {
        ans = (3*m12-m22)/(8*(m12-m22))
            - 1./4*(sqr(m12)*std::log(m12/m22))/sqr(m12-m22) + 1./4*std::log(q2/m22);
    }
  } else if (is_close(m1, m2, EPSTOL)) {
     ans = (3*m32-m12)/(8*(m32-m12)) - 1./4*(sqr(m32)*std::log(m32/m12))/sqr(m32-m12)
         + 1./4*std::log(q2/m12);
  } else if (is_close(m1, m3, EPSTOL)) {
     ans = (3*m22-m12)/(8*(m22-m12)) - 1./4*(sqr(m22)*std::log(m22/m12))/sqr(m22-m12)
         + 1./4*std::log(q2/m12);
  } else {
     ans = 3./8 - 1./4*sqr(m12)*std::log(m12/m32)/((m12-m22)*(m12-m32))
         - 1./4*sqr(m22)*std::log(m22/m32)/((m22-m12)*(m22-m32)) + 1./4*std::log(q2/m32);
  }

  return ans;
}

} // namespace softsusy
