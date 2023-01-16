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

#include "amm_loop_functions.hpp"
#include "Li2.hpp"
#include <cmath>
#include <complex>

namespace flexiblesusy {

/**
 * Barr-Zee 2-loop function with fermion loop (arXiv:1502.04199 Eq 26).
 *
 * @param m squared mass ratio
*/
double BarrZeeLoopFPS(double m)
{
   const double pi = 3.1415926535897932;

   // @todo(alex): check stability for small values
   if (m == 0) {
      return 0;
   }

   if (m == 0.25) {
      return 0.69314718055994531; // Log[2]
   }

   double r1, theta1, r2, theta2;

   if (m <= 0.25) {
      const double y = std::sqrt(1 - 4*m);
      r1 = 1 - (1 - y)/(2*m);
      r2 = 1 - (1 + y)/(2*m);
      if (r1 > 0) {
         theta1 = 0;
      } else {
         r1 *= -1;
         theta1 = pi;
      }
      if (r2 > 0) {
         theta2 = 0;
      } else {
         r2 *= -1;
         theta2 = pi;
      }
      return m / y * (Li2(std::polar(r1, theta1)).real() - Li2(std::polar(r2, theta2)).real());
   } else {
      const double y = std::sqrt(-1 + 4*m);
      const double real = 1 - 1/(2*m);
      r1 = r2 = 1;

      theta1 = std::atan2(y/(2*m), real);
      theta2 = std::atan2(-y/(2*m), real);

      return m / y * (Li2(std::polar(r1, theta1)).imag() - Li2(std::polar(r2, theta2)).imag());
   }
}

/**
 * Barr-Zee 2-loop function with fermion loop (arXiv:1502.04199 Eq 25).
 *
 * @param m squared mass ratio.
 */
double BarrZeeLoopFS(double m)
{
    if (m == 0) {
        return 0;
    }

    return (2*m - 1) * BarrZeeLoopFPS(m) - m * (2 + std::log(m));
}

/**
 * Barr-Zee 2-loop function with scalar loop (arXiv:1502.04199 Eq 27).
 *
 * @param m squared mass ratio
 */
double BarrZeeLoopS(double m)
{
   return 1 + std::log(m)/2 - BarrZeeLoopFPS(m);
}

/**
 * Barr-Zee 2-loop function with vector loop (arXiv:1502.04199 Eq 28).
 *
 * @param m squared mass ratio
 */
double BarrZeeLoopV(double m)
{
   if (m == 0.25) {
      return 4.75;
   }

   double j;

   double r1, theta1, r2, theta2;

   if (m < 0.25) {
      const double y = std::sqrt(1.0-4.0*m);
      r1 = 2 / (1-y);
      theta1 = std::atan2(0, r1);
      r1 = std::abs(r1);
      r2 = 2 / (1+y);
      theta2 = std::atan2(0, r2);
      r2 = std::abs(r2);

      j = 1/y * (std::log(std::abs(y-1) / (y+1)) * std::log(m) + Li2(std::polar(r1, theta1)).real() - Li2(std::polar(r2, theta2)).real());
   }
   else {
      const double y = std::sqrt(-1.0+4.0*m);
      r1 = r2 = std::sqrt(1/m);
      const double real = 1/(2*m);
      const double imag = y/(2*m);
      theta1 = std::atan2(imag, real);
      theta2 = std::atan2(-imag, real);
      j = 1/y * ( (std::atan2(y, -1) - std::atan2(y, 1)) * std::log(m) + Li2(std::polar(r1, theta1)).imag() - Li2(std::polar(r2, theta2)).imag());
   }

   return BarrZeeLoopS(m) + 15.0/2.0 * m * (2.0 + std::log(m)) + m/2 * (19 - 12*m) * j - 9*m * BarrZeeLoopFPS(m);
}

} // namespace flexiblesusy
