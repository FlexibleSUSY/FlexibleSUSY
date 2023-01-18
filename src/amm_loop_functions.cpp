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
#include "error.hpp"
#include "Cl2.hpp"
#include "Li2.hpp"
#include <cmath>
#include <complex>

namespace flexiblesusy {

/**
 * Barr-Zee 2-loop function with fermion loop (arXiv:1502.04199 Eq 26).
 *
 * @param z squared mass ratio
*/
double BarrZeeLoopFPS(double z)
{
   if (z < 0) {
      throw OutOfBoundsError("BarrZeeLoopFPS: argument must not be negative.");
   } else if (z == 0) {
      return 0;
   } else if (z < 0.25) {
      const double y = std::sqrt(1 - 4*z); // 0 < y < 1
      const double c = -4.9348022005446793; // -Pi^2/2
      const double q = (1 + y)/(1 - y);
      const double lq = std::log(q);
      return z/y*(2*Li2(1 + q) - lq*(std::log(z) - 0.5*lq) + c);
   } else if (z == 0.25) {
      return 0.69314718055994531; // Log[4]
   }

   // z > 0.25
   const double y = std::sqrt(-1 + 4*z);
   const double theta = std::atan2(y, 2*z - 1);
   return 2*z/y*Cl2(theta);
}

/**
 * Barr-Zee 2-loop function with fermion loop (arXiv:1502.04199 Eq 25).
 *
 * @param z squared mass ratio.
 */
double BarrZeeLoopFS(double z)
{
   if (z < 0) {
      throw OutOfBoundsError("BarrZeeLoopFS: argument must not be negative.");
   } else if (z == 0) {
      return 0;
   }

   return (2*z - 1)*BarrZeeLoopFPS(z) - z*(2 + std::log(z));
}

/**
 * Barr-Zee 2-loop function with scalar loop (arXiv:1502.04199 Eq 27).
 *
 * @param z squared mass ratio
 */
double BarrZeeLoopS(double z)
{
   if (z < 0) {
      throw OutOfBoundsError("BarrZeeLoopS: argument must not be negative.");
   } else if (z == 0) {
      return 0;
   }

   return 1 + 0.5*std::log(z) - BarrZeeLoopFPS(z);
}

/**
 * Barr-Zee 2-loop function with vector loop (arXiv:1502.04199 Eq 28).
 *
 * @param z squared mass ratio
 */
double BarrZeeLoopV(double z)
{
   if (z < 0) {
      throw OutOfBoundsError("BarrZeeLoopV: argument must not be negative.");
   } else if (z == 0) {
      return 0;
   } else if (z == 0.25) {
      return 4.75;
   }

   double j, r1, theta1, r2, theta2;

   if (z < 0.25) {
      const double y = std::sqrt(1.0-4.0*z);
      r1 = 2 / (1-y);
      theta1 = std::atan2(0, r1);
      r1 = std::abs(r1);
      r2 = 2 / (1+y);
      theta2 = std::atan2(0, r2);
      r2 = std::abs(r2);

      j = 1/y * (std::log(std::abs(y-1) / (y+1)) * std::log(z) + Li2(std::polar(r1, theta1)).real() - Li2(std::polar(r2, theta2)).real());
   }
   else {
      const double y = std::sqrt(-1.0+4.0*z);
      r1 = r2 = std::sqrt(1/z);
      const double real = 1/(2*z);
      const double imag = y/(2*z);
      theta1 = std::atan2(imag, real);
      theta2 = std::atan2(-imag, real);
      j = 1/y * ( (std::atan2(y, -1) - std::atan2(y, 1)) * std::log(z) + Li2(std::polar(r1, theta1)).imag() - Li2(std::polar(r2, theta2)).imag());
   }

   return BarrZeeLoopS(z) + 15.0/2.0 * z * (2.0 + std::log(z)) + z/2 * (19 - 12*z) * j - 9*z * BarrZeeLoopFPS(z);
}

} // namespace flexiblesusy
