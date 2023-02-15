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
namespace amm_loop_functions {
namespace two_loop {

/**
 * Barr-Zee 2-loop function with fermion loop and pseudoscalar and
 * photon mediators (arXiv:1502.04199 Eq 26).
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
 * Barr-Zee 2-loop function with fermion loop and scalar and photon
 * mediators (arXiv:1502.04199 Eq 25).
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
 * Barr-Zee 2-loop function with fermion loop and pseudoscalar and Z
 * boson mediators.
 *
 * @param x squared mass ratio (mf/ms)^2.
 * @param y squared mass ratio (mf/mz)^2.
 */
double BarrZeeLoopFPSZ(double x, double y)
{
   if (x < 0 || y < 0) {
      throw OutOfBoundsError("BarrZeeLoopFPSZ: arguments must not be negative.");
   } else if (y == 0) {
      return 0;
   } else if (x == y) {
      if (x == 0.25) {
         return -0.29543145370663021; // (1 - Log[16])/6
      }
      return (BarrZeeLoopFS(x) + 2*x)/(1 - 4*x);
   }

   return y*(BarrZeeLoopFPS(x) - BarrZeeLoopFPS(y))/(y - x);
}

/**
 * Barr-Zee 2-loop function with fermion loop and scalar and Z boson
 * mediators.
 *
 * @param x squared mass ratio (mf/ms)^2.
 * @param y squared mass ratio (mf/mz)^2.
 */
double BarrZeeLoopFSZ(double x, double y)
{
   if (x < 0 || y < 0) {
      throw OutOfBoundsError("BarrZeeLoopFSZ: arguments must not be negative.");
   } else if (x == y) {
      return 0; // @todo(alex) improve limit
   }

   return y*(BarrZeeLoopFS(x) - BarrZeeLoopFS(y))/(y - x);
}

/**
 * Barr-Zee 2-loop function with scalar loop and scalar and photon
 * mediators (arXiv:1502.04199 Eq 27).
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
 * Barr-Zee 2-loop function with vector loop and scalar and photon
 * mediators (arXiv:1502.04199 Eq 28).
 *
 * @param z squared mass ratio
 */
double BarrZeeLoopV(double z)
{
   if (z < 0) {
      throw OutOfBoundsError("BarrZeeLoopV: argument must not be negative.");
   } else if (z == 0) {
      return 0; // actually -inf; return 0 to avoid propagation of inf
   }

   return (17./2 - 15*z)*BarrZeeLoopFPS(z) + (0.5 + 7.5*z)*(2 + std::log(z));
}

} // namespace two_loop
} // namespace amm_loop_functions
} // namespace flexiblesusy
