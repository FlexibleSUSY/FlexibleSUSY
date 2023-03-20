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
#include <limits>

namespace flexiblesusy {
namespace amm_loop_functions {
namespace two_loop {

namespace {

void sort(double& x, double& y)
{
   if (x > y) {
      std::swap(x, y);
   }
}

} // anonymous namespace

/**
 * Barr-Zee 2-loop function with fermion loop and pseudoscalar and
 * photon mediators (arXiv:1502.04199 Eq 26).
 *
 * @param z squared mass ratio
*/
double BarrZeeLoopFP(double z)
{
   if (z < 0) {
      throw OutOfBoundsError("BarrZeeLoopFP: argument must not be negative.");
   } else if (z == 0) {
      return 0;
   } else if (z < std::numeric_limits<double>::epsilon()) {
      constexpr double pi26 = 1.6449340668482264; // Pi^2/6
      const double lz = std::log(z);
      return z*(pi26 + 0.5*lz*lz);
   } else if (z < 0.25) {
      const double y = std::sqrt(1 - 4*z); // 0 < y < 1
      constexpr double c = -4.9348022005446793; // -Pi^2/2
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
   } if (z > 1e2) {
      const double lz = std::log(z);
      const double iz = 1/z;
      return (-13./18 - 1./3*lz) + iz*(-26./300 - 15./300*lz
         + iz*(-673./44100 - 420./44100*lz + iz*(-971./317520 - 630./317520*lz)));
   }

   return (2*z - 1)*BarrZeeLoopFP(z) - z*(2 + std::log(z));
}

/**
 * Barr-Zee 2-loop function with fermion loop and pseudoscalar and Z
 * boson mediators.
 *
 * @param x squared mass ratio (mf/ms)^2.
 * @param y squared mass ratio (mf/mz)^2.
 */
double BarrZeeLoopFPZ(double x, double y)
{
   if (x < 0 || y < 0) {
      throw OutOfBoundsError("BarrZeeLoopFPZ: arguments must not be negative.");
   }

   sort(x, y);

   constexpr double eps = 1e-8;

   if (x == 0 || y == 0) {
      return 0;
   } else if (std::abs(1 - x/y) < eps) {
      if (std::abs(x - 0.25) < eps) {
         // -(1 + 2*Log[2])/6 + O(x - 1/4)
         return -0.39771572685331510 - 0.87269032593060833*(x - 0.25);
      }
      return x*(2*BarrZeeLoopFP(x) + std::log(x))/(1 - 4*x);
   }

   return (y*BarrZeeLoopFP(x) - x*BarrZeeLoopFP(y))/(x - y);
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
   } else if (y == 0) {
      return 0;
   } else if (std::abs(1 - x/y) < 1e-7) {
      if (std::abs(x - 0.25) < 1e-7) {
         const double d = x - 0.25;
         // 2/3*(Log[2] - 1) + O(x - 1/4)
         return -0.20456854629336979 + d*(0.30903548889591250 - 0.67527189418658333*d);
      } else if (x >= 1e3) {
         const double ix = 1/x;
         const double lx = std::log(x);
         return -1./3 + ix*(11./300 + 1./20*lx + ix*(463./22050 + 2./105*lx
            + ix*(761./105840 + 1./168*lx + ix*(8707./4.002075e6 + 2./1155*lx))));
      }
      return (x*(3 - 12*x - 2*(-1 + 3*x)*std::log(x))
              + (1 + 6*x*(-1 + 2*x))*BarrZeeLoopFP(x)
         )/(4*x - 1);
   }

   return y*(BarrZeeLoopFS(x) - BarrZeeLoopFS(y))/(x - y);
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

   return 1 + 0.5*std::log(z) - BarrZeeLoopFP(z);
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
   } if (z >= 1e2) {
      const double lz = std::log(z);
      const double iz = 1/z;
      return 89./12 + 42./12*lz + iz*(284./360 + 165./360*lz
         + iz*(6199./44100 + 3885./44100*lz
         + iz*(30017./1.0584e6 + 19530./1.0584e6*lz
         + iz*(83351./1.37214e7 + 55440./1.37214e7*lz
         + iz*(34978051./2.597186592e10 + 23603580./2.597186592e10*lz)))));
   }

   return (17./2 - 15*z)*BarrZeeLoopFP(z) + (0.5 + 7.5*z)*(2 + std::log(z));
}

} // namespace two_loop
} // namespace amm_loop_functions
} // namespace flexiblesusy
