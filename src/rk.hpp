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

/**
 * @file rk.hpp
 * @brief Integration of ODEs by Runge-Kutta
 * @author Ben Allanach, Alexander Voigt
 *
 * The implementation of the Runge-Kutta routines have been derived
 * from SOFTSUSY [hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305].
 */

#ifndef RK_H
#define RK_H

#include <algorithm>
#include <cmath>

#include "logger.hpp"
#include "error.hpp"

namespace flexiblesusy {

namespace runge_kutta {

namespace {
/// Returns |a| with sign of b in front
inline double sign(double a, double b) noexcept {
   return b >= 0 ? std::fabs(a) : -std::fabs(a);
}
} // anonymous namespace

/// A single step of Runge Kutta (5th order), input:
/// y and dydx (derivative of y), x is independent variable. yout is value
/// after step. derivs is a user-supplied function
template <typename ArrayType, typename Derivs>
void rungeKuttaStep(const ArrayType& y, const ArrayType& dydx, double x,
		    double h, ArrayType& yout, ArrayType& yerr, Derivs derivs)
{
   const double a2 = 0.2;
   const double a3 = 0.3;
   const double a4 = 0.6;
   const double a5 = 1.0;
   const double a6 = 0.875;
   const double b21 = 0.2;
   const double b31 = 3.0 / 40.0;
   const double b32 = 9.0 / 40.0;
   const double b41 = 0.3;
   const double b42 = -0.9;
   const double b43 = 1.2;
   const double b51 = -11.0 / 54.0;
   const double b52 = 2.5;
   const double b53 = -70.0 / 27.0;
   const double b54 = 35.0 / 27.0;
   const double b61 = 1631.0 / 55296.0;
   const double b62 = 175.0 / 512.0;
   const double b63 = 575.0 / 13824.0;
   const double b64 = 44275.0 / 110592.0;
   const double b65 = 253.0 / 4096.0;
   const double c1 = 37.0 / 378.0;
   const double c3 = 250.0 / 621.0;
   const double c4 = 125.0 / 594.0;
   const double c6 = 512.0 / 1771.0;
   const double dc5 = -277.00 / 14336.0;
   const double dc1 = c1 - 2825.0 / 27648.0;
   const double dc3 = c3 - 18575.0 / 48384.0;
   const double dc4 = c4 - 13525.0 / 55296.0;
   const double dc6 = c6 - 0.25;

   ArrayType ytemp = b21 * h * dydx + y;
   const ArrayType ak2 = derivs(x + a2 * h, ytemp);

   // Allowing piece-wise calculating of ytemp for speed reasons
   ytemp = y + h * (b31 * dydx + b32 * ak2);
   const ArrayType ak3 = derivs(x + a3 * h, ytemp);

   ytemp = y + h * (b41 * dydx + b42 * ak2 + b43 * ak3);
   const ArrayType ak4 = derivs(x+a4*h,ytemp);

   ytemp = y + h * (b51 * dydx + b52 * ak2 + b53 * ak3 + b54 * ak4);
   const ArrayType ak5 = derivs(x + a5 * h, ytemp);

   ytemp = y + h * (b61 * dydx + b62 * ak2 + b63 * ak3 + b64 * ak4 + b65 * ak5);
   const ArrayType ak6 = derivs(x + a6 * h, ytemp);

   yout = y + h * (c1 * dydx + c3 * ak3 + c4 * ak4 + c6 * ak6);
   yerr = h * (dc1 * dydx + dc3 * ak3 + dc4 * ak4 + dc5 * ak5 + dc6 * ak6);
}

/// organises the variable step-size for Runge-Kutta evolution
template <typename ArrayType, typename Derivs>
double odeStepper(ArrayType& y, const ArrayType& dydx, double& x, double htry,
                  double eps, const ArrayType& yscal, Derivs derivs,
                  int& max_step_dir)
{
   const double SAFETY = 0.9;
   const double PGROW = -0.2;
   const double PSHRNK = -0.25;
   const double ERRCON = 1.89e-4;
   const int n = y.size();
   double errmax;
   double h = htry;
   ArrayType yerr(n);
   ArrayType ytemp(n);

   for (;;) {
      rungeKuttaStep(y, dydx, x, h, ytemp, yerr, derivs);
      errmax = (yerr / yscal).abs().maxCoeff(&max_step_dir);
      errmax  /= eps;
      if (!std::isfinite(errmax)) {
#ifdef ENABLE_VERBOSE
         ERROR("odeStepper: non-perturbative running at Q = "
               << std::exp(x) << " GeV of parameter y(" << max_step_dir
               << ") = " << y(max_step_dir) << ", dy(" << max_step_dir
               << ")/dx = " << dydx(max_step_dir));
#endif
         throw NonPerturbativeRunningError(std::exp(x), max_step_dir, y(max_step_dir));
      }
      if (errmax <= 1.0) {
         break;
      }
      const double htemp = SAFETY * h * std::pow(errmax, PSHRNK);
      h = (h >= 0.0 ? std::max(htemp, 0.1 * h) : std::min(htemp, 0.1 * h));
      if (x + h == x) {
#ifdef ENABLE_VERBOSE
         ERROR("At Q = " << std::exp(x) << " GeV "
               "stepsize underflow in odeStepper in parameter y("
               << max_step_dir << ") = " << y(max_step_dir) << ", dy("
               << max_step_dir << ")/dx = " << dydx(max_step_dir));
#endif
         throw NonPerturbativeRunningError(std::exp(x), max_step_dir, y(max_step_dir));
      }
   }
   x += h;
   y = ytemp;

   return errmax > ERRCON ? SAFETY * h * std::pow(errmax,PGROW) : 5.0 * h;
}

/// Organises integration of 1st order system of ODEs
template <typename ArrayType, typename Derivs,
          typename Stepper = decltype(runge_kutta::odeStepper<ArrayType,Derivs>)>
void integrateOdes(ArrayType& ystart, double from, double to, double eps,
                   double h1, double hmin, Derivs derivs,
                   Stepper rkqs = runge_kutta::odeStepper<ArrayType,Derivs>, int max_steps = 400)
{
   const int nvar = ystart.size();
   const double TINY = 1.0e-16;
   double x = from;
   double h = sign(h1, to - from);
   ArrayType yscal(nvar);
   ArrayType y(ystart);
   ArrayType dydx;
   int max_step_dir = 0;

   for (int nstp = 0; nstp < max_steps; ++nstp) {
      dydx = derivs(x, y);
      yscal = y.abs() + (dydx * h).abs() + TINY;
      if ((x + h - to) * (x + h - from) > 0.0) {
         h = to - x;
      }

      const double hnext = rkqs(y, dydx, x, h, eps, yscal, derivs, max_step_dir);

      if ((x - to) * (to - from) >= 0.0) {
         ystart = y;
         return;
      }

      h = hnext;

      if (std::fabs(hnext) <= hmin) {
         break;
      }
   }

#ifdef ENABLE_VERBOSE
   ERROR("Bailed out of rk.cpp:too many steps in integrateOdes\n"
         "********** Q = " << std::exp(x) << " *********");
   ERROR("max step in direction of " << max_step_dir);
   for (int i = 0; i < nvar; i++)
      ERROR("y(" << i << ") = " << y(i) << " dydx(" << i <<
            ") = " << dydx(i));
#endif

   throw NonPerturbativeRunningError(std::exp(x), max_step_dir, y(max_step_dir));
}

} // namespace runge_kutta

} // namespace flexiblesusy

#endif // RK_H
