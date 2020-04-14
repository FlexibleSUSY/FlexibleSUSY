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
 * @file betafunction.cpp
 * @brief contains implemenation of class Beta_function
 */

#include "betafunction.hpp"
#include "basic_rk_integrator.hpp"
#include "error.hpp"

#include <cmath>

namespace flexiblesusy {

Beta_function::Beta_function()
   : integrator(runge_kutta::Basic_rk_integrator<Eigen::ArrayXd>())
{
}

void Beta_function::reset()
{
   num_pars = 0;
   loops = 0;
   thresholds = 0;
   scale = 0.0;
   tolerance = 1.e-4;
   min_tolerance = 1.0e-11;
   zero_threshold = 1.e-11;
   integrator = runge_kutta::Basic_rk_integrator<Eigen::ArrayXd>();
}

/**
 * Runs parameter objects of the generated models from current scale
 * to the scale x2 passed as in an argument.
 *
 * @param x2 renormalization scale to run parameters to
 * @param eps RG running precision
 */
void Beta_function::run_to(double x2, double eps)
{
   const double tol = get_tolerance(eps);
   run(scale, x2, tol);
}

/**
 * Runs parameter objects of the generated models from scale x1,
 * passed as first argument, to scale x2 passed as second argument.
 *
 * @param x1 renormalization scale to start RG running from
 * @param x2 renormalization scale to run parameters to
 * @param eps RG running precision
 */
void Beta_function::run(double x1, double x2, double eps)
{
   if (get_loops() > 0) {
      const double tol = get_tolerance(eps);

      if (std::fabs(x1) < tol)
         throw NonPerturbativeRunningError(x1);
      if (std::fabs(x2) < tol)
         throw NonPerturbativeRunningError(x2);

      Eigen::ArrayXd y(get());

      Derivs derivs = [this] (double x, const Eigen::ArrayXd& y) {
         return derivatives(x, y);
      };

      call_rk(x1, x2, y, derivs, tol);
      set(y);
   }

   set_scale(x2);
}

/**
 * Takes logarithm of renormalisation scale as first argument and
 * parameters of RGE passed in as an Eigen::ArrayXd object of dynamic
 * size in the second argument.  Returns the beta functions as an
 * Eigen::ArrayXd object.
 *
 * @param x logarithm of renormalization scale to calculate beta functions at
 * @param y array of model parameters
 *
 * @return array of beta functions
 */
Eigen::ArrayXd Beta_function::derivatives(double x, const Eigen::ArrayXd& y)
{
   set_scale(exp(x));
   set(y);
   return beta();
}

/**
 * Calls Runge Kutta routine, passing the start and end scales
 * as first and second argument respectively.  Parameters are passed
 * as third argument and the derivatives method as the fourth
 * argument.  The precison is passed as the last argument or if not
 * passed it is set to default value tolerance which is a private data
 * member of the class set to 1e-4 in the constructor.
 *
 * @param x1 renormalization scale to start RG running from
 * @param x2 renormalization scale to run parameters to
 * @param v array of model parameters
 * @param derivs function which calculates the derivatives (beta functions)
 * @param eps RG running precision
 */
void Beta_function::call_rk(double x1, double x2, Eigen::ArrayXd & v,
                            Derivs derivs, double eps)
{
   if (fabs(x1 - x2) < min_tolerance)
      return;

   const double start = std::log(fabs(x1));
   const double end = std::log(fabs(x2));
   const double tol = get_tolerance(eps);

   integrator(start, end, v, derivs, tol);

   set_scale(x2);
}

/**
 * Helper function, which returns the RG running precision, determined
 * from the argument eps.  If eps is less than zero, the default
 * running precision is returned.  Otherwise, if eps is positive and
 * less than the minimum allowed running precision, the minimum
 * allowed running precision is returned.  Otherwise, the value of eps
 * is returned.
 *
 * @param eps value to determine the RG running precision from
 *
 * @return RG running precision
 */
double Beta_function::get_tolerance(double eps)
{
   double tol;
   if (eps < 0.0)
      tol = tolerance;
   else if (eps < min_tolerance)
      tol = min_tolerance;
   else
      tol = eps;

   return tol;
}

} // namespace flexiblesusy
