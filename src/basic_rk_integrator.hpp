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
 * @file basic_rk_integrator.hpp
 * @brief Integration of ODEs by Runge-Kutta
 */

#ifndef BASIC_RK_INTEGRATOR_H
#define BASIC_RK_INTEGRATOR_H

#include <functional>

#include "rk.hpp"

namespace flexiblesusy {

namespace runge_kutta {

/**
 * @class Basic_rk_stepper
 * @brief Class to carry out a 5th order Runge-Kutta step
 *
 * @tparam StateType type of parameters vector
 * @tparam Derivs type of object returning the values of the derivatives
 */
template <typename StateType, typename Derivs>
class Basic_rk_stepper {
public:
   /// @brief Carries out a variable step-size Runge-Kutta step
   double step(StateType&, const StateType&, double&, double,
               double, const StateType&, Derivs, int&) const;
private:
   /// @brief Carries out a single 5th order Runge-Kutta step
   void runge_kutta_step(const StateType&, const StateType&, double,
                         double, StateType&, StateType&, Derivs) const;
};

/**
 * The step is calculated using the given fixed step-size.  In addition
 * to returning the estimate for the parameters at the next step, an
 * estimate for the error is also returned.
 *
 * @param[in] y current values of the parameters
 * @param[in] dydx current values of the parameter derivatives
 * @param[in] x current value of the independent variable
 * @param[in] h step-size to use
 * @param[out] yout updated values of the parameters
 * @param[out] yerr estimated truncation error
 * @param[in] derivs function calculating the derivatives
 */
template <typename StateType, typename Derivs>
void Basic_rk_stepper<StateType, Derivs>::runge_kutta_step(
   const StateType& y, const StateType& dydx, double x,
   double h, StateType& yout, StateType& yerr, Derivs derivs) const
{
   rungeKuttaStep(y, dydx, x, h, yout, yerr, derivs);
}

/**
 * The initial values of the independent and dependent variables
 * are updated to their new values after calling this function, i.e.,
 * the vector \c y contains the approximate values of the dependent
 * variables after carrying out the step, and \c x contains the
 * new value of the independent variable.
 *
 * @param[inout] y current values of the parameters
 * @param[in] dydx current values of the parameter derivatives
 * @param[inout] x current value of the independent variable
 * @param[in] htry initial step-size to try
 * @param[in] eps desired error tolerance
 * @param[in] yscal vector of scale values for fraction errors
 * @param[in] derivs function calculating the derivatives
 * @param[out] max_step_dir parameter with largest estimated error
 * @return estimated next step-size to use
 */
template <typename StateType, typename Derivs>
double Basic_rk_stepper<StateType,Derivs>::step(
   StateType& y, const StateType& dydx, double& x, double htry,
   double eps, const StateType& yscal, Derivs derivs,
   int& max_step_dir) const
{
   return odeStepper(y, dydx, x, htry, eps, yscal, derivs, max_step_dir);
}

/**
 * @class Basic_rk_integrator
 * @brief Class for integrating a system of first order ODEs
 *
 * @tparam StateType type of parameters vector
 * @tparam Derivs type of object returning the values of the derivatives
 * @tparam Stepper type of object implementing Runge-Kutta step
 */
template <typename StateType,
          typename Derivs
          = std::function<StateType(double, const StateType&)>,
          typename Stepper = Basic_rk_stepper<StateType,Derivs> >
class Basic_rk_integrator {
public:
   /// @brief Integrates the system over an interval
   void operator()(double start, double end, StateType& ystart,
                   Derivs derivs, double tolerance) const;

   /// @brief Sets the maximum number of allowed steps in the integration
   /// @param s maximum number of steps to allow
   void set_max_steps(int s) { max_steps = s; }

   /// @brief Returns the maximum number of allowed steps in the integration
   /// @return maximum number of steps to allow
   int get_max_steps() const { return max_steps; }

private:
   int max_steps{400}; ///< Maximum number of steps in integration
   Stepper stepper{};  ///< Stepper to provide a Runge-Kutta step
};

/**
 * The vector of the initial values of the parameters is
 * updated so that after calling this function, this vector contains
 * the updated values of the parameters at the end-point of the
 * integration.
 *
 * @param[in] start initial value of the independent variable
 * @param[in] end final value of the independent variable
 * @param[inout] ystart initial values of the parameters
 * @param[in] derivs function calculating the derivatives
 * @param[in] tolerance desired accuracy to use in integration step
 */
template <typename StateType, typename Derivs, typename Stepper>
void Basic_rk_integrator<StateType, Derivs, Stepper>::operator()(
   double start, double end, StateType& ystart, Derivs derivs,
   double tolerance) const
{
   const double guess = (start - end) * 0.1; // first step size
   const double hmin = (start - end) * tolerance * 1.0e-5;
   const auto rkqs = [this] (
      StateType& y, const StateType& dydx, double& x, double htry,
      double eps, const StateType& yscal, Derivs derivs,
      int& max_step_dir) -> double {
      return this->stepper.step(y, dydx, x, htry, eps,
                                yscal, derivs, max_step_dir);
   };

   integrateOdes(ystart, start, end, tolerance, guess, hmin,
                 derivs, rkqs, max_steps);
}

} // namespace runge_kutta

} // namespace flexiblesusy

#endif
