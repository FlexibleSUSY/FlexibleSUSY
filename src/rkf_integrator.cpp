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
 * @file rkf_integrator.cpp
 * @brief Implementation of the RKF_integrator class
 */

#include "rkf_integrator.hpp"
#include "logger.hpp"
#include "wrappers.hpp"

#include "config.h"

#ifdef ENABLE_ODEINT

#include <boost/numeric/odeint/algebra/algebra_dispatcher.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#include <boost/numeric/odeint/external/eigen/eigen_algebra_dispatcher.hpp>
#include <boost/numeric/odeint/external/eigen/eigen_resize.hpp>
#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>
#include <boost/numeric/odeint/stepper/generation.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_fehlberg78.hpp>

namespace boost {
namespace numeric {
namespace odeint {

template <typename Scalar, int R, int C, int O, int MR, int MC>
struct vector_space_norm_inf<Eigen::Array<Scalar,R,C,O,MR,MC> > {
   using result_type = typename Eigen::NumTraits<Scalar>::Real;
   result_type operator()(const Eigen::Array<Scalar,R,C,O,MR,MC>& x) const
      {
         return x.cwiseAbs().maxCoeff();
      }
};

} // namespace odeint
} // namespace numeric
} // namespace boost

namespace flexiblesusy {

namespace runge_kutta {

/**
 * The vector of the initial values of the parameters is
 * updated so that after calling this function, this vector contains
 * the updated values of the parameters at the end-point of the
 * integration.
 *
 * @param[in] start initial value of the independent variable
 * @param[in] end final value of the independent variable
 * @param[inout] pars initial values of the parameters
 * @param[in] derivs function calculating the derivatives
 * @param[in] tol desired accuracy to use in integration step
 */
void RKF_integrator::operator()(double start, double end,
                                Eigen::ArrayXd& pars, const Derivs& derivs,
                                double tol) const
{
   using state_type = Eigen::ArrayXd;
   using stepper_type = boost::numeric::odeint::runge_kutta_fehlberg78<
      state_type, double, state_type, double,
      boost::numeric::odeint::vector_space_algebra
      >;

   const double guess = (end - start) * 0.1; // first step size
   const auto derivatives = [derivs] (const state_type& y, state_type& dydt, double t) -> void {
      dydt = derivs(t, y);
   };

   const auto stepper = boost::numeric::odeint::make_controlled(tol, tol, stepper_type());
   boost::numeric::odeint::integrate_adaptive(
      stepper, derivatives, pars, start, end, guess, RKF_observer());
}

void RKF_integrator::RKF_observer::operator()(const Eigen::ArrayXd& state, double t) const
{
   if (!IsFinite(state)) {
      int max_step_dir = 0;
      for (int i = 0; i < state.size(); ++i) {
         if (!IsFinite(state(i))) {
            max_step_dir = i;
            break;
         }
      }
#ifdef ENABLE_VERBOSE
      ERROR("RKF_integrator: non-perturbative running at Q = "
            << std::exp(t) << " GeV of parameter y(" << max_step_dir
            << ") = " << state(max_step_dir));
#endif
      throw NonPerturbativeRunningError(std::exp(t), max_step_dir, state(max_step_dir));
   }
}

} // namespace runge_kutta

} // namespace flexiblesusy

#else

namespace flexiblesusy {

namespace runge_kutta {

/**
 * The vector of the initial values of the parameters is
 * updated so that after calling this function, this vector contains
 * the updated values of the parameters at the end-point of the
 * integration.
 *
 * @param[in] start initial value of the independent variable
 * @param[in] end final value of the independent variable
 * @param[inout] pars initial values of the parameters
 * @param[in] derivs function calculating the derivatives
 * @param[in] tol desired accuracy to use in integration step
 */
void RKF_integrator::operator()(double, double, Eigen::ArrayXd&, const Derivs&,
                                double) const
{
   throw DisabledOdeintError("Cannot call operator(), because odeint support is disabled.");
}

} // namespace runge_kutta

} // namespace flexiblesusy

#endif
