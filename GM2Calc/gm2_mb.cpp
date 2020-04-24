// ====================================================================
// This file is part of GM2Calc.
//
// GM2Calc is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// GM2Calc is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GM2Calc.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#include "gm2_mb.hpp"
#include "gm2_log.hpp"
#include "gm2_numerics.hpp"

#include <cmath>
#include <iostream>
#include <boost/math/tools/roots.hpp>

/**
 * \file gm2_mb.cpp
 *
 * Contains functions necessary to calculate the mb(Q) in the DR-bar
 * scheme.
 */

namespace gm2calc {

namespace mb {

namespace {
   const double Pi = 3.14159265358979323846;
} // anonymous namespace

/**
 * Calculates the strong coupling constant \f$\alpha_s(Q)\f$ in the
 * Standard Model with 5 active quark flavours in the MS-bar scheme at
 * the scale \f$Q\f$ using Eq (9) from arxiv:hep-ph/0207126 .
 *
 * @param scale renormalization scale \f$Q\f$
 * @param lambda_qcd \f$\Lambda_{\text{QCD}}\f$
 *
 * @return \f$\alpha_s(Q)\f$ MS-bar in the SM w/ 5 active flavours
 */
double calculate_alpha_s_SM5_at(double scale, double lambda_qcd) {
   const double t = std::log(sqr(scale/lambda_qcd));
   const double logt = std::log(t);

   return 12.*Pi/(23.*t) * (
      1 - 348./529. * logt/t
      + sqr(348./529.)/sqr(t) * (sqr(logt - 0.5) - 78073./242208.)
      );
}

/**
 * Calculates \f$\Lambda_{\text{QCD}}\f$ from a given
 * \f$\alpha_s(Q)\f$, where \f$\alpha_s(Q)\f$ is defined in the MS-bar
 * scheme in the SM with 5 active quark flavours.  A root finding
 * algorithm is used to solve Eq. (9) from arxiv:hep-ph/0207126
 * numerically for \f$\Lambda_{\text{QCD}}\f$.
 *
 * @param alpha \f$\alpha_s(Q)\f$ MS-bar, SM w/ 5 active quark flavours
 * @param scale renormalization scale Q
 * @param lambda_qcd_min minimum \f$\Lambda_{\text{QCD}}\f$
 * @param lambda_qcd_max maximum \f$\Lambda_{\text{QCD}}\f$
 * @param precision_goal accuracy goal
 * @param max_iterations maximum number of iterations
 *
 * @return \f$\Lambda_{\text{QCD}}\f$
 */
double calculate_lambda_qcd(double alpha, double scale,
                            double lambda_qcd_min = 0.001,
                            double lambda_qcd_max = 10.,
                            double precision_goal = 1e-10,
                            unsigned max_iterations = 1000)
{
   boost::uintmax_t it = max_iterations;

   // difference between goal alpha and alpha calculated using the
   // current lambda_qcd
   auto Difference_alpha = [alpha, scale](double lambda) -> double {
      const double alpha_current = calculate_alpha_s_SM5_at(scale, lambda);
      return alpha - alpha_current;
   };

   // stopping criterion, given two brackets a, b
   auto Stop_crit = [precision_goal](double a, double b) -> bool {
      return std::abs(a - b) < precision_goal;
   };

   double lambda_qcd = 0.217; // Nf = 5, PDG

   // find the root
   try {
      const std::pair<double,double> root =
         boost::math::tools::toms748_solve(Difference_alpha, lambda_qcd_min,
                                           lambda_qcd_max, Stop_crit, it);

      lambda_qcd = 0.5 * (root.first + root.second);
   } catch (const std::exception& e) {
      WARNING("Could not determine lambda_QCD: " << e.what()
              << ".  Using lambda_QCD = " << lambda_qcd);
   }

   if (it >= max_iterations) {
      const double precision = std::abs(Difference_alpha(lambda_qcd));
      WARNING("Calculation of Lambda_QCD did not converge"
              " (reached accuracy: " << precision <<
              ", accuracy goal: " << precision_goal <<
              ", max. iterations: " << max_iterations << ")");
   }

   return lambda_qcd;
}

/**
 * Calculates \f$F_b(\mu)\f$, Eq (5) of arxiv:hep-ph/0207126 .
 *
 * @param alpha strong coupling constant at the scale \f$\mu\f$
 *
 * @return \f$F_b(\mu)\f$, Eq (5) of arxiv:hep-ph/0207126
 */
double Fb(double alpha) {
   const double as = alpha / Pi;

   return std::pow(23.*as/6., 12./23.) * (
      1 + 3731./3174. * as + 1.500706 * as * as
      );
}

/**
 * MS-bar to DR-bar conversion for mb, Eq (11) of arxiv:hep-ph/0207126 .
 *
 * @param alpha strong coupling constant
 *
 * @return MS-bar to DR-bar conversion factor Eq (11)
 */
double conversion_mb_MSbar_to_DRbar(double alpha) {
   const double as = alpha / Pi;

   return 1. - as/3. - 29./72. * as * as;
}

} // namespace mb

/**
 * Calculates mb(Q) in the DR-bar scheme in the SM w/ 5 active quark
 * flavours using the approach described in arxiv:hep-ph/0207126 .
 *
 * @param mb_mb mb(mb) MS-bar in SM w/ 5 active quark flavours
 * @param alpha_s alpha_s MS-bar in SM w/ 5 quark flavours at scale Q
 * @param scale renormalization scale Q
 *
 * @return mb(Q) DR-bar in the SM w/ 5 quarks
 */
double calculate_mb_SM5_DRbar(
   double mb_mb, double alpha_s, double scale)
{
   // determine Lambda_QCD
   const double lambda_qcd = mb::calculate_lambda_qcd(alpha_s, scale);

   // calculate alpha_s(mb)
   const double alpha_s_mb = mb::calculate_alpha_s_SM5_at(mb_mb, lambda_qcd);

   // run mb to destination scale
   // Here alpha_s must be given at the destination scale `scale'.
   const double mb = mb_mb * mb::Fb(alpha_s) / mb::Fb(alpha_s_mb);

   // DR-bar conversion
   const double mb_DRbar = mb * mb::conversion_mb_MSbar_to_DRbar(alpha_s);

   return mb_DRbar;
}

} // namespace gm2calc
