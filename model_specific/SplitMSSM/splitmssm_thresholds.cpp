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

#include "splitmssm_thresholds.hpp"
#include "logger.hpp"
#include "threshold_loop_functions.hpp"
#include "numerics2.hpp"
#include <cmath>

#define ASSERT_NON_ZERO(p,fun)                  \
   do {                                         \
      if (is_zero(p))                           \
         FATAL(#fun ": " #p " is zero!");       \
   } while (false)

namespace flexiblesusy {
namespace splitmssm_thresholds {

namespace {

constexpr double Pi      = 3.141592653589793;
constexpr double oneLoop = 6.332573977646110963e-03; // 1/(4 Pi)^2

double sqr(double x) noexcept { return x*x; }
double pow3(double x) noexcept { return x*x*x; }
double pow4(double x) noexcept { return sqr(sqr(x)); }
double pow5(double x) noexcept { return pow4(x)*x; }
double abs_sqrt(double x) noexcept { return std::sqrt(std::abs(x)); }

} // anonymous namespace

std::ostream& operator<<(std::ostream& ostr, const Parameters& parameters)
{
   ostr
      << "scale = " << parameters.scale << ", "
      << "g1 = " << parameters.g1 << ", "
      << "g2 = " << parameters.g2 << ", "
      << "g3 = " << parameters.g3 << ", "
      << "yt = " << parameters.gt << ", "
      << "At = " << parameters.At << ", "
      << "mu = " << parameters.mu << ", "
      << "mA = " << parameters.mA << ", "
      << "m1 = " << parameters.m1 << ", "
      << "m2 = " << parameters.m2 << ", "
      << "tan_beta = " << parameters.tan_beta << ", "
      << "mq2(0,0) = " << parameters.mq2(0,0) << ", "
      << "mq2(1,1) = " << parameters.mq2(1,1) << ", "
      << "mq2(2,2) = " << parameters.mq2(2,2) << ", "
      << "mu2(0,0) = " << parameters.mu2(0,0) << ", "
      << "mu2(1,1) = " << parameters.mu2(1,1) << ", "
      << "mu2(2,2) = " << parameters.mu2(2,2) << ", "
      << "md2(0,0) = " << parameters.md2(0,0) << ", "
      << "md2(1,1) = " << parameters.md2(1,1) << ", "
      << "md2(2,2) = " << parameters.md2(2,2) << ", "
      << "ml2(0,0) = " << parameters.ml2(0,0) << ", "
      << "ml2(1,1) = " << parameters.ml2(1,1) << ", "
      << "ml2(2,2) = " << parameters.ml2(2,2) << ", "
      << "me2(0,0) = " << parameters.me2(0,0) << ", "
      << "me2(1,1) = " << parameters.me2(1,1) << ", "
      << "me2(2,2) = " << parameters.me2(2,2) << '\n';

   return ostr;
}

/**
 * Returns lambda at the tree-level, first term on the rhs of Eq (8)
 * of arXiv:1407.4081.
 *
 * @param parameters running parameters
 *
 * @return first term on the rhs of Eq. (8) of arXiv:1407.4081
 */
double lambda_tree_level(const Parameters& parameters)
{
   const double g1 = parameters.g1;
   const double g2 = parameters.g2;
   const double tan_beta = parameters.tan_beta;
   const double beta = std::atan(tan_beta);
   const double cos_2beta = std::cos(2*beta);
   const double lambda = 0.25 * (sqr(g2) + 0.6*sqr(g1)) * sqr(cos_2beta);

   return lambda;
}

/**
 * Returns \f$\tilde{g}_{1u}\f$ at the tree-level, first term on the
 * rhs of Eq (17) of arXiv:1407.4081.
 *
 * @param parameters running parameters
 *
 * @return first term on the rhs of Eq. (17) of arXiv:1407.4081
 */
double gYu_tree_level(const Parameters& parameters)
{
   const double g1 = parameters.g1;
   const double tan_beta = parameters.tan_beta;
   const double beta = std::atan(tan_beta);
   const double sin_beta = std::sin(beta);
   const double gYu = sqrt(0.6) * g1 * sin_beta;

   return gYu;
}

/**
 * Returns \f$\tilde{g}_{1d}\f$ at the tree-level, first term on the
 * rhs of Eq (18) of arXiv:1407.4081.
 *
 * @param parameters running parameters
 *
 * @return first term on the rhs of Eq. (18) of arXiv:1407.4081
 */
double gYd_tree_level(const Parameters& parameters)
{
   const double g1 = parameters.g1;
   const double tan_beta = parameters.tan_beta;
   const double beta = std::atan(tan_beta);
   const double cos_beta = std::cos(beta);
   const double gYd = sqrt(0.6) * g1 * cos_beta;

   return gYd;
}

/**
 * Returns \f$\tilde{g}_{2u}\f$ at the tree-level, first term on the
 * rhs of Eq (15) of arXiv:1407.4081.
 *
 * @param parameters running parameters
 *
 * @return first term on the rhs of Eq. (15) of arXiv:1407.4081
 */
double g2u_tree_level(const Parameters& parameters)
{
   const double g2 = parameters.g2;
   const double tan_beta = parameters.tan_beta;
   const double beta = std::atan(tan_beta);
   const double sin_beta = std::sin(beta);
   const double g2u = g2 * sin_beta;

   return g2u;
}

/**
 * Returns \f$\tilde{g}_{2d}\f$ at the tree-level, first term on the
 * rhs of Eq (16) of arXiv:1407.4081.
 *
 * @param parameters running parameters
 *
 * @return first term on the rhs of Eq. (16) of arXiv:1407.4081
 */
double g2d_tree_level(const Parameters& parameters)
{
   const double g2 = parameters.g2;
   const double tan_beta = parameters.tan_beta;
   const double beta = std::atan(tan_beta);
   const double cos_beta = std::cos(beta);
   const double g2d = g2 * cos_beta;

   return g2d;
}

/**
 * Finite 1-loop MS-bar--DR-bar transition counterterm for quartic
 * Higgs coupling (arXiv:1407.4081, Eq. (9)).
 *
 * @param parameters running parameters
 *
 * @return 1-loop MS-bar--DR-bar transition CT for quartic Higgs
 * coupling
 */
double delta_lambda_1loop_reg(const Parameters& parameters)
{
   const double g1 = parameters.g1;
   const double g2 = parameters.g2;
   const double tan_beta = parameters.tan_beta;
   const double beta = std::atan(tan_beta);
   const double cos_2beta = std::cos(2*beta);

   const double delta_lambda_1loop_reg
      = - 0.09*pow4(g1) - 0.3*sqr(g1)*sqr(g2)
      - (0.75 - sqr(cos_2beta)/6.) * pow4(g2);

   return delta_lambda_1loop_reg * oneLoop;
}

/**
 * Finite 1-loop threshold correction for quartic Higgs coupling from
 * sfermions (arXiv:1407.4081, Eq. (10)).
 *
 * @param parameters running parameters
 *
 * @return 1-loop threshold coupling for quartic Higgs coupling from
 * sfermions
 */
double delta_lambda_1loop_phi(const Parameters& parameters)
{
   using namespace threshold_loop_functions;

   const double g1 = parameters.g1;
   const double g2 = parameters.g2;
   const double g12 = sqr(g1);
   const double g22 = sqr(g2);
   const double g14 = sqr(g12);
   const double g24 = sqr(g22);
   const double tan_beta = parameters.tan_beta;
   const double beta = std::atan(tan_beta);
   const double cos_2beta = std::cos(2*beta);
   const double cos2_2beta = sqr(cos_2beta);
   const double cos_4beta = std::cos(4*beta);
   const double cos_8beta = std::cos(8*beta);
   const double sin_4beta = std::sin(4*beta);
   const double gt = parameters.gt;
   const double gt2 = sqr(gt);
   const double gt4 = sqr(gt2);
   const double scale2 = sqr(parameters.scale);
   const Eigen::Matrix<double,3,3> mq2(parameters.mq2);
   const Eigen::Matrix<double,3,3> mu2(parameters.mu2);
   const Eigen::Matrix<double,3,3> md2(parameters.md2);
   const Eigen::Matrix<double,3,3> ml2(parameters.ml2);
   const Eigen::Matrix<double,3,3> me2(parameters.me2);
   const double xQU = abs_sqrt(mq2(2,2)/mu2(2,2));
   const double At = parameters.At;
   const double mu = parameters.mu;
   const double xt = At - mu/tan_beta;
   const double xtt = sqr(xt)/abs_sqrt(mq2(2,2)*mu2(2,2));
   const double mA2 = sqr(parameters.mA);

   ASSERT_NON_ZERO(mq2(0,0), delta_lambda_1loop_phi);
   ASSERT_NON_ZERO(mq2(1,1), delta_lambda_1loop_phi);
   ASSERT_NON_ZERO(mq2(2,2), delta_lambda_1loop_phi);
   ASSERT_NON_ZERO(mu2(0,0), delta_lambda_1loop_phi);
   ASSERT_NON_ZERO(mu2(1,1), delta_lambda_1loop_phi);
   ASSERT_NON_ZERO(mu2(2,2), delta_lambda_1loop_phi);
   ASSERT_NON_ZERO(md2(0,0), delta_lambda_1loop_phi);
   ASSERT_NON_ZERO(md2(1,1), delta_lambda_1loop_phi);
   ASSERT_NON_ZERO(md2(2,2), delta_lambda_1loop_phi);
   ASSERT_NON_ZERO(ml2(0,0), delta_lambda_1loop_phi);
   ASSERT_NON_ZERO(ml2(1,1), delta_lambda_1loop_phi);
   ASSERT_NON_ZERO(ml2(2,2), delta_lambda_1loop_phi);
   ASSERT_NON_ZERO(me2(0,0), delta_lambda_1loop_phi);
   ASSERT_NON_ZERO(me2(1,1), delta_lambda_1loop_phi);
   ASSERT_NON_ZERO(me2(2,2), delta_lambda_1loop_phi);
   ASSERT_NON_ZERO(mA2     , delta_lambda_1loop_phi);
   ASSERT_NON_ZERO(tan_beta, delta_lambda_1loop_phi);
   ASSERT_NON_ZERO(scale2  , delta_lambda_1loop_phi);

   const double delta_lambda_1loop_phi =
      3*gt2*(gt2+0.5*(g22-g12/5.)*cos_2beta)*std::log(mq2(2,2)/scale2)
      +3*gt2*(gt2+0.4*g12*cos_2beta)*std::log(mu2(2,2)/scale2)
      +cos2_2beta/300*(
         3*(g14+25*g24)*(
            +std::log(mq2(0,0)/scale2)
            +std::log(mq2(1,1)/scale2)
            +std::log(mq2(2,2)/scale2)
         )
         +24*g14*(
            +std::log(mu2(0,0)/scale2)
            +std::log(mu2(1,1)/scale2)
            +std::log(mu2(2,2)/scale2)
         )
         +6*g14*(
            +std::log(md2(0,0)/scale2)
            +std::log(md2(1,1)/scale2)
            +std::log(md2(2,2)/scale2)
         )
         +(9*g14+25*g24)*(
            +std::log(ml2(0,0)/scale2)
            +std::log(ml2(1,1)/scale2)
            +std::log(ml2(2,2)/scale2)
         )
         +18*g14*(
            +std::log(me2(0,0)/scale2)
            +std::log(me2(1,1)/scale2)
            +std::log(me2(2,2)/scale2)
         )
      )
      +1./4800.*(261*g14+630*g12*g22+1325*g24
                 -4*cos_4beta*(9*g14+90*g12*g22+175*g24)
                 -9*cos_8beta*sqr(3*g12+5*g22))*std::log(mA2/scale2)
      -3./16.*sqr(3./5.*g12+g22)*sqr(sin_4beta)
      +6*gt4*xtt*(F1(xQU)-xtt/12*F2(xQU))
      +3./4.*gt2*xtt*cos_2beta*(3./5.*g12*F3(xQU)+g22*F4(xQU))
      -0.25*gt2*xtt*cos2_2beta*(3./5.*g12+g22)*F5(xQU);

   return delta_lambda_1loop_phi * oneLoop;
}

namespace {
   /**
    * Split-MSSM contribution to the beta-function of lambda
    * (arXiv:1407.4081, Eq. (12)).
    */
   double beta_lambda(double lambda, double gYu, double gYd, double g2u, double g2d)
   {
      const double gYu2 = sqr(gYu);
      const double gYd2 = sqr(gYd);
      const double g2u2 = sqr(g2u);
      const double g2d2 = sqr(g2d);
      const double gYu4 = sqr(gYu2);
      const double gYd4 = sqr(gYd2);
      const double g2u4 = sqr(g2u2);
      const double g2d4 = sqr(g2d2);

      const double beta =
         2*lambda*(gYd2 + gYu2 + 3*g2d2 + 3*g2u2)
         -gYd4 - gYu4 - 5*g2d4 - 5*g2u4
         -4*gYd*gYu*g2d*g2u
         -2*(gYd2+g2u2)*(gYu2+g2d2);

      return beta;
   }
} // anonymous namespace

/**
 * Finite 1-loop threshold correction for quartic Higgs coupling from
 * Higgsinos and gauginos (arXiv:1407.4081, Eq. (11)).
 *
 * The couplings \f$\tilde{g}_{1u}\f$, \f$\tilde{g}_{1d}\f$,
 * \f$\tilde{g}_{2u}\f$, \f$\tilde{g}_{2d}\f$, \f$\lambda\f$ are set
 * to their tree-level values.
 *
 * @param parameters running parameters
 *
 * @return 1-loop threshold correction for quartic Higgs coupling from
 * Higgsinos and gauginos
 */
double delta_lambda_1loop_chi_1(const Parameters& parameters)
{
   const double scale = parameters.scale;
   const double mu = parameters.mu;
   const double gYu = gYu_tree_level(parameters);
   const double gYd = gYd_tree_level(parameters);
   const double g2u = g2u_tree_level(parameters);
   const double g2d = g2d_tree_level(parameters);
   const double m1 = parameters.m1;
   const double m2 = parameters.m2;
   const double lambda = lambda_tree_level(parameters);

   ASSERT_NON_ZERO(scale, delta_lambda_1loop_chi_1);
   ASSERT_NON_ZERO(mu   , delta_lambda_1loop_chi_1);
   ASSERT_NON_ZERO(m1   , delta_lambda_1loop_chi_1);
   ASSERT_NON_ZERO(m2   , delta_lambda_1loop_chi_1);

   const double delta_lambda =
      delta_lambda_1loop_chi_1(scale,mu,lambda,gYu,gYd,g2u,g2d,m1,m2);

   return delta_lambda;
}

/**
 * Finite 1-loop threshold correction for quartic Higgs coupling from
 * Higgsinos and gauginos (arXiv:1407.4081, Eq. (11)).
 *
 * @param scale renormalization scale
 * @param mu bilinear Higgsino coupling
 * @param lambda quartic Higgs coupling
 * @param gYu Higgs-Up-type-Higgsino-bino coupling
 * @param gYd Higgs-Down-type-Higgsino-bino coupling
 * @param g2u Higgs-Up-type-Higgsino-wino coupling
 * @param g2d Higgs-Down-type-Higgsino-wino coupling
 * @param m1 bino mass parameter
 * @param m2 wino mass parameter
 *
 * @return r.h.s. of Eq. (11), including the loop factor
 */
double delta_lambda_1loop_chi_1(
   double scale,
   double mu,
   double lambda,
   double gYu,
   double gYd,
   double g2u,
   double g2d,
   double m1,
   double m2
)
{
   using namespace threshold_loop_functions;

   const double gYu2 = sqr(gYu);
   const double gYd2 = sqr(gYd);
   const double g2u2 = sqr(g2u);
   const double g2d2 = sqr(g2d);
   const double gYu4 = sqr(gYu2);
   const double gYd4 = sqr(gYd2);
   const double g2u4 = sqr(g2u2);
   const double g2d4 = sqr(g2d2);
   const double r1 = m1 / mu;
   const double r2 = m2 / mu;
   const double mu2 = sqr(mu);
   const double scale2 = sqr(scale);

   ASSERT_NON_ZERO(scale, delta_lambda_1loop_chi_1);
   ASSERT_NON_ZERO(mu   , delta_lambda_1loop_chi_1);
   ASSERT_NON_ZERO(m1   , delta_lambda_1loop_chi_1);
   ASSERT_NON_ZERO(m2   , delta_lambda_1loop_chi_1);

   const double delta_lambda =
      0.5*beta_lambda(lambda,gYu,gYd,g2u,g2d)*std::log(mu2/scale2)
      -7./12.*f1(r1)*(gYd4+gYu4)
      -9./4.*f2(r2)*(g2d4+g2u4)
      -3./2.*f3(r1)*gYd2*gYu2
      -7./2.*f4(r2)*g2d2*g2u2
      -8./3.*f5(r1,r2)*gYd*gYu*g2d*g2u
      -7./6.*f6(r1,r2)*(gYd2*g2d2+gYu2*g2u2)
      -1./6.*f7(r1,r2)*(gYd2*g2u2+gYu2*g2d2)
      -4./3.*f8(r1,r2)*(gYd*g2u+gYu*g2d)*(gYd*g2d+gYu*g2u)
      +2./3.*f(r1)*gYd*gYu*(lambda-2*(gYd2+gYu2))
      +2*f(r2)*g2d*g2u*(lambda-2*(g2d2+g2u2))
      +1./3.*g(r1)*lambda*(gYd2+gYu2)
      +g(r2)*lambda*(g2d2+g2u2);

   return oneLoop * delta_lambda;
}

/**
 * Finite 1-loop threshold correction for quartic Higgs coupling from
 * Higgsinos and gauginos (arXiv:1407.4081, Eq. (13)).
 *
 * @param parameters input parameters
 *
 * @return 1-loop threshold correction for quartic Higgs coupling from
 * Higgsinos and gauginos
 */
double delta_lambda_1loop_chi_2(const Parameters& parameters)
{
   const double scale = parameters.scale;
   const double mu = parameters.mu;
   const double m2 = parameters.m2;
   const double g1 = parameters.g1;
   const double g2 = parameters.g2;
   const double tan_beta = parameters.tan_beta;

   const double delta_lambda =
      delta_lambda_1loop_chi_2(scale,mu,m2,g1,g2,tan_beta);

   return delta_lambda;
}

/**
 * Finite 1-loop threshold correction for quartic Higgs coupling from
 * Higgsinos and gauginos (arXiv:1407.4081, Eq. (13)).
 *
 * @param scale renormalization scale
 * @param mu bilinear Higgsino coupling
 * @param m2 wino mass parameter
 * @param g1 running EFT gauge coupling \f$g_1\f$ (GUT normalized)
 * @param g2 running EFT gauge coupling \f$g_1\f$
 * @param tan_beta \f$\tan\beta\f$
 *
 * @return r.h.s. of Eq. (13), including the loop factor
 */
double delta_lambda_1loop_chi_2(
   double scale,
   double mu,
   double m2,
   double g1,
   double g2,
   double tan_beta)
{
   const double scale2 = sqr(scale);
   const double mu2 = sqr(mu);
   const double m22 = sqr(m2);
   const double g14 = pow4(g1);
   const double g24 = pow4(g2);
   const double beta = std::atan(tan_beta);
   const double cos2_2beta = sqr(std::cos(2*beta));

   ASSERT_NON_ZERO(scale, delta_lambda_1loop_chi_2);
   ASSERT_NON_ZERO(mu   , delta_lambda_1loop_chi_2);
   ASSERT_NON_ZERO(m2   , delta_lambda_1loop_chi_2);

   const double delta_lambda =
      -1./6. * cos2_2beta * (
         +2.*g24*std::log(m22/scale2)
         +(9./25.*g14 + g24)*std::log(mu2/scale2)
         );

   return oneLoop * delta_lambda;
}

/**
 * Finite 2-loop threshold correction for quartic Higgs coupling from
 * sfermions (arXiv:1407.4081, Eq. (34)).
 *
 * @param parameters input parameters
 *
 * @return 2-loop threshold coupling for quartic Higgs coupling from
 * sfermions (arXiv:1407.4081, Eq. (34)).
 */
double delta_lambda_2loop_phi(const Parameters& parameters)
{
   const double Pi4 = pow4(Pi);
   const double g3 = parameters.g3;
   const double g32 = sqr(g3);
   const double gt = parameters.gt;
   const double gt4 = pow4(gt);
   const double scale2 = sqr(parameters.scale);
   const Eigen::Matrix<double,3,3> mq2(parameters.mq2);
   const Eigen::Matrix<double,3,3> mu2(parameters.mu2);
   const double xQU = abs_sqrt(mq2(2,2)/mu2(2,2));
   const double xQU2 = sqr(xQU);
   const double At = parameters.At;
   const double mu = parameters.mu;
   const double tan_beta = parameters.tan_beta;
   const double xt = At - mu/tan_beta;
   const double xtt = sqr(xt)/abs_sqrt(mq2(2,2)*mu2(2,2));

   ASSERT_NON_ZERO(mq2(2,2), delta_lambda_2loop_phi);
   ASSERT_NON_ZERO(mu2(2,2), delta_lambda_2loop_phi);
   ASSERT_NON_ZERO(scale2  , delta_lambda_2loop_phi);
   ASSERT_NON_ZERO(tan_beta, delta_lambda_2loop_phi);

   double delta_lambda_2loop_phi;

   if (is_equal(xQU, 1.)) {
      delta_lambda_2loop_phi =
         3. - 2*xtt + sqr(xtt)/6.;
   } else {
      delta_lambda_2loop_phi =
         3.
         +4*std::log(xQU)
         +8*sqr(std::log(xQU))
         +6*sqr(std::log(mq2(2,2)/scale2))
         -4*(1+3*std::log(xQU))*std::log(mq2(2,2)/scale2)
         +xtt*(
            +(12*xQU*std::log(xQU))/(xQU2-1)*(2*std::log(mq2(2,2)/scale2)-1)
            -(16*xQU*(xQU2-2)*sqr(std::log(xQU)))/sqr(xQU2-1)
         )
         +sqr(xtt)*(
            +(6*xQU2*(5+xQU2)*std::log(xQU))/pow3(xQU2-1)
            +(4*xQU2*(pow4(xQU)-4*xQU2-5)*sqr(std::log(xQU)))/pow4(xQU2-1)
            -(10*xQU2)/sqr(xQU2-1)
            +(12*xQU2)/sqr(xQU2-1)*(1-(xQU2+1)/(xQU2-1)*std::log(xQU))*std::log(mq2(2,2)/scale2)
         );
   }

   return -(g32*gt4)/(32*Pi4) * delta_lambda_2loop_phi;
}

/**
 * Simplified, finite 2-loop threshold correction for quartic Higgs
 * coupling from sfermions (arXiv:1407.4081, Eq. (36)) in the
 * High-Scale SUSY scenario.
 *
 * @param parameters input parameters
 *
 * @return 2-loop threshold coupling for quartic Higgs coupling from
 * sfermions (arXiv:1407.4081, Eq. (36))
 */
double delta_lambda_2loop_phi_HSS(const Parameters& parameters)
{
   const double Pi4 = pow4(Pi);
   const double g32 = sqr(parameters.g3);
   const double gt4 = pow4(parameters.gt);
   const double scale = parameters.scale;
   const double At = parameters.At;
   const double mu = parameters.mu;
   const double tan_beta = parameters.tan_beta;
   const double xt = At - mu/tan_beta;
   const double r = xt / scale;

   ASSERT_NON_ZERO(scale   , delta_lambda_2loop_phi_HSS);
   ASSERT_NON_ZERO(tan_beta, delta_lambda_2loop_phi_HSS);

   const double delta_lambda_2loop_phi_HSS =
      -12.*r - 6.*sqr(r) + 14.*pow3(r) + 0.5*pow4(r) - pow5(r);

   return (g32*gt4)/(96.*Pi4) * delta_lambda_2loop_phi_HSS;
}

/**
 * Returns \f$\tilde{g}_{1u}\f$ at the one-loop level, all except the
 * first term on the rhs of Eq (17) of arXiv:1407.4081.
 *
 * @param parameters running parameters
 *
 * @return one-loop term on the rhs of Eq. (17) of arXiv:1407.4081
 */
double delta_gYu_1loop(const Parameters& parameters)
{
   const double g1 = parameters.g1;
   const double g2 = parameters.g2;
   const double g12 = sqr(g1);
   const double g22 = sqr(g2);
   const double tan_beta = parameters.tan_beta;
   const double beta = std::atan(tan_beta);
   const double cos2_beta = sqr(std::cos(beta));
   const double sin_beta = std::sin(beta);
   const double sin2_beta = sqr(sin_beta);
   const double gt = parameters.gt;
   const double gt2 = sqr(gt);
   const double scale2 = sqr(parameters.scale);
   const Eigen::Matrix<double,3,3> mq2(parameters.mq2);
   const Eigen::Matrix<double,3,3> mu2(parameters.mu2);
   const Eigen::Matrix<double,3,3> md2(parameters.md2);
   const Eigen::Matrix<double,3,3> ml2(parameters.ml2);
   const Eigen::Matrix<double,3,3> me2(parameters.me2);
   const double mA2 = sqr(parameters.mA);

   ASSERT_NON_ZERO(mq2(0,0), delta_gYu_1loop);
   ASSERT_NON_ZERO(mq2(1,1), delta_gYu_1loop);
   ASSERT_NON_ZERO(mq2(2,2), delta_gYu_1loop);
   ASSERT_NON_ZERO(mu2(0,0), delta_gYu_1loop);
   ASSERT_NON_ZERO(mu2(1,1), delta_gYu_1loop);
   ASSERT_NON_ZERO(mu2(2,2), delta_gYu_1loop);
   ASSERT_NON_ZERO(md2(0,0), delta_gYu_1loop);
   ASSERT_NON_ZERO(md2(1,1), delta_gYu_1loop);
   ASSERT_NON_ZERO(md2(2,2), delta_gYu_1loop);
   ASSERT_NON_ZERO(ml2(0,0), delta_gYu_1loop);
   ASSERT_NON_ZERO(ml2(1,1), delta_gYu_1loop);
   ASSERT_NON_ZERO(ml2(2,2), delta_gYu_1loop);
   ASSERT_NON_ZERO(me2(0,0), delta_gYu_1loop);
   ASSERT_NON_ZERO(me2(1,1), delta_gYu_1loop);
   ASSERT_NON_ZERO(me2(2,2), delta_gYu_1loop);
   ASSERT_NON_ZERO(mA2     , delta_gYu_1loop);
   ASSERT_NON_ZERO(tan_beta, delta_gYu_1loop);
   ASSERT_NON_ZERO(scale2  , delta_gYu_1loop);

   const double delta_gYu_1loop =
      (3*g22/16)*(-2+7*cos2_beta)
      +(3*g12/80)*(-44+7*cos2_beta)
      +(9*gt2)/(4*sin2_beta)
      +(4*g12-9*(g12+5*g22)*cos2_beta)/40*std::log(mA2/scale2)
      +(g12/10)*(
         +std::log(ml2(0,0)/scale2)+std::log(ml2(1,1)/scale2)+std::log(ml2(2,2)/scale2)
         +2*std::log(me2(0,0)/scale2)+2*std::log(me2(1,1)/scale2)+2*std::log(me2(2,2)/scale2)
      )
      +(g12/30)*(
         +std::log(mq2(0,0)/scale2)+std::log(mq2(1,1)/scale2)+std::log(mq2(2,2)/scale2)
         +8*std::log(mu2(0,0)/scale2)+8*std::log(mu2(1,1)/scale2)+8*std::log(mu2(2,2)/scale2)
         +2*std::log(md2(0,0)/scale2)+2*std::log(md2(1,1)/scale2)+2*std::log(md2(2,2)/scale2)
      )
      +gt2/(4*sin2_beta)*(7*std::log(mq2(2,2)/scale2)-13*std::log(mu2(2,2)/scale2));

   return delta_gYu_1loop * oneLoop * g1 * std::sqrt(0.6) * sin_beta;
}

/**
 * Returns \f$\tilde{g}_{1d}\f$ at the one-loop level, all except the
 * first term on the rhs of Eq (18) of arXiv:1407.4081.
 *
 * @param parameters running parameters
 *
 * @return one-loop term on the rhs of Eq. (18) of arXiv:1407.4081
 */
double delta_gYd_1loop(const Parameters& parameters)
{
   const double g1 = parameters.g1;
   const double g2 = parameters.g2;
   const double g12 = sqr(g1);
   const double g22 = sqr(g2);
   const double tan_beta = parameters.tan_beta;
   const double beta = std::atan(tan_beta);
   const double cos_beta = std::cos(beta);
   const double sin_beta = std::sin(beta);
   const double sin2_beta = sqr(sin_beta);
   const double scale2 = sqr(parameters.scale);
   const Eigen::Matrix<double,3,3> mq2(parameters.mq2);
   const Eigen::Matrix<double,3,3> mu2(parameters.mu2);
   const Eigen::Matrix<double,3,3> md2(parameters.md2);
   const Eigen::Matrix<double,3,3> ml2(parameters.ml2);
   const Eigen::Matrix<double,3,3> me2(parameters.me2);
   const double mA2 = sqr(parameters.mA);

   ASSERT_NON_ZERO(mq2(0,0), delta_gYd_1loop);
   ASSERT_NON_ZERO(mq2(1,1), delta_gYd_1loop);
   ASSERT_NON_ZERO(mq2(2,2), delta_gYd_1loop);
   ASSERT_NON_ZERO(mu2(0,0), delta_gYd_1loop);
   ASSERT_NON_ZERO(mu2(1,1), delta_gYd_1loop);
   ASSERT_NON_ZERO(mu2(2,2), delta_gYd_1loop);
   ASSERT_NON_ZERO(md2(0,0), delta_gYd_1loop);
   ASSERT_NON_ZERO(md2(1,1), delta_gYd_1loop);
   ASSERT_NON_ZERO(md2(2,2), delta_gYd_1loop);
   ASSERT_NON_ZERO(ml2(0,0), delta_gYd_1loop);
   ASSERT_NON_ZERO(ml2(1,1), delta_gYd_1loop);
   ASSERT_NON_ZERO(ml2(2,2), delta_gYd_1loop);
   ASSERT_NON_ZERO(me2(0,0), delta_gYd_1loop);
   ASSERT_NON_ZERO(me2(1,1), delta_gYd_1loop);
   ASSERT_NON_ZERO(me2(2,2), delta_gYd_1loop);
   ASSERT_NON_ZERO(mA2     , delta_gYd_1loop);
   ASSERT_NON_ZERO(tan_beta, delta_gYd_1loop);
   ASSERT_NON_ZERO(scale2  , delta_gYd_1loop);

   const double delta_gYd_1loop =
      (3*g22/16)*(-2+7*sin2_beta)
      +(3*g12/80)*(-44+7*sin2_beta)
      +(4*g12-9*(g12+5*g22)*sin2_beta)/40*std::log(mA2/scale2)
      +(g12/10)*(
         +std::log(ml2(0,0)/scale2)+std::log(ml2(1,1)/scale2)+std::log(ml2(2,2)/scale2)
         +2*std::log(me2(0,0)/scale2)+2*std::log(me2(1,1)/scale2)+2*std::log(me2(2,2)/scale2)
      )
      +(g12/30)*(
         +std::log(mq2(0,0)/scale2)+std::log(mq2(1,1)/scale2)+std::log(mq2(2,2)/scale2)
         +8*std::log(mu2(0,0)/scale2)+8*std::log(mu2(1,1)/scale2)+8*std::log(mu2(2,2)/scale2)
         +2*std::log(md2(0,0)/scale2)+2*std::log(md2(1,1)/scale2)+2*std::log(md2(2,2)/scale2)
      );

   return delta_gYd_1loop * oneLoop * g1 * std::sqrt(0.6) * cos_beta;
}

/**
 * Returns \f$\tilde{g}_{2u}\f$ at the one-loop level, all except the
 * first term on the rhs of Eq (15) of arXiv:1407.4081.
 *
 * @param parameters running parameters
 *
 * @return one-loop term on the rhs of Eq. (15) of arXiv:1407.4081
 */
double delta_g2u_1loop(const Parameters& parameters)
{
   const double g1 = parameters.g1;
   const double g2 = parameters.g2;
   const double g12 = sqr(g1);
   const double g22 = sqr(g2);
   const double tan_beta = parameters.tan_beta;
   const double beta = std::atan(tan_beta);
   const double cos2_beta = sqr(std::cos(beta));
   const double sin_beta = std::sin(beta);
   const double sin2_beta = sqr(sin_beta);
   const double gt = parameters.gt;
   const double gt2 = sqr(gt);
   const double scale2 = sqr(parameters.scale);
   const Eigen::Matrix<double,3,3> mq2(parameters.mq2);
   const Eigen::Matrix<double,3,3> mu2(parameters.mu2);
   const Eigen::Matrix<double,3,3> ml2(parameters.ml2);
   const double mA2 = sqr(parameters.mA);

   ASSERT_NON_ZERO(mq2(0,0), delta_g2u_1loop);
   ASSERT_NON_ZERO(mq2(1,1), delta_g2u_1loop);
   ASSERT_NON_ZERO(mq2(2,2), delta_g2u_1loop);
   ASSERT_NON_ZERO(mu2(0,0), delta_g2u_1loop);
   ASSERT_NON_ZERO(mu2(1,1), delta_g2u_1loop);
   ASSERT_NON_ZERO(mu2(2,2), delta_g2u_1loop);
   ASSERT_NON_ZERO(ml2(0,0), delta_g2u_1loop);
   ASSERT_NON_ZERO(ml2(1,1), delta_g2u_1loop);
   ASSERT_NON_ZERO(ml2(2,2), delta_g2u_1loop);
   ASSERT_NON_ZERO(mA2     , delta_g2u_1loop);
   ASSERT_NON_ZERO(tan_beta, delta_g2u_1loop);
   ASSERT_NON_ZERO(scale2  , delta_g2u_1loop);

   const double delta_g2u_1loop =
      -g22*(2./3.+11./16.*cos2_beta)
      +3*g12/80*(-2+7*cos2_beta)
      +9*gt2/(4*sin2_beta)
      +(20*g22+3*(-9*g12+35*g22)*cos2_beta)/120*std::log(mA2/scale2)
      +g22/6*(
         +std::log(ml2(0,0)/scale2)
         +std::log(ml2(1,1)/scale2)
         +std::log(ml2(2,2)/scale2)
      )
      +g22/2*(
         +std::log(mq2(0,0)/scale2)
         +std::log(mq2(1,1)/scale2)
         +std::log(mq2(2,2)/scale2)
      )
      -0.75*gt2/sin2_beta*(3*std::log(mq2(2,2)/scale2)-std::log(mu2(2,2)/scale2));

   return delta_g2u_1loop * oneLoop * g2 * sin_beta;
}

/**
 * Returns \f$\tilde{g}_{2d}\f$ at the one-loop level, all except the
 * first term on the rhs of Eq (16) of arXiv:1407.4081.
 *
 * @param parameters running parameters
 *
 * @return one-loop term on the rhs of Eq. (16) of arXiv:1407.4081
 */
double delta_g2d_1loop(const Parameters& parameters)
{
   const double g1 = parameters.g1;
   const double g2 = parameters.g2;
   const double g12 = sqr(g1);
   const double g22 = sqr(g2);
   const double tan_beta = parameters.tan_beta;
   const double beta = std::atan(tan_beta);
   const double cos_beta = std::cos(beta);
   const double sin_beta = std::sin(beta);
   const double sin2_beta = sqr(sin_beta);
   const double scale2 = sqr(parameters.scale);
   const Eigen::Matrix<double,3,3> mq2(parameters.mq2);
   const Eigen::Matrix<double,3,3> ml2(parameters.ml2);
   const double mA2 = sqr(parameters.mA);

   ASSERT_NON_ZERO(mq2(0,0), delta_g2d_1loop);
   ASSERT_NON_ZERO(mq2(1,1), delta_g2d_1loop);
   ASSERT_NON_ZERO(mq2(2,2), delta_g2d_1loop);
   ASSERT_NON_ZERO(ml2(0,0), delta_g2d_1loop);
   ASSERT_NON_ZERO(ml2(1,1), delta_g2d_1loop);
   ASSERT_NON_ZERO(ml2(2,2), delta_g2d_1loop);
   ASSERT_NON_ZERO(mA2     , delta_g2d_1loop);
   ASSERT_NON_ZERO(tan_beta, delta_g2d_1loop);
   ASSERT_NON_ZERO(scale2  , delta_g2d_1loop);

   const double delta_g2d_1loop =
      -g22*(2./3.+11./16.*sin2_beta)
      +(3*g12/80)*(-2+7*sin2_beta)
      +(g22/2)*std::log(mq2(0,0)/scale2)
      +(g22/2)*std::log(mq2(1,1)/scale2)
      +(g22/2)*std::log(mq2(2,2)/scale2)
      +(20*g22+3*(-9*g12+35*g22)*sin2_beta)/120*std::log(mA2/scale2)
      +(g22/6)*std::log(ml2(0,0)/scale2)
      +(g22/6)*std::log(ml2(1,1)/scale2)
      +(g22/6)*std::log(ml2(2,2)/scale2);

   return delta_g2d_1loop * oneLoop * g2 * cos_beta;
}

/**
 * Finite 1-loop threshold correction for top Yukawa coupling from
 * Higgsinos and gauginos (arXiv:1407.4081, Eq. (25)).
 *
 * @param scale renormalization scale
 * @param mu bilinear Higgsino coupling
 * @param gYu Higgs-Up-type-Higgsino-bino coupling
 * @param gYd Higgs-Down-type-Higgsino-bino coupling
 * @param g2u Higgs-Up-type-Higgsino-wino coupling
 * @param g2d Higgs-Down-type-Higgsino-wino coupling
 * @param m1 bino mass parameter
 * @param m2 wino mass parameter
 *
 * @return 1-loop threshold correction for the top Yukawa coupling
 * from Higgsinos and gauginos
 */
double delta_gt_1loop_chi(
   double scale, double mu, double gYu, double gYd,
   double g2u, double g2d, double m1, double m2)
{
   using namespace threshold_loop_functions;

   const double gYu2 = sqr(gYu);
   const double gYd2 = sqr(gYd);
   const double g2u2 = sqr(g2u);
   const double g2d2 = sqr(g2d);
   const double r1 = m1 / mu;
   const double r2 = m2 / mu;
   const double mu2 = sqr(mu);
   const double scale2 = sqr(scale);

   ASSERT_NON_ZERO(scale, delta_gt_1loop_chi);
   ASSERT_NON_ZERO(mu   , delta_gt_1loop_chi);
   ASSERT_NON_ZERO(m1   , delta_gt_1loop_chi);
   ASSERT_NON_ZERO(m2   , delta_gt_1loop_chi);

   const double delta_gt_chi =
      -1./6.*gYu*gYd*f(r1)
      -1./12.*(gYu2+gYd2)*(g(r1)+3*std::log(mu2/scale2))
      -0.5*g2u*g2d*f(r2)
      -0.25*(g2u2+g2d2)*(g(r2)+3*std::log(mu2/scale2));

   return oneLoop * delta_gt_chi;
}

namespace {
   /// arXiv:1502.06525, Eq. (44)
   double G(double x2, double q2)
   {
      if (is_zero(x2))
         return 0.;

      ASSERT_NON_ZERO(q2, G);

      return x2 * (std::log(x2/q2) - 1.);
   }

   /// arXiv:1502.06525, Eq. (A8)
   double G1(double m12, double m22, double q2)
   {
      ASSERT_NON_ZERO(q2, G1);

      if (is_equal(m12, m22)) {
         ASSERT_NON_ZERO(m12, G1);
         return std::log(m12/q2);
      }

      return (G(m12,q2) - G(m22,q2)) / (m12 - m22);
   }

   /// arXiv:1502.06525, Eq. (A13)
   double G2(double m12, double m22, double q2)
   {
      ASSERT_NON_ZERO(q2, G2);

      if (is_equal(m12, m22)) {
         ASSERT_NON_ZERO(m12, G2);
         return 2*m12*std::log(m12/q2) - m12;
      }

      return (m12*G(m12,q2) - m22*G(m22,q2)) / (m12 - m22);
   }
} // anonymous namespace

/**
 * Calculates the finite 1-loop threshold correction \f$\Delta m^2\f$
 * (note the sign!) for the Higgs mass parameter \f$m^2\f$ from
 * Higgsinos and gauginos (arXiv:1502.06525, Eq. (45)).  Only the
 * contribution from lines 4 and 5 are calculated.  The other lines
 * correspond to sfermions and the heavy Higgs and are ignored.
 *
 * @param parameters input parameters
 *
 * @return 1-loop threshold correction for Higgs mass parameter from
 * Higgsinos and gauginos
 */
double delta_m2_1loop_chi(const Parameters& parameters)
{
   const double g2 = parameters.g2;
   const double gY = parameters.g1 * std::sqrt(0.6);
   const double m1 = parameters.m1;
   const double m2 = parameters.m2;
   const double mu = parameters.mu;
   const double scale = parameters.scale;
   const double g22 = sqr(g2);
   const double gY2 = sqr(gY);
   const double m12 = sqr(m1);
   const double m22 = sqr(m2);
   const double mu2 = sqr(mu);
   const double scale2 = sqr(scale);

   const double minus_delta_m2 =
      - 6*g22*G2(m22,mu2,scale2)
      - 2*gY2*G2(m12,mu2,scale2)
      - 12*g22*m2*mu*G1(m22,mu2,scale2)
      - 4*gY2*m1*mu*G1(m12,mu2,scale2);

   return - minus_delta_m2 * 0.5 * oneLoop;
}

} // namespace splitmssm_thresholds
} // namespace flexiblesusy
