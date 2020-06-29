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

#include "ckm.hpp"
#include "error.hpp"
#include "ew_input.hpp"

#include <cmath>
#include <limits>

namespace flexiblesusy {

namespace {

bool is_zero(double x) noexcept
{
   return std::abs(x) <= std::numeric_limits<double>::epsilon();
}

double sqr(double x) { return x*x; }

double pow3(double x) { return x*x*x; }

double pow4(double x) { return sqr(sqr(x)); }

int sign(double x) { return x >= 0.0 ? 1 : -1; }

} // anonymous namespace

void CKM_parameters::reset_to_diagonal()
{
   theta_12 = 0.;
   theta_13 = 0.;
   theta_23 = 0.;
   delta    = 0.;
}

void CKM_parameters::reset_to_observation()
{
   theta_12 = Electroweak_constants::CKM_THETA12;
   theta_13 = Electroweak_constants::CKM_THETA13;
   theta_23 = Electroweak_constants::CKM_THETA23;
   delta    = Electroweak_constants::CKM_DELTA;
}

/**
 * Calculates V_CKM angles from Wolfenstein parameters (see
 * hep-ph/0406184)
 */
void CKM_parameters::set_from_wolfenstein(double lambdaW, double aCkm,
                                          double rhobar, double etabar)
{
   if (std::abs(lambdaW) > 1.) throw SetupError("Error: Wolfenstein lambda out of range!");
   if (std::abs(aCkm)    > 1.) throw SetupError("Error: Wolfenstein A parameter out of range!");
   if (std::abs(rhobar)  > 1.) throw SetupError("Error: Wolfenstein rho-bar parameter out of range!");
   if (std::abs(etabar)  > 1.) throw SetupError("Error: Wolfenstein eta-bar parameter out of range!");

   theta_12 = std::asin(lambdaW);
   theta_23 = std::asin(aCkm * sqr(lambdaW));

   const double lambdaW3 = pow3(lambdaW);
   const double lambdaW4 = pow4(lambdaW);

   const std::complex<double> rpe(rhobar, etabar);
   const std::complex<double> V13conj = aCkm * lambdaW3 * rpe
      * std::sqrt(1.0 - sqr(aCkm) * lambdaW4) /
      std::sqrt(1.0 - sqr(lambdaW)) / (1.0 - sqr(aCkm) * lambdaW4 * rpe);

   if (std::isfinite(std::real(V13conj)) && std::isfinite(std::imag(V13conj))) {
      theta_13 = std::asin(std::abs(V13conj));
      delta = std::arg(V13conj);
   }
}

/**
 * Calculates Wolfenstein parameters from V_CKM angles (see
 * hep-ph/0406184)
 */
void CKM_parameters::get_wolfenstein(double& lambdaW, double& aCkm,
                                     double& rhobar, double& etabar) const
{
   const double sin_12 = std::sin(theta_12);
   const double sin_13 = std::sin(theta_13);
   const double sin_23 = std::sin(theta_23);

   // Eq. (11.4) from PDG
   lambdaW  = sin_12;
   aCkm     = sin_23 / sqr(lambdaW);

   if (!std::isfinite(aCkm)) {
      aCkm = 0.0;
   }

   const double c = std::sqrt((1.0 - sqr(sin_23)) / (1.0 - sqr(lambdaW)));
   const std::complex<double> eid(std::polar(1.0, delta));
   const std::complex<double> r(sin_13 * eid /
      (c * aCkm * pow3(lambdaW) + sin_13 * eid * sqr(sin_23)));

   rhobar = std::isfinite(std::real(r)) ? std::real(r) : 0.;
   etabar = std::isfinite(std::imag(r)) ? std::imag(r) : 0.;
}

Eigen::Matrix<double,3,3> CKM_parameters::get_real_ckm() const
{
   const std::complex<double> eID(std::polar(1.0, delta));
   const double s12 = std::sin(theta_12);
   const double s13 = std::sin(theta_13);
   const double s23 = std::sin(theta_23);
   const double c12 = std::cos(theta_12);
   const double c13 = std::cos(theta_13);
   const double c23 = std::cos(theta_23);

   // set phase factor e^(i delta) to +1 or -1 depending on the sign
   // of its real part
   const int pf = sign(std::real(eID));

   Eigen::Matrix<double,3,3> ckm_matrix;
   ckm_matrix(0, 0) = c12 * c13;
   ckm_matrix(0, 1) = s12 * c13;
   ckm_matrix(0, 2) = pf * s13;
   ckm_matrix(1, 0) = -s12 * c23 - pf * c12 * s23 * s13;
   ckm_matrix(1, 1) = c12 * c23 - pf * s12 * s23 * s13;
   ckm_matrix(1, 2) = s23 * c13;
   ckm_matrix(2, 0) = s12 * s23 - pf * c12 * c23 * s13;
   ckm_matrix(2, 1) = -c12 * s23 - pf * s12 * c23 * s13;
   ckm_matrix(2, 2) = c23 * c13;

   return ckm_matrix;
}

Eigen::Matrix<std::complex<double>,3,3> CKM_parameters::get_complex_ckm() const
{
   const std::complex<double> eID(std::polar(1.0, delta));
   const double s12 = std::sin(theta_12);
   const double s13 = std::sin(theta_13);
   const double s23 = std::sin(theta_23);
   const double c12 = std::cos(theta_12);
   const double c13 = std::cos(theta_13);
   const double c23 = std::cos(theta_23);

   Eigen::Matrix<std::complex<double>,3,3> ckm_matrix;
   ckm_matrix(0, 0) = c12 * c13;
   ckm_matrix(0, 1) = s12 * c13;
   ckm_matrix(0, 2) = s13 / eID;
   ckm_matrix(1, 0) = -s12 * c23 - c12 * s23 * s13 * eID;
   ckm_matrix(1, 1) = c12 * c23 - s12 * s23 * s13 * eID;
   ckm_matrix(1, 2) = s23 * c13;
   ckm_matrix(2, 0) = s12 * s23 - c12 * c23 * s13 * eID;
   ckm_matrix(2, 1) = -c12 * s23 - s12 * c23 * s13 * eID;
   ckm_matrix(2, 2) = c23 * c13;

   return ckm_matrix;
}

void CKM_parameters::to_pdg_convention(Eigen::Matrix<double,3,3>& Vu,
                                       Eigen::Matrix<double,3,3>& Vd,
                                       Eigen::Matrix<double,3,3>& Uu,
                                       Eigen::Matrix<double,3,3>& Ud)
{
   Eigen::Matrix<double,3,3> ckm(Vu*Vd.adjoint());
   to_pdg_convention(ckm, Vu, Vd, Uu, Ud);
}

void CKM_parameters::to_pdg_convention(Eigen::Matrix<double,3,3>& ckm,
                                       Eigen::Matrix<double,3,3>& Vu,
                                       Eigen::Matrix<double,3,3>& Vd,
                                       Eigen::Matrix<double,3,3>& Uu,
                                       Eigen::Matrix<double,3,3>& Ud)
{
   Eigen::Matrix<double,3,3> signs_U(Eigen::Matrix<double,3,3>::Identity()),
      signs_D(Eigen::Matrix<double,3,3>::Identity());

   // make diagonal elements positive
   for (int i = 0; i < 3; i++) {
      if (ckm(i, i) < 0.) {
         signs_U(i, i) = -1.;
         for (int j = 0; j < 3; j++) {
            ckm(i, j) *= -1;
         }
      }
   }

   // make 12 element positive while keeping diagonal elements positive
   if (ckm(0, 1) < 0.) {
      signs_D(1, 1) = -1.;
      signs_U(1, 1) *= -1;
      for (int j = 0; j < 3; j++) {
         ckm(1, j) *= -1;
         ckm(j, 1) *= -1;
      }
   }

   // make 23 element positive while keeping diagonal elements positive
   if (ckm(1, 2) < 0.) {
      signs_D(2, 2) = -1.;
      signs_U(2, 2) *= -1;
      for (int j = 0; j < 3; j++) {
         ckm(2, j) *= -1;
         ckm(j, 2) *= -1;
      }
   }

   Vu = signs_U * Vu;
   Uu = signs_U * Uu;
   Vd = signs_D * Vd;
   Ud = signs_D * Ud;
}

void CKM_parameters::to_pdg_convention(Eigen::Matrix<std::complex<double>,3,3>& Vu,
                                       Eigen::Matrix<std::complex<double>,3,3>& Vd,
                                       Eigen::Matrix<std::complex<double>,3,3>& Uu,
                                       Eigen::Matrix<std::complex<double>,3,3>& Ud)
{
   Eigen::Matrix<std::complex<double>,3,3> ckm(Vu*Vd.adjoint());
   to_pdg_convention(ckm, Vu, Vd, Uu, Ud);
}

namespace {

template<typename T>
std::complex<T> phase(const std::complex<T>& z) noexcept
{
   T r = std::abs(z);
   return r == 0 ? 1 : z/r;
}

void calc_phase_factors
(const Eigen::Matrix<std::complex<double>,3,3>& ckm,
 const std::complex<double>& p,
 std::complex<double>& o,
 Eigen::DiagonalMatrix<std::complex<double>,3>& l,
 Eigen::DiagonalMatrix<std::complex<double>,3>& r)
{
   o = std::conj(phase(p * ckm(0,2)));
   l.diagonal().bottomRightCorner<2,1>() = (o * ckm.bottomRightCorner<2,1>()).
      unaryExpr([] (const std::complex<double>& z) { return phase<double>(z); }).conjugate();
   r.diagonal().topLeftCorner<2,1>() = (o * ckm.topLeftCorner<1,2>()).
      unaryExpr([] (const std::complex<double>& z) { return phase<double>(z); }).adjoint();
}

/// restrict sin or cos to interval [-1,1]
double sanitize_hypot(double sc) noexcept
{
   if (sc < -1.) { sc = -1.0; }
   if (sc > 1.) { sc = 1.0; }
   return sc;
}

} // anonymous namespace

void CKM_parameters::to_pdg_convention(Eigen::Matrix<std::complex<double>,3,3>& ckm,
                                       Eigen::Matrix<std::complex<double>,3,3>& Vu,
                                       Eigen::Matrix<std::complex<double>,3,3>& Vd,
                                       Eigen::Matrix<std::complex<double>,3,3>& Uu,
                                       Eigen::Matrix<std::complex<double>,3,3>& Ud)
{
   std::complex<double> o;
   Eigen::DiagonalMatrix<std::complex<double>,3> l(1,1,1), r(1,1,1);

   const double s13 = sanitize_hypot(std::abs(ckm(0,2)));
   const double c13_sq = 1. - sqr(s13);

   if (is_zero(c13_sq)) {
      o = std::conj(phase(ckm(0,2)));
      r.diagonal().block<2,1>(0,0) = (o * ckm.block<1,2>(1,0)).unaryExpr(
         [] (const std::complex<double>& z) { return phase<double>(z); }).adjoint();
      l.diagonal()[2] = std::conj(phase(o * r.diagonal()[1] * ckm(2,1)));
   } else {
      const double c13 = std::sqrt(c13_sq);
      const double s12 = sanitize_hypot(std::abs(ckm(0,1)) / c13);
      const double c12 = std::sqrt(1 - sqr(s12));
      const double s23 = sanitize_hypot(std::abs(ckm(1,2)) / c13);
      const double c23 = std::sqrt(1 - sqr(s23));
      const double side1 = s12*s23;
      const double side2 = c12*c23*s13;
      const double cosdelta = sanitize_hypot(
         is_zero(side1*side2) ? 1 // delta is removable
         : (sqr(side1)+sqr(side2)-std::norm(ckm(2,0))) / (2*side1*side2));
      const double sindelta = std::sqrt(1 - sqr(cosdelta));
      const std::complex<double> p(cosdelta, sindelta); // Exp[I delta]
      calc_phase_factors(ckm, p, o, l, r);
      const Eigen::Array<double,2,2>
         imagBL{(o * l * ckm * r).bottomLeftCorner<2,2>().imag()};
      if (!((imagBL <= 0).all() || (imagBL >= 0).all())) {
         calc_phase_factors(ckm, std::conj(p), o, l, r);
      }
   }

   Vu.transpose() *= l * o;
   Vd.transpose() *= r.diagonal().conjugate().asDiagonal();
   Uu.transpose() *= (l * o).diagonal().conjugate().asDiagonal();
   Ud.transpose() *= r;
   ckm = Vu * Vd.adjoint();
}

} // namespace flexiblesusy
