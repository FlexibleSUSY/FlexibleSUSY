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

#include "pmns.hpp"
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

int sign(double x) { return x >= 0.0 ? 1 : -1; }

} // anonymous namespace

void PMNS_parameters::reset_to_diagonal()
{
   theta_12 = 0.;
   theta_13 = 0.;
   theta_23 = 0.;
   delta    = 0.;
   alpha_1  = 0.;
   alpha_2  = 0.;
}

void PMNS_parameters::reset_to_observation()
{
   theta_12 = Electroweak_constants::PMNS_THETA12;
   theta_13 = Electroweak_constants::PMNS_THETA13;
   theta_23 = Electroweak_constants::PMNS_THETA23;
   delta    = Electroweak_constants::PMNS_DELTA;
   alpha_1  = Electroweak_constants::PMNS_ALPHA1;
   alpha_2  = Electroweak_constants::PMNS_ALPHA2;
}

Eigen::Matrix<double,3,3> PMNS_parameters::get_real_pmns() const
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

   Eigen::Matrix<double,3,3> pmns_matrix;
   pmns_matrix(0, 0) = c12 * c13;
   pmns_matrix(0, 1) = s12 * c13;
   pmns_matrix(0, 2) = pf * s13;
   pmns_matrix(1, 0) = -s12 * c23 - pf * c12 * s23 * s13;
   pmns_matrix(1, 1) = c12 * c23 - pf * s12 * s23 * s13;
   pmns_matrix(1, 2) = s23 * c13;
   pmns_matrix(2, 0) = s12 * s23 - pf * c12 * c23 * s13;
   pmns_matrix(2, 1) = -c12 * s23 - pf * s12 * c23 * s13;
   pmns_matrix(2, 2) = c23 * c13;

   return pmns_matrix;
}

Eigen::Matrix<std::complex<double>,3,3> PMNS_parameters::get_complex_pmns() const
{
   const std::complex<double> eID(std::polar(1.0, delta));
   const std::complex<double> eIAlpha1(std::polar(1.0, 0.5*alpha_1));
   const std::complex<double> eIAlpha2(std::polar(1.0, 0.5*alpha_2));
   const double s12 = std::sin(theta_12);
   const double s13 = std::sin(theta_13);
   const double s23 = std::sin(theta_23);
   const double c12 = std::cos(theta_12);
   const double c13 = std::cos(theta_13);
   const double c23 = std::cos(theta_23);

   Eigen::Matrix<std::complex<double>,3,3> pmns_matrix;
   pmns_matrix(0, 0) = c12 * c13 * eIAlpha1;
   pmns_matrix(0, 1) = s12 * c13 * eIAlpha2;
   pmns_matrix(0, 2) = s13 / eID;
   pmns_matrix(1, 0) = (-s12 * c23 - c12 * s23 * s13 * eID) * eIAlpha1;
   pmns_matrix(1, 1) = (c12 * c23 - s12 * s23 * s13 * eID) * eIAlpha2;
   pmns_matrix(1, 2) = s23 * c13;
   pmns_matrix(2, 0) = (s12 * s23 - c12 * c23 * s13 * eID) * eIAlpha1;
   pmns_matrix(2, 1) = (-c12 * s23 - s12 * c23 * s13 * eID) * eIAlpha2;
   pmns_matrix(2, 2) = c23 * c13;

   return pmns_matrix;
}

void PMNS_parameters::to_pdg_convention(Eigen::Matrix<double,3,3>& Vv,
                                        Eigen::Matrix<double,3,3>& Ve,
                                        Eigen::Matrix<double,3,3>& Ue)
{
   Eigen::Matrix<double,3,3> pmns(Ve*Vv.adjoint());
   to_pdg_convention(pmns, Vv, Ve, Ue);
}

void PMNS_parameters::to_pdg_convention(Eigen::Matrix<double,3,3>& pmns,
                                        Eigen::Matrix<double,3,3>& Vv,
                                        Eigen::Matrix<double,3,3>& Ve,
                                        Eigen::Matrix<double,3,3>& Ue)
{
   Eigen::Matrix<double,3,3> signs_E(Eigen::Matrix<double,3,3>::Identity());
   Eigen::Matrix<double,3,3> signs_V(Eigen::Matrix<double,3,3>::Identity());

   // make 33 element positive
   if (pmns(2, 2) < 0.) {
      signs_E(2, 2) = -1.;
      for (int j = 0; j < 3; ++j) {
         pmns(2, j) *= -1.;
      }
   }

   // make 23 element positive
   if (pmns(1, 2) < 0.) {
      signs_V(2, 2) = -1;
      signs_E(2, 2) *= -1;
      for (int j = 0; j < 3; ++j) {
         pmns(2, j) *= -1;
         pmns(j, 2) *= -1;
      }
   }

   Ve = signs_E * Ve;
   Ue = signs_E * Ue;
   Vv = signs_V * Vv;
}

void PMNS_parameters::to_pdg_convention(
   Eigen::Matrix<std::complex<double>,3,3>& Vv,
   Eigen::Matrix<std::complex<double>,3,3>& Ve,
   Eigen::Matrix<std::complex<double>,3,3>& Ue)
{
   Eigen::Matrix<std::complex<double>,3,3> pmns(Ve*Vv.adjoint());
   to_pdg_convention(pmns, Vv, Ve, Ue);
}

namespace {

template<typename T>
std::complex<T> phase(const std::complex<T>& z) noexcept
{
   T r = std::abs(z);
   return r == 0 ? 1 : z/r;
}

void calc_phase_factors(
   const Eigen::Matrix<std::complex<double>,3,3>& pmns,
   const std::complex<double>& p,
   std::complex<double>& o,
   Eigen::DiagonalMatrix<std::complex<double>,3>& l)
{
   o = std::conj(phase(p * pmns(0,2)));
   l.diagonal().bottomRightCorner<2,1>() = (o * pmns.bottomRightCorner<2,1>()).
      unaryExpr([] (const std::complex<double>& z) { return phase<double>(z); }).conjugate();
}

/// restrict sin or cos to interval [-1,1]
double sanitize_hypot(double sc) noexcept
{
   if (sc < -1.) { sc = -1.0; }
   if (sc > 1.) { sc = 1.0; }
   return sc;
}

} // anonymous namespace

void PMNS_parameters::to_pdg_convention(
   Eigen::Matrix<std::complex<double>,3,3>& pmns,
   Eigen::Matrix<std::complex<double>,3,3>& Vv,
   Eigen::Matrix<std::complex<double>,3,3>& Ve,
   Eigen::Matrix<std::complex<double>,3,3>& Ue)
{
   std::complex<double> o;
   Eigen::DiagonalMatrix<std::complex<double>,3> l(1,1,1);

   const double s13 = sanitize_hypot(std::abs(pmns(0,2)));
   const double c13_sq = 1. - sqr(s13);

   if (is_zero(c13_sq)) {
      o = std::conj(phase(pmns(0,2)));
      const auto rel_phase = std::sqrt(phase(pmns(1,0) * pmns(2,1)));
      const auto p = std::conj(rel_phase * rel_phase)
         * phase(pmns(1,0) * pmns(1,1));
      l.diagonal()[1] = std::conj(o * rel_phase);
      l.diagonal()[2] = std::conj(o * rel_phase) * p;
   } else {
      const double c13 = std::sqrt(c13_sq);
      const double s12 = sanitize_hypot(std::abs(pmns(0,1)) / c13);
      const double c12 = std::sqrt(1 - sqr(s12));
      const double s23 = sanitize_hypot(std::abs(pmns(1,2)) / c13);
      const double c23 = std::sqrt(1 - sqr(s23));
      const double jcp = std::imag(pmns(1,2) * std::conj(pmns(0,2))
                                   * pmns(0,1) * std::conj(pmns(1,1)));
      const double side1 = s12*s23;
      const double side2 = c12*c23*s13;
      const double cosdelta = sanitize_hypot(
         is_zero(jcp) ? 1 // delta is removable
         : (sqr(side1) + sqr(side2) - std::norm(pmns(2,0))) / (2*side1*side2));
      const double sindelta = std::sqrt(1 - sqr(cosdelta));
      const std::complex<double> p(cosdelta, sindelta); // Exp[I delta]
      calc_phase_factors(pmns, p, o, l);

      const auto a1 = phase(o*pmns(0,0));
      const auto a2 = phase(o*pmns(0,1));

      Eigen::Matrix<std::complex<double>,2,2> pmnsBL{
         (o * l * pmns).bottomLeftCorner<2,2>()};
      Eigen::Array<std::complex<double>,2,2> maj_phases;
      maj_phases << a1, a2, a1, a2;
      const Eigen::Array<double,2,2> imagBL{
         (maj_phases.conjugate() * pmnsBL.array()).imag()};
      if (!((imagBL <= 0).all() || (imagBL >= 0).all())) {
         calc_phase_factors(pmns, std::conj(p), o, l);
      }
   }

   Ve.transpose() *= l * o;
   Ue.transpose() *= (l * o).diagonal().conjugate().asDiagonal();
   pmns = Ve * Vv.adjoint();
}

} // namespace flexiblesusy
