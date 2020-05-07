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
 * @file sfermions.cpp
 * @brief finding mass eigenstates and mixing of sfermions in absence of 
 *        family mixing, where we have a 2 by 2 mass matrix.
 */

#include "sfermions.hpp"
#include "linalg2.hpp"

#include <cmath>

#include <Eigen/Core>

namespace flexiblesusy {

namespace {

double conj(double x) noexcept { return x; }

double sqr(double x) noexcept { return x*x; }

int sign(double x) noexcept { return x >= 0.0 ? 1 : -1; }

} // anonymous namespace

namespace sfermions {

static const double oneOverRoot2 = 0.7071067811865475; // 1/sqrt(2.)

const double Isospin[NUMBER_OF_MSSM_SPARTICLES] = {
   0.5, -0.5, 0.5, -0.5
};

const double Hypercharge_left[NUMBER_OF_MSSM_SPARTICLES] = {
   1./3., 1./3., -1., -1.
};

const double Hypercharge_right[NUMBER_OF_MSSM_SPARTICLES] = {
   -4./3., 2./3., 0., 2.
};

/**
 * Obtains 2 x 2 mass matrix using input parameters in first argument
 * and diagonalises it.  Fills the second argument with the (squared)
 * mass eigenvalues and returns the mixing angle.
 */
double diagonalize_sfermions_2x2(const Mass_data& pars, double& msf1, double& msf2)
{
   const double ml2    = pars.ml2;
   const double mr2    = pars.mr2;
   const double yf     = pars.yf;
   const double vd     = pars.vd;
   const double vu     = pars.vu;
   const double gY     = pars.gY;
   const double g2     = pars.g2;
   const double Tyf    = pars.Tyf;
   const double mu     = pars.mu;
   const double T3     = pars.T3;
   const double Yl     = pars.Yl;
   const double Yr     = pars.Yr;
   const double vev2   = 0.25 * (sqr(vd) - sqr(vu));
   Eigen::Matrix<double,2,2> mass_matrix;

   // fill sfermion phi in mass matix in basis (phi_L phi_R)
   if (sign(T3) > 0) {
      mass_matrix(0,0) = ml2 + 0.5 * sqr(yf) * sqr(vu)
         + (T3 * sqr(g2) - 0.5 * Yl * sqr(gY)) * vev2;
      mass_matrix(0,1) = oneOverRoot2 * (vu*conj(Tyf) - vd*conj(yf)*mu);
      mass_matrix(1,0) = conj(mass_matrix(0,1));
      mass_matrix(1,1) = mr2 + 0.5 * sqr(yf) * sqr(vu)
         - 0.5 * Yr * sqr(gY) * vev2;
   } else {
      mass_matrix(0,0) = ml2 + 0.5 * sqr(yf) * sqr(vd)
         + (T3 * sqr(g2) - 0.5 * Yl * sqr(gY)) * vev2;
      mass_matrix(0,1) = oneOverRoot2 * (vd*conj(Tyf) - vu*conj(yf)*mu);
      mass_matrix(1,0) = conj(mass_matrix(0,1));
      mass_matrix(1,1) = mr2 + 0.5 * sqr(yf) * sqr(vd)
         - 0.5 * Yr * sqr(gY) * vev2;
   }

   Eigen::Array<double,2,1> msf;
   Eigen::Matrix<double, 2, 2> Zf;
   diagonalize_hermitian(mass_matrix, msf, Zf);

   double theta;

   if (sign(Zf(0,0)) == sign(Zf(1,1))) {
      theta = std::acos(std::abs(Zf(0,0)));
   } else {
      theta = std::acos(std::abs(Zf(0,1)));
      Zf.col(0).swap(Zf.col(1));
      std::swap(msf(0), msf(1));
   }

   msf1 = msf(0);
   msf2 = msf(1);
   theta = sign(mass_matrix(0,1) / (mass_matrix(0,0) - mass_matrix(1,1)))
      * std::abs(theta);

   return theta;
}

} // namespace sfermions
} // namespace flexiblesusy
