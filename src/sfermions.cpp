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
#include "wrappers.hpp"
#include "logger.hpp"

namespace flexiblesusy {
namespace sfermions {

static const double oneOverRoot2 = 0.7071067811865475; // 1/sqrt(2.)

const std::array<double, NUMBER_OF_MSSM_SPARTICLES> Isospin = {
   0.5, -0.5, 0.5, -0.5
};

const std::array<double, NUMBER_OF_MSSM_SPARTICLES> Hypercharge_left = {
   1./3., 1./3., -1., -1.
};

const std::array<double, NUMBER_OF_MSSM_SPARTICLES> Hypercharge_right = {
   -4./3., 2./3., 0., 2.
};

/**
 * Obtains 2 x 2 mass matrix using input parameters in first argument 
 * and diagonalises it.  Fills the second argument with the eigenvalues
 * and returns the mixing angle.
 */ 
double diagonalize_sfermions_2x2(const Mass_data& pars,
                                 Eigen::Array<double,2,1>& msf)
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
   const double vev2   = 0.25 * (Sqr(vd) - Sqr(vu));
   Eigen::Matrix<double,2,2> mass_matrix;
   /// fill sfermion phi in mass matix in basis (phi_L phi_R)
   if (Sign(T3) > 0) {
      mass_matrix(0,0) = ml2 + 0.5 * AbsSqr(yf) * Sqr(vu)
         + (T3 * Sqr(g2) - 0.5 * Yl * Sqr(gY)) * vev2;
      mass_matrix(0,1) = oneOverRoot2 * (vu*Conj(Tyf) - vd*Conj(yf)*mu);
      mass_matrix(1,0) = Conj(mass_matrix(0,1));
      mass_matrix(1,1) = mr2 + 0.5 * AbsSqr(yf) * Sqr(vu)
         - 0.5 * Yr * Sqr(gY) * vev2;
   } else {
      mass_matrix(0,0) = ml2 + 0.5 * AbsSqr(yf) * Sqr(vd)
         + (T3 * Sqr(g2) - 0.5 * Yl * Sqr(gY)) * vev2;
      mass_matrix(0,1) = oneOverRoot2 * (vd*Conj(Tyf) - vu*Conj(yf)*mu);
      mass_matrix(1,0) = Conj(mass_matrix(0,1));
      mass_matrix(1,1) = mr2 + 0.5 * AbsSqr(yf) * Sqr(vd)
         - 0.5 * Yr * Sqr(gY) * vev2;
   }

   Eigen::Matrix<double, 2, 2> Zf;
   diagonalize_hermitian(mass_matrix, msf, Zf);

#ifdef ENABLE_VERBOSE
   if (msf.minCoeff() < 0.)
      WARNING("diagonalize_sfermions_2x2: sfermion tachyon");
#endif

   msf = AbsSqrt(msf);

   double theta;

   if (Sign(Zf(0,0)) == Sign(Zf(1,1))) {
      theta = ArcCos(Abs(Zf(0,0)));
   } else {
      theta = ArcCos(Abs(Zf(0,1)));
      Zf.col(0).swap(Zf.col(1));
      std::swap(msf(0), msf(1));
   }

   theta = Sign(mass_matrix(0,1) / (mass_matrix(0,0) - mass_matrix(1,1)))
      * Abs(theta);

   return theta;
}

} // namespace sfermions
} // namespace flexiblesusy
