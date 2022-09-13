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

#ifndef MIXINGS_H
#define MIXINGS_H

#include "numerics2.hpp"
#include <Eigen/Core>

namespace flexiblesusy {

void convert_symmetric_fermion_mixings_to_slha(double&,
                                               Eigen::Matrix<double, 1, 1>&);

template<int N>
void convert_symmetric_fermion_mixings_to_slha(Eigen::Array<double, N, 1>&,
                                               Eigen::Matrix<double, N, N>&)
{
}

void convert_symmetric_fermion_mixings_to_slha(double&,
                                               Eigen::Matrix<std::complex<double>, 1, 1>&);

/**
 * Converts the given vector of masses and the corresponding (complex)
 * mixing matrix to SLHA convention: Matrix rows with non-zero
 * imaginary parts are multiplied by i and the corresponding mass
 * eigenvalue is multiplied by -1.  As a result the mixing matrix will
 * be real and the mass eigenvalues might be positive or negative.  It
 * is assumed that these mixings result from diagonalizing a symmetric
 * fermion mass matrix in the convention of Haber and Kane,
 * Phys. Rept. 117 (1985) 75-263.  This conversion makes sense only if
 * the original symmetric mass matrix is real-valued.
 *
 * @param m vector of masses
 * @param z mixing matrix
 */
template<int N>
void convert_symmetric_fermion_mixings_to_slha(Eigen::Array<double, N, 1>& m,
                                               Eigen::Matrix<std::complex<double>, N, N>& z)
{
   for (int i = 0; i < N; i++) {
      // check if i'th row contains non-zero imaginary parts
      if (!is_zero(z.row(i).imag().cwiseAbs().maxCoeff())) {
         z.row(i) *= std::complex<double>(0.0,1.0);
         m(i) *= -1;
      }
   }
}

void convert_symmetric_fermion_mixings_to_hk(double&,
                                             Eigen::Matrix<double, 1, 1>&);

template<int N>
void convert_symmetric_fermion_mixings_to_hk(Eigen::Array<double, N, 1>&,
                                             Eigen::Matrix<double, N, N>&)
{
}

void convert_symmetric_fermion_mixings_to_hk(double&,
                                             Eigen::Matrix<std::complex<double>, 1, 1>&);

/**
 * Converts the given vector of masses and the corresponding (real)
 * mixing matrix to Haber-Kane convention (Phys. Rept. 117 (1985)
 * 75-263): Masses are positive and mixing matrices can be complex.
 *
 * @param m vector of masses
 * @param z mixing matrix
 */
template<int N>
void convert_symmetric_fermion_mixings_to_hk(Eigen::Array<double, N, 1>& m,
                                             Eigen::Matrix<std::complex<double>, N, N>& z)
{
   for (int i = 0; i < N; i++) {
      if (m(i) < 0.) {
         z.row(i) *= std::complex<double>(0.0,1.0);
         m(i) *= -1;
      }
   }
}

} // namespace flexiblesusy

#endif
