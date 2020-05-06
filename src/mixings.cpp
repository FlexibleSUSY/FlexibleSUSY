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

#include "mixings.hpp"
#include "logger.hpp"
#include "numerics2.hpp"
#include "config.h"

namespace flexiblesusy {

void convert_symmetric_fermion_mixings_to_slha(double& /*unused*/,
                                               Eigen::Matrix<double, 1, 1>& /*unused*/)
{
}

/**
 * @param m mass
 * @param z 1x1 mixing matrix
 */
void convert_symmetric_fermion_mixings_to_slha(double& m,
                                               Eigen::Matrix<std::complex<double>, 1, 1>& z)
{
   // check if 1st row contains non-zero imaginary parts
   if (!is_zero(std::abs(std::imag(z(0,0))))) {
      z(0,0) *= std::complex<double>(0.0,1.0);
      m *= -1;
#if defined(ENABLE_VERBOSE) || defined(ENABLE_DEBUG)
      if (!is_zero(std::abs(std::imag(z(0,0))))) {
         WARNING("Element (0,0) of the following fermion mixing matrix"
                 " contains entries which have non-zero real and imaginary"
                 " parts:\nZ = " << z);
      }
#endif
   }
}

void convert_symmetric_fermion_mixings_to_hk(double& /*unused*/,
                                             Eigen::Matrix<double, 1, 1>& /*unused*/)
{
}

/**
 * @param m mass
 * @param z 1x1 mixing matrix
 */
void convert_symmetric_fermion_mixings_to_hk(double& m,
                                             Eigen::Matrix<std::complex<double>, 1, 1>& z)
{
   if (m < 0.) {
      z(0,0) *= std::complex<double>(0.0,1.0);
      m *= -1;
   }
}

} // namespace flexiblesusy
