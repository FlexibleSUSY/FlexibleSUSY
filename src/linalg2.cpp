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

#include "linalg2.hpp"

namespace flexiblesusy {

template
void fs_diagonalize_hermitian(const Eigen::Matrix<double, 2, 2>&,
                              Eigen::Array<double, 2, 1>&,
                              Eigen::Matrix<double, 2, 2>&);

template
void fs_diagonalize_hermitian(const Eigen::Matrix<double, 3, 3>&,
                              Eigen::Array<double, 3, 1>&,
                              Eigen::Matrix<double, 3, 3>&);

template
void fs_diagonalize_hermitian(const Eigen::Matrix<double, 6, 6>&,
                              Eigen::Array<double, 6, 1>&,
                              Eigen::Matrix<double, 6, 6>&);

template
void fs_diagonalize_symmetric(const Eigen::Matrix<double, 4, 4>&,
                              Eigen::Array<double, 4, 1>&,
                              Eigen::Matrix<std::complex<double>, 4, 4>&);

template
void fs_svd(const Eigen::Matrix<double, 2, 2>&,
            Eigen::Array<double, 2, 1>&,
            Eigen::Matrix<double, 2, 2>&,
            Eigen::Matrix<double, 2, 2>&);

template
void fs_svd(const Eigen::Matrix<double, 2, 2>&,
            Eigen::Array<double, 2, 1>&,
            Eigen::Matrix<std::complex<double>, 2, 2>&,
            Eigen::Matrix<std::complex<double>, 2, 2>&);

template
void fs_svd(const Eigen::Matrix<double, 3, 3>&,
            Eigen::Array<double, 3, 1>&,
            Eigen::Matrix<double, 3, 3>&,
            Eigen::Matrix<double, 3, 3>&);

template
void fs_svd(const Eigen::Matrix<double, 3, 3>&,
            Eigen::Array<double, 3, 1>&,
            Eigen::Matrix<std::complex<double>, 3, 3>&,
            Eigen::Matrix<std::complex<double>, 3, 3>&);

} // namespace flexiblesusy
