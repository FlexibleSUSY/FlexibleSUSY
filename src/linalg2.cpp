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

#define INSTANTIATE_FS_HERMITIAN(N,Mtype,Ztype)                         \
   template                                                             \
   void fs_diagonalize_hermitian(const Eigen::Matrix<Mtype, N, N>&,     \
                                 Eigen::Array<double, N, 1>&,           \
                                 Eigen::Matrix<Ztype, N, N>&);

#define INSTANTIATE_FS_SYMMETRIC(N,Mtype,Ztype)                         \
   template                                                             \
   void fs_diagonalize_symmetric(const Eigen::Matrix<Mtype, N, N>&,     \
                                 Eigen::Array<double, N, 1>&,           \
                                 Eigen::Matrix<Ztype, N, N>&);

#define INSTANTIATE_FS_SVD(N,Mtype,Ztype)                               \
   template                                                             \
   void fs_svd(const Eigen::Matrix<Mtype, N, N>&,                       \
               Eigen::Array<double, N, 1>&,                             \
               Eigen::Matrix<Ztype, N, N>&,                             \
               Eigen::Matrix<Ztype, N, N>&);

INSTANTIATE_FS_HERMITIAN(2,double,double)
INSTANTIATE_FS_HERMITIAN(3,double,double)
INSTANTIATE_FS_HERMITIAN(6,double,double)

INSTANTIATE_FS_SYMMETRIC(4,double,std::complex<double>)

INSTANTIATE_FS_SVD(2,double,double)
INSTANTIATE_FS_SVD(2,double,std::complex<double>)
INSTANTIATE_FS_SVD(3,double,double)
INSTANTIATE_FS_SVD(3,double,std::complex<double>)

} // namespace flexiblesusy
