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

#ifndef EIGEN_TENSOR_H
#define EIGEN_TENSOR_H

#include <Eigen/Core>

// disable Eigen tensor module as it does not compiler with clang 5.0
#if 0 && EIGEN_VERSION_AT_LEAST(3,3,0)

#include <unsupported/Eigen/CXX11/Tensor>

namespace flexiblesusy {

template <class Scalar, int N, int M, int K>
class ZeroTensor3 : public Eigen::TensorFixedSize<Scalar, Eigen::Sizes<N,M,K> > {
public:
   using Base = Eigen::TensorFixedSize<Scalar, Eigen::Sizes<N,M,K>>;
   using Base::setZero;

   ZeroTensor3() : Base() {
      setZero();
   }
};

template <class Scalar, int N, int M, int K, int L>
class ZeroTensor4 : public Eigen::TensorFixedSize<Scalar, Eigen::Sizes<N,M,K,L> > {
public:
   using Base = Eigen::TensorFixedSize<Scalar, Eigen::Sizes<N,M,K,L>>;
   using Base::setZero;

   ZeroTensor4() : Base() {
      setZero();
   }
};

} // namespace flexiblesusy

#else

namespace flexiblesusy {

template <class Scalar, int N, int M, int K>
class ZeroTensor3;

template <class Scalar, int N, int M, int K, int L>
class ZeroTensor4;

} // namespace flexiblesusy

#endif

#endif
