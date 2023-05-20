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

#ifndef CONVERSION_H
#define CONVERSION_H

#include "linalg.h"
#include <Eigen/Core>

namespace flexiblesusy {

Eigen::ArrayXd ToEigenArray(const softsusy::DoubleVector&);
Eigen::ArrayXd ToEigenArray(double);
std::valarray<double> ToValarray(const softsusy::DoubleVector&);
std::valarray<double> ToValarray(double);
Eigen::MatrixXd ToEigenMatrix(const softsusy::DoubleMatrix&);
Eigen::MatrixXcd ToEigenMatrix(const softsusy::ComplexMatrix&);

template<class Derived>
softsusy::DoubleVector ToDoubleVector(const Eigen::ArrayBase<Derived>& a)
{
   softsusy::DoubleVector v(a.rows());
   for (int i = 0; i < a.rows(); i++)
      v(i + 1) = a(i);
   return v;
}

template<class Derived>
softsusy::ComplexMatrix ToComplexMatrix(const Eigen::MatrixBase<Derived>& m)
{
   const int r = m.rows();
   const int c = m.cols();
   softsusy::ComplexMatrix result(r,c);

   for (int i = 0; i < r; i++)
      for (int k = 0; k < c; k++)
         result(i+1, k+1) = m(i,k);

   return result;
}

template<class Derived>
softsusy::DoubleMatrix ToDoubleMatrix(const Eigen::MatrixBase<Derived>& m)
{
   const int r = m.rows();
   const int c = m.cols();
   softsusy::DoubleMatrix result(r,c);

   for (int i = 0; i < r; i++)
      for (int k = 0; k < c; k++)
         result(i+1, k+1) = m(i,k);

   return result;
}

} // namespace flexiblesusy

#endif
