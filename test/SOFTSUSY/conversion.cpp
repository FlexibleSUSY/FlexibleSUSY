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

#include "conversion.hpp"

namespace flexiblesusy {

Eigen::MatrixXd ToEigen(const softsusy::DoubleMatrix& m)
{
   const int cols = m.displayCols(), rows = m.displayRows();
   Eigen::MatrixXd eig(rows,cols);
   for (int i = 1; i <= rows; ++i)
      for (int j = 1; j <= rows; ++j)
         eig(i-1,j-1) = m(i,j);
   return eig;
}

softsusy::DoubleMatrix ToSoftsusy(const Eigen::MatrixXd& m)
{
   const int cols = m.cols(), rows = m.rows();
   softsusy::DoubleMatrix soft(rows,cols);
   for (int i = 1; i <= rows; ++i)
      for (int j = 1; j <= rows; ++j)
         soft(i,j) = m(i-1,j-1);
   return soft;
}

softsusy::DoubleVector ToSoftsusy(const Eigen::VectorXd& m)
{
   const int cols = m.rows();
   softsusy::DoubleVector soft(cols);
   for (int j = 1; j <= cols; ++j)
      soft(j) = m(j-1);
   return soft;
}

Eigen::ArrayXd ToEigenArray(const softsusy::DoubleVector& v)
{
   Eigen::ArrayXd a(v.size());
   for (int i = v.displayStart(); i <= v.displayEnd(); i++)
      a(i - v.displayStart()) = v(i);
   return a;
}

Eigen::ArrayXd ToEigenArray(double v)
{
   Eigen::ArrayXd a(1);
   a(0) = v;
   return a;
}

std::valarray<double> ToValarray(const softsusy::DoubleVector& v)
{
   std::valarray<double> a(v.size());
   for (int i = v.displayStart(); i <= v.displayEnd(); i++)
      a[i - v.displayStart()] = v(i);
   return a;
}

std::valarray<double> ToValarray(double v)
{
   std::valarray<double> a(1);
   a[0] = v;
   return a;
}

softsusy::DoubleVector ToDoubleVector(const Eigen::ArrayXd& a)
{
   softsusy::DoubleVector v(a.rows());
   for (int i = 0; i < a.rows(); i++)
      v(i + 1) = a(i);
   return v;
}

Eigen::MatrixXd ToEigenMatrix(const softsusy::DoubleMatrix& m)
{
   const int r = m.displayRows();
   const int c = m.displayCols();
   Eigen::MatrixXd result(r,c);

   for (int i = 1; i <= r; i++)
      for (int k = 1; k <= c; k++)
         result(i-1, k-1) = m(i,k);

   return result;
}

Eigen::MatrixXcd ToEigenMatrix(const softsusy::ComplexMatrix& m)
{
   const int r = m.displayRows();
   const int c = m.displayCols();
   Eigen::MatrixXcd result(r,c);

   for (int i = 1; i <= r; i++)
      for (int k = 1; k <= c; k++)
         result(i-1, k-1) = m(i,k);

   return result;
}

softsusy::DoubleMatrix ToDoubleMatrix(const Eigen::MatrixXd& m)
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
