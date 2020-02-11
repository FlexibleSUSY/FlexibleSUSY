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

#ifndef GSL_UTILS_H
#define GSL_UTILS_H

#include "error.hpp"
#include "gsl_vector.hpp"

#include <gsl/gsl_vector.h>
#include <Eigen/Core>

namespace flexiblesusy {

/// Returns true if GSL_vector contains only finite elements, false otherwise
bool is_finite(const GSL_vector&);
/// Returns true if GSL vector contains only finite elements, false otherwise
bool is_finite(const gsl_vector*);
Eigen::ArrayXd to_eigen_array(const gsl_vector*);
Eigen::ArrayXd to_eigen_array(const GSL_vector&);
Eigen::VectorXd to_eigen_vector(const gsl_vector*);
Eigen::VectorXd to_eigen_vector(const GSL_vector&);
GSL_vector to_GSL_vector(const gsl_vector*);

template <typename Derived>
GSL_vector to_GSL_vector(const Eigen::DenseBase<Derived>& v)
{
   using Index_t = typename Derived::Index;
   GSL_vector v2(v.rows());

   for (Index_t i = 0; i < v.rows(); i++) {
      v2[i] = v(i);
   }

   return v2;
}

template <int Size>
Eigen::Matrix<double,Size,1> to_eigen_vector_fixed(const gsl_vector* v)
{
   if (Size != v->size) {
      throw OutOfBoundsError("Size of GSL_vector does not match size of Eigen vector.");
   }

   using Result_t = Eigen::Matrix<double,Size,1>;
   using Index_t = typename Result_t::Index;
   Result_t result;

   for (Index_t i = 0; i < Size; i++) {
      result(i) = gsl_vector_get(v, i);
   }

   return result;
}

template <int Size>
Eigen::Matrix<double,Size,1> to_eigen_vector_fixed(const GSL_vector& v)
{
   if (Size != v.size()) {
      throw OutOfBoundsError(
         "Size of GSL_vector does not match size of Eigen vector.");
   }

   using Result_t = Eigen::Matrix<double,Size,1>;
   using Index_t = typename Result_t::Index;
   Result_t result;

   for (Index_t i = 0; i < Size; i++) {
      result(i) = v[i];
   }

   return result;
}

/**
 * Copies values from an Eigen array/matrix to a GSL vector.
 *
 * @param src Eigen array/matrix
 * @param dst GSL vector
 */
template <typename Derived>
void copy(const Eigen::DenseBase<Derived>& src, gsl_vector* dst)
{
   using Index_t = typename Derived::Index;
   const auto dim = src.rows();

   if (dst == nullptr) {
      throw SetupError("Destination GSL_vector is Null.");
   }

   if (static_cast<decltype(dst->size)>(dim) != dst->size) {
      throw OutOfBoundsError("Size of destination GSL_vector does not match "
                             "size of source Eigen vector.");
   }

   for (Index_t i = 0; i < dim; i++) {
      gsl_vector_set(dst, i, src(i));
   }
}

} // namespace flexiblesusy

#endif
