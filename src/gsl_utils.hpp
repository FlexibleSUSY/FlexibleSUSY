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

/// Returns true if GSL vector contains only finite elements, false otherwise
bool is_finite(const gsl_vector*);

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
Eigen::Matrix<double,Size,1> to_eigen_vector(const gsl_vector* v)
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

} // namespace flexiblesusy

#endif
