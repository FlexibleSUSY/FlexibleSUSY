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

#ifndef GSL_VECTOR_H
#define GSL_VECTOR_H

#include "gsl.hpp"
#include <gsl/gsl_vector.h>
#include <initializer_list>
#include <iosfwd>
#include <cstddef>

namespace flexiblesusy {

class GSL_vector {
public:
   GSL_vector() = default;
   GSL_vector(std::size_t);
   GSL_vector(const gsl_vector*);        ///< create copy of given gsl_vector
   GSL_vector(const GSL_vector&);
   GSL_vector(GSL_vector&&) noexcept;
   GSL_vector(std::initializer_list<double>);
   ~GSL_vector() noexcept;

   GSL_vector& operator=(const GSL_vector&);
   GSL_vector& operator=(GSL_vector&&) noexcept;
   double& operator[](std::size_t);      ///< element read/write access
   double operator[](std::size_t) const; ///< element read access
   double& operator()(std::size_t);      ///< element read/write access w/ range check
   double operator()(std::size_t) const; ///< element read access w/ range check

   bool empty() const noexcept;             ///< check if empty
   const gsl_vector* raw() const noexcept;  ///< get raw pointer
   gsl_vector* raw() noexcept;              ///< get raw pointer
   gsl::owner<gsl_vector>* release() noexcept; ///< release raw pointer from this object
   void set_all(double) noexcept;           ///< set all elemets to same value
   std::size_t size() const noexcept;       ///< number of elements

private:
   gsl::owner<gsl_vector>* vec{nullptr};    ///< raw gsl_vector

   void assign(const gsl_vector*);          ///< assign from gsl_vector
   void move_assign(GSL_vector&&) noexcept; ///< move assign
   void range_check(std::size_t) const;
};

double* begin(GSL_vector&); ///< iterator to begin of GSL_vector
double* end(GSL_vector&);   ///< iterator to end of GSL_vector

const double* cbegin(const GSL_vector&); ///< const iterator to begin of GSL_vector
const double* cend(const GSL_vector&);   ///< const iterator to end of GSL_vector

/// Returns true if GSL_vector contains only finite elements, false otherwise
bool is_finite(const GSL_vector&);

std::ostream& operator<<(std::ostream&, const GSL_vector&);

} // namespace flexiblesusy

#endif
