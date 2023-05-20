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

#include "gsl_vector.hpp"
#include "error.hpp"

#include <cmath>
#include <iostream>
#include <string>
#include <utility>

namespace flexiblesusy {

GSL_vector::GSL_vector(std::size_t size)
{
   if (size == 0u) {
      return;
   }

   vec = gsl_vector_calloc(size);

   if (vec == nullptr) {
      throw OutOfMemoryError(
         "Allocation of GSL_vector of size " + std::to_string(size)
         + " failed.");
   }
}

/**
 * Creates a new GSL_vector by copying the content of the given
 * gsl_vector .
 *
 * @param v gsl_vector to copy elements from
 */
GSL_vector::GSL_vector(const gsl_vector* v)
{
   assign(v);
}

GSL_vector::GSL_vector(const GSL_vector& other)
{
   vec = gsl_vector_alloc(other.size());

   if (vec == nullptr) {
      throw OutOfMemoryError(
         "Allocation of GSL_vector of size " + std::to_string(other.size())
         + " failed.");
   }

   gsl_vector_memcpy(vec, other.vec);
}

GSL_vector::GSL_vector(GSL_vector&& other) noexcept
{
   move_assign(std::move(other));
}

GSL_vector::GSL_vector(std::initializer_list<double> list)
{
   if (list.size() == 0) {
      return;
   }

   vec = gsl_vector_alloc(list.size());

   if (vec == nullptr) {
      throw OutOfMemoryError(
         "Allocation of GSL_vector of size " + std::to_string(list.size())
         + " failed.");
   }

   std::copy(list.begin(), list.end(), gsl_vector_ptr(vec, 0));
}

GSL_vector::~GSL_vector() noexcept
{
   gsl_vector_free(vec);
}

/**
 * Creates new GSL_vector with the content of the given pointer.
 *
 * @param other gsl_vector whose elements are copied
 */
void GSL_vector::assign(const gsl_vector* other)
{
   if (other == nullptr) {
      gsl_vector_free(vec);
      vec = nullptr;
      return;
   }

   // avoid free and alloc if other has same size
   if (size() != other->size) {
      gsl_vector_free(vec);
      vec = gsl_vector_alloc(other->size);

      if (vec == nullptr) {
         throw OutOfMemoryError(
            "Allocation of GSL_vector of size " + std::to_string(other->size)
            + " failed.");
      }
   }

   gsl_vector_memcpy(vec, other);
}

bool GSL_vector::empty() const noexcept
{
   return size() == 0;
}

GSL_vector& GSL_vector::operator=(const GSL_vector& rhs)
{
   if (this != &rhs) {
      assign(rhs.vec);
   }

   return *this;
}

GSL_vector& GSL_vector::operator=(GSL_vector&& rhs) noexcept
{
   if (this != &rhs) {
      gsl_vector_free(vec);
      move_assign(std::move(rhs));
   }

   return *this;
}

double& GSL_vector::operator[](std::size_t n)
{
   return *gsl_vector_ptr(vec, n);
}

double GSL_vector::operator[](std::size_t n) const
{
   return gsl_vector_get(vec, n);
}

double& GSL_vector::operator()(std::size_t n)
{
   range_check(n);
   return operator[](n);
}

double GSL_vector::operator()(std::size_t n) const
{
   range_check(n);
   return operator[](n);
}

std::size_t GSL_vector::size() const noexcept
{
   if (vec == nullptr) {
      return 0;
   }
   return vec->size;
}

/**
 * Releases the encapsulated gsl_vector from this object.  The pointer
 * to the gsl_vector is returned.  After this method has been called,
 * this object will no longer delete the the gsl_vector .
 *
 * @return pointer to gsl_vector
 */
gsl::owner<gsl_vector>* GSL_vector::release() noexcept
{
   gsl_vector* raw = vec;
   vec = nullptr;
   return raw;
}

const gsl_vector* GSL_vector::raw() const noexcept
{
   return vec;
}

gsl_vector* GSL_vector::raw() noexcept
{
   return vec;
}

void GSL_vector::set_all(double value) noexcept
{
   if (vec != nullptr) {
      gsl_vector_set_all(vec, value);
   }
}

std::ostream& operator<<(std::ostream& ostr, const GSL_vector& vec)
{
   ostr << "(";

   for (std::size_t i = 0; i < vec.size(); i++) {
      ostr << vec[i];
      if (i < vec.size() - 1) {
         ostr << ", ";
      }
   }

   ostr << ")";

   return ostr;
}

void GSL_vector::move_assign(GSL_vector&& other) noexcept
{
   vec = other.vec;
   other.vec = nullptr;
}

void GSL_vector::range_check(std::size_t n) const
{
   if (vec == nullptr) {
      throw OutOfBoundsError(
         "GSL_vector::operator[]: index " + std::to_string(n)
         + " out of range for vector of size 0.");
   }

   if (n >= size()) {
      throw OutOfBoundsError(
         "GSL_vector::operator[]: index " + std::to_string(n)
         + " out of range for vector of size " + std::to_string(size()) + ".");
   }
}

double* begin(GSL_vector& v)
{
   if (!v.empty()) {
      return gsl_vector_ptr(v.raw(), 0);
   }

   return nullptr;
}

double* end(GSL_vector& v)
{
   if (!v.empty()) {
      double* last = gsl_vector_ptr(v.raw(), v.size() - 1);
      return ++last;
   }

   return nullptr;
}

const double* cbegin(const GSL_vector& v)
{
   if (!v.empty()) {
      return gsl_vector_const_ptr(v.raw(), 0);
   }

   return nullptr;
}

const double* cend(const GSL_vector& v)
{
   if (!v.empty()) {
      const double* last = gsl_vector_const_ptr(v.raw(), v.size() - 1);
      return ++last;
   }

   return nullptr;
}

/**
 * Returns true if GSL_vector contains only finite elements (neither
 * nan nor inf), false otherwise.
 *
 * @param v GSL vector
 * @return true if vector contains only finite elements, false otherwise.
 */
bool is_finite(const GSL_vector& v)
{
   for (std::size_t i = 0; i < v.size(); i++) {
      if (!std::isfinite(v[i])) {
         return false;
      }
   }

   return true;
}

} // namespace flexiblesusy
