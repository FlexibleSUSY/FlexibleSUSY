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

#ifndef EIGEN_UTILS_H
#define EIGEN_UTILS_H

#include <Eigen/Core>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <limits>

namespace flexiblesusy {

template <typename Derived>
int closest_index(double mass, const Eigen::ArrayBase<Derived>& v)
{
   int pos;
   typename Derived::PlainObject tmp;
   tmp.setConstant(mass);

   (v - tmp).abs().minCoeff(&pos);

   return pos;
}

template <class BinaryOp, class Derived>
Derived binary_map(
   const Eigen::DenseBase<Derived>& a, const Eigen::DenseBase<Derived>& b, BinaryOp op)
{
   typename Derived::PlainObject result(a.rows(), b.cols());

   assert(a.rows() == b.rows());
   assert(a.cols() == b.cols());

   for (int k = 0; k < a.cols(); k++)
      for (int i = 0; i < a.rows(); i++)
         result(i,k) = op(a(i,k), b(i,k));

   return result;
}

/**
 * Divides a by b element wise.  If the quotient is not finite, it is
 * set to zero.
 *
 * @param a numerator
 * @param b denominator
 */
template <class Derived>
Derived div_safe(
   const Eigen::DenseBase<Derived>& a, const Eigen::DenseBase<Derived>& b)
{
   using Scalar = typename Derived::Scalar;

   return binary_map(a, b, [](Scalar x, Scalar y){
         const Scalar q = x / y;
         return std::isfinite(q) ? q : Scalar{};
      });
}

/**
 * Returns true when two objects are element-wise equal up to a
 * relative precision \a eps.
 *
 * @param a first object
 * @param b second object
 * @param eps maximum relative difference
 *
 * @return true when both objects are equal, false otherwise
 */
template <class Derived>
bool is_equal_rel(const Eigen::PlainObjectBase<Derived>& a,
                  const Eigen::PlainObjectBase<Derived>& b, double eps)
{
   if (a.rows() != b.rows() || a.cols() != b.cols()) {
      return false;
   }

   for (decltype(a.size()) i = 0; i < a.size(); i++) {
      const auto ai = a.data()[i];
      const auto bi = b.data()[i];
      const auto max = std::max(std::abs(ai), std::abs(bi));
      if (std::abs(ai - bi) >= eps*(1 + max)) {
         return false;
      }
   }

   return true;
}

/**
 * The element of v, which is closest to mass, is moved to the
 * position idx.
 *
 * @param idx new index of the mass eigenvalue
 * @param mass mass to compare against
 * @param v vector of masses
 * @param z corresponding mixing matrix
 */
template <typename DerivedArray, typename DerivedMatrix>
void move_goldstone_to(int idx, double mass, Eigen::ArrayBase<DerivedArray>& v,
                       Eigen::MatrixBase<DerivedMatrix>& z)
{
   int pos = closest_index(mass, v);
   if (pos == idx)
      return;

   const int sign = (idx - pos) < 0 ? -1 : 1;
   int steps = std::abs(idx - pos);

   // now we shuffle the states
   while (steps--) {
      const int new_pos = pos + sign;
      v.row(new_pos).swap(v.row(pos));
      z.row(new_pos).swap(z.row(pos));
      pos = new_pos;
   }
}

/**
 * Normalize each element of the given real matrix to be within the
 * interval [min, max].  Values < min are set to min.  Values > max
 * are set to max.
 *
 * @param m matrix
 * @param min minimum
 * @param max maximum
 */
template <int M, int N>
void normalize_to_interval(Eigen::Matrix<double,M,N>& m, double min = -1., double max = 1.)
{
   auto data = m.data();
   const auto size = m.size();

   for (int i = 0; i < size; i++) {
      if (data[i] < min)
         data[i] = min;
      else if (data[i] > max)
         data[i] = max;
   }
}

/**
 * Normalize each element of the given complex matrix to have a
 * magnitude within the interval [0, max].  If the magnitude of a
 * matrix element is > max, then the magnitude is set to max.  The
 * phase angles are not modified.
 *
 * @param m matrix
 * @param max_mag maximum magnitude
 */
template <int M, int N>
void normalize_to_interval(Eigen::Matrix<std::complex<double>,M,N>& m, double max_mag = 1.)
{
   auto data = m.data();
   const auto size = m.size();

   for (int i = 0; i < size; i++) {
      if (std::abs(data[i]) > max_mag)
         data[i] = std::polar(max_mag, std::arg(data[i]));
   }
}

namespace {
template <class T>
struct Is_not_finite {
   bool operator()(T x) const noexcept { return !std::isfinite(x); }
};
} // anonymous namespace

/**
 * Returns all elements from src, which are not close to the elements
 * in cmp.  The returned vector will have the length (src.size() -
 * cmp.size()).
 *
 * @param src source vector
 * @param cmp vector with elements to compare against
 * @return vector with elements of src not close to cmp
 */
template<class Real, int Nsrc, int Ncmp>
Eigen::Array<Real,Nsrc - Ncmp,1> remove_if_equal(
   const Eigen::Array<Real,Nsrc,1>& src,
   const Eigen::Array<Real,Ncmp,1>& cmp)
{
   Eigen::Array<Real,Nsrc,1> non_equal(src);
   Eigen::Array<Real,Nsrc - Ncmp,1> dst;

   for (int i = 0; i < Ncmp; i++) {
      const int idx = closest_index(cmp(i), non_equal);
      non_equal(idx) = std::numeric_limits<double>::infinity();
   }

   std::remove_copy_if(non_equal.data(), non_equal.data() + Nsrc,
                       dst.data(), Is_not_finite<Real>());

   return dst;
}

/**
 * @brief reorders vector v according to ordering in vector v2
 * @param v vector with elementes to be reordered
 * @param v2 vector with reference ordering
 */
template<class Real, int N>
void reorder_vector(
   Eigen::Array<Real,N,1>& v,
   const Eigen::Array<Real,N,1>& v2)
{
   Eigen::PermutationMatrix<N> p;
   p.setIdentity();
   std::sort(p.indices().data(), p.indices().data() + p.indices().size(),
             [&v2] (int i, int j) { return v2[i] < v2[j]; });

#if EIGEN_VERSION_AT_LEAST(3,1,4)
   v.matrix().transpose() *= p.inverse();
#else
   Eigen::Map<Eigen::Matrix<Real,N,1> >(v.data()).transpose() *= p.inverse();
#endif
}

/**
 * @brief reorders vector v according to ordering of diagonal elements in mass_matrix
 * @param v vector with elementes to be reordered
 * @param matrix matrix with diagonal elements with reference ordering
 */
template<class Derived>
void reorder_vector(
   Eigen::Array<double,Eigen::MatrixBase<Derived>::RowsAtCompileTime,1>& v,
   const Eigen::MatrixBase<Derived>& matrix)
{
   reorder_vector(v, matrix.diagonal().array().eval());
}

/// sorts an Eigen array
template<int N>
void sort(Eigen::Array<double, N, 1>& v)
{
   std::sort(v.data(), v.data() + v.size(),
             [] (double a, double b) { return std::abs(a) < std::abs(b); });
}

} // namespace flexiblesusy

#endif
