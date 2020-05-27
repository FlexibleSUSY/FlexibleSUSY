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

#ifndef NUMERICS_HPP
#define NUMERICS_HPP

#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <cstddef>
#include <cstdlib>
#include <type_traits>

namespace flexiblesusy {

namespace detail {

template <class T> struct is_complex : public std::false_type {};
template <class T> struct is_complex<const T > : public is_complex<T>{};
template <class T> struct is_complex<volatile const T > : public is_complex<T>{};
template <class T> struct is_complex<volatile T > : public is_complex<T>{};
template <class T> struct is_complex<std::complex<T> > : public std::true_type{};

} // namespace detail


/// compares a number for being close to zero
template <typename T>
typename std::enable_if<std::is_unsigned<T>::value, bool>::type
is_zero(T a, T prec = std::numeric_limits<T>::epsilon()) noexcept
{
   return a <= prec;
}


/// compares a number for being close to zero
template <typename T>
typename std::enable_if<!std::is_unsigned<T>::value &&
                        !detail::is_complex<T>::value, bool>::type
is_zero(T a, T prec = std::numeric_limits<T>::epsilon()) noexcept
{
   return std::abs(a) <= prec;
}


/// compares a complex number for being close to zero
template <typename T>
typename std::enable_if<detail::is_complex<T>::value, bool>::type
is_zero(T a,
        typename T::value_type prec
        = std::numeric_limits<typename T::value_type>::epsilon()) noexcept
{
   return is_zero(a.real(), prec) && is_zero(a.imag(), prec);
}


/// compares two numbers for (absolute) equality
template <typename T>
bool is_equal(T a, T b, T prec = std::numeric_limits<T>::epsilon()) noexcept
{
   return is_zero(a - b, prec);
}


/// compares two complex numbers for (absolute) equality
template <typename T>
bool is_equal(std::complex<T> a, std::complex<T> b,
              T prec = std::numeric_limits<T>::epsilon()) noexcept
{
   return is_equal(a.real(), b.real(), prec) &&
          is_equal(a.imag(), b.imag(), prec);
}


/// compares two numbers for relative equality
template <typename T>
typename std::enable_if<!std::is_unsigned<T>::value, bool>::type
is_equal_fraction(T a, T b, T prec = std::numeric_limits<T>::epsilon()) noexcept
{
   const T max = std::max(std::abs(a), std::abs(b));
   return is_zero(a - b, max*prec);
}


/// compares two numbers for relative equality, treating numbers with
/// small differences as equal
template <typename T>
typename std::enable_if<!std::is_unsigned<T>::value, bool>::type
is_equal_rel(T a, T b, T prec = std::numeric_limits<T>::epsilon()) noexcept
{
   if (is_equal(a, b, std::numeric_limits<T>::epsilon())) {
      return true;
   }

   const T min = std::min(std::abs(a), std::abs(b));

   if (min < std::numeric_limits<T>::epsilon()) {
      return is_equal(a, b, prec);
   }

   const T max = std::max(std::abs(a), std::abs(b));

   return is_equal(a, b, prec*max);
}


/// compares two numbers for relative equality, treating numbers with
/// small differences as equal
template <typename T>
typename std::enable_if<std::is_unsigned<T>::value, bool>::type
is_equal_rel(T a, T b, T prec = std::numeric_limits<T>::epsilon()) noexcept
{
   using ST = typename std::make_signed<T>::type;
   const auto sa = static_cast<ST>(a);
   const auto sb = static_cast<ST>(b);
   const auto sprec = static_cast<ST>(prec);

   return is_equal_rel(sa, sb, sprec);
}


/// returns true if all array elements are finite
bool is_finite(const double*, long) noexcept;


/// returns true if all array elements are finite
template <typename T, std::size_t N>
bool is_finite(const T (&v)[N]) noexcept
{
   return is_finite(&v[0], N);
}


/// returns true if all array elements are finite
template <typename T, std::size_t N>
bool is_finite(const std::array<T, N>& v) noexcept
{
   return is_finite(&v[0], N);
}


/// fast implementation of complex logarithm
template <class T>
std::complex<T> fast_log(const std::complex<T>& z) noexcept
{
   const T rz = std::real(z);
   const T iz = std::imag(z);

   return std::complex<T>(0.5*std::log(rz*rz + iz*iz), std::atan2(iz, rz));
}

} // namespace flexiblesusy

#endif
