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

#ifndef TEST_COMPLEX_EQUALITY_H
#define TEST_COMPLEX_EQUALITY_H

#include <boost/test/unit_test.hpp>
#include <complex>

#define TEST_COMPLEX_EQUALITY(a, b)                     \
   do {                                                 \
      const std::complex<double> aa(a), bb(b);          \
      BOOST_TEST(std::real(aa) == std::real(bb));       \
      BOOST_TEST(std::imag(aa) == std::imag(bb));       \
   } while (false)

#define TEST_COMPLEX_CLOSE_FRACTION(a, b, eps) do {                     \
      BOOST_CHECK_CLOSE_FRACTION(std::real(a), std::real(b), eps);      \
      BOOST_CHECK_CLOSE_FRACTION(std::imag(a), std::imag(b), eps);      \
   } while (false)

#endif
