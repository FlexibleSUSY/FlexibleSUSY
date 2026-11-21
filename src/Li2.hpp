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

#ifndef LI2_H
#define LI2_H

#include <complex>

namespace flexiblesusy {

/// real polylogarithm with n=2 (dilogarithm) with single precision
float Li2(float) noexcept;

/// real polylogarithm with n=2 (dilogarithm) with double precision
double Li2(double) noexcept;

/// real polylogarithm with n=2 (dilogarithm) with long double precision
long double Li2(long double) noexcept;

/// complex polylogarithm with n=2 (dilogarithm) with single precision
std::complex<float> Li2(const std::complex<float>&) noexcept;

/// complex polylogarithm with n=2 (dilogarithm) with double precision
std::complex<double> Li2(const std::complex<double>&) noexcept;

/// complex polylogarithm with n=2 (dilogarithm) with long double precision
std::complex<long double> Li2(const std::complex<long double>&) noexcept;

} // namespace flexiblesusy

#endif
