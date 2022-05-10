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

/// real dilogarithm
double Li2(double) noexcept;

/// real dilogarithm
long double Li2(long double) noexcept;

/// complex dilogarithm
std::complex<double> Li2(const std::complex<double>&) noexcept;

/// complex dilogarithm
std::complex<long double> Li2(const std::complex<long double>&) noexcept;

/// Clausen function Cl_2(x)
double Cl2(double) noexcept;

/// Clausen function Cl_2(x)
long double Cl2(long double) noexcept;

} // namespace flexiblesusy

#undef noexcept

#endif
