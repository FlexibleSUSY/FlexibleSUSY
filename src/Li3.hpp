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

#ifndef LI3_H
#define LI3_H

#include <complex>

namespace flexiblesusy {

/// real polylogarithm with n=3 (trilogarithm)
double Li3(double) noexcept;

/// complex polylogarithm with n=3 (trilogarithm)
std::complex<double> Li3(const std::complex<double>&) noexcept;

/// complex polylogarithm with n=3 (trilogarithm) with long double precision
std::complex<long double> Li3(const std::complex<long double>&) noexcept;

} // namespace flexiblesusy

#endif
