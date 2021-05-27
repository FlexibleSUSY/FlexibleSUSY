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

#ifndef FS_LI4_H
#define FS_LI4_H

#include <complex>

namespace flexiblesusy {

/// complex polylogarithm with n=4
std::complex<double> Li4(const std::complex<double>&) noexcept;

/// complex polylogarithm with n=4 with long double precision
std::complex<long double> Li4(const std::complex<long double>&) noexcept;

} // namespace flexiblesusy

#endif
