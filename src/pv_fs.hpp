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

/**
 * @file PV.hpp
 *
 * @brief Declaration of real Passarino-Veltman loop functions with
 * squared arguments.
 */

#pragma once

namespace flexiblesusy {

/// A0 Passarino-Veltman function
double a0(double m2, double q2) noexcept;
/// B0 Passarino-Veltman function
double b0(double p2, double m12, double m22, double q2) noexcept;
/// B0(s,x,x,q2) Passarino-Veltman function
double b0xx(double p2, double m2, double q2) noexcept;
/// derivative of B0 Passarino-Veltman function w.r.t. p^2, for p^2 = 0
double d1_b0(double m12, double m22) noexcept;
/// F0 Passarino-Veltman function
double f0(double p2, double m12, double m22, double q2) noexcept;
/// F0 Passarino-Veltman function
double g0(double p2, double m12, double m22, double q2) noexcept;

} // namespace flexiblesusy
