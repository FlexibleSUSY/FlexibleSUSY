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

/** \file numerics.h
   - Project:     SOFTSUSY
   - Author:      Ben Allanach, Alexander Voigt
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

   \brief loop functions
*/

#ifndef NUMERICS_H
#define NUMERICS_H

namespace softsusy {

double a0(double m, double q) noexcept;
double b0(double p, double m1, double m2, double q) noexcept;
double b1(double p, double m1, double m2, double q) noexcept;
double b22(double p,  double m1, double m2, double q) noexcept;
double c0(double m1, double m2, double m3) noexcept;
double c00(double m1, double m2, double m3, double q) noexcept;
double d27(double m1, double m2, double m3, double m4) noexcept;
double d0(double m1, double m2, double m3, double m4) noexcept;
double ffn(double p, double m1, double m2, double q) noexcept;
double gfn(double p, double m1, double m2, double q) noexcept;
double hfn(double p, double m1, double m2, double q) noexcept;
double b22bar(double p, double m1, double m2, double q) noexcept;

double d1_b0(double p2, double m2a, double m2b) noexcept;

} // namespace softsusy

#endif
