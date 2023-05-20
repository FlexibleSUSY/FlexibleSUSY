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

#ifndef SFERMIONS_H
#define SFERMIONS_H

namespace flexiblesusy {
namespace sfermions {

enum Sparticles {
   up = 0,
   down = 1,
   neutrino = 2,
   electron = 3,
   NUMBER_OF_MSSM_SPARTICLES
};

extern const double Isospin[NUMBER_OF_MSSM_SPARTICLES];
extern const double Hypercharge_left[NUMBER_OF_MSSM_SPARTICLES];
extern const double Hypercharge_right[NUMBER_OF_MSSM_SPARTICLES];

/**
 * parameters needed to fill 2 x 2 sfermion mass matrix
 */
struct Mass_data {
   double ml2{};      ///< soft mass of left-handed sfermion
   double mr2{};      ///< soft mass of right-handed sfermion
   double yf{};       ///< Yukawa coupling
   double vd{}, vu{}; ///< Higgs VEVs
   double gY{}, g2{}; ///< gauge couplings (not GUT normalized)
   double Tyf{};      ///< trilinear coupling
   double mu{};       ///< Superpotential parameter
   double T3{};       ///< weak isospin
   double Yl{};       ///< Hypercharge of left-handed sfermion
   double Yr{};       ///< Hypercharge of right-handed sfermion
};

/// returns mixing angle and sets (squared) mass eigenvalues
double diagonalize_sfermions_2x2(const Mass_data&, double&, double&);

} // namespace sfermions
} // namespace flexiblesusy

#endif
