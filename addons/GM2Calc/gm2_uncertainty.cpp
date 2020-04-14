// ====================================================================
// This file is part of GM2Calc.
//
// GM2Calc is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// GM2Calc is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GM2Calc.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#include "gm2calc/gm2_uncertainty.hpp"
#include "gm2calc/gm2_1loop.hpp"
#include "gm2calc/gm2_2loop.hpp"

#include <cmath>

/**
 * \file gm2_uncertainty.cpp
 *
 * Contains functions necessary to calculate theory uncertainty
 * associated with amu.
 */

namespace gm2calc {

/**
 * Calculates uncertainty associated with amu(0-loop) including
 * tan(beta) resummation.
 *
 * The estimated uncertainty is the magnitude amu(1-loop) (including
 * tan(beta) resummation).
 *
 * @param model model parameters
 *
 * @return uncertainty for amu(0-loop) w/ tan(beta) resummation
 */
double calculate_uncertainty_amu_0loop(const MSSMNoFV_onshell& model)
{
   const double amu_1L = calculate_amu_1loop(model);

   return std::abs(amu_1L);
}

/**
 * Calculates uncertainty associated with amu(1-loop) including
 * tan(beta) resummation.
 *
 * The estimated uncertainty is the sum of magnitude amu(2-loop)
 * (including tan(beta) resummation) and the 2-loop uncertainty,
 * calculated by calculate_uncertainty_amu_2loop().
 *
 * @param model model parameters
 *
 * @return uncertainty for amu(1-loop) w/ tan(beta) resummation
 */
double calculate_uncertainty_amu_1loop(const MSSMNoFV_onshell& model)
{
   const double amu_2L = calculate_amu_2loop(model);
   const double delta_amu_2L = calculate_uncertainty_amu_2loop(model);

   return std::abs(amu_2L) + std::abs(delta_amu_2L);
}

/**
 * Calculates uncertainty associated with amu(2-loop) using Eq (4).
 *
 * Eq. (4) takes into account the unknown two-loop contributions and
 * the employed approximation for the 2L(a) contributions.
 *
 * @param model model parameters
 *
 * @return uncertainty for amu(2-loop)
 */
double calculate_uncertainty_amu_2loop(const MSSMNoFV_onshell& model)
{
   const double amu_2La_Cha = amu2LaCha(model);
   const double amu_2La_Sferm = amu2LaSferm(model);

   return 2.3e-10 + 0.3 * (std::abs(amu_2La_Cha) + std::abs(amu_2La_Sferm));
}

} // namespace gm2calc
