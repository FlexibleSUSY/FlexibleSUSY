/* ====================================================================
 * This file is part of GM2Calc.
 *
 * GM2Calc is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * GM2Calc is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GM2Calc.  If not, see
 * <http://www.gnu.org/licenses/>.
 * ==================================================================== */

#ifndef GM2_UNCERTAINTY_H
#define GM2_UNCERTAINTY_H

#include "MSSMNoFV_onshell.h"

#ifdef __cplusplus
extern "C" {
#endif

/** calculates uncertainty for amu(0-loop) in the MSSMNoFV w/ tan(beta) resummation */
double gm2calc_mssmnofv_calculate_uncertainty_amu_0loop(const MSSMNoFV_onshell*);

/** calculates uncertainty for amu(1-loop) in the MSSMNoFV w/ tan(beta) resummation */
double gm2calc_mssmnofv_calculate_uncertainty_amu_1loop(const MSSMNoFV_onshell*);

/** calculates uncertainty for amu(2-loop) in the MSSMNoFV */
double gm2calc_mssmnofv_calculate_uncertainty_amu_2loop(const MSSMNoFV_onshell*);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
