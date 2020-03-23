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

#ifndef GM2_1LOOP_H
#define GM2_1LOOP_H

#include "MSSMNoFV_onshell.h"

/**
 * @file gm2_1loop.h
 * @brief contains declarations of C interface functions for 1-loop calculation
 *
 * This file contains the declarations for the C interface functions
 * used to calculate \f$a_\mu\f$ at the 1-loop level.
 */

#ifdef __cplusplus
extern "C" {
#endif

/** calculates full 1-loop SUSY contributions to (g-2) in the MSSMNoFV (w/ tan(beta) resummation) */
double gm2calc_mssmnofv_calculate_amu_1loop(const MSSMNoFV_onshell*);

/** calculates full 1-loop SUSY contributions to (g-2) in the MSSMNoFV (no tan(beta) resummation) */
double gm2calc_mssmnofv_calculate_amu_1loop_non_tan_beta_resummed(const MSSMNoFV_onshell*);

/* === routines for individual 1-loop contributions === */

/** 1-loop neutralino contribution */
double gm2calc_mssmnofv_amu1LChi0(const MSSMNoFV_onshell*);

/** 1-loop chargino contribution */
double gm2calc_mssmnofv_amu1LChipm(const MSSMNoFV_onshell*);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
