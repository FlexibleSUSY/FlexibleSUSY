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

#ifndef GM2_2LOOP_H
#define GM2_2LOOP_H

#include "MSSMNoFV_onshell.h"

/**
 * @file gm2_2loop.h
 * @brief contains declarations of C interface functions for 2-loop calculation
 *
 * This file contains the declarations for the C interface functions
 * used to calculate \f$a_\mu\f$ at the 2-loop level.
 */

#ifdef __cplusplus
extern "C" {
#endif

/** calculates best 2-loop SUSY contributions to a_mu in the MSSMNoFV (with tan(beta) resummation) */
double gm2calc_mssmnofv_calculate_amu_2loop(const MSSMNoFV_onshell*);

/** calculates best 2-loop SUSY contributions to a_mu in the MSSMNoFV (no tan(beta) resummation) */
double gm2calc_mssmnofv_calculate_amu_2loop_non_tan_beta_resummed(const MSSMNoFV_onshell*);

/* === 2-loop fermion/sfermion approximations === */

double gm2calc_mssmnofv_amu2LFSfapprox(const MSSMNoFV_onshell*);
double gm2calc_mssmnofv_amu2LFSfapprox_non_tan_beta_resummed(const MSSMNoFV_onshell*);

/* === photonic 2-loop corrections === */

double gm2calc_mssmnofv_amu2LChipmPhotonic(const MSSMNoFV_onshell*);
double gm2calc_mssmnofv_amu2LChi0Photonic(const MSSMNoFV_onshell*);

/* === SUSY 2L(a) diagrams === */

double gm2calc_mssmnofv_amu2LaSferm(const MSSMNoFV_onshell*);
double gm2calc_mssmnofv_amu2LaCha(const MSSMNoFV_onshell*);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
