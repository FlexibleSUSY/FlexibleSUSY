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

#ifndef GM2_CONSTANTS_H
#define GM2_CONSTANTS_H

/* fine-structure constant in the Thompson limit (Q = 0) from PDG (2019) */
#define GM2CALC_ALPHA_EM_THOMPSON (1.0/137.035999084)

/* quark and lepton contributions to the on-shell renormalized photon
   vacuum polarization */
#define GM2CALC_DELTA_ALPHA_EM_MZ (                       \
      + 0.0314979  /* leptonic */                         \
      - 0.00007180 /* top */                              \
      + 0.027611   /* hadronic, arXiv:1802.02995 */       \
)

/* fine-structure constant at Q = MZ */
#define GM2CALC_ALPHA_EM_MZ (                                           \
      GM2CALC_ALPHA_EM_THOMPSON / (1.0 - GM2CALC_DELTA_ALPHA_EM_MZ)     \
)

/* strong coupling alpha_s at Q = MZ */
#define GM2CALC_ALPHA_S_MZ 0.1184

/* W boson pole mass */
#define GM2CALC_MW 80.385

/* Z boson pole mass */
#define GM2CALC_MZ 91.1876

/* Top quark pole mass */
#define GM2CALC_MT 173.34

/* SM MS-bar bottom quark mass mb at the scale mb */
#define GM2CALC_MBMB 4.18

/* Electron pole mass */
#define GM2CALC_ME 0.000510998928

/* Muon pole mass */
#define GM2CALC_MM 0.1056583715

/* Tau lepton pole mass */
#define GM2CALC_ML 1.777

#endif
