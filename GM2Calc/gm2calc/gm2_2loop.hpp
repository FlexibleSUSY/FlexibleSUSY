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

#ifndef GM2_2LOOP_HPP
#define GM2_2LOOP_HPP

namespace gm2calc {

class MSSMNoFV_onshell;

/// calculates best 2-loop SUSY contributions to a_mu in the MSSM (with tan(beta) resummation)
double calculate_amu_2loop(const MSSMNoFV_onshell&);

/// calculates best 2-loop SUSY contributions to a_mu in the MSSM (no tan(beta) resummation)
double calculate_amu_2loop_non_tan_beta_resummed(const MSSMNoFV_onshell&);

// === routines for individual 2-loop contributions ===

/// 2-loop fermion/sfermion contribution (approximation)
double amu2LFSfapprox(const MSSMNoFV_onshell&);

/// 2-loop fermion/sfermion contribution (approximation) w/o tan(beta) resummation
double amu2LFSfapprox_non_tan_beta_resummed(const MSSMNoFV_onshell&);

/// 2-loop photonic chargino contribution
double amu2LChipmPhotonic(const MSSMNoFV_onshell&);

/// 2-loop photonic neutralino contribution
double amu2LChi0Photonic(const MSSMNoFV_onshell&);

/// 2-loop 2L(a) sfermion contribution
double amu2LaSferm(const MSSMNoFV_onshell&);

/// 2-loop 2L(a) chargino/neutralino contribution
double amu2LaCha(const MSSMNoFV_onshell&);

} // namespace gm2calc

#endif
