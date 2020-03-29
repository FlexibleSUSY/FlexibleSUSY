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

#ifndef GM2_2LOOP_HELPERS_HPP
#define GM2_2LOOP_HELPERS_HPP

#include <Eigen/Core>

namespace gm2calc {

class MSSMNoFV_onshell;

// === 2-loop fermion/sfermion approximations ===

// routines for individual fermion/sfermion 2-loop contributions for
// approximations

double amu2LWHnu(const MSSMNoFV_onshell&);
double amu2LWHmuL(const MSSMNoFV_onshell&);
double amu2LBHmuL(const MSSMNoFV_onshell&);
double amu2LBHmuR(const MSSMNoFV_onshell&);
double amu2LBmuLmuR(const MSSMNoFV_onshell&);

// routines for sub-expressions

double log_scale(const MSSMNoFV_onshell&);

double delta_g1(const MSSMNoFV_onshell&);
double delta_g2(const MSSMNoFV_onshell&);
double delta_yuk_higgsino(const MSSMNoFV_onshell&);
double delta_yuk_bino_higgsino(const MSSMNoFV_onshell&);
double delta_yuk_wino_higgsino(const MSSMNoFV_onshell&);
double delta_tan_beta(const MSSMNoFV_onshell&);

// === SUSY 2L(a) diagrams ===

// routines for sub-expressions (with tan(beta) resummation)

double tan_alpha(const MSSMNoFV_onshell&);
Eigen::Matrix<std::complex<double>,3,3> lambda_mu_cha(const MSSMNoFV_onshell&);
Eigen::Matrix<std::complex<double>,2,2> lambda_stop(const MSSMNoFV_onshell&);
Eigen::Matrix<std::complex<double>,2,2> lambda_sbot(const MSSMNoFV_onshell&);
Eigen::Matrix<std::complex<double>,2,2> lambda_stau(const MSSMNoFV_onshell&);

} // namespace gm2calc

#endif
