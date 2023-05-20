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

#ifndef THRESHOLD_LOOP_FUNCTIONS_H
#define THRESHOLD_LOOP_FUNCTIONS_H

#define TCF(n) threshold_loop_functions::F ## n
#define TCf(n) threshold_loop_functions::f ## n
#define TCfth(n) threshold_loop_functions::fth ## n
#define TCf0 threshold_loop_functions::f
#define TCg0 threshold_loop_functions::g
#define TCIabc threshold_loop_functions::Iabc
#define TCB0 threshold_loop_functions::B0
#define TCDB0 threshold_loop_functions::DB0
#define TCC0 threshold_loop_functions::C0
#define TCD0 threshold_loop_functions::D0
#define TCD2t threshold_loop_functions::D2t
#define TCD4t threshold_loop_functions::D4t
#define TCW threshold_loop_functions::W
#define TDelta threshold_loop_functions::delta_xyz
#define TPhi threshold_loop_functions::phi_xyz
#define TCD1F(n) threshold_loop_functions::D1F ## n
#define TCD2F(n) threshold_loop_functions::D2F ## n
#define TCD1f(n) threshold_loop_functions::D1f ## n
#define TCD10f(n) threshold_loop_functions::D10f ## n
#define TCD01f(n) threshold_loop_functions::D01f ## n
#define TCD1f0 threshold_loop_functions::D1f
#define TCD1g0 threshold_loop_functions::D1g

#include "cextensions.hpp"

namespace flexiblesusy {
namespace threshold_loop_functions {

// loop functions from arXiv:1407.4081

#define TCFATTR noexcept ATTR(const)

double F1(double) TCFATTR;
double F2(double) TCFATTR;
double F3(double) TCFATTR;
double F4(double) TCFATTR;
double F5(double) TCFATTR;
double F6(double) TCFATTR;
double F7(double) TCFATTR;
double F8(double, double) TCFATTR;
double F9(double, double) TCFATTR;

double f(double) TCFATTR;
double g(double) TCFATTR;

double f1(double) TCFATTR;
double f2(double) TCFATTR;
double f3(double) TCFATTR;
double f4(double) TCFATTR;
double f5(double, double) TCFATTR;
double f6(double, double) TCFATTR;
double f7(double, double) TCFATTR;
double f8(double, double) TCFATTR;

// 2-loop threshold function fth[y] from MhEFT-1.1
double fth1(double) TCFATTR;
double fth2(double) TCFATTR;
double fth3(double) TCFATTR;

// first derivatives

double D1F1(double) TCFATTR;
double D1F2(double) TCFATTR;
double D1F3(double) TCFATTR;
double D1F4(double) TCFATTR;
double D1F5(double) TCFATTR;
double D1F6(double) TCFATTR;
double D1F7(double) TCFATTR;
double D1f(double) TCFATTR;
double D1g(double) TCFATTR;
double D1f1(double) TCFATTR;
double D1f2(double) TCFATTR;
double D1f3(double) TCFATTR;
double D1f4(double) TCFATTR;
double D10f5(double, double) TCFATTR;
double D01f5(double, double) TCFATTR;
double D10f6(double, double) TCFATTR;
double D01f6(double, double) TCFATTR;
double D10f7(double, double) TCFATTR;
double D01f7(double, double) TCFATTR;
double D10f8(double, double) TCFATTR;
double D01f8(double, double) TCFATTR;

// second derivatives

double D2F1(double) TCFATTR;
double D2F2(double) TCFATTR;
double D2F3(double) TCFATTR;
double D2F4(double) TCFATTR;
double D2F5(double) TCFATTR;
double D2F6(double) TCFATTR;
double D2F7(double) TCFATTR;

/// \f$I_{abc}(a,b,c)\f$ (arguments are interpreted as unsquared)
double Iabc(double, double, double) TCFATTR;

/// \f$Delta_{xyz}(x,y,z)\f$ (arguments are interpreted as squared masses)
double delta_xyz(double, double, double) TCFATTR;

/// \f$phi_{xyz}(x,y,z)\f$ (arguments are interpreted as squared masses)
double phi_xyz(double, double, double) TCFATTR;

/// \f$B_0(p=0,m_1,m_2,Q)\f$ (arguments are interpreted as unsquared)
double B0(double, double, double) TCFATTR;

/// \f$B_0'(p=0,m_1,m_2)\f$ (arguments are interpreted as unsquared)
double DB0(double, double) TCFATTR;

/// \f$C_0(p=0,m_1,m_2,m_3)\f$ (arguments are interpreted as unsquared)
double C0(double, double, double) TCFATTR;

/// \f$D_0(p=0,m_1,m_2,m_3,m_4)\f$ (arguments are interpreted as unsquared)
double D0(double, double, double, double) TCFATTR;

/// \f$\tilde{D}_2(m_1,m_2,m_3,m_4)\f$ (arguments are interpreted as unsquared)
double D2t(double, double, double, double) TCFATTR;

/// \f$\tilde{D}_4(m_1,m_2,m_3,m_4,Q)\f$ (arguments are interpreted as unsquared)
double D4t(double, double, double, double, double) TCFATTR;

/// \f$Q(m_1,m_2,Q)\f$ (arguments are interpreted as unsquared)
double W(double, double, double) TCFATTR;

} // namespace threshold_loop_functions
} // namespace flexiblesusy

#endif
