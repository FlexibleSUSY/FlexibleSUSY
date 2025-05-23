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

#ifndef AMM_LOOP_FUNCTIONS_H
#define AMM_LOOP_FUNCTIONS_H

namespace flexiblesusy {
namespace amm_loop_functions {
namespace two_loop {

/// Barr-Zee 2-loop function with fermion loop
double BarrZeeLoopFP(double);
/// Barr-Zee 2-loop function with fermion loop
double BarrZeeLoopFS(double);
/// Barr-Zee 2-loop function with fermion loop
double BarrZeeLoopFPZ(double, double);
/// Barr-Zee 2-loop function with fermion loop
double BarrZeeLoopFSZ(double, double);
/// Barr-Zee 2-loop function with scalar loop
double BarrZeeLoopS(double);
/// Barr-Zee 2-loop function with vector
double BarrZeeLoopV(double);

} // namespace two_loop
} // namespace amm_loop_functions
} // namespace flexiblesusy

#endif
