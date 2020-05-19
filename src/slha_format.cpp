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

#include "slha_format.hpp"

namespace flexiblesusy {

/// SLHA line formatter for the MASS block entries
const char * const mass_formatter = " %9d   %16.8E   # %s\n";
/// SLHA line formatter for the mixing matrix entries  = NMIX, UMIX, VMIX, ...
const char * const mixing_matrix_formatter = " %2d %2d   %16.8E   # %s\n";
/// SLHA line formatter for vector entries
const char * const vector_formatter = " %5d   %16.8E   # %s\n";
/// SLHA number formatter
const char * const number_formatter = "         %16.8E   # %s\n";
/// SLHA line formatter for entries with three indices
const char * const tensor_formatter = " %8d %8d %8d   %16.8E   # %s\n";
/// SLHA scale formatter
const char * const scale_formatter = "%9.8E";
/// SLHA line formatter for the one-element entries  = HMIX, GAUGE, MSOFT, ...
const char * const single_element_formatter = " %5d   %16.8E   # %s\n";
/// SLHA line formatter for the SPINFO block entries
const char * const spinfo_formatter = " %5d   %s\n";

} // namespace flexiblesusy
