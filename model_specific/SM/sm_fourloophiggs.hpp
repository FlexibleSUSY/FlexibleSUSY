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

#ifndef SM_FOURLOOPHIGGS_H
#define SM_FOURLOOPHIGGS_H

namespace flexiblesusy {
namespace sm_fourloophiggs {

/// SM Higgs mass 4-loop contribution, only O(alpha_t alpha_s^3)
double delta_mh_4loop_at_as_as_as_sm(
   double scale, double mt, double yt, double g3);

} // namespace sm_fourloophiggs
} // namespace flexiblesusy

#endif
