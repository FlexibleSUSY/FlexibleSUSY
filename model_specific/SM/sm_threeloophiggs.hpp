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

#ifndef SM_THREELOOPHIGGS_H
#define SM_THREELOOPHIGGS_H

namespace flexiblesusy {
namespace sm_threeloophiggs {

/// SM Higgs self-energy 3-loop, only O(alpha_t alpha_s^2)
double delta_mh_3loop_at_as_as_sm(
   double scale, double mt, double yt, double g3);

/// SM Higgs self-energy 3-loop, only O(alpha_t^2 alpha_s)
double delta_mh_3loop_at_at_as_sm(
   double scale, double mt, double yt, double g3);

/// SM Higgs self-energy 3-loop, only O(alpha_t^3)
double delta_mh_3loop_at_at_at_sm(
   double scale, double mt, double yt, double mh);

} // namespace sm_threeloophiggs
} // namespace flexiblesusy

#endif
