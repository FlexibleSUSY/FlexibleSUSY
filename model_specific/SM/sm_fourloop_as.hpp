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

#ifndef SM_TWO_LOOP_AS_H
#define SM_TWO_LOOP_AS_H

#include <iosfwd>

namespace flexiblesusy {
namespace sm_fourloop_as {

struct Parameters {
    double as{};    ///< SM(5) strong coupling MS-bar
    double mt{};    ///< SM top mass MS-bar
    double Q{};     ///< renormalization scale
};

/// 1-loop O(alpha_s) contributions to Delta alpha_s [hep-ph/0004189]
double delta_alpha_s_1loop_as(const Parameters&);

/// 2-loop O(alpha_s^2) contributions to Delta alpha_s [hep-ph/0004189]
double delta_alpha_s_2loop_as_as(const Parameters&);

/// 3-loop O(alpha_s^3) contributions to Delta alpha_s [hep-ph/0004189]
double delta_alpha_s_3loop_as_as_as(const Parameters&);

/// 4-loop O(alpha_s^4) contributions to Delta alpha_s [hep-ph/0512060]
double delta_alpha_s_4loop_as_as_as_as(const Parameters&);

/// calculate alpha_s(SM(6)) from alpha_s(SM(5)) Eq (14) [hep-ph/0512060]
double calc_alpha_s(const Parameters&, int);

/// calculate alpha_s(SM(6)) from alpha_s(SM(5)) alternatively
double calc_alpha_s_alternative(const Parameters&, int);

std::ostream& operator<<(std::ostream&, const Parameters&);

} // namespace sm_fourloop_as
} // namespace flexiblesusy

#endif
