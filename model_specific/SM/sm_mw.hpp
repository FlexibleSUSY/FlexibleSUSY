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

#ifndef SM_MW_H
#define SM_MW_H

#include <utility>

namespace flexiblesusy {
namespace sm_mw {

/// returns W pole mass in the SM and corresponding uncertainty, using fit formula (MS-bar scheme)
std::pair<double, double> calculate_mw_pole_SM_fit_MSbar(double mh, double mt, double as, double da5had) noexcept;

/// returns W pole mass in the SM and corresponding uncertainty, using fit formula (OS scheme)
std::pair<double, double> calculate_mw_pole_SM_fit_OS(double mz, double mh, double mt, double as, double Da) noexcept;

} // namespace sm_mw
} // namespace flexiblesusy

#endif
