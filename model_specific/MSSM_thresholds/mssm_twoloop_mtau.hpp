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

// This file has been generated at Sat 22 Jul 2017 11:36:23
// with the script "tau_to_cpp.m".

#ifndef MSSM_TWO_LOOP_MTAU_H
#define MSSM_TWO_LOOP_MTAU_H

#include <iosfwd>

namespace flexiblesusy {
namespace mssm_twoloop_mtau {

using Real = long double;

struct Parameters {
    Real yt{};     ///< MSSM top Yukawa coupling DR-bar
    Real yb{};     ///< MSSM bottom Yukawa coupling DR-bar
    Real ytau{};   ///< MSSM tau Yukawa coupling DR-bar
    Real mt{};     ///< MSSM top mass DR-bar
    Real mb{};     ///< MSSM bottom mass DR-bar
    Real mtau{};   ///< MSSM tau mass DR-bar
    Real mst1{};   ///< MSSM light stop mass DR-bar
    Real mst2{};   ///< MSSM heavy stop mass DR-bar
    Real msb1{};   ///< MSSM light sbottom mass DR-bar
    Real msb2{};   ///< MSSM heavy sbottom mass DR-bar
    Real mstau1{}; ///< MSSM light stau mass DR-bar
    Real mstau2{}; ///< MSSM heavy stau mass DR-bar
    Real msntau{}; ///< MSSM tau sneutrino mass DR-bar
    Real xt{};     ///< MSSM stop mixing parameter DR-bar
    Real xb{};     ///< MSSM sbottom mixing parameter DR-bar
    Real xtau{};   ///< MSSM stau mixing parameter DR-bar
    Real mw{};     ///< MSSM W boson mass DR-bar
    Real mz{};     ///< MSSM Z boson mass DR-bar
    Real mh{};     ///< MSSM light CP-even Higgs mass DR-bar
    Real mH{};     ///< MSSM heavy CP-even Higgs mass DR-bar
    Real mC{};     ///< MSSM charged Higgs mass DR-bar
    Real mA{};     ///< MSSM CP-odd Higgs mass DR-bar
    Real mu{};     ///< MSSM mu superpotential parameter DR-bar
    Real tb{};     ///< MSSM tan(beta) DR-bar
    Real Q{};      ///< renormalization scale
};

/// 2-loop contribution to mtau O(alpha_tau^2) [0912.4652]
double delta_mtau_2loop_atau_atau(const Parameters&);

/// 2-loop contribution to mtau O(alpha_tau*alpha_t) [0912.4652]
double delta_mtau_2loop_atau_at(const Parameters&);

/// 2-loop contribution to mtau O(alpha_tau*alpha_b) [0912.4652]
double delta_mtau_2loop_atau_ab(const Parameters&);

std::ostream& operator<<(std::ostream&, const Parameters&);

} // namespace mssm_twoloop_mtau
} // namespace flexiblesusy

#endif
