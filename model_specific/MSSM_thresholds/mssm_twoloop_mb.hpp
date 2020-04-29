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

// This file has been generated at Sat 18 Apr 2020 01:20:19
// with the script "bquark_to_cpp.m".

#ifndef MSSM_TWO_LOOP_SQCD_MB_H
#define MSSM_TWO_LOOP_SQCD_MB_H

namespace flexiblesusy {
namespace mssm_twoloop_mb {

struct Parameters {
    Parameters() = default;
    Parameters(double g3_, double mt_, double mb_, double mg_,
               double mst1_, double mst2_,
               double msb1_, double msb2_,
               double msusy_,
               double xt_, double xb_, double Q_)
       : g3(g3_), mt(mt_), mb(mb_), mg(mg_)
       , mst1(mst1_), mst2(mst2_)
       , msb1(msb1_), msb2(msb2_)
       , msusy(msusy_)
       , xt(xt_), xb(xb_), Q(Q_)
       {}

    double g3{};    ///< MSSM strong gauge coupling DR-bar
    double mt{};    ///< MSSM top mass DR-bar
    double mb{};    ///< SM   bottom mass DR-bar
    double mg{};    ///< MSSM gluino mass DR-bar
    double mst1{};  ///< MSSM light stop mass DR-bar
    double mst2{};  ///< MSSM heavy stop mass DR-bar
    double msb1{};  ///< MSSM light sbottom mass DR-bar
    double msb2{};  ///< MSSM heavy sbottom mass DR-bar
    double msusy{}; ///< MSSM light squark masses DR-bar
    double xt{};    ///< MSSM sbottom mixing parameter DR-bar
    double xb{};    ///< MSSM stop mixing parameter DR-bar
    double Q{};     ///< renormalization scale
};

/// 2-loop full SQCD contributions to mb [arXiv:0707.0650]
double delta_mb_2loop(const Parameters&);

} // namespace mssm_twoloop_mb
} // namespace flexiblesusy

#endif
