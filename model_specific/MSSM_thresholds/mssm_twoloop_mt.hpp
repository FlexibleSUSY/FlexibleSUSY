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

// This file has been generated at Wed 17 Jun 2020 11:16:10
// with the script "tquark_to_cpp.m".

#ifndef MSSM_TWO_LOOP_SQCD_MT_H
#define MSSM_TWO_LOOP_SQCD_MT_H

namespace flexiblesusy {
namespace mssm_twoloop_mt {

struct Parameters {
    Parameters() = default;
    Parameters(double g3_, double mt_, double mg_, double mst1_,
               double mst2_, double msusy_, double xt_, double Q_)
       : g3(g3_), mt(mt_), mg(mg_), mst1(mst1_)
       , mst2(mst2_), msusy(msusy_), xt(xt_), Q(Q_)
       {}

    double g3{};    ///< MSSM strong gauge coupling DR-bar
    double mt{};    ///< MSSM top mass DR-bar
    double mg{};    ///< MSSM gluino mass DR-bar
    double mst1{};  ///< MSSM light stop mass DR-bar
    double mst2{};  ///< MSSM heavy stop mass DR-bar
    double msusy{}; ///< MSSM SUSY particle mass scale
    double xt{};    ///< MSSM stop mixing parameter DR-bar
    double Q{};     ///< renormalization scale
};

/// 1-loop QCD contributions to Delta Mt over mt [hep-ph/0210258]
double dMt_over_mt_1loop_qcd(const Parameters&);
/// 1-loop SUSY contributions to Delta Mt over mt [hep-ph/0210258]
double dMt_over_mt_1loop_susy(const Parameters&);
/// 1-loop full SQCD contributions to Delta Mt over mt [hep-ph/0210258]
double dMt_over_mt_1loop(const Parameters&);

/// 2-loop QCD contributions to Delta Mt over mt [hep-ph/0507139]
double dMt_over_mt_2loop_qcd(const Parameters&);
/// 2-loop SUSY contributions to Delta Mt over mt [hep-ph/0507139]
double dMt_over_mt_2loop_susy(const Parameters&);
/// 2-loop full SQCD contributions to Delta Mt over mt [hep-ph/0507139]
double dMt_over_mt_2loop(const Parameters&);

} // namespace mssm_twoloop_mt
} // namespace flexiblesusy

#endif
