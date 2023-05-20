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

#ifndef NMSSM_TWOLOOPHIGGS_H
#define NMSSM_TWOLOOPHIGGS_H

#include <Eigen/Core>

/**
 * @file nmssm_twoloophiggs.hpp
 * @brief function declarations for 2-loop NMSSM Higgs self-energies
 *        and tadpoles
 *
 * Notation:
 *
 * mt2    : squared DR-bar top mass in the NMSSM
 * mb2    : squared DR-bar bottom mass in the NMSSM
 * mtau2  : squared DR-bar tau mass in the NMSSM
 * mg     : DR-bar gluino mass in the NMSSM
 * mA2    : squared DR-bar CP-odd Higgs mass in the NMSSM
 * mst12  : squared DR-bar lightest stop mass
 * mst22  : squared DR-bar heaviest stop mass
 * msb12  : squared DR-bar lightest sbottom mass
 * msb22  : squared DR-bar heaviest sbottom mass
 * mstau12: squared DR-bar lightest stau mass
 * mstau22: squared DR-bar heaviest stau mass
 *
 * sxt    : sine of DR-bar stop mixing angle in the NMSSM
 * cxt    : cosine of DR-bar stop mixing angle in the NMSSM
 * sxb    : sine of DR-bar sbottom mixing angle in the NMSSM
 * cxb    : cosine of DR-bar sbottom mixing angle in the NMSSM
 * sintau : sine of DR-bar stau mixing angle in the NMSSM
 * costau : cosine of DR-bar stau mixing angle in the NMSSM
 *
 * gs     : DR-bar strong gauge coupling g3 in the NMSSM
 * mu     : DR-bar mu-parameter in the NMSSM (arXiv:0907.4682)
 * tanb   : DR-bar tan(beta) = vu/vd in the NMSSM
 * cotb   : DR-bar 1/tan(beta) in the NMSSM
 * vev2   : squared DR-bar vev^2 = (vu^2 + vd^2) in the NMSSM
 * lam    : DR-bar lambda in the NNMSSM
 * svev   : DR-bar singlet VEV = mu_eff / lam
 *
 * scheme : DR-bar scheme (0) or on-shell scheme (1)
 */

namespace flexiblesusy {
namespace nmssm_twoloophiggs {

Eigen::Matrix<double, 3, 1> tadpole_higgs_2loop_at_as_nmssm(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2,
   double mu, double tanb, double vev2, double gs, double svev);

Eigen::Matrix<double, 3, 1> tadpole_higgs_2loop_ab_as_nmssm(
   double mb2, double mg, double msb12, double msb22,
   double sxb, double cxb, double scale2,
   double mu, double cotb, double vev2, double gs, double svev);

Eigen::Matrix<double, 3, 3> self_energy_higgs_2loop_at_as_nmssm(
   double rmt, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double tanb, double vev2,
   double lam, double svev, double as, double mu);

Eigen::Matrix<double, 3, 3> self_energy_higgs_2loop_ab_as_nmssm(
   double rmb, double mg, double msb12, double msb22,
   double sxb, double cxb, double scale2, double cotb, double vev2,
   double lam, double svev, double as, double mu);

Eigen::Matrix<double, 3, 3> self_energy_pseudoscalar_2loop_at_as_nmssm(
   double rmt, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double tanb, double vev2,
   double lam, double svev, double as, double mu);

Eigen::Matrix<double, 3, 3> self_energy_pseudoscalar_2loop_ab_as_nmssm(
   double rmb, double mg, double msb12, double msb22,
   double sxb, double cxb, double scale2, double cotb, double vev2,
   double lam, double svev, double as, double mu);

Eigen::Matrix<double, 3, 3> self_energy_higgs_2loop_at_as_nmssm_with_tadpoles(
   double rmt, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double tanb, double vev2,
   double lam, double svev, double as);

Eigen::Matrix<double, 3, 3> self_energy_higgs_2loop_ab_as_nmssm_with_tadpoles(
   double rmb, double mg, double msb12, double msb22,
   double sxb, double cxb, double scale2, double cotb, double vev2,
   double lam, double svev, double as);

Eigen::Matrix<double, 3, 3> self_energy_pseudoscalar_2loop_at_as_nmssm_with_tadpoles(
   double rmt, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double tanb, double vev2,
   double lam, double svev, double as);

Eigen::Matrix<double, 3, 3> self_energy_pseudoscalar_2loop_ab_as_nmssm_with_tadpoles(
   double rmb, double mg, double msb12, double msb22,
   double sxb, double cxb, double scale2, double cotb, double vev2,
   double lam, double svev, double as);

} // namespace nmssm_twoloophiggs
} // namespace flexiblesusy

#endif
