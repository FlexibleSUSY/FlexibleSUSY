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

#ifndef MSSM_TWOLOOPHIGGS_H
#define MSSM_TWOLOOPHIGGS_H

#include <Eigen/Core>

/**
 * @file mssm_twoloophiggs.hpp
 * @brief function declarations for 2-loop MSSM Higgs self-energies
 *        and tadpoles
 *
 * Notation:
 *
 * mt2    : squared DR-bar top mass in the MSSM
 * mb2    : squared DR-bar bottom mass in the MSSM
 * mtau2  : squared DR-bar tau mass in the MSSM
 * mg     : DR-bar gluino mass in the MSSM
 * mA2    : squared DR-bar CP-odd Higgs mass in the MSSM
 * mst12  : squared DR-bar lightest stop mass
 * mst22  : squared DR-bar heaviest stop mass
 * msb12  : squared DR-bar lightest sbottom mass
 * msb22  : squared DR-bar heaviest sbottom mass
 * mstau12: squared DR-bar lightest stau mass
 * mstau22: squared DR-bar heaviest stau mass
 * msv2   : squared DR-bar tau sneutrino mass
 *
 * sxt    : sine of DR-bar stop mixing angle in the MSSM
 * cxt    : cosine of DR-bar stop mixing angle in the MSSM
 * sxb    : sine of DR-bar sbottom mixing angle in the MSSM
 * cxb    : cosine of DR-bar sbottom mixing angle in the MSSM
 * sintau : sine of DR-bar stau mixing angle in the MSSM
 * costau : cosine of DR-bar stau mixing angle in the MSSM
 *
 * gs     : DR-bar strong gauge coupling g3 in the MSSM
 * mu     : DR-bar mu-parameter in the MSSM (arXiv:0907.4682)
 * tanb   : DR-bar tan(beta) = vu/vd in the MSSM
 * cotb   : DR-bar 1/tan(beta) in the MSSM
 * vev2   : squared DR-bar vev^2 = (vu^2 + vd^2) in the MSSM
 *
 * scheme : DR-bar scheme (0) or on-shell scheme (1)
 */

namespace flexiblesusy {
namespace mssm_twoloophiggs {

// tadpoles

Eigen::Matrix<double, 2, 1> tadpole_higgs_2loop_at_as_mssm(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2,
   double mu, double tanb, double vev2, double gs);

Eigen::Matrix<double, 2, 1> tadpole_higgs_2loop_at_at_mssm(
   double mt2, double mb2, double mA2, double mst12,
   double mst22, double msb12, double msb22,
   double sxt, double cxt, double sxb, double cxb,
   double scale2, double mu, double tanb, double vev2);

Eigen::Matrix<double, 2, 1> tadpole_higgs_2loop_ab_as_mssm(
   double mb2, double mg, double msb12, double msb22,
   double sxb, double cxb, double scale2,
   double mu, double cotb, double vev2, double gs);

Eigen::Matrix<double, 2, 1> tadpole_higgs_2loop_atau_atau_mssm(
   double mtau2, double mA2, double msv2, double mstau12,
   double mstau22, double sintau, double costau, double scale2,
   double mu, double tanb, double vev2);

// self-energies

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_at_as_mssm(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double mu,
   double tanb, double vev2, double gs, int scheme = 0);

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_at_at_mssm(
   double mt2, double mb2, double mA2, double mst12,
   double mst22, double msb12, double msb22,
   double sxt, double cxt, double sxb, double cxb,
   double scale2, double mu, double tanb, double vev2);

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_ab_as_mssm(
   double mb2, double mg, double msb12, double msb22,
   double sxb, double cxb, double scale2, double mu,
   double cotb, double vev2, double gs, int scheme = 0);

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_atau_atau_mssm(
   double mtau2, double mA2, double msv2, double mstau12,
   double mstau22, double sintau, double costau, double scale2,
   double mu, double tanb, double vev2, int scheme = 0);


Eigen::Matrix<double, 2, 2> self_energy_pseudoscalar_2loop_at_as_mssm(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double mu,
   double tanb, double vev2, double gs);

Eigen::Matrix<double, 2, 2> self_energy_pseudoscalar_2loop_at_at_mssm(
   double mt2, double mb2, double mA2, double mst12,
   double mst22, double msb12, double msb22,
   double sxt, double cxt, double sxb, double cxb,
   double scale2, double mu, double tanb, double vev2);

Eigen::Matrix<double, 2, 2> self_energy_pseudoscalar_2loop_ab_as_mssm(
   double mb2, double mg, double msb12, double msb22,
   double sxb, double cxb, double scale2, double mu,
   double cotb, double vev2, double gs);

Eigen::Matrix<double, 2, 2> self_energy_pseudoscalar_2loop_atau_atau_mssm(
   double mtau2, double mA2, double msv2, double mstau12,
   double mstau22, double sintau, double costau, double scale2,
   double mu, double tanb, double vev2);

// self-energies with tadpoles added

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_at_as_mssm_with_tadpoles(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double mu,
   double tanb, double vev2, double gs, int scheme = 0);

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_at_at_mssm_with_tadpoles(
   double mt2, double mb2, double mA2, double mst12,
   double mst22, double msb12, double msb22,
   double sxt, double cxt, double sxb, double cxb,
   double scale2, double mu, double tanb, double vev2);

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_ab_as_mssm_with_tadpoles(
   double mb2, double mg, double msb12, double msb22,
   double sxb, double cxb, double scale2, double mu,
   double cotb, double vev2, double gs, int scheme = 0);

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_atau_atau_mssm_with_tadpoles(
   double mtau2, double mA2, double msv2, double mstau12,
   double mstau22, double sintau, double costau, double scale2,
   double mu, double tanb, double vev2, int scheme = 0);


double self_energy_pseudoscalar_2loop_at_as_mssm_with_tadpoles(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double mu,
   double tanb, double vev2, double gs);

double self_energy_pseudoscalar_2loop_at_at_mssm_with_tadpoles(
   double mt2, double mb2, double mA2, double mst12,
   double mst22, double msb12, double msb22,
   double sxt, double cxt, double sxb, double cxb,
   double scale2, double mu, double tanb, double vev2);

double self_energy_pseudoscalar_2loop_ab_as_mssm_with_tadpoles(
   double mb2, double mg, double msb12, double msb22,
   double sxb, double cxb, double scale2, double mu,
   double cotb, double vev2, double gs);

double self_energy_pseudoscalar_2loop_atau_atau_mssm_with_tadpoles(
   double mtau2, double mA2, double msv2, double mstau12,
   double mstau22, double sintau, double costau, double scale2,
   double mu, double tanb, double vev2);

} // namespace mssm_twoloophiggs
} // namespace flexiblesusy

#endif
