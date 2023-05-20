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

#ifndef SM_TWOLOOPHIGGS_H
#define SM_TWOLOOPHIGGS_H

namespace flexiblesusy {
namespace sm_twoloophiggs {

/// SM Higgs 1-loop contribution
double delta_mh_1loop_sm(
   double p, double scale, double mt, double yt,
   double v, double gY, double g2, double lambda);

/// SM Higgs 1-loop contribution, only O(alpha_t)
double delta_mh_1loop_at_sm(
   double p, double scale, double mt, double yt);

/// SM Higgs self-energy 2-loop, only O(alpha_t alpha_s)
double self_energy_higgs_2loop_at_as_sm(
   double p2, double scale, double mt, double yt, double g3);

/// SM Higgs self-energy 2-loop, only O(alpha_b alpha_s)
double self_energy_higgs_2loop_ab_as_sm(
   double p2, double scale, double mb, double yb, double g3);

/// SM Higgs self-energy 2-loop, only O((alpha_b + alpha_t)^2)
double self_energy_higgs_2loop_at_at_sm(
   double p2, double scale, double mt, double yt, double mb);

/// SM Higgs self-energy 2-loop, only O(alpha_tau^2)
double self_energy_higgs_2loop_atau_atau_sm(
   double p2, double scale, double mtau, double ytau);

/// SM Higgs tadpole 1-loop, only O(alpha_t)
double tadpole_higgs_1loop_at_sm(
   double scale, double mt, double yt);

/// SM Higgs tadpole 2-loop, only O(alpha_t alpha_s)
double tadpole_higgs_2loop_at_as_sm(
   double scale, double mt, double yt, double g3);

/// SM Higgs tadpole 2-loop, only O(alpha_b alpha_s)
double tadpole_higgs_2loop_ab_as_sm(
   double scale, double mb, double yb, double g3);

/// SM Higgs tadpole 2-loop, only O((alpha_b + alpha_t)^2)
double tadpole_higgs_2loop_at_at_sm(
   double scale, double mt, double yt, double mb);

/// SM Higgs tadpole 2-loop, only O(alpha_tau^2)
double tadpole_higgs_2loop_atau_atau_sm(
   double scale, double mtau, double ytau);

/// SM Higgs 2-loop contribution, only O(alpha_t alpha_s)
double delta_mh_2loop_at_as_sm(
   double p2, double scale, double mt, double yt, double g3);

/// SM Higgs 2-loop contribution, only O(alpha_t alpha_s)
double delta_mh_2loop_ab_as_sm(
   double p2, double scale, double mb, double yb, double g3);

/// SM Higgs 2-loop contribution, only O((alpha_b + alpha_t)^2)
double delta_mh_2loop_at_at_sm(
   double p2, double scale, double mt, double yt, double mb);

/// SM Higgs 2-loop contribution, only O(alpha_tau^2)
double delta_mh_2loop_atau_atau_sm(
   double p2, double scale, double mtau, double ytau);

/// SM Higgs 1-loop contribution from SUSYHD 1.0.2
double delta_mh_1loop_sm_SUSYHD(
   double vev, double Mt, double mh, double MW, double MZ, double Q);

/// SM Higgs 2-loop contribution from SUSYHD 1.0.2
double delta_mh_2loop_sm_SUSYHD(
   double vev, double Mt, double Mh, double g3);

} // namespace sm_twoloophiggs
} // namespace flexiblesusy

#endif
