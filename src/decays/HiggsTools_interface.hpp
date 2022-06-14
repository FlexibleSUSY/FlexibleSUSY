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


/**
 * @file HiggsTools_interface.hpp
 *
 * @brief contains interface to HiggsTools
 */

#ifndef HIGGSTOOLS_INTERFACE_H
#define HIGGSTOOLS_INTERFACE_H

#include "Higgs/Predictions.hpp"
#include "Higgs/Bounds.hpp"
#include "Higgs/Signals.hpp"

#include <array>
#include <algorithm>
#include <iostream>

#include "model_specific/SM/decays/standard_model_decays.hpp"

namespace flexiblesusy {

template <typename T>
double get_width_from_table(T const& decay_table, int inPDG, std::array<int, 2> const& outPDGs) {
   for (const auto& particle: decay_table) {
      const int pdg = particle.get_particle_id();
      if (pdg != inPDG) continue;
      for (const auto& decay: particle) {
         std::vector<int> final_state = decay.second.get_final_state_particle_ids();
         if (std::is_permutation(outPDGs.begin(), outPDGs.end(), final_state.begin(), final_state.end())) {
            return decay.second.get_width();
         }
      }
   }
   return 0.;
}

template <typename Decay_table>
void call_HiggsTools(
   Decay_table const& decay_table,
   Physical_input const& physical_input,
   softsusy::QedQcd const& qedqcd,
   FlexibleDecay_settings const& flexibledecay_settings) {

   standard_model::Standard_model sm {};
   sm.initialise_from_input(qedqcd);
   sm.solve_ewsb_tree_level();
   sm.calculate_DRbar_masses();
   sm.calculate_pole_masses();

   flexiblesusy::Standard_model_decays sm_decays(sm, qedqcd, physical_input, flexibledecay_settings);
   sm_decays.calculate_decays();

   const auto sm_decay_table = sm_decays.get_decay_table();
   for (const auto& fs: {std::array<int, 2>{21, 21}, {5,-5}, {22,22}, {15,-15}, {24, -24}, {23, 23}, {22, 22}, {22, 23}}) {
      std::cout << "{" << fs.at(0) << "," << fs.at(1) << "}: " << std::sqrt(get_width_from_table(decay_table, 25, fs)/get_width_from_table(sm_decay_table, 25, fs)) << std::endl;
   }

   auto pred = Higgs::Predictions{};
   // set model predictions on pred
   namespace HP = Higgs::predictions;
   auto& s = pred.addParticle(HP::BsmParticle("h1", HP::ECharge::neutral));
   s.setMass(125);


   // whether to calculate the ggH cross-section in terms of the effective top and bottom Yukawa couplings
   // or by rescaling the SM-like ggH XS by the squared of the effective gg coupling (no effects from colored BSM particles are taken into account)
   static constexpr bool calcggH = false;
   // whether to calculate the H->gaga decay width in terms of the effective couplings
   // or by rescaling the SM-like H->gaga decay by the squared of the effective gamgam coupling (no effects from charged BSM particles are taken into account).
   static constexpr bool calcHgamgam = false;

   auto effc = HP::NeutralEffectiveCouplings{
            2., 10., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2.};
   effectiveCouplingInput(s, effc, HP::ReferenceModel::SMHiggs, calcggH, calcHgamgam);

   auto& s2 = pred.addParticle(HP::BsmParticle("h2", HP::ECharge::neutral));
   s.setMass(255);
   auto effc2 = HP::NeutralEffectiveCouplings{
            20., 10., 20., 20., 20., 20., 20., 20., 20., 20., 20., 20., 20., 20.};
   effectiveCouplingInput(s2, effc2, HP::ReferenceModel::SMHiggs, calcggH, calcHgamgam);

   auto bounds = Higgs::Bounds {"/fs_dependencies/hbdataset"};
   auto hbResult = bounds(pred);
   // std::cout << hbResult << std::endl;
   // std::cout << "All applied limits: obsRatio (expRatio)\n";
   // for (const auto &al : hbResult.appliedLimits) {
   //   std::cout << al.limit()->id() << " " << al.limit()->processDesc()
   //             << ": " << al.obsRatio() << " (" << al.expRatio() << ")"
   //             << std::endl;
   //}

   const auto signals = Higgs::Signals {"/fs_dependencies/hsdataset"};
   auto resultHS = signals(pred);
   //std::cout << "\n HiggsSignals chisq: " << resultHS << " from "
   //           << signals.observableCount() << " observables" << std::endl;
}

} // flexiblesusy

#endif
