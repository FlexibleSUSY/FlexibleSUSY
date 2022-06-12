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

#include <array>
#include <algorithm>
#include <iostream>

#include "model_specific/SM/decays/standard_model_decays.hpp"
// #include "Higgs/Predictions.hpp"
// #include "Higgs/Bounds.hpp"
// #include "Higgs/Signals.hpp"

namespace flexiblesusy {

template <typename T>
double get_width_from_table(T const& decay_table, int inPDG, std::array<int, 2> const& outPDGs) {
   for (const auto& particle: decay_table) {
      const int pdg = particle.get_particle_id();
      if (pdg != inPDG) continue;
      for (const auto& decay: particle) {
         std::vector<int> final_state = decay.second.get_final_state_particle_ids();
         //if (decay.second.get_final_state_particle_ids()[0] != outPDGs.at(0) || decay.second.get_final_state_particle_ids()[1] != outPDGs.at(1)) continue;
         std::array<int, 2> sorted_outPDGs = outPDGs;
         std::sort(sorted_outPDGs.begin(), sorted_outPDGs.end(), std::greater<int>());
         std::sort(final_state.begin(), final_state.end(), std::greater<int>());
         if (std::equal(sorted_outPDGs.begin(), sorted_outPDGs.end(), final_state.begin())) {
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

   flexiblesusy::SM_decays sm_decays(sm, qedqcd, physical_input, flexibledecay_settings);
   sm_decays.calculate_decays();

   const auto sm_decay_table = sm_decays.get_decay_table();
   for (const auto& fs: {std::array<int, 2>{21, 21}, {5,-5}, {22,22}, {15,-15}}) {
      std::cout << "{" << fs.at(0) << "," << fs.at(1) << "}: " << std::sqrt(get_width_from_table(decay_table, 25, fs)/get_width_from_table(sm_decay_table, 25, fs)) << std::endl;
   }


   // auto bounds = Higgs::Bounds{"/path/to/HBDataSet"};
   // auto pred = Higgs::Predictions{};
   // set model predictions on pred

   // auto hbResult = bounds(pred);
   // std::cout << hbResult << std::endl;
}

} // flexiblesusy

#endif
