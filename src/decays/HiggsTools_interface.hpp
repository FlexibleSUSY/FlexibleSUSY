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
#include <fstream>
#include <tuple>

#include "decays/standard_model_decays.hpp"
#include "spectrum_generator_settings.hpp"
#include "minimizer.hpp"

namespace flexiblesusy {

std::complex<double> get_width_from_table(std::vector<std::tuple<std::string, int, int, double, std::complex<double>>> const& decay_table, std::array<int, 2> const& outPDGs) {
      for (const auto& decay: decay_table) {
         std::vector<int> final_state = {std::get<1>(decay), std::get<2>(decay)};
         if (std::is_permutation(outPDGs.begin(), outPDGs.end(), final_state.begin(), final_state.end())) {
            return std::get<4>(decay);
         }
      }
   return 0.;
}

void call_HiggsTools(
   std::vector<std::tuple<std::string, int, int, double, std::complex<double>>> const& bsm_input,
   Physical_input const& physical_input,
   softsusy::QedQcd const& qedqcd,
   Spectrum_generator_settings const& spectrum_generator_settings,
   FlexibleDecay_settings const& flexibledecay_settings) {

   std::set<std::string> Higgses;
   for (const auto& el : bsm_input) {
      Higgses.insert(std::get<0>(el));
   }

   auto pred = Higgs::Predictions();
   namespace HP = Higgs::predictions;
   // whether to calculate the ggH cross-section in terms of the effective top and bottom Yukawa couplings
   // or by rescaling the SM-like ggH XS by the squared of the effective gg coupling (no effects from colored BSM particles are taken into account)
   static constexpr bool calcggH = false;
   // whether to calculate the H->gaga decay width in terms of the effective couplings
   // or by rescaling the SM-like H->gaga decay by the squared of the effective gamgam coupling (no effects from charged BSM particles are taken into account).
   static constexpr bool calcHgamgam = false;

   for (auto& el : Higgses) {

      std::vector<std::tuple<std::string, int, int, double, std::complex<double>>> data;
      std::copy_if(
         bsm_input.begin(), bsm_input.end(), std::back_inserter(data),
         [&el](auto x) {return std::get<0>(x) == el;}
      );

      auto& s = pred.addParticle(HP::BsmParticle(el, HP::ECharge::neutral));
      const double mass = std::get<3>(data.at(0));
      s.setMass(mass);

      // create a SM equivalent to the BSM model, with mhSM == mass
      standard_model::Standard_model sm {};
      sm.initialise_from_input(qedqcd);
      sm.set_physical_input(physical_input);
      sm.set_pole_mass_loop_order(static_cast<int>(spectrum_generator_settings.get(Spectrum_generator_settings::pole_mass_loop_order)));
      sm.set_ewsb_loop_order(static_cast<int>(spectrum_generator_settings.get(Spectrum_generator_settings::ewsb_loop_order)));
      sm.set_precision(spectrum_generator_settings.get(Spectrum_generator_settings::precision));
      sm.set_threshold_corrections(spectrum_generator_settings.get_threshold_corrections());
      sm.set_loop_corrections(spectrum_generator_settings.get_loop_corrections());
      sm.set_loops(static_cast<int>(spectrum_generator_settings.get(Spectrum_generator_settings::beta_loop_order)));
      auto match_Higgs_mass = [&sm, mass](Eigen::Matrix<double,1,1> const& x) {
         sm.set_Lambdax(x[0]);
         sm.calculate_DRbar_masses(); // internaly calls solve_ewsb_tree_level()
         sm.solve_ewsb();
         sm.calculate_Mhh_pole();
         return std::abs(sm.get_physical().Mhh - mass);
      };
      Minimizer<1> minimizer(match_Higgs_mass, 100, 1.0e-5);
      Eigen::Matrix<double,1,1> start;
      start << 10;
      const int status = minimizer.minimize(start);
      const auto minimum_point = minimizer.get_solution();
      // finally, after fixing lambda to a value giving mhSM == mass,
      // compute all the other masses
      sm.calculate_pole_masses();

      std::cout << el << " matched masses: " << sm.get_physical().Mhh << ' ' << mass << ' ' << minimum_point[0] << ' '<< sm.get_physical().MVWp << std::endl;

      // calculate decays in the SM equivalent
      flexiblesusy::Standard_model_decays sm_decays(sm, qedqcd, physical_input, flexibledecay_settings);
      sm_decays.calculate_decays();
      const auto sm_input = sm_decays.get_higgstools_input();

      auto f = [&sm_input, &bsm_input](std::array<int, 2> const& a) {
         std::cout << "coup  " << a.at(0) << a.at(1) << ' ' << (std::abs(get_width_from_table(sm_input, a)) > 0 ? get_width_from_table(bsm_input, a)/get_width_from_table(sm_input, a) : 0.) << std::endl;
         return std::abs(get_width_from_table(sm_input, a)) > 0 ? get_width_from_table(bsm_input, a)/get_width_from_table(sm_input, a).real() : 0.;
      };

      auto effc = HP::NeutralEffectiveCouplings {};
      // quarks
      effc.dd = f({-1, 1});
      effc.uu = f({-2, 2});
      effc.ss = f({-3, 3});
      effc.cc = f({-4, 4});
      effc.bb = f({-5, 5});
      effc.tt = f({-6, 6});
      // leptons
      effc.ee =     f({-11, 11});
      effc.mumu =   f({-13, 13});
      effc.tautau = f({-15, 15});
      // gauge bosons
      effc.WW =     f({-24, 24}).real();
      effc.ZZ =     f({ 23, 23}).real();
      effc.Zgam =   f({ 23, 22}).real();
      effc.gamgam = f({ 22, 22}).real();
      effc.gg =     f({ 21, 21}).real();

      effectiveCouplingInput(s, effc, HP::ReferenceModel::SMHiggsEW, calcggH, calcHgamgam);
      // @todo: set total width to the one computed by FD as HiggsTools doesn't calculate
      // BSM decays of Higgs, e.g. H -> Ah Z
      // s.setTotalWidth(2.);
   }

   // HiggsBounds
   auto bounds = Higgs::Bounds {"/fs_dependencies/hbdataset"};
   auto hbResult = bounds(pred);
   std::ofstream hb_output("HiggsBounds.out");
   hb_output << hbResult;
   hb_output.close();
   std::cout << "All applied limits: obsRatio (expRatio)\n";
   for (const auto &al : hbResult.appliedLimits) {
      std::cout << al.limit()->id() << " " << al.limit()->processDesc()
                << ": " << al.obsRatio() << " (" << al.expRatio() << ")"
                << std::endl;
   }

   // HiggsSignals
   const auto signals = Higgs::Signals {"/fs_dependencies/hsdataset"};
   const double hs_chisq = signals(pred);
   std::cout << "\nHiggsSignals chisq: " << hs_chisq << " from "
              << signals.observableCount() << " observables" << std::endl;
}

} // flexiblesusy

#endif
