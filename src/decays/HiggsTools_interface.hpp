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

#include "model_specific/SM/decays/standard_model_decays.hpp"
#include "minimizer.hpp"

namespace flexiblesusy {

std::complex<double> get_width_from_table(std::vector<std::tuple<std::string, int, int, double, std::complex<double>>> const& decay_table, std::array<int, 2> const& outPDGs) {
      for (const auto& decay: decay_table) {
         std::vector<int> final_state = {std::get<1>(decay), std::get<2>(decay)};
         if (std::is_permutation(outPDGs.begin(), outPDGs.end(), final_state.begin(), final_state.end())) {
            return std::get<3>(decay);
         }
      }
   return 0.;
}

void call_HiggsTools(
   std::vector<std::tuple<std::string, int, int, double, std::complex<double>>> const& bsm_input,
   Physical_input const& physical_input,
   softsusy::QedQcd const& qedqcd,
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

      auto predicate = [&el](auto x) {return std::get<0>(x) == el;};
      std::vector<std::tuple<std::string, int, int, double, std::complex<double>>> data;
      std::copy_if(bsm_input.begin(), bsm_input.end(), std::back_inserter(data), predicate);

      auto& s = pred.addParticle(HP::BsmParticle(el, HP::ECharge::neutral));
      const double mass = std::get<3>(data.at(0));
      s.setMass(mass);

      standard_model::Standard_model sm {};
      sm.initialise_from_input(qedqcd);
      auto match_Higgs_mass = [&sm, mass](Eigen::Matrix<double,1,1> const& x) {
         sm.set_Lambdax(x[0]);
         sm.solve_ewsb_tree_level();
         sm.calculate_DRbar_masses();
         sm.calculate_pole_masses();
         return std::abs(sm.get_physical().Mhh - mass);
      };
      Minimizer<1> minimizer(match_Higgs_mass, 100, 1.0e-5);
      Eigen::Matrix<double,1,1> start;
      start << 10;
      const int status = minimizer.minimize(start);
      const auto minimum_point = minimizer.get_solution();
      std::cout << "matched masses: " << sm.get_physical().Mhh << ' ' << mass << std::endl;
    
   //BOOST_CHECK_EQUAL(status, GSL_SUCCESS);    
   //BOOST_CHECK_SMALL(minimizer.get_minimum_value(), 1.0e-5);    
   //BOOST_CHECK_CLOSE_FRACTION(minimum_point(0), 5.0, 1.0e-5);    
   //BOOST_CHECK_CLOSE_FRACTION(minimum_point(1), 1.0, 1.0e-5);

      flexiblesusy::Standard_model_decays sm_decays(sm, qedqcd, physical_input, flexibledecay_settings);
      sm_decays.calculate_decays();

      const auto sm_input = sm_decays.get_higgstools_input();

      auto effc = HP::NeutralEffectiveCouplings {};
      auto f = [&sm_input, &bsm_input](std::array<int, 2> const& a) {
         return std::abs(get_width_from_table(sm_input, a)) > 0 ? get_width_from_table(bsm_input, a)/get_width_from_table(sm_input, a) : 0.;
      };
      // quarks
      effc.uu = f({-2, 2});
      effc.dd = f({-1, 1});
      effc.cc = f({-4, 4});
      effc.ss = f({-3, 3});
      effc.tt = f({-6, 6});
      effc.bb = f({-5, 5});
      // leptons
      effc.ee = f({-11, 11});
      effc.mumu = f({-13, 13});
      effc.tautau = f({-15, 15});
      // gauge bosons
      effc.WW = f({-24, 24}).real();
      effc.ZZ = f({ 23, 23}).real();
      effc.Zgam = f({ 23, 22}).real();
      effc.gamgam = f({ 22, 22}).real();
      effc.gg = f({ 21, 21}).real();

      effectiveCouplingInput(s, effc, HP::ReferenceModel::SMHiggs, calcggH, calcHgamgam);
   }

   auto bounds = Higgs::Bounds {"/run/media/scratch/Pobrane/hbdataset-master"};
   auto hbResult = bounds(pred);
   std::ofstream hb_output;
   hb_output.open("HiggsBounds.out");
   hb_output << hbResult;
   hb_output.close();
   std::cout << "All applied limits: obsRatio (expRatio)\n";
   for (const auto &al : hbResult.appliedLimits) {
      std::cout << al.limit()->id() << " " << al.limit()->processDesc()
                << ": " << al.obsRatio() << " (" << al.expRatio() << ")"
                << std::endl;
   }

   const auto signals = Higgs::Signals {"/run/media/scratch/Pobrane/hsdataset-main"};
   const double hs_chisq = signals(pred);
   std::cout << "\nHiggsSignals chisq: " << hs_chisq << " from "
              << signals.observableCount() << " observables" << std::endl;
}

} // flexiblesusy

#endif
