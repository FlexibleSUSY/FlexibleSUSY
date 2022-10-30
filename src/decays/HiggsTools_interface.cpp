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
 * @file HiggsTools_interface.cpp
 *
 * @brief contains interface to HiggsTools
 */

#include "Higgs/Predictions.hpp"
#include "Higgs/Bounds.hpp"
#include "Higgs/Signals.hpp"

#include "HiggsTools_interface.hpp"
#include "decays/standard_model_decays.hpp"
#include "minimizer.hpp"

#include <algorithm>
#include <iostream>
#include <fstream>

namespace flexiblesusy {

std::pair<int, double> call_HiggsTools(
   EffectiveCoupling_list const& bsm_input,
   std::vector<SingleChargedHiggsInput> const& bsm_input2,
   Physical_input const& physical_input,
   softsusy::QedQcd const& qedqcd,
   Spectrum_generator_settings const& spectrum_generator_settings,
   FlexibleDecay_settings const& flexibledecay_settings) {

   auto pred = Higgs::Predictions();
   namespace HP = Higgs::predictions;
   // whether to calculate the ggH cross-section in terms of the effective top and bottom Yukawa couplings
   // or by rescaling the SM-like ggH XS by the squared of the effective gg coupling (no effects from colored BSM particles are taken into account)
   static constexpr bool calcggH = false;
   // whether to calculate the H->gaga decay width in terms of the effective couplings
   // or by rescaling the SM-like H->gaga decay by the squared of the effective gamgam coupling (no effects from charged BSM particles are taken into account).
   static constexpr bool calcHgamgam = false;

   for (auto const& el : bsm_input) {

      auto& s = pred.addParticle(HP::BsmParticle(el.particle, HP::ECharge::neutral));
      const double mass = el.mass;
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
      // set SM λ such that mhSM == mass
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

      if (sm.get_physical().Mhh > 0) {
         // calculate decays in the SM equivalent
         flexiblesusy::Standard_model_decays sm_decays(sm, qedqcd, physical_input, flexibledecay_settings);
         sm_decays.calculate_decays();
         const auto sm_input = sm_decays.get_higgstools_input();

         auto effc = HP::NeutralEffectiveCouplings {};
         // quarks
         effc.dd = std::abs(sm_input[0].dd) > 0 ? el.dd/sm_input[0].dd.real() : 0.;
         effc.uu = std::abs(sm_input[0].uu) > 0 ? el.uu/sm_input[0].uu.real() : 0.;
         effc.ss = std::abs(sm_input[0].ss) > 0 ? el.ss/sm_input[0].ss.real() : 0.;
         effc.cc = std::abs(sm_input[0].cc) > 0 ? el.cc/sm_input[0].cc.real() : 0.;
         effc.bb = std::abs(sm_input[0].bb) > 0 ? el.bb/sm_input[0].bb.real() : 0.;
         effc.tt = std::abs(sm_input[0].tt) > 0 ? el.tt/sm_input[0].tt.real() : 0.;
         // leptons
         effc.ee = std::abs(sm_input[0].ee) > 0 ? el.ee/sm_input[0].ee.real() : 0.;
         effc.mumu = std::abs(sm_input[0].mumu) > 0 ? el.mumu/sm_input[0].mumu.real() : 0.;
         effc.tautau = std::abs(sm_input[0].tautau) > 0 ? el.tautau/sm_input[0].tautau.real() : 0.;
         // gauge bosons
         effc.WW = std::abs(sm_input[0].WW) > 0 ? el.WW/sm_input[0].WW : 0.;
         effc.ZZ = std::abs(sm_input[0].ZZ) > 0 ? el.ZZ/sm_input[0].ZZ : 0.;
         effc.gamgam = std::abs(sm_input[0].gamgam) > 0 ? el.gamgam/sm_input[0].gamgam : 0.;
         effc.Zgam = std::abs(sm_input[0].Zgam) > 0 ? el.Zgam/sm_input[0].Zgam : 0.;
         effc.gg = std::abs(sm_input[0].gg) > 0 ? el.gg/sm_input[0].gg : 0.;

         effectiveCouplingInput(s, effc, HP::ReferenceModel::SMHiggsEW, calcggH, calcHgamgam);
         // effective coupligs are defined as sqrt(Gamma CP-even) + I sqrt(Gamma CP-odd)
         // so taking a norm gives a total partial width
         s.setDecayWidth(HP::Decay::emu,   std::norm(el.emu));
         s.setDecayWidth(HP::Decay::etau,  std::norm(el.etau));
         s.setDecayWidth(HP::Decay::mutau, std::norm(el.mutau));
         // set total width to the one computed by FD as HiggsTools doesn't calculate
         // some decays of Higgs at all, e.g. H -> Ah Z
         if (el.width > s.totalWidth()) {
            s.setDecayWidth("NP", "NP", el.width - s.totalWidth());
         }
      }
   }

   for (auto const& el : bsm_input2) {
      auto& Hpm = pred.addParticle(HP::BsmParticle(el.particle, HP::ECharge::single));
      // set mass
      Hpm.setMass(el.mass);

      // set H± decays
      if (el.width > 0) {
         Hpm.setTotalWidth(el.width);
         Hpm.setBr(HP::Decay::enu, el.brenu);
         Hpm.setBr(HP::Decay::munu, el.brmunu);
         Hpm.setBr(HP::Decay::taunu, el.brtaunu);
         Hpm.setBr(HP::Decay::tb, el.brtb);
         Hpm.setBr(HP::Decay::cs, el.brcs);
         Hpm.setBr(HP::Decay::WZ, el.brWZ);
         Hpm.setBr(HP::Decay::Wgam, el.brWgam);
      }

      // production
      // pp->ttbar, t(tbar)->Hp b
      double ppHpmtb_xsec = HP::EffectiveCouplingCxns::ppHpmtb(HP::Collider::LHC13, el.mass, el.cHpmtbR, el.cHpmtbL, el.brtHpb);
      Hpm.setCxn(HP::Collider::LHC13, HP::Production::Hpmtb, ppHpmtb_xsec);
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
              << signals.observableCount() << " observables" << ' ' << "(χ^2/ndf=" << hs_chisq/signals.observableCount() << ")" << std::endl;

   return {signals.observableCount(), hs_chisq};
}

} // flexiblesusy
