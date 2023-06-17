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

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <functional>
#include <random>

#include <gsl/gsl_min.h>
#include <gsl/gsl_errno.h>

namespace flexiblesusy {

std::tuple<int, double, std::vector<std::tuple<int, double, double, std::string>>> call_HiggsTools(
   EffectiveCoupling_list const& bsm_input,
   std::vector<SingleChargedHiggsInput> const& bsm_input2,
   Physical_input const& physical_input,
   softsusy::QedQcd const& qedqcd,
   Spectrum_generator_settings const& spectrum_generator_settings,
   FlexibleDecay_settings const& flexibledecay_settings,
   std::string const& higgsbounds_dataset, std::string const& higgssignals_dataset) {

   auto pred = Higgs::Predictions();
   namespace HP = Higgs::predictions;
   // whether to calculate the ggH cross-section in terms of the effective top and bottom Yukawa couplings
   // or by rescaling the SM-like ggH XS by the squared of the effective gg coupling (no effects from colored BSM particles are taken into account)
   static constexpr bool calcggH = false;
   // whether to calculate the H->gaga decay width in terms of the effective couplings
   // or by rescaling the SM-like H->gaga decay by the squared of the effective gamgam coupling (no effects from charged BSM particles are taken into account).
   static constexpr bool calcHgamgam = false;

   for (auto const& el : bsm_input) {

      const double mass = el.mass;
      // in the SM, λ = (mh/v)^2/2
      // for mh > 700 GeV this gives λ > 4
      // it probably makes no sense to use coupling strengh modifiers in this case so we skip those particles
      if (mass > 700) continue;

      auto& s = pred.addParticle(HP::BsmParticle(el.particle, HP::ECharge::neutral, static_cast<HP::CP>(el.CP)));
      s.setMass(mass);
      s.setMassUnc(0.03*mass); // set mass uncertainty to 3%

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
      auto match_Higgs_mass = [&sm, mass](double x) {
         sm.set_Lambdax(x);
         sm.calculate_DRbar_masses(); // internaly calls solve_ewsb_tree_level()
         sm.solve_ewsb();
         sm.calculate_Mhh_pole();
         return std::abs(sm.get_physical().Mhh - mass);
      };

      int status, iter = 0;
      static constexpr int max_iter = 50;
      const gsl_min_fminimizer_type *T;
      gsl_min_fminimizer *sGSL;
      // find λ in range [0, 5]
      double a = 0.0001, b = 5;
      // initial guess for the location of minimum: λ=(mass/v)^2/2
      double m = 0.5*Sqr(mass/247);

      // hack to pass lambda-function to GSL
      std::function<double(double)> f = std::bind(match_Higgs_mass, std::placeholders::_1);
      gsl_function F = {
         [](double d, void* vf) -> double {
            auto& f = *static_cast<std::function<double(double)>*>(vf);
            return f(d);
         },
         &f
      };

      // checked on a single point in the MRSSM2:
      //    brent seems faster and more accurate than quad_golden
      T = gsl_min_fminimizer_brent;
      sGSL = gsl_min_fminimizer_alloc (T);

      // gsl_min_fminimizer_set expects f(m) < f(a) && f(m) < f(b),
      // otherwise it returns GSL_EINVAL status.
      // In this case we randomly try different m from [max(a, 0.01m), min(b, 100m)]
      // until status != GSL_EINVAL
      gsl_error_handler_t * _error_handler = gsl_set_error_handler_off();
      status = gsl_min_fminimizer_set (sGSL, &F, m, a, b);
      if (status == GSL_EINVAL) {
         std::random_device rd;
         std::mt19937 gen(rd());
         std::uniform_real_distribution<> dis(std::max(a,1e-2*m), std::min(b,1e+2*m));
         do {
            m = dis(gen);
            status = gsl_min_fminimizer_set (sGSL, &F, m, a, b);
         } while (status == GSL_EINVAL);
      }
      gsl_set_error_handler (_error_handler);

      static constexpr double mass_precision = 1e-4;

      do
      {
           iter++;
           status = gsl_min_fminimizer_iterate (sGSL);

           m = gsl_min_fminimizer_x_minimum (sGSL);
           a = gsl_min_fminimizer_x_lower (sGSL);
           b = gsl_min_fminimizer_x_upper (sGSL);
      }
      while (std::abs(1. - sm.get_physical().Mhh/mass) > mass_precision && iter < max_iter);

      gsl_min_fminimizer_free (sGSL);

      if (const double diff = std::abs(1. - sm.get_physical().Mhh/mass); diff > mass_precision) {
         throw std::runtime_error("Higgstools interface: Cannot find a SM equivalent of " + el.particle +
                     " after " + std::to_string(iter+1) + "/" + std::to_string(max_iter) + " iterations. "
                     "Mass difference: " + std::to_string(mass) + " GeV (BSM) vs " +
                     std::to_string(sm.get_physical().Mhh) + " GeV (SM) for λSM = " + std::to_string(m) +
                     ". Difference: " + std::to_string(100*diff) + "%. ");
      }

      sm.calculate_pole_masses();

      if (sm.get_physical().Mhh > 0) {
         // calculate decays in the SM equivalent
         flexiblesusy::Standard_model_decays sm_decays(sm, qedqcd, physical_input, flexibledecay_settings);
         sm_decays.calculate_decays();
         const auto sm_input = sm_decays.get_higgstools_input();

         auto effc = HP::NeutralEffectiveCouplings {};

         // fermion channels are given as complex numbers
         // we normalize to real part of SM coupling
         // quarks
         effc.dd = std::abs(sm_input[0].dd) > 0 ? el.dd/sm_input[0].dd.real() : 0.;
         effc.uu = std::abs(sm_input[0].uu) > 0 ? el.uu/sm_input[0].uu.real() : 0.;
         effc.ss = std::abs(sm_input[0].ss) > 0 ? el.ss/sm_input[0].ss.real() : 0.;
         effc.cc = std::abs(sm_input[0].cc) > 0 ? el.cc/sm_input[0].cc.real() : 0.;
         effc.bb = std::abs(sm_input[0].bb) > 0 ? el.bb/sm_input[0].bb.real() : 0.;
         effc.tt = std::abs(sm_input[0].tt) > 0 ? el.tt/sm_input[0].tt.real() : 0.;
         // leptons
         effc.ee = std::abs(sm_input[0].ee)         > 0 ? el.ee/sm_input[0].ee.real()         : 0.;
         effc.mumu = std::abs(sm_input[0].mumu)     > 0 ? el.mumu/sm_input[0].mumu.real()     : 0.;
         effc.tautau = std::abs(sm_input[0].tautau) > 0 ? el.tautau/sm_input[0].tautau.real() : 0.;

         // gauge bosons
         effc.WW = std::abs(sm_input[0].WW) > 0         ? el.WW/sm_input[0].WW         : 0.;
         effc.ZZ = std::abs(sm_input[0].ZZ) > 0         ? el.ZZ/sm_input[0].ZZ         : 0.;
         effc.gamgam = std::abs(sm_input[0].gamgam) > 0 ? el.gamgam/sm_input[0].gamgam : 0.;
         effc.Zgam = std::abs(sm_input[0].Zgam) > 0     ? el.Zgam/sm_input[0].Zgam     : 0.;
         effc.gg = std::abs(sm_input[0].gg) > 0         ? el.gg/sm_input[0].gg         : 0.;

         effectiveCouplingInput(
            s, effc,
            HP::ReferenceModel::SMHiggsInterp,
            calcggH, calcHgamgam
         );

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

   /*
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
   */

   // HiggsBounds
   if (higgsbounds_dataset.empty()) {
      throw SetupError("Need to specify location of HiggsBounds database");
   }
   else if (!std::filesystem::exists(higgsbounds_dataset)) {
      throw SetupError("No HiggsBounds database found at " + higgsbounds_dataset);
   }
   auto bounds = Higgs::Bounds {higgsbounds_dataset};
   auto hbResult = bounds(pred);
   //std::ofstream hb_output("HiggsBounds.out");
   //hb_output << hbResult;
   //hb_output.close();
   // mimics the behavious of << operator
   std::vector<std::tuple<int, double, double, std::string>> hb_return {};
   for (const auto &[p, lim] : hbResult.selectedLimits) {
      auto found = std::find_if(std::begin(bsm_input), std::end(bsm_input), [&p](auto const& el) { return el.particle==p; });
      hb_return.push_back({found->pdgid, lim.obsRatio(), lim.expRatio(), lim.limit()->to_string()});
   }


   // HiggsSignals
   if (higgssignals_dataset.empty()) {
      throw SetupError("Need to specify location of HiggsSignals database");
   }
   else if (!std::filesystem::exists(higgssignals_dataset)) {
      throw SetupError("No HiggsSignals database found at " + higgssignals_dataset);
   }
   const auto signals = Higgs::Signals {higgssignals_dataset};
   //for (const auto &m : signals.measurements()) {
   //   std::cout << m.reference() << " " << m(pred) << std::endl;
   //}
   const double hs_chisq = signals(pred);

   return {signals.observableCount(), hs_chisq, hb_return};
}

} // flexiblesusy
