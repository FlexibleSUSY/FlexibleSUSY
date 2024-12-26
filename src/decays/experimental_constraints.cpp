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

#include "experimental_constraints.hpp"

#ifdef ENABLE_HIGGSTOOLS
#include "Higgs/Predictions.hpp"
#include "Higgs/Bounds.hpp"
#include "Higgs/Signals.hpp"
namespace HP = Higgs::predictions;
#endif

#ifdef ENABLE_LILITH
#include <Python.h>
#include "lilith.h"
#include "lilith.c"
#endif

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

namespace {

// Whether to calculate the X->gaga and X->gg decay widths in terms of the
// effective tree-level couplings (calcHgamgam=calcggH=true) or by rescaling
// the X->gaga and X->gg decay by the squared of the effective gamgam and gg
// couplings (calcHgamgam=calcggH=false). True would make sense only if we
// would not compute those loop-induced decays ourself.
constexpr bool calcggH     = false;
constexpr bool calcHgamgam = false;

// relative BSM Higgs-like state mass uncertainty
// example: 0.03 means 3% uncertainty
constexpr double relMassError = 0.03;

#ifdef ENABLE_HIGGSTOOLS
// Ref. model for computing brs and xsections on the HiggsTools side
constexpr auto refModel = HP::ReferenceModel::SMHiggsInterp;

double minChi2SM_hs(const double mhSM, std::string const& higgssignals_dataset) {
   const auto signals = Higgs::Signals {higgssignals_dataset};

   auto pred = Higgs::Predictions();
   auto& s = pred.addParticle(HP::BsmParticle("hSM", HP::ECharge::neutral, HP::CP::even));
   auto effc = HP::scaledSMlikeEffCouplings(1.0);
   s.setMass(mhSM);
   s.setMassUnc(0.);
   effectiveCouplingInput(
      s, effc,
      refModel,
      calcggH, calcHgamgam
   );
   return signals(pred);
}

void print_effc(double mass, HP::NeutralEffectiveCouplings const& effC) {
   std::cout << "Effective couplings for particle of mass " << mass << '\n';
   std::cout << "dd     " << effC.dd << std::endl;
   std::cout << "uu     " << effC.uu << std::endl;
   std::cout << "ss     " << effC.ss << std::endl;
   std::cout << "cc     " << effC.cc << std::endl;
   std::cout << "bb     " << effC.bb << std::endl;
   std::cout << "tt     " << effC.tt << std::endl;
   std::cout << "ee     " << effC.ee << std::endl;
   std::cout << "mumu   " << effC.mumu << std::endl;
   std::cout << "tautau " << effC.tautau << std::endl;

   std::cout << "WW     " << effC.WW << std::endl;
   std::cout << "ZZ     " << effC.ZZ << std::endl;
   std::cout << "gamgam " << effC.gamgam << std::endl;
   std::cout << "Zgam   " << effC.Zgam << std::endl;
   std::cout << "gg     " << effC.gg << std::endl;
   std::cout << "lam     " << effC.lam << std::endl;
}
#endif

#ifdef ENABLE_LILITH
double minChi2SM_lilith(const double mhSM) {

   Py_Initialize();
   char experimental_input[] = "";
   // Creating an object of the class Lilith: lilithcalc
   PyObject* lilithcalc2 = initialize_lilith(experimental_input);

   char XMLinputstring[6000] = "";
   char buffer[100];

   sprintf(buffer,"<?xml version=\"1.0\"?>\n");
   strcat(XMLinputstring, buffer);
   sprintf(buffer,"<lilithinput>\n");
   strcat(XMLinputstring, buffer);

   constexpr double cSM = 1.;
   constexpr double BRinv = 0.;
   constexpr double BRund = 0.;

   sprintf(buffer,"<reducedcouplings>\n");
   strcat(XMLinputstring, buffer);

   sprintf(buffer,"<mass>%f</mass>\n", mhSM);
   strcat(XMLinputstring, buffer);

   sprintf(buffer,"<C to=\"gammagamma\">%f</C>\n", cSM);
   strcat(XMLinputstring, buffer);
   sprintf(buffer,"<C to=\"Zgamma\">%f</C>\n", cSM);
   strcat(XMLinputstring, buffer);
   sprintf(buffer,"<C to=\"gg\">%f</C>\n", cSM);
   strcat(XMLinputstring, buffer);

   sprintf(buffer,"<C to=\"ZZ\">%f</C>\n", cSM);
   strcat(XMLinputstring, buffer);
   sprintf(buffer,"<C to=\"WW\">%f</C>\n", cSM);
   strcat(XMLinputstring, buffer);

   sprintf(buffer,"<C to=\"tt\" part=\"re\">%f</C>\n", cSM);
   strcat(XMLinputstring, buffer);
   sprintf(buffer,"<C to=\"cc\" part=\"re\">%f</C>\n", cSM);
   strcat(XMLinputstring, buffer);
   sprintf(buffer,"<C to=\"bb\" part=\"re\">%f</C>\n", cSM);
   strcat(XMLinputstring, buffer);
   sprintf(buffer,"<C to=\"tautau\" part=\"re\">%f</C>\n", cSM);
   strcat(XMLinputstring, buffer);

   sprintf(buffer,"<extraBR>\n");
   strcat(XMLinputstring, buffer);
   sprintf(buffer,"<BR to=\"invisible\">%f</BR>\n", BRinv);
   strcat(XMLinputstring, buffer);
   sprintf(buffer,"<BR to=\"undetected\">%f</BR>\n", BRund);
   strcat(XMLinputstring, buffer);
   sprintf(buffer,"</extraBR>\n");
   strcat(XMLinputstring, buffer);
   sprintf(buffer,"</reducedcouplings>\n");
   strcat(XMLinputstring, buffer);

   sprintf(buffer,"</lilithinput>\n");
   strcat(XMLinputstring, buffer);

   // Reading user input XML string
   lilith_readuserinput(lilithcalc2, XMLinputstring);

   // Getting -2LogL
   const double my_likelihood = lilith_computelikelihood(lilithcalc2);

   Py_Finalize();
   return my_likelihood;
}
#endif

void set_sm_settings_matching_bsm(
   standard_model::Standard_model& sm,
   Spectrum_generator_settings const& spectrum_generator_settings)
{
  sm.set_pole_mass_loop_order(static_cast<int>(spectrum_generator_settings.get(Spectrum_generator_settings::pole_mass_loop_order)));
  sm.set_ewsb_loop_order(static_cast<int>(spectrum_generator_settings.get(Spectrum_generator_settings::ewsb_loop_order)));
  sm.set_precision(spectrum_generator_settings.get(Spectrum_generator_settings::precision));
  sm.set_threshold_corrections(spectrum_generator_settings.get_threshold_corrections());
  sm.set_loop_corrections(spectrum_generator_settings.get_loop_corrections());
  sm.set_loops(static_cast<int>(spectrum_generator_settings.get(Spectrum_generator_settings::beta_loop_order)));

  Loop_corrections loop_corrections_;
  loop_corrections_.higgs_at_as = spectrum_generator_settings.get(Spectrum_generator_settings::higgs_2loop_correction_at_as);
  loop_corrections_.higgs_ab_as = spectrum_generator_settings.get(Spectrum_generator_settings::higgs_2loop_correction_ab_as);
  loop_corrections_.higgs_at_at = spectrum_generator_settings.get(Spectrum_generator_settings::higgs_2loop_correction_at_at);
  loop_corrections_.higgs_atau_atau = spectrum_generator_settings.get(Spectrum_generator_settings::higgs_2loop_correction_atau_atau);
  loop_corrections_.top_qcd = spectrum_generator_settings.get(Spectrum_generator_settings::top_pole_qcd_corrections);
  loop_corrections_.higgs_at_as_as = spectrum_generator_settings.get(Spectrum_generator_settings::higgs_3loop_correction_at_as2);
  loop_corrections_.higgs_ab_as_as = spectrum_generator_settings.get(Spectrum_generator_settings::higgs_3loop_correction_ab_as2);
  loop_corrections_.higgs_at_at_as = spectrum_generator_settings.get(Spectrum_generator_settings::higgs_3loop_correction_at2_as);
  loop_corrections_.higgs_at_at_at = spectrum_generator_settings.get(Spectrum_generator_settings::higgs_3loop_correction_at3);
  loop_corrections_.higgs_at_as_as_as = spectrum_generator_settings.get(Spectrum_generator_settings::higgs_4loop_correction_at_as3);
  sm.set_loop_corrections(loop_corrections_);

  sm.set_threshold_corrections(spectrum_generator_settings.get_threshold_corrections());
}

void set_sm_lambda_to_match_bsm_mh(standard_model::Standard_model& sm, double mass, std::string const& name) {

      // set SM λ such that mhSM == mass
      auto match_Higgs_mass = [&sm, mass](double x) {
         sm.set_Lambdax(x);
         sm.calculate_DRbar_masses(); // internaly calls solve_ewsb_tree_level()
         sm.solve_ewsb();
         sm.calculate_Mhh_pole();
         return std::abs(sm.get_physical().Mhh - mass);
      };

      int status, iter = 0;
      static constexpr int max_iter = 100;
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
         int iterCount = 0;
         do {
            m = dis(gen);
            status = gsl_min_fminimizer_set (sGSL, &F, m, a, b);
            iterCount++;
         } while (status == GSL_EINVAL && iterCount < 100);
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
         throw std::runtime_error("Normalized Higgs effective couplings: Cannot find a SM equivalent of " + name +
                     " after " + std::to_string(iter+1) + "/" + std::to_string(max_iter) + " iterations. "
                     "Mass difference: " + std::to_string(mass) + " GeV (BSM) vs " +
                     std::to_string(sm.get_physical().Mhh) + " GeV (SM) for λSM = " + std::to_string(m) +
                     ". Difference: " + std::to_string(100*diff) + "%. ");
      }


}
} // anonymous

EffectiveCoupling_list get_normalized_effective_couplings(
   EffectiveCoupling_list const& bsm_input,
   Physical_input const& physical_input,
   softsusy::QedQcd const& qedqcd,
   Spectrum_generator_settings const& spectrum_generator_settings,
   FlexibleDecay_settings const& flexibledecay_settings)
{
   // make sure we don't compute input for the EFFHIGGSCOUPLINGS block in the
   // built in SM
   auto flexibledecay_settings_ = flexibledecay_settings;

   EffectiveCoupling_list _bsm_input;
   for (auto const& el : bsm_input) {

      const double mass = el.mass;
      // in the SM, λ = (mh/v)^2/2
      // for mh > 700 GeV this gives λ > 4
      // it probably makes no sense to use coupling strengh modifiers in this case so we skip those particles
      // On the other hand there's a problem with finding a SM equivalent of very light states
      if (mass > 650 || mass < 1) continue;

      // create a SM equivalent to the BSM model, with mhSM == mass
      standard_model::Standard_model sm {};
      set_sm_settings_matching_bsm(sm, spectrum_generator_settings);
      sm.set_physical_input(physical_input);
      sm.initialise_from_input(qedqcd);
      set_sm_lambda_to_match_bsm_mh(sm, mass, el.particle);
      sm.calculate_pole_masses();

      if (sm.get_physical().Mhh > 0) {
         // calculate decays in the SM equivalent
         flexiblesusy::Standard_model_decays sm_decays(sm, qedqcd, physical_input, flexibledecay_settings_);
         sm_decays.calculate_decays();
         const auto sm_input = sm_decays.get_neutral_higgs_effc();

         // fermion channels are given as complex numbers
         // we normalize to real part of SM coupling
         // quarks
         NeutralHiggsEffectiveCouplings _coups {el};
         _coups.width_sm = sm_input[0].width;
         _coups.dd.second = std::abs(sm_input[0].dd.second) > 0 ? el.dd.second/sm_input[0].dd.second.real() : 0.;
         _coups.uu.second = std::abs(sm_input[0].uu.second) > 0 ? el.uu.second/sm_input[0].uu.second.real() : 0.;
         _coups.ss.second = std::abs(sm_input[0].ss.second) > 0 ? el.ss.second/sm_input[0].ss.second.real() : 0.;
         _coups.cc.second = std::abs(sm_input[0].cc.second) > 0 ? el.cc.second/sm_input[0].cc.second.real() : 0.;
         _coups.bb.second = std::abs(sm_input[0].bb.second) > 0 ? el.bb.second/sm_input[0].bb.second.real() : 0.;
         using namespace std::complex_literals;
         _coups.tt.second = std::abs(sm_input[0].tt.second) > 0 ? (std::abs(el.tt.second.real()) + 1i*std::abs(el.tt.second.imag()))/std::abs(sm_input[0].tt.second.real()) : 0.;
         // leptons
         _coups.ee.second = std::abs(sm_input[0].ee.second)         > 0 ? el.ee.second/sm_input[0].ee.second.real()         : 0.;
         _coups.mumu.second = std::abs(sm_input[0].mumu.second)     > 0 ? el.mumu.second/sm_input[0].mumu.second.real()     : 0.;
         _coups.tautau.second = std::abs(sm_input[0].tautau.second) > 0 ? el.tautau.second/sm_input[0].tautau.second.real() : 0.;

         // gauge bosons
         _coups.WW.second = std::abs(sm_input[0].WW.second)         > 0 ? el.WW.second/sm_input[0].WW.second         : 0.;
         _coups.ZZ.second = std::abs(sm_input[0].ZZ.second)         > 0 ? el.ZZ.second/sm_input[0].ZZ.second         : 0.;
         _coups.gamgam.second = std::abs(sm_input[0].gamgam.second) > 0 ? el.gamgam.second/sm_input[0].gamgam.second : 0.;
         _coups.Zgam.second = std::abs(sm_input[0].Zgam.second)     > 0 ? el.Zgam.second/sm_input[0].Zgam.second     : 0.;
         _coups.gg.second = std::abs(sm_input[0].gg.second)         > 0 ? el.gg.second/sm_input[0].gg.second         : 0.;

         _coups.lam = std::abs(sm_input[0].lam)       > 0 ? el.lam/sm_input[0].lam       : 0.;
         _bsm_input.push_back(std::move(_coups));
      }
   }

   return _bsm_input;
}

#ifdef ENABLE_HIGGSTOOLS
std::tuple<SignalResult, std::vector<std::tuple<int, double, double, std::string>>> call_higgstools(
   EffectiveCoupling_list const& bsm_input,
   Physical_input const& physical_input,
   std::string const& higgsbounds_dataset, std::string const& higgssignals_dataset) {

   // check location of databases
   // HiggsBounds
   if (higgsbounds_dataset.empty()) {
      throw SetupError("Need to specify location of HiggsBounds database");
   }
   else if (!std::filesystem::exists(higgsbounds_dataset)) {
      throw SetupError("No HiggsBounds database found at " + higgsbounds_dataset);
   }
   // HiggsSignals
   if (higgssignals_dataset.empty()) {
      throw SetupError("Need to specify location of HiggsSignals database");
   }
   else if (!std::filesystem::exists(higgssignals_dataset)) {
      throw SetupError("No HiggsSignals database found at " + higgssignals_dataset);
   }

   auto pred = Higgs::Predictions();
   namespace HP = Higgs::predictions;

   for (auto const& el : bsm_input) {
      auto effc = HP::NeutralEffectiveCouplings {};
      auto& s = pred.addParticle(HP::BsmParticle(el.particle, HP::ECharge::neutral, static_cast<HP::CP>(el.CP)));
      s.setMass(el.mass);
      s.setMassUnc(relMassError*el.mass); // set mass uncertainty to 3%

      // quarks
      effc.dd = el.dd.second;
      effc.uu = el.uu.second;
      effc.ss = el.ss.second;
      effc.cc = el.cc.second;
      effc.bb = el.bb.second;
      effc.tt = el.tt.second;

      // leptons
      effc.ee = el.ee.second;
      effc.mumu = el.mumu.second;
      effc.tautau = el.tautau.second;

      // gauge bosons
      effc.WW = el.WW.second;
      effc.ZZ = el.ZZ.second;
      effc.gamgam = el.gamgam.second;
      effc.Zgam = el.Zgam.second;
      effc.gg = el.gg.second;

      effc.lam = el.lam;

      effectiveCouplingInput(s, effc, refModel, calcggH, calcHgamgam);

      // Effective coupligs below are defined as sqrt(Gamma CP-even) + I sqrt(Gamma CP-odd)
      // (note that this is different than couplings like gg, WW etc)
      // so taking a norm gives a total partial width
      s.setDecayWidth(HP::Decay::emu,   std::norm(el.emu.second));
      s.setDecayWidth(HP::Decay::etau,  std::norm(el.etau.second));
      s.setDecayWidth(HP::Decay::mutau, std::norm(el.mutau.second));

      // Higgs to LSP decay (if model contains one)
      s.setDecayWidth("Inv", "Inv", el.invWidth);

      // all remaining partial widths
      s.setDecayWidth("Undetected", "Undetected", el.get_undetected_width());
   }

   auto bounds = Higgs::Bounds {higgsbounds_dataset};
   auto hbResult = bounds(pred);
   std::vector<std::tuple<int, double, double, std::string>> hb_return {};
   for (auto const& _hb: hbResult.selectedLimits) {
      auto found = std::find_if(
         std::begin(bsm_input), std::end(bsm_input),
         [&_hb](auto const& el) { return el.particle==_hb.first; }
      );
      hb_return.push_back({found->pdgid, _hb.second.obsRatio(), _hb.second.expRatio(), _hb.second.limit()->to_string()});
   }

   const auto signals = Higgs::Signals {higgssignals_dataset};
   //for (const auto &m : signals.measurements()) {
   //   std::cout << m.reference() << " " << m(pred) << std::endl;
   //}
   const double hs_chisq = signals(pred);

   const double mhSMref = physical_input.get(Physical_input::mh_pole);

   auto smChi2 = minChi2SM_hs(mhSMref, higgssignals_dataset);

   return {{signals.observableCount(), mhSMref, hs_chisq, smChi2}, hb_return};
}
#endif

#ifdef ENABLE_LILITH
std::optional<SignalResult> call_lilith(
   EffectiveCoupling_list const& bsm_input,
   Physical_input const& physical_input,
   std::string const& lilith_db) {

   // Lilith requires Higgs mass to be in range [123, 128]
   bool higgs_in_range = false;
   for (auto const& el : bsm_input) {
      const double mh = el.mass;
      if (mh > 123.0 && mh < 128.0) {
         higgs_in_range = true || higgs_in_range;
      }
   }
   if (!higgs_in_range) {
      return {};
   }

   Py_Initialize();
   // Creating an object of the class Lilith: lilithcalc
   PyObject* lilithcalc = initialize_lilith(const_cast<char*>(lilith_db.c_str()));

   char XMLinputstring[6000]="";
   char buffer[100];

   sprintf(buffer,"<?xml version=\"1.0\"?>\n");
   strcat(XMLinputstring, buffer);
   sprintf(buffer,"<lilithinput>\n");
   strcat(XMLinputstring, buffer);

   for (auto const& el : bsm_input) {

      const double mh = el.mass;
      if (mh < 123.0 || mh > 128.0) continue;

      const double BRinv = el.invWidth/el.width;
      const double BRund = el.get_undetected_width()/el.width;

      sprintf(buffer,"<reducedcouplings>\n");
      strcat(XMLinputstring, buffer);

      sprintf(buffer,"<mass>%f</mass>\n", mh);
      strcat(XMLinputstring, buffer);

      // massless gauge bosons
      sprintf(buffer,"<C to=\"gammagamma\">%f</C>\n", el.gamgam.second);
      strcat(XMLinputstring, buffer);
      sprintf(buffer,"<C to=\"Zgamma\">%f</C>\n", el.Zgam.second);
      strcat(XMLinputstring, buffer);
      // the same reduce coupling for production and decay for gluons
      sprintf(buffer,"<C to=\"gg\">%f</C>\n", el.gg.second);
      strcat(XMLinputstring, buffer);

      // massive gauge bosons
      sprintf(buffer,"<C to=\"ZZ\">%f</C>\n", el.ZZ.second);
      strcat(XMLinputstring, buffer);
      sprintf(buffer,"<C to=\"WW\">%f</C>\n", el.WW.second);
      strcat(XMLinputstring, buffer);

      // fermions
      // tt
      sprintf(buffer,"<C to=\"tt\" part=\"re\">%f</C>\n", std::real(el.tt.second));
      strcat(XMLinputstring, buffer);
      sprintf(buffer,"<C to=\"tt\" part=\"im\">%f</C>\n", std::imag(el.tt.second));
      strcat(XMLinputstring, buffer);
      // cc
      sprintf(buffer,"<C to=\"cc\" part=\"re\">%f</C>\n", std::real(el.cc.second));
      strcat(XMLinputstring, buffer);
      sprintf(buffer,"<C to=\"cc\" part=\"im\">%f</C>\n", std::imag(el.cc.second));
      strcat(XMLinputstring, buffer);
      // bb
      sprintf(buffer,"<C to=\"bb\" part=\"re\">%f</C>\n", std::real(el.bb.second));
      strcat(XMLinputstring, buffer);
      sprintf(buffer,"<C to=\"bb\" part=\"im\">%f</C>\n", std::imag(el.bb.second));
      strcat(XMLinputstring, buffer);
      // tautau
      sprintf(buffer,"<C to=\"tautau\" part=\"re\">%f</C>\n", std::real(el.tautau.second));
      strcat(XMLinputstring, buffer);
      sprintf(buffer,"<C to=\"tautau\" part=\"im\">%f</C>\n", std::imag(el.tautau.second));
      strcat(XMLinputstring, buffer);

      sprintf(buffer,"<extraBR>\n");
      strcat(XMLinputstring, buffer);
         sprintf(buffer,"<BR to=\"invisible\">%f</BR>\n", BRinv);
         strcat(XMLinputstring, buffer);
         sprintf(buffer,"<BR to=\"undetected\">%f</BR>\n", BRund);
         strcat(XMLinputstring, buffer);
      sprintf(buffer,"</extraBR>\n");
      strcat(XMLinputstring, buffer);

      // in the BEST-QCD mode, only the real part of the coupling is taken into account (see 1502.04138)
      sprintf(buffer,"<precision>%s</precision>\n", el.CP == 1 ? "BEST-QCD" : "LO");
      strcat(XMLinputstring, buffer);
      sprintf(buffer,"</reducedcouplings>\n");
      strcat(XMLinputstring, buffer);
    }

    sprintf(buffer,"</lilithinput>\n");
    strcat(XMLinputstring, buffer);

    // reading user input XML string
    lilith_readuserinput(lilithcalc, XMLinputstring);

    // getting -2*log(L)
    const double my_likelihood = lilith_computelikelihood(lilithcalc);

    // getting ndf
    const std::size_t exp_ndf = static_cast<std::size_t>(lilith_exp_ndf(lilithcalc));

    const double mhSMref = physical_input.get(Physical_input::mh_pole);

    const double sm_likelihood = minChi2SM_lilith(mhSMref);

    Py_Finalize();

    const SignalResult res {exp_ndf, mhSMref, my_likelihood, sm_likelihood};
    return res;
}
#endif

} // flexiblesusy
