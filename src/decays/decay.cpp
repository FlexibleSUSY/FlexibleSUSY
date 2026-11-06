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

#include "decay.hpp"
#include "error.hpp"
#include "wrappers.hpp"

#include <boost/functional/hash.hpp>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <sstream>

#include "string_utils.hpp"

namespace flexiblesusy {

namespace {

template <class Container>
std::size_t hash_pid_list(int pid_in, Container pids_out)
{
   Container sorted(pids_out);
   std::sort(sorted.begin(), sorted.end());

   boost::hash<int> hash_pid;
   auto seed = hash_pid(pid_in);
   boost::hash_range(seed, sorted.begin(), sorted.end());

   return seed;
}

} // anonymous namespace

static std::array<int, 2> pdg_id_pairs[PDG_id_pairs::NUMBER_OF_PDG_IDS] = {
   {-1, 1},
   {-2, 2},
   {-3, 3},
   {-4, 4},
   {-5, 5},
   {-6, 6},
   {-11, 11},
   {-13, 13},
   {-15, 15},
   {-24, 24},
   {23, 23},
   {22, 22},
   {22, 23},
   {21, 21},
   {-11, 13},
   {-11, 15},
   {-13, 15}
};

std::size_t hash_decay(const Decay& decay)
{
   int pid_in = decay.get_initial_particle_id();
   const auto& pids_out = decay.get_final_state_particle_ids();
   return hash_pid_list(pid_in, pids_out);
}

Decays_list::Decays_list(int initial_pdg_)
   : initial_pdg(initial_pdg_)
{
}

void Decays_list::clear()
{
   decays.clear();
   total_width = 0.;
}

const Decay& Decays_list::get_decay(
   std::initializer_list<int> product_pdgs) const
{
   const Decay decay(initial_pdg, product_pdgs, 0., std::string());
   const auto decay_hash = hash_decay(decay);

   const auto pos = decays.find(decay_hash);

   if (pos == std::end(decays)) {
      std::ostringstream msg;
      msg << "Decay of particle " << initial_pdg
          << " into particles {"
          << concat(product_pdgs.begin(), product_pdgs.end(), ", ")
          << "} not found\n";

      throw OutOfBoundsError(msg.str());
   }

   return pos->second;
}

/**
 * Sort decays of every particle according to their width
 */
std::vector<Decay> sort_decays_list(const Decays_list& decays_list) {
   std::vector<Decay> decays_list_as_vector;
   decays_list_as_vector.reserve(decays_list.size());
   for (const auto& el : decays_list) {
      decays_list_as_vector.push_back(el.second);
   }

   std::sort(
      decays_list_as_vector.begin(),
      decays_list_as_vector.end(),
      [](const auto& d1, const auto& d2) {
         return d1.get_width() > d2.get_width();
      }
   );

   return decays_list_as_vector;
}

std::string strip_field_namespace(std::string const& s) {
   std::string result = s.substr(s.find_last_of(':')+1);
   if (s.find("bar") != std::string::npos) {
      result.pop_back();
      return "bar" + result;
   } else if (s.find("conj") != std::string::npos) {
      result.pop_back();
      return "conj" + result;
   } else {
      return result;
   }
}

namespace {

constexpr double hVV_4body_rel_tolerance = 1e-3;
constexpr std::size_t hVV_4body_max_intervals = 200;

struct HVV_4body_parameters {
   double mH{};
   double mV_sq{};
   double mV_GammaV{};
   gsl_integration_workspace* inner_workspace{nullptr};
   double q1_sq{};

   double t_from_q_sq(double q_sq) const { return std::atan((q_sq - mV_sq)/mV_GammaV); }
   double q_sq_from_t(double t) const { return mV_sq + mV_GammaV*std::tan(t); }
};

double hVV_4body_inner(double t2, void* params)
{
   const auto* p = static_cast<const HVV_4body_parameters*>(params);
   const double mH_sq = Sqr(p->mH);
   const double q2_sq = p->q_sq_from_t(t2);
   const double kl = KallenLambda(1., p->q1_sq/mH_sq, q2_sq/mH_sq);

   return kl > 0. ? std::sqrt(kl)*(kl + 12.*p->q1_sq*q2_sq/Power4(p->mH)) : 0.;
}

double hVV_4body_outer(double t1, void* params)
{
   auto* p = static_cast<HVV_4body_parameters*>(params);
   p->q1_sq = p->q_sq_from_t(t1);

   if (p->q1_sq <= 0. || p->q1_sq >= Sqr(p->mH)) {
      return 0.;
   }

   gsl_function inner;
   inner.function = &hVV_4body_inner;
   inner.params = params;

   double result = 0.;
   double error = 0.;
   gsl_integration_qag(
      &inner, p->t_from_q_sq(0.), p->t_from_q_sq(Sqr(p->mH - std::sqrt(p->q1_sq))),
      0., hVV_4body_rel_tolerance/3., hVV_4body_max_intervals, GSL_INTEG_GAUSS21,
      p->inner_workspace, &result, &error);

   return result;
}

} // anonymous namespace

/**
 * Integral over the virtualities of both vector bosons in H -> V* V*.
 *
 * The substitution q^2 = mV^2 + mV GammaV tan(t) turns each Breit-Wigner
 * factor times dq^2 into dt, which removes the resonance peaks from the
 * integrand.
 *
 * @param mH mass of the decaying Higgs boson
 * @param mV mass of the vector boson
 * @param GammaV width of the vector boson
 *
 * @return value of the integral and estimate of its absolute error
 */
HVV_4body_result hVV_4body(double mH, double mV, double GammaV)
{
   using Workspace = std::unique_ptr<gsl_integration_workspace, decltype(&gsl_integration_workspace_free)>;

   const Workspace outer_workspace(
      gsl_integration_workspace_alloc(hVV_4body_max_intervals), &gsl_integration_workspace_free);
   const Workspace inner_workspace(
      gsl_integration_workspace_alloc(hVV_4body_max_intervals), &gsl_integration_workspace_free);

   HVV_4body_parameters params{mH, Sqr(mV), mV*GammaV, inner_workspace.get()};

   gsl_function outer;
   outer.function = &hVV_4body_outer;
   outer.params = &params;

   HVV_4body_result result;

   // By default, GSL aborts if the requested precision is not reached.
   // The caller checks the error estimate instead.
   gsl_error_handler_t* error_handler = gsl_set_error_handler_off();
   const int status = gsl_integration_qag(
      &outer, params.t_from_q_sq(0.), params.t_from_q_sq(Sqr(mH)),
      0., hVV_4body_rel_tolerance, hVV_4body_max_intervals, GSL_INTEG_GAUSS21,
      outer_workspace.get(), &result.value, &result.error);
   gsl_set_error_handler(error_handler);

   if (status != GSL_SUCCESS) {
      result.error = std::max(result.error, std::abs(result.value));
   }

   return result;
}

void EffectiveCoupling_list::set_invisible_width(std::string const& p, double c) {
      auto found = std::find_if(
         std::begin(effective_coupling_list), std::end(effective_coupling_list),
         [&p](NeutralHiggsEffectiveCouplings const& effC) {return effC.particle == p;}
      );
      if (found == std::end(effective_coupling_list)) {
         auto effC = NeutralHiggsEffectiveCouplings {};
         effC.particle = p;
         effC.invWidth = c;
         effective_coupling_list.push_back(std::move(effC));
      }
      else {
         found->invWidth = c;
      }
}

void EffectiveCoupling_list::add_coupling(std::string const& p, std::array<int, 2> const& fs, std::pair<std::string, double> const& c, double partialW) {
      auto found = std::find_if(
         std::begin(effective_coupling_list), std::end(effective_coupling_list),
         [&p](NeutralHiggsEffectiveCouplings const& effC) {return effC.particle == p;}
      );
      auto are_the_same = [](std::array<int, 2> const& a1, std::array<int, 2> const& a2) {return std::is_permutation(a1.begin(), a1.end(), a2.begin(), a2.end());};
      if (found == std::end(effective_coupling_list)) {
         auto effC = NeutralHiggsEffectiveCouplings {};
         effC.particle = p;
         if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::gg])) {
            effC.gg = c;
            effC.undetectedWidth -= partialW;
         }
         else if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::gamgam])) {
            effC.gamgam = c;
            effC.undetectedWidth -= partialW;
         }
         else if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::ZZ])) {
            effC.ZZ = c;
            effC.undetectedWidth -= partialW;
         }
         else if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::WW])) {
            effC.WW = c;
            effC.undetectedWidth -= partialW;
         }
         else if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::Zgam])) {
            effC.Zgam = c;
            effC.undetectedWidth -= partialW;
         }
         effective_coupling_list.push_back(std::move(effC));
      }
      else {
         if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::gg])) {
            found->gg = c;
            found->undetectedWidth -= partialW;
         }
         else if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::gamgam])) {
            found->gamgam = c;
            found->undetectedWidth -= partialW;
         }
         else if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::ZZ])) {
            found->ZZ = c;
            found->undetectedWidth -= partialW;
         }
         else if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::WW])) {
            found->WW = c;
            found->undetectedWidth -= partialW;
         }
         else if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::Zgam])) {
            found->Zgam = c;
            found->undetectedWidth -= partialW;
         }
      }
   }

void EffectiveCoupling_list::add_coupling(std::string const& p, std::array<int, 2> const& fs, std::pair<std::string, std::complex<double>> const& c, double partialW) {

   auto found = std::find_if(
      std::begin(effective_coupling_list), std::end(effective_coupling_list),
      [&p](NeutralHiggsEffectiveCouplings const& effC) {return effC.particle == p;}
   );

   auto are_the_same = [](std::array<int, 2> const& a1, std::array<int, 2> const& a2) {
         return std::is_permutation(a1.begin(), a1.end(), a2.begin(), a2.end());
   };
   auto negate = [](std::array<int, 2> in) {
      std::for_each(in.begin(), in.end(), [](int& el){el *= -1; });
      return in;
   };

   if (found == std::end(effective_coupling_list)) {
      auto effC = NeutralHiggsEffectiveCouplings {};
      effC.particle = p;
      if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::dd])) {
         effC.dd = c;
         effC.undetectedWidth -= partialW;
      }
      else if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::uu])) {
         effC.uu = c;
         effC.undetectedWidth -= partialW;
      }
      else if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::ss])) {
         effC.ss = c;
         effC.undetectedWidth -= partialW;
      }
      else if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::cc])) {
         effC.cc = c;
         effC.undetectedWidth -= partialW;
      }
      else if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::bb])) {
         effC.bb = c;
         effC.undetectedWidth -= partialW;
      }
      else if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::tt])) {
         effC.tt = c;
         effC.undetectedWidth -= partialW;
      }
      else if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::ee])) {
         effC.ee = c;
         effC.undetectedWidth -= partialW;
      }
      else if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::mumu])) {
         effC.mumu = c;
         effC.undetectedWidth -= partialW;
      }
      else if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::tautau])) {
         effC.tautau = c;
         effC.undetectedWidth -= partialW;
      }
      else if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::emu]) || are_the_same(fs, negate(pdg_id_pairs[PDG_id_pairs::emu]))) {
         effC.emu.first = c.first;
         effC.emu.second += c.second;
         effC.undetectedWidth -= partialW;
      }
      else if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::etau]) || are_the_same(fs, negate(pdg_id_pairs[PDG_id_pairs::etau]))) {
         effC.etau.first = c.first;
         effC.etau.second += c.second;
         effC.undetectedWidth -= partialW;
      }
      else if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::mutau]) || are_the_same(fs, negate(pdg_id_pairs[PDG_id_pairs::mutau]))) {
         effC.mutau.first = c.first;
         effC.mutau.second = c.second;
         effC.undetectedWidth -= partialW;
      }
      else {
         WARNING("HiggsTools interface warning: trying to add an unknown decay channel {" + std::to_string(fs[0]) + ", " + std::to_string(fs[1]) + "}");
      }
      effective_coupling_list.push_back(std::move(effC));
   }
   else {
      if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::dd])) {
         found->dd = c;
         found->undetectedWidth -= partialW;
      }
      else if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::uu])) {
         found->uu = c;
         found->undetectedWidth -= partialW;
      }
      else if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::ss])) {
         found->ss = c;
         found->undetectedWidth -= partialW;
      }
      else if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::cc])) {
         found->cc = c;
         found->undetectedWidth -= partialW;
      }
      else if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::bb])) {
         found->bb = c;
         found->undetectedWidth -= partialW;
      }
      else if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::tt])) {
         found->tt = c;
         found->undetectedWidth -= partialW;
      }
      else if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::ee])) {
         found->ee = c;
         found->undetectedWidth -= partialW;
      }
      else if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::mumu])) {
         found->mumu = c;
         found->undetectedWidth -= partialW;
      }
      else if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::tautau])) {
         found->tautau = c;
         found->undetectedWidth -= partialW;
      }
      else if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::emu]) || are_the_same(fs, negate(pdg_id_pairs[PDG_id_pairs::emu]))) {
         found->emu.first = c.first;
         found->emu.second += c.second;
         found->undetectedWidth -= partialW;
      }
      else if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::etau]) || are_the_same(fs, negate(pdg_id_pairs[PDG_id_pairs::etau]))) {
         found->etau.first = c.first;
         found->etau.second += c.second;
         found->undetectedWidth -= partialW;
      }
      else if (are_the_same(fs, pdg_id_pairs[PDG_id_pairs::mutau]) || are_the_same(fs, negate(pdg_id_pairs[PDG_id_pairs::mutau]))) {
         found->mutau.first = c.first;
         found->mutau.second += c.second;
         found->undetectedWidth -= partialW;
      }
      else {
         WARNING("HiggsTools interface warning: trying to add an unknown decay channel {" + std::to_string(fs[0]) + ", " + std::to_string(fs[1]) + "}");
      }
   }
}

} // namespace flexiblesusy
