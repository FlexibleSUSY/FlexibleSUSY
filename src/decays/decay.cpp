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

#include <boost/functional/hash.hpp>

#include <algorithm>
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

std::size_t hash_decay(const Decay& decay)
{
   int pid_in = decay.get_initial_particle_id();
   const auto& pids_out = decay.get_final_state_particle_ids();
   return hash_pid_list(pid_in, pids_out);
}

Decay::Decay(
   int pid_in_, std::initializer_list<int> pids_out_, double width_, std::string const& proc_string_)
   : pid_in(pid_in_)
   , pids_out(pids_out_)
   , width(width_)
   , proc_string(proc_string_)
{
   std::sort(pids_out.begin(), pids_out.end());
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

void Decays_list::set_decay(double width, std::initializer_list<int> pids_out, std::string const& proc_string)
{
   const Decay decay(initial_pdg, pids_out, width, proc_string);
   const auto decay_hash = hash_decay(decay);

   const auto pos = decays.find(decay_hash);
   if (pos != std::end(decays)) {
      total_width -= pos->second.get_width();
      pos->second.set_width(width);
   } else {
      decays.insert(pos, std::make_pair(decay_hash, std::move(decay)));
   }

   // some channels give small negative withs
   // we later check if for channels with width < 0
   // |width/total_width| < threshold
   // for that it makes more sense to calculate total_width
   // form sum of |width|
   total_width += std::abs(width);
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

double hVV_4body(double *q2, size_t dim, void *params)
{
  (void)(dim); /* avoid unused parameter warnings */
  struct my_f_params * fp = (struct my_f_params *)params;
  const double mHOS = fp->mHOS;
  if (q2[1] > Sqr(mHOS - std::sqrt(q2[0]))) return 0.;
  const double mVOS = fp->mVOS;
  const double GammaV = fp->GammaV;
  const double kl = KallenLambda(1., q2[0]/Sqr(mHOS), q2[1]/Sqr(mHOS));
  return
     mVOS*GammaV/(Sqr(q2[0] - Sqr(mVOS)) + Sqr(mVOS*GammaV))
     * mVOS*GammaV/(Sqr(q2[1] - Sqr(mVOS)) + Sqr(mVOS*GammaV))
     * std::sqrt(kl)*(kl + 12.*q2[0]*q2[1]/Power4(mHOS));
}

void EffectiveCoupling_list::add_coupling(std::string const& p, std::array<int, 2> const& fs, double c) {
      std::cout << "double: " << fs.at(0) << ' ' << fs.at(1) << std::endl;
      auto found = std::find_if(
         std::begin(effective_coupling_list), std::end(effective_coupling_list),
         [&p](NeutralHiggsEffectiveCouplings const& effC) {return effC.particle == p;}
      );
      auto are_the_same = [](std::array<int, 2> const& a1, std::array<int, 2> const& a2) {return std::is_permutation(a1.begin(), a1.end(), a2.begin(), a2.end());};
      if (found == std::end(effective_coupling_list)) {
         auto effC = NeutralHiggsEffectiveCouplings {};
         effC.particle = p;
         if (are_the_same(fs, {21, 21})) {
            effC.gg = c;
         }
         else if (are_the_same(fs, {22, 22})) {
            effC.gamgam = c;
         }
         else if (are_the_same(fs, {23, 23})) {
            effC.ZZ = c;
         }
         else if (are_the_same(fs, {-24, 24})) {
            effC.WW = c;
         }
         else if (are_the_same(fs, {22, 23})) {
            effC.Zgam = c;
         }
         effective_coupling_list.push_back(std::move(effC));
      }
      else {
         if (are_the_same(fs, {21, 21})) {
            found->gg = c;
         }
         else if (are_the_same(fs, {22, 22})) {
            found->gamgam = c;
         }
         else if (are_the_same(fs, {23, 23})) {
            found->ZZ = c;
         }
         else if (are_the_same(fs, {-24, 24})) {
            found->WW = c;
         }
         else if (are_the_same(fs, {22, 23})) {
            found->Zgam = c;
         }
      }
   }
void EffectiveCoupling_list::add_coupling(std::string const& p, std::array<int, 2> const& fs, std::complex<double> c) {
   std::cout << "complex: " << fs.at(0) << ' ' << fs.at(1) << std::endl;

   auto found = std::find_if(
      std::begin(effective_coupling_list), std::end(effective_coupling_list),
      [&p](NeutralHiggsEffectiveCouplings const& effC) {return effC.particle == p;}
   );

   auto are_the_same = [](std::array<int, 2> const& a1, std::array<int, 2> const& a2) {
         return std::is_permutation(a1.begin(), a1.end(), a2.begin(), a2.end());
   };

   if (found == std::end(effective_coupling_list)) {
      auto effC = NeutralHiggsEffectiveCouplings {};
      effC.particle = p;
      if (are_the_same(fs, {-1, 1})) {
         effC.dd = c;
      }
      else if (are_the_same(fs, {-2, 2})) {
         effC.uu = c;
      }
      else if (are_the_same(fs, {-3, 3})) {
         effC.ss = c;
      }
      else if (are_the_same(fs, {-4, 4})) {
         effC.cc = c;
      }
      else if (are_the_same(fs, {-5, 5})) {
         effC.bb = c;
      }
      else if (are_the_same(fs, {-6, 6})) {
         effC.tt = c;
      }
      else if (are_the_same(fs, {-11, 11})) {
         effC.ee = c;
      }
      else if (are_the_same(fs, {-13, 13})) {
         effC.mumu = c;
      }
      else if (are_the_same(fs, {-15, 15})) {
         effC.tautau = c;
      }
      else if (are_the_same(fs, {-11, 13}) || are_the_same(fs, {11, -13})) {
         effC.emu += c;
      }
      else if (are_the_same(fs, {-11, 15}) || are_the_same(fs, {11, -15})) {
         effC.etau += c;
      }
      else if (are_the_same(fs, {-13, 15}) || are_the_same(fs, {13, -15})) {
         effC.mutau += c;
      }
      effective_coupling_list.push_back(std::move(effC));
   }
   else {
      if (are_the_same(fs, {-1, 1})) {
         found->dd = c;
      }
      else if (are_the_same(fs, {-2, 2})) {
         found->uu = c;
      }
      else if (are_the_same(fs, {-3, 3})) {
         found->ss = c;
      }
      else if (are_the_same(fs, {-4, 4})) {
         found->cc = c;
      }
      else if (are_the_same(fs, {-5, 5})) {
         found->bb = c;
      }
      else if (are_the_same(fs, {-6, 6})) {
         found->tt = c;
      }
      else if (are_the_same(fs, {-11, 11})) {
         found->ee = c;
      }
      else if (are_the_same(fs, {-13, 13})) {
         found->mumu = c;
      }
      else if (are_the_same(fs, {-15, 15})) {
         found->tautau = c;
      }
      else if (are_the_same(fs, {-11, 13}) || are_the_same(fs, {11, -13})) {
         found->emu += c;
      }
      else if (are_the_same(fs, {-11, 15}) || are_the_same(fs, {11, -15})) {
         found->etau += c;
      }
      else if (are_the_same(fs, {-13, 15}) || are_the_same(fs, {13, -15})) {
         found->mutau += c;
      }
   }
}

} // namespace flexiblesusy
