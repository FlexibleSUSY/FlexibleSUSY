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

void Decays_list::set_decay(double width, std::initializer_list<int> pids_out, std::string proc_string)
{
   const Decay decay(initial_pdg, pids_out, width, proc_string);
   const auto decay_hash = hash_decay(decay);

   const auto pos = decays.find(decay_hash);
   if (pos != std::end(decays)) {
      total_width -= pos->second.get_width();
      pos->second.set_width(width);
   } else {
      decays.insert(pos, std::make_pair(decay_hash, decay));
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
          << " into particles {";

      std::size_t count = 0;
      const std::size_t n_final = product_pdgs.size();
      for (int pdg : product_pdgs) {
         count++;
         msg << pdg;
         if (++count != n_final) {
            msg << ", ";
         }
      }

      msg << "} not found\n";

      throw OutOfBoundsError(msg.str());
   }

   return pos->second;
}

} // namespace flexiblesusy
