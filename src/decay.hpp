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

#ifndef DECAY_H
#define DECAY_H

#include <initializer_list>
#include <unordered_map>
#include <vector>

namespace flexiblesusy {

class Decay {
public:
   Decay(int, std::initializer_list<int>, double);
   ~Decay() = default;
   Decay(const Decay&) = default;
   Decay(Decay&&) = default;
   Decay& operator=(const Decay&) = default;
   Decay& operator=(Decay&&) = default;

   int get_initial_particle_id() const { return pid_in; }
   const std::vector<int>& get_final_state_particle_ids() const {
      return pids_out;
   }
   std::size_t get_final_state_size() const { return pids_out.size(); }

   double get_width() const { return width; }
   void set_width(double w) { width = w; }

private:
   int pid_in{0};
   std::vector<int> pids_out{};
   double width{0.};
};

std::size_t hash_decay(const Decay& decay);

class Decays_list {
private:
   using List_type = std::unordered_map<std::size_t, Decay>;
public:
   using iterator = List_type::iterator;
   using const_iterator = List_type::const_iterator;

   explicit Decays_list(int);
   ~Decays_list() = default;
   Decays_list(const Decays_list&) = default;
   Decays_list(Decays_list&&) = default;
   Decays_list& operator=(const Decays_list&) = default;
   Decays_list& operator=(Decays_list&&) = default;

   iterator begin() noexcept { return decays.begin(); }
   const_iterator begin() const noexcept { return decays.begin(); }
   const_iterator cbegin() const noexcept { return decays.cbegin(); }
   iterator end() noexcept { return decays.end(); }
   const_iterator end() const noexcept { return decays.end(); }
   const_iterator cend() const noexcept { return decays.end(); }

   std::size_t size() const noexcept { return decays.size(); }

   void clear();
   void set_decay(double width, std::initializer_list<int> products);
   int get_particle_id() const { return initial_pdg; }
   const Decay& get_decay(std::initializer_list<int> products) const;
   double get_total_width() const { return total_width; }

private:
   int initial_pdg{0};
   std::unordered_map<std::size_t, Decay> decays{};
   double total_width{0.};
};

} // namespace flexiblesusy

#endif
