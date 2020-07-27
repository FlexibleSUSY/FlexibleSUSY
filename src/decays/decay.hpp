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
#include <map>
#include <vector>
#include <string>

#include <boost/core/demangle.hpp>

namespace flexiblesusy {

class Decay {
public:
   Decay(int, std::initializer_list<int>, double, std::string const&);
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
   std::string get_proc_string() const { return proc_string; }
   void set_width(double w) { width = w; }

private:
   int pid_in{0};
   std::vector<int> pids_out{};
   double width{0.};
   std::string proc_string;
};

std::size_t hash_decay(const Decay& decay);

class Decays_list {
private:
   /* map is slower than unordered_map but will preserve order of etries */
   using List_type = std::map<std::size_t, Decay>;
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
   // @todo
   // the string cannot be taken by const&
   // if you do this, the printed to the screen total_width of Higgs boson
   // in the SM is shifted from 4.07 MeV 3.98 MeV
   // none of this makes any sense
   // this happend with clang 10.0.0, but not with gcc 9.3.1
   void set_decay(double width, std::initializer_list<int> products, std::string);
   int get_particle_id() const { return initial_pdg; }
   const Decay& get_decay(std::initializer_list<int> products) const;
   double get_total_width() const { return total_width; }

private:
   int initial_pdg{0};
   List_type decays{};
   double total_width{0.};
};

std::string strip_field_namespace(std::string const&);

template<typename FieldIn, typename FieldOut1, typename FieldOut2>
std::string create_process_string(
      std::array<int, FieldIn::numberOfFieldIndices> const in,
      std::array<int, FieldOut1::numberOfFieldIndices> const out1,
      std::array<int, FieldOut2::numberOfFieldIndices> const out2) {

   auto vector_to_idx = [](auto v) {
      if (v.empty()) {
         return std::string();
      }
      else {
         return "(" + std::to_string(v[0]) + ")";
      }
   };

   using boost::core::demangle;
   std::string process_string =
         strip_field_namespace(demangle(typeid(FieldIn).name())) + vector_to_idx(in)
         + "->{" +
         strip_field_namespace(demangle(typeid(FieldOut1).name())) + vector_to_idx(out1) + "," +
         strip_field_namespace(demangle(typeid(FieldOut2).name())) + vector_to_idx(out2) +
         "}";

   return process_string;
}

} // namespace flexiblesusy

#endif
