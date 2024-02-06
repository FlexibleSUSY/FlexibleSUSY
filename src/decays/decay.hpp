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
#include <boost/range/algorithm/equal.hpp>

#include "always_false.hpp"
#include "cxx_qft/fields.hpp"

namespace flexiblesusy {

class Decay {
public:

template <typename T>
   Decay(
      int pid_in_, std::initializer_list<int> pids_out_, double width_, T&& proc_string_)
      : pid_in{pid_in_}
      , pids_out{pids_out_}
      , width{width_}
      , proc_string{std::forward<T>(proc_string_)}
   {
      std::sort(pids_out.begin(), pids_out.end());
   }
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
   /* map is slower than unordered_map but will preserve order of entries */
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

   template <typename T>
   void set_decay(double width, std::initializer_list<int> pids_out, T&& proc_string)
   {
      const Decay decay(initial_pdg, pids_out, width, std::forward<T>(proc_string));
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
      // as the sum of |width|
      total_width += std::abs(width);
   }

   int get_particle_id() const { return initial_pdg; }
   const Decay& get_decay(std::initializer_list<int> products) const;
   double get_total_width() const { return total_width; }

private:
   int initial_pdg{0};
   List_type decays{};
   double total_width{0.};
};

/// sort decays w.r.t. their width
std::vector<Decay> sort_decays_list(const Decays_list&);

std::string strip_field_namespace(std::string const&);

template<typename Field>
std::string field_as_string(std::array<int, Field::numberOfFieldIndices> const& idx) {
   const std::string field = strip_field_namespace(boost::core::demangle(typeid(Field).name()));
   if constexpr (Field::numberOfFieldIndices == 0) {
      return field;
   }
   else if constexpr (Field::numberOfFieldIndices == 1) {
      // in the output we count particles from 1 (not 0)
      return field + "(" + std::to_string(idx[0]+1) + ")";
   }
   else {
      static_assert(always_false<Field>, "Field is expected to have 0 or 1 index");
   }
}

template<typename FieldIn, typename FieldOut1, typename FieldOut2>
std::string create_process_string(
      std::array<int, FieldIn::numberOfFieldIndices> const in,
      std::array<int, FieldOut1::numberOfFieldIndices> const out1,
      std::array<int, FieldOut2::numberOfFieldIndices> const out2) {

   std::string process_string =
      field_as_string<FieldIn>(in) + " -> " +
      field_as_string<FieldOut1>(out1) + " " + field_as_string<FieldOut2>(out2);

   return process_string;
}

// returns a squared color generator for a 3 point amplitude with FieldIn, FieldOut1 and FieldOut2
// averaged over inital state colors
// the generator is guessed from color representations of FieldIn, FieldOut1 and FieldOut2
// This is not a bulletproof solution and might fail in general but is enough for
// decays of color singlets

template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr double squared_color_generator() noexcept {
   if constexpr (cxx_diagrams::fields::is_singlet_v<FieldIn>) {
      // 1 -> 1, 1
      if constexpr (cxx_diagrams::fields::is_singlet_v<FieldOut1> && cxx_diagrams::fields::is_singlet_v<FieldOut2>) {
         return 1.;
      }
      // 1 -> 3, -3
      else if constexpr (
         (cxx_diagrams::fields::is_triplet_v<FieldOut1> && cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>)
         || (cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> && cxx_diagrams::fields::is_triplet_v<FieldOut2>)
      ) {
         return 3.;
      }
      // 1 -> 8, 8
      else if constexpr (cxx_diagrams::fields::is_octet_v<FieldOut1> && cxx_diagrams::fields::is_octet_v<FieldOut2>) {
         return 8.;
      }
      else {
         static_assert(always_false<FieldIn, FieldOut1, FieldOut2>, "Unknow colour structure in decay");
      }
   }
   else if constexpr (cxx_diagrams::fields::is_triplet_v<FieldIn>) {
      // 3 -> 1, 3
      if constexpr (
         (cxx_diagrams::fields::is_triplet_v<FieldOut1> && cxx_diagrams::fields::is_singlet_v<FieldOut2>)
         || (cxx_diagrams::fields::is_singlet_v<FieldOut1> && cxx_diagrams::fields::is_triplet_v<FieldOut2>)
      ) {
         return 1.;
      }
      // 3 -> 3, 8
      else if constexpr (
         (cxx_diagrams::fields::is_triplet_v<FieldIn> && cxx_diagrams::fields::is_octet_v<FieldOut2>)
         || (cxx_diagrams::fields::is_octet_v<FieldOut1> && cxx_diagrams::fields::is_triplet_v<FieldOut2>)
      ) {
         return 4.;
      }
      else {
         static_assert(always_false<FieldIn, FieldOut1, FieldOut2>, "Unknow colour structure in decay");
      }
   }
   else if constexpr (cxx_diagrams::fields::is_anti_triplet_v<FieldIn>) {
      // -3 -> 1, -3
      if constexpr (
         (cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> && cxx_diagrams::fields::is_singlet_v<FieldOut2>)
         || (cxx_diagrams::fields::is_singlet_v<FieldOut1> && cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>)
      ) {
         return 1.;
      }
      // -3 -> -3, 8
      else if constexpr (
         (cxx_diagrams::fields::is_anti_triplet_v<FieldIn> && cxx_diagrams::fields::is_octet_v<FieldOut2>)
         || (cxx_diagrams::fields::is_octet_v<FieldOut1> && cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>)
      ) {
         return 4.;
      }
      else {
         static_assert(always_false<FieldIn, FieldOut1, FieldOut2>, "Unknow colour structure in decay");
      }
   }
   else if constexpr (cxx_diagrams::fields::is_octet_v<FieldIn>) {
      // 8 -> 1, 8
      if constexpr (
         (cxx_diagrams::fields::is_octet_v<FieldOut1> && cxx_diagrams::fields::is_singlet_v<FieldOut2>)
         || (cxx_diagrams::fields::is_singlet_v<FieldOut1> && cxx_diagrams::fields::is_octet_v<FieldOut2>)
      ) {
         return 1.;
      }
      else if constexpr (
         (cxx_diagrams::fields::is_triplet_v<FieldOut1> && cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>)
         || (cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> && cxx_diagrams::fields::is_triplet_v<FieldOut2>)
      ) {
         return 1./2.;
      }
      // 8 -> 8, 8 with identical particles in the final state
      // because of symmetry of the final state it must be proportional to d^2
      else if constexpr (
         cxx_diagrams::fields::is_octet_v<FieldOut1> && cxx_diagrams::fields::is_octet_v<FieldOut2> && std::is_same_v<FieldOut1, FieldOut2>
      ) {
         // color:   d^2 = (2 (4 - 5 Nc^2 + Nc^4) TR)/Nc = 40/3
         // average: 1/8
         return 40/24.;
      }
      // 8 -> 8, 8 with differnt particles in the final state
      else if constexpr (
         cxx_diagrams::fields::is_octet_v<FieldOut1> && cxx_diagrams::fields::is_octet_v<FieldOut2> && !std::is_same_v<FieldOut1, FieldOut2>
      ) {
         // color:   f^2 = 2 Nc (-1 + Nc^2) TR = 24
         // average: 1/8
         return 3.;
      }
      else {
         static_assert(always_false<FieldIn, FieldOut1, FieldOut2>, "Unknow colour structure in decay");
      }
   }
}

template <typename Field1, typename Field2>
double final_state_symmetry_factor(typename cxx_diagrams::field_indices<Field1>::type const& idx1,
                            typename cxx_diagrams::field_indices<Field2>::type const& idx2)
{
   if constexpr (!std::is_same_v<Field1, Field2>) {
      return 1.;
   }
   else {
      if (boost::range::equal(idx1, idx2)) {
         return 0.5;
      }
      else {
         return 1.;
      }
   }
}

// utility functions for H->V*V*
double hVV_4body(double *q2, size_t dim, void *params);
struct hVV_4body_params {
   double mHOS {};
   double mVOS {};
   double GammaV {};
};

} // namespace flexiblesusy

#endif
