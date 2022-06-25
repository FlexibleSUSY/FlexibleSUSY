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

#include <algorithm>
#include <initializer_list>
#include <map>
#include <string>
#include <vector>

#include <boost/core/demangle.hpp>

#include "wrappers.hpp"
#include "cxx_qft/fields.hpp"
#include "cxx_qft/vertices.hpp"

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

struct NeutralHiggsEffectiveCouplings {
   std::string particle = "";
   double mass = 0.;
   double width = 0.;
   std::complex<double> dd = 0.;
   std::complex<double> uu = 0.;
   std::complex<double> ss = 0.;
   std::complex<double> cc = 0.;
   std::complex<double> bb = 0.;
   std::complex<double> tt = 0.;
   std::complex<double> ee = 0.;
   std::complex<double> mumu = 0.;
   std::complex<double> tautau = 0.;
   double WW = 0.;
   double ZZ = 0.;
   double Zgam = 0.;
   double gamgam = 0.;
   double gg = 0.;
};

class EffectiveCoupling_list {
public:
   EffectiveCoupling_list() = default;
   ~EffectiveCoupling_list() = default;
   EffectiveCoupling_list(const EffectiveCoupling_list&) = default;
   EffectiveCoupling_list(EffectiveCoupling_list&&) = default;
   EffectiveCoupling_list& operator=(const EffectiveCoupling_list&) = default;
   EffectiveCoupling_list& operator=(EffectiveCoupling_list&&) = default;

   std::vector<NeutralHiggsEffectiveCouplings>::iterator begin() noexcept { return effective_coupling_list.begin(); }
   std::vector<NeutralHiggsEffectiveCouplings>::const_iterator begin() const noexcept { return effective_coupling_list.begin(); }
   std::vector<NeutralHiggsEffectiveCouplings>::const_iterator cbegin() const noexcept { return effective_coupling_list.cbegin(); }
   std::vector<NeutralHiggsEffectiveCouplings>::iterator end() noexcept { return effective_coupling_list.end(); }
   std::vector<NeutralHiggsEffectiveCouplings>::const_iterator end() const noexcept { return effective_coupling_list.end(); }
   std::vector<NeutralHiggsEffectiveCouplings>::const_iterator cend() const noexcept { return effective_coupling_list.end(); }
   NeutralHiggsEffectiveCouplings const& operator[](int index) const {
      return effective_coupling_list[index];
   }

   void add_coupling(std::string p, std::array<int, 2> fs, double c) {
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
         effective_coupling_list.push_back(effC);
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
   void add_coupling(std::string p, std::array<int, 2> fs, std::complex<double> c) {
      auto found = std::find_if(
         std::begin(effective_coupling_list), std::end(effective_coupling_list),
         [&p](NeutralHiggsEffectiveCouplings const& effC) {return effC.particle == p;}
      );
      auto are_the_same = [](std::array<int, 2> const& a1, std::array<int, 2> const& a2) {return std::is_permutation(a1.begin(), a1.end(), a2.begin(), a2.end());};
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
         effective_coupling_list.push_back(effC);
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
      }
   }
   std::size_t size() const noexcept { return effective_coupling_list.size(); }
private:
   std::vector<NeutralHiggsEffectiveCouplings> effective_coupling_list {};
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
   void set_decay(double width, std::initializer_list<int> products, std::string const&);
   int get_particle_id() const { return initial_pdg; }
   const Decay& get_decay(std::initializer_list<int> products) const;
   double get_total_width() const { return total_width; }

private:
   int initial_pdg{0};
   List_type decays{};
   double total_width{0.};
};

std::string strip_field_namespace(std::string const&);

template<typename Field>
std::string field_as_string(std::array<int, Field::numberOfFieldIndices> const& idx) {
   auto vector_to_idx = [](auto v) {
      if (v.empty()) {
         return std::string();
      }
      else {
         return "(" + std::to_string(v[0]) + ")";
      }
   };

   using boost::core::demangle;
   return strip_field_namespace(demangle(typeid(Field).name())) + vector_to_idx(idx);
}

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

// returns a squared color generator for a 3 point amplitude with FieldIn, FieldOut1 and FieldOut2
// averaged over inital state colors
// the generator is guessed from color representations of FieldIn, FieldOut1 and FieldOut2
// This is not a bulletproof solution and might fail in general but is enough for
// decays of color singlets

// 1 -> 1, 1
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
cxx_diagrams::fields::is_singlet_v<FieldIn> &&
cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
cxx_diagrams::fields::is_singlet_v<FieldOut2>, double>
squared_color_generator() {return 1.;}

// 1 -> 3, 3bar
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
   cxx_diagrams::fields::is_singlet_v<FieldIn>
   &&
   (
   (cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
   cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>)
   ||
   (cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
   cxx_diagrams::fields::is_triplet_v<FieldOut2>)
   ), double
>
squared_color_generator() {return 3.;}

// 1 -> 8, 8
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
cxx_diagrams::fields::is_singlet_v<FieldIn> &&
cxx_diagrams::fields::is_octet_v<FieldOut1> &&
cxx_diagrams::fields::is_octet_v<FieldOut2>, double>
squared_color_generator() {return 8.;}

// 3 -> 3, 1; 3bar -> 3bar, 1
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
(
cxx_diagrams::fields::is_triplet_v<FieldIn> &&
((cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
cxx_diagrams::fields::is_triplet_v<FieldOut2>))
)
||
(
cxx_diagrams::fields::is_anti_triplet_v<FieldIn> &&
((cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>))
)
, double>
squared_color_generator() {return 1.;}

// 3 -> 3, 8; 3bar -> 3bar, 8
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
(
cxx_diagrams::fields::is_triplet_v<FieldIn> &&
((cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
cxx_diagrams::fields::is_octet_v<FieldOut2>) ||
(cxx_diagrams::fields::is_octet_v<FieldOut1> &&
cxx_diagrams::fields::is_triplet_v<FieldOut2>))
)
||
(
cxx_diagrams::fields::is_anti_triplet_v<FieldIn> &&
((cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
cxx_diagrams::fields::is_octet_v<FieldOut2>) ||
(cxx_diagrams::fields::is_octet_v<FieldOut1> &&
cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>))
)
, double>
squared_color_generator() {return 4.;}

// 8 -> 8, 1
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
cxx_diagrams::fields::is_octet_v<FieldIn> &&
((cxx_diagrams::fields::is_octet_v<FieldOut1> &&
cxx_diagrams::fields::is_singlet_v<FieldOut2>) ||
(cxx_diagrams::fields::is_singlet_v<FieldOut1> &&
cxx_diagrams::fields::is_octet_v<FieldOut2>))
, double>
squared_color_generator() {return 1.;}

// 8 -> 3, 3bar
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
   cxx_diagrams::fields::is_octet_v<FieldIn>
   &&
   (
   (cxx_diagrams::fields::is_triplet_v<FieldOut1> &&
   cxx_diagrams::fields::is_anti_triplet_v<FieldOut2>)
   ||
   (cxx_diagrams::fields::is_anti_triplet_v<FieldOut1> &&
   cxx_diagrams::fields::is_triplet_v<FieldOut2>)
   ), double
>
squared_color_generator() {return 1./2.;}

// 8 -> 8, 8 with identical particles in the final state
// because of symmetry of the final state it must be proportional to d^2
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
cxx_diagrams::fields::is_octet_v<FieldIn> &&
cxx_diagrams::fields::is_octet_v<FieldOut1> &&
cxx_diagrams::fields::is_octet_v<FieldOut2> &&
std::is_same<FieldOut1, FieldOut2>::value
, double>
// color:   d^2 = (2 (4 - 5 Nc^2 + Nc^4) TR)/Nc = 40/3
// average: 1/8
squared_color_generator() {return 40/24.;}

// 8 -> 8, 8 with differnt particles in the final state
template<typename FieldIn, typename FieldOut1, typename FieldOut2>
constexpr
std::enable_if_t<
cxx_diagrams::fields::is_octet_v<FieldIn> &&
cxx_diagrams::fields::is_octet_v<FieldOut1> &&
cxx_diagrams::fields::is_octet_v<FieldOut2> &&
!std::is_same<FieldOut1, FieldOut2>::value
, double>
// color:   f^2 = 2 Nc (-1 + Nc^2) TR = 24
// average: 1/8
squared_color_generator() {return 3.;}

template <typename Field1, typename Field2>
constexpr std::enable_if_t<!std::is_same<Field1, Field2>::value, double>
final_state_symmetry_factor(typename cxx_diagrams::field_indices<Field1>::type const& idx1,
                            typename cxx_diagrams::field_indices<Field2>::type const& idx2)
{
   return 1.;
}

template <typename Field1, typename Field2>
std::enable_if_t<std::is_same<Field1, Field2>::value, double>
final_state_symmetry_factor(typename cxx_diagrams::field_indices<Field1>::type const& idx1,
                            typename cxx_diagrams::field_indices<Field2>::type const& idx2)
{
   if (std::equal(idx1.begin(), idx1.end(), idx2.begin(), idx2.end())) {
      return 0.5;
   }
   else {
      return 1.;
   }
}

struct my_f_params {double mHOS; double mVOS; double GammaV;};
double hVV_4body(double *q2, size_t dim, void *params);

} // namespace flexiblesusy

#endif
