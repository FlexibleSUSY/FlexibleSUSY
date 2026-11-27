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
#include "always_false.hpp"
#include "cxx_qft/fields.hpp"
#include "cxx_qft/vertices.hpp"

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

/// CP property of a neutral Higgs boson
enum class CP_property {
   even,     ///< CP-even
   odd,      ///< CP-odd
   undefined ///< not a CP eigenstate (CP violation)
};

class NeutralHiggsEffectiveCouplings {
public:
   std::string particle {};
   int pdgid {};
   double mass {};
   double width {};
   double width_sm {};
   double invWidth {0.};
   CP_property CP{CP_property::undefined};
   std::pair<std::string, std::complex<double>> dd {};
   std::pair<std::string, std::complex<double>> uu {};
   std::pair<std::string, std::complex<double>> ss {};
   std::pair<std::string, std::complex<double>> cc {};
   std::pair<std::string, std::complex<double>> bb {};
   std::pair<std::string, std::complex<double>> tt {};
   std::pair<std::string, std::complex<double>> ee {};
   std::pair<std::string, std::complex<double>> mumu {};
   std::pair<std::string, std::complex<double>> tautau {};
   // sums of the partial widths of both charge combinations, in GeV
   double emu_width {0.};
   double etau_width {0.};
   double mutau_width {0.};
   std::pair<std::string, double> WW {};
   std::pair<std::string, double> ZZ {};
   std::pair<std::string, double> Zgam {};
   std::pair<std::string, double> gamgam {};
   std::pair<std::string, double> gg = {};
   double lam {};
   double undetectedWidth {0.};

   double get_undetected_width() const { return undetectedWidth; }
   /// adds the total width minus the invisible width to the
   /// undetected width (from which the widths of the channels passed
   /// to HiggsTools and Lilith were subtracted); a negative result is
   /// set to 0, with a warning if it is larger than rounding errors
   void calculate_undetected_br();
};

class EffectiveCoupling_list {
public:
   EffectiveCoupling_list() = default;
   ~EffectiveCoupling_list() = default;

   std::vector<NeutralHiggsEffectiveCouplings>::iterator begin() noexcept { return effective_coupling_list.begin(); }
   std::vector<NeutralHiggsEffectiveCouplings>::const_iterator begin() const noexcept { return effective_coupling_list.begin(); }
   std::vector<NeutralHiggsEffectiveCouplings>::const_iterator cbegin() const noexcept { return effective_coupling_list.cbegin(); }
   std::vector<NeutralHiggsEffectiveCouplings>::iterator end() noexcept { return effective_coupling_list.end(); }
   std::vector<NeutralHiggsEffectiveCouplings>::const_iterator end() const noexcept { return effective_coupling_list.end(); }
   std::vector<NeutralHiggsEffectiveCouplings>::const_iterator cend() const noexcept { return effective_coupling_list.end(); }

   NeutralHiggsEffectiveCouplings const& operator[](int index) const {
      return effective_coupling_list[index];
   }

   /// adds the effective coupling and the partial width of a decay of
   /// particle p into the final state fs (PDG codes); for the lepton
   /// flavour violating final states (e mu, e tau, mu tau) only the
   /// partial width is used, as HiggsTools takes widths for them
   void add_coupling(std::string const& p, std::array<int, 2> const& fs,
                     std::pair<std::string, double> const&, double partialW);
   void add_coupling(std::string const& p, std::array<int, 2> const& fs,
                     std::pair<std::string, std::complex<double>> const&,
                     double partialW);
   void set_invisible_width(std::string const& p, double);
   void push_back(NeutralHiggsEffectiveCouplings&& el) { effective_coupling_list.push_back(el); };

   std::size_t size() const noexcept { return effective_coupling_list.size(); }

private:
   std::vector<NeutralHiggsEffectiveCouplings> effective_coupling_list {};
};

std::size_t hash_decay(const Decay& decay);

/**
 * Whether a channel has a negative partial width that is more than
 * rounding noise, i.e. whose branching ratio (relative to the total
 * width, which is the sum of the absolute partial widths) is below
 * -1e-11.  Such channels are unphysical and are left out of the
 * output.
 */
bool has_unphysical_negative_width(const Decay&, double total_width);

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
         decays.insert(pos, std::make_pair(decay_hash, std::move(decay)));
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

   /**
    * Sets the partial widths of all channels with an unphysical
    * negative width (see has_unphysical_negative_width()) to zero and
    * removes them from the total width.
    *
    * @return a warning for each such channel
    */
   std::vector<std::string> remove_unphysical_negative_widths();

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
   if constexpr (cxx_diagrams::fields::is_singlet<FieldIn>) {
      // 1 -> 1, 1
      if constexpr (cxx_diagrams::fields::is_singlet<FieldOut1> && cxx_diagrams::fields::is_singlet<FieldOut2>) {
         return 1.;
      }
      // 1 -> 3, -3
      else if constexpr (
         (cxx_diagrams::fields::is_triplet<FieldOut1> && cxx_diagrams::fields::is_anti_triplet<FieldOut2>)
         || (cxx_diagrams::fields::is_anti_triplet<FieldOut1> && cxx_diagrams::fields::is_triplet<FieldOut2>)
      ) {
         return 3.;
      }
      // 1 -> 8, 8
      else if constexpr (cxx_diagrams::fields::is_octet<FieldOut1> && cxx_diagrams::fields::is_octet<FieldOut2>) {
         return 8.;
      }
      else {
         static_assert(always_false<FieldIn, FieldOut1, FieldOut2>, "Unknow colour structure in decay");
      }
   }
   else if constexpr (cxx_diagrams::fields::is_triplet<FieldIn>) {
      // 3 -> 1, 3
      if constexpr (
         (cxx_diagrams::fields::is_triplet<FieldOut1> && cxx_diagrams::fields::is_singlet<FieldOut2>)
         || (cxx_diagrams::fields::is_singlet<FieldOut1> && cxx_diagrams::fields::is_triplet<FieldOut2>)
      ) {
         return 1.;
      }
      // 3 -> 3, 8
      else if constexpr (
         (cxx_diagrams::fields::is_triplet<FieldIn> && cxx_diagrams::fields::is_octet<FieldOut2>)
         || (cxx_diagrams::fields::is_octet<FieldOut1> && cxx_diagrams::fields::is_triplet<FieldOut2>)
      ) {
         return 16/3.;
      }
      else {
         static_assert(always_false<FieldIn, FieldOut1, FieldOut2>, "Unknow colour structure in decay");
      }
   }
   else if constexpr (cxx_diagrams::fields::is_anti_triplet<FieldIn>) {
      // -3 -> 1, -3
      if constexpr (
         (cxx_diagrams::fields::is_anti_triplet<FieldOut1> && cxx_diagrams::fields::is_singlet<FieldOut2>)
         || (cxx_diagrams::fields::is_singlet<FieldOut1> && cxx_diagrams::fields::is_anti_triplet<FieldOut2>)
      ) {
         return 1.;
      }
      // -3 -> -3, 8
      else if constexpr (
         (cxx_diagrams::fields::is_anti_triplet<FieldIn> && cxx_diagrams::fields::is_octet<FieldOut2>)
         || (cxx_diagrams::fields::is_octet<FieldOut1> && cxx_diagrams::fields::is_anti_triplet<FieldOut2>)
      ) {
         return 16/3.;
      }
      else {
         static_assert(always_false<FieldIn, FieldOut1, FieldOut2>, "Unknow colour structure in decay");
      }
   }
   else if constexpr (cxx_diagrams::fields::is_octet<FieldIn>) {
      // 8 -> 1, 8
      if constexpr (
         (cxx_diagrams::fields::is_octet<FieldOut1> && cxx_diagrams::fields::is_singlet<FieldOut2>)
         || (cxx_diagrams::fields::is_singlet<FieldOut1> && cxx_diagrams::fields::is_octet<FieldOut2>)
      ) {
         return 1.;
      }
      else if constexpr (
         (cxx_diagrams::fields::is_triplet<FieldOut1> && cxx_diagrams::fields::is_anti_triplet<FieldOut2>)
         || (cxx_diagrams::fields::is_anti_triplet<FieldOut1> && cxx_diagrams::fields::is_triplet<FieldOut2>)
      ) {
         return 2.;
      }
      // 8 -> 8, 8 with identical particles in the final state
      // because of symmetry of the final state it must be proportional to d^2
      else if constexpr (
         cxx_diagrams::fields::is_octet<FieldOut1> && cxx_diagrams::fields::is_octet<FieldOut2> && std::is_same_v<FieldOut1, FieldOut2>
      ) {
         // color:   d^2 = (2 (4 - 5 Nc^2 + Nc^4) TR)/Nc = 40/3
         // average: 1/8
         return 40/24.;
      }
      // 8 -> 8, 8 with differnt particles in the final state
      else if constexpr (
         cxx_diagrams::fields::is_octet<FieldOut1> && cxx_diagrams::fields::is_octet<FieldOut2> && !std::is_same_v<FieldOut1, FieldOut2>
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
      if (idx1 == idx2) {
         return 0.5;
      }
      else {
         return 1.;
      }
   }
}

// utility functions for H->V*V*
struct HVV_4body_result {
   double value {};
   double error {};
};
HVV_4body_result hVV_4body(double mH, double mV, double GammaV);

/**
 * Width of H -> VV (V = W, Z) for x = mV^2/mH^2, combined from the double
 * off-shell (4-body), single off-shell (3-body) and on-shell (2-body) results
 * according to FlexibleDecay[4] (offshell_VV_decays).
 *
 * The results are switched where they are equal, so that the width is
 * continuous at the thresholds:
 *
 * - The on-shell width vanishes at the VV threshold, the double off-shell
 *   width doesn't.  Above the threshold, the larger of the two is used
 *   (offshell_VV_decays != 0).  They cross at x ~ 0.246 for W and Z, i.e.
 *   1.2 (W) and 1.4 (Z) GeV above the threshold, and the on-shell width is
 *   the larger one from there on, so the double off-shell width is only
 *   calculated for x > 0.24.
 * - The single off-shell width (offshell_VV_decays = 1) treats one vector
 *   boson as on-shell with zero width, so it diverges at the VV threshold.
 *   For x < 0.3 the smaller of the single and double off-shell widths is
 *   used.  They cross at x ~ 0.272 for W and Z.
 *
 * Below the single vector threshold (x >= 1) the double off-shell width is
 * used if offshell_VV_decays != 0.  With offshell_VV_decays = 1 the width
 * jumps there from the double off-shell width to the single off-shell one,
 * which vanishes at x = 1.  This step is small: 1.3e-6 GeV (W) and 1.0e-6
 * GeV (Z) for an SM-like coupling.
 *
 * @param x mV^2/mH^2
 * @param offshell_VV_decays value of FlexibleDecay[4]
 * @param four_body function returning the double off-shell width and the
 *    absolute error of its integration
 * @param three_body function returning the single off-shell width
 * @param two_body function returning the on-shell width
 *
 * @return width, and the absolute integration error if the double off-shell
 *    width is used (0 otherwise)
 */
template <typename FourBody, typename ThreeBody, typename TwoBody>
HVV_4body_result hVV_width(double x, int offshell_VV_decays,
                           const FourBody& four_body,
                           const ThreeBody& three_body, const TwoBody& two_body)
{
   static constexpr double x_min_four_body = 0.24;
   static constexpr double x_max_three_body = 0.3;

   if (4. * x < 1.) {
      const double on_shell = two_body();
      if (offshell_VV_decays != 0 && x > x_min_four_body) {
         const HVV_4body_result double_off_shell = four_body();
         if (double_off_shell.value > on_shell) {
            return double_off_shell;
         }
      }
      return {on_shell, 0.};
   }

   if (offshell_VV_decays == 0) {
      return {};
   }

   if (offshell_VV_decays == 2 || x >= 1.) {
      return four_body();
   }

   const double single_off_shell = three_body();
   if (x < x_max_three_body) {
      const HVV_4body_result double_off_shell = four_body();
      // the negation also replaces a NaN at x = 1/4
      if (!(single_off_shell <= double_off_shell.value)) {
         return double_off_shell;
      }
   }
   return {single_off_shell, 0.};
}

enum PDG_id_pairs: int {
   dd = 0,
   uu = 1,
   ss = 2,
   cc = 3,
   bb = 4,
   tt = 5,
   ee = 6,
   mumu = 7,
   tautau = 8,
   WW = 9,
   ZZ = 10,
   gamgam = 11,
   Zgam = 12,
   gg = 13,
   emu = 14,
   etau = 15,
   mutau = 16,
   NUMBER_OF_PDG_IDS
};

} // namespace flexiblesusy

#endif
