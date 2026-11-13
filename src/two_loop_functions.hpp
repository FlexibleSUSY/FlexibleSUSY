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

#ifndef TWO_LOOP_FUNCTIONS
#define TWO_LOOP_FUNCTIONS

#include "cextensions.hpp"

#include <array>
#include <cstddef>
#include <functional>
#include <unordered_map>

#ifndef PVATTR
#define PVATTR noexcept ATTR(const)
#endif

namespace flexiblesusy {
namespace two_loop_functions {

/// renormalization scheme of the 2-loop functions
enum class Scheme { DRbar, MSbar };

// tadpole functions

double TfSS(double, double, double, double) PVATTR;
double TfSSS(double, double, double, double) PVATTR;
double TfSSSS(double, double, double, double, double) PVATTR;
double TfSSFF(double, double, double, double, double) PVATTR;
double TfSSFbFb(double, double, double, double, double) PVATTR;
double TfFFFbS(double, double, double, double, double) PVATTR;
double TfFFbFS(double, double, double, double, double) PVATTR;
double TfFbFbFbS(double, double, double, double, double) PVATTR;
double TfSV(double, double) PVATTR;
double TfFV(double, double, Scheme = Scheme::DRbar) PVATTR;

// self-energy functions

double WfSSSS(double, double, double, double, double, double) PVATTR;
double XfSSS(double, double, double, double, double) PVATTR;
double YfSSSS(double, double, double, double, double, double) PVATTR;
double ZfSSSS(double, double, double, double, double, double) PVATTR;
double SfSSS(double, double, double, double, double) PVATTR;
double UfSSSS(double, double, double, double, double, double) PVATTR;
double VfSSSSS(double, double, double, double, double, double, double) PVATTR;
double MfSSSSS(double, double, double, double, double, double, double) PVATTR;
double WfSSSV(double, double, double) PVATTR;
double MfSSSSV(double, double, double, double) PVATTR;
double WfSSFF(double, double, double, double, double, double) PVATTR;
double WfSSFbFb(double, double, double, double, double, double) PVATTR;
double MfFbFbFbFbS(double, double, double, double, double, double, double) PVATTR;
double MfFFbFbFS(double, double, double, double, double, double, double) PVATTR;
double MfFFbFFbS(double, double, double, double, double, double, double) PVATTR;
double MfFFFbFbS(double, double, double, double, double, double, double) PVATTR;
double MfFFFFS(double, double, double, double, double, double, double) PVATTR;
double MfSFbSFbFb(double, double, double, double, double, double, double) PVATTR;
double MfSFSFbF(double, double, double, double, double, double, double) PVATTR;
double MfSFSFFb(double, double, double, double, double, double, double) PVATTR;
double VfSSSFbFb(double, double, double, double, double, double, double) PVATTR;
double VfSSSFF(double, double, double, double, double, double, double) PVATTR;
double VfFbFbFbFbS(double, double, double, double, double, double, double) PVATTR;
double VfFbFFbFS(double, double, double, double, double, double, double) PVATTR;
double VfFbFFFbS(double, double, double, double, double, double, double) PVATTR;
double VfFFbFbFS(double, double, double, double, double, double, double) PVATTR;
double VfFFFbFbS(double, double, double, double, double, double, double) PVATTR;
double VfFFFFS(double, double, double, double, double, double, double) PVATTR;
double GfFFV(double, double, double, double, Scheme = Scheme::DRbar) PVATTR;
double GfFbFbV(double, double, double, double, Scheme = Scheme::DRbar) PVATTR;

/**
 * Cache for the values of a pure function of N double arguments.
 *
 * The generated 2-loop self-energies call the same loop functions
 * with identical arguments many times (once for each external index
 * combination), so the loop function values are cached.  The key
 * must contain all arguments, including the renormalization scale.
 * The cache is cleared when it reaches max_size entries.
 */
template <std::size_t N>
class Memo {
public:
   using Key = std::array<double, N>;

   template <class F>
   double operator()(const Key& key, F f)
   {
      const auto it = cache.find(key);
      if (it != cache.end()) {
         return it->second;
      }
      if (cache.size() >= max_size) {
         cache.clear();
      }
      const double value = f();
      cache.emplace(key, value);
      return value;
   }

private:
   struct Hash {
      std::size_t operator()(const Key& key) const noexcept
      {
         std::size_t seed = 0;
         for (const double k : key) {
            seed ^= std::hash<double>{}(k) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
         }
         return seed;
      }
   };

   static constexpr std::size_t max_size = std::size_t(1) << 18;
   std::unordered_map<Key, double, Hash> cache;
};

} // namespace two_loop_functions
} // namespace flexiblesusy

#endif // TWO_LOOP_FUNCTIONS
