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
 * @file threshold_corrections.hpp
 *
 * @brief contains struct for selection of low-energy threshold
 * correction loop orders
 */

#ifndef THRESHOLD_CORRECTIONS_H
#define THRESHOLD_CORRECTIONS_H

#include <cstdint>
#include <iosfwd>

namespace flexiblesusy {

struct Threshold_corrections {
   using Flags_t = std::int_fast32_t;

   enum class Positions : int {
      alpha_em, sin_theta_w, alpha_s, mz, mw, mh, mt, mb, mtau
   };

   Threshold_corrections() = default;
   Threshold_corrections(Flags_t);

   int alpha_em{1};
   int sin_theta_w{2};
   int alpha_s{4};
   int mz{1};
   int mw{1};
   int mh{1};
   int mt{4};
   int mb{2};
   int mtau{1};

   void set(Flags_t);   ///< sets all values to digits in given a Flags_t
   Flags_t get() const; ///< returns all value in a Flags_t
};

std::ostream& operator<<(std::ostream&, const Threshold_corrections&);

} // namespace flexiblesusy

#endif
