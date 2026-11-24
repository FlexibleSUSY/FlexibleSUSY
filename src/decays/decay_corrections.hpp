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
 * @file decay_corrections.hpp
 *
 * @brief contains struct for selection of decay correction
 * loop orders
 */

#ifndef DECAY_CORRECTIONS_H
#define DECAY_CORRECTIONS_H

#include <cstdint>

namespace flexiblesusy {

struct Decay_corrections {
   using Flags_t = std::int_fast32_t;

   enum class Positions : int {
      bsm, sm
   };

   Decay_corrections() = default;
   Decay_corrections(Flags_t);

   int sm{4};
   int bsm{1};

   void set(Flags_t);   ///< sets all values to digits in given a Flags_t
   Flags_t get() const; ///< returns all value in a Flags_t
};

} // namespace flexiblesusy

#endif
