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

#ifndef FUNCTORS_H
#define FUNCTORS_H

#include <cmath>
#include <limits>

namespace flexiblesusy {

/**
 * @class Chop
 * @brief function object, whose operator() returns zero if the value
 * it is applied to is smaller than the given threshold
 * @tparam real scalar type
 */
template <typename RealScalar>
struct Chop {
   Chop() = default;
   Chop(RealScalar threshold_) : threshold(std::abs(threshold_)) {}
   const RealScalar operator()(const RealScalar& x) const noexcept {
      return std::abs(x) < threshold ? RealScalar{} : x;
   }
   RealScalar threshold{std::numeric_limits<RealScalar>::epsilon()};
};

} // namespace flexiblesusy

#endif
