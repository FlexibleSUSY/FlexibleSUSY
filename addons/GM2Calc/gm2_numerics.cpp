// ====================================================================
// This file is part of GM2Calc.
//
// GM2Calc is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// GM2Calc is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GM2Calc.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#include "gm2_numerics.hpp"

namespace gm2calc {

double abs_sqrt(double x) {
   return std::sqrt(std::abs(x));
}

int sign(double x) noexcept { return x < 0 ? -1 : 1; }

double signed_sqr(double x) { return sign(x) * x * x; }

double signed_abs_sqrt(double x) {
   return sign(x) * std::sqrt(std::abs(x));
}

} // namespace gm2calc
