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

#ifndef SCAN_HPP
#define SCAN_HPP

#include <cstddef>
#include <vector>

namespace flexiblesusy {

/// returns range of floating point values between start (included) and stop (excluded)
std::vector<double> float_range(double start, double stop,
                                std::size_t number_of_steps);

/// returns range of floating point values between start (included) and stop (excluded) with logarithmic spacing
std::vector<double> float_range_log(double start, double stop,
                                    std::size_t number_of_steps);

/// returns range of (number_of_divisions + 1) floating point values between (including) start and (including) stop
std::vector<double> subdivide(double start, double stop,
                              std::size_t number_of_divisions);

/// returns range of (number_of_divisions + 1) floating point values between (including) start and (including) stop with logarithmic spacing
std::vector<double> subdivide_log(double start, double stop,
                                  std::size_t number_of_divisions);

} // namespace flexiblesusy

#endif
