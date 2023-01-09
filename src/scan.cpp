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

#include "scan.hpp"
#include "error.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>

/**
 * @file scan.cpp
 * @brief contains helper functions and classes for parameter scans
 */

namespace flexiblesusy {

/**
 * Returns vector with number_of_steps floating point values between
 * start and stop.  The endpoint is excluded.
 *
 * @param start smallest value
 * @param stop largest value (excluded)
 * @param number_of_steps number of values
 *
 * @return vector of floating point values
 */
std::vector<double> float_range(double start, double stop,
                                std::size_t number_of_steps)
{
   const double step_size = (stop - start) / number_of_steps;
   std::vector<double> result(number_of_steps);

   for (std::size_t i = 0; i < number_of_steps; ++i) {
      result[i] = start + i*step_size;
   }

   return result;
}

/**
 * Returns vector with number_of_steps floating point values between
 * start and stop.  The values are logarithmically distributed,
 * i.e. the difference of the logarithm of two consecutive points is
 * constant.  The endpoint is excluded.
 *
 * @param start smallest value
 * @param stop largest value (excluded)
 * @param number_of_steps number of values
 *
 * @return vector of floating point values
 */
std::vector<double> float_range_log(double start, double stop,
                                    std::size_t number_of_steps)
{
   if (start <= 0 || stop <= 0) {
      throw flexiblesusy::OutOfBoundsError("float_range_log: interval boundaries must be > 0.");
   }

   const double log_start = std::log(start);
   const double log_stop = std::log(stop);
   auto result = float_range(log_start, log_stop, number_of_steps);

   std::transform(result.begin(), result.end(), result.begin(),
                  [] (double x) { return std::exp(x); });

   return result;
}

/**
 * Returns vector with (number_of_divisions + 1) floating point values
 * between (including) start and (including) stop.
 *
 * @param start smallest value (included)
 * @param stop largest value (included)
 * @param number_of_divisions number of divisions
 *
 * @return vector of floating point values
 */
std::vector<double> subdivide(double start, double stop,
                              std::size_t number_of_divisions)
{
   if (number_of_divisions == 0) {
      throw flexiblesusy::OutOfBoundsError("subdivide: number_of_divisions must be > 0.");
   }

   auto result = float_range(start, stop, number_of_divisions);
   result.push_back(stop);

   return result;
}

/**
 * Returns vector with (number_of_divisions + 1) floating point values
 * between (including) start and (including) stop with logarithmic
 * spacing.
 *
 * @param start smallest value (included)
 * @param stop largest value (included)
 * @param number_of_divisions number of divisions
 *
 * @return vector of floating point values
 */
std::vector<double> subdivide_log(double start, double stop,
                                  std::size_t number_of_divisions)
{
   if (number_of_divisions == 0) {
      throw flexiblesusy::OutOfBoundsError("subdivide_log: number_of_divisions must be > 0.");
   }

   auto result = float_range_log(start, stop, number_of_divisions);
   result.push_back(stop);

   return result;
}

} // namespace flexiblesusy
