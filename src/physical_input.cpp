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

#include "physical_input.hpp"
#include "error.hpp"
#include <algorithm>

namespace flexiblesusy {

/**
 * Default constructor
 *
 * Calls reset() to initialize all physical input parameters to their
 * default values.
 */
Physical_input::Physical_input()
{
   reset();
}

double Physical_input::get(Input o) const
{
   return values.at(o);
}

Eigen::ArrayXd Physical_input::get() const
{
   Eigen::ArrayXd vec(values.size());

   std::copy(values.cbegin(), values.cend(), vec.data());

   return vec;
}

const std::array<std::string, Physical_input::NUMBER_OF_INPUT_PARAMETERS>& Physical_input::get_names()
{
   static const std::array<std::string, NUMBER_OF_INPUT_PARAMETERS> names = {
      "alpha_em(0)",
      "mh_pole"
   };
   return names;
}

void Physical_input::set(Input o, double value)
{
   values.at(o) = value;
}

void Physical_input::set(const Eigen::ArrayXd& vec)
{
   if (vec.size() != static_cast<decltype(vec.size())>(values.size()))
      throw SetupError("Parameters array has wrong size");

   std::copy(vec.data(), vec.data() + vec.size(), values.begin());
}

/**
 * Resets all physical input parameters to their default values.
 *
 * | enum                             | possible values              | default value   |
 * |----------------------------------|------------------------------|-----------------|
 * | alpha_em_0                       | any positive double          | 1/137.035999074 |
 * | mh_pole                          | any positive double          | 125.09          |
 */
void Physical_input::reset()
{
   values[alpha_em_0] = 1./137.035999074;
   values[mh_pole] = 125.09;
}

} // namespace flexiblesusy
