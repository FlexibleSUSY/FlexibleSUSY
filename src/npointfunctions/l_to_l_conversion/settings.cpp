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

#include "settings.hpp"
#include "error.hpp"
#include "string_format.hpp"

#include <cmath>
#include <iostream>
#include <string>

namespace flexiblesusy {

namespace {
const std::array<std::string, LToLConversion_settings::NUMBER_OF_OPTIONS> descriptions = {
   "include tensor contribution",
   "include gluonic contribution",
};

bool is_integer(double value)
{
   double intpart;
   return std::modf(value, &intpart) == 0.0;
}

void assert_bool(double value, const char* quantity)
{
   if (value != 0.0 && value != 1.0) {
      throw SetupError(std::string(quantity) + " must either 0 or 1");
   }
}

void assert_integer(double value, const char* quantity)
{
   if (!is_integer(value)) {
      throw SetupError(std::string(quantity) + " must be an integer");
   }
}

void assert_ge(double value, double lower_bound, const char* quantity)
{
   if (value < lower_bound) {
      throw SetupError(std::string(quantity) +
                       " must be greater than or equal to " +
                       flexiblesusy::to_string(lower_bound));
   }
}

void assert_le(double value, double upper_bound, const char* quantity)
{
   if (value > upper_bound) {
      throw SetupError(std::string(quantity) +
                       " must be lower than or equal to " +
                       flexiblesusy::to_string(upper_bound));
   }
}

} // anonymous namespace

/**
 * Default constructor
 *
 * Calls reset() to initialize all spectrum generator settings to
 * their default values.
 */
LToLConversion_settings::LToLConversion_settings()
{
   reset();
}

double LToLConversion_settings::get(Settings o) const
{
   return values.at(o);
}

LToLConversion_settings::Settings_t LToLConversion_settings::get() const
{
   Settings_t s(&values[0]);
   return s;
}

std::string LToLConversion_settings::get_description(Settings o) const
{
   return descriptions.at(o);
}

void LToLConversion_settings::set(Settings o, double value)
{
   switch (o) {
   case include_tensor_contribution: // 0 [bool]
      assert_bool(value, descriptions.at(o).c_str());
      break;
   case include_gluonic_contribution: // 1 [bool]
      assert_bool(value, descriptions.at(o).c_str());
      break;
   /*case min_br_to_print: // 1 [double >= 0 and <= 1]
      assert_ge(value, 0, descriptions.at(o).c_str());
      assert_le(value, 1, descriptions.at(o).c_str());
      break;*/
   default:
      break;
   }

   values.at(o) = value;
}

void LToLConversion_settings::set(const LToLConversion_settings::Settings_t& s)
{
   std::copy(s.data(), s.data() + s.size(), values.begin());
}

/**
 * Resets all spectrum generator settings to their defaults.
 *
 * | enum                             | possible values   | default value   |
 * |----------------------------------|-------------------|-----------------|
 * | include_tensor_contribution      | 0 (no) or 1 (yes) | 1 (= enabled)   |
 * | include_gluonic_contribution     | 0 (no) or 1 (yes) | 1 (= enabled)   |
 */
void LToLConversion_settings::reset()
{
   values[include_tensor_contribution]  = 1.0;
   values[include_gluonic_contribution] = 1.0;
}
bool is_integer(double value)
{
   double intpart;
   return std::modf(value, &intpart) == 0.0;
}

}
