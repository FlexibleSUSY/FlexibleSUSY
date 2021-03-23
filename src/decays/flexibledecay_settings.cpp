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

#include "flexibledecay_settings.hpp"
#include "error.hpp"
#include "string_format.hpp"

#include <cmath>
#include <iostream>
#include <string>

namespace flexiblesusy {

namespace {
const std::array<std::string, FlexibleDecay_settings::NUMBER_OF_OPTIONS> descriptions = {
   "calculate particle decays",
   "higher order corrections in decays",
   "off-shell decays into VV pair",
   "min. br to print"
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
FlexibleDecay_settings::FlexibleDecay_settings()
{
   reset();
}

double FlexibleDecay_settings::get(Settings o) const
{
   return values.at(o);
}

FlexibleDecay_settings::Settings_t FlexibleDecay_settings::get() const
{
   Settings_t s(&values[0]);
   return s;
}

std::string FlexibleDecay_settings::get_description(Settings o) const
{
   return descriptions.at(o);
}

void FlexibleDecay_settings::set(Settings o, double value)
{
   switch (o) {
   case calculate_decays: // 0 [bool]
      assert_bool(value, descriptions.at(o).c_str());
      break;
   case include_higher_order_corrections: // 1 [int >= 0 and <= 1]
      assert_integer(value, descriptions.at(o).c_str());
      break;
   case offshell_VV_decays: // 2 [int >= 0 and <= 2]
      assert_integer(value, descriptions.at(o).c_str());
      assert_ge(value, 0, descriptions.at(o).c_str());
      assert_le(value, 2, descriptions.at(o).c_str());
      break;
   case min_br_to_print: // 2 [int >= 0 and <= 2]
      assert_ge(value, 0, descriptions.at(o).c_str());
      assert_le(value, 1, descriptions.at(o).c_str());
      break;
   default:
      break;
   }

   values.at(o) = value;
}

void FlexibleDecay_settings::set(const FlexibleDecay_settings::Settings_t& s)
{
   std::copy(s.data(), s.data() + s.size(), values.begin());
}

/**
 * Resets all spectrum generator settings to their defaults.
 *
 * | enum                             | possible values                                 | default value   |
 * |----------------------------------|-------------------------------------------------|-----------------|
 * | calculate_decays                 | 0 (no) or 1 (yes)                               | 1 (= enabled)   |
 * | include_higher_order_corrections | 0 (no) or 1 (yes)                               | 1 (= enabled)   |
 */
void FlexibleDecay_settings::reset()
{
   values[calculate_decays]                 = 1.;
   values[include_higher_order_corrections] = 1.;
   values[offshell_VV_decays]               = 2.;
   values[min_br_to_print]                  = 1e-5;
}
bool is_integer(double value)
{
   double intpart;
   return std::modf(value, &intpart) == 0.0;
}

}
