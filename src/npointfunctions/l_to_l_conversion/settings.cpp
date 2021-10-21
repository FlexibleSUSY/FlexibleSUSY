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
   "scalar form factor G^S_{p,u} without masses",
   "scalar form factor G^S_{n,u} without masses",
   "scalar form factor G^S_{p,d} without masses",
   "scalar form factor G^S_{n,d} without masses",
   "scalar form factor G^S_{p,s} without masses",
   "scalar form factor G^T_{n,s} without masses",
   "vector form factor G^V_{p,u}",
   "vector form factor G^V_{n,u}",
   "vector form factor G^V_{p,d}",
   "vector form factor G^V_{n,d}",
   "tensor form factor G^T_{p,u}",
   "tensor form factor G^T_{n,u}",
   "tensor form factor G^T_{p,d}",
   "tensor form factor G^T_{n,d}",
   "tensor form factor G^T_{p,s}",
   "tensor form factor G^T_{n,s}"
};

void assert_bool(double value, const char* quantity)
{
   if (value != 0.0 && value != 1.0) {
      throw SetupError(std::string(quantity) + " must be either 0 or 1");
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
   case include_tensor_contribution: // [bool]
      assert_bool(value, descriptions.at(o).c_str());
      break;
   case include_gluonic_contribution: // [bool]
      assert_bool(value, descriptions.at(o).c_str());
      break;
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
 * Resets all settings to their defaults.
 *
 * | enum                             | possible values   | default value     |
 * |----------------------------------|-------------------|-------------------|
 * | include_tensor_contribution      | 0 (no) or 1 (yes) | 0 (= disabled)    |
 * | include_gluonic_contribution     | 0 (no) or 1 (yes) | 0 (= disabled)    |
 * | other coefficients               | any double        | see function body |
 */
void LToLConversion_settings::reset()
{
   values[include_tensor_contribution]  = 0.0;
   values[include_gluonic_contribution] = 0.0;
   values[scalar_pu] = 0.021;
   values[scalar_nu] = 0.019;
   values[scalar_pd] = 0.041;
   values[scalar_nd] = 0.045;
   values[scalar_ps] = 0.043;
   values[scalar_ns] = 0.043;
   values[vector_pu] = 2.0;
   values[vector_nu] = 1.0;
   values[vector_pd] = 1.0;
   values[vector_nd] = 2.0;
   values[tensor_pu] = 0.77;
   values[tensor_nu] = -0.23;
   values[tensor_pd] = -0.23;
   values[tensor_nd] = 0.77;
   values[tensor_ps] = 0.008;
   values[tensor_ns] = 0.008;
}

} // namespace flexiblesusy
