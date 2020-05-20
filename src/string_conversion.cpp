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

#include "string_conversion.hpp"
#include "error.hpp"
#include <string>

namespace flexiblesusy {

int to_int(const std::string& str)
{
   int i = 0;

   try {
      i = std::stoi(str);
   } catch (std::exception& e) {
      throw ReadError(e.what());
   }

   return i;
}

double to_double(const std::string& str)
{
   double d = 0.0;

   try {
      d = std::stod(str);
   } catch (std::exception& e) {
      throw ReadError(e.what());
   }

   return d;
}

} // namespace flexiblesusy
