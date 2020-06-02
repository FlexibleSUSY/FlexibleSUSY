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

#include "string_utils.hpp"

namespace flexiblesusy {

std::string concat(const std::vector<std::string>& strings)
{
   std::string result;

   for (const auto& s: strings)
      result += s;

   return result;
}


std::string concat(const std::vector<std::string>& strings, const std::string& separator)
{
   std::string result;

   for (auto it = strings.cbegin(), end = strings.end(); it != end; ++it) {
      if (it != strings.cbegin()) {
         result += separator;
      }
      result += *it;
   }

   return result;
}


std::string concat(const std::vector<std::string>& strings, char separator)
{
   return concat(strings, std::string(1, separator));
}

} // namespace flexiblesusy
