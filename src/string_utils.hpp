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

#ifndef STRING_UTILS_H
#define STRING_UTILS_H

#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>

namespace flexiblesusy {

/// concatenate strings
std::string concat(const std::vector<std::string>&);

/// concatenate strings with separator
template <typename T>
std::string concat(const std::vector<std::string>& strings, const T& separator)
{
   std::string result;

   for (const auto& s: strings)
      result += s + separator;

   return result;
}

template <typename T>
std::string to_string(T a)
{
   return boost::lexical_cast<std::string>(a);
}

} // namespace flexiblesusy

#endif
