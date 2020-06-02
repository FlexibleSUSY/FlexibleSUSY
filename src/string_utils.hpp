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
#include <sstream>
#include <vector>

namespace flexiblesusy {

/// concatenate strings
std::string concat(const std::vector<std::string>&);


/// concatenate strings with separator
std::string concat(const std::vector<std::string>&, const std::string&);


/// concatenate strings with separator
std::string concat(const std::vector<std::string>&, char);


/// concatenate strings with separator
template <class InputIterator>
std::string concat(InputIterator start, InputIterator end, const std::string& separator)
{
   std::ostringstream result;

   for (InputIterator it = start; it != end; ++it) {
      if (it != start) {
         result << separator;
      }
      result << *it;
   }

   return result.str();
}


/// concatenate strings with separator
template <class InputIterator>
std::string concat(InputIterator start, InputIterator end, char separator)
{
   return concat(start, end, std::string(1, separator));
}

} // namespace flexiblesusy

#endif
