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
#include <cerrno>
#include <climits>
#include <cstdlib>

namespace flexiblesusy {

int to_int(const char* s)
{
   char* end = nullptr;
   errno = 0;

   const long l = std::strtol(s, &end, 10);

   if ((errno == ERANGE && l == LONG_MAX) || l > INT_MAX) {
      errno = 0;
      throw ReadError("range overflow occurred in conversion to int");
   }
   if ((errno == ERANGE && l == LONG_MIN) || l < INT_MIN) {
      errno = 0;
      throw ReadError("range underflow occurred in conversion to int");
   }
   if (*s == '\0' || *end != '\0') {
      errno = 0;
      throw ReadError("cannot convert string to int");
   }

   errno = 0;

   return static_cast<int>(l);
}


long to_long(const char* s)
{
   char* end = nullptr;
   errno = 0;

   const long l = std::strtol(s, &end, 10);

   if (errno == ERANGE && l == LONG_MAX) {
      errno = 0;
      throw ReadError("range overflow occurred in conversion to long");
   }
   if (errno == ERANGE && l == LONG_MIN) {
      errno = 0;
      throw ReadError("range underflow occurred in conversion to long");
   }
   if (*s == '\0' || *end != '\0') {
      errno = 0;
      throw ReadError("cannot convert string to long");
   }

   errno = 0;

   return l;
}


double to_double(const char* s)
{
   char* end = nullptr;
   errno = 0;

   const double d = std::strtod(s, &end);

   if (errno == ERANGE) {
      errno = 0;
      throw ReadError("range error occurred in conversion to double");
   }
   if (*s == '\0' || *end != '\0') {
      errno = 0;
      throw ReadError("cannot convert string to double");
   }

   errno = 0;

   return d;
}

} // namespace flexiblesusy
