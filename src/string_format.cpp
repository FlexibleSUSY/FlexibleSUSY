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

#include "string_format.hpp"
#include <boost/lexical_cast.hpp>

namespace flexiblesusy {

#define DEFINE_TO_STRING(type)                     \
   std::string to_string(type a)                   \
   {                                               \
      return boost::lexical_cast<std::string>(a);  \
   }

DEFINE_TO_STRING(char              );
DEFINE_TO_STRING(unsigned char     );
DEFINE_TO_STRING(unsigned short    );
DEFINE_TO_STRING(unsigned int      );
DEFINE_TO_STRING(unsigned long     );
DEFINE_TO_STRING(unsigned long long);
DEFINE_TO_STRING(signed char       );
DEFINE_TO_STRING(signed short      );
DEFINE_TO_STRING(signed int        );
DEFINE_TO_STRING(signed long       );
DEFINE_TO_STRING(signed long long  );
DEFINE_TO_STRING(double            );

std::string to_string(const std::complex<double>& a)
{
   return "(" + to_string(std::real(a)) + "," + to_string(std::imag(a)) + ")";
}

} // namespace flexiblesusy
