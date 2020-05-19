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

#ifndef STRING_FORMAT_H
#define STRING_FORMAT_H

#include <complex>
#include <string>

namespace flexiblesusy {

std::string to_string(char);
std::string to_string(unsigned char);
std::string to_string(unsigned short);
std::string to_string(unsigned int);
std::string to_string(unsigned long);
std::string to_string(unsigned long long);
std::string to_string(signed char);
std::string to_string(signed short);
std::string to_string(signed int);
std::string to_string(signed long);
std::string to_string(signed long long);
std::string to_string(double);
std::string to_string(const std::complex<double>&);

} // namespace flexiblesusy

#endif
