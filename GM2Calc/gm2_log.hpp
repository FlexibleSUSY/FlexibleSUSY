// ====================================================================
// This file is part of GM2Calc.
//
// GM2Calc is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// GM2Calc is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GM2Calc.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#ifndef GM2_LOG_HPP
#define GM2_LOG_HPP

#include <iostream>

#define ERROR(message)                                                  \
   do { std::cerr << "Error: " << message << std::endl; } while (0)

#define VERBOSE(message)                                                \
   do { std::cerr << message << std::endl; } while (0)

#define WARNING(message)                                                \
   do { std::cerr << "Warning: " << message << std::endl; } while (0)

#endif
