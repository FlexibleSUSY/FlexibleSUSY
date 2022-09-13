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

#pragma once
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

namespace flexiblesusy {
namespace test {

const char PATH_SEPARATOR =
#ifdef _WIN32
   '\\';
#else
   '/';
#endif

/**
 * Reads real numbers from a file line by line into a vector of vectors of T.
 *
 * @param filename file name
 * @tparam T data type
 *
 * @return vector of vectors of doubles.
 */
template <typename T>
std::vector<std::vector<T>>
read_from_file(const std::string& filename)
{
   std::vector<std::vector<T>> data;
   std::string line;
   std::ifstream fstr(filename);

   while (std::getline(fstr, line)) {
      std::istringstream iss(line);
      data.emplace_back(std::vector<T>{std::istream_iterator<T>(iss),
                                       std::istream_iterator<T>()});
   }

   return data;
}

} // namespace test
} // namespace flexiblesusy
