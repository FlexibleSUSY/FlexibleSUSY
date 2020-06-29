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

#include "observable_problems_format.hpp"

#include <iostream>
#include <iterator>

namespace flexiblesusy {

std::ostream& operator<<(std::ostream& ostr, const Observable_problems& op)
{
   ostr << "Observable problems:\n";
   copy_problem_strings(op, std::ostream_iterator<std::string>(ostr, "\n"));
   return ostr;
}

} // namespace flexiblesusy
