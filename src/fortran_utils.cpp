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

#include "fortran_utils.hpp"
#include <unistd.h>
#include <iostream>

extern "C" {
void flush_impl();
}

namespace flexiblesusy
{
namespace futils
{

void swap() noexcept
{
   std::cout << std::flush;
   std::cerr << std::flush;
   int stdout_copy = dup(STDOUT_FILENO);
   dup2(STDERR_FILENO, STDOUT_FILENO);
   dup2(stdout_copy, STDERR_FILENO);
   close(stdout_copy);
}

void flush() noexcept { flush_impl(); }

} // namespace futils
} // namespace flexiblesusy
