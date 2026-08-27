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

#include "gsl_min_fminimizer.hpp"
#include "error.hpp"
#include <string>

namespace flexiblesusy {

GSL_min_fminimizer::GSL_min_fminimizer(
   const gsl_min_fminimizer_type* type,
   gsl_function* f, double start,
   double end)
{
   solver = gsl_min_fminimizer_alloc(type);

   if (!solver) {
      throw OutOfMemoryError(
         std::string("Cannot allocate gsl_multimin_fminimizer ") +
         gsl_min_fminimizer_name(solver));
   }

   gsl_min_fminimizer_set(solver, f, 0.5*(end-start), start, end);
}

} // namespace flexiblesusy
