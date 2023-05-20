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

#include "gsl_multimin_fminimizer.hpp"
#include "logger.hpp"
#include "error.hpp"
#include <string>

namespace flexiblesusy {

GSL_multimin_fminimizer::GSL_multimin_fminimizer(
   const gsl_multimin_fminimizer_type* type, std::size_t dim,
   gsl_multimin_function* f, const GSL_vector& start,
   const GSL_vector& step_size)
{
   solver = gsl_multimin_fminimizer_alloc(type, dim);

   if (!solver) {
      throw OutOfMemoryError(
         std::string("Cannot allocate gsl_multimin_fminimizer ") +
         gsl_multimin_fminimizer_name(solver));
   }

   gsl_multimin_fminimizer_set(solver, f, start.raw(), step_size.raw());
}

GSL_multimin_fminimizer::~GSL_multimin_fminimizer() noexcept
{
   gsl_multimin_fminimizer_free(solver);
}

GSL_vector GSL_multimin_fminimizer::get_minimum_point() const
{
   return solver->x;
}

double GSL_multimin_fminimizer::get_minimum_value() const
{
   return solver->fval;
}

int GSL_multimin_fminimizer::iterate()
{
   return gsl_multimin_fminimizer_iterate(solver);
}

void GSL_multimin_fminimizer::print_state(std::size_t iteration) const
{
   VERBOSE_MSG("\t\t\tIteration " << iteration
               << ": x = " << GSL_vector(solver->x)
               << ", f(x) = " << solver->fval);
}

int GSL_multimin_fminimizer::test_residual(double precision) const noexcept
{
   const double size = gsl_multimin_fminimizer_size(solver);
   return gsl_multimin_test_size(size, precision);
}

} // namespace flexiblesusy
