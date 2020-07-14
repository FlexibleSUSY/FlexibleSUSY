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

#include "gsl_multiroot_fsolver.hpp"
#include "logger.hpp"
#include "error.hpp"
#include <string>

namespace flexiblesusy {

GSL_multiroot_fsolver::GSL_multiroot_fsolver(
   const gsl_multiroot_fsolver_type* type, std::size_t dim,
   gsl_multiroot_function* f, const GSL_vector& start)
{
   solver = gsl_multiroot_fsolver_alloc(type, dim);

   if (!solver) {
      throw OutOfMemoryError(
         std::string("Cannot allocate gsl_multiroot_fsolver ") +
         gsl_multiroot_fsolver_name(solver));
   }

   gsl_multiroot_fsolver_set(solver, f, start.raw());
}

GSL_multiroot_fsolver::~GSL_multiroot_fsolver() noexcept
{
   gsl_multiroot_fsolver_free(solver);
}

GSL_vector GSL_multiroot_fsolver::get_root() const { return solver->x; }

int GSL_multiroot_fsolver::iterate()
{
   return gsl_multiroot_fsolver_iterate(solver);
}

void GSL_multiroot_fsolver::print_state(std::size_t iteration) const
{
   VERBOSE_MSG("\t\t\tIteration " << iteration
                                  << ": x = " << GSL_vector(solver->x)
                                  << ", f(x) = " << GSL_vector(solver->f));
}

int GSL_multiroot_fsolver::test_residual(double precision) const noexcept
{
   return gsl_multiroot_test_residual(solver->f, precision);
}

} // namespace flexiblesusy
