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

#ifndef GSL_MULTIROOT_FSOLVER_H
#define GSL_MULTIROOT_FSOLVER_H

#include "gsl_vector.hpp"
#include <gsl/gsl_multiroots.h>

namespace flexiblesusy {

class GSL_multiroot_fsolver
{
public:
   GSL_multiroot_fsolver(const gsl_multiroot_fsolver_type* type, std::size_t dim,
                         gsl_multiroot_function* f, const GSL_vector& start);
   GSL_multiroot_fsolver(const GSL_multiroot_fsolver&) = delete;
   GSL_multiroot_fsolver(GSL_multiroot_fsolver&&) = delete;
   ~GSL_multiroot_fsolver() noexcept;
   GSL_multiroot_fsolver& operator=(const GSL_multiroot_fsolver&) = delete;
   GSL_multiroot_fsolver& operator=(GSL_multiroot_fsolver&&) = delete;

   GSL_vector get_root() const;
   int iterate();
   void print_state(std::size_t iteration) const;
   int test_residual(double precision) const noexcept;

private:
   gsl_multiroot_fsolver* solver = nullptr;
};

} // namespace flexiblesusy

#endif
