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

#ifndef GSL_MULTIMIN_FMINIMIZER_H
#define GSL_MULTIMIN_FMINIMIZER_H

#include "gsl_vector.hpp"
#include <gsl/gsl_multimin.h>

namespace flexiblesusy {

/**
 * RAII wrapper for gsl_multimin_fminimizer
 */
class GSL_multimin_fminimizer
{
public:
   GSL_multimin_fminimizer(const gsl_multimin_fminimizer_type* type, std::size_t dim,
                           gsl_multimin_function* f, const GSL_vector& start,
                           const GSL_vector& step_size);
   GSL_multimin_fminimizer(const GSL_multimin_fminimizer&) = delete;
   GSL_multimin_fminimizer(GSL_multimin_fminimizer&&) = delete;
   ~GSL_multimin_fminimizer() noexcept;
   GSL_multimin_fminimizer& operator=(const GSL_multimin_fminimizer&) = delete;
   GSL_multimin_fminimizer& operator=(GSL_multimin_fminimizer&&) = delete;

   GSL_vector get_minimum_point() const;
   double get_minimum_value() const;
   int iterate();
   void print_state(std::size_t iteration) const;
   int test_residual(double precision) const noexcept;

private:
   gsl_multimin_fminimizer* solver = nullptr;
};

} // namespace flexiblesusy

#endif
