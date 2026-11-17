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

#ifndef GSL_MIN_FMINIMIZER_H
#define GSL_MIN_FMINIMIZER_H

#include "gsl_vector.hpp"
#include <gsl/gsl_min.h>

namespace flexiblesusy {

class GSL_min_fminimizer
{
public:
   GSL_min_fminimizer(const gsl_min_fminimizer_type* type,
                      gsl_function* f, double start,
                      double end);
   GSL_min_fminimizer(const GSL_min_fminimizer&) = delete;
   GSL_min_fminimizer(GSL_min_fminimizer&&) = delete;
   ~GSL_min_fminimizer() noexcept;
   GSL_min_fminimizer& operator=(const GSL_min_fminimizer&) = delete;
   GSL_min_fminimizer& operator=(GSL_min_fminimizer&&) = delete;

   GSL_vector get_minimum_point() const;
   double get_minimum_value() const;
   int iterate();
   void print_state(std::size_t iteration) const;
   int test_residual(double precision) const noexcept;

private:
   gsl_min_fminimizer* solver = nullptr;
};

} // namespace flexiblesusy

#endif
