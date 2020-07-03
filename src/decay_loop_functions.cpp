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

#include "decay_loop_functions.hpp"
#include <cmath>
#include <complex>
#include <limits>

namespace flexiblesusy {

namespace {

/// Eq.(2.31) of hep-ph/0503172
double RT_general(double x) noexcept
{
   const std::complex<double> z(x, 0.0);

   return std::real(
      3.*(1. - 8.*z + 20.*z*z)/std::sqrt(4.*z - 1.)*std::acos((3.*z - 1.)/(2.*std::pow(z, 3./2.)))
      - (1. - z)/(2.*z)*(2. - 13.*z + 47.*z*z)
      - 3./2. * (1. - 6.*z + 4.*z*z)*std::log(z)
   );
}

} // anonymous namespace

/// Eq.(2.31) of hep-ph/0503172, including edge cases
double RT(double x) noexcept
{
   if (x < 0.25) {
      return RT_general(x);
   } else if (x == 0.25) {
      return std::numeric_limits<double>::quiet_NaN();
   } else if (x < 0.95) {
      return RT_general(x);
   } else if (x < 1) {
      const double d = x - 1;
      const double d2 = d*d;
      const double d4 = d2*d2;
      const double d5 = d4*d;
      return d5*(-3./10. + d*(13./20. + d*(-15./14. + d*(447./280. + d*(-95./42. + 523.*d/168)))));
   } else if (x == 1) {
      return 0;
   } else if (x < 1.01) {
      const double d = x - 1;
      return d*(39 + d*(75.0/2.0 + d*(-6 + 5.0/4.0*d)));
   }
   return RT_general(x);
}

} // namespace flexiblesusy
