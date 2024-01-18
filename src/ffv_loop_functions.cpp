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

#include "ffv_loop_functions.hpp"
#include "numerics2.hpp"
#include "wrappers.hpp"
#include <cmath>

namespace flexiblesusy {
namespace ffv_loop_functions {
namespace one_loop {

/**
* @defgroup LoopFunctions 1-loop functions
* @brief The 1-loop functions necessary for the Fe_I -> Fe_J gamma 1-loop calculations.
*
* These are FA(), FB() as specified in arXiv:0808.1819
*/

// function from eq. 15 of hep-ph/9510309
double FA(double r)
{
   if (is_zero(1.0 - r)) {
      return 1.5;
   } else if (is_zero(r)) {
      return 2.0;
   } else {
      return (2.0 - 9.0 * r + 18.0 * r * r - 11.0 * r * r * r +
              6.0 * r * r * r * std::log(r)) /
             Power4(1.0 - r);
   }
}

// function from eq. 16 of hep-ph/9510309
double FB(double r)
{
   const double y = r - 1.0;
   if (is_zero(r)) {
      return 2.0;
   } else if (std::abs(y) < 0.23) {
      // error around x=1 is <= 10^-12 on an intel i7
      return (1.0000000000000000000 -
              0.4000000000000000000  * y +
              0.2000000000000000000  * y * y -
              0.11428571428571428571 * y * y * y +
              0.07142857142857142857 * y * y * y * y -
              0.04761904761904761905 * y * y * y * y * y +
              0.03333333333333333333 * y * y * y * y * y * y -
              0.02424242424242424242 * y * y * y * y * y * y * y +
              0.0181818181818181818  * y * y * y * y * y * y * y * y -
              0.01398601398601398601 * y * y * y * y * y * y * y * y * y +
              0.01098901098901098901 * y * y * y * y * y * y * y * y * y * y -
              0.0087912087912087912  * y * y * y * y * y * y * y * y * y * y * y +
              0.00714285714285714286 * y * y * y * y * y * y * y * y * y * y * y * y -
              0.0058823529411764706  * y * y * y * y * y * y * y * y * y * y * y * y * y +
              0.0049019607843137255  * y * y * y * y * y * y * y * y * y * y * y * y * y * y -
              0.0041279669762641899  * y * y * y * y * y * y * y * y * y * y * y * y * y * y * y);
   } else {
      return 2. *
             (1.0 - 6.0 * r + 3.0 * r * r + 2.0 * r * r * r -
              6.0 * r * r * std::log(r)) /
             Power4(1.0 - r);
   }
}

// function from eq. 16 of hep-ph/9510309
double FC(double r)
{
   const double y = r - 1.0;
   if (is_zero(r)) {
      return 3.0;
   } else if (std::isinf(r)) {
      return 0.0;
   } else if (std::abs(y) < 0.185) {
      // error around x=1 is <= 10^-13 on an intel i7
      return (1.0000000000000000000 -
              0.50000000000000000000 * y +
              0.30000000000000000000 * y * y -
              0.2000000000000000000  * y * y * y +
              0.14285714285714285714 * y * y * y * y -
              0.10714285714285714286 * y * y * y * y * y +
              0.08333333333333333333 * y * y * y * y * y * y -
              0.06666666666666666667 * y * y * y * y * y * y * y +
              0.05454545454545454545 * y * y * y * y * y * y * y * y -
              0.0454545454545454545  * y * y * y * y * y * y * y * y * y +
              0.0384615384615384615  * y * y * y * y * y * y * y * y * y * y -
              0.03296703296703296703 * y * y * y * y * y * y * y * y * y * y * y +
              0.0285714285714285714  * y * y * y * y * y * y * y * y * y * y * y * y -
              0.02500000000000000000 * y * y * y * y * y * y * y * y * y * y * y * y * y +
              0.0220588235294117647  * y * y * y * y * y * y * y * y * y * y * y * y * y * y -
              0.0196078431372549020  * y * y * y * y * y * y * y * y * y * y * y * y * y * y * y);
   } else {
      return 3. * (1.0 - r * r + 2.0 * r * std::log(r)) / Power3(1.0 - r);
   }
}

// function from eq. 17 of hep-ph/9510309
double FD(double r)
{
   if (std::isinf(r)) {
      return 0.0;
   } else if (is_zero(1.0 - r)) {
      return -4.5;
   } else {
      return (16.0 - 45.0 * r + 36.0 * r * r - 7.0 * r * r * r +
              6.0 * (2.0 - 3.0 * r) * std::log(r)) /
             Power4(1.0 - r);
   }
}

// function from eq. 18 of hep-ph/9510309
double FE(double r)
{
   const double y = r - 1.0;
   if (is_zero(r)) {
      return 4.0;
   } else if (std::abs(y) < 0.21) {
      // error around x=1 is <= 10^-12 on an intel i7
      return (1.0000000000000000000 -
              0.60000000000000000000  * y +
              0.40000000000000000000  * y * y -
              0.28571428571428571429  * y * y * y +
              0.21428571428571428571  * y * y * y * y -
              0.16666666666666666667  * y * y * y * y * y +
              0.13333333333333333333  * y * y * y * y * y * y -
              0.10909090909090909091  * y * y * y * y * y * y * y +
              0.090909090909090909091 * y * y * y * y * y * y * y * y -
              0.076923076923076923077 * y * y * y * y * y * y * y * y * y +
              0.065934065934065934066 * y * y * y * y * y * y * y * y * y * y -
              0.057142857142857142857 * y * y * y * y * y * y * y * y * y * y * y +
              0.050000000000000000000 * y * y * y * y * y * y * y * y * y * y * y * y -
              0.044117647058823529412 * y * y * y * y * y * y * y * y * y * y * y * y * y +
              0.039215686274509803922 * y * y * y * y * y * y * y * y * y * y * y * y * y * y -
              0.035087719298245614035 * y * y * y * y * y * y * y * y * y * y * y * y * y * y * y);
   } else {
      return 2. *
             (2.0 + 3.0 * r - 6.0 * r * r + r * r * r + 6.0 * r * std::log(r)) /
             Power4(1.0 - r);
   }
}

// function from eq. 18 of hep-ph/9510309
double FF(double r)
{
   const double y = r - 1.0;
   if (std::abs(y) < 0.155) {
      // error around x=1 is <= 10^-13 on an intel i7
      return (1.0 -
              0.75 * y +
              0.6 * y * y -
              0.50000000000000000000 * y * y * y +
              0.4285714285714285714  * y * y * y * y -
              0.37500000000000000000 * y * y * y * y * y +
              0.33333333333333333333 * y * y * y * y * y * y -
              0.3000000000000000000  * y * y * y * y * y * y * y +
              0.2727272727272727273  * y * y * y * y * y * y * y * y -
              0.2500000000000000000  * y * y * y * y * y * y * y * y * y +
              0.23076923076923076923 * y * y * y * y * y * y * y * y * y * y -
              0.21428571428571428571 * y * y * y * y * y * y * y * y * y * y * y +
              0.2000000000000000000  * y * y * y * y * y * y * y * y * y * y * y * y -
              0.1875000000000000000  * y * y * y * y * y * y * y * y * y * y * y * y * y +
              0.1764705882352941176  * y * y * y * y * y * y * y * y * y * y * y * y * y * y -
              0.16666666666666666667 * y * y * y * y * y * y * y * y * y * y * y * y * y * y * y);
   } else {
      return 3. / 2. * (-3.0 + 4.0 * r - r * r - 2.0 * std::log(r)) /
             Power3(1.0 - r);
   }
}

//  One loop function for the form factor F1 which is the coefficient
//  of ubar(pp,mj) \gamma^\mu u(p,mi), for the FFV_VVF contribution.
double FH(double r)
{
   if (is_zero(r)) {
      return 41.0;
   } else if (is_zero(1.0 - r)){
      return 27.0;
   } else {
      return (41.0 - 207.0 * r + 279.0 * r * r - 113.0 * r * r * r + 6.0 * r * r * (-15.0 + 13.0 * r)
              * std::log(r)) / Power4(1.0 - r);
      //return (5.0 - 27.0 * r + 27.0 * r * r - 5.0 * r * r * r + 6.0 * (r - 3.0) * r * r * std::log(r)) /
      //        Power4(1.0 - r);
   }
}

//  One loop function for the form factor F1 which is the coefficient
//  of ubar(pp,mj) \gamma^\mu u(p,mi), for the FFV_VVF contribution.
double FI(double r)
{
   if (is_zero(r)) {
      return 1.0;
   } else if (is_zero(1.0 - r)){
      return 0.5;
   } else {
      return (1.0 - r + r * std::log(r)) / Power2(1.0 - r);
   }
}

// function from eq. 78 of hep-ph/2104.03691
double FJ(double r)
{
   if (is_zero(r)) {
      return 7.0;
   } else if (is_zero(1.0 - r)){
      return 5.0;
   } else {
      return (7.0 - 33.0 * r + 57.0 * r * r - 31.0 * r * r * r + 6.0 * r * r * (3.0 * r - 1.0)
                   * std::log(r)) / Power4(1.0 - r);
   }
}

// function from eq. 79 of hep-ph/2104.03691
double FK(double r)
{
   if (is_zero(r)) {
      return 1.0;
   } else if (is_zero(1.0 - r)){
      return 2.0 / 3.0;
   } else {
      return (1.0 - 4.0 * r + 3.0 * r * r - 2 * r * r * std::log(r)) / Power3(1.0 - r);
   }
}

//  One loop function for the form factor F1 which is the coefficient
//  of ubar(pp,mj) \gamma^\mu u(p,mi), for the FFV_FFV contribution.
double FL(double r)
{
   if (is_zero(1.0 - r)){
      return 10.5;
   } else if (std::isinf(r)) {
      return 0.0;
   } else {
      return (2.0 + 27.0 * r - 54.0 * r * r + 25.0 * r * r * r - 6.0 * (2.0 - 9.0 * r + 6.0 * r * r)
                 * std::log(r)) / Power4(1.0 - r);
   }
}

// function from eq. 80 of hep-ph/2104.03691
double FM(double r)
{
   if (is_zero(r)) {
      return 4.0;
   } else if (is_zero(1.0 - r)){
      return 1.5;
   } else if (std::isinf(r)) {
      return 0.0;
   } else {
      return (4.0 - 9.0 * r + 5.0 * r * r * r + 6.0 * (1.0 - 2.0 * r) * r * std::log(r)) /
              Power4(1.0 - r);
   }
}

// function from eq. 81 of hep-ph/2104.03691
double FN(double r, double s)
{
   if (is_zero(r - s)) {
      if (is_zero(r)) {
         return 1.0 / 2.0;
      } else if (is_zero(1.0 - r)){
         return 1.0 / 3.0;
      } else {
         return (1.0 - 4.0 * r + 3 * r * r - 2 * r * r * std::log(r)) / 2.0 / Power3(1.0 - r);
      }
   } else {
      return s * (r * Power2(r - s) * std::log(r) + (r - 1.0) * ((r - s) * (s - 1.0) - (r - 1.0)
                              * r * std::log(r/s))) / Power2(1.0 - r) / Power2(r - s) / (s - 1.0);
   }
}

} // namespace one_loop
} // namespace ffv_loop_functions
} // namespace flexiblesusy
