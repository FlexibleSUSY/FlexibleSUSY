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

#include "threshold_loop_functions.hpp"
#include "dilog.hpp"
#include "loop_libraries/loop_library.hpp"
#include "logger.hpp"
#include "numerics.h"

#include <cmath>
#include <limits>

namespace flexiblesusy {
namespace threshold_loop_functions {

namespace {
   template <typename T> T sqr(T x) noexcept { return x*x; }
   template <typename T> T cube(T x) noexcept { return x*x*x; }
   template <typename T> T quad(T x) noexcept { return x*x*x*x; }
   template <typename T> T pow5(T x) noexcept { return x*x*x*x*x; }
   template <typename T> T pow6(T x) noexcept { return x*x*x*x*x*x; }
   template <typename T> T pow7(T x) noexcept { return x*x*x*x*x*x*x; }
   template <typename T> T pow8(T x) noexcept { return x*x*x*x*x*x*x*x; }
   template <typename T> T pow9(T x) noexcept { return x*x*x*x*x*x*x*x*x; }

   template <typename T>
   bool is_zero(T a, T prec) noexcept
   {
      return std::fabs(a) < prec;
   }

   template <typename T>
   bool is_equal(T a, T b, T prec) noexcept
   {
      return is_zero(a - b, prec);
   }

   template <typename T>
   bool is_equal_rel(T a, T b, T prec) noexcept
   {
      if (is_equal(a, b, std::numeric_limits<T>::epsilon())) {
         return true;
      }

      const T min = std::min(std::abs(a), std::abs(b));

      if (min < std::numeric_limits<T>::epsilon()) {
         return is_equal(a, b, prec);
      }

      const double max = std::max(std::fabs(a), std::fabs(b));

      return is_equal(a, b, prec*max);
   }

   bool is_close(double a, double b, double prec) noexcept
   {
      const double max = std::max(std::abs(a), std::abs(b));
      return is_zero(a - b, prec*(1.0 + max));
   }

   double logx(double x) noexcept
   {
      if (is_zero(x, 10.0*std::numeric_limits<double>::epsilon())) {
         return 0.;
      }

      return std::log(x);
   }

   double xlogx(double x) noexcept
   {
      if (is_zero(x, 1e-14)) {
         return 0.;
      }

      return x * std::log(x);
   }

} // anonymous namespace

double F1(double x) noexcept
{
   const double eps = std::numeric_limits<double>::epsilon();

   if (is_zero(x, eps)) {
      return 0.;
   }

   if (is_equal(x, 1., 0.01)) {
      const double d = x - 1;
      const double d2 = sqr(d);
      return 1 + d2*(-1/6. + d*(1/6. - 2/15.*d));
   }

   const double x2 = sqr(x);

   return x*std::log(x2)/(x2 - 1);
}

double F2(double x) noexcept
{
   const double eps = std::numeric_limits<double>::epsilon();

   if (is_zero(x, eps)) {
      return 0.;
   }

   if (is_equal(x, 1., 0.01)) {
      const double d = (x - 1)*(x + 1);
      const double d2 = sqr(d);
      return 1 + d2*(-0.1 + d*(0.1 - 3/35.*d));
   }

   const double x2 = sqr(x);

   return 6*x2*(2 - 2*x2 + (1 + x2)*std::log(x2))/cube(x2 - 1);
}

double F3(double x) noexcept
{
   const double eps = std::numeric_limits<double>::epsilon();

   if (is_zero(x, eps)) {
      return 0.;
   }

   if (is_equal(x, 1., 0.01)) {
      const double d = x - 1;
      return 1 + d*(5/9. + d*(-4/9. + d*(2/9. - 7/90.*d)));
   }

   const double x2 = sqr(x);

   return 2*x*(5*(1 - x2) + (1 + 4*x2)*std::log(x2))/(3*sqr(x2 - 1));
}

double F4(double x) noexcept
{
   const double eps = std::numeric_limits<double>::epsilon();

   if (is_zero(x, eps)) {
      return 0.;
   }

   if (is_equal(x, 1., 0.01)) {
      const double d = x - 1;
      const double d2 = sqr(d);
      return 1 + d*(-1/3. + d2*(2/15. - d/6.));
   }

   const double x2 = sqr(x);

   return 2*x*(x2 - 1 - std::log(x2))/sqr(x2 - 1);
}

double F5(double x) noexcept
{
   const double eps = std::numeric_limits<double>::epsilon();

   if (is_zero(x, eps)) {
      return 0.;
   }

   if (is_equal(x, 1., 0.01)) {
      const double d = x - 1;
      const double d2 = sqr(d);
      return 1 + d2*(-0.3 + d*(0.3 - 3/14.*d));
   }

   if (is_equal(x, -1., 0.01)) {
      const double d = x + 1;
      const double d2 = sqr(d);
      return -1 + d2*(0.3 + d*(0.3 + 3/14.*d));
   }

   const double x2 = sqr(x);
   const double x4 = sqr(x2);

   return 3*x*(1 - x4 + 2*x2*std::log(x2))/cube(1 - x2);
}

double F6(double x) noexcept
{
   const double eps = std::numeric_limits<double>::epsilon();

   if (is_zero(x, eps)) {
      return -0.75;
   }

   const double x2 = sqr(x);

   if (is_equal(x2, 1., 0.01)) {
      const double d = (x - 1)*(x + 1);
      return d*(1/3. + d*(-1/8. + d*(1/15. + d/24.)));
   }

   return (x2 - 3)/(4*(1 - x2)) + x2*(x2 - 2)/(2*sqr(1 - x2))*std::log(x2);
}

double F7(double x) noexcept
{
   const double eps = std::numeric_limits<double>::epsilon();

   if (is_zero(x, eps)) {
      return -1.5;
   }

   const double x2 = sqr(x);

   if (is_equal(x2, 1., 0.01)) {
      const double d = (x - 1)*(x + 1);
      return 1 + d*(3/2. + d*(-9/20. + d*(0.2 - 3/28.*d)));
   }

   const double x4 = sqr(x2);

   return -3*(x4 - 6*x2 + 1)/(2*sqr(x2 - 1))
      + 3*x4*(x2 - 3)/cube(x2 - 1)*std::log(x2);
}

/// F8(x1,x2) in the limit x1 -> 1 and x2 -> 1
static double F8_1_1(double x1, double x2) noexcept
{
   const double d1 = (x1 - 1)*(x1 + 1);
   const double d2 = (x2 - 1)*(x2 + 1);
   const double d12 = sqr(d1);
   const double d22 = sqr(d2);
   const double d13 = d1*d12;
   const double d23 = d2*d22;
   const double d14 = sqr(d12);
   const double d24 = sqr(d22);

   return 1 + 2/3.*(d1 + d2) + 1/6.*(-d12 - d1*d2 - d22)
      + (d13 + d12*d2 + d1*d22 + d23)/15.
      + (-d14 - d13*d2 - d12*d22 - d1*d23 - d24)/30.;
}

/// F8(x1,x2) in the limit x1 -> 1, x2 != 1
static double F8_1_x2(double x1, double x2) noexcept
{
   const double x22 = sqr(x2);
   const double d = (x1 - 1)*(x1 + 1);

   if (is_zero(x2, 0.01)) {
      return 2*x22 + d*(1 - x22 + d*(-1/3. + 2/3.*x22));
   }

   const double d2 = sqr(d);
   const double d3 = d*d2;
   const double d4 = sqr(d2);
   const double y = (x2 - 1)*(x2 + 1);
   const double y2 = sqr(y);
   const double y3 = y*y2;
   const double y4 = sqr(y2);
   const double y5 = y2*y3;
   const double y6 = sqr(y3);
   const double lx22 = std::log(x22);

   return
      - 2 + 2*(1 + x22*(-1 + lx22*x22))/y2
      + d*(-1 + x22*(4 + x22*(-3 + 2*lx22)))/y3
      + d2*(-1 + x22*(6 + x22*(-3 + 6*lx22 + -2*x22)))/(3*y4)
      + d3*(-1 + x22*(8 + x22*(12*lx22 + x22*(-8 + x22))))/(6*y5)
      + d4*(-3 + x22*(30 + x22*(20 + 60*lx22 + x22*(-60 + x22*(15 - 2*x22)))))/(30*y6);
}

/// F8(x1,x2) in the limit x1 -> 0, x2 != 0
static double F8_0_x2(double x1, double x2) noexcept
{
   const double x12 = sqr(x1);
   const double x22 = sqr(x2);
   const double lx22 = std::log(x22);

   return 2*(1 - x22 + (x12 + x22)*lx22)/(-1 + x22);
}

// F8(x1,x2) in the limit x1 -> x2, x2 != 1
static double F8_x1_x2(double x1, double x2) noexcept
{
   const double x22 = sqr(x2);
   const double d = (x1 - x2)*(x1 + x2);
   const double d2 = sqr(d);
   const double d3 = d*d2;
   const double y = (x2 - 1)*(x2 + 1);
   const double y2 = sqr(y);
   const double y3 = y*y2;
   const double y4 = sqr(y2);
   const double y5 = y2*y3;
   const double lx22 = std::log(x22);

   return
      + 2*((-2 + x22)*x22*lx22 + y)/y2
      + d*(3 + x22*(-4 + x22) + 2*lx22)/y3
      + d2*(-2 + x22*(-3 - 6*lx22 + x22*(6 - x22)))/(3*x22*y4)
      + d3*(-1 + x22*(8 + x22*(12*lx22 + x22*(-8 + x22))))/(6*x22*x22*y5);
}

double F8(double x1, double x2) noexcept
{
   const double eps = 10.*std::numeric_limits<double>::epsilon();

   if (is_zero(x1, eps) && is_zero(x2, eps)) {
      return -2.;
   }

   const double ax1 = std::fabs(x1);
   const double ax2 = std::fabs(x2);

   if (is_equal(ax1, 1., 0.01) && is_equal(ax2, 1., 0.01)) {
      return F8_1_1(x1, x2);
   }

   if (is_equal(ax1, 1., 0.01)) {
      return F8_1_x2(x1, x2);
   }

   if (is_equal(ax2, 1., 0.01)) {
      return F8_1_x2(x2, x1);
   }

   if (is_zero(x1, 0.00001)) {
      return F8_0_x2(x1, x2);
   }

   if (is_zero(x2, 0.00001)) {
      return F8_0_x2(x2, x1);
   }

   if (is_equal(ax1, ax2, 0.00001)) {
      return F8_x1_x2(x1, x2);
   }

   const double x12 = sqr(x1);
   const double x22 = sqr(x2);
   const double x14 = sqr(x12);
   const double x24 = sqr(x22);

   return -2. + 2./(x12 - x22)*(
      + x14/(x12 - 1)*std::log(x12)
      - x24/(x22 - 1)*std::log(x22));
}

/// F9(x1,x2) in the limit x1 -> 1 and x2 -> 1
static double F9_1_1(double x1, double x2) noexcept
{
   const double d1 = (x1 - 1)*(x1 + 1);
   const double d2 = (x2 - 1)*(x2 + 1);
   const double d12 = sqr(d1);
   const double d22 = sqr(d2);
   const double d13 = d1*d12;
   const double d23 = d2*d22;
   const double d14 = sqr(d12);
   const double d24 = sqr(d22);

   return 1
      + 1/3.*(-d2 - d1)
      + 1/6.*(d12 + d1*d2 + d22)
      + 1/10.*(-d13 - d12*d2 - d1*d22 - d23)
      + 1/15.*(d14 + d13*d2 + d12*d22 + d1*d23 + d24);
}

/// F9(x1,x2) in the limit x1 -> 1
static double F9_1_x2(double x1, double x2) noexcept
{
   const double x22 = sqr(x2);
   const double d = (x1 - 1)*(x1 + 1);
   const double d2 = sqr(d);
   const double d3 = d*d2;
   const double d4 = sqr(d2);
   const double y = (x2 - 1)*(x2 + 1);
   const double y2 = sqr(y);
   const double y3 = y*y2;
   const double y4 = sqr(y2);
   const double y5 = y2*y3;
   const double y6 = sqr(y3);
   const double lx22 = std::log(x22);

   if (is_zero(x2, 0.01)) {
      return 2 + d;
   }

   return
      + 2*(1 + x22*(-1 + lx22))/y2
      + d*(1 + x22*(2*lx22 - x22))/y3
      + d2*(2 + x22*(3 + 6*lx22 + x22*(-6 + x22)))/(3*y4)
      + d3*(3 + x22*(10 + 12*lx22 + x22*(-18 + x22*(6 - x22))))/(6*y5)
      + d4*(12 + x22*(65 + 60*lx22 + x22*(-120 + x22*(60 + x22*(-20 + 3*x22)))))/(30*y6);
}

/// F9(x1,x2) in the limit x1 -> 0, x2 != 1, x2 != 0
static double F9_0_x2(double /* x1 */, double x2) noexcept
{
   const double x22 = sqr(x2);

   return 2*std::log(x22)/(-1 + x22);
}

/// F9(x1,x2) in the limit x1 -> x2, x2 != 0
static double F9_x1_x2(double x1, double x2) noexcept
{
   const double x22 = sqr(x2);
   const double d = (x1 - x2)*(x1 + x2);
   const double d2 = sqr(d);
   const double d3 = d*d2;
   const double y = (x2 - 1)*(x2 + 1);
   const double y2 = sqr(y);
   const double y3 = y*y2;
   const double y4 = sqr(y2);
   const double y5 = y2*y3;
   const double lx22 = std::log(x22);

   return
      + 2*(y - lx22)/y2
      + d*(1 + x22*(2*lx22 - x22))/(x22*y3)
      + d2*(1 + x22*(-6 + x22*(3 - 6*lx22 + 2*x22)))/(3*x22*x22*y4)
      + d3*(1 + x22*(-6 + x22*(18 + x22*(-10 + 12*lx22 - 3*x22))))/(6*x22*x22*x22*y5);
}

double F9(double x1, double x2) noexcept
{
   const double ax1 = std::fabs(x1);
   const double ax2 = std::fabs(x2);

   if (is_equal(ax1, 1., 0.01) && is_equal(ax2, 1., 0.01)) {
      return F9_1_1(x1, x2);
   }

   if (is_equal(ax1, 1., 0.01)) {
      return F9_1_x2(x1, x2);
   }

   if (is_equal(ax2, 1., 0.01)) {
      return F9_1_x2(x2, x1);
   }

   if (is_zero(x1, 0.0001)) {
      return F9_0_x2(x1, x2);
   }

   if (is_zero(x2, 0.0001)) {
      return F9_0_x2(x2, x1);
   }

   if (is_equal(ax1, ax2, 0.00001)) {
      return F9_x1_x2(x1, x2);
   }

   const double x12 = sqr(x1);
   const double x22 = sqr(x2);

   return 2/(x12 - x22)*(
      + x12/(x12 - 1)*std::log(x12)
      - x22/(x22 - 1)*std::log(x22));
}

double f(double r) noexcept
{
   return F5(r);
}

double g(double r) noexcept
{
   return F7(r);
}

double f1(double r) noexcept
{
   const double r2 = sqr(r);

   if (is_zero(r, 0.01)) {
      return 18./7.*r2;
   }

   const double d = r2 - 1;

   if (is_equal(std::fabs(r), 1., 0.01)) {
      return 1 + d*(4/7. + d*(-13/70. + d*(3/35. - 23/490.*d)));
   }

   const double d2 = sqr(d);
   const double d3 = d*d2;

   return 6*(r2 + 3)*r2/(7*d2) + 6*(r2 - 5)*sqr(r2)*std::log(r2)/(7*d3);
}

double f2(double r) noexcept
{
   const double r2 = sqr(r);

   if (is_zero(r, 0.01)) {
      return 22./9.*r2;
   }

   const double d = r2 - 1;

   if (is_equal(std::fabs(r), 1., 0.01)) {
      return 1 + d*(16/27. + d*(-49/270. + d*(11/135. - 83/1890.*d)));
   }

   const double d2 = sqr(d);
   const double d3 = d*d2;

   return 2*(r2 + 11)*r2/(9*d2) + 2*(5*r2 - 17)*sqr(r2)*std::log(r2)/(9*d3);
}

double f3(double r) noexcept
{
   if (is_zero(r, 1e-6)) {
      return 4./3.;
   }

   const double r2 = sqr(r);
   const double d = r2 - 1;

   if (is_equal(std::fabs(r), 1., 0.01)) {
      return 1 + d*(2/9. + d*(1/90. + d*(-2/45. + 29/630.*d)));
   }

   const double r4 = sqr(r2);
   const double d2 = sqr(d);
   const double d3 = d*d2;

   return 2*(r4 + 9*r2 + 2)/(3*d2) + 2*(r4 - 7*r2 - 6)*r2*std::log(r2)/(3*d3);
}

double f4(double r) noexcept
{
   if (is_zero(r, 1e-6)) {
      return 12./7.;
   }

   const double r2 = sqr(r);
   const double d = r2 - 1;

   if (is_equal(std::fabs(r), 1., 0.01)) {
      return 1 + d*(2/21. + d*(13/210. + d*(-8/105. + 101/1470.*d)));
   }

   const double r4 = sqr(r2);
   const double d2 = sqr(d);
   const double d3 = d*d2;

   return 2*(5*r4 + 25*r2 + 6)/(7*d2) + 2*(r4 - 19*r2 - 18)*r2*std::log(r2)/(7*d3);
}

/// f5(r1,r2) in the limit r1 -> 0, r2 -> 0
static double f5_0_0(double r1, double r2) noexcept
{
   if (is_zero(r1, 1e-6) && is_zero(r2, 1e-6)) {
      return 0.75;
   }

   if (std::abs(r1) > std::abs(r2)) {
      std::swap(r1, r2);
   }

   const double r12 = sqr(r1);
   const double r22 = sqr(r2);
   const double r2lr2 = xlogx(r2);
   const double lr2 = logx(r2);

   // expansion of f5 in
   //   r1 up to (including) r1^2
   //   r2 up to (including) r1^3
   return 0.75*(
      + 1
      + 2*r1*(r2 + r2lr2)
      + 2*r22 + 2*r2*r2lr2
      + r12*(2 + 2*lr2 + r2*(2*r2 + 6*r2lr2))
      + r1*r22*(6*r2lr2 + 2*r2)
   );
}

/// f5(r1,r2) in the limit r1 -> 1 and r2 -> 1
static double f5_1_1(double r1, double r2) noexcept
{
   const double d1 = r1 - 1;
   const double d2 = r2 - 1;
   const double d12 = sqr(d1);
   const double d22 = sqr(d2);
   const double d13 = d1*d12;
   const double d23 = d2*d22;
   const double d14 = sqr(d12);
   const double d24 = sqr(d22);

   return 1
      + 3./8.*(d1 + d2)
      + 1./40.*(d12 + d1*d2 + d22)
      + 1./20.*(-d13 - d12*d2 - d1*d22 - d23)
      + 1./280.*(9*(d14 + d13*d2 + d12*d22 + d1*d23 + d24));
}

/// f5(r1,r2) in the limit r1 -> 1, r2 != 1
static double f5_1_r2(double r1, double r2) noexcept
{
   if (is_zero(r2, 0.01)) {
      const double d = r1 - 1;
      return 0.75*(1 + d*(1./3. + d/6. + r2/3.));
   }

   const double d1 = r1 - 1;
   const double d2 = r2 - 1;
   const double d12 = sqr(d1);
   const double d13 = d12*d1;
   const double d14 = d13*d1;
   const double d22 = sqr(d2);
   const double d23 = d22*d2;
   const double d24 = d23*d2;
   const double d25 = d24*d2;
   const double d26 = d25*d2;
   const double d27 = d26*d2;
   const double y = 1 + r2;
   const double y2 = sqr(y);
   const double r22 = sqr(r2);
   const double lr22 = std::log(r22);

   return
      + (-3 + r2*(3 + r2*(6 + r2*(3*lr22 + r2*(-3 + (-3 + 3*lr22)*r2)))))/(4.*d23*y2)
      + d1*(1 + r2*(-1 + r2*(-2 + r2*(8 + 3*lr22 + r2*(1 + (-7 + 3*lr22)*r2)))))/(4.*d24*y2)
      + d12*(-1 + r2*(4 + r2*(-1 + r2*(4 + 6*lr22 + r2*(5 + (-8 + 6*lr22 - 3*r2)*r2)))))/(8.*d25*y2)
      + d13*(-4 + r2*(17 + r2*(-4 + r2*(25 + 30*lr22 + r2*(20 + r2*(-41 + 30*lr22 + (-12 - r2)*r2))))))/(40.*d26*y2)
      + d14*(-2 + r2*(9 + r2*(4 + r2*(33 + 30*lr22 + r22*(-33 + 30*lr22 + r2*(-4 + r2*(-9 + 2*r2)))))))/(40.*d27*y2);
}

/// f5(r1,r2) in the limit r1 -> 0
static double f5_0_r2(double r1, double r2) noexcept
{
   const double r22 = sqr(r2);
   const double lr22 = std::log(r22);

   // expansion in terms of r1 around r1 = 0 up to maximum possible
   // power such that no IR-divergent terms appear (e.g. log(r1))
   const double res =
      (1 + r22*(lr22 + (-1 + lr22)*r22)
       + r1*(r1*(2 + lr22 + (-2 + lr22)*r22) + r2*(2 + lr22 + (-2 + lr22)*r22)))/sqr(-1 + r22);

   return 0.75 * res;
}

/// f5(r1,r2) in the limit r1 -> r2
static double f5_r1_r2(double r1, double r2) noexcept
{
   const double d = r1 - r2;
   const double d2 = sqr(d);
   const double d3 = d2*d;
   const double r22 = sqr(r2);
   const double y = r22 - 1;
   const double y3 = sqr(y)*y;
   const double y4 = y3*y;
   const double y5 = y4*y;
   const double y6 = y5*y;
   const double lr22 = std::log(r22);

   return 0.75*(
      + (-1 + r22*(-5 - 3*lr22 + r22*(5 - 6*lr22 + (1 + lr22)*r22)))/y3
      + d*r2*(11 + 3*lr22 + r22*(3 + 18*lr22 + r22*(-15 + 3*lr22 + r22)))/y4
      + d2*(-17 - 3*lr22 + r22*(-116 - 75*lr22 + r22*(90 - 105*lr22 + r22*(44 - 9*lr22 - r22))))/(3.*y5)
      + d3*(3 + r22*(273 + 90*lr22 + r22*(314 + 510*lr22 + r22*(-498 + 342*lr22 + r22*(-93 + 18*lr22 + r22)))))/(6.*r2*y6)
   );
}

double f5(double r1, double r2) noexcept
{
   const double eps_zero = 1e-4;
   const double eps_one = 1e-2;

   if (is_zero(r1, eps_zero) && is_zero(r2, eps_zero)) {
      return f5_0_0(r1, r2);
   }

   if (is_equal(r1, 1., eps_one) && is_equal(r2, 1., eps_one)) {
      return f5_1_1(r1, r2);
   }

   if (is_equal(r1, -1., eps_one) && is_equal(r2, -1., eps_one)) {
      return f5_1_1(-r1, -r2);
   }

   if (is_equal(r1, 1., eps_one)) {
      return f5_1_r2(r1, r2);
   }

   if (is_equal(r2, 1., 0.01)) {
      return f5_1_r2(r2, r1);
   }

   if (is_zero(r1, eps_zero)) {
      return f5_0_r2(r1, r2);
   }

   if (is_zero(r2, eps_zero)) {
      return f5_0_r2(r2, r1);
   }

   if (is_equal(r1, r2, 0.0001)) {
      return f5_r1_r2(r2, r1);
   }

   const double r12 = sqr(r1);
   const double r22 = sqr(r2);

   const double result
      = (1 + sqr(r1 + r2) - r12*r22)/((r12 - 1)*(r22 - 1))
      + (cube(r1)*(r12 + 1)*std::log(r12))/(sqr(r12 - 1)*(r1 - r2))
      - (cube(r2)*(r22 + 1)*std::log(r22))/((r1 - r2)*sqr(r22 - 1));

   return 0.75 * result;
}

/// f6(r1,r2) in the limit r1 -> 0 and r2 -> 0
static double f6_0_0(double r1, double r2) noexcept
{
   if (is_zero(r1, 1e-10) && is_zero(r2, 1e-10)) {
      return 0.;
   }

   const double r22 = sqr(r2);
   const double lr22 = std::log(r22);

   const double res =
      r22*(1 + (1 + lr22)*r22)
      + r1*(r2*(1 + (1 + lr22)*r22)
      + r1*(1 + r22*(1 + lr22 + (1 + 2*lr22)*r22)
      + r1*(r2*(1 + lr22 + (1 + 2*lr22)*r22)
      + r1*(1 + lr22 + r22*(1 + 2*lr22 + (1 + 3*lr22)*r22)))));

   return 6./7.*res;
}

/// f6(r1,r2) in the limit r1 -> 1 and r2 -> 1
static double f6_1_1(double r1, double r2) noexcept
{
   const double d1 = r1 - 1;
   const double d2 = r2 - 1;

   return
      1 + d2*(4./7 + d2*(-2./35 + (-1./70 + 9.*d2/490)*d2)) +
      d1*(4./7 + d2*(-2./35 + d2*(-1./70 + (9./490 - 3.*d2/245)*d2)) +
      d1*(-2./35 + d2*(-1./70 + d2*(9./490 + (-3./245 + d2/147)*d2)) +
      d1*(-1./70 + d2*(9./490 + d2*(-3./245 + (1./147 - d2/294)*d2)) +
      d1*(9./490 + d2*(-3./245 + d2*(1./147 + (-1./294 + 5.*d2/3234)*d2))))));
}

/// f6(r1,r2) in the limit r1 -> 0 and r2 -> 1
static double f6_0_1(double r1, double r2) noexcept
{
   const double d2 = r2 - 1;
   const double d22 = sqr(d2);

   return
      + 3./7 + d2*(4./7 + (-2./35 + 3.*d2/70)*d22)
      + r1*(3./7 + d2*(1./7 + d2*(-1./7 + (3./35 - 3.*d2/70)*d2))
      + r1*(3./7 + d2*(-2./7 + d2*(1./7 + (-2./35 + d2/70)*d2))
      + r1*(-3./7 + d2*(1./7 + (-2./35 + d2/14)*d22)
      + r1*(-3./7 + d2*(4./7 + d2*(-4./7 + (18./35 - 31.*d2/70)*d2))))));
}

/// f6(r1,r2) in the limit r1 -> 1
static double f6_1_r2(double r1, double r2) noexcept
{
   if (is_zero(r2, 1e-4)) {
      return f6_0_1(r2, r1);
   }

   const double r22 = sqr(r2);
   const double lr22 = std::log(r22);
   const double y = 1 + r2;
   const double y2 = sqr(y);
   const double z = r2 - 1;
   const double z2 = sqr(z);
   const double z3 = z2*z;
   const double z4 = z3*z;
   const double z5 = z4*z;
   const double z6 = z5*z;
   const double z7 = z6*z;
   const double d = r1 - 1;
   const double d2 = sqr(d);
   const double d3 = d2*d;
   const double d4 = d3*d;

   return
      (-3 + r22*(6 + r2*(6 + r2*(-3 + r2*(-6 + 6*lr22)))))/(7.*z3*y2)
      + d*(4 + r2*(-7 + r2*(-8 + r2*(20 + r2*(4 + r2*(-13 + 6*lr22))))))/(7.*z4*y2)
      + d2*r2*(1 + r2*(-4 + r2*(4 + r2*(8 + r2*(-5 - 4*r2 + 6*lr22)))))/(7.*z5*y2)
      + d3*(-2 + r2*(11 + r2*(-22 + r2*(10 + r2*(50 + r2*(-23 + r2*(-26 + 2*r2) + 30*lr22))))))/(35.*z6*y2)
      + d4*(-3 + r2*(18 + r2*(-40 + r2*(24 + r2*(90 + r2*(-42 + r2*(-48 + r22) + 60*lr22))))))/(70.*z7*y2);
}

/// f6(r1,r2) in the limit r1 -> 0
static double f6_0_r2(double r1, double r2) noexcept
{
   const double r22 = sqr(r2);
   const double lr22 = std::log(r22);

   return
      6./7.*(r22*(1 + (-1 + lr22)*r22)
             + r1*(r2*(1 + (-1 + lr22)*r22)
             + r1*(1 + (-1 + lr22)*r22
             + r1*(r1*(1 + lr22 - r22) + r2*(1 + lr22 - r22)))))/sqr(-1 + r22);
}

// f6(r1,r2) in the limit r1 -> r2
static double f6_r1_r2(double r1, double r2) noexcept
{
   const double d = r1 - r2;
   const double d2 = sqr(d);
   const double d3 = d2*d;
   const double d4 = d3*d;
   const double r22 = sqr(r2);
   const double y = r22 - 1;
   const double y2 = sqr(y);
   const double y3 = y2*y;
   const double y4 = y3*y;
   const double y5 = y4*y;
   const double y6 = y5*y;
   const double y7 = y6*y;
   const double lr22 = std::log(r22);

   return
      6./7. * (
          r22*(-3 + r22*(2 - 5*lr22 + (1 + lr22)*r22))/y3
          + d*r2*(3 + r22*(7 + 10*lr22 + r22*(-11 + 2*lr22 + r22)))/y4
          + d2*(-3 + r22*(-62 - 30*lr22 + r22*(36 - 60*lr22 + r22*(30 - 6*lr22 - r22))))/(3.*y5)
          + d3*r2*(107 + 30*lr22 + r22*(206 + 240*lr22 + r22*(-252 + 198*lr22 + r22*(-62 + 12*lr22 + r22))))/(6.*y6)
          + d4*(-167 - 30*lr22 + r22*(-2215 - 1050*lr22 + r22*(-510 - 3150*lr22 + r22*(2570 - 1470*lr22 + r22*(325 - 60*lr22 - 3*r22)))))/(30.*y7)
      );
}

double f6(double r1, double r2) noexcept
{
   const double eps_zero = 1e-4;

   if (is_zero(r1, eps_zero) && is_zero(r2, eps_zero)) {
      return f6_0_0(r1, r2);
   }

   if (is_equal(r1, 1., 0.02) && is_equal(r2, 1., 0.02)) {
      return f6_1_1(r1, r2);
   }

   if (is_equal(r1, -1., 0.02) && is_equal(r2, -1., 0.02)) {
      return f6_1_1(-r1, -r2);
   }

   if (is_equal(r1, 1., 0.001)) {
      return f6_1_r2(r1, r2);
   }

   if (is_equal(r2, 1., 0.001)) {
      return f6_1_r2(r2, r1);
   }

   if (is_zero(r1, eps_zero)) {
      return f6_0_r2(r1, r2);
   }

   if (is_zero(r2, eps_zero)) {
      return f6_0_r2(r2, r1);
   }

   if (is_equal(r1, r2, 1e-4)) {
      return f6_r1_r2(r2, r1);
   }

   const double r12 = sqr(r1);
   const double r22 = sqr(r2);

   const double result
      = (r12 + r22 + r1*r2 - r12*r22)/((r12 - 1)*(r22 - 1))
      + (pow5(r1)*std::log(r12))/(sqr(r12 - 1)*(r1 - r2))
      - (pow5(r2)*std::log(r22))/((r1 - r2)*sqr(r22 - 1));

   return 6./7. * result;
}

/// f7(r1,r2) in the limit r1 -> 0, r2 -> 0
static double f7_0_0(double r1, double r2) noexcept
{
   if (is_zero(r1, 1e-6) && is_zero(r2, 1e-6)) {
      return 6.;
   }

   if (std::abs(r1) > std::abs(r2)) {
      std::swap(r1, r2);
   }

   const double r22 = sqr(r2);
   const double lr22 = std::log(r22);

   const double res =
      1 + (1 + lr22)*r22
      + r1*((1 + lr22)*r2 + r1*(1 + lr22 + (1 + 2*lr22)*r22));

   return 6*res;
}

/// f7(r1,r2) in the limit r1 -> 1 and r2 -> 1
static double f7_1_1(double r1, double r2) noexcept
{
   const double d1 = r1 - 1;
   const double d2 = r2 - 1;

   return
      1 + d2*(-1 + d2*(3./5 + (-3./10 + 9.*d2/70)*d2))
      + d1*(-1 + d2*(3./5 + d2*(-3./10 + (9./70 - 3.*d2/70)*d2))
      + d1*(3./5 + d2*(-3./10 + d2*(9./70 + (-3./70 + d2/210)*d2))
      + d1*(-3./10 + d2*(9./70 + d2*(-3./70 + (1./210 + d2/105)*d2))
      + d1*(9./70 + d2*(-3./70 + d2*(1./210 + (1./105 - d2/77)*d2))))));
}

/// f7(r1,r2) in the limit r1 -> 0 and r2 -> 1
static double f7_0_1(double r1, double r2) noexcept
{
   const double d = r2 - 1;
   const double r12 = sqr(r1);

   return 3 + r1*(-4 - 3*r1 + r2)
      + d*(-2 + 4*r12
      + d*(1 - 4*r12
      + d*(-0.4 + r1*(-0.4 + (18*r1)/5.)
      + d*(0.1 + (0.5 - (31*r1)/10.)*r1))));
}

/// f7(r1,r2) in the limit r1 -> 1
static double f7_1_r2(double r1, double r2) noexcept
{
   if (is_zero(r2, 0.0001)) {
      return f7_0_1(r2, r1);
   }

   const double d = r1 - 1;
   const double d2 = sqr(d);
   const double d3 = d2*d;
   const double d4 = d3*d;
   const double y = r2 - 1;
   const double y2 = sqr(y);
   const double y3 = y2*y;
   const double y4 = y3*y;
   const double y5 = y4*y;
   const double y6 = y5*y;
   const double y7 = y6*y;
   const double z = 1 + r2;
   const double z2 = sqr(z);
   const double r22 = sqr(r2);
   const double lr22 = std::log(r22);

   return
      (-3 + r2*(6 + r2*(6 + (-6 + 6*lr22 - 3*r2)*r2)))/(y3*z2)
      + d*(-2 + r2*(5 + r2*(4 + r2*(-4 + 6*lr22 + (-2 - r2)*r2))))/(y4*z2)
      + d2*(-1 + r2*(3 + r2*(3 + r2*(6*lr22 + r2*(-3 + (-3 + r2)*r2)))))/(y5*z2)
      + d3*(-2 + r2*(6 + r2*(18 + r2*(15 + 30*lr22 + r2*(-30 + r2*(-18 + (14 - 3*r2)*r2))))))/(5.*y6*z2)
      + d4*(-1 + r22*(48 + r2*(42 + 60*lr22 + r2*(-90 + r2*(-24 + r2*(40 + r2*(-18 + 3*r2)))))))/(10.*y7*z2);
}

/// f7(r1,r2) in the limit r1 -> 0
static double f7_0_r2(double r1, double r2) noexcept
{
   const double r12 = sqr(r1);
   const double r22 = sqr(r2);
   const double lr22 = std::log(r22);
   const double r12lr12 = xlogx(r12);

   return
      6*(
         + r22*(1 + (-1 + lr22)*r22)
         + r1*(r2*(-r12lr12 + r22*(1 + lr22 + 2*r12lr12 + (-1 - r12lr12)*r22))
         + r1*(-r12lr12 + r22*(1 + lr22 + 2*r12lr12 + (-1 - r12lr12)*r22)
         + r1*(r1*(lr22 + r22*(1 - r22)) + r2*(lr22 + r22*(1 - r22)))))
         )/(r22*sqr(-1 + r22));
}

/// f7(r1,r2) in the limit r1 -> r2
static double f7_r1_r2(double r1, double r2) noexcept
{
   const double d = r1 - r2;
   const double d2 = sqr(d);
   const double d3 = d2*d;
   const double d4 = d3*d;
   const double r22 = sqr(r2);
   const double y = r22 - 1;
   const double y2 = sqr(y);
   const double y3 = y2*y;
   const double y4 = y3*y;
   const double y5 = y4*y;
   const double y6 = y5*y;
   const double y7 = y6*y;
   const double lr22 = std::log(r22);

   const double res =
      (-1 + r22*(-2 - 3*lr22 + (3 - lr22)*r22))/y3
      + d*r2*(8 + 3*lr22 + r22*(-4 + 8*lr22 + (-4 + lr22)*r22))/y4
      + d2*(-14 - 3*lr22 + r22*(-54 - 45*lr22 + r22*(54 - 45*lr22 + (14 - 3*lr22)*r22)))/(3.*y5)
      + d3*(3 + r22*(166 + 60*lr22 + r22*(108 + 270*lr22 + r22*(-246 + 144*lr22 + (-31 + 6*lr22)*r22))))/(6.*r2*y6)
      + d4*(3 + r22*(-325 - 60*lr22 + r22*(-2570 - 1470*lr22 + r22*(510 - 3150*lr22 + r22*(2215 - 1050*lr22 + (167 - 30*lr22)*r22)))))/(30.*r22*y7);

   return 6*res;
}

double f7(double r1, double r2) noexcept
{
   const double eps_zero = 1e-4;
   const double eps_one = 2e-2;

   if (is_zero(r1, eps_zero) && is_zero(r2, eps_zero)) {
      return f7_0_0(r1, r2);
   }

   if (is_equal(r1, 1., eps_one) && is_equal(r2, 1., eps_one)) {
      return f7_1_1(r1, r2);
   }

   if (is_equal(r1, -1., eps_one) && is_equal(r2, -1., eps_one)) {
      return f7_1_1(-r1, -r2);
   }

   if (is_equal(r1, 1., eps_one)) {
      return f7_1_r2(r1, r2);
   }

   if (is_equal(r2, 1., eps_one)) {
      return f7_1_r2(r2, r1);
   }

   if (is_zero(r1, eps_zero)) {
      return f7_0_r2(r1, r2);
   }

   if (is_zero(r2, eps_zero)) {
      return f7_0_r2(r2, r1);
   }

   if (is_equal(r1, r2, 0.0001)) {
      return f7_r1_r2(r2, r1);
   }

   const double r12 = sqr(r1);
   const double r22 = sqr(r2);

   const double result
      = (1 + r1*r2)/((r12 - 1)*(r22 - 1))
      + (cube(r1)*std::log(r12))/(sqr(r12 - 1)*(r1 - r2))
      - (cube(r2)*std::log(r22))/((r1 - r2)*sqr(r22 - 1));

   return 6. * result;
}

/// f8(r1,r2) in the limit r1 -> 0 and r2 -> 0
static double f8_0_0(double r1, double r2) noexcept
{
   if (is_zero(r1, 1e-10) && is_zero(r2, 1e-10)) {
      return 0.;
   }

   if (std::abs(r1) > std::abs(r2)) {
      std::swap(r1, r2);
   }

   const double r22 = sqr(r2);
   const double lr22 = logx(r22);
   const double r2lr22 = r2*lr22;
   const double r22lr22 = r2*r2lr22;

   if (is_zero(r1, 1e-8) && is_zero(r2, 1e-8)) {
      return
         (3*r2)/2.
         + r1*(1.5 + (3*r22)/2. + (3*r22lr22)/2.
         + r1*((3*r2)/2. + (3*r2lr22)/2.));
   }

   return
      r2*(1.5 + (3*r22)/2. + (3*r22lr22)/2.)
      + r1*(1.5 + (3*r22)/2. + (3*r22lr22)/2.
      + r1*(r2*(1.5 + (3*r22)/2. + 3*r22lr22)
      + r1*(1.5 + (3*lr22)/2. + (3*r22)/2. + 3*r22lr22) + (3*r2lr22)/2.));
}

/// f8(r1,r2) in the limit r1 -> 1 and r2 -> 1
static double f8_1_1(double r1, double r2) noexcept
{
   const double d1 = r1 - 1;
   const double d2 = r2 - 1;
   const double d22 = sqr(d2);

   return
      1 + d22*(-1./10 + (3./40 - 3.*d2/70)*d2)
      + d1*(d2*(-1./10 + d2*(3./40 + (-3./70 + 3.*d2/140)*d2))
      + d1*(-1./10 + d2*(3./40 + d2*(-3./70 + (3./140 - d2/105)*d2))
      + d1*(3./40 + d2*(-3./70 + d2*(3./140 + (-1./105 + d2/280)*d2))
      + d1*(-3./70 + d2*(3./140 + d2*(-1./105 + (1./280 - d2/1155)*d2))))));
}

/// f8(r1,r2) in the limit r1 -> 0 and r2 -> 1
static double f8_0_1(double r1, double r2) noexcept
{
   const double r12 = sqr(r1);
   const double d2 = r2 - 1;

   return 0.75 + r1*(0.75 + r1*(-0.75 + r1*(-1.75 + r2)))
      + d2*(0.25 + (-0.5 + r1/4.)*r1
      + d2*(-0.25 + r1*(0.25 - r12)
      + d2*(0.15 + r1*(-0.1 + (-0.1 + (9*r1)/10.)*r1))));
}

/// f8(r1,r2) in the limit r1 -> 1
static double f8_1_r2(double r1, double r2) noexcept
{
   if (is_zero(r2, 1e-4)) {
      return f8_0_1(r2, r1);
   }

   const double d1 = r1 - 1;
   const double d12 = sqr(d1);
   const double d13 = d12*d1;
   const double d14 = d13*d1;
   const double y = 1 + r2;
   const double y2 = sqr(y);
   const double d2 = r2 - 1;
   const double d22 = sqr(d2);
   const double d23 = d22*d2;
   const double d24 = d23*d2;
   const double d25 = d24*d2;
   const double d26 = d25*d2;
   const double d27 = d26*d2;
   const double r22 = sqr(r2);
   const double lr22 = std::log(r22);

   return
      (-3 + r22*(12 + (-9 + 6*lr22)*r22))/(4.*d23*y2)
      + d1*(1 + r2*(-4 + r2*(4 + r2*(8 + (-5 + 6*lr22 - 4*r2)*r2))))/(4.*d24*y2)
      + d12*(1 + r2*(-4 + r2*(4 + r2*(8 + (-5 + 6*lr22 - 4*r2)*r2))))/(4.*d25*y2)
      + d13*(3 + r2*(-14 + r2*(18 + r2*(30 + r2*(-15 + 30*lr22 + r2*(-18 + r2*(-6 + 2*r2)))))))/(20.*d26*y2)
      + d14*(3 + r2*(-16 + r2*(24 + r2*(48 + r2*(60*lr22 + r2*(-48 + r2*(-24 + (16 - 3*r2)*r2)))))))/(40.*d27*y2);
}

/// f8(r1,r2) in the limit r1 -> 0
static double f8_0_r2(double r1, double r2) noexcept
{
   const double r12 = sqr(r1);
   const double r22 = sqr(r2);
   const double lr22 = std::log(r22);
   const double r12lr12 = xlogx(r12);

   return
      (r22*(3 + (-3 + 3*lr22)*r22)
      + r1*(r2*(3 + (-3 + 3*lr22)*r22)
      + r1*(-3*r12lr12 + r22*(3 + 3*lr22 + 6*r12lr12 + (-3 - 3*r12lr12)*r22)
      + r1*(r2*(3 + 3*lr22 - 3*r22) + r1*(3*lr22 + r22*(3 - 3*r22))))))/(2*r2*sqr(-1 + r22));
}

/// f8(r1,r2) in the limit r1 -> r2
static double f8_r1_r2(double r1, double r2) noexcept
{
   const double d = r1 - r2;
   const double d2 = sqr(d);
   const double d3 = d2*d;
   const double d4 = d3*d;
   const double r22 = sqr(r2);
   const double y = r22 - 1;
   const double y2 = sqr(y);
   const double y3 = y2*y;
   const double y4 = y3*y;
   const double y5 = y4*y;
   const double y6 = y5*y;
   const double y7 = y6*y;
   const double lr22 = std::log(r22);

   return
      (r2*(-3 + r22*(-6*lr22 + 3*r22)))/y3
      + d*(3 + r22*(27 + 18*lr22 + r22*(-27 + 18*lr22 - 3*r22)))/(2.*y4)
      + d2*r2*(-19 - 6*lr22 + r22*(-9 - 30*lr22 + r22*(27 - 12*lr22 + r22)))/y5
      + d3*(31 + 6*lr22 + r22*(246 + 144*lr22 + r22*(-108 + 270*lr22 + r22*(-166 + 60*lr22 - 3*r22))))/(4.*y6)
      + d4*(-3 + r22*(-285 - 90*lr22 + r22*(-570 - 630*lr22 + r22*(570 - 630*lr22 + r22*(285 - 90*lr22 + 3*r22)))))/(5.*r2*y7);
}

double f8(double r1, double r2) noexcept
{
   const double eps_zero = 5e-5;
   const double eps_one = 1e-2;

   if (is_zero(r1, eps_zero) && is_zero(r2, eps_zero)) {
      return f8_0_0(r1, r2);
   }

   if (is_equal(r1, 1., eps_one) && is_equal(r2, 1., eps_one)) {
      return f8_1_1(r1, r2);
   }

   if (is_equal(r1, -1., eps_one) && is_equal(r2, -1., eps_one)) {
      return -1.;
   }

   if (is_equal(r1, 1., eps_one)) {
      return f8_1_r2(r1, r2);
   }

   if (is_equal(r2, 1., eps_one)) {
      return f8_1_r2(r2, r1);
   }

   if (is_zero(r1, eps_zero)) {
      return f8_0_r2(r1, r2);
   }

   if (is_zero(r2, eps_zero)) {
      return f8_0_r2(r2, r1);
   }

   if (is_equal(r1, r2, 0.0001)) {
      return f8_r1_r2(r2, r1);
   }

   const double r12 = sqr(r1);
   const double r22 = sqr(r2);

   const double result
      = (r1 + r2)/((r12 - 1)*(r22 - 1))
      + (quad(r1)*std::log(r12))/(sqr(r12 - 1)*(r1 - r2))
      - (quad(r2)*std::log(r22))/((r1 - r2)*sqr(r22 - 1));

   return 1.5 * result;
}

double fth1(double y) noexcept
{
   const double eps = 10.*std::numeric_limits<double>::epsilon();

   if (is_zero(y, eps)) {
      return 0.;
   }

   if (is_equal(std::abs(y), 1.0, eps)) {
      return -1.;
   }

   if (!is_zero(y, eps) && !is_equal(std::abs(y), 1.0, eps)) {
      const double y2 = sqr(y);

      return y2*std::log(y2) / (1. - y2);
   }

   return 0.;
}

double fth2(double y) noexcept
{
   const double eps = 10.*std::numeric_limits<double>::epsilon();

   if (is_zero(y, eps)) {
      return 0.;
   }

   if (is_equal(std::abs(y), 1.0, eps)) {
      return 0.5;
   }

   if (!is_zero(y, eps) && !is_equal(std::abs(y), 1.0, eps)) {
      const double y2 = sqr(y);

      return (1. + y2*std::log(y2) / (1. - y2)) / (1 - y2);
   }

   return 0.;
}

double fth3(double y) noexcept
{
   const double eps = 10.*std::numeric_limits<double>::epsilon();
   const double z2 = 1.644934066848226;

   if (is_zero(y, eps)) {
      return z2;
   }

   if (is_equal(std::abs(y), 1.0, eps)) {
      return -9./4.;
   }

   if (!is_zero(y, eps) && !is_equal(std::abs(y), 1.0, eps)) {
      const double y2 = sqr(y);
      const double y4 = sqr(y2);
      const double ly = std::log(y2);

      return (-1. + 2*y2 + 2*y4)
         *(-z2 - y2*ly + ly*std::log(1. - y2) + dilog(y2))
         / sqr(1 - y2);
   }

   return 0.;
}

/// First derivative of F1
double D1F1(double x) noexcept
{
   if (is_equal(x, 1., 0.0001)) {
      return (5 - 8*x + 3*sqr(x))/6.;
   }

   return (2*(-1 + sqr(x)) - std::log(sqr(x))*(1 + sqr(x)))/sqr(-1 + sqr(x));
}

/// First derivative of F2
double D1F2(double x) noexcept
{
   if (is_equal(x, 1., 0.001)) {
      return (133 - 326*x - 138*cube(x) + 25*quad(x) + 306*sqr(x))/35.;
   }

   return (-12*x*(3 - 3*quad(x) + std::log(sqr(x))*(1 + quad(x) + 4*sqr(x)))) /
      quad(-1 + sqr(x));
}

/// First derivative of F3
double D1F3(double x) noexcept
{
   if (is_equal(x, 1., 0.0001)) {
      return (1541 - 2048*x - 256*cube(x) + 15*quad(x) + 1098*sqr(x))/630.;
   }

   return (-2*(7 - 13*quad(x) + 6*sqr(x) + std::log(sqr(x))*(1 + 4*quad(x)
      + 15*sqr(x))))/(3.*cube(-1 + sqr(x)));
}

/// First derivative of F4
double D1F4(double x) noexcept
{
   if (is_equal(x, 1., 0.0001)) {
      return -0.3333333333333333 - (2*cube(-1 + x))/3. + (11*quad(-1 + x))/14.
         + (2*sqr(-1 + x))/5.;
   }

   return (-2*(-3 + quad(x) + 2*sqr(x)) + std::log(sqr(x))*(2 + 6*sqr(x))) /
         cube(-1 + sqr(x));
}

/// First derivative of F5
double D1F5(double x) noexcept
{
   if (is_equal(x, 1., 0.001)) {
      return (3*(70 - 176*x - 80*cube(x) + 15*quad(x) + 171*sqr(x)))/70.;
   }

   return (-3*(-1 + pow6(x) + 9*quad(x) - 9*sqr(x) - 6*std::log(sqr(x))*(quad(x)
      + sqr(x))))/quad(-1 + sqr(x));
}

/// First derivative of F6
double D1F6(double x) noexcept
{
   if (is_equal(x, 1., 0.0001)) {
      return (204 - 11*x + 87*cube(x) - 20*quad(x) - 120*sqr(x))/210.;
   }

   return (x*(3 + 2*std::log(sqr(x)) + quad(x) - 4*sqr(x)))/cube(-1 + sqr(x));
}

/// First derivative of F7
double D1F7(double x) noexcept
{
   if (is_equal(x, 1., 0.001)) {
      return (-3*(-14 - 80*x - 51*cube(x) + 10*quad(x) + 100*sqr(x)))/35.;
   }

   return (6*x*(2 + pow6(x) - 6*quad(x) + 3*sqr(x)
                + 6*std::log(sqr(x))*sqr(x)))/quad(-1 + sqr(x));
}

/// First derivative of f
double D1f(double x) noexcept
{
   if (is_equal(x, 1., 0.001)) {
      return (3*(70 - 176*x - 80*cube(x) + 15*quad(x) + 171*sqr(x)))/70.;
   }

   return (-3*(-1 + pow6(x) + 9*quad(x) - 9*sqr(x) - 6*std::log(sqr(x))*(quad(x)
      + sqr(x))))/quad(-1 + sqr(x));
}

/// First derivative of g
double D1g(double x) noexcept
{
   if (is_equal(x, 1., 0.001)) {
      return (-3*(-14 - 80*x - 51*cube(x) + 10*quad(x) + 100*sqr(x)))/35.;
   }

   return (6*x*(2 + pow6(x) - 6*quad(x) + 3*sqr(x)
      + 6*std::log(sqr(x))*sqr(x)))/quad(-1 + sqr(x));
}

/// First derivative of f1
double D1f1(double x) noexcept
{
   if (is_equal(x, 1., 0.01)) {
      return (-2*(-71 - 315*x - 225*cube(x) + 45*quad(x) + 426*sqr(x)))/245.;
   }

   return (12*x*(3 + pow6(x) - 11*quad(x) + 7*sqr(x)
      + 2*std::log(sqr(x))*sqr(x)*(5 + sqr(x))))/(7.*quad(-1 + sqr(x)));
}

/// First derivative of f2
double D1f2(double x) noexcept
{
   if (is_equal(x, 1., 0.01)) {
      return (-2*(-239 - 1275*x - 837*cube(x) + 165*quad(x) + 1626*sqr(x)))/945.;
   }

   return (4*x*(11 + 5*pow6(x) - 35*quad(x) + 19*sqr(x)
      + 2*std::log(sqr(x))*sqr(x)*(17 + sqr(x))))/(9.*quad(-1 + sqr(x)));
}

/// First derivative of f3
double D1f3(double x) noexcept
{
   if (is_equal(x, 1., 0.001)) {
      return (-2*(386 - 1143*x - 495*cube(x) + 90*quad(x) + 1092*sqr(x)))/315.;
   }

   return (4*x*(19 + pow6(x) - 19*quad(x) - sqr(x)
      + std::log(sqr(x))*(6 + 4*quad(x) + 26*sqr(x))))/(3.*quad(-1 + sqr(x)));
}

/// First derivative of f4
double D1f4(double x) noexcept
{
   if (is_equal(x, 1., 0.001)) {
      return (-2*(1184 - 3099*x - 1323*cube(x) + 240*quad(x) + 2928*sqr(x)))/735.;
   }

   return (4*x*(55 + pow6(x) - 55*quad(x) - sqr(x)
      + 2*std::log(sqr(x))*(9 + 8*quad(x) + 37*sqr(x))))/(7.*quad(-1 + sqr(x)));
}

/// First derivative of f5 w.r.t. 1st argument
double D10f5(double x, double y) noexcept
{
   if (is_equal(x, 1., 0.001) && is_equal(y, 1., 0.001))
      return (-117 + 306*y - 91*sqr(y) - 3*sqr(x)*(55 - 36*y + 9*sqr(y))
         + x*(450 - 344*y + 90*sqr(y)))/560.;

   if (is_equal(x, 1., 0.001))
      return (30*cube(y)*std::log(sqr(y))*(1 + sqr(y))*(6 + 2*x*(-4 + y) - 4*y +
         3*sqr(x) + sqr(y)) - (-1 + sqr(y))*(-12 + 71*y + 306*cube(y) +
         43*pow5(y) - 164*quad(y) + 3*sqr(x)*(-4 + 17*y + 42*cube(y) + pow5(y)
         + 12*quad(y) - 8*sqr(y)) - 64*sqr(y) + 2*x*(17 - 76*y - 176*cube(y) +
         12*pow5(y) - 11*quad(y) + 54*sqr(y))))/(40.*pow6(-1 + y)*sqr(1 + y));

   if (is_equal(y, 1., 0.001))
      return (-6*std::log(sqr(x))*sqr(x)*(y*pow5(x) + 3*cube(x)*(-5 + sqr(y)) + 3*(3 -
         3*y + sqr(y)) + quad(x)*(6 - 3*y + 2*sqr(y)) + x*(-3 - 5*y +
         3*sqr(y)) + sqr(x)*(23 - 24*y + 9*sqr(y))) + (-1 + sqr(x))*(3 - 4*y +
         12*pow6(x) + sqr(y) + x*(1 - 4*y + 3*sqr(y)) + pow5(x)*(-35 + 20*y +
         3*sqr(y)) + 2*sqr(x)*(69 - 68*y + 23*sqr(y)) + cube(x)*(-74 - 40*y +
         30*sqr(y)) + quad(x)*(75 - 76*y + 37*sqr(y))))/(8.*cube(1 +
         x)*pow6(-1 + x));

   if (is_equal(x, y, 0.001))
      return (6*y*std::log(sqr(y))*(y + 9*cube(y) + 205*pow5(y) + 247*pow7(y) +
         18*pow9(y) - 2*x*(-1 + 203*pow6(y) + 12*pow8(y) + 245*quad(y) +
         21*sqr(y)) + 3*y*sqr(x)*(15 + 3*pow6(y) + 57*quad(y) + 85*sqr(y))) +
         (-1 + sqr(y))*(3*sqr(x)*(-3 - 92*pow6(y) + pow8(y) - 590*quad(y) -
         276*sqr(y)) + (-7 - 548*pow6(y) + 13*pow8(y) - 2022*quad(y) -
         316*sqr(y))*sqr(y) + 2*x*y*(-25 + 364*pow6(y) - 5*pow8(y) +
         1950*quad(y) + 596*sqr(y))))/(8.*y*pow6(-1 + sqr(y)));

   return (3*((-1 + sqr(x))*(-2*(cube(y) + cube(y)*quad(x) + 2*y*sqr(x)*(-2 +
      sqr(y)) - 3*cube(x)*(-1 + sqr(y)) - pow5(x)*(-1 + sqr(y)))*(-1 +
      sqr(y)) + cube(y)*std::log(sqr(y))*(1 + sqr(y))*sqr(-1 + sqr(x))) -
      std::log(sqr(x))*sqr(x)*(2*x - 3*y + 6*cube(x) + y*quad(x) -
      6*y*sqr(x))*sqr(-1 + sqr(y))))/(4.*cube(-1 + sqr(x))*sqr(x -
      y)*sqr(-1 + sqr(y)));
}

/// First derivative of f5 w.r.t. 2nd argument
double D01f5(double x, double y) noexcept
{
   if (is_equal(x, 1., 0.001) && is_equal(y, 1., 0.001))
      return (sqr(x)*(-91 + 90*y - 27*sqr(y)) + 2*x*(153 - 172*y + 54*sqr(y))
         - 3*(39 - 150*y + 55*sqr(y)))/560.;

   if (is_equal(x, 1., 0.001))
      return (-6*std::log(sqr(y))*sqr(y)*(9 - 3*y - 15*cube(y) + 6*quad(y) + x*(-9 -
         5*y + pow5(y) - 3*quad(y) - 24*sqr(y)) + 23*sqr(y) + sqr(x)*(3 + 3*y
         + 3*cube(y) + 2*quad(y) + 9*sqr(y))) + (-1 + sqr(y))*(3 + y -
         74*cube(y) - 35*pow5(y) + 12*pow6(y) + 75*quad(y) + 4*x*(-1 - y -
         10*cube(y) + 5*pow5(y) - 19*quad(y) - 34*sqr(y)) + 138*sqr(y) +
         sqr(x)*(1 + 3*y + 30*cube(y) + 3*pow5(y) + 37*quad(y) +
         46*sqr(y))))/(8.*cube(1 + y)*pow6(-1 + y));

   if (is_equal(y, 1., 0.001))
      return (30*cube(x)*std::log(sqr(x))*(1 + sqr(x))*(6 + 2*x*(-2 + y) - 8*y +
         sqr(x) + 3*sqr(y)) - (-1 + sqr(x))*(pow5(x)*(43 + 24*y + 3*sqr(y)) -
         4*sqr(x)*(16 - 27*y + 6*sqr(y)) - 2*(6 - 17*y + 6*sqr(y)) +
         2*quad(x)*(-82 - 11*y + 18*sqr(y)) + x*(71 - 152*y + 51*sqr(y)) +
         2*cube(x)*(153 - 176*y + 63*sqr(y))))/(40.*pow6(-1 + x)*sqr(1 + x));

   if (is_equal(x, y, 0.001))
      return ((-1 + sqr(y))*(sqr(x)*(-3 - 92*pow6(y) + pow8(y) - 590*quad(y) -
         276*sqr(y)) - 4*x*y*(7 - 68*pow6(y) + pow8(y) - 340*quad(y) -
         80*sqr(y)) + sqr(y)*(-35 - 276*pow6(y) + 9*pow8(y) - 662*quad(y) +
         4*sqr(y))) + 6*y*std::log(sqr(y))*(y*(2 + 101*pow6(y) + 9*pow8(y) +
         45*quad(y) + 3*sqr(y)) - x*(-1 + 146*pow6(y) + 9*pow8(y) +
         160*quad(y) + 6*sqr(y)) + y*sqr(x)*(15 + 3*pow6(y) + 57*quad(y) +
         85*sqr(y))))/(8.*y*pow6(-1 + sqr(y)));

   return (3*(cube(x)*cube(-1 + sqr(y))*std::log(sqr(x))*(1 + sqr(x)) - (-1 +
      sqr(x))*(std::log(sqr(y))*sqr(y)*(-2*(y + 3*cube(y)) + 2*(y +
      3*cube(y))*sqr(x) + cube(x)*(-3 + quad(y) - 6*sqr(y)) + x*(3 -
      quad(y) + 6*sqr(y))) + 2*(-1 + sqr(y))*(-4*x*sqr(y) + cube(y)*(3 +
      sqr(y)) - cube(y)*sqr(x)*(3 + sqr(y)) + cube(x)*sqr(1 +
      sqr(y))))))/(4.*cube(-1 + sqr(y))*sqr(x - y)*sqr(-1 + sqr(x)));
}

/// First derivative of f6 w.r.t. 1st argument
double D10f6(double x, double y) noexcept
{
   if (is_equal(x, 1., 0.001) && is_equal(y, 1., 0.001))
      return (259 + 99*y - 43*sqr(y) - 3*sqr(x)*(22 - 21*y + 6*sqr(y))
         + 2*x*(54 - 88*y + 27*sqr(y)))/490.;

   if (is_equal(x, 1., 0.001))
      return (30*std::log(sqr(y))*pow5(y)*(6 + 2*x*(-4 + y) - 4*y + 3*sqr(x) + sqr(y))
         + (-1 + sqr(y))*(-14 + 32*y - 223*cube(y) - 19*pow5(y) + 82*quad(y) +
         52*sqr(y) + sqr(x)*(6 - 33*y - 63*cube(y) + 6*pow5(y) - 78*quad(y) +
         72*sqr(y)) - 2*x*(6 - 38*y - 108*cube(y) + 26*pow5(y) - 73*quad(y) +
         97*sqr(y))))/(35.*pow6(-1 + y)*sqr(1 + y));

   if (is_equal(y, 1., 0.001))
      return (-6*std::log(sqr(x))*quad(x)*(y*cube(x) + 5*(3 - 3*y + sqr(y)) + 3*x*(-3
         - y + sqr(y)) + sqr(x)*(4 - 3*y + 2*sqr(y))) + (-1 + sqr(x))*(-1 -
         3*y + 12*pow6(x) + x*(-9 + 17*y - 5*sqr(y)) + sqr(y) + pow5(x)*(-36 +
         17*y + 4*sqr(y)) + sqr(x)*(47 - 30*y + 7*sqr(y)) + cube(x)*(-9 - 46*y
         + 19*sqr(y)) + quad(x)*(56 - 75*y + 34*sqr(y))))/(7.*cube(1 +
         x)*pow6(-1 + x));

   if (is_equal(x, y, 0.001))
      return (6*y*std::log(sqr(y))*(3*sqr(x)*(5 + 2*pow6(y) + 33*quad(y) + 40*sqr(y))
         + sqr(y)*(5 + 12*pow6(y) + 141*quad(y) + 82*sqr(y)) - 2*x*y*(5 +
         8*pow6(y) + 117*quad(y) + 110*sqr(y))) + (-1 +
         sqr(y))*(3*y*sqr(x)*(-107 + pow6(y) - 61*quad(y) - 313*sqr(y)) +
         y*(-6 - 375*pow6(y) + 13*pow8(y) - 975*quad(y) - 97*sqr(y)) + x*(-12
         + 486*pow6(y) - 10*pow8(y) + 2022*quad(y) + 394*sqr(y))))/(7.*pow6(-1
         + sqr(y)));

   return (6*((-1 + sqr(x))*(-((-1 + sqr(y))*(cube(y) + y*sqr(x)*(-3 + sqr(y))
      - 2*cube(x)*(-1 + sqr(y)) - 2*pow5(x)*(-1 + sqr(y)) + y*quad(x)*(-1 +
      2*sqr(y)))) + std::log(sqr(y))*pow5(y)*sqr(-1 + sqr(x))) -
      std::log(sqr(x))*quad(x)*(4*x - 5*y + y*sqr(x))*sqr(-1 +
      sqr(y))))/(7.*cube(-1 + sqr(x))*sqr(x - y)*sqr(-1 + sqr(y)));
}

/// First derivative of f6 w.r.t. 2nd argument
double D01f6(double x, double y) noexcept
{
   if (is_equal(x, 1., 0.001) && is_equal(y, 1., 0.001))
      return (259 + 108*y + sqr(x)*(-43 + 54*y - 18*sqr(y)) - 66*sqr(y)
         + x*(99 - 176*y + 63*sqr(y)))/490.;

   if (is_equal(x, 1., 0.001))
      return (-6*std::log(sqr(y))*quad(y)*(15 - 9*y + x*(-15 - 3*y + cube(y) -
         3*sqr(y)) + 4*sqr(y) + sqr(x)*(5 + 3*y + 2*sqr(y))) + (-1 +
         sqr(y))*(-1 - 9*y - 9*cube(y) - 36*pow5(y) + 12*pow6(y) + 56*quad(y)
         + x*(-3 + 17*y - 46*cube(y) + 17*pow5(y) - 75*quad(y) - 30*sqr(y)) +
         47*sqr(y) + sqr(x)*(1 - 5*y + 19*cube(y) + 4*pow5(y) + 34*quad(y) +
         7*sqr(y))))/(7.*cube(1 + y)*pow6(-1 + y));

   if (is_equal(y, 1., 0.001))
      return (30*std::log(sqr(x))*pow5(x)*(6 + 2*x*(-2 + y) - 8*y + sqr(x) + 3*sqr(y))
         + (-1 + sqr(x))*(quad(x)*(82 + 146*y - 78*sqr(y)) + cube(x)*(-223 +
         216*y - 63*sqr(y)) + x*(32 + 76*y - 33*sqr(y)) + 2*(-7 - 6*y +
         3*sqr(y)) + pow5(x)*(-19 - 52*y + 6*sqr(y)) + 2*sqr(x)*(26 - 97*y +
         36*sqr(y))))/(35.*pow6(-1 + x)*sqr(1 + x));

   if (is_equal(x, y, 0.001))
      return (6*std::log(sqr(y))*(cube(y)*(5 + 6*pow6(y) + 57*quad(y) + 12*sqr(y)) +
         y*sqr(x)*(5 + 2*pow6(y) + 33*quad(y) + 40*sqr(y)) - 2*x*quad(y)*(35 +
         3*quad(y) + 42*sqr(y))) + (-1 + sqr(y))*(y*sqr(x)*(-107 + pow6(y) -
         61*quad(y) - 313*sqr(y)) + y*(-12 - 193*pow6(y) + 9*pow8(y) -
         277*quad(y) - 7*sqr(y)) + x*(-6 + 182*pow6(y) - 4*pow8(y) +
         698*quad(y) + 90*sqr(y))))/(7.*pow6(-1 + sqr(y)));

   return (6*(cube(-1 + sqr(y))*std::log(sqr(x))*pow5(x) - (-1 +
      sqr(x))*(std::log(sqr(y))*quad(y)*(-1 + sqr(x))*(4*y + x*(-5 + sqr(y))) +
      (-1 + sqr(y))*(2*(cube(y) + pow5(y)) - 2*(cube(y) + pow5(y))*sqr(x) -
      x*sqr(y)*(3 + sqr(y)) + cube(x)*(1 + 2*quad(y) +
      sqr(y))))))/(7.*cube(-1 + sqr(y))*sqr(x - y)*sqr(-1 + sqr(x)));
}

/// First derivative of f7 w.r.t. 1st argument
double D10f7(double x, double y) noexcept
{
   if (is_equal(x, 1., 0.001) && is_equal(y, 1., 0.001))
      return (-376 + 207*y - 48*sqr(y) - 9*sqr(x)*(11 - 5*y + sqr(y))
         + 6*x*(57 - 28*y + 6*sqr(y)))/70.;

   if (is_equal(x, 1., 0.001))
      return (30*cube(y)*std::log(sqr(y))*(6 + 2*x*(-4 + y) - 4*y + 3*sqr(x) + sqr(y))
         - (-1 + sqr(y))*(-26 + 103*y + 83*cube(y) + 24*pow5(y) - 82*quad(y) -
         12*sqr(y) + 3*sqr(x)*(-2 + 6*y + 21*cube(y) + 3*pow5(y) - 14*quad(y)
         + 16*sqr(y)) - 2*x*(-11 + 38*y + 68*cube(y) + 14*pow5(y) - 62*quad(y)
         + 43*sqr(y))))/(5.*pow6(-1 + y)*sqr(1 + y));

   if (is_equal(y, 1., 0.001))
      return -(((-1 + sqr(x))*(-4 + y + sqr(x)*(-91 + 106*y - 39*sqr(y)) +
         cube(x)*(65 - 6*y - 11*sqr(y)) + x*(-10 + 21*y - 8*sqr(y)) +
         quad(x)*(-19 + y - 3*sqr(y)) + pow5(x)*(-1 - 3*y + sqr(y))) +
         6*std::log(sqr(x))*sqr(x)*(3*(-2 + y)*cube(x) + 2*quad(x) + 3*(3 - 3*y +
         sqr(y)) + x*(-3 - 5*y + 3*sqr(y)) + sqr(x)*(8 - 9*y +
         4*sqr(y))))/(cube(1 + x)*pow6(-1 + x)));

   if (is_equal(x, y, 0.001))
      return (6*y*std::log(sqr(y))*(y + 4*cube(y) + 123*pow5(y) + 106*pow7(y) +
         6*pow9(y) - 2*x*(-1 + 86*pow6(y) + 4*pow8(y) + 135*quad(y) +
         16*sqr(y)) + 3*y*sqr(x)*(10 + pow6(y) + 24*quad(y) + 45*sqr(y))) -
         (-1 + sqr(y))*(1047*pow6(y) + 173*pow8(y) + 219*quad(y) + sqr(y) -
         2*x*y*(-19 + 121*pow6(y) + 939*quad(y) + 399*sqr(y)) + sqr(x)*(9 +
         93*pow6(y) + 831*quad(y) + 507*sqr(y))))/(y*pow6(-1 + sqr(y)));

   return (6*((-1 + sqr(x))*(-((-1 + sqr(y))*(cube(y) + y*quad(x) -
      4*cube(x)*(-1 + sqr(y)) + y*sqr(x)*(-5 + 3*sqr(y)))) +
      cube(y)*std::log(sqr(y))*sqr(-1 + sqr(x))) - std::log(sqr(x))*sqr(x)*(2*x - 3*y
      + 2*cube(x) - y*sqr(x))*sqr(-1 + sqr(y))))/(cube(-1 + sqr(x))*sqr(x -
      y)*sqr(-1 + sqr(y)));
}

/// First derivative of f7 w.r.t. 2nd argument
double D01f7(double x, double y) noexcept
{
   if (is_equal(x, 1., 0.001) && is_equal(y, 1., 0.001))
      return (-376 + 342*y + sqr(x)*(-48 + 36*y - 9*sqr(y)) - 99*sqr(y)
         + 3*x*(69 - 56*y + 15*sqr(y)))/70.;

   if (is_equal(x, 1., 0.001))
      return -((6*std::log(sqr(y))*sqr(y)*(9 - 3*y - 6*cube(y) + 2*quad(y) + x*(-9 -
         5*y + 3*cube(y) - 9*sqr(y)) + 8*sqr(y) + sqr(x)*(3 + 3*y + 4*sqr(y)))
         + (-1 + sqr(y))*(-4 - 10*y + 65*cube(y) - pow5(y) - 19*quad(y) +
         y*sqr(x)*(-8 - 39*y - 3*cube(y) + quad(y) - 11*sqr(y)) - 91*sqr(y) +
         x*(1 + 21*y - 6*cube(y) - 3*pow5(y) + quad(y) + 106*sqr(y))))/(cube(1
         + y)*pow6(-1 + y)));

   if (is_equal(y, 1., 0.001))
      return (30*cube(x)*std::log(sqr(x))*(6 + 2*x*(-2 + y) - 8*y + sqr(x) + 3*sqr(y))
         - (-1 + sqr(x))*(-26 + 22*y - 6*sqr(y) + pow5(x)*(24 - 28*y +
         9*sqr(y)) + x*(103 - 76*y + 18*sqr(y)) - 2*quad(x)*(41 - 62*y +
         21*sqr(y)) + 2*sqr(x)*(-6 - 43*y + 24*sqr(y)) + cube(x)*(83 - 136*y +
         63*sqr(y))))/(5.*pow6(-1 + x)*sqr(1 + x));

   if (is_equal(x, y, 0.001))
      return (6*y*std::log(sqr(y))*(y*(2 + 44*pow6(y) + 3*pow8(y) + 33*quad(y) -
         2*sqr(y)) - x*(-1 + 62*pow6(y) + 3*pow8(y) + 90*quad(y) + 6*sqr(y)) +
         y*sqr(x)*(10 + pow6(y) + 24*quad(y) + 45*sqr(y))) - (-1 +
         sqr(y))*((23 + 83*pow6(y) + 385*quad(y) - 11*sqr(y))*sqr(y) -
         2*x*y*(-11 + 45*pow6(y) + 331*quad(y) + 115*sqr(y)) + sqr(x)*(3 +
         31*pow6(y) + 277*quad(y) + 169*sqr(y))))/(y*pow6(-1 + sqr(y)));

   return (6*(cube(x)*cube(-1 + sqr(y))*std::log(sqr(x)) + (-1 +
      sqr(x))*(std::log(sqr(y))*sqr(y)*(2*(y + cube(y)) - 2*(y + cube(y))*sqr(x)
      - x*(3 + sqr(y)) + cube(x)*(3 + sqr(y))) - (-1 + sqr(y))*(4*cube(y) -
      4*cube(y)*sqr(x) + x*(-5 + sqr(y))*sqr(y) + cube(x)*(1 +
      3*sqr(y))))))/(cube(-1 + sqr(y))*sqr(x - y)*sqr(-1 + sqr(x)));
}

/// First derivative of f8 w.r.t. 1st argument
double D10f8(double x, double y) noexcept
{
   if (is_equal(x, 1., 0.001) && is_equal(y, 1., 0.001))
      return (288 - 232*y + x*(-356 + 234*y - 60*sqr(y)) + 63*sqr(y)
         + 9*sqr(x)*(13 - 8*y + 2*sqr(y)))/280.;

   if (is_equal(x, 1., 0.001))
      return ((-24 + 122*y + 12*cube(y) - 14*pow5(y) + 37*quad(y) - 2*x*(-14 +
         67*y - 43*cube(y) + 6*pow5(y) + 2*quad(y) - 108*sqr(y)) +
         3*sqr(x)*(-3 + 14*y - 16*cube(y) + 2*pow5(y) - 6*quad(y) - 21*sqr(y))
         - 223*sqr(y))*(-1 + sqr(y)) + 30*std::log(sqr(y))*quad(y)*(6 + 2*x*(-4 +
         y) - 4*y + 3*sqr(x) + sqr(y)))/(20.*pow6(-1 + y)*sqr(1 + y));

   if (is_equal(y, 1., 0.001))
      return (-6*cube(x)*std::log(sqr(x))*((-3 + 2*y)*cube(x) + quad(x) + 4*(3 - 3*y +
         sqr(y)) + 3*sqr(x)*(2 - 2*y + sqr(y)) + x*(-6 - 4*y + 3*sqr(y))) +
         (-1 + sqr(x))*(-6 + 4*y + (17 + 4*y)*pow5(x) - sqr(y) + x*(26 - 14*y
         + 3*sqr(y)) + quad(x)*(-51 + 10*y + 8*sqr(y)) + sqr(x)*(3 - 26*y +
         11*sqr(y)) + cube(x)*(71 - 98*y + 39*sqr(y))))/(4.*cube(1 +
         x)*pow6(-1 + x));

   if (is_equal(x, y, 0.001))
      return (6*std::log(sqr(y))*(153*pow6(y) + 52*pow8(y) + 34*quad(y) + sqr(y) +
         3*sqr(x)*(1 + 10*pow6(y) + 45*quad(y) + 24*sqr(y)) - 2*x*y*(-1 +
         38*pow6(y) + 147*quad(y) + 56*sqr(y))) - (-1 + sqr(y))*(6 +
         771*pow6(y) + 23*pow8(y) + 651*quad(y) - 11*sqr(y) - 2*x*y*(17 +
         13*pow6(y) + 615*quad(y) + 795*sqr(y)) + sqr(x)*(93 + 9*pow6(y) +
         507*quad(y) + 831*sqr(y))))/(4.*pow6(-1 + sqr(y)));

   return (-3*(-((-1 + sqr(x))*(-((-1 + sqr(y))*(sqr(x)*(1 - 3*sqr(y)) +
      quad(x)*(3 - 2*sqr(y)) + 2*x*y*(-1 + sqr(y)) + 2*y*cube(x)*(-1 +
      sqr(y)) + sqr(y))) + std::log(sqr(y))*quad(y)*sqr(-1 + sqr(x)))) +
      cube(x)*(3*x - 4*y + cube(x))*std::log(sqr(x))*sqr(-1 +
      sqr(y))))/(2.*cube(-1 + sqr(x))*sqr(x - y)*sqr(-1 + sqr(y)));
}

/// First derivative of f8 w.r.t. 2nd argument
double D01f8(double x, double y) noexcept
{
   if (is_equal(x, 1., 0.001) && is_equal(y, 1., 0.001))
      return (288 - 356*y + x*(-232 + 234*y - 72*sqr(y)) + 117*sqr(y)
         + 3*sqr(x)*(21 - 20*y + 6*sqr(y)))/280.;

   if (is_equal(x, 1., 0.001))
      return (-6*cube(y)*std::log(sqr(y))*(12 - 6*y - 3*cube(y) + quad(y) + 2*x*(-6 -
         2*y + cube(y) - 3*sqr(y)) + 6*sqr(y) + sqr(x)*(4 + 3*y + 3*sqr(y))) +
         (-1 + sqr(y))*(-6 + 26*y + 71*cube(y) + 17*pow5(y) - 51*quad(y) +
         2*x*(2 - 7*y - 49*cube(y) + 2*pow5(y) + 5*quad(y) - 13*sqr(y)) +
         3*sqr(y) + sqr(x)*(-1 + 3*y + 39*cube(y) + 8*quad(y) +
         11*sqr(y))))/(4.*cube(1 + y)*pow6(-1 + y));

   if (is_equal(y, 1., 0.001))
      return (30*std::log(sqr(x))*quad(x)*(6 + 2*x*(-2 + y) - 8*y + sqr(x) + 3*sqr(y))
         + (-1 + sqr(x))*(-24 + 28*y + sqr(x)*(-223 + 216*y - 63*sqr(y)) +
         cube(x)*(12 + 86*y - 48*sqr(y)) + quad(x)*(37 - 4*y - 18*sqr(y)) -
         9*sqr(y) + 2*pow5(x)*(-7 - 6*y + 3*sqr(y)) + 2*x*(61 - 67*y +
         21*sqr(y))))/(20.*pow6(-1 + x)*sqr(1 + x));

   if (is_equal(x, y, 0.001))
      return (6*std::log(sqr(y))*(sqr(y)*(3 + 24*pow6(y) + 51*quad(y) + 2*sqr(y)) -
         2*x*y*(-1 + 14*pow6(y) + 51*quad(y) + 16*sqr(y)) + sqr(x)*(1 +
         10*pow6(y) + 45*quad(y) + 24*sqr(y))) - (-1 + sqr(y))*(6 +
         325*pow6(y) + 13*pow8(y) + 133*quad(y) + 3*sqr(y) - 2*x*y*(-7 +
         5*pow6(y) + 223*quad(y) + 259*sqr(y)) + sqr(x)*(31 + 3*pow6(y) +
         169*quad(y) + 277*sqr(y))))/(4.*pow6(-1 + sqr(y)));

   return (-3*(-(cube(-1 + sqr(y))*std::log(sqr(x))*quad(x)) + (-1 + sqr(x))*((-1 +
      sqr(y))*(-2*x*(y + cube(y)) + 2*cube(x)*(y + cube(y)) + 3*quad(y) +
      sqr(x)*(1 - 2*quad(y) - 3*sqr(y)) + sqr(y)) +
      cube(y)*std::log(sqr(y))*(4*x - 4*cube(x) - y*(3 + sqr(y)) + y*sqr(x)*(3 +
      sqr(y))))))/(2.*cube(-1 + sqr(y))*sqr(x - y)*sqr(-1 + sqr(x)));
}

/// Second derivative of F1
double D2F1(double x) noexcept
{
   if (is_equal(x, 1., 0.001)) {
      return -1.3333333333333333 + x + 2*cube(-1 + x) - (31*quad(-1 + x))/14.
         - (8*sqr(-1 + x))/5.;
   }

   return (2 - 6*quad(x) + 4*sqr(x) + 2*std::log(sqr(x))*sqr(x)*(3 +
      sqr(x)))/(x*cube(-1 + sqr(x)));
}

/// Second derivative of F2
double D2F2(double x) noexcept
{
   if (is_equal(x, 1., 0.01)) {
      return (-381 + 832*x + 320*cube(x) - 55*quad(x) - 744*sqr(x))/35.;
   }

   return (12*(5 - 11*pow6(x) - 21*quad(x) + 27*sqr(x) + std::log(sqr(x))*(1 +
      3*pow6(x) + 25*quad(x) + 19*sqr(x))))/pow5(-1 + sqr(x));
}

/// Second derivative of F3
double D2F3(double x) noexcept
{
   if (is_equal(x, 1., 0.001)) {
      return (2*(-392 + 69*x - 465*cube(x) + 120*quad(x) + 528*sqr(x)))/315.;
   }

   return (4*(1 - 17*pow6(x) - 25*quad(x) + 41*sqr(x) +
      2*std::log(sqr(x))*sqr(x)*(9 + 2*quad(x) + 19*sqr(x))))/(3.*x*quad(-1 +
      sqr(x)));
}

/// Second derivative of F4
double D2F4(double x) noexcept
{
   if (is_equal(x, 1., 0.001)) {
      return (-2*(174 - 529*x - 335*cube(x) + 70*quad(x) + 620*sqr(x)))/35.;
   }

   return (4*(-1 + pow6(x) + 9*quad(x) - 9*sqr(x) - 6*std::log(sqr(x))*(quad(x) +
      sqr(x))))/(x*quad(-1 + sqr(x)));
}

/// Second derivative of F5
double D2F5(double x) noexcept
{
   if (is_equal(x, 1., 0.01)) {
      return (-334 + 793*x + 370*cube(x) - 70*quad(x) - 780*sqr(x))/35.;
   }

   return (6*x*(-19 + pow6(x) + 27*quad(x) - 9*sqr(x) - 6*std::log(sqr(x))*(1 +
      2*quad(x) + 5*sqr(x))))/pow5(-1 + sqr(x));
}

/// Second derivative of F6
double D2F6(double x) noexcept
{
   if (is_equal(x, 1., 0.001)) {
      return (109 - 720*x - 560*cube(x) + 120*quad(x) + 981*sqr(x))/210.;
   }

   return (-7 - pow6(x) + 7*quad(x) + sqr(x) - 2*std::log(sqr(x))*(1 +
      5*sqr(x)))/quad(-1 + sqr(x));
}

/// Second derivative of F7
double D2F7(double x) noexcept
{
   if (is_equal(x, 1., 0.001)) {
      return 10 - (208*x)/7. - 16*cube(x) + (22*quad(x))/7. + (1119*sqr(x))/35.;
   }

   return (-6*(2 + pow8(x) - 11*pow6(x) - 27*quad(x) + 35*sqr(x) +
      6*std::log(sqr(x))*sqr(x)*(3 + 5*sqr(x))))/pow5(-1 + sqr(x));
}

namespace {

/// I2abc(a,a,a), squared arguments, a != 0
double I2aaa(double a, double b, double c) noexcept {
   const double ba = b - a;
   const double ca = c - a;
   const double a2 = sqr(a);
   const double a3 = a2*a;

   return 0.5/a + (-ba - ca)/(6.0*a2) + (sqr(ba) + ba*ca + sqr(ca))/(12.0*a3);
}

/// I2abc(a,a,c), squared arguments, a != c
double I2aac(double a, double b, double c) noexcept {
   const double ba = b - a;
   const double ac = a - c;
   const double a2 = sqr(a);
   const double a3 = a2*a;
   const double c2 = sqr(c);
   const double c3 = c2*c;
   const double ac2 = sqr(ac);
   const double ac3 = ac2*ac;
   const double ac4 = ac2*ac2;
   const double lac = std::log(a/c);

   return (ac - c*lac)/ac2
      + ba*(-a2 + c2 + 2*a*c*lac)/(2.0*a*ac3)
      + sqr(ba)*((2*a3 + 3*a2*c - 6*a*c2 + c3 - 6*a2*c*lac)/(6.*a2*ac4));
}

/// I2abc(a,a,0), squared arguments, a != 0
double I2aa0(double a, double b) noexcept {
   const double a2 = sqr(a);
   const double a3 = a2*a;
   const double ba = b - a;
   const double ba2 = sqr(ba);

   return 1.0/a - ba/(2.0*a2) + ba2/(3.0*a3);
}

/// I2abc(0,b,c), squared arguments, b != c
double I20bc(double b, double c) noexcept {
   return std::log(b/c)/(b - c);
}

} // anonymous namespace

double Iabc(double a, double b, double c) noexcept {
   const double eps = 10.0*std::numeric_limits<double>::epsilon();

   if ((is_zero(a, eps) && is_zero(b, eps) && is_zero(c, eps)) ||
       (is_zero(a, eps) && is_zero(b, eps)) ||
       (is_zero(a, eps) && is_zero(c, eps)) ||
       (is_zero(b, eps) && is_zero(c, eps))) {
      return 0.0;
   }

   const double a2 = sqr(a);
   const double b2 = sqr(b);
   const double c2 = sqr(c);
   const double eps_eq = 0.001;

   if (is_close(a2, b2, eps_eq) && is_close(a2, c2, eps_eq)) {
      return I2aaa(a2, b2, c2);
   }

   if (is_close(a2, b2, eps_eq)) {
      if (is_zero(c, eps)) {
         return I2aa0(a2, b2);
      }
      return I2aac(a2, b2, c2);
   }

   if (is_close(b2, c2, eps_eq)) {
      if (is_zero(a, eps)) {
         return I2aa0(b2, c2);
      }
      return I2aac(b2, c2, a2);
   }

   if (is_close(a2, c2, eps_eq)) {
      if (is_zero(b, eps)) {
         return I2aa0(a2, c2);
      }
      return I2aac(a2, c2, b2);
   }

   if (is_zero(a, eps)) {
      return I20bc(b2, c2);
   }

   if (is_zero(b, eps)) {
      return I20bc(c2, a2);
   }

   if (is_zero(c, eps)) {
      return I20bc(a2, b2);
   }

   return (+ a2 * b2 * std::log(a2/b2)
           + b2 * c2 * std::log(b2/c2)
           + c2 * a2 * std::log(c2/a2))
           / ((a2 - b2) * (b2 - c2) * (a2 - c2));
}

/// Delta function from hep-ph/0907.47682v1
double delta_xyz(double x, double y, double z) noexcept
{
   return sqr(x)+sqr(y)+sqr(z)-2*(x*y+x*z+y*z);
}

namespace {
   /// lambda^2(u,v)
   double lambda_2(double u, double v) noexcept
   {
      return sqr(1 - u - v) - 4*u*v;
   }

   /// u < 1 && v < 1, lambda^2(u,v) > 0; note: phi_pos(u,v) = phi_pos(v,u)
   double phi_pos(double u, double v) noexcept
   {
      const double eps = 1.0e-7;

      if (is_equal(u, 1.0, eps) && is_equal(v, 1.0, eps)) {
         return 2.343907238689459;
      }

      const double Pi = 3.141592653589793;
      const auto lambda = std::sqrt(lambda_2(u,v));

      if (is_equal(u, v, eps)) {
         return (-(sqr(std::log(u)))
                 + 2*sqr(std::log((1 - lambda)/2.))
                 - 4*dilog((1 - lambda)/2.)
                 + sqr(Pi)/3.)/lambda;
      }

      return (-(std::log(u)*std::log(v))
              + 2*std::log((1 - lambda + u - v)/2.)*std::log((1 - lambda - u + v)/2.)
              - 2*dilog((1 - lambda + u - v)/2.)
              - 2*dilog((1 - lambda - u + v)/2.)
              + sqr(Pi)/3.)/lambda;
   }

   /// lambda^2(u,v) < 0, u = 1
   double phi_neg_1v(double v, double lambda) noexcept
   {
      return 2*(+ clausen_2(2*std::acos((2 - v)/2))
                + 2*clausen_2(2*std::acos(0.5*std::sqrt(v))))/lambda;
   }

   /// lambda^2(u,v) < 0; note: phi_neg(u,v) = phi_neg(v,u)
   double phi_neg(double u, double v) noexcept
   {
      const double eps = 1.0e-7;

      if (is_equal(u, 1.0, eps) && is_equal(v, 1.0, eps)) {
         return 2.343907238689459;
      }

      const auto lambda = std::sqrt(-lambda_2(u,v));

      if (is_equal(u, 1.0, eps)) {
         return phi_neg_1v(v, lambda);
      }

      if (is_equal(v, 1.0, eps)) {
         return phi_neg_1v(u, lambda);
      }

      if (is_equal(u, v, eps)) {
         return 2*(2*clausen_2(2*std::acos(1/(2.*std::sqrt(u))))
                   + clausen_2(2*std::acos((-1 + 2*u)/(2.*std::abs(u)))))/lambda;
      }

      return 2*(+ clausen_2(2*std::acos((1 + u - v)/(2.*std::sqrt(u))))
                + clausen_2(2*std::acos((1 - u + v)/(2.*std::sqrt(v))))
                + clausen_2(2*std::acos((-1 + u + v)/(2.*std::sqrt(u*v)))))/lambda;
   }

   /**
    * Phi(u,v) with u = x/z, v = y/z.
    *
    * The following identities hold:
    * Phi(u,v) = Phi(v,u) = Phi(1/u,v/u)/u = Phi(1/v,u/v)/v
    */
   double phi_uv(double u, double v) noexcept
   {
      const auto lambda = lambda_2(u,v);

      if (is_zero(lambda, std::numeric_limits<double>::epsilon())) {
         // phi_uv is always multiplied by lambda.  So, in order to
         // avoid nans if lambda == 0, we simply return 0
         return 0.0;
      }

      if (lambda > 0.) {
         if (u <= 1 && v <= 1) {
            return phi_pos(u,v);
         }
         if (u >= 1 && v/u <= 1) {
            return phi_pos(1./u,v/u)/u;
         }
         // v >= 1 && u/v <= 1
         return phi_pos(1./v,u/v)/v;
      }

      return phi_neg(u,v);
   }
} // anonymous namespace

/**
 * \f$\Phi(x,y,z)\f$ function.  The arguments x, y and z are
 * interpreted as squared masses.
 *
 * Davydychev and Tausk, Nucl. Phys. B397 (1993) 23
 *
 * @param x squared mass
 * @param y squared mass
 * @param z squared mass
 *
 * @return \f$\Phi(x,y,z)\f$
 */
double phi_xyz(double x, double y, double z) noexcept
{
   const auto u = x/z, v = y/z;
   return phi_uv(u,v);
}

/**
 * B0 function for zero momentum, arxiv:0901.2065 Eq (130)
 *
 * @param m1 \f$m_1\f$ (not squared)
 * @param m2 \f$m_2\f$ (not squared)
 * @param scale \f$Q\f$ (not squared)
 *
 * @return \f$B_0(p=0,m_1,m_2,Q)\f$
 */
double B0(double m1, double m2, double scale) noexcept
{
   return Loop_library::get().B0(0, m1*m1, m2*m2, scale*scale).real();
}

/**
 * B0' function for zero momentum, arxiv:0901.2065 Eq (130)
 *
 * @param m1 \f$m_1\f$ (not squared)
 * @param m2 \f$m_2\f$ (not squared)
 *
 * @return \f$B_0'(p=0,m_1,m_2)\f$
 */
double DB0(double m1, double m2) noexcept
{
   const double eps = 10.*std::numeric_limits<double>::epsilon();
   const double m12 = sqr(m1);
   const double m14 = sqr(m12);
   const double m22 = sqr(m2);
   const double m24 = sqr(m22);

   if (is_zero(m12, eps) || is_zero(m22, eps)) {
      return 0.;
   }

   if (is_equal_rel(m12, m22, 1e-3)) {
      return 1./(6. * m22);
   }

   return (m14 - m24 + 2*m12*m22*std::log(m22/m12))/
      (2*cube(m12 - m22));
}

/**
 * C0 function for zero momentum, arxiv:0901.2065 Eq (130)
 *
 * \f$C_0(0,m_1,m_2,m_3) = -I_{abc}(m_1,m_2,m_3)\f$
 *
 * @param m1 \f$m_1\f$ (not squared)
 * @param m2 \f$m_2\f$ (not squared)
 * @param m3 \f$m_3\f$ (not squared)
 *
 * @return \f$C_0(p=0,m_1,m_2,m_3)\f$
 */
double C0(double m1, double m2, double m3) noexcept
{
   return softsusy::c0(m1, m2, m3);
}

double D0(double m1, double m2, double m3, double m4) noexcept
{
   return softsusy::d0(m1,m2,m3,m4);
}

/**
 * \f$\tilde{D}_2\f$ for zero momentum, arxiv:0901.2065 Eq (131)
 *
 * @param m1  \f$m_1\f$ (not squared)
 * @param m2  \f$m_2\f$ (not squared)
 * @param m3  \f$m_3\f$ (not squared)
 * @param m4  \f$m_4\f$ (not squared)
 *
 * @return \f$\tilde{D}_2(m_1,m_2,m_3,m_4)\f$
 */
double D2t(double m1, double m2, double m3, double m4) noexcept
{
   return C0(m2, m3, m4) + m1*m1 * D0(m1, m2, m3, m4);
}

/**
 * \f$\tilde{D}_4\f$ for zero momentum, arxiv:0901.2065 Eq (131)
 *
 * @param m1  \f$m_1\f$ (not squared)
 * @param m2  \f$m_2\f$ (not squared)
 * @param m3  \f$m_3\f$ (not squared)
 * @param m4  \f$m_4\f$ (not squared)
 * @param scale \f$Q\f$ (not squared)
 *
 * @return \f$\tilde{D}_4(m_1,m_2,m_3,m_4,Q)\f$
 */
double D4t(double m1, double m2, double m3, double m4, double scale) noexcept
{
   return B0(m3, m4, scale) + (m1*m1 + m2*m2) * C0(m2, m3, m4)
      + quad(m1) * D0(m1, m2, m3, m4);
}

/**
 * \f$W\f$ for zero momentum, arxiv:0901.2065 Eq (130)
 *
 * @param m1 \f$m_1\f$ (not squared)
 * @param m2 \f$m_2\f$ (not squared)
 * @param scale \f$Q\f$ (not squared)
 *
 * @return \f$W(m_1,m_2,Q)\f$
 */
double W(double m1, double m2, double scale) noexcept
{
   const double eps = 10.*std::numeric_limits<double>::epsilon();
   const double m12 = sqr(m1);
   const double m14 = sqr(m12);
   const double m22 = sqr(m2);
   const double m24 = sqr(m22);
   const double m26 = m24 * m22;
   const double Q2  = sqr(scale);

   if (is_zero(m12, eps) || is_zero(m22, eps) || is_zero(Q2, eps)) {
      return 0.;
   }

   if (is_equal_rel(m12, m22, 1e-3)) {
      return 2./3. - 2. * std::log(Q2/m22);
   }

   return (- 2*std::log(Q2/m12)
           - std::log(m22/m12)*(2*m26 - 6*m12*m24)/cube(m12 - m22)
           - (m14 - 6*m22*m12 + m24)/sqr(m12 - m22));
}

} // namespace threshold_loop_functions
} // namespace flexiblesusy
