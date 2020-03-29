// ====================================================================
// This file is part of GM2Calc.
//
// GM2Calc is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// GM2Calc is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GM2Calc.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#include "gm2_ffunctions.hpp"
#include "gm2_dilog.hpp"
#include "gm2_log.hpp"

#include <cmath>
#include <complex>
#include <iostream>
#include <limits>

namespace gm2calc {

namespace {
   const double eps = 10.0*std::numeric_limits<double>::epsilon();

   /// returns number squared
   template <typename T> T sqr(T x) { return x*x; }

   /// returns number cubed
   template <typename T> T pow3(T x) { return x*x*x; }

   /// returns number to the power 4
   template <typename T> T pow4(T x) { return x*x*x*x; }

   bool is_zero(double a, double prec)
   {
      return std::fabs(a) < prec;
   }

   bool is_equal(double a, double b, double prec)
   {
      const double max = std::max(std::abs(a), std::abs(b));
      return is_zero(a - b, prec*(1.0 + max));
   }

} // anonymous namespace

double F1C(double x) {
   if (is_zero(x, eps)) {
      return 4.0;
   }

   const double d = x - 1.0;

   if (is_equal(x, 1.0, 0.03)) {
      return 1.0 + d*(-0.6 + d*(0.4 + d*(-2.0/7.0
         + d*(3.0/14.0 + d*(-1.0/6.0
         + 2.0/15.0*d)))));
   }

   return 2.0/pow4(d)*(2.0 + x*(3.0 + 6.0*std::log(x) + x*(-6.0 + x)));
}

double F2C(double x) {
   if (is_zero(x, eps)) {
      return 0.0;
   }

   if (is_equal(x, 1.0, 0.03)) {
      const double d = x - 1.0;

      return 1.0 + d*(-0.75 + d*(0.6 + d*(-0.5 + d*(3.0/7.0
         + d*(-0.375 + 1.0/3.0*d)))));
   }

   return 3.0/(2.0*pow3(1.0 - x))*(-3.0 - 2.0*std::log(x) + x*(4.0 - x));
}

double F3C(double x) {
   const double d = x - 1.0;

   if (is_equal(x, 1.0, 0.03)) {
      return 1.0
         + d*(1059.0/1175.0
         + d*(-4313.0/3525.0
         + d*(70701.0/57575.0
         + d*(-265541.0/230300.0
         + d*(+48919.0/46060.0
         - 80755.0/82908.0*d)))));
   }

   const double lx = std::log(x);
   const double x2 = sqr(x);

   return 4.0/(141.0*pow4(d)) * (
      + (1.0 - x) * (151.0 * x2 - 335.0 * x + 592.0)
      + 6.0 * (21.0 * pow3(x) - 108.0 * x2 - 93.0 * x + 50.0) * lx
      - 54.0 * x * (x2 - 2.0 * x - 2.0) * sqr(lx)
      - 108.0 * x * (x2 - 2.0 * x + 12.0) * dilog(1.0 - x)
      );
}

double F4C(double x) {
   if (is_zero(x, eps)) {
      return 0.0;
   }

   if (is_equal(x, 1.0, 0.03)) {
      const double d = x - 1.0;

      return 1.0
         + d*(-45.0/122.0
         + d*(941.0/6100.0
         + d*(-17.0/305.0
         + d*(+282.0/74725.0
         + d*(+177.0/6832.0 - 47021.0/1076040.0*d)))));
   }

   const double lx = std::log(x);
   const double x2 = sqr(x);

   return -9.0/(122.0 * pow3(1.0 - x)) * (
      + 8.0 * (x2 - 3.0 * x + 2.0)
      + (11.0 * x2 - 40.0 * x + 5.0) * lx
      - 2.0 * (x2 - 2.0 * x - 2.0) * sqr(lx)
      - 4.0 * (x2 - 2.0 * x + 9.0) * dilog(1.0 - x)
      );
}

double F1N(double x) {
   if (is_zero(x, eps)) {
      return 2.0;
   }

   const double d = x - 1.0;

   if (is_equal(x, 1.0, 0.03)) {
      return 1.0 + d*(-0.4 + d*(0.2 + d*(-4.0/35.0
         + d*(1.0/14.0 + d*(-1.0/21.0 + 1.0/30.0*d)))));
   }

   return 2.0/pow4(d)*(1.0 + x*(-6.0 + x*(+3.0 - 6.0 * std::log(x) + 2.0 * x)));
}

double F2N(double x) {
   if (is_zero(x, eps)) {
      return 3.0;
   }

   if (is_equal(x, 1.0, 0.04)) {
      const double d = x - 1.0;

      return 1. + d*(-0.5 + d*(0.3 + d*(-0.2
         + d*(1.0/7.0 + d*(-3.0/28.0 + 1.0/12.0*d)))));
   }

   return 3.0/pow3(1.0 - x) * (1.0 + x*(2.0 * std::log(x) - x));
}

double F3N(double x) {
   if (is_zero(x, eps)) {
      return 8.0/105.0;
   }

   const double d = x - 1.0;

   if (is_equal(x, 1.0, 0.03)) {
      return 1.0 + d*(76/875.0 + d*(-431/2625.0 + d*(5858/42875.0
         + d*(-3561/34300.0 + d*(23/294.0 - 4381/73500.0*d)))));
   }

   const double x2 = sqr(x);

   return 4.0/(105.0 * pow4(d)) * (
      + (1.0 - x) * (- 97.0 * x2 - 529.0 * x + 2.0)
      + 6.0 * x2 * (13.0 * x + 81.0) * std::log(x)
      + 108.0 * x * (7.0 * x + 4.0) * dilog(1.0 - x)
      );
}

double F4N(double x) {
   const double PI = 3.14159265358979323846;
   const double PI2 = PI * PI;

   if (is_zero(x, eps)) {
      return -3.0/4.0*(-9.0 + PI2);
   }

   if (is_equal(x, 1.0, 0.03)) {
      const double d = x - 1.0;

      return 1.0 + sqr(d)*(-111.0/800.0 + d*(59.0/400.0 + d*(-129.0/980.0
         + d*(177.0/1568.0 - 775.0/8064.0*d))));
   }

   return -2.25/pow3(1.0 - x) * (
      + (x + 3.0) * (x * std::log(x) + x - 1.0)
      + (6.0 * x + 2.0) * dilog(1.0 - x)
      );
}

namespace {

/// Fb(1,1)
double Fb11(double x, double y) {
   const double x1 = x - 1.0;
   const double y1 = y - 1.0;

   return
      + 1.0/12.0 + (-0.05 + y1/30.)*y1
      + x1*(-0.05 + (1.0/30.0 - y1/42.0)*y1
      + x1*(1.0/30.0 + (-1.0/42.0 + y1/56.0)*y1));
}

/// Fb(x,1), x != 1, x != 0
double Fb1(double x, double y) {
   const double x1 = x - 1.0;
   const double y1 = y - 1.0;
   const double lx = std::log(x);
   const double x14 = pow4(x1);
   const double x15 = x14*x1;
   const double x16 = x15*x1;

   return (2.0 + x*(3.0 + 6.0*lx + x*(-6.0 + x)))/(6.0*x14)
      + y1*(3.0 + x*(10.0 + 12.0*lx + x*(-18.0 + x*(6.0 - x))))/(12.0*x15)
      + sqr(y1)*(12.0 + x*(65.0 + 60.0*lx + x*(-120.0 + x*(60.0 + x*(-20.0 + 3.0*x)))))/(60.*x16);
}

/// Fb(x,x), x != 1, x != 0
double Fbx(double x, double y) {
   const double x1 = x - 1.0;
   const double d = y - x;
   const double lx = std::log(x);
   const double x14 = pow4(x1);
   const double x15 = x14*x1;
   const double x16 = x15*x1;

   return (-5.0 - 2.0*lx + x*(4.0 - 4.0*lx + x))/(2.0*x14)
      - d*(-1.0 + x*(-9.0 - 6.0*lx + x*(9.0 - 6.0*lx + x)))/(2.0*x15*x)
      - sqr(d)*(-1.0 + x*(12.0 + x*(36.0 + 36.0*lx + x*(-44.0 + 24.0*lx - 3.0*x))))/(6.*x16*sqr(x));
}

} // anonymous namespace

double Fb(double x, double y) {
   if ((is_zero(x, eps) && is_zero(y, eps)) || is_zero(x, eps) ||
       is_zero(y, eps)) {
      return 0.0;
   }

   if (is_equal(x, 1.0, 0.01) && is_equal(y, 1.0, 0.01)) {
      return Fb11(x, y);
   }

   if (is_equal(x, 1.0, 0.01)) {
      return Fb1(y, x);
   }

   if (is_equal(y, 1.0, 0.01)) {
      return Fb1(x, y);
   }

   if (is_equal(x, y, 0.01)) {
      return Fbx(x, y);
   }

   return - (G4(x) - G4(y)) / (x - y);
}

namespace {

/// Fa(1,1)
double Fa11(double x, double y) {
   const double x1 = x - 1.0;
   const double y1 = y - 1.0;

   return
      0.25 + (-0.2 + y1/6.)*y1
      + x1*(-0.2 + (1.0/6.0 - y1/7.0)*y1)
      + sqr(x1)*(1.0/6.0 + (-1.0/7.0 + y1/8.0)*y1);
}

/// Fa(x,1), x != 1, x != 0
double Fa1(double x, double y) {
   const double x1 = x - 1.0;
   const double y1 = y - 1.0;
   const double lx = std::log(x);
   const double x14 = pow4(x1);
   const double x15 = x14*x1;
   const double x16 = x15*x1;

   return (-11.0 - 6.0*lx + x*(18.0 + x*(-9.0 + 2.0*x)))/(6.0*x14)
      + y1*(-25.0 - 12.0*lx + x*(48.0 + x*(-36.0 + x*(16.0 - 3.0*x))))/(12.0*x15)
      + sqr(y1)*(-137.0 - 60.0*lx + x*(300.0 + x*(-300.0 + x*(200.0 + x*(-75.0 + 12.0*x)))))/(60.0*x16);
}

/// Fa(x,x), x != 1, x != 0
double Fax(double x, double y) {
   const double x1 = x - 1.0;
   const double d = y - x;
   const double lx = std::log(x);
   const double x14 = pow4(x1);
   const double x15 = x14*x1;
   const double x16 = x15*x1;
   const double x2 = sqr(x);
   const double x3 = x2*x;

   return (2.0 + x*(3.0 + 6.0*lx + x*(-6.0 + x)))/(2.0*x14*x)
      - d*(-1.0 + x*(8.0 + x*(12.0*lx + x*(-8.0 + x))))/(2.0*x15*x2)
      - sqr(d)*(-2.0 + x*(15.0 + x*(-60.0 + x*(20.0 - 60.0*lx + x*(30.0 - 3.0*x)))))/(6.0*x16*x3);
}

} // anonymous namespace

double Fa(double x, double y) {
   if ((is_zero(x, eps) && is_zero(y, eps)) || is_zero(x, eps) ||
       is_zero(y, eps)) {
      return 0.0;
   }

   if (is_equal(x, 1.0, 0.001) && is_equal(y, 1.0, 0.001)) {
      return Fa11(x,y);
   }

   if (is_equal(x, 1.0, 0.001)) {
      return Fa1(y,x);
   }

   if (is_equal(y, 1.0, 0.001)) {
      return Fa1(x,y);
   }

   if (is_equal(x, y, 0.001)) {
      return Fax(x,y);
   }

   return - (G3(x) - G3(y)) / (x - y);
}

double G3(double x) {
   if (is_equal(x, 1.0, 0.01)) {
      const double d = x - 1.0;
      return 1.0/3.0 + d*(-0.25 + d*(0.2 + (-1.0/6.0 + d/7.0)*d));
   }

   return 1.0/(2.0*pow3(x - 1.0))*((x - 1.0)*(x - 3.0) + 2.0*std::log(x));
}

double G4(double x) {
   if (is_equal(x, 1.0, 0.01)) {
      const double d = x - 1.0;
      return 1.0/6.0 + d*(-1.0/12.0 + d*(0.05 + (-1.0/30.0 + d/42.0)*d));
   }

   return 1.0/(2.0*pow3(x - 1.0))*((x - 1.0)*(x + 1.0) - 2.0*x*std::log(x));
}

namespace {

/// I2abc(a,a,a), squared arguments, a != 0
double I2aaa(double a, double b, double c) {
   const double ba = b - a;
   const double ca = c - a;
   const double a2 = sqr(a);
   const double a3 = a2*a;

   return 0.5/a + (-ba - ca)/(6.0*a2) + (sqr(ba) + ba*ca + sqr(ca))/(12.0*a3);
}

/// I2abc(a,a,c), squared arguments, a != c
double I2aac(double a, double b, double c) {
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
double I2aa0(double a, double b) {
   const double a2 = sqr(a);
   const double a3 = a2*a;
   const double ba = b - a;
   const double ba2 = sqr(ba);

   return 1.0/a - ba/(2.0*a2) + ba2/(3.0*a3);
}

/// I2abc(0,b,c), squared arguments, b != c
double I20bc(double b, double c) {
   return std::log(b/c)/(b - c);
}

} // anonymous namespace

double Iabc(double a, double b, double c) {
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

   if (is_equal(a2, b2, eps_eq) && is_equal(a2, c2, eps_eq)) {
      return I2aaa(a2, b2, c2);
   }

   if (is_equal(a2, b2, eps_eq)) {
      if (is_zero(c, eps)) {
         return I2aa0(a2, b2);
      }
      return I2aac(a2, b2, c2);
   }

   if (is_equal(b2, c2, eps_eq)) {
      if (is_zero(a, eps)) {
         return I2aa0(b2, c2);
      }
      return I2aac(b2, c2, a2);
   }

   if (is_equal(a2, c2, eps_eq)) {
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

/**
 * Calculates \f$f_{PS}(z)\f$, Eq (70) arXiv:hep-ph/0609168
 */
double f_PS(double z) {
   double result = 0.0;

   if (z < 0.25) {
      const double y = std::sqrt(1. - 4. * z);
      result = 2.0*z/y*(dilog(1.0 - 0.5*(1.0 - y)/z) - dilog(1.0 - 0.5*(1.0 + y)/z));
   } else {
      const std::complex<double> y = std::sqrt(std::complex<double>(1.0 - 4.0*z, 0.0));
      const std::complex<double> zc(z, 0.0);
      result = std::real(2.0*zc/y*(dilog(1.0 - 0.5*(1.0 - y)/zc) - dilog(1.0 - 0.5*(1.0 + y)/zc)));
   }

   return result;
}

/**
 * Calculates \f$f_S(z)\f$, Eq (71) arXiv:hep-ph/0609168
 */
double f_S(double z) {
   if (z < 0.0) {
      ERROR("f_S: z must not be negative!");
   }

   return (2.0*z - 1.0)*f_PS(z) - 2.0*z*(2.0 + std::log(z));
}

/**
 * Calculates \f$f_{\tilde{f}}(z)\f$, Eq (72) arXiv:hep-ph/0609168
 */
double f_sferm(double z) {
   if (z < 0.0) {
      ERROR("f_sferm: z must not be negative!");
   }

   return 0.5*z*(2.0 + std::log(z) - f_PS(z));
}

} // namespace gm2calc
