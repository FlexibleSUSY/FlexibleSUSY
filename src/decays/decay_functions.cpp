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

#include <cmath>
#include <limits>

#include "decay_functions.hpp"
#include "dilog.hpp"
#include "trilog.hpp"
#include "Li4.hpp"
#include "numerics2.hpp"
#include "wrappers.hpp"

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

/// Eq.(2.7) of hep-ph/0503173
double calc_A(double b) noexcept
{
   const double log_b {std::log(b)};
   const double log_ratio {std::log((1 + b) / (1 - b))};

   return (1 + b * b) *
             (4. * dilog((1 - b) / (1 + b)) + 2. * dilog(-(1 - b) / (1 + b)) -
              3. * log_ratio * std::log(2 / (1 + b)) -
              2. * log_ratio * log_b) -
          3. * b * std::log(4. / (1 - b * b)) - 4. * b * log_b;
}

} // anonymous namespace

/// Eq.(2.6) of hep-ph/0503173
double calc_DeltaH(double b) noexcept
{
   const double b2 = b*b;

   return calc_A(b)/b
      + 1./(16*b2*b) * (3 + b2*(34 - 13*b2)) * std::log((1+b)/(1-b))
      + 3./(8.*b2)*(7*b2 - 1);
}

// Eq.(2.6 of) hep-ph/0503173
double calc_DeltaAH(double b) noexcept
{
   const double b2 = b*b;

   return calc_A(b)/b + 1./(16*b) * (19 + b2*(2 + 3*b2))
      * std::log((1+b)/(1-b)) + 3./8.*(7 - b2);
}

/// Eq.(2.11) of hep-ph/0503173, 2-loop and higher order
double calc_Deltaqq(double alpha_s_red, double Nf) noexcept
{
   // order alpha_s_red^1 is taken into account with mass dependence somewhere else
   return
      + alpha_s_red*(0. * 17./3.
      + alpha_s_red*(35.94 - 1.36*Nf
      + alpha_s_red*(164.14 + Nf*(-25.77 + 0.259*Nf)
      + alpha_s_red*(39.34 + Nf*(-220.9 + Nf*(9.685 - 0.0205*Nf))))));
}

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

std::complex<double> f(double tau) noexcept {
   if (tau <= 1) {
      return Sqr(std::asin(std::sqrt(tau)));
   }
   else {
      using namespace std::literals::complex_literals;
      return -0.25*Sqr(std::log((1.+std::sqrt(1.+1./tau))/(1.+std::sqrt(1.-1./tau))) - 1i*Pi);
   }
}

std::complex<double> taum1fprime(double tau) noexcept {
   if (is_zero(tau-1)) {
      return 0.;
   }
   else if (tau < 1) {
      return (tau-1.)*std::asin(std::sqrt(tau))/std::sqrt(tau - Sqr(tau));
   }
   else {
      using namespace std::literals::complex_literals;
      return
(tau-1.)*(1i*std::log(-1. - 1i*Pi + 2.*tau + 2.*std::sqrt((-1. + tau)*tau)))/
   (-2i*std::sqrt((-1. + tau)*tau) + 
     2*Pi*(2*(-1 + tau)*tau + std::sqrt((-1 + tau)*tau) - 2*std::sqrt((-1 + tau)*Power3(tau))));
   }
}

std::complex<double> fprime(double tau) noexcept {
   if (tau < 1) {
      return std::asin(std::sqrt(tau))/std::sqrt(tau - Sqr(tau));
   }
   else {
      using namespace std::literals::complex_literals;
      return
(1i*std::log(-1. - 1i*Pi + 2.*tau + 2.*std::sqrt((-1. + tau)*tau)))/
   (-2i*std::sqrt((-1. + tau)*tau) + 
     2*Pi*(2*(-1 + tau)*tau + std::sqrt((-1 + tau)*tau) - 2*std::sqrt((-1 + tau)*Power3(tau))));
   }
}

// eq. 2.5, 2.7 & 2.8 of https://arxiv.org/pdf/hep-ph/0509189.pdf
std::complex<double> delta_hAA_2loopQCD_for_quark_loop(double mH, double mq, double mu) noexcept
{
   const double tau = Sqr(mH/(2.*mq));
   const std::complex<double> z = is_zero(tau) ? 1. : (std::sqrt(std::complex<double>(1. - 1./tau))-1.)/(std::sqrt(std::complex<double>(1.-1./tau))+1.);
   const std::complex<double> ln = std::log(z);
   const std::complex<double> li2p = dilog(z);
   const std::complex<double> li2m = dilog(-z);
   const std::complex<double> li3p = trilog(z);
   const std::complex<double> li3m = trilog(-z);
   const std::complex<double> li4p = Li4(z);
   const std::complex<double> li4m = Li4(-z);

   const std::complex<double> F0H = 3./(2.*Sqr(tau))*(tau+(tau-1.)*f(tau));

   const std::complex<double> F0HC1H =
      -z*(1. + z + Sqr(z) + Cube(z))/Power5(1.-z) * (108.*li4p + 144.*li4m - 64.*li3p*ln
            - 64.*li3m*ln + 14.*li2p*Sqr(ln) + 8.*li2m*Sqr(ln) + 1./12.*Power4(ln)
            + 4.*zeta2*Sqr(ln) + 16.*zeta3*ln + 18.*zeta4)
      + z*Sqr(1+z)/Power4(1.-z) * (-32.*li3m + 16.*li2m*ln - 4.*zeta2*ln)
      - 4.*z*(7.-2*z+7*Sqr(z))/Power4(1.-z)*li3p + 8.*z*(3.-2*z+3*Sqr(z))/Power4(1.-z)*li2p*ln
      + 2.*z*(5.-6.*z+5.*Sqr(z))/Power4(1.-z)*std::log(1.-z)*Sqr(ln) + z*(3.+25.*z-7.*Sqr(z)+3.*Cube(z))/(3.*Power5(1.-z))*Cube(ln)
      + 4.*z*(1.-14.*z+Sqr(z))/Power4(1.-z)*zeta3 + 12.*Sqr(z)/Power4(1.-z)*Sqr(ln) - 12.*z*(1.+z)/Power3(1.-z)*ln - 20.*z/Sqr(1.-z);
   const std::complex<double> F0HC2H = 3./Sqr(tau)*(tau+(tau-2.)*f(tau) - tau*taum1fprime(tau));

   return (F0HC1H + F0HC2H*std::log(Sqr(mu/mq)))/F0H;
}

// eq. 2.17, 2.19 & 2.20 of https://arxiv.org/pdf/hep-ph/0509189.pdf
std::complex<double> delta_AhAA_2loopQCD_for_quark_loop(double mAh, double mq, double mu) noexcept {
   const double tau = Sqr(mAh/(2.*mq));
   const std::complex<double> z = is_zero(tau) ? 1. : (std::sqrt(std::complex<double>(1. - 1./tau))-1.)/(std::sqrt(std::complex<double>(1.-1./tau))+1.);
   const std::complex<double> ln = std::log(z);
   const std::complex<double> li2p = dilog(z);
   const std::complex<double> li2m = dilog(-z);
   const std::complex<double> li3p = trilog(z);
   const std::complex<double> li3m = trilog(-z);
   const std::complex<double> li4p = Li4(z);
   const std::complex<double> li4m = Li4(-z);

   const std::complex<double> F0A = f(tau)/tau;

   const std::complex<double> F0AC1A =
      - z*(1.+Sqr(z))/(Cube(1.-z)*(1.+z))*(72.*li4p + 96.*li4m - 128./3.*(li3p + li3m)*ln
            + 28./3.*li2p*Sqr(ln) + 16./3.*li2m*Sqr(ln) + 1./18.*Power4(ln)
            + 8./3.*zeta2*Sqr(ln) + 32./3.*zeta3*ln + 12.*zeta4)
      + z/Sqr(1.-z)*(-56./.3*li3p - 64./3.*li3m + 16.*li2p*ln
            + 32./3*li2m*ln + 20./3.*std::log(1.-z)*Sqr(ln) - 8./3.*zeta2*ln + 8./3.*zeta3)
      + 2.*z*(1.+z)/(3.*Cube(1.-z))*Cube(ln);
   const std::complex<double> F0AC2A = 2./tau*(f(tau) - tau*fprime(tau));

   return (F0AC1A + F0AC2A*std::log(Sqr(mu/mq)))/F0A;
}

/* 2-loop QCD corrections to H->gamma gamma amplitude through scalar color triplet loop
 * from hep-ph/0611266 */

/* this functions could be improved to cover all range in r by linking
 * to CHAPLIN library (https://arxiv.org/abs/1106.5739) */
std::complex<double> delta_hAA_2loopQCD_for_squark_loop(double mH, double msq, double mu) noexcept {
   const std::complex<double> r = Sqr(mH/msq);
   if (std::abs(r) < 0.7) {
      // eq. A.3
      const std::complex<double> F01l = -1./3. - 2./45*r - 1./140.*Sqr(r) - 2./1575*Cube(r) - 1./4158*Power4(r);
      const std::complex<double> F02la = -3./4. - 29./216.*r - 4973./226800.*Sqr(r) - 3137./882000.*Cube(r) - 1180367./2095632000.*Power4(r);
      const std::complex<double> F02lb = -1./4. - 1./15.*r - 9./560.*Sqr(r) - 2./525.*Cube(r) - 5./5544.*Power4(r);
      const std::complex<double> F02lc = - 3/.4*F01l;
      return (F02la + (F02lb + F02lc)*std::log(Sqr(msq/mu)))/F01l;
   }
   else if (std::abs(r) > 1.5) {
      // eq. A.6
      const std::complex<double> F01l = 4./r + Sqr(2.*std::log(-r)/r);
      const std::complex<double> F02la =
         (14.-3*std::log(-r))/r
         + (6. - 72./5*Sqr(zeta2) - 24.*zeta3 - 8.*(1.-zeta2-4.*zeta3)*std::log(-r)
               +(17.-8*zeta2)*Sqr(std::log(-r)) - 5./3.*Cube(std::log(-r)) - 1./6.*Power4(std::log(-r))
            )/Sqr(r);
      const std::complex<double> F02lb =
         (6.*std::log(-r) - 3.*Sqr(std::log(-r)))/Sqr(r);
      const std::complex<double> F02lc = - 3/.4*F01l;
      return (F02la + (F02lb + F02lc)*std::log(Sqr(msq/mu)))/F01l;
   }
   else {
      return 0.;
   }
}

std::complex<double> delta_AhAA_2loopQCD_for_squark_loop(double mAH, double msq, double mu) noexcept {
   const double r = Sqr(mAH/msq);
   if (r < 0.7) {
      return 0.;
   }
   else if (r > 1.5) {
      return 0.;
   }
   else {
      return 0.;
   }
}

// eq. 6 of https://arxiv.org/pdf/1109.5304.pdf
std::complex<double> hgg_SM_loop_function(double x) noexcept {
   if(is_zero(1.-x)) {
      return Sqr(M_PI/2.);
   }
   else if (x > 1) {
      return Sqr(std::asin(1./std::sqrt(x)));
   }
   else {
      return -0.25*Sqr(std::log((1.+std::sqrt(1.-x))/(1.-std::sqrt(1.-x))) - M_PI*std::complex<double>(0.,1.));
   }
}

unsigned int number_of_active_flavours(double m) noexcept
{
   /*
   if (m < qedqcd.displayMbMb(){
      return 4.0;
   } else if (m < qedqcd.displayPoleMt()) {
      return 5.0;
   } else {
      return 6.0;
   }
   */
   return 5;
}

} // namespace flexiblesusy
