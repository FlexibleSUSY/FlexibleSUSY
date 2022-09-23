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

#include "sm_twoloop_mt.hpp"

#include "config.h"

#ifdef ENABLE_TSIL

#include <tsil_cpp.h>
#include <cmath>

/**
 * @file sm_twoloop_mt.cpp
 * @brief contains 2-loop corrections to mt(MS-bar,SM) from [arXiv:1604.01134]
 * @author Harun Acaroglu
 */

namespace flexiblesusy {
namespace sm_twoloop_mt {

namespace {

const auto zeta3 = 1.2020569031595942853997381615114L;
const auto sqrt2 = 1.4142135623730950488016887242097L;
const auto PI = 3.1415926535897932384626433832795L;
const auto PI2 = PI*PI;
const auto k  = 1.0L/(16.0L*PI2);
const auto k2 = k*k;

#ifdef TSIL_SIZE_DOUBLE

/*
 * The following operators need to be defined, if TSIL was compiled
 * with precision -DTSIL_SIZE_DOUBLE instead of -DTSIL_SIZE_LONG.
 */

std::complex<double> operator*(long double a, std::complex<double> b)
{
   return {static_cast<double>(a) * b.real(),
           static_cast<double>(a) * b.imag()};
}

std::complex<double> operator*(std::complex<double> b, long double a)
{
   return {static_cast<double>(a) * b.real(),
           static_cast<double>(a) * b.imag()};
}

std::complex<double> operator/(std::complex<double> b, long double a)
{
   return {b.real() * (1./static_cast<double>(a)),
           b.imag() * (1./static_cast<double>(a))};
}

std::complex<double> operator+(long double a, std::complex<double> b)
{
   return static_cast<double>(a) + b;
}

std::complex<double> operator+(std::complex<double> b, long double a)
{
   return static_cast<double>(a) + b;
}

std::complex<double> operator-(long double a, std::complex<double> b)
{
   return static_cast<double>(a) - b;
}

std::complex<double> operator-(std::complex<double> b, long double a)
{
   return b - static_cast<double>(a);
}

#endif

/**
 * Returns \f$\frac{1}{t} \cdot \delta^{(1)}_{\text{QCD}}\f$ in Eq.(2.1)
 from [arxiv:1604.01134]
 *
 * @param t squared running top mass \f$m^2_t\f$
 * @param qq renormalization scale \f$Q^2\f$
 *
 * @return \f$\frac{1}{t} \cdot \delta^{(1)}_{\text{QCD}}\f$
 */
TSIL_COMPLEXCPP delta1QCD(TSIL_REAL t,TSIL_REAL qq)
{
   const TSIL_COMPLEXCPP a = 32.0L/3.0L
      -8.0L*(TSIL_A_(t,qq)*(1.0L/t)+1.0L);

   return a;
}

/**
 * Returns gaugeless limit of the scalar part of the top self-energy
 \f$\frac{\Sigma^{(s)}_t}{(4 \pi)^2}\f$ in Eq.(7) from [arxiv:1710.03760]
 *
 * @param t squared running top mass \f$m^2_t\f$
 * @param h squared running Higgs mass \f$m^2_h\f$
 * @param yt running top Yukawa coupling \f$y_t\f$
 * @param s external momentum invariant \f$p^2\f$
 * @param qq renormalization scale \f$Q^2\f$
 *
 * \f$\frac{\Sigma^{(s)}}{(4 \pi)^2}|_{g=g'=0}\f$
 */
TSIL_COMPLEXCPP SigmaS(TSIL_REAL t, TSIL_REAL h, TSIL_REAL yt,
   TSIL_REAL s, TSIL_REAL qq)
{
   const auto mt = std::sqrt(t);
   const auto v = mt / yt * sqrt2;
   const auto v2 = v*v;
   const auto t32 = t*std::sqrt(t);
   const auto Bht = TSIL_B_(h,t,s,qq);
   const auto Bt0 = TSIL_B_(t,0.0L,s,qq);

   const TSIL_COMPLEXCPP a = ((Bht - Bt0)*t32)/v2;

   return a;
}

/**
 * Returns the derivative of SigmaS with respect to the external
 momentum invariant s
 *
 * @param t squared running top mass \f$m^2_t\f$
 * @param h squared running Higgs mass \f$m^2_h\f$
 * @param yt running top Yukawa coupling \f$y_t\f$
 * @param s external m,omentum invariant \f$p^2\f$
 * @param qq renormalization scale \f$Q^2\f$
 *
 * \f$\frac{1}{(4 \pi)^2} \cdot \frac{d\Sigma^{(s)}|_{g=g'=0}}{ds}\f$
 */
TSIL_COMPLEXCPP dSigmaSds(TSIL_REAL t, TSIL_REAL h, TSIL_REAL yt,
   TSIL_REAL s, TSIL_REAL qq)
{
   const auto mt = std::sqrt(t);
   const auto v = mt / yt * sqrt2;
   const auto v2 = v*v;
   const auto t32 = t*std::sqrt(t);
   const auto dBhtds = TSIL_dBds_(h,t,s,qq);
   const auto dBt0ds = TSIL_dBds_(t,0.0L,s,qq);

   const TSIL_COMPLEXCPP a = ((dBhtds - dBt0ds)*t32)/v2;

   return a;
}

/**
 * Returns gaugeless limit of the right-handed part of the top
 self-energy \f$\frac{\Sigma^{(R)}_t}{(4 \pi)^2}\f$ in Eq.(7)
 from [arxiv:1710.03760]
 *
 * @param t squared running top mass \f$m^2_t\f$
 * @param h squared running Higgs mass \f$m^2_h\f$
 * @param yt running top Yukawa coupling \f$y_t\f$
 * @param s external m,omentum invariant \f$p^2\f$
 * @param qq renormalization scale \f$Q^2\f$
 *
 * \f$\frac{\Sigma^{(R)}_t}{(4 \pi)^2}|_{g=g'=0}\f$
 */
TSIL_COMPLEXCPP SigmaR(TSIL_REAL t, TSIL_REAL h, TSIL_REAL yt,
   TSIL_REAL s, TSIL_REAL qq)
{
   const auto mt = std::sqrt(t);
   const auto v = mt / yt * sqrt2;
   const auto v2 = v*v;
   const auto At = TSIL_A_(t,qq);
   const auto Ah = TSIL_A_(h,qq);
   const auto Bht = TSIL_B_(h,t,s,qq);
   const auto Bt0 = TSIL_B_(t,0.0L,s,qq);

   const TSIL_COMPLEXCPP a = (t*(-Ah + 2.0L*At - Bht*h + Bht*s + Bt0*s
      + Bht*t + Bt0*t))/(4.0L*s*v2);

   return a;
}

/**
 * Returns the derivative of SigmaR with respect to the external
 momentum invariant s
 *
 * @param t squared running top mass \f$m^2_t\f$
 * @param h squared running Higgs mass \f$m^2_h\f$
 * @param yt running top Yukawa coupling \f$y_t\f$
 * @param s external m,omentum invariant \f$p^2\f$
 * @param qq renormalization scale \f$Q^2\f$
 *
 * \f$\frac{1}{(4 \pi)^2} \cdot \frac{d\Sigma^{(R)}|_{g=g'=0}}{ds}\f$
 */
TSIL_COMPLEXCPP dSigmaRds(TSIL_REAL t, TSIL_REAL h, TSIL_REAL yt,
   TSIL_REAL s, TSIL_REAL qq)
{
   const auto mt = std::sqrt(t);
   const auto v = mt / yt * sqrt2;
   const auto v2 = v*v;
   const auto s2 = s*s;
   const auto At = TSIL_A_(t,qq);
   const auto Ah = TSIL_A_(h,qq);
   const auto Bht = TSIL_B_(h,t,s,qq);
   const auto Bt0 = TSIL_B_(t,0.0L,s,qq);
   const auto dBhtds = TSIL_dBds_(h,t,s,qq);
   const auto dBt0ds = TSIL_dBds_(t,0.0L,s,qq);

   const TSIL_COMPLEXCPP a = (t*(Ah - 2.0L*At + Bht*h - Bht*t - Bt0*t
      + s*(-h + s + t)*dBhtds + s*(s + t)*dBt0ds))/(4.0L*s2*v2);

   return a;
}

/**
 * Returns gaugeless limit of the left-handed part of the top
 self-energy \f$\frac{\Sigma^{(L)}_t}{(4 \pi)^2}\f$ in Eq.(7)
 from [arxiv:1710.03760]
 *
 * @param t squared running top mass \f$m^2_t\f$
 * @param h squared running Higgs mass \f$m^2_h\f$
 * @param yt running top Yukawa coupling \f$y_t\f$
 * @param s external m,omentum invariant \f$p^2\f$
 * @param qq renormalization scale \f$Q^2\f$
 *
 * \f$\frac{\Sigma^{(L)}_t}{(4 \pi)^2}|_{g=g'=0}\f$
 */
TSIL_COMPLEXCPP SigmaL(TSIL_REAL t, TSIL_REAL h, TSIL_REAL yt,
   TSIL_REAL s, TSIL_REAL qq)
{
   const auto mt = std::sqrt(t);
   const auto v = mt / yt * sqrt2;
   const auto v2 = v*v;
   const auto At = TSIL_A_(t,qq);
   const auto Ah = TSIL_A_(h,qq);
   const auto Bht = TSIL_B_(h,t,s,qq);
   const auto Bt0 = TSIL_B_(t,0.0L,s,qq);
   const auto B00 = TSIL_B_(0.0L,0.0L,s,qq);

   const TSIL_COMPLEXCPP a = (t*(-Ah + 2.0L*At - Bht*h + 2.0L*B00*s
      + Bht*s + Bt0*s + Bht*t + Bt0*t))/(4.0L*s*v2);

   return a;
}

/**
 * Returns the derivative of SigmaL with respect to the external
 momentum invariant s
 *
 * @param t squared running top mass \f$m^2_t\f$
 * @param h squared running Higgs mass \f$m^2_h\f$
 * @param yt running top Yukawa coupling \f$y_t\f$
 * @param s external m,omentum invariant \f$p^2\f$
 * @param qq renormalization scale \f$Q^2\f$
 *
 * \f$\frac{1}{(4 \pi)^2} \cdot \frac{d\Sigma^{(L)}|_{g=g'=0}}{ds}\f$
 */
TSIL_COMPLEXCPP dSigmaLds(TSIL_REAL t, TSIL_REAL h, TSIL_REAL yt,
   TSIL_REAL s, TSIL_REAL qq)
{
   const auto mt = std::sqrt(t);
   const auto v = mt / yt * sqrt2;
   const auto v2 = v*v;
   const auto s2 = s*s;
   const auto At = TSIL_A_(t,qq);
   const auto Ah = TSIL_A_(h,qq);
   const auto Bht = TSIL_B_(h,t,s,qq);
   const auto Bt0 = TSIL_B_(t,0.0L,s,qq);
   const auto dBhtds = TSIL_dBds_(h,t,s,qq);
   const auto dBt0ds = TSIL_dBds_(t,0.0L,s,qq);
   const auto dB00ds = TSIL_dBds_(0.0L,0.0L,s,qq);

   const TSIL_COMPLEXCPP a = (t*(Ah - 2.0L*At + Bht*h - Bht*t - Bt0*t
      + 2.0L*s2*dB00ds + s*(-h + s + t)*dBhtds + s2*dBt0ds
      + s*t*dBt0ds))/(4.0L*s2*v2);

   return a;
}

/**
 * Returns \f$\frac{1}{t} \cdot \delta^{(2)}_{\text{QCD}}\f$ in Eq.(2.1)
 from [arxiv:1604.01134]
 *
 * @param t squared running top mass \f$m^2_t\f$
 * @param qq renormalization scale \f$Q^2\f$
 *
 * @return \f$\frac{1}{t} \cdot \delta^{(2)}_{\text{QCD}}\f$
 */
TSIL_COMPLEXCPP delta2QCD(TSIL_REAL t, TSIL_REAL qq)
{
   const auto t2 = t*t;
   const auto At = TSIL_A_(t,qq);
   const auto At2 = At*At;

   const TSIL_COMPLEXCPP a = 112.55555555555556L + (16.0L*PI2)/9.0L
   + (60.0L*At2)*(1.0L/t2) - 84.0L*At*(1.0L/t)
   + (32.0L*PI2*std::log(2.0L))/9.0L - (16.0L*zeta3)/3.0L;

   return a;
}

/**
 * Returns gaugeless limit of \f$\frac{1}{t} \cdot
 \delta^{(2)}_{\text{mixed}}\f$ in Eq.(2.1) from [arxiv:1604.01134]
 *
 * @param t squared running top mass \f$m^2_t\f$
 * @param h squared running Higgs mass \f$m^2_h\f$
 * @param yt running top Yukawa coupling \f$y_t\f$
 * @param qq renormalization scale \f$Q^2\f$
 *
 * @return \f$\frac{1}{t} \cdot \delta^{(2)}_{\text{mixed}}|_{g=g'=0}\f$
 */
TSIL_COMPLEXCPP delta2mixed(TSIL_REAL t, TSIL_REAL h, TSIL_REAL yt,
   TSIL_REAL qq)
{
   TSIL_DATA           result;
   TSIL_SetParameters(&result,0.0L,t,t,h,t,qq);
   TSIL_Evaluate(&result,t);
   const auto Uthtt = TSIL_GetFunction_(&result, "Uyuzv");
   const auto M0ttht = TSIL_GetFunction_(&result, "M");

   TSIL_COMPLEXCPP     Tbar0ht,Th0t,M00t00;
   TSIL_Tbaranalytic_(0.0L,h,t,t,qq,&Tbar0ht);
   TSIL_Tanalytic_(h,0.0L,t,t,qq,&Th0t);
   TSIL_Manalytic_(0.0L,0.0L,t,0.0L,0.0L,t,&M00t00);

   const auto mt = std::sqrt(t);
   const auto v = mt / yt * sqrt2;
   const auto v2 = v*v;
   const auto t2 = t*t;
   const auto h2 = h*h;
   const auto h3 = h2*h;
   const auto At = TSIL_A_(t,qq);
   const auto Ah = TSIL_A_(h,qq);
   const auto At2 = At*At;
   const auto Bht = TSIL_B_(h,t,t,qq);
   const auto Ihtt = TSIL_I2_(h,t,t,qq);

   const TSIL_COMPLEXCPP a = (-24.0L*At2*h*(h - t)
      + 24.0L*At*t*(Ah*(3.0L*h + 5.0L*t) + h*((5.0L + 2.0L*Bht)*h
      + 2.0L*(-12.0L + Bht)*t)) + 2.0L*t*(-60.0L*Ah*t2
      - 2.0L*h2*(3.0L*Ihtt + t*(-6.0L - 51.0L*Bht - 8.0L*PI2
      + 72.0L*M0ttht*t + 12.0L*Tbar0ht - 24.0L*Th0t))
      - 6.0L*h3*(1.0L + Bht - 4.0L*M0ttht*t + Th0t) + h*t*(-54.0L*Ah
      + 84.0L*Ihtt + t*(333.0L - 288.0L*Bht - 61.0L*PI2
      + 24.0L*M00t00*t + 192.0L*M0ttht*t + 96.0L*Tbar0ht - 60.0L*Th0t
      - 96.0L*Uthtt))))/(9.0L*h*t2*v2);

   return a;
}

/**
 * Returns gaugeless limit of \f$\frac{1}{t} \cdot
 \delta^{(2)}_{\text{non-QCD}}\f$ in Eq.(2.1) from [arxiv:1604.01134]
 *
 * @param t squared running top mass \f$m^2_t\f$
 * @param h squared running Higgs mass \f$m^2_h\f$
 * @param yt running top Yukawa coupling \f$y_t\f$
 * @param qq renormalization scale \f$Q^2\f$
 *
 * @return \f$\frac{1}{t} \cdot \delta^{(2)}_{\text{non-QCD}}|_{g=g'=0}\f$
 */
TSIL_COMPLEXCPP delta2Higgs(TSIL_REAL t, TSIL_REAL h, TSIL_REAL yt,
   TSIL_REAL qq)
{
   TSIL_DATA           result1, result2, result3, result4;
   const TSIL_REAL     EPS = 1.0L;

   TSIL_SetParameters(&result1,h,0.0L,t,t,0.0L,qq);
   TSIL_Evaluate(&result1,t);
   const auto Mh0tt0 = TSIL_GetFunction_(&result1, "M");

   TSIL_SetParameters(&result2,h,h,t,t,h,qq);
   TSIL_Evaluate(&result2,t);
   const auto Thht = TSIL_GetFunction_(&result2, "Tvyz");
   const auto Uthhh = TSIL_GetFunction_(&result2, "Uzxyv");
   const auto Mhhtth = TSIL_GetFunction_(&result2, "M");

   TSIL_SetParameters(&result3,h,t,t,h,t,qq);
   TSIL_Evaluate(&result3,t);
   const auto Uhtht = TSIL_GetFunction_(&result3, "Uxzuv");
   const auto Uthtt = TSIL_GetFunction_(&result3, "Uzxyv");
   const auto Mhttht = TSIL_GetFunction_(&result3, "M");

   TSIL_SetParameters(&result4,t,t,0.0L,0.0L,h,qq);
   TSIL_Evaluate(&result4,t);
   const auto Mtt00h = TSIL_GetFunction_(&result4, "M");

   TSIL_COMPLEXCPP   S0h0, Tbar000, Tbar00h, Tbar0ht, Tbar0t0, Th00,
                     Th0t, Tht0, Tth0, U0000, U000h, U00h0, U000t,
                     U0tht, U0tt0, Uht00, Uhtt0, Ut0h0, Uth00, M00t00,
                     M0t0h0a,M0t0h0b;
   TSIL_Sanalytic_(0.0L,h,0.0L,t,qq,&S0h0);
   TSIL_Tbaranalytic_(0.0L,0.0L,0.0L,t,qq,&Tbar000);
   TSIL_Tbaranalytic_(0.0L,0.0L,h,t,qq,&Tbar00h);
   TSIL_Tbaranalytic_(0.0L,0.0L,h,t,qq,&Tbar00h);
   TSIL_Tbaranalytic_(0.0L,h,t,t,qq,&Tbar0ht);
   TSIL_Tbaranalytic_(0.0L,t,0.0L,t,qq,&Tbar0t0);
   TSIL_Tanalytic_(h,0.0L,0.0L,t,qq,&Th00);
   TSIL_Tanalytic_(h,0.0L,t,t,qq,&Th0t);
   TSIL_Tanalytic_(h,t,0.0L,t,qq,&Tht0);
   TSIL_Tanalytic_(t,h,0.0L,t,qq,&Tth0);
   TSIL_Uanalytic_(0.0L,0.0L,0.0L,0.0L,t,qq,&U0000);
   TSIL_Uanalytic_(0.0L,0.0L,0.0L,h,t,qq,&U000h);
   TSIL_Uanalytic_(0.0L,0.0L,h,0.0L,t,qq,&U00h0);
   TSIL_Uanalytic_(0.0L,0.0L,0.0L,t,t,qq,&U000t);
   TSIL_Uanalytic_(0.0L,t,h,t,t,qq,&U0tht);
   TSIL_Uanalytic_(0.0L,t,t,0.0L,t,qq,&U0tt0);
   TSIL_Uanalytic_(h,t,0.0L,0.0L,t,qq,&Uht00);
   TSIL_Uanalytic_(h,t,t,0.0L,t,qq,&Uhtt0);
   TSIL_Uanalytic_(t,0.0L,h,0.0L,t,qq,&Ut0h0);
   TSIL_Uanalytic_(t,h,0.0L,0.0L,t,qq,&Uth00);
   TSIL_Manalytic_(0.0L,0.0L,t,0.0L,0.0L,t,&M00t00);

   TSIL_Manalytic_(0.0L,t,0.0L,h,0.0L,t-EPS,&M0t0h0a);
   TSIL_Manalytic_(0.0L,t,0.0L,h,0.0L,t+EPS,&M0t0h0b);
   const auto M0t0h0 = (M0t0h0a+M0t0h0b)/2.0L;

   const auto mt = std::sqrt(t);
   const auto v = mt / yt * sqrt2;

   const auto q = 0.5L*std::log(qq);
   const auto q2 = q*q;
   const auto v4 = v*v*v*v;
   const auto h2 = h*h;
   const auto h3 = h2*h;
   const auto h4 = h2*h2;
   const auto t2 = t*t;
   const auto t3 = t2*t;
   const auto t4 = t2*t2;
   const auto Logt = std::log(t);
   const auto Logh = std::log(h);
   const auto Logtmh = std::log(t-h);
   const auto Log1mth = std::log(t/h-1.0L);
   const auto Logt2 = Logt*Logt;
   const auto Logh2 = Logh*Logh;
   const auto Ah = TSIL_A_(h,qq);
   const auto At = TSIL_A_(t,qq);
   const auto At2 = At*At;
   const auto Ah2 = Ah*Ah;
   const auto B00 = TSIL_B_(0.0L,0.0L,t,qq);
   const auto Bt0 = TSIL_B_(t,0.0L,t,qq);
   const auto B0h = TSIL_B_(0.0L,h,t,qq);
   const auto Bht = TSIL_B_(h,t,t,qq);
   const auto B002 = B00*B00;
   const auto Bht2 = Bht*Bht;
   const auto Ih00 = TSIL_I2_(h,0.0L,0.0L,qq);
   const auto I0t0 = TSIL_I2_(0.0L,t,0.0L,qq);
   const auto Ihhh = TSIL_I2_(h,h,h,qq);
   const auto Ihtt = TSIL_I2_(h,t,t,qq);
   const auto Diloga = TSIL_Dilog_(t/h);
   const auto Dilogb = TSIL_Dilog_(t/(t-h));

   const TSIL_COMPLEXCPP a = (-72.0L*Bht*h3 + 12.0L*h3*PI2
      - 336.0L*h2*t + 528.0L*Bht*h2*t - 60.0L*Bht2*h2*t
      + 168.0L*h*Ih00*t - 288.0L*h*Ihhh*t + 96.0L*h*Ihtt*t
      - 52.0L*h2*PI2*t + 96.0L*h2*q*t + 96.0L*B0h*h2*q*t
      - 336.0L*h*S0h0*t + 258.0L*h*t2 - 768.0L*Bht*h*t2
      - 72.0L*B00*Bht*h*t2 + 336.0L*Bht2*h*t2 + 96.0L*Bt0*h*t2
      - 96.0L*Bht*Bt0*h*t2 + 72.0L*Ih00*t2 + 336.0L*I0t0*t2
      - 960.0L*Ihtt*t2 - 96.0L*h2*Mh0tt0*t2 + 432.0L*h2*Mhhtth*t2
      + 48.0L*h2*Mhttht*t2 - 48.0L*h2*Mtt00h*t2 + 64.0L*h*PI2*t2
      - 288.0L*h*q*t2 - 96.0L*B00*h*q*t2 + 96.0L*B0h*h*q*t2
      + 288.0L*Bt0*h*q*t2 + 72.0L*S0h0*t2 + 1131.0L*t3 - 444.0L*B00*t3
      + 96.0L*B002*t3 - 2112.0L*Bht*t3 + 48.0L*B00*Bht*t3
      - 384.0L*Bht2*t3 + (1728.0L*Ihtt*t3)/h - 192.0L*h*M0t0h0*t3
      - 1152.0L*h*Mhhtth*t3 + 384.0L*h*Mhttht*t3 - 96.0L*PI2*t3
      - 144.0L*q*t3 + 144.0L*Bt0*q*t3 + (1152.0L*t4)/h
      + (2304.0L*Bht*t4)/h - 768.0L*Mhttht*t4
      + (36.0L*Ah2*(h2 - 6.0L*h*t - 5.0L*t2))/h + (96.0L*At2*(h2
      - 18.0L*h*t + 42.0L*t2))/h - (24.0L*At*(6.0L*(-1.0L + Bht)*h3
      - 4.0L*h2*(-7.0L + 8.0L*Bht - Bt0 + 3.0L*q)*t + h*(-27.0L
      - 10.0L*B00 + 64.0L*Bht + 18.0L*q)*t2 + 48.0L*(3.0L
      - 2.0L*Bht)*t3))/h - (24.0L*Ah*(3.0L*h3 + h2*(-35.0L + 2.0L*B00
      - 4.0L*Bht + Bt0 - 4.0L*q)*t + (-25.0L - 4.0L*B00
      + 4.0L*Bht)*h*t2 + (9.0L + 4.0L*B00 + 32.0L*Bht)*t3
      + At*(7.0L*h2 - 8.0L*h*t + 39.0L*t2)))/h + 84.0L*t3*Tbar000
      + 48.0L*h2*t*Tbar00h + 48.0L*h*t2*Tbar00h + 144.0L*h*t2*Tbar0t0
      + 72.0L*t3*Tbar0t0 - 120.0L*h2*t*Th00 + 168.0L*h*t2*Th00
      - 96.0L*t3*Th00 - 72.0L*h3*Th0t + 384.0L*h2*t*Th0t
      - 384.0L*h*t2*Th0t + 408.0L*h*t2*Thht - 408.0L*t3*Thht
      + 48.0L*h2*t*Th0t + 24.0L*h*t2*Th0t + 24.0L*h2*t*Tth0
      + 48.0L*h*t2*Tth0 - 144.0L*t3*Tth0 - 12.0L*t3*U0000
      + 72.0L*h*t2*U000h + 576.0L*t3*U000t - 48.0L*h2*t*U000h
      + 144.0L*h*t2*U000h - 96.0L*h*t2*U0tht + 96.0L*h*t2*U0tt0
      - 120.0L*h*t2*Uht00 + 240.0L*t3*Uht00 - 72.0L*h2*t*Uhtht
      + 384.0L*t3*Uhtht - 96.0L*h*t2*Uhtt0 - 24.0L*h2*t*Ut0h0
      - 144.0L*h2*t*Uth00 + 528.0L*h*t2*Uth00 - 288.0L*h2*t*Uthhh
      + 432.0L*h*t2*Uthhh - 864.0L*h*t2*Uthtt + 2688.0L*t3*Uthtt
      + (4.0L*(3.0L*h4*(7.0L - 12.0L*PI2) + 3.0L*h3*(22.0L
      + 19.0L*PI2)*t + h2*(33.0L - 23.0L*PI2)*t2 + 2.0L*h*(13.0L*PI2
      - 72.0L*(1.0L+2.0L*q+q2))*t3 - 24.0L*(PI2 - 6.0L*(1.0L+2.0L*q
      + q2))*t4 + 3.0L*h*(h - t)*t*(h + 2.0L*t)*Logh2 - 6.0L*t2*(-2.0L*h
      + 3.0L*t)*(-h + 4.0L*t)*Logt + 36.0L*(h - t)*t3*Logt2
      + 3.0L*h*Logtmh*(2.0L*t*(-5.0L*h2 + h*t + t2)
      - (-h + t)*(-4.0L*h2 + h*t + 2.0L*t2)*Logtmh)
      + 6.0L*h*Logh*(h*t*(h + 5.0L*t) + (-h + t)*(2.0L*t2*Logt
      + h*(-2.0L*h + t)*Logtmh)) + 18.0L*h4*Log1mth
      + 6.0L*h*t*(-h + t)*((-h - 2.0L*t)*Diloga
      - (h - 2.0L*t)*Dilogb)))/(h - t))/(48.0L*t*v4);

   return a;
}

TSIL_COMPLEXCPP deltamt1QCD(TSIL_REAL g3, TSIL_REAL t, TSIL_REAL qq)
{
   const auto g32= g3*g3;

   const auto a = -0.5L*k*delta1QCD(t,qq)*g32;

   return a;
}

TSIL_COMPLEXCPP Sigma1S(TSIL_REAL t, TSIL_REAL h, TSIL_REAL yt,
   TSIL_REAL T, TSIL_REAL qq)
{
   const auto a = SigmaS(t, h, yt, T, qq)*k;

   return a;
}

TSIL_COMPLEXCPP Sigma1R(TSIL_REAL t, TSIL_REAL h, TSIL_REAL yt,
   TSIL_REAL T, TSIL_REAL qq)
{
   const auto a = SigmaR(t, h, yt, T, qq)*k;

   return a;
}

TSIL_COMPLEXCPP Sigma1L(TSIL_REAL t, TSIL_REAL h, TSIL_REAL yt,
   TSIL_REAL T, TSIL_REAL qq)
{
   const auto a = SigmaL(t, h, yt, T, qq)*k;

   return a;
}

TSIL_COMPLEXCPP deltamt2QCD(TSIL_REAL g3, TSIL_REAL t, TSIL_REAL qq)
{
   const auto g34 = g3*g3*g3*g3;
   const auto d1 = delta1QCD(t,qq);

   const auto a = ((3.0L*d1*d1 - 4.0L*delta2QCD(t,qq))*g34*k2)/8.0L;

   return a;
}

TSIL_COMPLEXCPP deltamt2mixed(TSIL_REAL g3, TSIL_REAL t, TSIL_REAL h,
   TSIL_REAL yt, TSIL_REAL T, TSIL_REAL qq)
{
   const auto g32 = g3*g3;
   const auto d1 = delta1QCD(t, qq);

   const auto a = (g32*k2*(-delta2mixed(t, h, yt, qq)
      - 3.0L*d1*(SigmaL(t, h, yt, t, qq)
      + SigmaR(t, h, yt, t, qq))
      - 2.0L*d1*t*(dSigmaRds(t, h, yt, T, qq)
      + dSigmaLds(t, h, yt, T, qq))))/2.0L;

   return a;
}

TSIL_COMPLEXCPP deltamt2Higgs(TSIL_REAL t, TSIL_REAL h, TSIL_REAL yt,
   TSIL_REAL T, TSIL_REAL qq)
{
   const auto sl = SigmaL(t, h, yt, t, qq);
   const auto sr = SigmaR(t, h, yt, t, qq);

   const auto a = (k2*(-delta2Higgs(t, h, yt, qq)
      + 3.0L*sl*sl + 6.0L*sl*sr + 3.0L*sr*sr
      + 4.0L*(SigmaS(t, h, yt, t, qq)*std::sqrt(t)
      + (sr + sl)*t)*(dSigmaRds(t, h, yt, T, qq)
      + dSigmaLds(t, h, yt, T, qq))))/2.0L;

   return a;
}

TSIL_COMPLEXCPP Sigma2Smixed(TSIL_REAL g3, TSIL_REAL t, TSIL_REAL h,
   TSIL_REAL yt, TSIL_REAL T, TSIL_REAL qq)
{
   const auto g32 = g3*g3;

   const auto a = -(delta1QCD(t,qq)*g32*k2*(SigmaS(t, h, yt, t, qq)
      + t*dSigmaSds(t, h, yt, T, qq)));

   return a;
}

TSIL_COMPLEXCPP Sigma2SHiggs(TSIL_REAL t, TSIL_REAL h,
   TSIL_REAL yt, TSIL_REAL T, TSIL_REAL qq)
{
   const auto mt = std::sqrt(t);
   const auto ss = SigmaS(t, h, yt, t, qq);
   const auto sl = SigmaL(t, h, yt, t, qq);
   const auto sr = SigmaR(t, h, yt, t, qq);

   const auto a = (k2*(ss*(ss
      + 4.0L*(sl + sr)*mt)
      + 4.0L*(ss + (sl + sr)*mt)*t*dSigmaSds(t, h, yt, T, qq)))
      /(2.0L*mt);

   return a;
}

/*

TSIL_COMPLEXCPP mt_flexiblesusy(TSIL_REAL g3, TSIL_REAL t, TSIL_REAL h,
   TSIL_REAL yt, TSIL_REAL T, TSIL_REAL qq)
{
   const auto Mt = std::sqrt(T);

   const auto a = Mt + Sigma1S(t, h, yt, T, qq)
      + Sigma2Smixed(g3, t, h, yt, T, qq)
      + Sigma2SHiggs(t, h, yt, T, qq) + Mt*(Sigma1L(t, h, yt, T, qq)
      + Sigma1R(t, h, yt, T, qq) + deltamt1QCD(g3, t, qq)
      + deltamt2QCD(g3, t, qq) + deltamt2mixed(g3, t, h, yt, T, qq)
      + deltamt2Higgs(t, h, yt, T, qq));

   return a;
}

*/

TSIL_COMPLEXCPP deltamt2QCDspheno(TSIL_REAL g3, TSIL_REAL t,
   TSIL_REAL qq)
{
   const auto g34 = g3*g3*g3*g3;
   const auto d1 = delta1QCD(t,qq);

   const auto a = ((d1*d1 - 4.0L*delta2QCD(t,qq))*g34*k2)/8.0L;

   return a;
}

TSIL_COMPLEXCPP deltamt2mixedspheno(TSIL_REAL g3, TSIL_REAL t,
   TSIL_REAL h, TSIL_REAL yt, TSIL_REAL T, TSIL_REAL qq)
{
   const auto g32 = g3*g3;
   const auto mt = std::sqrt(t);
   const auto d1 = delta1QCD(t,qq);

   const auto a = -(g32*k2*(delta2mixed(t, h, yt, qq)
      + d1*SigmaL(t, h, yt, t, qq)
      + d1*SigmaR(t, h, yt, t, qq)
      + 2.0L*d1*t*(dSigmaLds(t, h, yt, T, qq)
      + dSigmaRds(t, h, yt, T, qq))
      + 2.0L*d1*mt*dSigmaSds(t, h, yt, T, qq)))/2.0L;

   return a;
}

TSIL_COMPLEXCPP deltamt2Higgsspheno(TSIL_REAL t, TSIL_REAL h,
   TSIL_REAL yt, TSIL_REAL T, TSIL_REAL qq)
{
   const auto mt = std::sqrt(t);
   const auto sl = SigmaL(t, h, yt, t, qq);
   const auto sr = SigmaR(t, h, yt, t, qq);

   const auto a = (k2*(-delta2Higgs(t, h, yt, qq)
      + (sl + sr)*(sl + sr) + 4.0L*(SigmaS(t, h, yt, t, qq)
      + (sl + sr)*mt)*(mt*(dSigmaLds(t, h, yt, T, qq)
      + dSigmaRds(t, h, yt, T, qq)) + dSigmaSds(t, h, yt, T, qq))))
   /2.0L;

   return a;
}

TSIL_COMPLEXCPP Sigma2Smixedspheno(TSIL_REAL g3, TSIL_REAL t,
   TSIL_REAL h, TSIL_REAL yt, TSIL_REAL qq)
{
   const auto g32 = g3*g3;

   const auto a = -(delta1QCD(t,qq)*g32*k2*SigmaS(t, h, yt, t, qq))
      /2.0L;

   return a;
}

TSIL_COMPLEXCPP Sigma2SHiggsspheno(TSIL_REAL t,
   TSIL_REAL h, TSIL_REAL yt, TSIL_REAL qq)
{
   const auto mt = std::sqrt(t);
   const auto ss = SigmaS(t, h, yt, t, qq);

   const auto a =
      (k2*ss*(2.0L*(SigmaL(t, h, yt, t, qq) + SigmaR(t, h, yt, t, qq))
      + ss/mt))/2.0L;

   return a;
}

/*

TSIL_COMPLEXCPP mt_spheno(TSIL_REAL g3, TSIL_REAL t, TSIL_REAL h,
   TSIL_REAL yt, TSIL_REAL T, TSIL_REAL qq)
{
   const auto Mt = std::sqrt(T);
   const auto mt = std::sqrt(t);

   const auto a = Mt + Sigma1S(t, h, yt, T, qq)
      + Sigma2Smixedspheno(g3, t, h, yt, qq)
      + Sigma2SHiggsspheno(t, h, yt, qq) + mt*(Sigma1L(t, h, yt, T, qq)
      + Sigma1R(t, h, yt, T, qq) + deltamt1QCD(g3, t, qq)
      + deltamt2QCDspheno(g3, t, qq)
      + deltamt2mixedspheno(g3, t, h, yt, T, qq)
      + deltamt2Higgsspheno(t, h, yt, T, qq));

   return a;
}

*/

} // anonymous namespace

/* ******************** 1-loop ******************** */

TSIL_REAL delta_Mt_1loop_as(TSIL_REAL g3, TSIL_REAL t, TSIL_REAL qq)
{
   return std::real(-deltamt1QCD(g3, t, qq));
}

TSIL_REAL delta_Mt_1loop_at(TSIL_REAL t, TSIL_REAL h, TSIL_REAL yt, TSIL_REAL qq)
{
   const auto mt = std::sqrt(t);
   const auto ss = Sigma1S(t, h, yt, t, qq);
   const auto sl = Sigma1L(t, h, yt, t, qq);
   const auto sr = Sigma1R(t, h, yt, t, qq);

   return -std::real(ss + mt*(sl + sr));
}

TSIL_REAL delta_mt_1loop_as(TSIL_REAL g3, TSIL_REAL t, TSIL_REAL qq)
{
   return std::real(deltamt1QCD(g3, t, qq));
}

TSIL_REAL delta_mt_1loop_at_S(TSIL_REAL yt, TSIL_REAL t, TSIL_REAL h, TSIL_REAL s, TSIL_REAL qq)
{
   return std::real(Sigma1S(t, h, yt, s, qq));
}

TSIL_REAL delta_mt_1loop_at_L(TSIL_REAL yt, TSIL_REAL t, TSIL_REAL h, TSIL_REAL s, TSIL_REAL qq)
{
   return std::real(Sigma1L(t, h, yt, s, qq));
}

TSIL_REAL delta_mt_1loop_at_R(TSIL_REAL yt, TSIL_REAL t, TSIL_REAL h, TSIL_REAL s, TSIL_REAL qq)
{
   return std::real(Sigma1R(t, h, yt, s, qq));
}

/* ******************** 2-loop ******************** */

TSIL_REAL delta_Mt_2loop_as_as(TSIL_REAL g3, TSIL_REAL t, TSIL_REAL qq)
{
   const auto g34 = g3*g3*g3*g3;
   const auto d1 = delta1QCD(t,qq);

   const auto a = -d1*d1/8.0L + 0.5L*delta2QCD(t,qq);

   return std::real(k2*g34*a);
}

/// FlexibleSUSY convention ///

TSIL_REAL delta_mt_2loop_as_as_flexiblesusy(
   TSIL_REAL g3, TSIL_REAL t, TSIL_REAL qq)
{
   return std::real(deltamt2QCD(g3, t, qq));
}

TSIL_REAL delta_mt_2loop_as_at_S_flexiblesusy(
   TSIL_REAL g3, TSIL_REAL yt, TSIL_REAL t, TSIL_REAL h, TSIL_REAL s, TSIL_REAL qq)
{
   return std::real(Sigma2Smixed(g3, t, h, yt, s, qq));
}

TSIL_REAL delta_mt_2loop_as_at_LR_flexiblesusy(
   TSIL_REAL g3, TSIL_REAL yt, TSIL_REAL t, TSIL_REAL h, TSIL_REAL s, TSIL_REAL qq)
{
   return std::real(deltamt2mixed(g3, t, h, yt, s, qq));
}

TSIL_REAL delta_mt_2loop_at_at_S_flexiblesusy(
   TSIL_REAL yt, TSIL_REAL t, TSIL_REAL h, TSIL_REAL s, TSIL_REAL qq)
{
   return std::real(Sigma2SHiggs(t, h, yt, s, qq));
}

TSIL_REAL delta_mt_2loop_at_at_LR_flexiblesusy(
   TSIL_REAL yt, TSIL_REAL t, TSIL_REAL h, TSIL_REAL s, TSIL_REAL qq)
{
   return std::real(deltamt2Higgs(t, h, yt, s, qq));
}

/// SPheno convention ///

TSIL_REAL delta_mt_2loop_as_as_spheno(
   TSIL_REAL g3, TSIL_REAL t, TSIL_REAL qq)
{
   return std::real(deltamt2QCDspheno(g3, t, qq));
}

TSIL_REAL delta_mt_2loop_as_at_S_spheno(
   TSIL_REAL g3, TSIL_REAL yt, TSIL_REAL t, TSIL_REAL h, TSIL_REAL /* s */, TSIL_REAL qq)
{
   return std::real(Sigma2Smixedspheno(g3, t, h, yt, qq));
}

TSIL_REAL delta_mt_2loop_as_at_LR_spheno(
   TSIL_REAL g3, TSIL_REAL yt, TSIL_REAL t, TSIL_REAL h, TSIL_REAL s, TSIL_REAL qq)
{
   return std::real(deltamt2mixedspheno(g3, t, h, yt, s, qq));
}

TSIL_REAL delta_mt_2loop_at_at_S_spheno(
   TSIL_REAL yt, TSIL_REAL t, TSIL_REAL h, TSIL_REAL /* s */, TSIL_REAL qq)
{
   return std::real(Sigma2SHiggsspheno(t, h, yt, qq));
}

TSIL_REAL delta_mt_2loop_at_at_LR_spheno(
   TSIL_REAL yt, TSIL_REAL t, TSIL_REAL h, TSIL_REAL s, TSIL_REAL qq)
{
   return std::real(deltamt2Higgsspheno(t, h, yt, s, qq));
}

} // namespace sm_twoloop_mt
} // namespace flexiblesusy

#else // ENABLE_TSIL

#include "error.hpp"

namespace flexiblesusy {
namespace sm_twoloop_mt {

#define ERROR_NO_TSIL(fun)                                      \
   do {                                                         \
      throw SetupError("TSIL is required to call " fun);        \
   } while (false)

/* ******************** 1-loop ******************** */

TSIL_REAL delta_mt_1loop_as(TSIL_REAL, TSIL_REAL, TSIL_REAL)
{
   ERROR_NO_TSIL("delta_mt_1loop_as");
}

TSIL_REAL delta_mt_1loop_at_S(TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL)
{
   ERROR_NO_TSIL("delta_mt_1loop_at_S");
}

TSIL_REAL delta_mt_1loop_at_L(TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL)
{
   ERROR_NO_TSIL("delta_mt_1loop_at_L");
}

TSIL_REAL delta_mt_1loop_at_R(TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL)
{
   ERROR_NO_TSIL("delta_mt_1loop_at_R");
}

/* ******************** 2-loop ******************** */

/// FlexibleSUSY convention ///

TSIL_REAL delta_mt_2loop_as_as_flexiblesusy(
   TSIL_REAL, TSIL_REAL, TSIL_REAL)
{
   ERROR_NO_TSIL("delta_mt_2loop_as_as_flexiblesusy");
}

TSIL_REAL delta_mt_2loop_as_at_S_flexiblesusy(
   TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL)
{
   ERROR_NO_TSIL("delta_mt_2loop_as_at_S_flexiblesusy");
}

TSIL_REAL delta_mt_2loop_as_at_LR_flexiblesusy(
   TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL)
{
   ERROR_NO_TSIL("delta_mt_2loop_as_at_LR_flexiblesusy");
}

TSIL_REAL delta_mt_2loop_at_at_S_flexiblesusy(
   TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL)
{
   ERROR_NO_TSIL("delta_mt_2loop_at_at_S_flexiblesusy");
}

TSIL_REAL delta_mt_2loop_at_at_LR_flexiblesusy(
   TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL)
{
   ERROR_NO_TSIL("delta_mt_2loop_at_at_LR_flexiblesusy");
}

/// SPheno convention ///

TSIL_REAL delta_mt_2loop_as_as_spheno(
   TSIL_REAL, TSIL_REAL, TSIL_REAL)
{
   ERROR_NO_TSIL("delta_mt_2loop_as_as_spheno");
}

TSIL_REAL delta_mt_2loop_as_at_S_spheno(
   TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL)
{
   ERROR_NO_TSIL("delta_mt_2loop_as_at_S_spheno");
}

TSIL_REAL delta_mt_2loop_as_at_LR_spheno(
   TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL)
{
   ERROR_NO_TSIL("delta_mt_2loop_as_at_LR_spheno");
}

TSIL_REAL delta_mt_2loop_at_at_S_spheno(
   TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL)
{
   ERROR_NO_TSIL("delta_mt_2loop_at_at_S_spheno");
}

TSIL_REAL delta_mt_2loop_at_at_LR_spheno(
   TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL, TSIL_REAL)
{
   ERROR_NO_TSIL("delta_mt_2loop_at_at_LR_spheno");
}

} // namespace sm_twoloop_mt
} // namespace flexiblesusy

#endif // ENABLE_TSIL
