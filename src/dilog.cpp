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

#include "dilog.hpp"
#include "complex.hpp"
#include <cfloat>
#include <cmath>
#include <limits>

namespace flexiblesusy {

namespace {

   template <typename T, int N>
   T horner(T x, const T (&c)[N]) noexcept
   {
      T p = c[N - 1];
      for (int i = N - 2; i >= 0; --i) {
         p = p*x + c[i];
      }
      return p;
   }

   template <int Nstart, typename T, int N>
   Complex<T> horner(const Complex<T>& z, const T (&coeffs)[N]) noexcept
   {
      static_assert(0 <= Nstart && Nstart < N && N >= 2, "invalid array bounds");

      const T r = z.re + z.re;
      const T s = z.re * z.re + z.im * z.im;
      T a = coeffs[N - 1], b = coeffs[N - 2];

      for (int i = N - 3; i >= Nstart; --i) {
         const T t = a;
         a = b + r * a;
         b = coeffs[i] - s * t;
      }

      return Complex<T>(z.re*a + b, z.im*a);
   }

} // anonymous namespace

/**
 * @brief Clausen function \f$\mathrm{Cl}_2(\theta) = \mathrm{Im}(\mathrm{Li}_2(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\mathrm{Cl}_2(\theta)\f$
 * @author Alexander Voigt
 * @note Implemented as economized Padé approximation.
 */
double clausen_2(double x) noexcept
{
   const double PI = 3.14159265358979324;
   const double PI2 = 2*PI, PIH = PI/2, PI28 = PI*PI/8;
   double sgn = 1;

   if (x < 0) {
      x = -x;
      sgn = -1;
   }

   if (x >= PI2) {
      x = std::fmod(x, PI2);
   }

   if (x > PI) {
      const double p0 = 6.28125;
      const double p1 = 0.0019353071795864769253;
      x = (p0 - x) + p1;
      sgn = -sgn;
   }

   if (x == 0 || x == PI) {
      return 0;
   }

   double h = 0;

   if (x < PIH) {
      const double P[] = {
         2.7951565822419270e-02, -8.8865360514541522e-04,
         6.8282348222485902e-06, -7.5276232403566808e-09
      };
      const double Q[] = {
         1.0000000000000000e+00, -3.6904397961160525e-02,
         3.7342870576106476e-04, -8.7460760866531179e-07
      };
      const double y = x*x;
      const double z = y - PI28;
      const double z2 = z*z;
      const double p = P[0] + z * P[1] + z2 * (P[2] + z * P[3]);
      const double q = Q[0] + z * Q[1] + z2 * (Q[2] + z * Q[3]);

      h = x*(1 - std::log(x) + y*p/q/2);
   } else {
      const double P[] = {
         6.4005702446195512e-01, -2.0641655351338783e-01,
         2.4175305223497718e-02, -1.2355955287855728e-03,
         2.5649833551291124e-05, -1.4783829128773320e-07
      };
      const double Q[] = {
         1.0000000000000000e+00, -2.5299102015666356e-01,
         2.2148751048467057e-02, -7.8183920462457496e-04,
         9.5432542196310670e-06, -1.8184302880448247e-08
      };
      const double y = PI - x;
      const double z = y*y - PI28;
      const double z2 = z*z;
      const double z4 = z2*z2;
      const double p = P[0] + z * P[1] + z2 * (P[2] + z * P[3]) +
         z4 * (P[4] + z * P[5]);
      const double q = Q[0] + z * Q[1] + z2 * (Q[2] + z * Q[3]) +
         z4 * (Q[4] + z * Q[5]);

      h = y*p/q;
   }

   return sgn*h;
}

/**
 * @brief Clausen function \f$\mathrm{Cl}_2(\theta) = \mathrm{Im}(\mathrm{Li}_2(e^{i\theta}))\f$ with long double precision
 * @param x real angle
 * @return \f$\mathrm{Cl}_2(\theta)\f$
 * @author Alexander Voigt
 *
 * Implemented as an economized Padé approximation with a maximum
 * error of approximately 3.26e-41 (for long double and quadruple
 * precision).
 */
long double clausen_2(long double x) noexcept
{
   const long double PI = 3.14159265358979323846264338327950288L;
   const long double PI2 = 2*PI, PIH = PI/2, PI28 = PI*PI/8;
   long double sgn = 1;

   if (x < 0) {
      x = -x;
      sgn = -1;
   }

   if (x >= PI2) {
      x = std::fmod(x, PI2);
   }

   if (x > PI) {
      const long double p0 = 6.28125L;
      const long double p1 = 0.0019353071795864769252867665590057683943L;
      x = (p0 - x) + p1;
      sgn = -sgn;
   }

   if (x == 0 || x == PI) {
      return 0;
   }

   long double h = 0;

   if (x < PIH) {
      const long double P[] = {
         2.795156582241927046412081735910646612854e-02L,
        -2.704528039782130831727668760352473119745e-03L,
         1.058576547177802928762582430994046913011e-04L,
        -2.147507975446829791077479076828450780219e-06L,
         2.401415296681270093111305488326496124531e-08L,
        -1.450571790543608936928129678333156785370e-10L,
         4.280534901040925211965221454555516657749e-13L,
        -4.792802237226483806823208684186867186935e-16L,
         8.883657381852830471176782778999368430017e-20L
      };
      const long double Q[] = {
         1.000000000000000000000000000000000000000e+00L,
        -1.018694323414614410071369720193716012304e-01L,
         4.248408782245281612900840467370146443889e-03L,
        -9.337710301347963985908084056584187570954e-05L,
         1.159379163193822053103946363960603543601e-06L,
        -8.083352720393357000801675734774176899515e-09L,
         2.949313240431512997069808854213308209519e-11L,
        -4.742700419624204182400715964695278593777e-14L,
         2.158380636740175386190809152807629331877e-17L
      };
      const long double y = x*x;
      const long double z = y - PI28;
      const long double z2 = z*z;
      const long double z4 = z2*z2;
      const long double z8 = z4*z4;
      const long double p = P[0] + z * P[1] + z2 * (P[2] + z * P[3]) +
         z4 * (P[4] + z * P[5] + z2 * (P[6] + z * P[7])) + z8 * P[8];
      const long double q = Q[0] + z * Q[1] + z2 * (Q[2] + z * Q[3]) +
         z4 * (Q[4] + z * Q[5] + z2 * (Q[6] + z * Q[7])) + z8 * Q[8];

      h = x*(1 - std::log(x) + y*p/q/2);
   } else {
      const long double P[] = {
         6.400570244619551220929428522356830562481e-01L,
        -4.651631624886004423703445967760673575997e-01L,
         1.487130845262105644024901814213146749895e-01L,
        -2.749665174801454303884783494225610407035e-02L,
         3.251522465413666561950482170352156048203e-03L,
        -2.567438381297475310848635518657180974512e-04L,
         1.372076105130164861564020074178493529151e-05L,
        -4.924179093673498700461153483531075799113e-07L,
         1.153267936031337440182387313169828395860e-08L,
        -1.667578310677508029208023423625588832295e-10L,
         1.348437292247918547169070120217729056878e-12L,
        -5.052245092698477071447850656280954693011e-15L,
         5.600638109466570497480415519182233229048e-18L,
      };
      const long double Q[] = {
         1.000000000000000000000000000000000000000e+00L,
        -6.572465772185054284667746526549393897676e-01L,
         1.886234634096976582977630140163583172173e-01L,
        -3.103347567899737687117030083178445406132e-02L,
         3.230860399291224478336071920154030050234e-03L,
        -2.216546569335921728108984951507368387512e-04L,
         1.011949972138985643994631167412420906088e-05L,
        -3.033400935206767852937290458763547850384e-07L,
         5.748454611964843057644023468691231929690e-09L,
        -6.408350048413952604351408631173781906861e-11L,
         3.678584366662951864267349037579031493395e-13L,
        -8.240439699357036167611014086997683837396e-16L,
         3.041049046123062158788159773779755292771e-19L
      };
      const long double y = PI - x;
      const long double z = y*y - PI28;
      const long double z2 = z*z;
      const long double z4 = z2*z2;
      const long double z8 = z4*z4;
      const long double p = P[0] + z * P[1] + z2 * (P[2] + z * P[3]) +
         z4 * (P[4] + z * P[5] + z2 * (P[6] + z * P[7])) +
         z8 * (P[8] + z * P[9] + z2 * (P[10] + z * P[11]) + z4 * P[12]);
      const long double q = Q[0] + z * Q[1] + z2 * (Q[2] + z * Q[3]) +
         z4 * (Q[4] + z * Q[5] + z2 * (Q[6] + z * Q[7])) +
         z8 * (Q[8] + z * Q[9] + z2 * (Q[10] + z * Q[11]) + z4 * Q[12]);

      h = y*p/q;
   }

   return sgn*h;
}

/**
 * @brief Real dilogarithm \f$\mathrm{Li}_2(x)\f$
 * @param x real argument
 * @return \f$\mathrm{Li}_2(x)\f$
 * @author Alexander Voigt
 *
 * Implemented as an economized Pade approximation with a
 * maximum error of 4.16e-18.
 */
double dilog(double x) noexcept
{
   const double PI = 3.1415926535897932;
   const double P[] = {
      1.0706105563309304277e+0,
     -4.5353562730201404017e+0,
      7.4819657596286408905e+0,
     -6.0516124315132409155e+0,
      2.4733515209909815443e+0,
     -4.6937565143754629578e-1,
      3.1608910440687221695e-2,
     -2.4630612614645039828e-4
   };
   const double Q[] = {
      1.0000000000000000000e+0,
     -4.5355682121856044935e+0,
      8.1790029773247428573e+0,
     -7.4634190853767468810e+0,
      3.6245392503925290187e+0,
     -8.9936784740041174897e-1,
      9.8554565816757007266e-2,
     -3.2116618742475189569e-3
   };

   double y = 0, r = 0, s = 1;

   // transform to [0, 1/2)
   if (x < -1) {
      const double l = std::log(1 - x);
      y = 1/(1 - x);
      r = -PI*PI/6 + l*(0.5*l - std::log(-x));
      s = 1;
   } else if (x == -1) {
      return -PI*PI/12;
   } else if (x < 0) {
      const double l = std::log1p(-x);
      y = x/(x - 1);
      r = -0.5*l*l;
      s = -1;
   } else if (x == 0) {
      return 0;
   } else if (x < 0.5) {
      y = x;
      r = 0;
      s = 1;
   } else if (x < 1) {
      y = 1 - x;
      r = PI*PI/6 - std::log(x)*std::log(1 - x);
      s = -1;
   } else if (x == 1) {
      return PI*PI/6;
   } else if (x < 2) {
      const double l = std::log(x);
      y = 1 - 1/x;
      r = PI*PI/6 - l*(std::log(1 - 1/x) + 0.5*l);
      s = 1;
   } else {
      const double l = std::log(x);
      y = 1/x;
      r = PI*PI/3 - 0.5*l*l;
      s = -1;
   }

   const double z = y - 0.25;
   const double z2 = z*z;
   const double z4 = z2*z2;

   const double p = P[0] + z * P[1] + z2 * (P[2] + z * P[3]) +
                    z4 * (P[4] + z * P[5] + z2 * (P[6] + z * P[7]));
   const double q = Q[0] + z * Q[1] + z2 * (Q[2] + z * Q[3]) +
                    z4 * (Q[4] + z * Q[5] + z2 * (Q[6] + z * Q[7]));

   return r + s*y*p/q;
}

/**
 * @brief Real dilogarithm \f$\mathrm{Li}_2(z)\f$ with long double precision
 * @param x real argument
 * @return \f$\mathrm{Li}_2(z)\f$
 * @author Alexander Voigt
 *
 * Implemented as an economized Pade approximation with a maximum
 * error of 2.13e-20 (long double precision) and 1.03e-38 (quadruple
 * precision), respectively.
 */
long double dilog(long double x) noexcept
{
   const long double PI  = 3.14159265358979323846264338327950288L;

#if LDBL_DIG <= 18
   const long double P[] = {
      1.07061055633093042767673531395124630e+0L,
     -5.25056559620492749887983310693176896e+0L,
      1.03934845791141763662532570563508185e+1L,
     -1.06275187429164237285280053453630651e+1L,
      5.95754800847361224707276004888482457e+0L,
     -1.78704147549824083632603474038547305e+0L,
      2.56952343145676978700222949739349644e-1L,
     -1.33237248124034497789318026957526440e-2L,
      7.91217309833196694976662068263629735e-5L
   };
   const long double Q[] = {
      1.00000000000000000000000000000000000e+0L,
     -5.20360694854541370154051736496901638e+0L,
      1.10984640257222420881180161591516337e+1L,
     -1.24997590867514516374467903875677930e+1L,
      7.97919868560471967115958363930214958e+0L,
     -2.87732383715218390800075864637472768e+0L,
      5.49210416881086355164851972523370137e-1L,
     -4.73366369162599860878254400521224717e-2L,
      1.23136575793833628711851523557950417e-3L
   };
#else
   const long double P[] = {
      1.0706105563309304276767353139512463033e+0L,
     -1.0976999936495889694316759019749736514e+1L,
      5.0926847456521671435598196295391041399e+1L,
     -1.4137651316583179225403308056234710093e+2L,
      2.6167438966024905838946321639772104746e+2L,
     -3.4058684964478756354428571563597134889e+2L,
      3.2038320769322335817541819891481680235e+2L,
     -2.2042602382067574460998180123431383994e+2L,
      1.1098617775200190748144566197068215051e+2L,
     -4.0511655788875306581280201485439376229e+1L,
      1.0505482055068989911823839054507173824e+1L,
     -1.8711701614474515075977773430379529915e+0L,
      2.1700009060344649725726270289240541652e-1L,
     -1.5030137964245010798524220531488795883e-2L,
      5.3441650657913964794668981233490297624e-4L,
     -7.0813883289232115572593407126124042503e-6L,
      9.9037749983459843207756017329825486520e-9L
   };
   const long double Q[] = {
      1.0000000000000000000000000000000000000e+0L,
     -1.0552362671556327276432317611336360654e+1L,
      5.0559576537049301822461383847757604795e+1L,
     -1.4554341156550484644439941130107249144e+2L,
      2.8070530313005385434963768448877956414e+2L,
     -3.8295927403204227055447294669165437162e+2L,
      3.8034123327297619837621395664129395552e+2L,
     -2.7878320883603471098333035075019703189e+2L,
      1.5127204204944666467295435798726892255e+2L,
     -6.0401294167088266781532200159750431265e+1L,
      1.7480717870136655299932069333598290551e+1L,
     -3.5731199220918363429317367913920234597e+0L,
      4.9537900060585746351429453828993042241e-1L,
     -4.3750415383481013167382471956794084257e-2L,
      2.2234568413141078158637017098109483650e-3L,
     -5.4149527323164210917629185927908003034e-5L,
      4.1629116823949062782848730828350635143e-7L
   };
#endif

   long double y = 0, r = 0, s = 1;

   // transform to [0, 1/2)
   if (x < -1) {
      const long double l = std::log(1 - x);
      y = 1/(1 - x);
      r = -PI*PI/6 + l*(0.5L*l - std::log(-x));
      s = 1;
   } else if (x == -1) {
      return -PI*PI/12;
   } else if (x < 0) {
      const long double l = std::log1p(-x);
      y = x/(x - 1);
      r = -0.5L*l*l;
      s = -1;
   } else if (x == 0) {
      return 0;
   } else if (x < 0.5L) {
      y = x;
      r = 0;
      s = 1;
   } else if (x < 1) {
      y = 1 - x;
      r = PI*PI/6 - std::log(x)*std::log(1 - x);
      s = -1;
   } else if (x == 1) {
      return PI*PI/6;
   } else if (x < 2) {
      const long double l = std::log(x);
      y = 1 - 1/x;
      r = PI*PI/6 - l*(std::log(1 - 1/x) + 0.5L*l);
      s = 1;
   } else {
      const long double l = std::log(x);
      y = 1/x;
      r = PI*PI/3 - 0.5L*l*l;
      s = -1;
   }

   const long double z = y - 0.25L;

   const long double p = horner(z, P);
   const long double q = horner(z, Q);

   return r + s*y*p/q;
}

/**
 * @brief Complex dilogarithm \f$\mathrm{Li}_2(z)\f$
 * @param z_ complex argument
 * @return \f$\mathrm{Li}_2(z)\f$
 * @note Implementation translated from SPheno to C++
 * @author Werner Porod
 * @note translated to C++ by Alexander Voigt
 */
std::complex<double> dilog(const std::complex<double>& z_) noexcept
{
   const double PI = 3.1415926535897932;
   const Complex<double> z = { std::real(z_), std::imag(z_) };

   // bf[1..N-1] are the even Bernoulli numbers / (2 n + 1)!
   // generated by: Table[BernoulliB[2 n]/(2 n + 1)!, {n, 1, 9}]
   const double bf[10] = {
      - 1.0/4.0,
      + 1.0/36.0,
      - 1.0/3600.0,
      + 1.0/211680.0,
      - 1.0/10886400.0,
      + 1.0/526901760.0,
      - 4.0647616451442255e-11,
      + 8.9216910204564526e-13,
      - 1.9939295860721076e-14,
      + 4.5189800296199182e-16
   };

   // special cases
   if (z.im == 0) {
      if (z.re <= 1) {
         return dilog(z.re);
      }
      // z.re > 1
      return { dilog(z.re), -PI*std::log(z.re) };
   }

   const double nz = norm_sqr(z);

   if (nz < std::numeric_limits<double>::epsilon()) {
      return z*(1.0 + 0.25*z);
   }

   Complex<double> u(0.0, 0.0), rest(0.0, 0.0);
   double sgn = 1;

   // transformation to |z|<1, Re(z)<=0.5
   if (z.re <= 0.5) {
      if (nz > 1) {
         const Complex<double> lz = log(-z);
         u = -log(1.0 - 1.0 / z);
         rest = -0.5*lz*lz - PI*PI/6;
         sgn = -1;
      } else { // nz <= 1
         u = -log(1.0 - z);
         rest = 0;
         sgn = 1;
      }
   } else { // z.re > 0.5
      if (nz <= 2*z.re) {
         u = -log(z);
         rest = u*log(1.0 - z) + PI*PI/6;
         sgn = -1;
      } else { // nz > 2*z.re
         const Complex<double> lz = log(-z);
         u = -log(1.0 - 1.0 / z);
         rest = -0.5*lz*lz - PI*PI/6;
         sgn = -1;
      }
   }

   const Complex<double> u2(u*u);

   return sgn*(u + u2*(bf[0] + u*horner<1>(u2, bf))) + rest;
}

/**
 * @brief Complex dilogarithm \f$\mathrm{Li}_2(z)\f$ with long double precision
 * @param z_ complex argument
 * @return \f$\mathrm{Li}_2(z)\f$
 * @note Implementation translated from SPheno to C++
 * @author Werner Porod
 * @note translated to C++ and extended to long double precision by Alexander Voigt
 */
std::complex<long double> dilog(const std::complex<long double>& z_) noexcept
{
   const long double PI = 3.14159265358979323846264338327950288L;
   const Complex<long double> z = { std::real(z_), std::imag(z_) };

   // bf[1..N-1] are the even Bernoulli numbers / (2 n + 1)!
   // generated by: Table[BernoulliB[2 n]/(2 n + 1)!, {n, 1, 22}]
   const long double bf[] = {
      -1.0L/4.0L                                 ,
       1.0L/36.0L                                ,
      -1.0L/3600.0L                              ,
       1.0L/211680.0L                            ,
      -1.0L/10886400.0L                          ,
       1.0L/526901760.0L                         ,
      -4.06476164514422552680590938629196667e-11L,
       8.92169102045645255521798731675274885e-13L,
      -1.99392958607210756872364434779378971e-14L,
       4.51898002961991819165047655285559323e-16L,
      -1.03565176121812470144834115422186567e-17L,
#if LDBL_DIG > 18
       2.39521862102618674574028374300098038e-19L,
      -5.58178587432500933628307450562541991e-21L,
       1.30915075541832128581230739918659230e-22L,
      -3.08741980242674029324227976486646243e-24L,
       7.31597565270220342035790560925214859e-26L,
      -1.74084565723400074098905514775970255e-27L,
       4.15763564461389971961789962077522667e-29L,
      -9.96214848828462210319400670245583885e-31L,
       2.39403442489616530052116798789374956e-32L,
      -5.76834735536739008429179316187765424e-34L,
       1.39317947964700797782788660391154833e-35L,
      -3.37212196548508947046847363525493096e-37L
#endif
   };

   // special cases
   if (z.im == 0) {
      if (z.re <= 1) {
         return dilog(z.re);
      }
      // z.re > 1
      return { dilog(z.re), -PI*std::log(z.re) };
   }

   const long double nz = norm_sqr(z);

   if (nz < std::numeric_limits<long double>::epsilon()) {
      return z*(1.0L + 0.25L*z);
   }

   Complex<long double> u(0.0L, 0.0L), rest(0.0L, 0.0L);
   long double sgn = 1;

   // transformation to |z|<1, Re(z)<=0.5
   if (z.re <= 0.5L) {
      if (nz > 1) {
         const Complex<long double> lz = log(-z);
         u = -log(1.0L - 1.0L/z);
         rest = -0.5L*lz*lz - PI*PI/6;
         sgn = -1;
      } else { // nz <= 1
         u = -log(1.0L - z);
         rest = 0;
         sgn = 1;
      }
   } else { // z.re > 0.5L
      if (nz <= 2*z.re) {
         u = -log(z);
         rest = u*log(1.0L - z) + PI*PI/6;
         sgn = -1;
      } else { // nz > 2*z.re
         const Complex<long double> lz = log(-z);
         u = -log(1.0L - 1.0L/z);
         rest = -0.5L*lz*lz - PI*PI/6;
         sgn = -1;
      }
   }

   const Complex<long double> u2(u*u);

   return sgn*(u + u2*(bf[0] + u*horner<1>(u2, bf))) + rest;
}

} // namespace flexiblesusy
