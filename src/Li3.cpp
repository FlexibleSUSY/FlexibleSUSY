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

#include "Li3.hpp"
#include "horner.hpp"
#include "log.hpp"
#include <cfloat>
#include <cmath>
#include <complex>

namespace flexiblesusy {

namespace {

   /// Li_3(x) for x in [-1,0]
   double li3_neg(double x) noexcept
   {
      const double cp[] = {
         0.9999999999999999795e+0, -2.0281801754117129576e+0,
         1.4364029887561718540e+0, -4.2240680435713030268e-1,
         4.7296746450884096877e-2, -1.3453536579918419568e-3
      };
      const double cq[] = {
         1.0000000000000000000e+0, -2.1531801754117049035e+0,
         1.6685134736461140517e+0, -5.6684857464584544310e-1,
         8.1999463370623961084e-2, -4.0756048502924149389e-3,
         3.4316398489103212699e-5
      };

      const double x2 = x*x;
      const double x4 = x2*x2;
      const double p = cp[0] + x*cp[1] + x2*(cp[2] + x*cp[3]) +
         x4*(cp[4] + x*cp[5]);
      const double q = cq[0] + x*cq[1] + x2*(cq[2] + x*cq[3]) +
         x4*(cq[4] + x*cq[5] + x2*cq[6]);

      return x*p/q;
   }

    /// Li_3(x) for x in [0,1/2]
    double li3_pos(double x) noexcept
    {
       const double cp[] = {
          0.9999999999999999893e+0, -2.5224717303769789628e+0,
          2.3204919140887894133e+0, -9.3980973288965037869e-1,
          1.5728950200990509052e-1, -7.5485193983677071129e-3
       };
       const double cq[] = {
          1.0000000000000000000e+0, -2.6474717303769836244e+0,
          2.6143888433492184741e+0, -1.1841788297857667038e+0,
          2.4184938524793651120e-1, -1.8220900115898156346e-2,
          2.4927971540017376759e-4
       };

       const double x2 = x*x;
       const double x4 = x2*x2;
       const double p = cp[0] + x*cp[1] + x2*(cp[2] + x*cp[3]) +
          x4*(cp[4] + x*cp[5]);
       const double q = cq[0] + x*cq[1] + x2*(cq[2] + x*cq[3]) +
          x4*(cq[4] + x*cq[5] + x2*cq[6]);

       return x*p/q;
   }

} // anonymous namespace

/**
 * @brief Real trilogarithm \f$\operatorname{Li}_3(x)\f$
 * @param x real argument
 * @return \f$\operatorname{Li}_3(x)\f$
 * @author Alexander Voigt
 */
double Li3(double x) noexcept
{
   const double zeta2 = 1.6449340668482264;
   const double zeta3 = 1.2020569031595943;

   // transformation to [-1,0] and [0,1/2]
   if (x < -1) {
      const double l = std::log(-x);
      return li3_neg(1/x) - l*(zeta2 + 1.0/6*l*l);
   } else if (x == -1) {
      return -0.75*zeta3;
   } else if (x < 0) {
      return li3_neg(x);
   } else if (x == 0) {
      return x;
   } else if (x < 0.5) {
      return li3_pos(x);
   } else if (x == 0.5) {
      return 0.53721319360804020;
   } else if (x < 1) {
      const double l = std::log(x);
      return -li3_neg(1 - 1/x) - li3_pos(1 - x)
         + zeta3 + l*(zeta2 + l*(-0.5*std::log1p(-x) + 1.0/6*l));
   } else if (x == 1) {
      return zeta3;
   } else if (x < 2) {
      const double l = std::log(x);
      return -li3_neg(1 - x) - li3_pos(1 - 1/x)
         + zeta3 + l*(zeta2 + l*(-0.5*std::log(x - 1) + 1.0/6*l));
   } else { // x >= 2.0
      const double l = std::log(x);
      return li3_pos(1/x) + l*(2*zeta2 - 1.0/6*l*l);
   }
}

/**
 * @brief Complex trilogarithm \f$\operatorname{Li}_3(z)\f$
 * @param z complex argument
 * @return \f$\operatorname{Li}_3(z)\f$
 * @author Alexander Voigt
 */
std::complex<double> Li3(const std::complex<double>& z) noexcept
{
   const double PI    = 3.1415926535897932;
   const double zeta2 = 1.6449340668482264;
   const double zeta3 = 1.2020569031595943;
   const double bf[18] = {
      1.0                   , -3.0/8.0               ,
      17.0/216.0            , -5.0/576.0             ,
      1.2962962962962963e-04,  8.1018518518518519e-05,
     -3.4193571608537595e-06, -1.3286564625850340e-06,
      8.6608717561098513e-08,  2.5260875955320400e-08,
     -2.1446944683640648e-09, -5.1401106220129789e-10,
      5.2495821146008294e-11,  1.0887754406636318e-11,
     -1.2779396094493695e-12, -2.3698241773087452e-13,
      3.1043578879654623e-14,  5.2617586299125061e-15
   };

   const double rz = std::real(z);
   const double iz = std::imag(z);

   if (iz == 0) {
      if (rz <= 1) {
         return { Li3(rz), iz };
      } else {
         const double l = std::log(rz);
         return { Li3(rz), -0.5*PI*l*l };
      }
   }

   const double nz  = std::abs(z);
   const double pz  = std::arg(z);
   const double lnz = std::log(nz);

   if (lnz*lnz + pz*pz < 1) { // |log(z)| < 1
      const std::complex<double> u(lnz, pz); // log(z)
      const auto u2 = u*u;
      const auto u4 = u2*u2;
      const auto u8 = u4*u4;
      const auto c0 = zeta3 + u*(zeta2 - u2/12.0);
      const auto c1 = 0.25 * (3.0 - 2.0*pos_log(-u));

      const double cs[7] = {
         -3.4722222222222222e-03, 1.1574074074074074e-05,
         -9.8418997228521038e-08, 1.1482216343327454e-09,
         -1.5815724990809166e-11, 2.4195009792525152e-13,
         -3.9828977769894877e-15
      };

      return
         c0 +
         c1*u2 +
         u4*(cs[0] + u2*cs[1]) +
         u8*(cs[2] + u2*cs[3] + u4*(cs[4] + u2*cs[5])) +
         u8*u8*cs[6];
   }

   std::complex<double> u(0.0, 0.0), rest(0.0, 0.0);

   if (nz <= 1) {
      u = -log1p(-z);
   } else { // nz > 1
      const double arg = pz > 0.0 ? pz - PI : pz + PI;
      const std::complex<double> lmz(lnz, arg); // log(-z)
      u = -log1p(-1.0/z);
      rest = -lmz*(lmz*lmz/6.0 + zeta2);
   }

   const auto u2 = u*u;
   const auto u4 = u2*u2;
   const auto u8 = u4*u4;

   return
      rest +
      u*bf[0] +
      u2*(bf[1] + u*bf[2]) +
      u4*(bf[3] + u*bf[4] + u2*(bf[5] + u*bf[6])) +
      u8*(bf[7] + u*bf[8] + u2*(bf[9] + u*bf[10]) +
          u4*(bf[11] + u*bf[12] + u2*(bf[13] + u*bf[14]))) +
      u8*u8*(bf[15] + u*bf[16] + u2*bf[17]);
}

/**
 * @brief Complex trilogarithm \f$\operatorname{Li}_3(z)\f$ with long double precision
 * @param z complex argument
 * @return \f$\operatorname{Li}_3(z)\f$
 * @author Alexander Voigt
 */
std::complex<long double> Li3(const std::complex<long double>& z) noexcept
{
   const long double PI    = 3.14159265358979323846264338327950288L;
   const long double zeta2 = 1.64493406684822643647241516664602519L;
   const long double zeta3 = 1.20205690315959428539973816151144999L;
   const long double bf[] = {
      1.0L,
     -3.0L/8.0L,
      17.0L/216.0L,
     -5.0L/576.0L,
      7.0L/54000.0L,
      7.0L/86400.0L,
     -3.41935716085375949321527552820069827e-06L,
     -1.32865646258503401360544217687074830e-06L,
      8.66087175610985134794658604182413706e-08L,
      2.52608759553203997648442092886537331e-08L,
     -2.14469446836406476093388507573649032e-09L,
     -5.14011062201297891533581769272004962e-10L,
      5.24958211460082943639408880855807284e-11L,
      1.08877544066363183753729715704249107e-11L,
     -1.27793960944936953055818317540722120e-12L,
     -2.36982417730874520997977788101244891e-13L,
      3.10435788796546229428475327046556211e-14L,
      5.26175862991250608413183925112250061e-15L,
     -7.53847954994926536599250143226771028e-16L,
     -1.18623225777522852530825009512459322e-16L,
      1.83169799654913833820892731212815349e-17L,
      2.70681710318373501514907347126169436e-18L,
#if LDBL_DIG > 18
     -4.45543389782963882643263099217632212e-19L,
     -6.23754849225569465036532224739838641e-20L,
      1.08515215348745349131365609968642833e-20L,
      1.44911748660360819307349049665275324e-21L,
     -2.64663397544589903347408911861443741e-22L,
     -3.38976534885101047219258165860814078e-23L,
      6.46404773360331088903253098219534234e-24L,
      7.97583448960241242420922272590502795e-25L,
     -1.58091787902874833559211176293826770e-25L,
     -1.88614997296228681931102253988531956e-26L,
      3.87155366384184733039971271888313319e-27L,
      4.48011750023456073048653898320511684e-28L,
     -9.49303387191183612641753676027699150e-29L,
     -1.06828138090773812240182143033807908e-29L,
      2.33044789361030518600785199019281371e-30L,
      2.55607757265197540805635698286695865e-31L,
     -5.72742160613725968447274458033057100e-32L,
     -6.13471321379642358258549296897773326e-33L,
      1.40908086040689448401268688489421700e-33L,
      1.47642223976665341443182801167106626e-34L,
     -3.47010516489959160555004020312910903e-35L,
     -3.56210662409746357967357370318293608e-36L,
      8.55369656823692105754731289124468101e-37L
#endif
   };

   const long double rz = std::real(z);
   const long double iz = std::imag(z);

   if (iz == 0) {
      if (rz == 0) {
         return { rz, iz };
      }
      if (rz == 1) {
         return { zeta3, iz };
      }
      if (rz == -1) {
         return { -0.75L*zeta3, iz };
      }
      if (rz == 0.5L) {
         return { 0.537213193608040200940623225594965827L, iz };
      }
   }

   const long double nz  = std::abs(z);
   const long double pz  = std::arg(z);
   const long double lnz = std::log(nz);

   if (lnz*lnz + pz*pz < 1) { // |log(z)| < 1
      const std::complex<long double> u(lnz, pz); // log(z)
      const auto u2 = u*u;
      const auto c0 = zeta3 + u*(zeta2 - u2/12.0L);
      const auto c1 = 0.25L * (3.0L - 2.0L*pos_log(-u));

      const long double cs[] = {
        -3.47222222222222222222222222222222222e-03L,
         1.15740740740740740740740740740740741e-05L,
        -9.84189972285210380448475686570924666e-08L,
         1.14822163433274544385655496766607878e-09L,
        -1.58157249908091658933409775160616911e-11L,
         2.41950097925251519452732701564998016e-13L,
        -3.98289777698948774786517290926462002e-15L,
         6.92336661830592905806820954095065870e-17L,
        -1.25527223044997727545846570912655367e-18L,
#if LDBL_DIG > 18
         2.35375400276846523056441171414060379e-20L,
        -4.53639890345868701844750708901700830e-22L,
         8.94516967039264316712031170773304472e-24L,
        -1.79828400469549627172020247141015426e-25L,
         3.67549976479373844433604733912674099e-27L,
        -7.62080797156479522953948500963765478e-29L,
         1.60004196436948597517376392257325602e-30L,
        -3.39676114756037558792312060520851852e-32L,
         7.28227228675776469531725636144432664e-34L,
        -1.57502264795800348718497893940378261e-35L,
         3.43354009248058933588797212016527038e-37L
#endif
      };

      return c0 + u2*(c1 + u2*horner<0>(u2, cs));
   }

   std::complex<long double> u(0.0L, 0.0L), rest(0.0L, 0.0L);

   if (nz <= 1) {
      u = -log1p(-z);
   } else { // nz > 1
      const long double arg = pz > 0.0 ? pz - PI : pz + PI;
      const std::complex<long double> lmz(lnz, arg); // log(-z)
      u = -log1p(-1.0L/z);
      rest = -lmz*(lmz*lmz/6.0L + zeta2);
   }

   return rest + u*horner<0>(u, bf);
}

} // namespace flexiblesusy
