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

#include "mssm_twoloophiggs.hpp"
#include "mssm_twoloophiggs.h"
#include "config.h"
#include "dilog.hpp"
#include <cmath>
#include <limits>
#include <utility>

#ifdef ENABLE_THREADS
   #include <mutex>
   #define LOCK_MUTEX() std::lock_guard<std::mutex> lg(mtx_mssm)
   namespace flexiblesusy {namespace mssm_twoloophiggs {
      static std::mutex mtx_mssm; /// locks MSSM fortran functions
   } // namespace mssm_twoloophiggs
} // namespace flexiblesusy
#else
   #define LOCK_MUTEX()
#endif

namespace flexiblesusy {
namespace mssm_twoloophiggs {

namespace {

template <typename T> T constexpr sqr(T a) { return a * a; }
template <typename T> T constexpr pow3(T a) { return a * a * a; }
template <typename T> T sqrtabs(T a) { return std::sqrt(std::abs(a)); }
template <typename T> T logabs(T x) { return std::log(std::abs(x)); }

double phi(double x, double y, double z)
{
   using std::log;

   const double u = x/z, v = y/z;
   const double lambda = sqrtabs(sqr(1 - u - v) - 4*u*v);
   const double xp = 0.5 * (1 + (u - v) - lambda);
   const double xm = 0.5 * (1 - (u - v) - lambda);

   return 1./lambda * (2*logabs(xp)*logabs(xm) - logabs(u)*logabs(v) -
                       2*(dilog(xp) + dilog(xm)) + M_PI*M_PI/3.);
}

/// First derivative of phi[t,T,g] w.r.t. T
double dphi_010(double t, double T, double g)
{
   using std::fabs;
   using std::sqrt;
   using std::log;
   using std::pow;

   constexpr double Pi2 = M_PI * M_PI;
   const double g2 = sqr(g);
   const double abbr = (-4*t*T)/g2 + sqr(1 - t/g - T/g);
   const double rabbr = sqrtabs(abbr);

   return ((g + t - T)*(Pi2 - 6*dilog((g - rabbr*g + t - T)/(2.*g)) -
      6*dilog((g - rabbr*g - t + T)/(2.*g)) -
      3*logabs(t/g)*logabs(T/g) + 6*logabs((g - rabbr*g + t -
      T)/(2.*g))*logabs((g - rabbr*g - t + T)/(2.*g))) + (3*rabbr*g* (
      rabbr*g*((-1 + rabbr)*g + t - T)*logabs(t/g) +
      2*T*(-2*g*logabs(4.) + (g + rabbr*g + t - T)*logabs((g - rabbr*g
      + t - T)/g) + (g + rabbr*g + t - T)*logabs((g + rabbr*g + t -
      T)/g) + g*logabs((g - rabbr*g - t + T)/g) - rabbr*g*logabs((g -
      rabbr*g - t + T)/g) - t*logabs((g - rabbr*g - t + T)/g) +
      T*logabs((g - rabbr*g - t + T)/g) + g*logabs((g + rabbr*g - t +
      T)/g) - rabbr*g*logabs((g + rabbr*g - t + T)/g) - t*logabs((g +
      rabbr*g - t + T)/g) + T*logabs((g + rabbr*g - t + T)/g)) ) ) /
      (T*(g - rabbr*g - t + T)))/(3.*pow(fabs(abbr),1.5)*g2);
}

/// First derivative of phi[g,t,T] w.r.t. T
double dphi_001(double g, double t, double T)
{
   using std::sqrt;

   const double Pi2 = 9.869604401089359;
   const double T2 = sqr(T);
   const double T3 = T2*T;
   const double x = sqrt(sqr(1 - g/T - t/T) - 4*g*t/T2);
   const double y = -(2*(g/T2 + t/T2)*(1 - g/T - t/T) + (8*g*t)/T3)/(2*x);
   const double ym = -g/T + t/T;
   const double yp = g/T - t/T;
   const double lgT = logabs(g/T);
   const double ltT = logabs(t/T);
   const double lxmym = logabs(0.5*(1 - x + ym));
   const double lxmyp = logabs(0.5*(1 - x + yp));
   const double lxpym = logabs(0.5*(1 + x + ym));
   const double lxpyp = logabs(0.5*(1 + x + yp));
   const double li2xym = dilog(0.5*(1 - x + ym));
   const double li2xyp = dilog(0.5*(1 - x + yp));

   return ((t*(t - T) - g*(2*t + T) + sqr(g))*
       (-6*li2xym - 6*li2xyp - 3*lgT*ltT + 6*lxmym*lxmyp + Pi2)
       + 3*T*(lgT*T + ltT*T +
              +2*lxmyp*(g - t + y*T2)/(1 - x + ym)
              +2*lxmym*(-g + t + y*T2)/(1 - x + yp)
              -2*lxpyp*(g - t + y*T2)/(-1 + x - ym)
              -2*lxpym*(-g + t + y*T2)/(-1 + x - yp)
              )*sqr(x))/
      (3.*pow3(T)*pow3(x));
}

double calc_At(double mt2, double mst12, double mst22,
   double sxt, double cxt, double mu, double tanb)
{
   const double s2t = 2*cxt*sxt;
   const double Xt = (mst12 - mst22)*s2t/2./sqrtabs(mt2);
   const double At = Xt - mu/tanb;

   return At;
}

/// limit st -> 0
Eigen::Matrix<double, 2, 1> tadpole_higgs_2loop_at_as_mssm_st_0(
   double mt2, double mg, double mst12, double mst22,
   double /* sxt */, double /* cxt */, double scale2,
   double mu, double tanb, double vev2, double gs)
{
   using std::atan;
   using std::sin;
   using std::cos;

   const double gs2 = sqr(gs);
   const double g = sqr(mg);
   const double q = scale2;
   const double q2 = sqr(q);
   const double t = mt2;
   const double T1 = mst12;
   const double T2 = mst22;
   const double g2 = sqr(g);
   const double v = sqrtabs(vev2);
   const double beta = std::atan(tanb);
   const double v2 = v * std::sin(beta);
   const double v1 = v * std::cos(beta);
   const double ltg = logabs(t/g);
   const double ltq = logabs(t/q);
   const double lgq = logabs(g/q);
   const double lT1q = logabs(T1/q);
   const double lT2q = logabs(T2/q);
   const double lgtq2 = logabs(g*t/q2);
   const double del1 = g2 + sqr(t) + sqr(T1) - 2*(g*t + g*T1 + t*T1);
   const double del2 = g2 + sqr(t) + sqr(T2) - 2*(g*t + g*T2 + t*T2);

   const double t1 =
      (16*mg*mu*(T1*T2*(g*(-lT1q + lT2q)*ltg + lT1q*ltg*t - lT2q*ltg*t + 5*T1 +
                        (-4 + lgtq2)*lT1q*T1 - lgq*ltq*T1 - 5*T2 + 4*lT2q*T2 -
                        lgtq2*lT2q*T2 + lgq*ltq*T2) + del1*T2*phi(g,t,T1) -
                 del2*T1*phi(g,t,T2)))/(T1*(T1 - T2)*T2*tanb*sqr(v1));

   const double t2 =
      (16*(-(T2*(del1*mg*mu - (g2 + g*(t - T1))*(T1 - T2)*tanb)*phi(g,t,T1)) +
           T1*(mg*mu*T2*(g*(lT1q - lT2q)*ltg + lT2q*ltg*t - 5*T1 + lgq*ltq*T1 -
                         lT1q*(ltg*t + (-4 + lgtq2)*T1) +
                         (5 + (-4 + lgtq2)*lT2q - lgq*ltq)*T2) +
               T2*(-T1 + T2)*(g*(2 - 2*lgq*ltq + lT1q*(-2 + lgq + ltq) +
                                 lT2q*(-2 + lgq + ltq)) +
                              2*(-5 + lT1q*(-1 + ltq) + lT2q*(-1 + ltq) - 3*(-2 + ltq)*ltq)*
                              t + (-4 + lT1q)*lT1q*T1 + (-4 + lT2q)*lT2q*T2 + 5*(T1 + T2))*
               tanb + (del2*mg*mu + (g2 + g*(t - T2))*(T1 - T2)*tanb)*phi(g,t,T2)
              )))/(T1*(T1 - T2)*T2*tanb*sqr(v2));

   Eigen::Matrix<double, 2, 1> result;
   result << t1, t2;

   const double k2 = 0.00004010149318236068; // 1/(4 Pi)^4
   const double pref = k2*mt2*gs2;

   return -result*pref;
}

/// limit st -> 0 and mst1 -> mst2
Eigen::Matrix<double, 2, 1> tadpole_higgs_2loop_at_as_mssm_st_0_mst1_eq_mst2(
   double mt2, double mg, double mst12, double /* mst22 */,
   double /* sxt */, double /* cxt */, double scale2,
   double mu, double tanb, double vev2, double gs)
{
   using std::atan;
   using std::sin;
   using std::cos;

   const double gs2 = sqr(gs);
   const double q = scale2;
   const double q2 = sqr(q);
   const double g = sqr(mg);
   const double g2 = sqr(g);
   const double t = mt2;
   const double T = mst12;
   const double Tsqr = sqr(T);
   const double ltg = logabs(t/g);
   const double lTq = logabs(T/q);
   const double ltq = logabs(t/q);
   const double lgq = logabs(g/q);
   const double lgtq2 = logabs(g*t/q2);
   const double v = sqrtabs(vev2);
   const double beta = std::atan(tanb);
   const double v2 = v * std::sin(beta);
   const double v1 = v * std::cos(beta);

   const double t1 =
      (16*mg*mu*(-((g2 - 2*g*t + sqr(t) - Tsqr)*phi(g,t,T)) +
                 T*(ltg*(-g + t) + (1 + lgtq2 - lgq*ltq - 4*lTq + lgtq2*lTq)*T +
                    (g2 - 2*t*T - 2*g*(t + T) + sqr(t) + Tsqr)*
                    dphi_001(g,t,T))))/(tanb*Tsqr*sqr(v1));

   const double t2 =
      (-16*(-((g2*(mg*mu + 2*T*tanb) + mg*mu*(sqr(t) - Tsqr) -
               2*g*(mg*mu*t - t*T*tanb + tanb*Tsqr))*phi(g,t,T)) +
            T*(ltg*mg*mu*(-g + t) +
               T*((1 + lgtq2 - lgq*ltq - 4*lTq + lgtq2*lTq)*mg*mu +
                  2*tanb*(g*(1 + (-2 + ltq)*lTq + lgq*(-ltq + lTq)) +
                          t*(-5 - 2*lTq + 2*ltq*(3 + lTq) - 3*sqr(ltq)) +
                          T*(5 - 4*lTq + sqr(lTq)))) +
               mg*mu*(g2 - 2*t*T - 2*g*(t + T) + sqr(t) + Tsqr)*
               dphi_001(g,t,T))))/(tanb*Tsqr*sqr(v2));

   Eigen::Matrix<double, 2, 1> result;
   result << t1, t2;

   const double k2 = 0.00004010149318236068; // 1/(4 Pi)^4
   const double pref = k2*mt2*gs2;

   return -result*pref;
}

/// Pietro Slavich implementation
Eigen::Matrix<double, 2, 1> tadpole_higgs_2loop_at_as_mssm_general(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2,
   double mu, double tanb, double vev2, double gs)
{
   Eigen::Matrix<double, 2, 1> result;

   ewsb2loop_(&mt2, &mg, &mst12, &mst22, &sxt, &cxt, &scale2,
              &mu, &tanb, &vev2, &gs, &result(0), &result(1));

   // workaround for intel or Eigen bug causing unexpected behaviour
   // of result.allFinite()
   if (!std::isfinite(result(0)) || !std::isfinite(result(1)))
       result.setZero();

   return -result;
}

/// limit st -> 0 and mst1 -> mst2
Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_at_as_mssm_with_tadpoles_st_0_mst1_eq_mst2(
   double mt2, double mg, double mst12, double /* mst22 */,
   double /* sxt */, double /* cxt */, double scale2, double mu,
   double tanb, double vev2, double gs, int /* scheme */)
{
   using std::atan;
   using std::sin;

   const double gs2 = sqr(gs);
   const double g = sqr(mg);
   const double g2 = sqr(g);
   const double q = scale2;
   const double t = mt2;
   const double T = mst12;
   const double t2 = sqr(t);
   const double t3 = t2*t;
   const double Tsqr = sqr(T);
   const double Tcub = Tsqr*T;
   const double ltg = logabs(t/g);
   const double lTg = logabs(T/g);
   const double lTq = logabs(T/q);
   const double ltq = logabs(t/q);
   const double lgq = logabs(g/q);
   const double lT2t2 = logabs(Tsqr/t2);
   const double del = g2 + t2 + Tsqr - 2*(g*t + g*T + t*T);
   const double sb = sin(atan(tanb));
   const double ht2 = 2./vev2*mt2/sqr(sb);

   Eigen::Matrix<double, 2, 2> result;

   result(0,0) = 0.;
   result(0,1) = (32*mg*mu*(-1 + ltq - T*dphi_010(t,T,g)))/T;
   result(1,0) = result(0,1);
   result(1,1) =
      (8*(8*del*mg*mu - 8*del*ltq*mg*mu + 8*del*g*tanb - 8*del*g*lgq*tanb +
          8*del*t*tanb - 8*del*lgq*t*tanb - 8*del*T*tanb +
          5*del*lT2t2*T*tanb + 8*del*ltg*T*tanb + 4*g2*ltg*T*tanb -
          4*del*lTg*T*tanb + 16*g2*lTg*T*tanb + 6*del*ltq*T*tanb +
          2*del*lTq*T*tanb + 40*g*ltg*t*T*tanb + 8*g*ltg*t2*tanb +
          12*ltg*T*t2*tanb - 8*ltg*t3*tanb - 4*ltg*tanb*Tcub +
          8*g*(g + t - T)*T*tanb*phi(t,T,g) + del*T*tanb*sqr(lT2t2) +
          8*del*T*tanb*sqr(ltq) - 8*del*T*tanb*sqr(lTq) +
          8*del*mg*mu*T*dphi_010(t,T,g)))/(del*T*tanb);

   const double k2 = 0.00004010149318236068; // 1/(4 Pi)^4
   const double pref = k2*ht2*mt2*gs2;

   return -result*pref;
}

/// Pietro Slavich implementation
Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_at_as_mssm_with_tadpoles_general(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double mu,
   double tanb, double vev2, double gs, int scheme)
{
   Eigen::Matrix<double, 2, 2> result;

   dszhiggs_(&mt2, &mg, &mst12, &mst22, &sxt, &cxt, &scale2, &mu,
             &tanb, &vev2, &gs, &scheme,
             &result(0,0), &result(1,1), &result(0,1));

   result(1,0) = result(0,1);

   return -result;
}

double self_energy_pseudoscalar_2loop_at_as_mssm_with_tadpoles_mst1_eq_mst2(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double mu,
   double tanb, double vev2, double gs)
{
   using std::atan;
   using std::log;
   using std::sin;

   constexpr double Pi2 = M_PI * M_PI;
   const double g = sqr(mg);
   const double g2 = sqr(g);
   const double q = scale2;
   const double q2 = sqr(scale2);
   const double t = mt2;
   const double T = mst12;
   const double sb = sin(atan(tanb));
   const double ht2 = 2./vev2*mt2/sqr(sb);
   const double At = calc_At(mt2, mst12, mst22, sxt, cxt, mu, tanb);

   const double result = (-2*(g*(2*At*g + 2*At*t - At*T + mg*T + mg*(g
      - t)*logabs(g/t) - At*T*logabs(g/q)*logabs(t/q) -
      mg*T*logabs(g/q)*logabs(t/q) - 4*mg*T*logabs(T/q) -
      2*At*T*sqr(logabs(T/q)) + logabs((g*t)/q2)*(-(At*(g + t - T)) +
      mg*T + (At + mg)*T*logabs(T/q))) - 2*(At + mg)*(g + t -
      T)*T*phi(t,T,g) + T*(At*(g2 + sqr(t - T) - 2*g*T) + mg*(g2 +
      sqr(t - T) - 2*g*(t + T)))*dphi_010(t,T,g)))/ (g*T);

   const double pref = 4*sqr(gs)/sqr(16*Pi2) * ht2*mu*(1./tanb + tanb);

   return -pref * result;
}

double self_energy_pseudoscalar_2loop_at_as_mssm_with_tadpoles_general(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double mu,
   double tanb, double vev2, double gs)
{
   double result;

   dszodd_(&mt2, &mg, &mst12, &mst22, &sxt, &cxt, &scale2, &mu,
           &tanb, &vev2, &gs, &result);

   return -result;
}

} // anonymous namespace

Eigen::Matrix<double, 2, 1> tadpole_higgs_2loop_at_as_mssm(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2,
   double mu, double tanb, double vev2, double gs)
{
   if (std::abs(sxt) < 1e-8) {
      if (std::abs((mst12 - mst22)/mst12) < 1e-6)
         return tadpole_higgs_2loop_at_as_mssm_st_0_mst1_eq_mst2(
            mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, gs);

      return tadpole_higgs_2loop_at_as_mssm_st_0(
         mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, gs);
   }

   return tadpole_higgs_2loop_at_as_mssm_general(
      mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, gs);
}

Eigen::Matrix<double, 2, 1> tadpole_higgs_2loop_at_at_mssm(
   double mt2, double mb2, double mA2, double mst12,
   double mst22, double msb12, double msb22,
   double sxt, double cxt, double sxb, double cxb,
   double scale2, double mu, double tanb, double vev2)
{
   Eigen::Matrix<double, 2, 1> result;

   {
      LOCK_MUTEX();

      ddstad_(&mt2, &mb2, &mA2, &mst12, &mst22, &msb12, &msb22,
              &sxt, &cxt, &sxb, &cxb, &scale2, &mu, &tanb, &vev2,
              &result(0), &result(1));
   }

   // workaround for intel or Eigen bug causing unexpected behaviour
   // of result.allFinite()
   if (!std::isfinite(result(0)) || !std::isfinite(result(1)))
       result.setZero();

   return -result;
}

Eigen::Matrix<double, 2, 1> tadpole_higgs_2loop_ab_as_mssm(
   double mb2, double mg, double msb12, double msb22,
   double sxb, double cxb, double scale2,
   double mu, double cotb, double vev2, double gs)
{
   Eigen::Matrix<double, 2, 1> result(tadpole_higgs_2loop_at_as_mssm(
      mb2, mg, msb12, msb22, sxb, cxb, scale2,
      mu, cotb, vev2, gs));

   std::swap(result(0), result(1));

   return result;
}

Eigen::Matrix<double, 2, 1> tadpole_higgs_2loop_atau_atau_mssm(
   double mtau2, double mA2, double msv2, double mstau12,
   double mstau22, double sintau, double costau, double scale2,
   double mu, double tanb, double vev2)
{
   Eigen::Matrix<double, 2, 1> result;

   tausqtad_(&mtau2, &mA2, &msv2, &mstau12, &mstau22, &sintau,
             &costau, &scale2, &mu, &tanb, &vev2, &result(0), &result(1));

   // workaround for intel or Eigen bug causing unexpected behaviour
   // of result.allFinite()
   if (!std::isfinite(result(0)) || !std::isfinite(result(1)))
       result.setZero();

   return -result;
}

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_at_as_mssm_with_tadpoles(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double mu,
   double tanb, double vev2, double gs, int scheme)
{
   if (std::abs((mst12 - mst22)/mst12) < 1e-8)
      return self_energy_higgs_2loop_at_as_mssm_with_tadpoles_st_0_mst1_eq_mst2(
         mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, gs, scheme);

   return self_energy_higgs_2loop_at_as_mssm_with_tadpoles_general(
      mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, gs, scheme);
}

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_at_at_mssm_with_tadpoles(
   double mt2, double mb2, double mA2, double mst12,
   double mst22, double msb12, double msb22,
   double sxt, double cxt, double sxb, double cxb,
   double scale2, double mu, double tanb, double vev2)
{
   Eigen::Matrix<double, 2, 2> result;

   {
      LOCK_MUTEX();

      ddshiggs_(&mt2, &mb2, &mA2, &mst12, &mst22, &msb12, &msb22,
                &sxt, &cxt, &sxb, &cxb, &scale2, &mu, &tanb, &vev2,
                &result(0,0), &result(0,1), &result(1,1));
   }

   result(1,0) = result(0,1);

   return -result;
}

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_ab_as_mssm_with_tadpoles(
   double mb2, double mg, double msb12, double msb22,
   double sxb, double cxb, double scale2, double mu,
   double cotb, double vev2, double gs, int scheme)
{
   Eigen::Matrix<double, 2, 2> result(self_energy_higgs_2loop_at_as_mssm_with_tadpoles(
      mb2, mg, msb12, msb22, sxb, cxb, scale2, mu,
      cotb, vev2, gs, scheme));

   std::swap(result(0,0), result(1,1));

   return result;
}

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_atau_atau_mssm_with_tadpoles(
   double mtau2, double mA2, double msv2, double mstau12,
   double mstau22, double sintau, double costau, double scale2,
   double mu, double tanb, double vev2, int scheme)
{
   Eigen::Matrix<double, 2, 2> result;

   tausqhiggs_(&mtau2, &mA2, &msv2, &mstau12, &mstau22, &sintau,
               &costau, &scale2, &mu, &tanb, &vev2, &scheme,
               &result(0,0), &result(1,1), &result(0,1));

   result(1,0) = result(0,1);

   return -result;
}

double self_energy_pseudoscalar_2loop_at_as_mssm_with_tadpoles(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double mu,
   double tanb, double vev2, double gs)
{
   if (std::abs((mst12 - mst22)/mst12) < 1e-8) {
      const double At = calc_At(mt2, mst12, mst22, sxt, cxt, mu, tanb);

      // if At = 0 => mu = 0 => dMA(2L) = 0
      if (std::abs(At) < std::numeric_limits<double>::epsilon())
         return 0.;

      return self_energy_pseudoscalar_2loop_at_as_mssm_with_tadpoles_mst1_eq_mst2(
         mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, gs);
   }

   return self_energy_pseudoscalar_2loop_at_as_mssm_with_tadpoles_general(
      mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, gs);
}

double self_energy_pseudoscalar_2loop_at_at_mssm_with_tadpoles(
   double mt2, double mb2, double mA2, double mst12,
   double mst22, double msb12, double msb22,
   double sxt, double cxt, double sxb, double cxb,
   double scale2, double mu, double tanb, double vev2)
{
   double result;

   {
      LOCK_MUTEX();

      ddsodd_(&mt2, &mb2, &mA2, &mst12, &mst22, &msb12, &msb22,
              &sxt, &cxt, &sxb, &cxb, &scale2, &mu, &tanb, &vev2, &result);
   }

   return -result;
}

double self_energy_pseudoscalar_2loop_ab_as_mssm_with_tadpoles(
   double mb2, double mg, double msb12, double msb22,
   double sxb, double cxb, double scale2, double mu,
   double cotb, double vev2, double gs)
{
   return self_energy_pseudoscalar_2loop_at_as_mssm_with_tadpoles(
      mb2, mg, msb12, msb22, sxb, cxb, scale2, mu,
      cotb, vev2, gs);
}

double self_energy_pseudoscalar_2loop_atau_atau_mssm_with_tadpoles(
   double mtau2, double mA2, double msv2, double mstau12,
   double mstau22, double sintau, double costau, double scale2,
   double mu, double tanb, double vev2)
{
   double result;

   tausqodd_(&mtau2, &mA2, &msv2, &mstau12, &mstau22, &sintau,
             &costau, &scale2, &mu, &tanb, &vev2, &result);

   return -result;
}

// self-energies without tadpoles

Eigen::Matrix<double, 2, 2> rotate_scalar(
   double self_energy, double tanb)
{
   const double tanb2 = sqr(tanb);
   const double sinb = tanb / sqrtabs(1. + tanb2);
   const double cosb = 1. / sqrtabs(1. + tanb2);

   Eigen::Matrix<double, 2, 2> result;

   result(0,0) = self_energy * sqr(sinb);
   result(0,1) = - self_energy * sinb * cosb;
   result(1,0) = result(0,1);
   result(1,1) = self_energy * sqr(cosb);

   return result;
}

Eigen::Matrix<double, 2, 2> subtract_mssm_tadpoles_scalar(
   double self_energy, const Eigen::Matrix<double, 2, 1>& tadpoles,
   double tanb)
{
   return rotate_scalar(self_energy, tanb) + Eigen::Matrix<double, 2, 2>(tadpoles.asDiagonal());
}

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_at_as_mssm(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double mu,
   double tanb, double vev2, double gs, int scheme)
{
   const Eigen::Matrix<double, 2, 2> result =
      self_energy_higgs_2loop_at_as_mssm_with_tadpoles(
         mt2, mg, mst12, mst22, sxt, cxt, scale2, mu,
         tanb, vev2, gs, scheme);

   const double dMA = self_energy_pseudoscalar_2loop_at_as_mssm_with_tadpoles(
      mt2, mg, mst12, mst22, sxt, cxt, scale2, mu,
      tanb, vev2, gs);

   const Eigen::Matrix<double, 2, 1> tadpoles =
      tadpole_higgs_2loop_at_as_mssm(
         mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, gs);

   const Eigen::Matrix<double, 2, 2> tM =
      subtract_mssm_tadpoles_scalar(dMA, tadpoles, tanb);

   return result + tM;
}

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_at_at_mssm(
   double mt2, double mb2, double mA2, double mst12,
   double mst22, double msb12, double msb22,
   double sxt, double cxt, double sxb, double cxb,
   double scale2, double mu, double tanb, double vev2)
{
   const Eigen::Matrix<double, 2, 2> result =
      self_energy_higgs_2loop_at_at_mssm_with_tadpoles(
         mt2, mb2, mA2, mst12, mst22, msb12, msb22,
         sxt, cxt, sxb, cxb, scale2, mu, tanb, vev2);

   const double dMA = self_energy_pseudoscalar_2loop_at_at_mssm_with_tadpoles(
      mt2, mb2, mA2, mst12, mst22, msb12, msb22,
      sxt, cxt, sxb, cxb, scale2, mu, tanb, vev2);

   const Eigen::Matrix<double, 2, 1> tadpoles =
      tadpole_higgs_2loop_at_at_mssm(
         mt2, mb2, mA2, mst12, mst22, msb12, msb22,
         sxt, cxt, sxb, cxb, scale2, mu, tanb, vev2);

   const Eigen::Matrix<double, 2, 2> tM =
      subtract_mssm_tadpoles_scalar(dMA, tadpoles, tanb);

   return result + tM;
}

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_ab_as_mssm(
   double mb2, double mg, double msb12, double msb22,
   double sxb, double cxb, double scale2, double mu,
   double cotb, double vev2, double gs, int scheme)
{
   const Eigen::Matrix<double, 2, 2> result =
      self_energy_higgs_2loop_ab_as_mssm_with_tadpoles(
         mb2, mg, msb12, msb22, sxb, cxb, scale2, mu,
         cotb, vev2, gs, scheme);

   const double dMA = self_energy_pseudoscalar_2loop_ab_as_mssm_with_tadpoles(
      mb2, mg, msb12, msb22, sxb, cxb, scale2, mu, cotb, vev2, gs);

   const Eigen::Matrix<double, 2, 1> tadpoles =
      tadpole_higgs_2loop_ab_as_mssm(
         mb2, mg, msb12, msb22, sxb, cxb, scale2, mu, cotb, vev2, gs);

   const Eigen::Matrix<double, 2, 2> tM =
      subtract_mssm_tadpoles_scalar(dMA, tadpoles, 1./cotb);

   return result + tM;
}

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_atau_atau_mssm(
   double mtau2, double mA2, double msv2, double mstau12,
   double mstau22, double sintau, double costau, double scale2,
   double mu, double tanb, double vev2, int scheme)
{
   const Eigen::Matrix<double, 2, 2> result =
      self_energy_higgs_2loop_atau_atau_mssm_with_tadpoles(
         mtau2, mA2, msv2, mstau12, mstau22, sintau, costau, scale2,
         mu, tanb, vev2, scheme);

   const double dMA = self_energy_pseudoscalar_2loop_atau_atau_mssm_with_tadpoles(
      mtau2, mA2, msv2, mstau12, mstau22, sintau, costau, scale2,
      mu, tanb, vev2);

   const Eigen::Matrix<double, 2, 1> tadpoles =
      tadpole_higgs_2loop_atau_atau_mssm(
         mtau2, mA2, msv2, mstau12, mstau22, sintau, costau, scale2,
         mu, tanb, vev2);

   const Eigen::Matrix<double, 2, 2> tM =
      subtract_mssm_tadpoles_scalar(dMA, tadpoles, tanb);

   return result + tM;
}

Eigen::Matrix<double, 2, 2> rotate_pseudoscalar(
   double self_energy, double tanb)
{
   const double tanb2 = sqr(tanb);
   const double sinb = tanb / sqrtabs(1. + tanb2);
   const double cosb = 1. / sqrtabs(1. + tanb2);

   Eigen::Matrix<double, 2, 2> result;

   // see hep-ph/0105096 Eq. (9)
   result(0,0) = self_energy * sqr(sinb);
   result(0,1) = self_energy * sinb * cosb;
   result(1,0) = result(0,1);
   result(1,1) = self_energy * sqr(cosb);

   return result;
}

Eigen::Matrix<double, 2, 2> subtract_mssm_tadpoles_pseudoscalar(
   double self_energy, const Eigen::Matrix<double, 2, 1>& tadpoles,
   double tanb)
{
   return rotate_pseudoscalar(self_energy, tanb) + Eigen::Matrix<double, 2, 2>(tadpoles.asDiagonal());
}

Eigen::Matrix<double, 2, 2> self_energy_pseudoscalar_2loop_at_as_mssm(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double mu,
   double tanb, double vev2, double gs)
{
   const double se = self_energy_pseudoscalar_2loop_at_as_mssm_with_tadpoles(
      mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, gs);

   const Eigen::Matrix<double, 2, 1> tadpoles =
      tadpole_higgs_2loop_at_as_mssm(
         mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, gs);

   return subtract_mssm_tadpoles_pseudoscalar(se, tadpoles, tanb);
}

Eigen::Matrix<double, 2, 2> self_energy_pseudoscalar_2loop_at_at_mssm(
   double mt2, double mb2, double mA2, double mst12,
   double mst22, double msb12, double msb22,
   double sxt, double cxt, double sxb, double cxb,
   double scale2, double mu, double tanb, double vev2)
{
   const double se = self_energy_pseudoscalar_2loop_at_at_mssm_with_tadpoles(
      mt2, mb2, mA2, mst12, mst22, msb12, msb22,
      sxt, cxt, sxb, cxb, scale2, mu, tanb, vev2);

   const Eigen::Matrix<double, 2, 1> tadpoles =
      tadpole_higgs_2loop_at_at_mssm(
         mt2, mb2, mA2, mst12, mst22, msb12, msb22,
         sxt, cxt, sxb, cxb, scale2, mu, tanb, vev2);

   return subtract_mssm_tadpoles_pseudoscalar(se, tadpoles, tanb);
}

Eigen::Matrix<double, 2, 2> self_energy_pseudoscalar_2loop_ab_as_mssm(
   double mb2, double mg, double msb12, double msb22,
   double sxb, double cxb, double scale2, double mu,
   double cotb, double vev2, double gs)
{
   const double se = self_energy_pseudoscalar_2loop_ab_as_mssm_with_tadpoles(
      mb2, mg, msb12, msb22, sxb, cxb, scale2, mu, cotb, vev2, gs);

   const Eigen::Matrix<double, 2, 1> tadpoles =
      tadpole_higgs_2loop_ab_as_mssm(
         mb2, mg, msb12, msb22, sxb, cxb, scale2, mu, cotb, vev2, gs);

   return subtract_mssm_tadpoles_pseudoscalar(se, tadpoles, 1./cotb);
}

Eigen::Matrix<double, 2, 2> self_energy_pseudoscalar_2loop_atau_atau_mssm(
   double mtau2, double mA2, double msv2, double mstau12,
   double mstau22, double sintau, double costau, double scale2,
   double mu, double tanb, double vev2)
{
   const double se = self_energy_pseudoscalar_2loop_atau_atau_mssm_with_tadpoles(
      mtau2, mA2, msv2, mstau12, mstau22, sintau, costau, scale2,
      mu, tanb, vev2);

   const Eigen::Matrix<double, 2, 1> tadpoles =
      tadpole_higgs_2loop_atau_atau_mssm(
         mtau2, mA2, msv2, mstau12, mstau22, sintau, costau, scale2,
         mu, tanb, vev2);

   return subtract_mssm_tadpoles_pseudoscalar(se, tadpoles, tanb);
}

} // namespace mssm_twoloophiggs
} // namespace flexiblesusy
