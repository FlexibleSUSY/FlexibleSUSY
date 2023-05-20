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

#include "nmssm_twoloophiggs.hpp"
#include "mssm_twoloophiggs.hpp"
#include "nmssm2loop.h"
#include <cmath>
#include <utility>

using namespace flexiblesusy::mssm_twoloophiggs;

namespace flexiblesusy {
namespace nmssm_twoloophiggs {

namespace {
template <typename T> T sqr(T a) { return a * a; }
} // namespace

Eigen::Matrix<double, 3, 1> tadpole_higgs_2loop_at_as_nmssm(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2,
   double mu, double tanb, double vev2, double gs, double svev)
{
   const double cosb = 1. / std::sqrt(1. + sqr(tanb));

   const Eigen::Matrix<double, 2, 1> t_mssm = tadpole_higgs_2loop_at_as_mssm(
      mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, gs);

   const double vd = std::sqrt(vev2) * cosb;
   const double vs = svev * std::sqrt(2.);

   Eigen::Matrix<double, 3, 1> result;
   result.head<2>() = t_mssm;
   // rescale T1 to get TS
   result(2) = t_mssm(0) * sqr(vd/vs);

   return result;
}

Eigen::Matrix<double, 3, 1> tadpole_higgs_2loop_ab_as_nmssm(
   double mb2, double mg, double msb12, double msb22,
   double sxb, double cxb, double scale2,
   double mu, double cotb, double vev2, double gs, double svev)
{
   const double tanb = 1./cotb;
   const double sinb = tanb / std::sqrt(1. + sqr(tanb));

   const Eigen::Matrix<double, 2, 1> t_mssm = tadpole_higgs_2loop_ab_as_mssm(
      mb2, mg, msb12, msb22, sxb, cxb, scale2, mu, cotb, vev2, gs);

   const double vu = std::sqrt(vev2) * sinb;
   const double vs = svev * std::sqrt(2.);

   Eigen::Matrix<double, 3, 1> result;
   result.head<2>() = t_mssm;
   // rescale T1 to get TS
   result(2) = t_mssm(0) * sqr(vu/vs);

   return result;
}

Eigen::Matrix<double, 3, 3> self_energy_higgs_2loop_at_as_nmssm(
   double rmt, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double tanb, double vev2,
   double lam, double svev, double as, double mu)
{
   const double gs = std::sqrt(as * 4. * M_PI);

   const Eigen::Matrix<double, 3, 3> se =
      self_energy_higgs_2loop_at_as_nmssm_with_tadpoles(
         rmt, mg, mst12, mst22, sxt, cxt, scale2, tanb, vev2, lam, svev, as);

   const Eigen::Matrix<double, 3, 3> tadpoles =
      tadpole_higgs_2loop_at_as_nmssm(
         sqr(rmt), mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, gs, svev).asDiagonal();

   return se + tadpoles;
}

Eigen::Matrix<double, 3, 3> self_energy_higgs_2loop_ab_as_nmssm(
   double rmb, double mg, double msb12, double msb22,
   double sxb, double cxb, double scale2, double cotb, double vev2,
   double lam, double svev, double as, double mu)
{
   const double gs = std::sqrt(as * 4. * M_PI);

   const Eigen::Matrix<double, 3, 3> se =
      self_energy_higgs_2loop_ab_as_nmssm_with_tadpoles(
         rmb, mg, msb12, msb22, sxb, cxb, scale2, cotb, vev2, lam, svev, as);

   const Eigen::Matrix<double, 3, 3> tadpoles =
      tadpole_higgs_2loop_ab_as_nmssm(
         sqr(rmb), mg, msb12, msb22, sxb, cxb, scale2, mu, cotb, vev2, gs, svev).asDiagonal();

   return se + tadpoles;
}

Eigen::Matrix<double, 3, 3> self_energy_pseudoscalar_2loop_at_as_nmssm(
   double rmt, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double tanb, double vev2,
   double lam, double svev, double as, double mu)
{
   const double gs = std::sqrt(as * 4. * M_PI);

   const Eigen::Matrix<double, 3, 3> se =
      self_energy_pseudoscalar_2loop_at_as_nmssm_with_tadpoles(
         rmt, mg, mst12, mst22, sxt, cxt, scale2, tanb, vev2,
         lam, svev, as);

   const Eigen::Matrix<double, 3, 3> tadpoles =
      tadpole_higgs_2loop_at_as_nmssm(
         sqr(rmt), mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, gs, svev).asDiagonal();

   return se + tadpoles;
}

Eigen::Matrix<double, 3, 3> self_energy_pseudoscalar_2loop_ab_as_nmssm(
   double rmb, double mg, double msb12, double msb22,
   double sxb, double cxb, double scale2, double cotb, double vev2,
   double lam, double svev, double as, double mu)
{
   const double gs = std::sqrt(as * 4. * M_PI);

   const Eigen::Matrix<double, 3, 3> se =
      self_energy_pseudoscalar_2loop_ab_as_nmssm_with_tadpoles(
         rmb, mg, msb12, msb22, sxb, cxb, scale2, cotb, vev2,
         lam, svev, as);

   const Eigen::Matrix<double, 3, 3> tadpoles =
      tadpole_higgs_2loop_ab_as_nmssm(
         sqr(rmb), mg, msb12, msb22, sxb, cxb, scale2, mu, cotb, vev2, gs, svev).asDiagonal();

   return se + tadpoles;
}

Eigen::Matrix<double, 3, 3> self_energy_higgs_2loop_at_as_nmssm_with_tadpoles(
   double rmt, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double tanb, double vev2,
   double lam, double svev, double as)
{
   int loop = 2;
   double vev = std::sqrt(0.5 * vev2);
   double DMS[3][3] = {{ 0. }}, DMP[3][3] = {{ 0. }};

   effpot_(&loop, &rmt, &mg, &mst12, &mst22, &sxt, &cxt,
           &scale2, &tanb, &vev, &lam, &svev, &as, &DMS, &DMP);

   Eigen::Matrix<double, 3, 3> result;
   result << DMS[0][0], DMS[0][1], DMS[0][2],
             DMS[1][0], DMS[1][1], DMS[1][2],
             DMS[2][0], DMS[2][1], DMS[2][2];

   return -result;
}

Eigen::Matrix<double, 3, 3> self_energy_higgs_2loop_ab_as_nmssm_with_tadpoles(
   double rmb, double mg, double msb12, double msb22,
   double sxb, double cxb, double scale2, double cotb, double vev2,
   double lam, double svev, double as)
{
   Eigen::Matrix<double, 3, 3> result =
      self_energy_higgs_2loop_at_as_nmssm_with_tadpoles(
         rmb, mg, msb12, msb22, sxb, cxb, scale2, cotb, vev2,
         lam, svev, as);

   // Make appropriate substitutions for elements following 0907.4682
   // bottom of page 9
   std::swap(result(0,0), result(1,1));
   std::swap(result(0,2), result(1,2));

   return result;
}

Eigen::Matrix<double, 3, 3> self_energy_pseudoscalar_2loop_at_as_nmssm_with_tadpoles(
   double rmt, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double tanb, double vev2,
   double lam, double svev, double as)
{
   int loop = 2;
   double vev = std::sqrt(0.5 * vev2);
   double DMS[3][3] = {{ 0. }}, DMP[3][3] = {{ 0. }};

   effpot_(&loop, &rmt, &mg, &mst12, &mst22, &sxt, &cxt,
           &scale2, &tanb, &vev, &lam, &svev, &as, &DMS, &DMP);

   Eigen::Matrix<double, 3, 3> result;
   result << DMP[0][0], DMP[0][1], DMP[0][2],
             DMP[1][0], DMP[1][1], DMP[1][2],
             DMP[2][0], DMP[2][1], DMP[2][2];

   return -result;
}

Eigen::Matrix<double, 3, 3> self_energy_pseudoscalar_2loop_ab_as_nmssm_with_tadpoles(
   double rmb, double mg, double msb12, double msb22,
   double sxb, double cxb, double scale2, double cotb, double vev2,
   double lam, double svev, double as)
{
   Eigen::Matrix<double, 3, 3> result =
      self_energy_pseudoscalar_2loop_at_as_nmssm_with_tadpoles(
         rmb, mg, msb12, msb22, sxb, cxb, scale2, cotb, vev2,
         lam, svev, as);

   // Make appropriate substitutions for elements following 0907.4682
   // bottom of page 9
   std::swap(result(0,0), result(1,1));
   std::swap(result(0,2), result(1,2));

   return result;
}

} // namespace nmssm_twoloophiggs
} // namespace flexiblesusy
