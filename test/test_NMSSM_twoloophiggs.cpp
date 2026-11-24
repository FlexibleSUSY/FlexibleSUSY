#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_NMSSM_twoloophiggs

#include <boost/test/unit_test.hpp>

#include "mssm_twoloophiggs.hpp"
#include "nmssm_twoloophiggs.hpp"

#include <Eigen/Core>

#include <cmath>
#include <numbers>

using namespace flexiblesusy;

namespace
{

struct Point {
   double tanb = 39.;
   double vev2 = 246. * 246.;
   double scale2 = 1000. * 1000.;
   double as = 0.09;
   double mg = 2200.;
   double mt = 150., mst12 = 1100. * 1100., mst22 = 1600. * 1600., sxt = 0.3;
   double mb = 2.4, msb12 = 1400. * 1400., msb22 = 1700. * 1700., sxb = 0.2;
   double lam = 0.3;
   double mu =
      -150.; ///< mu in the convention of the routines (= -mu_eff of SARAH)
   [[nodiscard]] double svev() const
   {
      return -mu / lam;
   } ///< v_S/sqrt(2), with its sign
};

Eigen::Matrix3d pseudoscalar_at_as(const Point& p)
{
   return nmssm_twoloophiggs::
      self_energy_pseudoscalar_2loop_at_as_nmssm_with_tadpoles(
         p.mt, p.mg, p.mst12, p.mst22, p.sxt, std::sqrt(1. - p.sxt * p.sxt),
         p.scale2, p.tanb, p.vev2, p.lam, p.svev(), p.as);
}

Eigen::Matrix3d pseudoscalar_ab_as(const Point& p)
{
   return nmssm_twoloophiggs::
      self_energy_pseudoscalar_2loop_ab_as_nmssm_with_tadpoles(
         p.mb, p.mg, p.msb12, p.msb22, p.sxb, std::sqrt(1. - p.sxb * p.sxb),
         p.scale2, 1. / p.tanb, p.vev2, p.lam, p.svev(), p.as);
}

Eigen::Matrix3d higgs_ab_as(const Point& p)
{
   return nmssm_twoloophiggs::self_energy_higgs_2loop_ab_as_nmssm_with_tadpoles(
      p.mb, p.mg, p.msb12, p.msb22, p.sxb, std::sqrt(1. - p.sxb * p.sxb),
      p.scale2, 1. / p.tanb, p.vev2, p.lam, p.svev(), p.as);
}

/// neutral Goldstone direction in the basis (Im H_d^0, Im H_u^0, Im S)
Eigen::Vector3d goldstone(double tanb)
{
   const double cb = 1. / std::sqrt(1. + tanb * tanb);
   return {cb, -tanb * cb, 0.};
}

/// heavy CP-odd doublet direction A_NSM in the same basis
Eigen::Vector3d a_nsm(double tanb)
{
   const double cb = 1. / std::sqrt(1. + tanb * tanb);
   return {tanb * cb, cb, 0.};
}

double rel_asymmetry(const Eigen::Matrix3d& m)
{
   return (m - m.transpose()).cwiseAbs().maxCoeff() / m.cwiseAbs().maxCoeff();
}

} // anonymous namespace

// The O(ab as) results are obtained from the O(at as) ones with H_u <-> H_d
// exchanged; the exchange must keep the matrices symmetric.
BOOST_AUTO_TEST_CASE(test_ab_as_symmetric)
{
   for (const double lam : {0.3, -0.3}) {
      Point p;
      p.lam = lam;
      BOOST_TEST_CONTEXT("lambda = " << lam)
      {
         BOOST_CHECK_SMALL(rel_asymmetry(pseudoscalar_ab_as(p)), 1e-14);
         BOOST_CHECK_SMALL(rel_asymmetry(higgs_ab_as(p)), 1e-14);
      }
   }
}

// The tadpole-subtracted CP-odd mass matrix derived from the effective
// potential must have the Goldstone boson as a zero mode.
BOOST_AUTO_TEST_CASE(test_goldstone_zero_mode)
{
   for (const double lam : {0.3, -0.3}) {
      for (const double mu : {150., -150.}) {
         Point p;
         p.lam = lam;
         p.mu = mu;
         const Eigen::Vector3d g = goldstone(p.tanb);
         const Eigen::Matrix3d mt = pseudoscalar_at_as(p);
         const Eigen::Matrix3d mb = pseudoscalar_ab_as(p);
         BOOST_TEST_CONTEXT("lambda = " << lam << ", mu = " << mu)
         {
            BOOST_CHECK_SMALL((mt * g).norm() / mt.cwiseAbs().maxCoeff(),
                              1e-10);
            BOOST_CHECK_SMALL((mb * g).norm() / mb.cwiseAbs().maxCoeff(),
                              1e-10);
            // the lower triangle (read by the eigenvalue solvers) alone
            const Eigen::Matrix3d mb_lower = mb.selfadjointView<Eigen::Lower>();
            BOOST_CHECK_SMALL((mb_lower * g).norm() / mb.cwiseAbs().maxCoeff(),
                              1e-10);
         }
      }
   }
}

// For lambda -> 0 at fixed mu, the doublet part must reproduce the MSSM
// result with the same mu. This fixes the argument svev = -mu/lambda
// (= v_S/sqrt(2) with its sign, as FlexibleSUSY passes mu = -mu_eff).
BOOST_AUTO_TEST_CASE(test_mssm_limit)
{
   for (const double lam : {1e-4, -1e-4}) {
      for (const double mu : {150., -150.}) {
         Point p;
         p.lam = lam;
         p.mu = mu;
         p.tanb = 5.;
         const double gs = std::sqrt(4. * std::numbers::pi * p.as);
         const Eigen::Vector3d a = a_nsm(p.tanb);

         const double dma2_t = mssm_twoloophiggs::
            self_energy_pseudoscalar_2loop_at_as_mssm_with_tadpoles(
               p.mt * p.mt, p.mg, p.mst12, p.mst22, p.sxt,
               std::sqrt(1. - p.sxt * p.sxt), p.scale2, p.mu, p.tanb, p.vev2,
               gs);
         const double dma2_b = mssm_twoloophiggs::
            self_energy_pseudoscalar_2loop_ab_as_mssm_with_tadpoles(
               p.mb * p.mb, p.mg, p.msb12, p.msb22, p.sxb,
               std::sqrt(1. - p.sxb * p.sxb), p.scale2, p.mu, 1. / p.tanb,
               p.vev2, gs);

         BOOST_TEST_CONTEXT("lambda = " << lam << ", mu = " << mu)
         {
            BOOST_CHECK_CLOSE_FRACTION(a.dot(pseudoscalar_at_as(p) * a), dma2_t,
                                       1e-6);
            BOOST_CHECK_CLOSE_FRACTION(a.dot(pseudoscalar_ab_as(p) * a), dma2_b,
                                       1e-6);
         }
      }
   }
}
