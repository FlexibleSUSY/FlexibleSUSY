#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SMSSM_singlet_shift

#include <boost/test/unit_test.hpp>

#include "SMSSM_input_parameters.hpp"
#include "SMSSM_mass_eigenstates.hpp"

#include <Eigen/Core>

#include <cmath>
#include <numbers>

using namespace flexiblesusy;

namespace
{

// parameters of the point in the frame with vS = vS0
constexpr double vS0 = 600.;
constexpr double lambda = 0.3, kappa = 0.2;
constexpr double MS0 = 100., L10 = 1e4, BMS0 = 1e4;
constexpr double TLambdax = lambda * 500., TKappa = kappa * -200.;
constexpr double ms2 = 1e5;

/**
 * Creates the SMSSM at the point above after the shift S -> S + delta
 * of the singlet superfield, which gives vS = vS0 - sqrt(2) delta,
 * and solves the 1-loop EWSB equations for mu, B[mu] and L[L1].
 *
 * The shift is an exact symmetry of the SMSSM: with
 *
 *   mu    -> mu + lambda delta
 *   MS    -> MS + 2 kappa delta
 *   L1    -> L1 + MS delta + kappa delta^2
 *   B[mu] -> B[mu] + T[lambda] delta
 *   B[MS] -> B[MS] + 2 T[kappa] delta
 *   L[L1] -> L[L1] + B[MS] delta + T[kappa] delta^2 + ms2 delta
 *
 * the Lagrangian expanded around the vacuum is unchanged.  mu, B[mu]
 * and L[L1] are fixed by the EWSB.
 */
SMSSM_mass_eigenstates make_model(double vS)
{
   const double delta = (vS0 - vS) / std::numbers::sqrt2;
   const double tanb = 5., vev = 246., cb = 1. / std::sqrt(1. + tanb * tanb);
   const double msf2 = 1e6, At = 2000.;

   SMSSM_input_parameters input;
   input.SignMu = 1; // sign of mu_eff = mu + lambda vS/sqrt(2)

   SMSSM_mass_eigenstates m(input);
   m.set_scale(1000.);
   m.set_loops(1);
   m.set_ewsb_iteration_precision(1e-12);
   m.set_g1(0.47);
   m.set_g2(0.64);
   m.set_g3(1.05);
   m.set_Yu(2, 2, 0.85);
   m.set_Yd(2, 2, 0.08);
   m.set_Ye(2, 2, 0.05);
   m.set_TYu(2, 2, 0.85 * At);
   m.set_TYd(2, 2, 0.08 * 1000.);
   m.set_TYe(2, 2, 0.05 * 1000.);
   for (int i = 0; i < 3; ++i) {
      m.set_mq2(i, i, msf2);
      m.set_mu2(i, i, msf2);
      m.set_md2(i, i, msf2);
      m.set_ml2(i, i, msf2);
      m.set_me2(i, i, msf2);
   }
   m.set_mHd2(2e5);
   m.set_mHu2(-1e5);
   m.set_ms2(ms2);
   m.set_MassB(500.);
   m.set_MassWB(1000.);
   m.set_MassG(2000.);
   m.set_vd(vev * cb);
   m.set_vu(vev * tanb * cb);
   m.set_Lambdax(lambda);
   m.set_Kappa(kappa);
   m.set_TLambdax(TLambdax);
   m.set_TKappa(TKappa);
   m.set_vS(vS);
   m.set_MS(MS0 + 2 * kappa * delta);
   m.set_L1(L10 + MS0 * delta + kappa * delta * delta);
   m.set_BMS(BMS0 + 2 * TKappa * delta);

   BOOST_REQUIRE_EQUAL(m.solve_ewsb(), 0);
   m.calculate_DRbar_masses();
   BOOST_REQUIRE(!m.get_problems().have_problem());

   m.calculate_Mhh_pole();
   m.calculate_MAh_pole();
   m.calculate_MHpm_pole();
   m.calculate_MChi_pole();
   m.calculate_MCha_pole();
   BOOST_REQUIRE(!m.get_problems().have_problem());

   return m;
}

double max_rel_diff(const Eigen::ArrayXd& a, const Eigen::ArrayXd& b)
{
   return ((a - b).abs() / a.abs().max(1.)).maxCoeff();
}

void check_same_spectrum(const SMSSM_mass_eigenstates& a,
                         const SMSSM_mass_eigenstates& b, double eps)
{
   const auto& pa = a.get_physical();
   const auto& pb = b.get_physical();

   BOOST_TEST_MESSAGE("vS = " << a.get_vS() << " vs. " << b.get_vS());
   BOOST_TEST_MESSAGE("Mhh: " << pa.Mhh.transpose() << " vs. "
                              << pb.Mhh.transpose());
   BOOST_TEST_MESSAGE("MAh: " << pa.MAh.transpose() << " vs. "
                              << pb.MAh.transpose());

   BOOST_CHECK_SMALL(max_rel_diff(pa.Mhh, pb.Mhh), eps);
   BOOST_CHECK_SMALL(max_rel_diff(pa.MAh, pb.MAh), eps);
   BOOST_CHECK_SMALL(max_rel_diff(pa.MHpm, pb.MHpm), eps);
   BOOST_CHECK_SMALL(max_rel_diff(pa.MChi, pb.MChi), eps);
   BOOST_CHECK_SMALL(max_rel_diff(pa.MCha, pb.MCha), eps);
}

} // anonymous namespace

// The pole masses must not depend on the frame of the singlet shift,
// also for small and negative vS.  Before, the running masses inside
// the loops were calculated with ms2 from the tree-level EWSB
// equations, which differs from the loop-level vacuum by t_S/vS.
BOOST_AUTO_TEST_CASE(test_singlet_shift_invariance)
{
   const auto m0 = make_model(vS0);

   for (const double vS : {300., 50., 20., -20., -300.}) {
      check_same_spectrum(m0, make_model(vS), 1e-8);
   }
}
