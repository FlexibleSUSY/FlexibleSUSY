#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_lowNMSSM_2loop_singlet_sign

#include <boost/test/unit_test.hpp>

#include "lowNMSSM_mass_eigenstates.hpp"

#include <Eigen/Core>

#include <cmath>
#include <numbers>

using namespace flexiblesusy;

namespace
{

lowNMSSM_mass_eigenstates make_model(double sign)
{
   const double tanb = 5., vev = 246., cb = 1. / std::sqrt(1. + tanb * tanb);
   const double lambda = sign * 0.3, kappa = sign * 0.2, mueff = 300.;
   const double Alambda = 600., Akappa = -100., At = 2000.;
   const double ms2 = 1e6;

   lowNMSSM_mass_eigenstates m;
   m.set_scale(1000.);
   m.set_loops(2);
   m.set_g1(0.47);
   m.set_g2(0.64);
   m.set_g3(1.05);
   m.set_Yu(2, 2, 0.85);
   m.set_Yd(2, 2, 0.08);
   m.set_Ye(2, 2, 0.05);
   m.set_TYu(2, 2, 0.85 * At);
   m.set_TYd(2, 2, 0.08 * 1000.);
   m.set_TYe(2, 2, 0.05 * 1000.);
   m.set_vd(vev * cb);
   m.set_vu(vev * tanb * cb);
   m.set_Lambdax(lambda);
   m.set_Kappa(kappa);
   m.set_vS(std::numbers::sqrt2 * mueff /
            lambda); // EffectiveMu = lambda vS/sqrt(2) > 0
   m.set_TLambdax(lambda * Alambda);
   m.set_TKappa(kappa * Akappa);
   for (int i = 0; i < 3; ++i) {
      m.set_mq2(i, i, ms2);
      m.set_mu2(i, i, ms2);
      m.set_md2(i, i, ms2);
      m.set_ml2(i, i, ms2);
      m.set_me2(i, i, ms2);
   }
   m.set_mHd2(1e5);
   m.set_mHu2(-1e5);
   m.set_ms2(1e5);
   m.set_MassB(500.);
   m.set_MassWB(1000.);
   m.set_MassG(2000.);
   m.calculate_DRbar_masses();
   return m;
}

double max_rel_diff(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b)
{
   return (a - b).cwiseAbs().maxCoeff() / a.cwiseAbs().maxCoeff();
}

} // anonymous namespace

// The NMSSM is invariant under S -> -S with (lambda, kappa, vS, TLambdax,
// TKappa) -> -(lambda, kappa, vS, TLambdax, TKappa). The 2-loop Higgs
// tadpoles and self-energies must transform accordingly: the singlet
// components change sign, the rest is unchanged.
BOOST_AUTO_TEST_CASE(test_singlet_sign_flip)
{
   const auto m = make_model(+1.);
   const auto f = make_model(-1.);

   Eigen::Matrix3d p = Eigen::Matrix3d::Identity();
   p(2, 2) = -1.;

   BOOST_CHECK_SMALL(
      max_rel_diff(m.tadpole_hh_2loop(), p * f.tadpole_hh_2loop()), 1e-10);
   BOOST_CHECK_SMALL(
      max_rel_diff(m.self_energy_hh_2loop(), p * f.self_energy_hh_2loop() * p),
      1e-10);
   BOOST_CHECK_SMALL(
      max_rel_diff(m.self_energy_Ah_2loop(), p * f.self_energy_Ah_2loop() * p),
      1e-10);
}
