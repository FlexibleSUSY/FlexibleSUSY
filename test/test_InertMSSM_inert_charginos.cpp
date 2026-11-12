#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_InertMSSM_inert_charginos

#include <boost/test/unit_test.hpp>

#include "InertMSSM_mass_eigenstates.hpp"

#include <Eigen/SVD>

#include <cmath>
#include <complex>

using namespace flexiblesusy;

namespace
{

InertMSSM_mass_eigenstates make_model()
{
   InertMSSM_mass_eigenstates m;
   m.set_g1(0.47);
   m.set_g2(0.65);
   m.set_vd(25.);
   m.set_vu(243.);

   // non-diagonal and non-symmetric, so that UIM != UIP
   Eigen::Matrix<double, 2, 2> mu_prime;
   mu_prime << 300., 150., -40., 700.;
   m.set_MuPrime(mu_prime);

   m.calculate_MVPVZ();
   m.calculate_MInertCha();

   return m;
}

} // anonymous namespace

// The charged inert higgsinos form two Dirac fermions of charge -1
// (issue #563), with the mass matrix that SARAH derives: mu'^T.
BOOST_AUTO_TEST_CASE(test_InertCha_mass_matrix)
{
   const auto m = make_model();

   BOOST_CHECK_EQUAL(m.get_MInertCha().size(), 2);

   const Eigen::Matrix<double, 2, 2> mm = m.get_mass_matrix_InertCha();
   const Eigen::Matrix<double, 2, 2> mu_prime_T = m.get_MuPrime().transpose();

   for (int i = 0; i < 2; i++) {
      for (int k = 0; k < 2; k++) {
         BOOST_CHECK_EQUAL(mm(i, k), mu_prime_T(i, k));
      }
   }
}

BOOST_AUTO_TEST_CASE(test_InertCha_masses)
{
   const auto m = make_model();

   // singular values of mu' in ascending order
   const Eigen::JacobiSVD<Eigen::Matrix<double, 2, 2>> svd(m.get_MuPrime());
   const auto sv = svd.singularValues();

   BOOST_CHECK_CLOSE_FRACTION(m.get_MInertCha(0), sv(1), 1e-14);
   BOOST_CHECK_CLOSE_FRACTION(m.get_MInertCha(1), sv(0), 1e-14);

   // UIM^* M UIP^dagger = diag(MInertCha)
   const Eigen::Matrix<std::complex<double>, 2, 2> diag =
      m.get_UIM().conjugate() *
      m.get_mass_matrix_InertCha().cast<std::complex<double>>() *
      m.get_UIP().adjoint();

   for (int i = 0; i < 2; i++) {
      for (int k = 0; k < 2; k++) {
         const double expected = i == k ? m.get_MInertCha(i) : 0.;
         BOOST_CHECK_SMALL(std::abs(diag(i, k) - expected), 1e-10);
      }
   }
}

// A Dirac fermion of charge -1 couples vectorially and diagonally to the
// photon, with strength e.
BOOST_AUTO_TEST_CASE(test_InertCha_photon_coupling)
{
   const auto m = make_model();

   const double gY = std::sqrt(3. / 5.) * m.get_g1();
   const double g2 = m.get_g2();
   const double e = gY * g2 / std::sqrt(gY * gY + g2 * g2);

   for (int i = 0; i < 2; i++) {
      for (int k = 0; k < 2; k++) {
         const double expected = i == k ? e : 0.;
         BOOST_CHECK_SMALL(m.CpbarInertChaInertChaVPPL(i, k) - expected, 1e-14);
         BOOST_CHECK_SMALL(m.CpbarInertChaInertChaVPPR(i, k) - expected, 1e-14);
      }
   }
}
