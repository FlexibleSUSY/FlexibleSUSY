
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_ckm

#include <boost/test/unit_test.hpp>

#include "config.h"

#include "test.hpp"
#include "ckm.hpp"
#include "wrappers.hpp"
#include "linalg2.hpp"
#include "random_matrix.hpp"

#include <algorithm>

#ifdef ENABLE_RANDOM
#include <random>
#endif

using namespace flexiblesusy;

bool is_orthogonal(const Eigen::Matrix<double,3,3>& m,
                   double max_deviation = std::numeric_limits<double>::epsilon())
{
   const Eigen::Matrix<double,3,3> unity(m * m.transpose());

   bool result = is_equal(unity, Eigen::Matrix<double,3,3>::Identity(), max_deviation);

   if (!result)
      BOOST_TEST_MESSAGE("matrix is not orthogonal!"
                    "\n   Matrix = " << m <<
                    "\n   Matrix * Matrix^T = " << unity <<
                    "\n   Maximum allowed max_deviation from unity = " << max_deviation);

   return result;
}

bool is_unitary(const Eigen::Matrix<double,3,3>& m,
                double max_deviation = std::numeric_limits<double>::epsilon())
{
   return is_orthogonal(m, max_deviation);
}

bool is_unitary(const Eigen::Matrix<std::complex<double>,3,3>& m,
                double max_deviation = std::numeric_limits<double>::epsilon())
{
   const Eigen::Matrix<std::complex<double>,3,3> unity(m * m.adjoint());

   bool result = is_equal(unity, Eigen::Matrix<std::complex<double>,3,3>::Identity(), max_deviation);

   if (!result)
      BOOST_TEST_MESSAGE("matrix is not unitary!"
                    "\n   Matrix = " << m <<
                    "\n   Matrix * Matrix^+ = " << unity <<
                    "\n   Maximum allowed max_deviation from unity = " << max_deviation);

   return result;
}

BOOST_AUTO_TEST_CASE( test_CKM_unitarity_from_wolfenstein )
{
   CKM_parameters ckm_pars;

   const double lambda = 0.2272;
   const double A      = 0.818;
   const double rhobar = 0.221;
   const double etabar = 0.34;

   ckm_pars.set_from_wolfenstein(lambda, A, rhobar, etabar);

   const Eigen::Matrix<double,3,3> ckm(ckm_pars.get_real_ckm());

   BOOST_CHECK(is_unitary(ckm));
}

BOOST_AUTO_TEST_CASE( test_CKM_unitarity_from_angles )
{
   CKM_parameters ckm_pars;

   ckm_pars.theta_12 = 0.229206;
   ckm_pars.theta_13 = 0.003960;
   ckm_pars.theta_23 = 0.042223;
   ckm_pars.delta    = 0.1;

   const Eigen::Matrix<double,3,3> ckm_real(ckm_pars.get_real_ckm());
   BOOST_CHECK(is_unitary(ckm_real));

   const Eigen::Matrix<std::complex<double>,3,3> ckm_cmpl(ckm_pars.get_complex_ckm());
   BOOST_CHECK(is_unitary(ckm_cmpl));
}

BOOST_AUTO_TEST_CASE( test_real_CKM_pdg_convention )
{
   // fermion mass matrices
   const Eigen::Matrix<double,3,3>
      mu(Eigen::Matrix<double,3,3>::Random()),
      md(Eigen::Matrix<double,3,3>::Random());

   // mass eigenvalues
   Eigen::Array<double,3,1> su, sd;

   // mixing matrices
   Eigen::Matrix<double,3,3> vu, vd, uu, ud;

   fs_svd(mu, su, uu, vu);
   fs_svd(md, sd, ud, vd);

   BOOST_CHECK(is_equal(mu, uu.transpose() * su.matrix().asDiagonal() * vu, 1.e-10));
   BOOST_CHECK(is_equal(md, ud.transpose() * sd.matrix().asDiagonal() * vd, 1.e-10));

   Eigen::Matrix<double,3,3> ckm(vu*vd.adjoint());

   // transpose in order to make (0,1), (1,2) negative
   ckm.transposeInPlace();
   BOOST_CHECK(is_unitary(ckm, 1.e-10));

   CKM_parameters::to_pdg_convention(ckm, vu, vd, uu, ud);

   BOOST_CHECK(ckm(0,0) > 0.);
   BOOST_CHECK(ckm(1,1) > 0.);
   BOOST_CHECK(ckm(2,2) > 0.);

   BOOST_CHECK(ckm(0,1) > 0.);
   // BOOST_CHECK(ckm(0,2) < 0.); // 13 element is not made positive
   BOOST_CHECK(ckm(1,2) > 0.);

   BOOST_CHECK(is_equal(mu, uu.transpose() * su.matrix().asDiagonal() * vu, 1.e-10));
   BOOST_CHECK(is_equal(md, ud.transpose() * sd.matrix().asDiagonal() * vd, 1.e-10));

   {
      // check that converted matrices are consistent
      const Eigen::Matrix<double,3,3> ckm_check(vu*vd.adjoint());
      for (int i = 0; i < 3; i++) {
         for (int k = 0; k < 3; k++) {
            // check for transposed equality here, because the CKM
            // matrix from above was transposed to generate negative
            // (0,1), (1,2) elements
            BOOST_CHECK_CLOSE_FRACTION(Re(ckm(i,k)), Re(ckm_check(k,i)), 1.0e-10);
            BOOST_CHECK_CLOSE_FRACTION(Im(ckm(i,k)), Im(ckm_check(k,i)), 1.0e-10);
         }
      }
   }
}

BOOST_AUTO_TEST_CASE( test_complex_CKM_pdg_convention )
{
   // fermion mass matrices
   const Eigen::Matrix<std::complex<double>,3,3>
      mu(Eigen::Matrix<std::complex<double>,3,3>::Random()),
      md(Eigen::Matrix<std::complex<double>,3,3>::Random());

   // mass eigenvalues
   Eigen::Array<double,3,1> su, sd;

   // mixing matrices
   Eigen::Matrix<std::complex<double>,3,3> vu, vd, uu, ud;

   fs_svd(mu, su, uu, vu);
   fs_svd(md, sd, ud, vd);

   BOOST_CHECK(is_equal(mu, uu.transpose() * su.matrix().asDiagonal() * vu, 1.e-10));
   BOOST_CHECK(is_equal(md, ud.transpose() * sd.matrix().asDiagonal() * vd, 1.e-10));

   Eigen::Matrix<std::complex<double>,3,3> ckm(vu*vd.adjoint());
   const Eigen::Matrix<double,3,3> ckm_squared_invariants(ckm.cwiseAbs2());

   CKM_parameters::to_pdg_convention(ckm, vu, vd, uu, ud);

   const Eigen::Matrix<double,3,3> ckm_squared_invariants_pdg(ckm.cwiseAbs2());

   BOOST_CHECK(is_unitary(ckm, 1.e-10));

   BOOST_CHECK(is_equal(mu, uu.transpose() * su.matrix().asDiagonal() * vu, 1.e-10));
   BOOST_CHECK(is_equal(md, ud.transpose() * sd.matrix().asDiagonal() * vd, 1.e-10));

   // check signs of sij & cij
   BOOST_CHECK_LT(std::abs(std::arg(ckm(0,0))), 1e-15);
   BOOST_CHECK_LT(std::abs(std::arg(ckm(0,1))), 1e-15);
   BOOST_CHECK_LT(std::abs(std::arg(ckm(1,2))), 1e-15);
   BOOST_CHECK_LT(std::abs(std::arg(ckm(2,2))), 1e-15);
   BOOST_CHECK((ckm.bottomLeftCorner<2,2>().imag().array() <= 0).all() ||
	       (ckm.bottomLeftCorner<2,2>().imag().array() >= 0).all());

   BOOST_CHECK(is_equal(ckm_squared_invariants, ckm_squared_invariants_pdg, 1.e-12));
}

#ifdef ENABLE_RANDOM

std::complex<double> random_phase()
{
   return std::polar(1., 2. * M_PI * rand() / RAND_MAX);
}

BOOST_AUTO_TEST_CASE( test_complex_CKM_pdg_convention_zero_c13 )
{
   Eigen::Matrix<std::complex<double>,2,2> ckmBL;
   Eigen::Matrix<std::complex<double>,3,3> vd;

   std::mt19937 generator;
   random_cue_matrix(ckmBL, generator);
   random_cue_matrix(vd, generator);

   Eigen::Matrix<std::complex<double>,3,3> ckm(
      Eigen::Matrix<std::complex<double>,3,3>::Zero());
   ckm(0,2) = random_phase();
   ckm.bottomLeftCorner<2,2>() = ckmBL;

   Eigen::Matrix<std::complex<double>,3,3> vu = ckm * vd;

   BOOST_REQUIRE(is_unitary(ckm, 1.e-14));
   BOOST_REQUIRE(is_unitary(vd, 1.e-14));
   BOOST_REQUIRE(is_unitary(vu, 1.e-14));
   BOOST_REQUIRE(is_equal(ckm, vu*vd.adjoint(), 1.e-15));

   // remaining mixing matrices
   Eigen::Matrix<std::complex<double>,3,3> uu;
   Eigen::Matrix<std::complex<double>,3,3> ud;

   random_cue_matrix(uu, generator);
   random_cue_matrix(ud, generator);

   BOOST_REQUIRE(is_unitary(uu, 1.e-14));
   BOOST_REQUIRE(is_unitary(ud, 1.e-14));

   // mass eigenvalues
   Eigen::Array<double,3,1> su(Eigen::Array<double,3,1>::Random().abs());
   Eigen::Array<double,3,1> sd(Eigen::Array<double,3,1>::Random().abs());
   std::sort(su.data(), su.data() + su.size());
   std::sort(sd.data(), sd.data() + sd.size());

   // mass matrices
   const Eigen::Matrix<std::complex<double>,3,3> mu(
      uu.transpose() * su.matrix().asDiagonal() * vu);
   const Eigen::Matrix<std::complex<double>,3,3> md(
      ud.transpose() * sd.matrix().asDiagonal() * vd);

   const Eigen::Matrix<double,3,3> ckm_squared_invariants(ckm.cwiseAbs2());

   CKM_parameters::to_pdg_convention(ckm, vu, vd, uu, ud);

   const Eigen::Matrix<double,3,3> ckm_squared_invariants_pdg(ckm.cwiseAbs2());

   BOOST_CHECK(is_unitary(ckm, 1.e-10));

   BOOST_CHECK(is_equal(mu, uu.transpose() * su.matrix().asDiagonal() * vu, 1.e-10));
   BOOST_CHECK(is_equal(md, ud.transpose() * sd.matrix().asDiagonal() * vd, 1.e-10));

   // check all elements are real
   BOOST_CHECK((ckm.imag().cwiseAbs().array() < 1.e-12).all());

   BOOST_CHECK(is_equal(ckm_squared_invariants, ckm_squared_invariants_pdg, 1.e-12));
}

// checks that converting a CKM matrix given in PDG convention leaves
// it unchanged
BOOST_AUTO_TEST_CASE( test_complex_CKM_pdg_convention_consistent )
{
   CKM_parameters ckm_pars;

   ckm_pars.theta_12 = 0.229206;
   ckm_pars.theta_13 = 0.003960;
   ckm_pars.theta_23 = 0.042223;
   ckm_pars.delta    = 0.1;

   const Eigen::Matrix<std::complex<double>,3,3> ckm_cmplx(ckm_pars.get_complex_ckm());
   Eigen::Matrix<std::complex<double>,3,3> ckm(ckm_cmplx);

   const Eigen::Matrix<std::complex<double>,3,3> md(Eigen::Matrix<std::complex<double>,3,3>::Random());
   Eigen::Array<double,3,1> sd;
   Eigen::Matrix<std::complex<double>,3,3> vd;
   Eigen::Matrix<std::complex<double>,3,3> ud;

   fs_svd(md, sd, ud, vd);

   Eigen::Array<double,3,1> su(Eigen::Array<double,3,1>::Random().abs());
   std::sort(su.data(), su.data() + su.size());

   std::mt19937 generator;
   Eigen::Matrix<std::complex<double>,3,3> vu(ckm*vd);
   Eigen::Matrix<std::complex<double>,3,3> uu;
   random_cue_matrix(uu, generator);

   const Eigen::Matrix<std::complex<double>,3,3> mu(
      uu.transpose() * su.matrix().asDiagonal() * vu);

   BOOST_REQUIRE(is_equal(ckm, vu*vd.adjoint(), 1.e-10));

   CKM_parameters::to_pdg_convention(ckm, vu, vd, uu, ud);

   BOOST_CHECK(is_equal(ckm, ckm_cmplx, 1.e-10));
}

#endif
