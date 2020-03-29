
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_pmns

#include <boost/test/unit_test.hpp>

#include "test.hpp"
#include "linalg2.hpp"
#include "pmns.hpp"
#include "random_matrix.hpp"
#include "wrappers.hpp"

using namespace flexiblesusy;

bool is_orthogonal(const Eigen::Matrix<double,3,3>& m,
                   double max_deviation = std::numeric_limits<double>::epsilon())
{
   const Eigen::Matrix<double,3,3> unity(m * m.transpose());

   bool result = is_equal(unity, Eigen::Matrix<double,3,3>::Identity(),
                          max_deviation);

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

bool is_symmetric(const Eigen::Matrix<double,3,3>& m,
                  double max_deviation = std::numeric_limits<double>::epsilon())
{
   bool result = is_equal(m, m.transpose(), max_deviation);

   if (!result)
      BOOST_TEST_MESSAGE("matrix is not symmetric!"
                         "\n   Matrix = " << m <<
                         "\n   Matrix^T = " << m.transpose() <<
                         "\n   Maximum allowed max_deviation from equality = " << max_deviation);

   return result;
}

bool is_symmetric(const Eigen::Matrix<std::complex<double>,3,3>& m,
                  double max_deviation = std::numeric_limits<double>::epsilon())
{
   bool result = is_equal(m, m.transpose(), max_deviation);

   if (!result)
      BOOST_TEST_MESSAGE("matrix is not symmetric!"
                         "\n   Matrix = " << m <<
                         "\n   Matrix^T = " << m.transpose() <<
                         "\n   Maximum allowed max_deviation from equality = " << max_deviation);

   return result;
}

double random_angle()
{
   return 2. * M_PI * rand() / RAND_MAX;
}

BOOST_AUTO_TEST_CASE( test_PMNS_unitarity_from_angles )
{
   PMNS_parameters pmns_pars;

   pmns_pars.theta_12 = 0.5764;
   pmns_pars.theta_13 = 0.1472;
   pmns_pars.theta_23 = 0.7101;
   pmns_pars.delta = 4.3354;

   pmns_pars.alpha_1 = 0.;
   pmns_pars.alpha_2 = 0.;

   const Eigen::Matrix<double,3,3> pmns_real_no_alpha(
      pmns_pars.get_real_pmns());
   BOOST_CHECK(is_unitary(pmns_real_no_alpha));

   const Eigen::Matrix<std::complex<double>,3,3> pmns_complex_no_alpha(
      pmns_pars.get_complex_pmns());
   BOOST_CHECK(is_unitary(pmns_complex_no_alpha));

   pmns_pars.alpha_1 = random_angle();
   pmns_pars.alpha_2 = random_angle();

   const Eigen::Matrix<double,3,3> pmns_real(pmns_pars.get_real_pmns());
   BOOST_CHECK(is_unitary(pmns_real));

   const Eigen::Matrix<std::complex<double>,3,3> pmns_complex(
      pmns_pars.get_complex_pmns());
   BOOST_CHECK(is_unitary(pmns_complex));
}

BOOST_AUTO_TEST_CASE( test_real_PMNS_pdg_convention )
{
   // fermion mass matrices
   const Eigen::Matrix<double,3,3> me(Eigen::Matrix<double,3,3>::Random());
   Eigen::Matrix<double,3,3> mv(Eigen::Matrix<double,3,3>::Random());

   Symmetrize(mv);

   BOOST_REQUIRE(is_symmetric(mv));

   // mass eigenvalues
   Eigen::Array<double,3,1> se;
   Eigen::Array<double,3,1> sv;

   // mixing matrices
   Eigen::Matrix<double,3,3> ve;
   Eigen::Matrix<double,3,3> ue;
   Eigen::Matrix<double,3,3> vv;

   fs_svd(me, se, ue, ve);
   fs_diagonalize_hermitian(mv, sv, vv);

   BOOST_CHECK(is_equal(me, ue.transpose() * se.matrix().asDiagonal() * ve,
                        1.e-10));
   BOOST_CHECK(is_equal(mv, vv.transpose() * sv.matrix().asDiagonal() * vv,
                        1.e-10));

   Eigen::Matrix<double,3,3> pmns(ve*vv.adjoint());

   // transpose in order to make (1, 2) negative
   pmns.transposeInPlace();
   BOOST_CHECK(is_unitary(pmns, 1.e-10));

   PMNS_parameters::to_pdg_convention(pmns, vv, ve, ue);

   BOOST_CHECK(pmns(1,2) >= 0.);
   BOOST_CHECK(pmns(2,2) >= 0.);

   BOOST_CHECK(is_equal(me, ue.transpose() * se.matrix().asDiagonal() * ve,
                        1.e-10));
   BOOST_CHECK(is_equal(mv, vv.transpose() * sv.matrix().asDiagonal() * vv,
                        1.e-10));

   {
      // check that converted matrices are consistent
      const Eigen::Matrix<double,3,3> pmns_check(ve*vv.adjoint());
      for (int i = 0; i < 3; i++) {
         for (int k = 0; k < 3; k++) {
            // check for transposed equality here, because the PMNS
            // matrix from above was transposed to generate a negative
            // (1,2) element
            BOOST_CHECK_CLOSE_FRACTION(Re(pmns(i,k)), Re(pmns_check(k,i)),
                                       1.0e-10);
            BOOST_CHECK_CLOSE_FRACTION(Im(pmns(i,k)), Im(pmns_check(k,i)),
                                       1.0e-10);
         }
      }
   }
}

BOOST_AUTO_TEST_CASE( test_complex_PMNS_pdg_convention )
{
   // fermion mass matrices
   const Eigen::Matrix<std::complex<double>,3,3> me(
      Eigen::Matrix<std::complex<double>,3,3>::Random());
   Eigen::Matrix<std::complex<double>,3,3> mv(
      Eigen::Matrix<std::complex<double>,3,3>::Random());

   Symmetrize(mv);

   BOOST_REQUIRE(is_symmetric(mv));

   // mass eigenvalues
   Eigen::Array<double,3,1> se;
   Eigen::Array<double,3,1> sv;

   // mixing matrices
   Eigen::Matrix<std::complex<double>,3,3> ve;
   Eigen::Matrix<std::complex<double>,3,3> ue;
   Eigen::Matrix<std::complex<double>,3,3> vv;

   fs_svd(me, se, ue, ve);
   fs_diagonalize_symmetric(mv, sv, vv);

   BOOST_CHECK(is_equal(me, ue.transpose() * se.matrix().asDiagonal() * ve,
                        1.e-10));
   BOOST_CHECK(is_equal(mv, vv.transpose() * sv.matrix().asDiagonal() * vv,
                        1.e-10));

   Eigen::Matrix<std::complex<double>,3,3> pmns(ve*vv.adjoint());

   PMNS_parameters::to_pdg_convention(pmns, vv, ve, ue);

   BOOST_CHECK(is_unitary(pmns, 1.e-10));

   BOOST_CHECK(is_equal(me, ue.transpose() * se.matrix().asDiagonal() * ve,
                        1.e-10));
   BOOST_CHECK(is_equal(mv, vv.transpose() * sv.matrix().asDiagonal() * vv,
                        1.e-10));

   BOOST_CHECK_LT(std::abs(std::arg(pmns(1,2))), 1e-15);
   BOOST_CHECK_LT(std::abs(std::arg(pmns(2,2))), 1e-15);

   const auto a1 = std::polar(1., std::arg(pmns(0,0)));
   const auto a2 = std::polar(1., std::arg(pmns(0,1)));

   const Eigen::Matrix<std::complex<double>,2,2> pmnsBL(
      pmns.bottomLeftCorner<2,2>());
   Eigen::Array<std::complex<double>,2,2> maj_phases_expected;
   maj_phases_expected << a1, a2, a1, a2;
   const Eigen::Array<std::complex<double>,2,2> no_maj_phases
      = maj_phases_expected.conjugate() * pmnsBL.array();
   const Eigen::Array<double,2,2> imagBL(no_maj_phases.imag());

   BOOST_CHECK((imagBL <= 0).all() || (imagBL >= 0).all());
}

#ifdef ENABLE_RANDOM

std::complex<double> random_phase()
{
   return std::polar(1., 2. * M_PI * rand() / RAND_MAX);
}

BOOST_AUTO_TEST_CASE( test_complex_PMNS_pdg_convention_zero_c13 )
{
   Eigen::Matrix<std::complex<double>,2,2> pmnsBL;
   Eigen::Matrix<std::complex<double>,3,3> vv;

   std::mt19937 generator;
   random_cue_matrix(pmnsBL, generator);
   random_cue_matrix(vv, generator);

   Eigen::Matrix<std::complex<double>,3,3> pmns(
      Eigen::Matrix<std::complex<double>,3,3>::Zero());
   pmns(0,2) = random_phase();
   pmns.bottomLeftCorner<2,2>() = pmnsBL;

   Eigen::Matrix<std::complex<double>,3,3> ve = pmns * vv;

   BOOST_REQUIRE(is_unitary(pmns, 1.e-15));
   BOOST_REQUIRE(is_unitary(vv, 1.e-15));
   BOOST_REQUIRE(is_unitary(ve, 1.e-15));
   BOOST_REQUIRE(is_equal(pmns, ve*vv.adjoint(), 1.e-15));

   // remaining mixing matrices
   Eigen::Matrix<std::complex<double>,3,3> ue;
   random_cue_matrix(ue, generator);

   BOOST_REQUIRE(is_unitary(ue, 1.e-15));

   // mass eigenvalues
   Eigen::Array<double,3,1> se(Eigen::Array<double,3,1>::Random().abs());
   Eigen::Array<double,3,1> sv(Eigen::Array<double,3,1>::Random().abs());
   std::sort(se.data(), se.data() + se.size());
   std::sort(sv.data(), sv.data() + sv.size());

   // mass matrices
   const Eigen::Matrix<std::complex<double>,3,3> me(
      ue.transpose() * se.matrix().asDiagonal() * ve);
   const Eigen::Matrix<std::complex<double>,3,3> mv(
      vv.transpose() * sv.matrix().asDiagonal() * vv);

   BOOST_REQUIRE(is_symmetric(mv, 1.e-15));

   const Eigen::Matrix<double,3,3> pmns_squared_invariants(pmns.cwiseAbs2());

   PMNS_parameters::to_pdg_convention(pmns, vv, ve, ue);

   const Eigen::Matrix<double,3,3> pmns_squared_invariants_pdg(
      pmns.cwiseAbs2());

   BOOST_CHECK(is_unitary(pmns, 1.e-10));

   BOOST_CHECK(is_equal(me, ue.transpose() * se.matrix().asDiagonal() * ve,
                        1.e-10));
   BOOST_CHECK(is_equal(mv, vv.transpose() * sv.matrix().asDiagonal() * vv,
                        1.e-10));

   // check CP-violating phase is eliminated
   BOOST_CHECK_LT(std::abs(std::arg(pmns(0,2))), 1.e-15);

   // check remaining Majorana phases are consistent (note: relative sign
   // between (1,0) and (2,0) elements)
   BOOST_CHECK(is_equal(std::arg(pmns(1,0)), std::arg(-pmns(2,0)), 1.e-10));
   BOOST_CHECK(is_equal(std::arg(pmns(1,1)), std::arg(pmns(2,1)), 1.e-10));

   BOOST_CHECK(is_equal(pmns_squared_invariants, pmns_squared_invariants_pdg,
                        1.e-12));
}

// checks that converting to a PMNS matrix given in PDG convention leaves
// it unchanged
BOOST_AUTO_TEST_CASE( test_complex_PMNS_pdg_convention_consistent_no_alpha )
{
   PMNS_parameters pmns_pars;

   pmns_pars.theta_12 = 0.5764;
   pmns_pars.theta_13 = 0.1472;
   pmns_pars.theta_23 = 0.7101;
   pmns_pars.delta = 4.3354;

   pmns_pars.alpha_1 = 0.;
   pmns_pars.alpha_2 = 0.;

   const Eigen::Matrix<std::complex<double>,3,3> pmns_complex_no_alpha(
      pmns_pars.get_complex_pmns());
   Eigen::Matrix<std::complex<double>,3,3> pmns(pmns_complex_no_alpha);

   Eigen::Matrix<std::complex<double>,3,3> mv(
      Eigen::Matrix<std::complex<double>,3,3>::Random());
   Symmetrize(mv);

   BOOST_REQUIRE(is_symmetric(mv, 1.e-15));

   Eigen::Array<double,3,1> sv;
   Eigen::Matrix<std::complex<double>,3,3> vv;

   fs_diagonalize_symmetric(mv, sv, vv);

   Eigen::Array<double,3,1> se(Eigen::Array<double,3,1>::Random().abs());
   std::sort(se.data(), se.data() + se.size());

   std::mt19937 generator;
   Eigen::Matrix<std::complex<double>,3,3> ve(pmns*vv);
   Eigen::Matrix<std::complex<double>,3,3> ue;
   random_cue_matrix(ue, generator);

   const Eigen::Matrix<std::complex<double>,3,3> me(
      ue.transpose() * se.matrix().asDiagonal() * ve);

   BOOST_REQUIRE(is_equal(pmns, ve*vv.adjoint(), 1.e-10));

   PMNS_parameters::to_pdg_convention(pmns, vv, ve, ue);

   BOOST_CHECK(is_equal(pmns, pmns_complex_no_alpha, 1.e-10));
}

BOOST_AUTO_TEST_CASE( test_complex_PMNS_pdg_convention_consistent_no_delta )
{
   PMNS_parameters pmns_pars;

   pmns_pars.theta_12 = 0.5764;
   pmns_pars.theta_13 = 0.1472;
   pmns_pars.theta_23 = 0.7101;
   pmns_pars.delta = 0.;

   pmns_pars.alpha_1 = random_angle();
   pmns_pars.alpha_2 = random_angle();

   const Eigen::Matrix<std::complex<double>,3,3> pmns_complex_no_delta(
      pmns_pars.get_complex_pmns());
   Eigen::Matrix<std::complex<double>,3,3> pmns(pmns_complex_no_delta);

   Eigen::Matrix<std::complex<double>,3,3> mv(
      Eigen::Matrix<std::complex<double>,3,3>::Random());
   Symmetrize(mv);

   BOOST_REQUIRE(is_symmetric(mv, 1.e-15));

   Eigen::Array<double,3,1> sv;
   Eigen::Matrix<std::complex<double>,3,3> vv;

   fs_diagonalize_symmetric(mv, sv, vv);

   Eigen::Array<double,3,1> se(Eigen::Array<double,3,1>::Random().abs());
   std::sort(se.data(), se.data() + se.size());

   std::mt19937 generator;
   Eigen::Matrix<std::complex<double>,3,3> ve(pmns*vv);
   Eigen::Matrix<std::complex<double>,3,3> ue;
   random_cue_matrix(ue, generator);

   const Eigen::Matrix<std::complex<double>,3,3> me(
      ue.transpose() * se.matrix().asDiagonal() * ve);

   BOOST_REQUIRE(is_equal(pmns, ve*vv.adjoint(), 1.e-10));

   PMNS_parameters::to_pdg_convention(pmns, vv, ve, ue);

   BOOST_CHECK(is_equal(pmns, pmns_complex_no_delta, 1.e-10));

   const auto a1 = std::polar(1., std::arg(pmns_complex_no_delta(0,0)));
   const auto a2 = std::polar(1., std::arg(pmns_complex_no_delta(0,1)));

   BOOST_CHECK_LT(std::abs(std::imag(std::conj(a1) * pmns(1,0))), 1.e-10);
   BOOST_CHECK_LT(std::abs(std::imag(std::conj(a2) * pmns(1,1))), 1.e-10);
   BOOST_CHECK_LT(std::abs(std::imag(std::conj(a1) * pmns(2,0))), 1.e-10);
   BOOST_CHECK_LT(std::abs(std::imag(std::conj(a2) * pmns(2,1))), 1.e-10);
}

BOOST_AUTO_TEST_CASE( test_complex_PMNS_pdg_convention_consistent )
{
   PMNS_parameters pmns_pars;

   pmns_pars.theta_12 = 0.5764;
   pmns_pars.theta_13 = 0.1472;
   pmns_pars.theta_23 = 0.7101;
   pmns_pars.delta = random_angle();

   pmns_pars.alpha_1 = random_angle();
   pmns_pars.alpha_2 = random_angle();

   const Eigen::Matrix<std::complex<double>,3,3> pmns_complex(
      pmns_pars.get_complex_pmns());
   Eigen::Matrix<std::complex<double>,3,3> pmns(pmns_complex);

   Eigen::Matrix<std::complex<double>,3,3> mv(
      Eigen::Matrix<std::complex<double>,3,3>::Random());
   Symmetrize(mv);

   BOOST_REQUIRE(is_symmetric(mv, 1.e-15));

   Eigen::Array<double,3,1> sv;
   Eigen::Matrix<std::complex<double>,3,3> vv;

   fs_diagonalize_symmetric(mv, sv, vv);

   Eigen::Array<double,3,1> se(Eigen::Array<double,3,1>::Random().abs());
   std::sort(se.data(), se.data() + se.size());

   std::mt19937 generator;
   Eigen::Matrix<std::complex<double>,3,3> ve(pmns*vv);
   Eigen::Matrix<std::complex<double>,3,3> ue;
   random_cue_matrix(ue, generator);

   const Eigen::Matrix<std::complex<double>,3,3> me(
      ue.transpose() * se.matrix().asDiagonal() * ve);

   BOOST_REQUIRE(is_equal(pmns, ve*vv.adjoint(), 1.e-10));

   PMNS_parameters::to_pdg_convention(pmns, vv, ve, ue);

   BOOST_CHECK(is_equal(pmns, pmns_complex, 1.e-10));
}

#endif
