#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_max_abs_re_eigenvalue

#include <boost/test/unit_test.hpp>

#include "unitarity.hpp"

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

using namespace flexiblesusy;

namespace {

double from_all_eigenvalues(Eigen::MatrixXcd const& m)
{
   return m.eigenvalues().real().cwiseAbs().maxCoeff();
}

/// random blocks of the given sizes, with rows and columns shuffled
Eigen::MatrixXcd make_shuffled_blocks(std::vector<int> const& sizes, bool symmetric, unsigned seed)
{
   std::mt19937 gen(seed);
   std::uniform_real_distribution<double> dist(-1., 1.);

   const int n = std::accumulate(sizes.begin(), sizes.end(), 0);
   std::vector<int> position(n);
   std::iota(position.begin(), position.end(), 0);
   std::shuffle(position.begin(), position.end(), gen);

   Eigen::MatrixXcd m = Eigen::MatrixXcd::Zero(n, n);

   int offset = 0;
   for (const int size: sizes) {
      for (int i = 0; i < size; ++i) {
         for (int j = symmetric ? i : 0; j < size; ++j) {
            const std::complex<double> el(dist(gen), dist(gen));
            m(position[offset + i], position[offset + j]) = el;
            if (symmetric) {
               m(position[offset + j], position[offset + i]) = el;
            }
         }
      }
      offset += size;
   }

   return m;
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(test_diagonal_matrix_uses_real_parts)
{
   Eigen::MatrixXcd m = Eigen::MatrixXcd::Zero(3, 3);
   m(0, 0) = std::complex<double>(0.1, 5.);
   m(1, 1) = std::complex<double>(-0.3, 0.);
   m(2, 2) = std::complex<double>(0.2, -7.);

   BOOST_CHECK_EQUAL(max_abs_re_eigenvalue(m), 0.3);
}

BOOST_AUTO_TEST_CASE(test_zero_matrix)
{
   BOOST_CHECK_EQUAL(max_abs_re_eigenvalue(Eigen::MatrixXcd::Zero(4, 4)), 0.);
}

BOOST_AUTO_TEST_CASE(test_blocks_found_through_one_triangle_only)
{
   // eigenvalues {2, 3} and -4, where 0 and 2 are connected by a single element
   Eigen::MatrixXcd m = Eigen::MatrixXcd::Zero(3, 3);
   m(0, 0) = 2.;
   m(0, 2) = 1.;
   m(2, 2) = 3.;
   m(1, 1) = -4.;

   BOOST_CHECK_CLOSE_FRACTION(max_abs_re_eigenvalue(m), 4., 1e-14);

   m(1, 1) = 1.;

   BOOST_CHECK_CLOSE_FRACTION(max_abs_re_eigenvalue(m), 3., 1e-14);
}

BOOST_AUTO_TEST_CASE(test_same_result_as_whole_matrix)
{
   const std::vector<int> sizes{1, 7, 2, 1, 12, 3, 5, 1};

   for (unsigned seed = 0; seed < 20; ++seed) {
      for (const bool symmetric: {true, false}) {
         const auto m = make_shuffled_blocks(sizes, symmetric, seed);

         BOOST_CHECK_CLOSE_FRACTION(max_abs_re_eigenvalue(m), from_all_eigenvalues(m), 1e-11);
      }
   }
}

BOOST_AUTO_TEST_CASE(test_single_block)
{
   const auto m = make_shuffled_blocks({30}, false, 42);

   BOOST_CHECK_CLOSE_FRACTION(max_abs_re_eigenvalue(m), from_all_eigenvalues(m), 1e-11);
}

BOOST_AUTO_TEST_CASE(test_non_finite_elements_give_nan)
{
   Eigen::MatrixXcd m = Eigen::MatrixXcd::Identity(3, 3);
   m(0, 1) = std::numeric_limits<double>::quiet_NaN();

   BOOST_CHECK(std::isnan(max_abs_re_eigenvalue(m)));

   m(0, 1) = std::numeric_limits<double>::infinity();

   BOOST_CHECK(std::isnan(max_abs_re_eigenvalue(m)));
}
