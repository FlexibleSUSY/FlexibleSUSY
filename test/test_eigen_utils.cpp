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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_eigen_utils

#include <boost/test/unit_test.hpp>
#include <typeinfo>
#include "eigen_utils.hpp"
#include <complex>

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE(test_div_safe_zero)
{
   Eigen::Array<double,3,3> v1, v2, result;
   v1.setZero();
   v2.setZero();
   result = div_safe(v1, v2);

   BOOST_CHECK(result.allFinite());

   for (int i = 0; i < result.rows(); i++)
      for (int k = 0; k < result.cols(); k++)
         BOOST_CHECK_EQUAL(result(i,k), 0.);
}

BOOST_AUTO_TEST_CASE(test_div_safe)
{
   Eigen::Array<double,3,1> v1, v2, result;
   v1 << 1, 2, 3;
   v2 << 0, 1, 2;
   result = div_safe(v1, v2);

   BOOST_CHECK_EQUAL(result(0), 0.);
   BOOST_CHECK_EQUAL(result(1), 2.);
   BOOST_CHECK_EQUAL(result(2), 3./2.);
}

BOOST_AUTO_TEST_CASE(test_div_safe_dynamic)
{
   Eigen::ArrayXd v1(3,1), v2(3,1), result(3,1);
   v1 << 1, 2, 3;
   v2 << 0, 1, 2;
   result = div_safe(v1, v2);

   BOOST_CHECK_EQUAL(result(0), 0.);
   BOOST_CHECK_EQUAL(result(1), 2.);
   BOOST_CHECK_EQUAL(result(2), 3./2.);
}

BOOST_AUTO_TEST_CASE(test_reorder_vector)
{
   Eigen::Array<double,3,1> vec1, vec2;
   vec1 << 1, 2, 3;
   vec2 << 3, 1, 2;

   reorder_vector(vec1, vec2);

   BOOST_CHECK_EQUAL(vec1(0), 3);
   BOOST_CHECK_EQUAL(vec1(1), 1);
   BOOST_CHECK_EQUAL(vec1(2), 2);
}

BOOST_AUTO_TEST_CASE(test_reorder_vector_abs)
{
   Eigen::Array<double,3,1> vec1, vec2;
   vec1 << 1, 2, 3;
   vec2 << 1, 2, -3;

   reorder_vector(vec1, vec2);

   BOOST_CHECK_EQUAL(vec1(0), 2);
   BOOST_CHECK_EQUAL(vec1(1), 3);
   BOOST_CHECK_EQUAL(vec1(2), 1);
}

BOOST_AUTO_TEST_CASE(test_reorder_vector_with_matrix)
{
   Eigen::Array<double,3,1> vec1;
   vec1 << 1, 2, 3;

   Eigen::Matrix<double,3,3> matrix;
   matrix << 3, 0, 0,
             0, 1, 0,
             0, 0, 2;

   reorder_vector(vec1, matrix);

   BOOST_CHECK_EQUAL(vec1(0), 3);
   BOOST_CHECK_EQUAL(vec1(1), 1);
   BOOST_CHECK_EQUAL(vec1(2), 2);
}

BOOST_AUTO_TEST_CASE(test_remove_if_equal)
{
   Eigen::Array<double,5,1> src;
   src << 1, 2, 3, 4, 5;

   Eigen::Array<double,2,1> cmp; // to be removed from src
   cmp << 2, 4;

   const Eigen::Array<double,3,1> dst = remove_if_equal(src, cmp);

   BOOST_CHECK_EQUAL(dst(0), 1);
   BOOST_CHECK_EQUAL(dst(1), 3);
   BOOST_CHECK_EQUAL(dst(2), 5);
}

BOOST_AUTO_TEST_CASE(test_remove_if_equal_double_indices)
{
   Eigen::Array<double,5,1> src;
   src << 1, 2, 2.3, 4, 5;

   Eigen::Array<double,2,1> cmp; // to be removed from src
   cmp << 2, 2.1;

   const Eigen::Array<double,3,1> dst = remove_if_equal(src, cmp);

   BOOST_CHECK_EQUAL(dst(0), 1);
   BOOST_CHECK_EQUAL(dst(1), 4);
   BOOST_CHECK_EQUAL(dst(2), 5);
}

BOOST_AUTO_TEST_CASE(test_normalize_to_interval_real)
{
   Eigen::Matrix<double,2,2> m;
   m << -2., -1., 0., 2;

   normalize_to_interval(m, -1., 1.);

   BOOST_CHECK_EQUAL(m(0,0), -1.);
   BOOST_CHECK_EQUAL(m(0,1), -1.);
   BOOST_CHECK_EQUAL(m(1,0),  0.);
   BOOST_CHECK_EQUAL(m(1,1),  1.);
}

BOOST_AUTO_TEST_CASE(test_normalize_to_interval_complex)
{
   Eigen::Matrix<std::complex<double>,2,2> m;
   m << -2., -1., 0., 2.;

   normalize_to_interval(m, 1.);

   BOOST_CHECK_CLOSE(std::real(m(0,0)), -1., 1e-15);
   BOOST_CHECK_CLOSE(std::real(m(0,1)), -1., 1e-15);
   BOOST_CHECK_SMALL(std::real(m(1,0)), 1e-15);
   BOOST_CHECK_CLOSE(std::real(m(1,1)),  1., 1e-15);
   BOOST_CHECK_SMALL(std::imag(m(0,0)), 1e-15);
   BOOST_CHECK_SMALL(std::imag(m(0,1)), 1e-15);
   BOOST_CHECK_SMALL(std::imag(m(1,0)), 1e-15);
   BOOST_CHECK_SMALL(std::imag(m(1,1)), 1e-15);

   m << std::polar(2.,1.), std::polar(1.,1.),
        std::polar(0.5,1.), std::polar(2.,-1.);

   normalize_to_interval(m, 1.);

   BOOST_CHECK_CLOSE(std::abs(m(0,0)),  1., 1e-15);
   BOOST_CHECK_CLOSE(std::abs(m(0,1)),  1., 1e-15);
   BOOST_CHECK_CLOSE(std::abs(m(1,0)), 0.5, 1e-15);
   BOOST_CHECK_CLOSE(std::abs(m(1,1)),  1., 1e-15);
   BOOST_CHECK_CLOSE(std::arg(m(0,0)),  1., 1e-15);
   BOOST_CHECK_CLOSE(std::arg(m(0,1)),  1., 1e-15);
   BOOST_CHECK_CLOSE(std::arg(m(1,0)),  1., 1e-15);
   BOOST_CHECK_CLOSE(std::arg(m(1,1)), -1., 1e-15);
}
