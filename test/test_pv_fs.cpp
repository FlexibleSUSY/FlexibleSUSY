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
#define BOOST_TEST_MODULE test_pv_fs

#include <boost/test/unit_test.hpp>

#include "pv_fs.hpp"
#include "read_data.hpp"
#include <algorithm>
#include <iomanip>
#include <limits>
#include <vector>


bool is_close(double a, double b, double eps)
{
   const double max = std::max(std::abs(a), std::abs(b));
   return std::abs(a - b) <= eps*max;
}


double rel_dev(double a, double b)
{
   const double max = std::max(std::abs(a), std::abs(b));
   return std::abs(a - b)/max;
}


using namespace flexiblesusy;

struct A_data {
   double m2{}, q2{}, a0{};
};

struct B_data {
   double p2{}, m12{}, m22{}, q2{}, b0{};
};

struct DB_data {
   double m12{}, m22{}, db0{};
};

std::ostream& operator<<(std::ostream& ostr, const A_data& a0)
{
   ostr << std::setprecision(std::numeric_limits<double>::digits10)
        << "A0(m2=" << a0.m2 << ", q2=" << a0.q2 << ") = " << a0.a0;
   return ostr;
}

std::ostream& operator<<(std::ostream& ostr, const B_data& b0)
{
   ostr << std::setprecision(std::numeric_limits<double>::digits10)
        << "B(p2=" << b0.p2 << ", m12=" << b0.m12
        << ", m22=" << b0.m22 << ", q2=" << b0.q2 << ") = " << b0.b0;
   return ostr;
}

std::ostream& operator<<(std::ostream& ostr, const DB_data& db0)
{
   ostr << std::setprecision(std::numeric_limits<double>::digits10)
        << "DB0(m12=" << db0.m12 << ", m22=" << db0.m22 << ") = " << db0.db0;
   return ostr;
}

/// tranform vector of elements of type A -> B
template <class B, class A, class F>
std::vector<B> fmap(F f, const std::vector<A>& in)
{
   std::vector<B> out;
   std::transform(in.cbegin(), in.cend(), std::back_inserter(out), f);
   return out;
}

/// read A function data
std::vector<A_data> read_a0(const std::string& filename)
{
   return fmap<A_data>(
      [](const auto& d) {
         return A_data{d.at(0), d.at(1), d.at(2)};
      },
      test::read_from_file<double>(filename));
}

/// read B function data
std::vector<B_data> read_b0(const std::string& filename)
{
   return fmap<B_data>(
      [](const auto& d) {
         return B_data{d.at(0), d.at(1), d.at(2), d.at(3), d.at(4)};
      },
      test::read_from_file<double>(filename));
}

/// read DB function data
std::vector<DB_data> read_db0(const std::string& filename)
{
   return fmap<DB_data>(
      [](const auto& d) {
         return DB_data{d.at(0), d.at(1), d.at(2)};
      },
      test::read_from_file<double>(filename));
}

BOOST_AUTO_TEST_CASE( test_ReA0_values )
{
   const auto filename = std::string(TEST_DATA_DIR) + test::PATH_SEPARATOR + "A0.dat";
   const auto data = read_a0(filename);
   const double eps = 1e-12;

   for (auto d: data) {
      const auto a0 = flexiblesusy::a0(d.m2, d.q2);
      BOOST_CHECK_MESSAGE(is_close(d.a0, a0, eps), "expected: " << d << ", observed: " << a0);
   }
}

BOOST_AUTO_TEST_CASE( test_ReB0xx_values )
{
   const auto filename = std::string(TEST_DATA_DIR) + test::PATH_SEPARATOR + "B0xx.dat";
   const auto data = read_b0(filename);
   const double eps = 1e-12;

   for (auto d: data) {
      const auto b0 = flexiblesusy::b0xx(d.p2, d.m12, d.q2);
      BOOST_CHECK_MESSAGE(is_close(d.b0, b0, eps), "expected: " << d << ", observed: " << b0);
   }
}

BOOST_AUTO_TEST_CASE( test_ReB0_values )
{
   const auto filename = std::string(TEST_DATA_DIR) + test::PATH_SEPARATOR + "B0.dat";
   const auto data = read_b0(filename);
   const double eps = 1e-12;

   for (auto d: data) {
      const auto b0 = flexiblesusy::b0(d.p2, d.m12, d.m22, d.q2);
      BOOST_CHECK_MESSAGE(is_close(d.b0, b0, eps), "expected: " << d << ", observed: " << b0);
   }
}

BOOST_AUTO_TEST_CASE( test_ReD1B0_values )
{
   const auto filename = std::string(TEST_DATA_DIR) + test::PATH_SEPARATOR + "DB0.dat";
   const auto data = read_db0(filename);
   const double eps = 1e-12;

   for (auto d: data) {
      const auto db0 = flexiblesusy::d1_b0(d.m12, d.m22);
      BOOST_CHECK_MESSAGE(is_close(d.db0, db0, eps), "expected: " << d << ", observed: " << db0);
   }
}

BOOST_AUTO_TEST_CASE( test_ReB22_values )
{
   const auto filename = std::string(TEST_DATA_DIR) + test::PATH_SEPARATOR + "B22.dat";
   const auto data = read_b0(filename);
   const double eps = 1e-12;

   for (auto d: data) {
      const auto b22 = flexiblesusy::b22(d.p2, d.m12, d.m22, d.q2);
      BOOST_CHECK_MESSAGE(is_close(d.b0, b22, eps), "expected: " << d << ", observed: " << b22
                          << ", rel. deviation: " << rel_dev(d.b0, b22));
   }
}
