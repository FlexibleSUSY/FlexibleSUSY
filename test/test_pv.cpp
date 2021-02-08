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

#include <cmath>
#include <iomanip>
#include <limits>
#include <vector>
#include "derivative.hpp"
#include "numerics.h"
#include "pv.hpp"
#include "read_data.hpp"
#include "stopwatch.hpp"
#include "wrappers.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_pv

#include <boost/test/unit_test.hpp>
#include "numerics2.hpp"
#include "rk.hpp"
#include <Eigen/Dense>

using namespace std;
using namespace flexiblesusy;
using namespace flexiblesusy::passarino_veltman;

struct B0_args {
   double p{}, m1{}, m2{}, q{};
};

constexpr double TOL = 1e-4;
constexpr double sqr(double a) noexcept { return a*a; }

double refnfn(double x, double p, double m1, double m2, double q) noexcept
{
  using flexiblesusy::fast_log;
  const static std::complex<double> iEpsilon(0.0, TOL * 1.0e-20);

  return std::real(x *
    fast_log(((1 - x) * sqr(m1) + x * sqr(m2)
              - x * (1 - x) * sqr(p) - iEpsilon) / sqr(q)));
}

double integrandThreshbnr(double x, double p, double m1, double m2, double q) noexcept
{
  return refnfn(x, p, m1, m2, q);
}

Eigen::Array<double,1,1> dd(double x, double p, double m1, double m2, double q) noexcept
{
  Eigen::Array<double,1,1> dydx;
  dydx(0) = -integrandThreshbnr(x, p, m1, m2, q);
  return dydx;
}

// Returns real part of integral
double bIntegral(double p, double m1, double m2, double q) noexcept
{
  using namespace flexiblesusy;

  const double from = 0.0, to = 1.0, guess = 0.1, hmin = TOL * 1.0e-5;
  const double eps = TOL * 1.0e-3;
  Eigen::Array<double,1,1> v;
  v(0) = 1.0;

  try {
     runge_kutta::integrateOdes(v, from, to, eps, guess, hmin,
                                [p, m1, m2, q] (double x, const Eigen::Array<double,1,1>&) {
                                   return dd(x, p, m1, m2, q);
                                });
  } catch (...) {
     ERROR("B1 integral did not converge.");
     v(0) = 1.;
  }

  return v(0) - 1.0;
}

BOOST_AUTO_TEST_CASE(test_real_parts)
{
    BOOST_CHECK_CLOSE_FRACTION(ReA0 (2,     1),  0.61370563888010943, 1e-14);
    BOOST_CHECK_CLOSE_FRACTION(ReB0 (2,3,4, 1), -1.14772143349266490, 1e-14);
    BOOST_CHECK_CLOSE_FRACTION(ReB1 (2,3,4, 1),  0.59926550299197479, 1e-14);
    BOOST_CHECK_CLOSE_FRACTION(ReB00(2,3,4, 1), -0.24981786818504056, 1e-14);
}

const double scale  = 100;
const double scale2 = Sqr(scale);
const double p  = 91.0;
const double p2 = Sqr(p);

BOOST_AUTO_TEST_CASE( test_ReA0 )
{
   BOOST_CHECK_EQUAL(ReA0(0., scale2), 0.);
}

BOOST_AUTO_TEST_CASE( test_ReB0 )
{
   BOOST_CHECK_EQUAL(ReB0(0., 0., 0., scale2), 0.);

   BOOST_CHECK_CLOSE(ReB0(p2, 0., 0., scale2), 2. - log(p2/scale2), 1.0e-10);
   BOOST_CHECK_CLOSE(ReB0(0., p2, 0., scale2), 1. - log(p2/scale2), 1.0e-10);
   BOOST_CHECK_EQUAL(ReB0(0., 0., p2, scale2), ReB0(0., p2, 0., scale2));

   BOOST_CHECK_CLOSE(ReB0(p2, p2, 0., scale2), 2. - log(p2/scale2), 0.005);
   BOOST_CHECK_CLOSE(ReB0(0., p2, p2, scale2), 0. - log(p2/scale2), 0.005);
   BOOST_CHECK_EQUAL(ReB0(p2, 0., p2, scale2), ReB0(p2, p2, 0., scale2));

   BOOST_CHECK_CLOSE(ReB0(1., 1., 2., 1.), -0.263943507355163, 1.0e-9);
}

struct A0_data {
   double m2{}, q2{}, a0{};
};

struct B0_data {
   double p2{}, m12{}, m22{}, q2{}, b0{};
};

struct DB0_data {
   double m12{}, m22{}, db0{};
};

std::ostream& operator<<(std::ostream& ostr, const A0_data& a0)
{
   ostr << std::setprecision(std::numeric_limits<double>::digits10)
        << "A0(m2=" << a0.m2 << ", q2=" << a0.q2 << ") = " << a0.a0;
   return ostr;
}

std::ostream& operator<<(std::ostream& ostr, const B0_data& b0)
{
   ostr << std::setprecision(std::numeric_limits<double>::digits10)
        << "B0(p2=" << b0.p2 << ", m12=" << b0.m12
        << ", m22=" << b0.m22 << ", q2=" << b0.q2 << ") = " << b0.b0;
   return ostr;
}

std::ostream& operator<<(std::ostream& ostr, const DB0_data& db0)
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

/// read A0 function data
std::vector<A0_data> read_a0(const std::string& filename)
{
   return fmap<A0_data>(
      [](const auto& d) {
         return A0_data{d.at(0), d.at(1), d.at(2)};
      },
      test::read_from_file<double>(filename));
}

/// read B0 function data
std::vector<B0_data> read_b0(const std::string& filename)
{
   return fmap<B0_data>(
      [](const auto& d) {
         return B0_data{d.at(0), d.at(1), d.at(2), d.at(3), d.at(4)};
      },
      test::read_from_file<double>(filename));
}

/// read DB0 function data
std::vector<DB0_data> read_db0(const std::string& filename)
{
   return fmap<DB0_data>(
      [](const auto& d) {
         return DB0_data{d.at(0), d.at(1), d.at(2)};
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
      BOOST_TEST_MESSAGE("expected: " << d);
      BOOST_TEST_MESSAGE("observed: " << a0);
      BOOST_CHECK_CLOSE_FRACTION(d.a0, a0, eps);
   }
}

BOOST_AUTO_TEST_CASE( test_ReD1B0_values )
{
   const auto filename = std::string(TEST_DATA_DIR) + test::PATH_SEPARATOR + "DB0.dat";
   const auto data = read_db0(filename);
   const double eps = 1e-12;

   for (auto d: data) {
      const auto db0 = flexiblesusy::d1_b0(d.m12, d.m22);
      BOOST_TEST_MESSAGE("expected: " << d);
      BOOST_TEST_MESSAGE("observed: " << db0);
      BOOST_CHECK_CLOSE_FRACTION(d.db0, db0, eps);
   }
}

BOOST_AUTO_TEST_CASE( test_ReD1B0 )
{
   BOOST_CHECK_EQUAL(ReD1B0(0., 0., 0.), 0.);

   BOOST_CHECK_CLOSE_FRACTION(ReD1B0(0., p2, 0.), 1./(2.*p2), 2e-6);
   BOOST_CHECK_EQUAL(ReD1B0(0., 0., p2), ReD1B0(0., p2, 0.));
   BOOST_CHECK_CLOSE_FRACTION(ReD1B0(p2, 0., p2), ReD1B0(p2, p2, 0.), 1e-12);
}

BOOST_AUTO_TEST_CASE( test_ReD1B0_numerical )
{
   const std::vector<B0_args> vals = {
      {1e-2, 1., 0., 1. },
      {1e-2, 1., 1., 1. },
      {1e-2, 1., 1., 10.},
      {1e-2, 2., 1., 1. },
      {1e-2, 1., 2., 1. }
   };

   for (const auto v: vals) {
      auto db0 = [v](double p2) { return ReB0(p2, v.m1, v.m2, v.q); };
      const auto v1 = ReD1B0(v.p, v.m1, v.m2);
      const auto v2 = derivative_central<3>(db0, v.p);

      BOOST_CHECK_CLOSE_FRACTION(v1, v2, 0.01);
   }
}

BOOST_AUTO_TEST_CASE( test_ReB1 )
{
   BOOST_CHECK_EQUAL(ReB1(0., 0., 0., scale2), 0.);

   BOOST_CHECK_CLOSE(ReB1(p2, 0,0, scale2), -0.5*ReB0(p2, 0,0, scale2), 1e-10);
   BOOST_CHECK_CLOSE(ReB1(0,p2,p2, scale2), -0.5*ReB0(0,p2,p2, scale2), 1e-10);

   BOOST_CHECK(ReB1(0,1e300,0,1e-10) != 0.);
}

BOOST_AUTO_TEST_CASE( test_ReB1_integral )
{
   const std::vector<B0_args> vals = {
      {0.  , 1.  , 0., 1.},
      {1e-1, 1.  , 0., 1.},
      {1e-2, 1.  , 0., 1.},
      {1e-3, 1.  , 0., 1.},
      {1e-4, 1.  , 0., 1.},
      {1e-5, 1.  , 0., 1.},
      {0.  , 1.  , 1e-15, 1.},
      {1e-1, 1.  , 1e-15, 1.},
      {1e-2, 1.  , 1e-15, 1.},
      {1e-3, 1.  , 1e-15, 1.},
      {1e-4, 1.  , 1e-15, 1.},
      {1e-5, 1.  , 1e-15, 1.},
      {0.  , 1e20, 0., 1.},
      {1e-1, 1e20, 0., 1.},
      {1e-2, 1e20, 0., 1.},
      {1e-3, 1e20, 0., 1.},
      {1e-4, 1e20, 0., 1.},
      {1e-5, 1e20, 0., 1.},
      {0.  , 0.  , 1., 1.},
      {1e-1, 0.  , 1., 1.},
      {1e-2, 0.  , 1., 1.},
      {1e-3, 0.  , 1., 1.},
      {1e-4, 0.  , 1., 1.},
      {1e-5, 0.  , 1., 1.},
      {0.  , 1e20, 1., 1.},
      {1e-1, 1e20, 1., 1.},
      {1e-2, 1e20, 1., 1.},
      {1e-3, 1e20, 1., 1.},
      {1e-4, 1e20, 1., 1.},
      {1e-5, 1e20, 1., 1.},
      {0.  , 0., 1e20, 1.},
      {1e-1, 0., 1e20, 1.},
      {1e-2, 0., 1e20, 1.},
      {1e-3, 0., 1e20, 1.},
      {1e-4, 0., 1e20, 1.},
      {1e-5, 0., 1e20, 1.},
      {1.  , 1.  , 1., 1.}
   };

   for (const auto v: vals) {
      const auto p = v.p;
      const auto m1 = v.m1;
      const auto m2 = v.m2;
      const auto q = v.q;

      const auto v1 = softsusy::b1(p,m1,m2,q);
      const auto v2 = bIntegral(p,m1,m2,q);
      BOOST_CHECK_CLOSE_FRACTION(v1, v2, 8e-5);
   }
}

BOOST_AUTO_TEST_CASE( test_ReB00 )
{
   BOOST_CHECK_EQUAL(ReB00(0., 0., 0., scale2), 0.);
}

BOOST_AUTO_TEST_CASE( test_ReB22 )
{
   BOOST_CHECK_EQUAL(ReB22(0., 0., 0., scale2), 0.);
}

BOOST_AUTO_TEST_CASE( test_ReH0 )
{
   BOOST_CHECK_EQUAL(ReH0(0., 0., 0., scale2), 0.);
}

BOOST_AUTO_TEST_CASE( test_ReF0 )
{
   BOOST_CHECK_EQUAL(ReF0(0., 0., 0., scale2), 0.);
}

BOOST_AUTO_TEST_CASE( test_ReD1F0_numerical )
{
   const std::vector<B0_args> vals = {
      {1e-2, 1.  , 0., 1. },
      {1e-2, 1.  , 1., 1. },
      {1e-2, 1.  , 1., 10.},
      {1e-2, 2.  , 1., 1. },
      {1e-2, 1.  , 2., 1. }
   };

   for (const auto v: vals) {
      auto df0 = [v](double p2) { return ReF0(p2, v.m1, v.m2, v.q); };
      const auto v1 = ReD1F0(v.p, v.m1, v.m2, v.q);
      const auto v2 = derivative_central<3>(df0, v.p);

      BOOST_CHECK_CLOSE_FRACTION(v1, v2, 0.01);
   }
}

BOOST_AUTO_TEST_CASE( test_ReG0 )
{
   BOOST_CHECK_EQUAL(ReG0(0., 0., 0., scale2), 0.);
}

BOOST_AUTO_TEST_CASE( test_ReD1G0_numerical )
{
   const std::vector<B0_args> vals = {
      {1e-2, 1.  , 0., 1. },
      {1e-2, 1.  , 1., 1. },
      {1e-2, 1.  , 1., 10.},
      {1e-2, 2.  , 1., 1. },
      {1e-2, 1.  , 2., 1. }
   };

   for (const auto v: vals) {
      auto dg0 = [v](double p2) { return ReG0(p2, v.m1, v.m2, v.q); };
      const auto v1 = ReD1G0(v.p, v.m1, v.m2, v.q);
      const auto v2 = derivative_central<3>(dg0, v.p);

      BOOST_CHECK_CLOSE_FRACTION(v1, v2, 0.01);
   }
}

#if defined(ENABLE_LOOPTOOLS) || defined(ENABLE_FFLITE)

const double eps = numeric_limits<double>::min();
const complex<double> neg_m2 = complex<double>(-0.1, -eps);

BOOST_AUTO_TEST_CASE(test_negative_m2)
{
    BOOST_CHECK_SMALL(abs(A0 (neg_m2,1) -
	complex<double>(-0.3302585092994046,-0.3141592653589793)), 1e-14);

    BOOST_CHECK_SMALL(abs(B0 (2,neg_m2,4, 1) -
	complex<double>(0.04191013584834996,0.14336849876951852)), 1e-14);
    BOOST_CHECK_SMALL(abs(B0 (2,3,neg_m2, 1) -
	complex<double>(0.49731133798115157,0.24955609428693504)), 1e-14);

    BOOST_CHECK_SMALL(abs(B1 (2,neg_m2,4, 1) -
	complex<double>(0.32573255511542315,-0.003271354485747599)), 1e-14);
    BOOST_CHECK_SMALL(abs(B1 (2,3,neg_m2, 1) -
	complex<double>(-0.6254665451021995,-0.2396442038760973)), 1e-14);

    BOOST_CHECK_SMALL(abs(B00(2,neg_m2,4, 1) -
	complex<double>(0.16595591599028536,-0.003633975888972291)), 1e-14);
    BOOST_CHECK_SMALL(abs(B00(2,3,neg_m2, 1) -
	complex<double>(0.2828439119832702,-0.006501356567577554)), 1e-14);
}

#endif
