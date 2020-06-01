
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_numerics

#include <boost/test/unit_test.hpp>
#include "numerics2.hpp"
#include "stopwatch.hpp"
#include <complex>
#include <cstdlib>

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE(test_is_equal_zero)
{
   short s = (short)0;
   int i = 0;
   long l = 0L;
   long long ll = 0LL;
   unsigned short us = 0;
   unsigned int ui = 0;
   unsigned long ul = 0UL;
   unsigned long long ull = 0ULL;
   std::complex<double> c(0.0, 0.0);

   BOOST_CHECK(is_zero(s));
   BOOST_CHECK(is_zero(i));
   BOOST_CHECK(is_zero(l));
   BOOST_CHECK(is_zero(ll));
   BOOST_CHECK(is_zero(us));
   BOOST_CHECK(is_zero(ui));
   BOOST_CHECK(is_zero(ul));
   BOOST_CHECK(is_zero(ull));
   BOOST_CHECK(is_zero(c));

   // BOOST_CHECK(is_equal(s, (short)0));
   BOOST_CHECK(is_equal(i, 0));
   BOOST_CHECK(is_equal(l, 0L));
   BOOST_CHECK(is_equal(ll, 0LL));
   // BOOST_CHECK(is_equal(us, (unsigned short)0));
   BOOST_CHECK(is_equal(ui, 0U));
   BOOST_CHECK(is_equal(ul, 0UL));
   BOOST_CHECK(is_equal(ull, 0ULL));
   BOOST_CHECK(is_equal(c, std::complex<double>(0.0, 0.0)));

   // BOOST_CHECK(is_equal_rel(s, (short)0));
   BOOST_CHECK(is_equal_rel(i, 0));
   BOOST_CHECK(is_equal_rel(l, 0L));
   BOOST_CHECK(is_equal_rel(ll, 0LL));
   // BOOST_CHECK(is_equal_rel(us, (unsigned short)0));
   BOOST_CHECK(is_equal_rel(ui, 0U));
   BOOST_CHECK(is_equal_rel(ul, 0UL));
   BOOST_CHECK(is_equal_rel(ull, 0ULL));
}


BOOST_AUTO_TEST_CASE(test_is_finite)
{
   const std::array<double, 4> sa = { 1.0, 2.0, 3.0, 4.0 };
   const double ca[4] = { 1.0, 2.0, 3.0, 4.0 };
   const double* pa = &ca[0];

   BOOST_CHECK(is_finite(sa));
   BOOST_CHECK(is_finite(ca));
   BOOST_CHECK(is_finite(pa, 4));
}


BOOST_AUTO_TEST_CASE(test_is_finite_inf)
{
   const double inf = std::numeric_limits<double>::infinity();
   const std::array<double, 4> sa = { 1.0, 2.0, 3.0, inf };
   const double ca[4] = { 1.0, 2.0, 3.0, inf };
   const double* pa = &ca[0];

   BOOST_CHECK(!is_finite(sa));
   BOOST_CHECK(!is_finite(ca));
   BOOST_CHECK(!is_finite(pa, 4));
}


BOOST_AUTO_TEST_CASE(test_is_finite_nan)
{
   const double nan = std::numeric_limits<double>::quiet_NaN();
   const std::array<double, 4> sa = { 1.0, 2.0, 3.0, nan };
   const double ca[4] = { 1.0, 2.0, 3.0, nan };
   const double* pa = &ca[0];

   BOOST_CHECK(!is_finite(sa));
   BOOST_CHECK(!is_finite(ca));
   BOOST_CHECK(!is_finite(pa, 4));
}


BOOST_AUTO_TEST_CASE(test_is_zero)
{
   const double epss[] = {
      1e-00, 1e-01, 1e-02, 1e-03, 1e-04,
      1e-05, 1e-06, 1e-07, 1e-08, 1e-09,
      1e-10, 1e-11, 1e-12, 1e-13, 1e-14,
      1e-15, 1e-16, 1e-17, 1e-18, std::numeric_limits<double>::epsilon()};

   const int E = sizeof(epss) / sizeof(epss[0]);

   const double nums[] = {1e-00, 1e-01, 1e-02, 1e-03, 1e-04, 1e-05, 1e-06,
                          1e-07, 1e-08, 1e-09, 1e-10, 1e-11, 1e-12, 1e-13,
                          1e-14, 1e-15, 1e-16, 1e-17, 1e-18};

   const int N = sizeof(nums) / sizeof(nums[0]);

   for (int e = 0; e < E; e++) {
      for (int n = 0; n < N; n++) {
         const double num = nums[n];
         const double eps = epss[n];

         if (std::abs(num) <= eps) {
            BOOST_CHECK(is_zero(num, eps));
            BOOST_CHECK(is_zero(-num, eps));
         } else {
            BOOST_CHECK(!is_zero(num, eps));
            BOOST_CHECK(!is_zero(-num, eps));
         }
      }
   }
}


BOOST_AUTO_TEST_CASE(test_is_equal_rel)
{
   double x = 1.0, eps = 0.5;
   BOOST_CHECK(is_equal_rel(x, x + eps, 1.1*eps));
   BOOST_CHECK(is_equal_rel(x, x - eps, 1.1*eps));

   x = 1.0, eps = 1e-10;
   BOOST_CHECK(is_equal_rel(x, x + eps, 1.1*eps));
   BOOST_CHECK(is_equal_rel(x, x - eps, 1.1*eps));

   x = 1.0, eps = std::numeric_limits<double>::epsilon();
   BOOST_CHECK(is_equal_rel(x, x + eps, 1.1*eps));
   BOOST_CHECK(is_equal_rel(x, x - eps, 1.1*eps));

   x = 0.0, eps = std::numeric_limits<double>::epsilon();
   BOOST_CHECK(is_equal_rel(x, x + eps, 1.1*eps));
   BOOST_CHECK(is_equal_rel(x, x - eps, 1.1*eps));

   eps = std::numeric_limits<double>::epsilon();
   BOOST_CHECK(is_equal_rel(eps, eps, 1.1*eps));

   // specific cases where numbers with small differences are treated as equal
   x = eps = std::numeric_limits<double>::epsilon();
   BOOST_CHECK(is_equal_rel(x, x + eps, 1.1*eps));
   BOOST_CHECK(is_equal_rel(x, x - eps, 1.1*eps));

   eps = std::numeric_limits<double>::epsilon();
   BOOST_CHECK(is_equal_rel(eps, 0.1*eps, 0.01*eps));
   BOOST_CHECK(is_equal_rel(1e-18, 1e-19, 1e-20));
}


BOOST_AUTO_TEST_CASE(test_is_equal_fraction)
{
   double x = 1.0, eps = 0.5;
   BOOST_CHECK(is_equal_fraction(x, x + eps, 1.1*eps));
   BOOST_CHECK(is_equal_fraction(x, x - eps, 1.1*eps));

   x = 1.0, eps = 1e-10;
   BOOST_CHECK(is_equal_fraction(x, x + eps, 1.1*eps));
   BOOST_CHECK(is_equal_fraction(x, x - eps, 1.1*eps));

   x = 1.0, eps = std::numeric_limits<double>::epsilon();
   BOOST_CHECK(is_equal_fraction(x, x + eps, 1.1*eps));
   BOOST_CHECK(is_equal_fraction(x, x - eps, 1.1*eps));

   x = 0.0, eps = std::numeric_limits<double>::epsilon();
   BOOST_CHECK(!is_equal_fraction(x, x + eps, 1.1*eps));
   BOOST_CHECK(!is_equal_fraction(x, x - eps, 1.1*eps));

   eps = std::numeric_limits<double>::epsilon();
   BOOST_CHECK(is_equal_fraction(eps, eps, 1.1*eps));

   x = eps = std::numeric_limits<double>::epsilon();
   BOOST_CHECK(!is_equal_fraction(x, x + eps, 1.1*eps));
   BOOST_CHECK(!is_equal_fraction(x, x - eps, 1.1*eps));

   eps = std::numeric_limits<double>::epsilon();
   BOOST_CHECK(!is_equal_fraction(eps, 0.1*eps, 0.01*eps));
   BOOST_CHECK(!is_equal_fraction(1e-18, 1e-19, 1e-20));
}


template <long N>
std::array<std::complex<double>, N> make_logs()
{
   std::array<std::complex<double>, N> vals;

   for (auto& v: vals)
      v = std::complex<double>(1.*std::rand()/RAND_MAX, 1.*std::rand()/RAND_MAX);

   return vals;
}

BOOST_AUTO_TEST_CASE(test_complex_log)
{
   const auto vals = make_logs<10000>();

   for (auto v: vals) {
      BOOST_CHECK_CLOSE_FRACTION(std::real(std::log(v)), std::real(fast_log(v)), 1e-9);
      BOOST_CHECK_CLOSE_FRACTION(std::imag(std::log(v)), std::imag(fast_log(v)), 1e-9);
   }
}

BOOST_AUTO_TEST_CASE(test_complex_log_bench)
{
   const auto vals = make_logs<10000>();

   Stopwatch sw;

   sw.start();
   for (auto v: vals)
      volatile auto res = std::log(v);
   sw.stop();
   const double t_log = sw.get_time_in_seconds();

   sw.start();
   for (auto v: vals)
      volatile auto res = fast_log(v);
   sw.stop();
   const double t_fslog = sw.get_time_in_seconds();

   BOOST_TEST_MESSAGE("time for log     : " << t_log);
   BOOST_TEST_MESSAGE("time for fast_log: " << t_fslog);

   // BOOST_CHECK_LT(t_fslog, t_log);
}
