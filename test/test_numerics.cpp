
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_numerics

#include <boost/test/unit_test.hpp>
#include "numerics2.hpp"
#include "stopwatch.hpp"
#include <complex>
#include <cstdlib>

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE(test_is_equal)
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

   // BOOST_CHECK(is_equal_rel(s, (short)0));
   BOOST_CHECK(is_equal_rel(i, 0));
   BOOST_CHECK(is_equal_rel(l, 0L));
   BOOST_CHECK(is_equal_rel(ll, 0LL));
   // BOOST_CHECK(is_equal_rel(us, (unsigned short)0));
   BOOST_CHECK(is_equal_rel(ui, 0U));
   BOOST_CHECK(is_equal_rel(ul, 0UL));
   BOOST_CHECK(is_equal_rel(ull, 0ULL));
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
