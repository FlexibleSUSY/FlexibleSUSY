#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_thread_pool

#include <boost/test/unit_test.hpp>

#include "thread_pool.hpp"
#include "stopwatch.hpp"

#include <chrono>
#include <thread>

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE(test_0_threads)
{
   int result = 0;

   {
      Thread_pool tp(0);
      tp.run_task([&result](){ result += 1; });
   }

   BOOST_CHECK_EQUAL(result, 1);
}

BOOST_AUTO_TEST_CASE(test_0_threads_packaged_task)
{
   int result = 0;

   Thread_pool tp(0);
   auto fut = tp.run_packaged_task([](){ return 1; });

   result += fut.get();

   BOOST_CHECK_EQUAL(result, 1);
}

BOOST_AUTO_TEST_CASE(test_1_thread)
{
   int result = 0;

   {
      Thread_pool tp(1);
      tp.run_task([&result](){ result += 1; });
   }

   BOOST_CHECK_EQUAL(result, 1);
}

BOOST_AUTO_TEST_CASE(test_1_thread_packaged_task)
{
   int result = 0;

   Thread_pool tp(1);
   auto fut = tp.run_packaged_task([](){ return 1; });

   result += fut.get();

   BOOST_CHECK_EQUAL(result, 1);
}

BOOST_AUTO_TEST_CASE(test_fill_array)
{
   std::array<double, 100> a{};

   {
      Thread_pool tp(std::thread::hardware_concurrency());
      for (std::size_t i = 0; i < a.size(); i++)
         tp.run_task([&a, i](){ a[i] = i; });
   }

   for (std::size_t i = 0; i < a.size(); i++)
      BOOST_CHECK_EQUAL(a[i], i);
}

BOOST_AUTO_TEST_CASE(test_fill_array_packaged_task)
{
   std::array<double, 100> a{};
   std::array<std::future<double>, 100> f{};

   Thread_pool tp(std::thread::hardware_concurrency());

   for (std::size_t i = 0; i < a.size(); i++)
      f[i] = tp.run_packaged_task([i](){ return 1.*i; });

   for (std::size_t i = 0; i < a.size(); i++)
      a[i] = f[i].get();

   for (std::size_t i = 0; i < a.size(); i++)
      BOOST_CHECK_EQUAL(a[i], i);
}

template <typename Array>
void fill_sequential(Array& a)
{
   for (typename Array::size_type i = 0; i < a.size(); i++) {
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
      a[i] = i;
   }
}

template <typename Array>
void fill_parallel(Array& a)
{
   Thread_pool tp(std::thread::hardware_concurrency());
   for (typename Array::size_type i = 0; i < a.size(); i++)
      tp.run_task([&a, i](){
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
            a[i] = i;
         });
}

template <typename F>
double measure_time(F&& f)
{
   flexiblesusy::Stopwatch s;
   s.start();
   f();
   s.stop();

   return s.get_time_in_seconds();
}

BOOST_AUTO_TEST_CASE(test_fill_array_benchmark)
{
   std::array<double, 10000> a{}, b{};

   const double time_parallel = measure_time([&a](){ fill_parallel(a); });
   const double time_sequential = measure_time([&b](){ fill_sequential(b); });

   for (std::size_t i = 0; i < a.size(); i++)
      BOOST_CHECK_EQUAL(a[i], b[i]);

   BOOST_TEST_MESSAGE("parallel array fill  : " << time_parallel << "s");
   BOOST_TEST_MESSAGE("sequential array fill: " << time_sequential << "s");

   BOOST_CHECK_LT(time_parallel, time_sequential);
}
