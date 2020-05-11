#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_sum

#include <boost/test/unit_test.hpp>
#include "sum.hpp"
#include "stopwatch.hpp"
#include <complex>
#include <Eigen/Core>

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE(test_SUM_int)
{
   BOOST_CHECK_EQUAL(SUM(i,1,0,1), 0);
   BOOST_CHECK_EQUAL(SUM(i,1,1,1), 1);
   BOOST_CHECK_EQUAL(SUM(i,1,2,1), 2);
   BOOST_CHECK_EQUAL(SUM(i,1,3,1), 3);

   BOOST_CHECK_EQUAL(SUM(i,1,0,i), 0);
   BOOST_CHECK_EQUAL(SUM(i,1,1,i), 1);
   BOOST_CHECK_EQUAL(SUM(i,1,2,i), 3);
   BOOST_CHECK_EQUAL(SUM(i,1,3,i), 6);

   BOOST_CHECK_EQUAL(SUM(i,-1,0,i), -1);
   BOOST_CHECK_EQUAL(SUM(i,-1,1,i), 0);
   BOOST_CHECK_EQUAL(SUM(i,-1,2,i), 2);
   BOOST_CHECK_EQUAL(SUM(i,-1,3,i), 5);

   BOOST_CHECK_EQUAL(SUM(i,1,0,2*i+1), 0);
   BOOST_CHECK_EQUAL(SUM(i,1,1,2*i+1), 3);
   BOOST_CHECK_EQUAL(SUM(i,1,2,2*i+1), 8);
   BOOST_CHECK_EQUAL(SUM(i,1,3,2*i+1), 15);

   BOOST_CHECK_EQUAL(SUM(i,1,2,SUM(k,1,2,k*i)), 9);
   BOOST_CHECK_EQUAL(SUM(i,1,2,i*SUM(k,1,2,k*i)), 15);
}

BOOST_AUTO_TEST_CASE(test_SUM_double)
{
   BOOST_CHECK_EQUAL(SUM(i,1,0,1.), 0.);
   BOOST_CHECK_EQUAL(SUM(i,1,1,1.), 1.);
   BOOST_CHECK_EQUAL(SUM(i,1,2,1.), 2.);
   BOOST_CHECK_EQUAL(SUM(i,1,3,1.), 3.);

   BOOST_CHECK_EQUAL(SUM(i,1,0,1.*i), 0.);
   BOOST_CHECK_EQUAL(SUM(i,1,1,1.*i), 1.);
   BOOST_CHECK_EQUAL(SUM(i,1,2,1.*i), 3.);
   BOOST_CHECK_EQUAL(SUM(i,1,3,1.*i), 6.);

   BOOST_CHECK_EQUAL(SUM(i,1,2,SUM(k,1,2,1.5*k*i)), 13.5);
   BOOST_CHECK_EQUAL(SUM(i,1,2,SUM(k,i,3*(i+1),1.5*k*i)), 163.5);
}

BOOST_AUTO_TEST_CASE(test_SUM_complex)
{
   const std::complex<double> sum_1 = SUM(i,1,3,std::complex<double>(i,i));
   const std::complex<double> sum_2 = SUM(i,1,2,SUM(k,1,2,std::complex<double>(k*i,k*i)));

   BOOST_CHECK_EQUAL(sum_1.real(), 6.);
   BOOST_CHECK_EQUAL(sum_1.imag(), 6.);
   BOOST_CHECK_EQUAL(sum_2.real(), 9.);
   BOOST_CHECK_EQUAL(sum_2.imag(), 9.);
}

BOOST_AUTO_TEST_CASE(test_SUM_benchmark_large_sum)
{
   constexpr std::size_t N = 10000;
   Stopwatch stopwatch;

   stopwatch.start();
   const std::size_t result_SUM = SUM(i,1,N,SUM(k,1,N,k*i));
   stopwatch.stop();
   const double time_SUM = stopwatch.get_time_in_seconds();

   stopwatch.start();
   std::size_t result_loop = 0;
   for (std::size_t i = 1; i <= N; i++)
      for (std::size_t k = 1; k <= N; k++)
         result_loop += k*i;
   stopwatch.stop();
   const double time_loop = stopwatch.get_time_in_seconds();

   BOOST_TEST_MESSAGE("SUM  took " << time_SUM  << " seconds");
   BOOST_TEST_MESSAGE("loop took " << time_loop << " seconds");

   BOOST_CHECK_EQUAL(result_loop, result_SUM);
   // BOOST_CHECK_LT(10*time_loop, time_SUM);
}

BOOST_AUTO_TEST_CASE(test_SUM_benchmark_small_sum)
{
   constexpr std::size_t N = 10;
   constexpr std::size_t M = 1000000;
   Stopwatch stopwatch;

   stopwatch.start();
   std::size_t result_SUM = 0;
   for (std::size_t i = 0; i < M; i++)
      result_SUM += SUM(i,1,N,SUM(k,1,N,k*i));
   stopwatch.stop();
   const double time_SUM = stopwatch.get_time_in_seconds();

   stopwatch.start();
   std::size_t result_loop = 0;
   for (std::size_t l = 0; l < M; l++)
      for (std::size_t i = 1; i <= N; i++)
         for (std::size_t k = 1; k <= N; k++)
            result_loop += k*i;
   stopwatch.stop();
   const double time_loop = stopwatch.get_time_in_seconds();

   BOOST_TEST_MESSAGE("SUM  took " << time_SUM  << " seconds");
   BOOST_TEST_MESSAGE("loop took " << time_loop << " seconds");

   BOOST_CHECK_EQUAL(result_loop, result_SUM);

   // If the following test passes, we should use the SUM() macro
   // instead of an explicit sum over Feynman diagrams.
   // BOOST_CHECK_LE(10*time_SUM, time_loop);
}

BOOST_AUTO_TEST_CASE(test_sum_eigen_vector)
{
#define UVec(i) (Eigen::Matrix<double,2,1>::Unit(i))

   Eigen::Matrix<double,2,1> v;
   v.setZero();

   v += SUM(i,0,1, 2*(i+1)*UVec(i));

   BOOST_CHECK_EQUAL(v(0), 2.);
   BOOST_CHECK_EQUAL(v(1), 4.);

#undef UVec
}
