#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_multiindex

#include <boost/test/unit_test.hpp>
#include <boost/mpl/pair.hpp>
#include <boost/mpl/vector_c.hpp>

#include "benchmark.hpp"
#include "multiindex.hpp"

template <int start, int end>
struct Vertex {
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, start>,
     boost::mpl::vector_c<int, end>
   >;
};


BOOST_AUTO_TEST_CASE( test_multiindex_empty )
{
   using namespace flexiblesusy;

   auto ir = index_range<Vertex<0,0> >();

   unsigned count = 0;

   for (const auto& i: ir)
      count++;

   BOOST_CHECK_EQUAL(count, 1u);
}


BOOST_AUTO_TEST_CASE( test_multiindex_3_generations )
{
   using namespace flexiblesusy;

   auto ir = index_range<Vertex<0,3> >();

   unsigned count = 0;

   for (const auto& i: ir)
      count++;

   BOOST_CHECK_EQUAL(count, 3u);
}


BOOST_AUTO_TEST_CASE( bench_multiindex_increment )
{
   using namespace flexiblesusy;

   constexpr unsigned N = 1'000'000;
   auto ir = index_range<Vertex<0,N>>();
   unsigned count = 0;

   Stopwatch sw;
   sw.start();

   for (const auto& i: ir)
      do_not_optimize(count++);

   sw.stop();

   BOOST_CHECK_EQUAL(count, N);
   BOOST_TEST_MESSAGE("index loop (N = " << N << ") took "
                      << sw.get_time_in_seconds() << " seconds");
}
