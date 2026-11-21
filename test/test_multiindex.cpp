#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_multiindex

#include <boost/mpl/pair.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/test/unit_test.hpp>

#include "benchmark.hpp"
#include "multiindex.hpp"

#include <array>
#include <cstddef>
#include <iterator>
#include <type_traits>
#include <vector>

using namespace flexiblesusy;

/// object with one index in [start, end)
template <int start, int end>
struct Vertex {
   using index_bounds = boost::mpl::pair<boost::mpl::vector_c<int, start>,
                                         boost::mpl::vector_c<int, end>>;
};

/// object with two indices in [start1, end1) x [start2, end2)
template <int start1, int end1, int start2, int end2>
struct Vertex2 {
   using index_bounds =
      boost::mpl::pair<boost::mpl::vector_c<int, start1, start2>,
                       boost::mpl::vector_c<int, end1, end2>>;
};

/// object without indices
struct VertexNoIndices {
   using index_bounds =
      boost::mpl::pair<boost::mpl::vector_c<int>, boost::mpl::vector_c<int>>;
};

template <class Range>
auto collect(const Range& range)
{
   std::vector<
      typename std::iterator_traits<decltype(range.begin())>::value_type>
      result;
   for (const auto& idx : range) {
      result.push_back(idx);
   }
   return result;
}

BOOST_AUTO_TEST_CASE(test_multiindex_no_indices)
{
   const auto combinations = collect(index_range<VertexNoIndices>());

   BOOST_CHECK_EQUAL(combinations.size(), 1U);
}

BOOST_AUTO_TEST_CASE(test_multiindex_empty_range)
{
   BOOST_CHECK(collect(index_range<Vertex<0, 0>>()).empty());
   BOOST_CHECK(collect(index_range<Vertex<2, 1>>()).empty());
   BOOST_CHECK(collect(index_range<Vertex2<0, 3, 1, 1>>()).empty());
}

BOOST_AUTO_TEST_CASE(test_multiindex_3_generations)
{
   const auto combinations = collect(index_range<Vertex<0, 3>>());

   BOOST_REQUIRE_EQUAL(combinations.size(), 3U);
   for (int i = 0; i < 3; ++i) {
      BOOST_CHECK_EQUAL(combinations[i][0], i);
   }
}

BOOST_AUTO_TEST_CASE(test_multiindex_lower_bound)
{
   const auto combinations = collect(index_range<Vertex<2, 5>>());

   BOOST_REQUIRE_EQUAL(combinations.size(), 3U);
   BOOST_CHECK_EQUAL(combinations[0][0], 2);
   BOOST_CHECK_EQUAL(combinations[1][0], 3);
   BOOST_CHECK_EQUAL(combinations[2][0], 4);
}

BOOST_AUTO_TEST_CASE(test_multiindex_two_indices)
{
   // first index runs fastest
   const std::vector<std::array<int, 2>> expected{{0, 1}, {1, 1}, {2, 1},
                                                  {0, 2}, {1, 2}, {2, 2}};

   const auto combinations = collect(index_range<Vertex2<0, 3, 1, 3>>());

   BOOST_REQUIRE_EQUAL(combinations.size(), expected.size());
   for (std::size_t i = 0; i < expected.size(); ++i) {
      BOOST_CHECK(combinations[i] == expected[i]);
   }
}

BOOST_AUTO_TEST_CASE(test_multiindex_iterator)
{
   const auto range = index_range<Vertex<0, 2>>();

   auto it = range.begin();
   BOOST_CHECK(it != range.end());
   BOOST_CHECK_EQUAL((*it)[0], 0);
   BOOST_CHECK_EQUAL(it->size(), 1U);

   const auto old = it++;
   BOOST_CHECK_EQUAL((*old)[0], 0);
   BOOST_CHECK_EQUAL((*it)[0], 1);

   ++it;
   BOOST_CHECK(it == range.end());
}

BOOST_AUTO_TEST_CASE(test_multiindex_same_bounds_same_type)
{
   struct OtherVertex : Vertex<0, 3> {};

   BOOST_CHECK((std::is_same_v<decltype(index_range<Vertex<0, 3>>()),
                               decltype(index_range<OtherVertex>())>));
}

BOOST_AUTO_TEST_CASE(bench_multiindex_increment)
{
   constexpr unsigned N = 1'000'000;
   const auto ir = index_range<Vertex<0, N>>();
   unsigned count = 0;

   Stopwatch sw;
   sw.start();

   for ([[maybe_unused]] const auto& i : ir) {
      do_not_optimize(count++);
   }

   sw.stop();

   BOOST_CHECK_EQUAL(count, N);
   BOOST_TEST_MESSAGE("index loop (N = " << N << ") took "
                                         << sw.get_time_in_seconds()
                                         << " seconds");
}
