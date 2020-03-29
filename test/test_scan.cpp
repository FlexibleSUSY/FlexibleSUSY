#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_scan

#include <boost/test/unit_test.hpp>
#include <algorithm>
#include <random>

/**
 * @file scan.hpp
 * @brief contains helper functions and classes for parameter scans
 */

namespace flexiblesusy {

/**
 * @class Uniform
 * @brief real uniform distribution
 *
 * Usage:
 * @code
 * double x = Uniform<>::dice(0., 1.); // uses std::minstd_rand
 * double y = Uniform<std::mt19937>::dice(0., 1.); // uses std::mt19937
 * @endcode
 */
template <typename Generator = std::minstd_rand, typename RealType = double>
class Uniform {
public:
   /// returns random number between start and stop
   static RealType dice(RealType start, RealType stop) {
      static Generator generator;
      static std::uniform_real_distribution<RealType> distribution(start, stop);
      return distribution(generator);
   }
   /// returns random number between 0. and 1.
   static RealType dice() {
      static Generator generator;
      static std::uniform_real_distribution<RealType> distribution;
      return distribution(generator);
   }
};

} // namespace flexiblesusy

BOOST_AUTO_TEST_CASE( uniform_distribution_0_1 )
{
   using namespace flexiblesusy;
   const std::size_t N = 1000;
   std::array<double, N> a{};

   std::generate(std::begin(a), std::end(a),
                 [] () { return Uniform<>::dice(); });

   BOOST_CHECK(*std::max_element(std::begin(a), std::end(a)) <= 1.);
   BOOST_CHECK(*std::min_element(std::begin(a), std::end(a)) >= 0.);
}

BOOST_AUTO_TEST_CASE( uniform_distribution_start_stop )
{
   using namespace flexiblesusy;
   const std::size_t N = 1000;
   const double start = 0., stop = 10.;
   std::array<double, N> a{};

   std::generate(std::begin(a), std::end(a),
                 [start, stop] () { return Uniform<>::dice(start, stop); });

   BOOST_CHECK(*std::max_element(std::begin(a), std::end(a)) <= stop);
   BOOST_CHECK(*std::min_element(std::begin(a), std::end(a)) >= start);
}
