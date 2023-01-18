#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_amm_loop_functions

#include <boost/test/unit_test.hpp>
#include "amm_loop_functions.hpp"
#include "read_data.hpp"

/**
 * Compares the Barr-Zee 2-loop function FPS (arXiv:1502.04199 Eq 26)
 * with generated data.
 *
 * Note: The data has been taken from GM2Calc, where Eq. (70) from
 * arXiv:hep-ph/0609168 is implemented, which differs by
 * Eq. (26) from arXiv:1502.04199 by a factor 2.
 */
BOOST_AUTO_TEST_CASE( test_FPS )
{
   using namespace flexiblesusy;
   const double eps = 1e-14;
   const auto data = test::read_from_file<double>(
      std::string(TEST_DATA_DIR) + test::PATH_SEPARATOR + "fPS.txt");

   for (const auto d: data) {
      BOOST_REQUIRE(d.size() == 2);
      const double x = d.at(0);
      const double expected = d.at(1);
      const double value = 2*BarrZeeLoopFPS(x);
      BOOST_CHECK_CLOSE_FRACTION(value, expected, eps);
   }
}

/**
 * Compares the Barr-Zee 2-loop function FS (arXiv:1502.04199 Eq 25)
 * with generated data.
 *
 * Note: The data has been taken from GM2Calc, where Eq. (71)
 * arXiv:hep-ph/0609168 from arXiv:hep-ph/0609168 is implemented,
 * which differs by Eq. (25) from arXiv:1502.04199 by a factor 2.
 */
BOOST_AUTO_TEST_CASE( test_FS )
{
   using namespace flexiblesusy;
   const double eps = 1e-14;
   const auto data = test::read_from_file<double>(
      std::string(TEST_DATA_DIR) + test::PATH_SEPARATOR + "fS.txt");

   for (const auto d: data) {
      BOOST_REQUIRE(d.size() == 2);
      const double x = d.at(0);
      const double expected = d.at(1);
      const double value = 2*BarrZeeLoopFS(x);
      BOOST_CHECK_CLOSE_FRACTION(value, expected, eps);
   }
}
