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
   using namespace amm_loop_functions::two_loop;
   const double eps = 1e-14;
   const auto data = test::read_from_file<double>(
      std::string(TEST_DATA_DIR) + test::PATH_SEPARATOR + "fPS.txt");

   for (const auto d: data) {
      BOOST_REQUIRE(d.size() == 2);
      const double x = d.at(0);
      const double expected = d.at(1);
      const double value = 2*BarrZeeLoopFPS(x); // note the factor 2
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
   using namespace amm_loop_functions::two_loop;
   const double eps = 1e-14;
   const auto data = test::read_from_file<double>(
      std::string(TEST_DATA_DIR) + test::PATH_SEPARATOR + "fS.txt");

   for (const auto d: data) {
      BOOST_REQUIRE(d.size() == 2);
      const double x = d.at(0);
      const double expected = d.at(1);
      const double value = 2*BarrZeeLoopFS(x); // note the factor 2
      BOOST_CHECK_CLOSE_FRACTION(value, expected, eps);
   }
}

/**
 * Compares the Barr-Zee 2-loop function S (arXiv:1502.04199 Eq 27)
 * with generated data.
 *
 * Note: The data has been taken from GM2Calc, where Eq. (72)
 * arXiv:hep-ph/0609168 from arXiv:hep-ph/0609168 is implemented,
 * which differs by Eq. (27) from arXiv:1502.04199 by a prefactor
 * given by the mass ratio passed to the loop function.
 */
BOOST_AUTO_TEST_CASE( test_S )
{
   using namespace flexiblesusy;
   using namespace amm_loop_functions::two_loop;
   const double eps = 1e-14;
   const auto data = test::read_from_file<double>(
      std::string(TEST_DATA_DIR) + test::PATH_SEPARATOR + "fsferm.txt");

   for (const auto d: data) {
      BOOST_REQUIRE(d.size() == 2);
      const double x = d.at(0);
      const double expected = d.at(1);
      const double value = x*BarrZeeLoopS(x); // note the factor x
      BOOST_CHECK_CLOSE_FRACTION(value, expected, eps);
   }
}

/**
 * Compares the Barr-Zee 2-loop function V (arXiv:1502.04199 Eq 28)
 * with generated data.
 */
BOOST_AUTO_TEST_CASE( test_V )
{
   using namespace flexiblesusy;
   using namespace amm_loop_functions::two_loop;
   const double eps = 1e-14;
   const auto data = test::read_from_file<double>(
      std::string(TEST_DATA_DIR) + test::PATH_SEPARATOR + "fv.txt");

   for (const auto d: data) {
      BOOST_REQUIRE(d.size() == 2);
      const double x = d.at(0);
      const double expected = d.at(1);
      const double value = BarrZeeLoopV(x);
      BOOST_CHECK_CLOSE_FRACTION(value, expected, eps);
   }
}

/**
 * Compares the Barr-Zee 2-loop function FPSZ with generated data.
 */
BOOST_AUTO_TEST_CASE( test_FPSZ )
{
   using namespace flexiblesusy;
   using namespace amm_loop_functions::two_loop;
   const double eps = 1e-13;
   const auto data = test::read_from_file<double>(
      std::string(TEST_DATA_DIR) + test::PATH_SEPARATOR + "FPSZ.txt");

   for (const auto d: data) {
      BOOST_REQUIRE(d.size() == 3);
      const double x = d.at(0);
      const double y = d.at(1);
      const double expected = d.at(2);
      const double value = 2*BarrZeeLoopFPSZ(x, y);
      BOOST_CHECK_CLOSE_FRACTION(value, expected, eps);
   }
}

/**
 * Compares the Barr-Zee 2-loop function FSZ with generated data.
 */
BOOST_AUTO_TEST_CASE( test_FSZ )
{
   using namespace flexiblesusy;
   using namespace amm_loop_functions::two_loop;
   const double eps = 1e-13;
   const auto data = test::read_from_file<double>(
      std::string(TEST_DATA_DIR) + test::PATH_SEPARATOR + "FSZ.txt");

   for (const auto d: data) {
      BOOST_REQUIRE(d.size() == 3);
      const double x = d.at(0);
      const double y = d.at(1);
      const double expected = d.at(2);
      const double value = 2*BarrZeeLoopFSZ(x, y);
      BOOST_CHECK_CLOSE_FRACTION(value, expected, eps);
   }
}
