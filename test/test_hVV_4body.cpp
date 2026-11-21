#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_hVV_4body

#include <boost/test/unit_test.hpp>

#include "decays/decay.hpp"

using namespace flexiblesusy;

namespace {

constexpr double mW = 80.385;
constexpr double GammaW = 2.085;
constexpr double mZ = 91.1876;
constexpr double GammaZ = 2.4952;

void check(double mH, double mV, double GammaV, double expected)
{
   const auto result = hVV_4body(mH, mV, GammaV);
   BOOST_CHECK_CLOSE_FRACTION(result.value, expected, 1e-6);
   // error estimate is below the threshold for the warning in the decay code
   BOOST_CHECK_LE(result.error, 1e-3*result.value);
}

} // anonymous namespace

// The reference values were obtained with a tolerance of 1e-9.  They agree
// with 1000000 point MISER integrations within the statistical errors.

BOOST_AUTO_TEST_CASE(test_below_single_vector_threshold)
{
   check(70., mW, GammaW, 3.330360444e-05);
}

BOOST_AUTO_TEST_CASE(test_between_thresholds)
{
   check(125.09, mW, GammaW, 0.01304931174);
   check(125.09, mZ, GammaZ, 0.003260272064);
   check(155., mW, GammaW, 0.1830731849);
}

BOOST_AUTO_TEST_CASE(test_close_to_pair_threshold)
{
   check(160.7, mW, GammaW, 0.6289885477);
   check(182., mZ, GammaZ, 0.5955580784);
}

BOOST_AUTO_TEST_CASE(test_deterministic)
{
   const auto first = hVV_4body(125.09, mW, GammaW);
   const auto second = hVV_4body(125.09, mW, GammaW);
   BOOST_CHECK_EQUAL(first.value, second.value);
}
