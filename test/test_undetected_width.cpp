#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_undetected_width

#include <boost/test/unit_test.hpp>

#include "decays/decay.hpp"

using namespace flexiblesusy;

namespace
{

/// undetected width for the given total and invisible widths, after
/// the widths of the channels passed to HiggsTools/Lilith were subtracted
double undetected_width(double total, double invisible, double subtracted)
{
   NeutralHiggsEffectiveCouplings effc;
   effc.particle = "hh(1)";
   effc.width = total;
   effc.invWidth = invisible;
   effc.undetectedWidth = -subtracted;
   effc.calculate_undetected_br();
   return effc.get_undetected_width();
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(test_positive_undetected_width_is_kept)
{
   BOOST_CHECK_CLOSE_FRACTION(undetected_width(4e-3, 1e-3, 2.5e-3), 0.5e-3,
                              1e-12);
}

BOOST_AUTO_TEST_CASE(test_rounding_level_negative_width_is_zero)
{
   BOOST_CHECK_EQUAL(undetected_width(4e-3, 1e-3, 3e-3 * (1 + 1e-14)), 0.);
}

BOOST_AUTO_TEST_CASE(test_larger_negative_width_is_zero)
{
   // e.g. a channel subtracted twice
   BOOST_CHECK_EQUAL(undetected_width(4e-3, 1e-3, 6e-3), 0.);
}

BOOST_AUTO_TEST_CASE(test_zero_total_width)
{
   BOOST_CHECK_EQUAL(undetected_width(0., 0., 1e-3), 0.);
   BOOST_CHECK_EQUAL(undetected_width(0., 0., 0.), 0.);
}
