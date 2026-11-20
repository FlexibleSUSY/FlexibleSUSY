#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_decays_list

#include <boost/test/unit_test.hpp>

#include "decays/decay.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE(test_total_width_is_sum_of_absolute_widths)
{
   Decays_list decays(25);
   decays.set_decay(10., {5, -5}, "h -> b bbar");
   decays.set_decay(-1., {23, 22}, "h -> Z gamma");

   BOOST_CHECK_CLOSE_FRACTION(decays.get_total_width(), 11., 1e-15);
}

BOOST_AUTO_TEST_CASE(test_total_width_after_replacing_positive_width)
{
   Decays_list decays(25);
   decays.set_decay(10., {5, -5}, "h -> b bbar");
   decays.set_decay(1., {23, 22}, "h -> Z gamma");
   decays.set_decay(2., {23, 22}, "h -> Z gamma");

   BOOST_CHECK_CLOSE_FRACTION(decays.get_total_width(), 12., 1e-15);
}

BOOST_AUTO_TEST_CASE(test_total_width_after_replacing_negative_width)
{
   Decays_list decays(25);
   decays.set_decay(10., {5, -5}, "h -> b bbar");
   decays.set_decay(-1., {23, 22}, "h -> Z gamma");
   decays.set_decay(-1., {23, 22}, "h -> Z gamma");

   BOOST_CHECK_CLOSE_FRACTION(decays.get_total_width(), 11., 1e-15);
}
