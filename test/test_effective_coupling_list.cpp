#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_effective_coupling_list

#include <boost/test/unit_test.hpp>

#include "decays/decay.hpp"

#include <array>
#include <complex>
#include <string>
#include <utility>

using namespace flexiblesusy;

namespace {

const std::pair<std::string, std::complex<double>> coupling {"", std::complex<double>(0.1, 0.2)};

// The coupling must not enter: both charge combinations are summed as widths.
EffectiveCoupling_list add_both_charges(int pdg1, int pdg2, double width1, double width2)
{
   EffectiveCoupling_list list;
   list.add_coupling("hh(1)", std::array<int, 2>{pdg1, -pdg2}, coupling, width1);
   list.add_coupling("hh(1)", std::array<int, 2>{-pdg1, pdg2}, coupling, width2);
   return list;
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(test_emu_width_is_sum_over_charges)
{
   const auto list = add_both_charges(11, 13, 0.3, 0.2);

   BOOST_CHECK_CLOSE_FRACTION(list[0].emu_width, 0.5, 1e-15);
   BOOST_CHECK_EQUAL(list[0].etau_width, 0.);
   BOOST_CHECK_EQUAL(list[0].mutau_width, 0.);
   BOOST_CHECK_CLOSE_FRACTION(list[0].get_undetected_width(), -0.5, 1e-15);
}

BOOST_AUTO_TEST_CASE(test_etau_width_is_sum_over_charges)
{
   const auto list = add_both_charges(11, 15, 0.3, 0.2);

   BOOST_CHECK_CLOSE_FRACTION(list[0].etau_width, 0.5, 1e-15);
   BOOST_CHECK_EQUAL(list[0].emu_width, 0.);
   BOOST_CHECK_EQUAL(list[0].mutau_width, 0.);
   BOOST_CHECK_CLOSE_FRACTION(list[0].get_undetected_width(), -0.5, 1e-15);
}

BOOST_AUTO_TEST_CASE(test_mutau_width_is_sum_over_charges)
{
   const auto list = add_both_charges(13, 15, 0.3, 0.2);

   BOOST_CHECK_CLOSE_FRACTION(list[0].mutau_width, 0.5, 1e-15);
   BOOST_CHECK_EQUAL(list[0].emu_width, 0.);
   BOOST_CHECK_EQUAL(list[0].etau_width, 0.);
   BOOST_CHECK_CLOSE_FRACTION(list[0].get_undetected_width(), -0.5, 1e-15);
}

BOOST_AUTO_TEST_CASE(test_widths_of_different_particles_are_separate)
{
   EffectiveCoupling_list list;
   list.add_coupling("hh(1)", std::array<int, 2>{11, -13}, coupling, 0.3);
   list.add_coupling("hh(2)", std::array<int, 2>{11, -13}, coupling, 0.7);

   BOOST_CHECK_CLOSE_FRACTION(list[0].emu_width, 0.3, 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(list[1].emu_width, 0.7, 1e-15);
}

BOOST_AUTO_TEST_CASE(test_flavour_diagonal_coupling_keeps_cp_parts)
{
   EffectiveCoupling_list list;
   list.add_coupling("hh(1)", std::array<int, 2>{-15, 15}, coupling, 0.3);

   BOOST_CHECK_EQUAL(list[0].tautau.second, coupling.second);
}
