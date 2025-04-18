#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SM_observable_problems

#include <boost/test/unit_test.hpp>

#include "test_SM.hpp"
#include "lowe.h"
#include "observable_problems_format_slha.hpp"
#include "physical_input.hpp"
#include "spectrum_generator_settings.hpp"
#include "SM_mass_eigenstates.hpp"
#include "SM_input_parameters.hpp"
#include "SM_observables.hpp"
#include "observables/l_to_l_conversion/settings.hpp"
#include "config.h"


BOOST_AUTO_TEST_CASE( test_non_perturbative_running )
{
   softsusy::QedQcd qedqcd;
   flexiblesusy::Physical_input physical_input;
   flexiblesusy::SM_input_parameters input;
   flexiblesusy::SM_mass_eigenstates sm;
   flexiblesusy::Spectrum_generator_settings settings;
   const double scale = 0.0;

   setup_SM_const(sm, input);

#if defined(ENABLE_FEYNARTS) && defined(ENABLE_FORMCALC)
   flexiblesusy::LToLConversion_settings ltolconversion_settings;
   const auto obs = flexiblesusy::calculate_observables(sm, qedqcd, ltolconversion_settings, physical_input, settings, scale);
#else
   const auto obs = flexiblesusy::calculate_observables(sm, qedqcd, physical_input, settings, scale);
#endif

   const auto op = obs.problems;

   // BOOST_CHECK(op.have_problem());
   // BOOST_TEST_MESSAGE(op);
}
