#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_spectrum_generator_settings

#include <boost/test/unit_test.hpp>
#include "spectrum_generator_settings.hpp"
#include "error.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_get )
{
   Spectrum_generator_settings s;
   const auto es = s.get();

   BOOST_REQUIRE_EQUAL(Spectrum_generator_settings::NUMBER_OF_OPTIONS, es.size());

   for (int i = 0; i < Spectrum_generator_settings::NUMBER_OF_OPTIONS; i++) {
      BOOST_CHECK_EQUAL(s.get(static_cast<Spectrum_generator_settings::Settings>(i)), es(i));
   }
}

BOOST_AUTO_TEST_CASE( test_set )
{
   const int N = Spectrum_generator_settings::NUMBER_OF_OPTIONS;
   Spectrum_generator_settings::Settings_t es = Spectrum_generator_settings::Settings_t::Random(N);

   BOOST_REQUIRE_EQUAL(Spectrum_generator_settings::NUMBER_OF_OPTIONS, es.size());

   Spectrum_generator_settings s;
   s.set(es);

   for (int i = 0; i < Spectrum_generator_settings::NUMBER_OF_OPTIONS; i++) {
      BOOST_CHECK_EQUAL(s.get(static_cast<Spectrum_generator_settings::Settings>(i)), es(i));
   }
}

BOOST_AUTO_TEST_CASE( test_set_invalid )
{
   Spectrum_generator_settings s;

   // 0 [double > 0]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::precision, 1e-4));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::precision,  0.0), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::precision, -1.0), flexiblesusy::Error);

   // 1 [int]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::max_iterations, 0));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::max_iterations, 100));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::max_iterations,  100.1), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::max_iterations, -100.0), flexiblesusy::Error);

   // 2 [int]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::solver, 0));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::solver, 1));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::solver, 2));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::solver, -1), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::solver, 0.1), flexiblesusy::Error);

   // 3 [bool]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::calculate_sm_masses, 0));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::calculate_sm_masses, 1));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::calculate_sm_masses, 2), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::calculate_sm_masses, -1), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::calculate_sm_masses, 0.1), flexiblesusy::Error);

   // 4 [int]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::pole_mass_loop_order, 0));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::pole_mass_loop_order, 1));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::pole_mass_loop_order, 2));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::pole_mass_loop_order, -1), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::pole_mass_loop_order, 0.1), flexiblesusy::Error);

   // 5 [int]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::ewsb_loop_order, 0));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::ewsb_loop_order, 1));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::ewsb_loop_order, 2));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::ewsb_loop_order, -1), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::ewsb_loop_order, 0.1), flexiblesusy::Error);

   // 6 [int]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::beta_loop_order, 0));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::beta_loop_order, 1));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::beta_loop_order, 2));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::beta_loop_order, -1), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::beta_loop_order, 0.1), flexiblesusy::Error);

   // 7 [int]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::threshold_corrections_loop_order, 0));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::threshold_corrections_loop_order, 1));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::threshold_corrections_loop_order, 2));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::threshold_corrections_loop_order, -1), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::threshold_corrections_loop_order, 0.1), flexiblesusy::Error);

   // 8 [bool]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::higgs_2loop_correction_at_as, 0));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::higgs_2loop_correction_at_as, 1));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::higgs_2loop_correction_at_as, 2), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::higgs_2loop_correction_at_as, -1), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::higgs_2loop_correction_at_as, 0.1), flexiblesusy::Error);

   // 9 [bool]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::higgs_2loop_correction_ab_as, 0));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::higgs_2loop_correction_ab_as, 1));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::higgs_2loop_correction_ab_as, 2), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::higgs_2loop_correction_ab_as, -1), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::higgs_2loop_correction_ab_as, 0.1), flexiblesusy::Error);

   // 10 [bool]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::higgs_2loop_correction_at_at, 0));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::higgs_2loop_correction_at_at, 1));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::higgs_2loop_correction_at_at, 2), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::higgs_2loop_correction_at_at, -1), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::higgs_2loop_correction_at_at, 0.1), flexiblesusy::Error);

   // 11 [bool]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::higgs_2loop_correction_atau_atau, 0));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::higgs_2loop_correction_atau_atau, 1));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::higgs_2loop_correction_atau_atau, 2), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::higgs_2loop_correction_atau_atau, -1), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::higgs_2loop_correction_atau_atau, 0.1), flexiblesusy::Error);

   // 12 [bool]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::force_output, 0));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::force_output, 1));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::force_output, 2), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::force_output, -1), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::force_output, 0.1), flexiblesusy::Error);

   // 13 [int]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::top_pole_qcd_corrections, 0));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::top_pole_qcd_corrections, 1));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::top_pole_qcd_corrections, 2));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::top_pole_qcd_corrections, -1), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::top_pole_qcd_corrections, 0.1), flexiblesusy::Error);

   // 14 [double > 0]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::beta_zero_threshold, 1e-4));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::beta_zero_threshold, 0.0));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::beta_zero_threshold, -1.0), flexiblesusy::Error);

   // 15 [bool]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::calculate_observables, 0));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::calculate_observables, 1));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::calculate_observables, 2), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::calculate_observables, -1), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::calculate_observables, 0.1), flexiblesusy::Error);

   // 16 [bool]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::force_positive_masses, 0));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::force_positive_masses, 1));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::force_positive_masses, 2), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::force_positive_masses, -1), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::force_positive_masses, 0.1), flexiblesusy::Error);

   // 17 [double > 0]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::pole_mass_scale, 1000));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::pole_mass_scale, 0.0));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::pole_mass_scale, -1.0), flexiblesusy::Error);

   // 18 [double > 0]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::eft_pole_mass_scale, 1000));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::eft_pole_mass_scale, 0.0));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::eft_pole_mass_scale, -1.0), flexiblesusy::Error);

   // 19 [double > 0]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::eft_matching_scale, 1000));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::eft_matching_scale, 0.0));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::eft_matching_scale, -1.0), flexiblesusy::Error);

   // 20 [int]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::eft_matching_loop_order_up, 0));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::eft_matching_loop_order_up, 1));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::eft_matching_loop_order_up, 2));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::eft_matching_loop_order_up, -1), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::eft_matching_loop_order_up, 0.1), flexiblesusy::Error);

   // 21 [int]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::eft_matching_loop_order_down, 0));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::eft_matching_loop_order_down, 1));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::eft_matching_loop_order_down, 2));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::eft_matching_loop_order_down, -1), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::eft_matching_loop_order_down, 0.1), flexiblesusy::Error);

   // 22 [int]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::eft_higgs_index, 0));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::eft_higgs_index, 1));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::eft_higgs_index, 2));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::eft_higgs_index, -1), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::eft_higgs_index, 0.1), flexiblesusy::Error);

   // 23 [bool]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::calculate_bsm_masses, 0));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::calculate_bsm_masses, 1));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::calculate_bsm_masses, 2), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::calculate_bsm_masses, -1), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::calculate_bsm_masses, 0.1), flexiblesusy::Error);

   // 24 [int]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::threshold_corrections, 0));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::threshold_corrections, 1));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::threshold_corrections, 2));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::threshold_corrections, -1), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::threshold_corrections, 0.1), flexiblesusy::Error);

   // 25 [int]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::higgs_3loop_ren_scheme_atb_as2, 0));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::higgs_3loop_ren_scheme_atb_as2, 1));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::higgs_3loop_ren_scheme_atb_as2, 2));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::higgs_3loop_ren_scheme_atb_as2, 3), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::higgs_3loop_ren_scheme_atb_as2, -1), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::higgs_3loop_ren_scheme_atb_as2, 0.1), flexiblesusy::Error);

   // 26 [bool]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::higgs_3loop_correction_at_as2, 0));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::higgs_3loop_correction_at_as2, 1));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::higgs_3loop_correction_at_as2, 2), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::higgs_3loop_correction_at_as2, -1), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::higgs_3loop_correction_at_as2, 0.1), flexiblesusy::Error);

   // 27 [bool]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::higgs_3loop_correction_ab_as2, 0));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::higgs_3loop_correction_ab_as2, 1));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::higgs_3loop_correction_ab_as2, 2), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::higgs_3loop_correction_ab_as2, -1), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::higgs_3loop_correction_ab_as2, 0.1), flexiblesusy::Error);

   // 28 [bool]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::higgs_3loop_correction_at2_as, 0));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::higgs_3loop_correction_at2_as, 1));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::higgs_3loop_correction_at2_as, 2), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::higgs_3loop_correction_at2_as, -1), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::higgs_3loop_correction_at2_as, 0.1), flexiblesusy::Error);

   // 29 [bool]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::higgs_3loop_correction_at3, 0));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::higgs_3loop_correction_at3, 1));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::higgs_3loop_correction_at3, 2), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::higgs_3loop_correction_at3, -1), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::higgs_3loop_correction_at3, 0.1), flexiblesusy::Error);

   // 30 [bool]
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::higgs_4loop_correction_at_as3, 0));
   BOOST_CHECK_NO_THROW(s.set(Spectrum_generator_settings::higgs_4loop_correction_at_as3, 1));
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::higgs_4loop_correction_at_as3, 2), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::higgs_4loop_correction_at_as3, -1), flexiblesusy::Error);
   BOOST_CHECK_THROW(s.set(Spectrum_generator_settings::higgs_4loop_correction_at_as3, 0.1), flexiblesusy::Error);

}
