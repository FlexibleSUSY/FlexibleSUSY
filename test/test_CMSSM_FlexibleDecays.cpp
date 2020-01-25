
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSM_FlexibleDecays

#include <iomanip>
#include <cmath>

#include <boost/test/unit_test.hpp>

#include "CMSSM_two_scale_spectrum_generator.hpp"
#include "CMSSM_slha_io.hpp"
#include "CMSSM_decays.hpp"
#include "spectrum_generator_settings.hpp"

// #include "test_CMSSM.hpp"
// #include "CMSSM_two_scale_model.hpp"
// #include "CMSSM_two_scale_high_scale_constraint.hpp"
// #include "CMSSM_spectrum_generator.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_CMSSM_FlexibleDecays )
{
   // flexiblesusy::CMSSM_slha_io slha_io {};
   softsusy::QedQcd qedqcd;
   //qedqcd.setAlpha(softsusy::ALPHA, 1./1.2794001e+02);
   //qedqcd.setAlphaEmInput(1.0 / 1.2794001e+02);
   //qedqcd.setFermiConstant(1.1663787e-5);
   //qedqcd.setPoleMZ(91.1876);
   //qedqcd.setMass(softsusy::mTau, 1.77682);
   //qedqcd.setPoleMtau(1.77682);

   // input.TanBeta = 10.;
   // input.m0 = 125.;
   // input.m12 = 200.;
   // input.SignMu = 1;
   // input.Azero = 0.;
   // point from Peter

   CMSSM_input_parameters input;
   input.TanBeta = 49.8794511629;
   input.m0 = 9000.6279265;
   input.m12 = 2256.47237065;
   input.SignMu = -1;
   input.Azero = 9206.07878542;

   // try {
      // slha_io.fill(qedqcd);
      // slha_io.fill(input);
   // } catch (const Error& error) {
      // ERROR(error.what());
   //   return EXIT_FAILURE;
   // }


   Spectrum_generator_settings spectrum_generator_settings;
   //spectrum_generator_settings.set(Spectrum_generator_settings::calculate_sm_masses, 1.);
   spectrum_generator_settings.set(Spectrum_generator_settings::pole_mass_loop_order, 2.);
   spectrum_generator_settings.set(Spectrum_generator_settings::ewsb_loop_order, 2.);
   spectrum_generator_settings.set(Spectrum_generator_settings::beta_loop_order, 2.);
   spectrum_generator_settings.set(Spectrum_generator_settings::threshold_corrections_loop_order, 2.);
   spectrum_generator_settings.set(Spectrum_generator_settings::higgs_2loop_correction_at_as, 1.);
   spectrum_generator_settings.set(Spectrum_generator_settings::higgs_2loop_correction_ab_as, 1.);
   spectrum_generator_settings.set(Spectrum_generator_settings::higgs_2loop_correction_at_at, 1.);
   spectrum_generator_settings.set(Spectrum_generator_settings::higgs_2loop_correction_atau_atau, 1.);


   CMSSM_spectrum_generator<Two_scale> spectrum_generator;
   /*
   spectrum_generator.set_parameter_output_scale(
      //slha_io.get_parameter_output_scale()
      1e+3
      );
      */
   spectrum_generator.set_settings(spectrum_generator_settings);

   // solve BVP
   spectrum_generator.run(qedqcd, input);

   // get your model object:
   auto model = spectrum_generator.get_models_slha();
   auto m = std::get<0>(model);

   if (m.get_problems().have_problem()) {
      std::ostringstream ostr;
      m.get_problems().print_problems(ostr);
      BOOST_FAIL(ostr.str());
   }
   // get the problems:
   const auto& problems = spectrum_generator.get_problems();

/*
   CMSSM<Two_scale> m;
   m.set_input_parameters(input);
   CMSSM_high_scale_constraint<Two_scale> constraint(&m);

   // setup_CMSSM_const(m, input);
   m.calculate_DRbar_masses();

   m.set_pole_mass_loop_order(1);

   */
   CMSSM_decays decays(m, qedqcd, input);

   // check consistency of model parameters with the ones injected into HDECAY
   BOOST_CHECK_CLOSE_FRACTION(m.get_physical().MFe[2], 1.763858808146417,
                              1e-15);
   BOOST_CHECK_CLOSE_FRACTION(m.get_physical().Mhh[0], 90.7454505909439, 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(m.get_physical().Mhh[1], 765.475139007419, 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(m.get_physical().MHpm[1], 769.717887768395,
                              1e-15);
   BOOST_CHECK_CLOSE_FRACTION(m.get_physical().MAh[1], 765.24323718682024,
                              1e-15);
   BOOST_CHECK_CLOSE_FRACTION(m.get_physical().MVZ, 90.260768810873827, 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(
      m.get_MVZ(),
      1. / 2. * std::sqrt(pow(m.get_vu(), 2) + pow(m.get_vd(), 2)) *
         std::sqrt(3. / 5. * pow(m.get_g1(), 2) + pow(m.get_g2(), 2)),
      1e-15);
      std::cout <<m.get_g1() << '\n';

   // sin of the tree-level Higgs mixing angle (sometimes called alpha)
   const auto sin_alpha = 0.10228157265640475;
   BOOST_CHECK_CLOSE_FRACTION(m.get_ZH(0,0), sin_alpha, 1e-14);
   
   // --------------------------- light Higgs -------------------------------
   // h -> b bbar
   BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_barFdFd(m, 0, 2, 2),
                              2.6118180765322455E-003, 2e-15);
   // h -> tau+ tau-
   BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_barFeFe(m, 0, 2, 2),
                              1.9569212629432192E-004, 5e-15);
   // h -> W+ W-
   BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VWmconjVWm(m, 0),
                              8.4705126919250480E-004, 1e-14);
   // h -> Z Z
   BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VZVZ(m, 0),
                              8.5120714342385694E-014, 1e-14);
   // h -> gluon gluo
   // BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VGVG(m), , 1e-15);
   // h -> gamma gamma
   // BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VPVP(m), , 1e-15);
   // h -> Z gamma
   // BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VZVP(m), , 1e-15);

   // --------------------------- heavy Higgs -------------------------------
   // h -> b bbar
   // BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_barFdFd(m, 1, 2, 2),
                              // 2.6118180765322455E-003, 2e-15);
   // h -> tau+ tau-
   // BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_barFeFe(m, 1, 2, 2),
                              // 2.6800741077127161E-004, 1e-15);
   // h -> W+ W-
   // BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VWmconjVWm(m, 1),
                              // 8.4705126919250480E-004, 1e-14);
   // h -> Z Z
   // BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VZVZ(m, 1),
                              // 9.4231401208556083E-005, 1e-14);
   // h -> gluon gluon
   // BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VGVG(m), , 1e-15);
   // h -> gamma gamma
   // BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VPVP(m), , 1e-15);
   // h -> Z gamma
   // BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VZVP(m), , 1e-15);

   // --------------------------- pseudoscalar Higgs ------------------------
   // h -> b bbar
   // BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_Ah_to_barFdFd(m, 0, 2, 2),
                              // 2.6118180765322455E-003, 2e-15);
   // h -> tau+ tau-
   // BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_Ah_to_barFeFe(m, 0, 2, 2),
                              // 2.6800741077127161E-004, 1e-15);
   // h -> W+ W-
   // BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_Ah_to_VWmconjVWm(m, 0),
                              // 8.4705126919250480E-004, 1e-14);
   // h -> Z Z
   // BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_Ah_to_VZVZ(m, 0),
                              // 9.4231401208556083E-005, 1e-14);
   // h -> gluon gluon
   // BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VGVG(m), , 1e-15);
   // h -> gamma gamma
   // BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VPVP(m), , 1e-15);
   // h -> Z gamma
   // BOOST_CHECK_CLOSE_FRACTION(decays.partial_width_hh_to_VZVP(m), , 1e-15);
}
