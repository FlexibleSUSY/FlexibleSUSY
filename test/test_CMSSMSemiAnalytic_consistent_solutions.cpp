#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSMSemiAnalytic_consistent_solutions

#include <boost/test/unit_test.hpp>

#define private public

#include "CMSSMSemiAnalytic_input_parameters.hpp"
#include "CMSSMSemiAnalytic_slha_io.hpp"
#include "CMSSMSemiAnalytic_semi_analytic_ewsb_solver.hpp"
#include "CMSSMSemiAnalytic_semi_analytic_high_scale_constraint.hpp"
#include "CMSSMSemiAnalytic_semi_analytic_low_scale_constraint.hpp"
#include "CMSSMSemiAnalytic_semi_analytic_soft_parameters_constraint.hpp"
#include "CMSSMSemiAnalytic_semi_analytic_spectrum_generator.hpp"
#include "CMSSMSemiAnalytic_semi_analytic_susy_convergence_tester.hpp"
#include "CMSSMSemiAnalytic_semi_analytic_susy_scale_constraint.hpp"

#include "CMSSM_input_parameters.hpp"
#include "CMSSM_slha_io.hpp"
#include "CMSSM_two_scale_high_scale_constraint.hpp"
#include "CMSSM_two_scale_low_scale_constraint.hpp"
#include "CMSSM_two_scale_spectrum_generator.hpp"
#include "CMSSM_two_scale_susy_scale_constraint.hpp"

#include "two_scale_running_precision.hpp"
#include "two_scale_solver.hpp"

using namespace flexiblesusy;

template <class NewModel, class OldModel>
NewModel copy_parameters_from_model(const OldModel& old_model)
{
   NewModel new_model;

   new_model.set_loops(old_model.get_loops());
   new_model.set_scale(old_model.get_scale());
   new_model.set_ewsb_loop_order(old_model.get_ewsb_loop_order());
   new_model.set_thresholds(old_model.get_thresholds());

   new_model.set_Yu(old_model.get_Yu());
   new_model.set_Yd(old_model.get_Yd());
   new_model.set_Ye(old_model.get_Ye());

   new_model.set_g1(old_model.get_g1());
   new_model.set_g2(old_model.get_g2());
   new_model.set_g3(old_model.get_g3());

   new_model.set_Mu(old_model.get_Mu());

   new_model.set_vd(old_model.get_vd());
   new_model.set_vu(old_model.get_vu());

   new_model.set_TYu(old_model.get_TYu());
   new_model.set_TYd(old_model.get_TYd());
   new_model.set_TYe(old_model.get_TYe());

   new_model.set_mHd2(old_model.get_mHd2());
   new_model.set_mHu2(old_model.get_mHu2());
   new_model.set_BMu(old_model.get_BMu());

   new_model.set_mq2(old_model.get_mq2());
   new_model.set_mu2(old_model.get_mu2());
   new_model.set_md2(old_model.get_md2());
   new_model.set_ml2(old_model.get_ml2());
   new_model.set_me2(old_model.get_me2());

   new_model.set_MassB(old_model.get_MassB());
   new_model.set_MassWB(old_model.get_MassWB());
   new_model.set_MassG(old_model.get_MassG());

   return new_model;
}

CMSSM_input_parameters initialize_two_scale_input(
   double m0Sq, const CMSSMSemiAnalytic_input_parameters& semi_analytic_input)
{
   CMSSM_input_parameters input;

   input.m0 = Sqrt(m0Sq);
   input.m12 = semi_analytic_input.m12;
   input.Azero = semi_analytic_input.Azero;
   input.TanBeta = semi_analytic_input.TanBeta;
   input.SignMu = Sign(semi_analytic_input.MuInput);

   return input;
}

CMSSM<Two_scale> initialize_two_scale_model(
   const CMSSMSemiAnalytic<Semi_analytic>& semi_analytic_model,
   const CMSSMSemiAnalytic_input_parameters& semi_analytic_input)
{
   CMSSM<Two_scale> two_scale_model
      = copy_parameters_from_model<CMSSM<Two_scale> >(semi_analytic_model);

   two_scale_model.set_input_parameters(
      initialize_two_scale_input(semi_analytic_model.get_m0Sq(),
                                 semi_analytic_input));

   two_scale_model.calculate_DRbar_masses();

   return two_scale_model;
}

CMSSMSemiAnalytic_input_parameters initialize_semi_analytic_input(
   double MuInput, const CMSSM_input_parameters& two_scale_input)
{
   CMSSMSemiAnalytic_input_parameters input;

   input.m12 = two_scale_input.m12;
   input.Azero = two_scale_input.Azero;
   input.TanBeta = two_scale_input.TanBeta;
   input.MuInput = MuInput;

   return input;
}

CMSSMSemiAnalytic<Semi_analytic> initialize_semi_analytic_model(
   const CMSSM<Two_scale>& two_scale_model,
   const CMSSM_input_parameters& two_scale_input,
   double Mu0, double BMu0, double high_scale)
{
   CMSSMSemiAnalytic<Semi_analytic> semi_analytic_model
      = copy_parameters_from_model<CMSSMSemiAnalytic<Semi_analytic> >(two_scale_model);

   semi_analytic_model.set_input_parameters(
      initialize_semi_analytic_input(two_scale_model.get_Mu(),
                                     two_scale_input));

   semi_analytic_model.set_m0Sq(Sqr(two_scale_input.m0));
   semi_analytic_model.set_BMu0(BMu0);
   semi_analytic_model.set_MuBV(Mu0);

   semi_analytic_model.calculate_semi_analytic_solutions(high_scale);
   semi_analytic_model.calculate_DRbar_masses();

   return semi_analytic_model;
}

template <class Model>
double max_susy_parameters_rel_diff(const Model& old_model, const Model& new_model)
{
   std::array<double, 33> diff{};

   diff[0] = MaxRelDiff(old_model.get_g1(), new_model.get_g1());
   diff[1] = MaxRelDiff(old_model.get_g2(), new_model.get_g2());
   diff[2] = MaxRelDiff(old_model.get_g3(), new_model.get_g3());
   diff[3] = MaxRelDiff(old_model.get_vd(), new_model.get_vd());
   diff[4] = MaxRelDiff(old_model.get_vu(), new_model.get_vu());
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 5] = MaxRelDiff(old_model.get_Yd(i,j), new_model.get_Yd(i,j));
      }
   }
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 14] = MaxRelDiff(old_model.get_Ye(i,j), new_model.get_Ye(i,j));
      }
   }
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 23] = MaxRelDiff(old_model.get_Yu(i,j), new_model.get_Yu(i,j));
      }
   }
   diff[32] = MaxRelDiff(old_model.get_Mu(), new_model.get_Mu());

   return *std::max_element(diff.cbegin(), diff.cend());
}

template <class Model>
double max_soft_parameters_rel_diff(const Model& old_model, const Model& new_model)
{
   std::array<double,78> diff{};

   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j] = MaxRelDiff(old_model.get_TYd(i,j), new_model.get_TYd(i,j));
      }
   }
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 9] = MaxRelDiff(old_model.get_TYe(i,j), new_model.get_TYe(i,j));
      }
   }
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 18] = MaxRelDiff(old_model.get_TYu(i,j), new_model.get_TYu(i,j));
      }
   }
   diff[27] = MaxRelDiff(old_model.get_BMu(), new_model.get_BMu());
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 28] = MaxRelDiff(old_model.get_mq2(i,j), new_model.get_mq2(i,j));
      }
   }
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 37] = MaxRelDiff(old_model.get_ml2(i,j), new_model.get_ml2(i,j));
      }
   }
   diff[46] = MaxRelDiff(old_model.get_mHd2(), new_model.get_mHd2());
   diff[47] = MaxRelDiff(old_model.get_mHu2(), new_model.get_mHu2());
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 48] = MaxRelDiff(old_model.get_md2(i,j), new_model.get_md2(i,j));
      }
   }
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 57] = MaxRelDiff(old_model.get_mu2(i,j), new_model.get_mu2(i,j));
      }
   }
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 66] = MaxRelDiff(old_model.get_me2(i,j), new_model.get_me2(i,j));
      }
   }
   diff[75] = MaxRelDiff(old_model.get_MassB(), new_model.get_MassB());
   diff[76] = MaxRelDiff(old_model.get_MassWB(), new_model.get_MassWB());
   diff[77] = MaxRelDiff(old_model.get_MassG(), new_model.get_MassG());

   return *std::max_element(diff.cbegin(), diff.cend());
}

template <class Model>
double max_mass_rel_diff(const Model& old_model, const Model& new_model)
{
   std::array<double, 34> diff{};

   diff[0] = MaxRelDiff(old_model.get_MGlu(), new_model.get_MGlu());
   for (int i = 0; i < 6; ++i) {
      diff[i + 1] = MaxRelDiff(old_model.get_MSd(i), new_model.get_MSd(i));
   }
   for (int i = 0; i < 3; ++i) {
      diff[i + 7] = MaxRelDiff(old_model.get_MSv(i), new_model.get_MSv(i));
   }
   for (int i = 0; i < 6; ++i) {
      diff[i + 10] = MaxRelDiff(old_model.get_MSu(i), new_model.get_MSu(i));
   }
   for (int i = 0; i < 6; ++i) {
      diff[i + 16] = MaxRelDiff(old_model.get_MSe(i), new_model.get_MSe(i));
   }
   for (int i = 0; i < 2; ++i) {
      diff[i + 22] = MaxRelDiff(old_model.get_Mhh(i), new_model.get_Mhh(i));
   }
   for (int i = 1; i < 2; ++i) {
      diff[i + 24] = MaxRelDiff(old_model.get_MAh(i), new_model.get_MAh(i));
   }
   for (int i = 1; i < 2; ++i) {
      diff[i + 26] = MaxRelDiff(old_model.get_MHpm(i), new_model.get_MHpm(i));
   }
   for (int i = 0; i < 4; ++i) {
      diff[i + 28] = MaxRelDiff(old_model.get_MChi(i), new_model.get_MChi(i));
   }
   for (int i = 0; i < 2; ++i) {
      diff[i + 32] = MaxRelDiff(old_model.get_MCha(i), new_model.get_MCha(i));
   }

   return *std::max_element(diff.cbegin(), diff.cend());
}

CMSSM<Two_scale> run_single_two_scale_iteration(const CMSSM<Two_scale>& model, const CMSSM_scales& scales)
{
   CMSSM<Two_scale> next_model(model);

   CMSSM_high_scale_constraint<Two_scale> high_scale_constraint;
   CMSSM_susy_scale_constraint<Two_scale> susy_scale_constraint;
   CMSSM_low_scale_constraint<Two_scale> low_scale_constraint;

   high_scale_constraint.set_model(&next_model);
   susy_scale_constraint.set_model(&next_model);
   low_scale_constraint.set_model(&next_model);

   softsusy::QedQcd qedqcd;
   qedqcd.to(qedqcd.displayPoleMZ());

   susy_scale_constraint.set_sm_parameters(qedqcd);
   low_scale_constraint.set_sm_parameters(qedqcd);

   high_scale_constraint.initialize();
   susy_scale_constraint.initialize();
   low_scale_constraint.initialize();

   high_scale_constraint.scale = scales.HighScale;
   susy_scale_constraint.scale = scales.SUSYScale;
   low_scale_constraint.scale = scales.LowScale;

   // apply constraints once
   next_model.run_to(scales.LowScale);
   low_scale_constraint.apply();

   next_model.run_to(scales.HighScale);
   high_scale_constraint.apply();

   next_model.run_to(scales.SUSYScale);
   susy_scale_constraint.apply();

   return next_model;
}

CMSSMSemiAnalytic<Semi_analytic> run_single_semi_analytic_iteration(
   const CMSSMSemiAnalytic<Semi_analytic>& model,
   const CMSSMSemiAnalytic_scales& scales,
   double precision)
{
   CMSSMSemiAnalytic<Semi_analytic> next_model(model);

   CMSSMSemiAnalytic_semi_analytic_solutions& solutions(
      next_model.get_semi_analytic_solutions());

   CMSSMSemiAnalytic_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_semi_analytic_solutions(&solutions);
   next_model.set_ewsb_solver(
      std::make_shared<CMSSMSemiAnalytic_ewsb_solver<Semi_analytic> >(ewsb_solver));

   CMSSMSemiAnalytic_high_scale_constraint<Semi_analytic> high_scale_constraint;
   CMSSMSemiAnalytic_susy_scale_constraint<Semi_analytic> susy_scale_constraint;
   CMSSMSemiAnalytic_low_scale_constraint<Semi_analytic> low_scale_constraint;
   CMSSMSemiAnalytic_soft_parameters_constraint<Semi_analytic> soft_constraint;

   high_scale_constraint.set_model(&next_model);
   susy_scale_constraint.set_model(&next_model);
   low_scale_constraint.set_model(&next_model);
   soft_constraint.set_model(&next_model);

   softsusy::QedQcd qedqcd;
   qedqcd.to(qedqcd.displayPoleMZ());

   susy_scale_constraint.set_sm_parameters(qedqcd);
   low_scale_constraint.set_sm_parameters(qedqcd);
   soft_constraint.set_sm_parameters(qedqcd);

   soft_constraint.set_boundary_scale(
      [&high_scale_constraint] () {
         return high_scale_constraint.get_scale(); });
   susy_scale_constraint.set_scale(
      [&soft_constraint] () {
         return soft_constraint.get_scale(); });

   high_scale_constraint.initialize();
   susy_scale_constraint.initialize();
   low_scale_constraint.initialize();
   soft_constraint.initialize();

   high_scale_constraint.scale = scales.HighScale;
   susy_scale_constraint.scale = scales.SUSYScale;
   low_scale_constraint.scale = scales.LowScale;
   soft_constraint.scale = scales.SUSYScale;

   CMSSMSemiAnalytic_susy_convergence_tester<Semi_analytic> convergence_tester(
      &next_model, precision);

   Two_scale_increasing_precision running_precision(10.0, precision);

   // apply constraints once
   RGFlow<Two_scale> inner_solver;
   inner_solver.reset();
   inner_solver.set_convergence_tester(&convergence_tester);
   inner_solver.set_running_precision(&running_precision);

   inner_solver.add(&low_scale_constraint, &next_model);
   inner_solver.add(&high_scale_constraint, &next_model);
   inner_solver.add(&susy_scale_constraint, &next_model);

   inner_solver.solve();

   next_model.run_to(scales.SUSYScale);
   soft_constraint.apply();

   return next_model;
}

std::vector<CMSSMSemiAnalytic_input_parameters> initialize_semi_analytic_inputs()
{
   std::vector<CMSSMSemiAnalytic_input_parameters> inputs{3};

   inputs[0].m12 = 660.;
   inputs[0].TanBeta = 40.;
   inputs[0].Azero = 0.;
   inputs[0].MuInput = 550.;

   inputs[1].m12 = 500.;
   inputs[1].TanBeta = 10.;
   inputs[1].Azero = 0.;
   inputs[1].MuInput = 623.4;

   inputs[2].m12 = 920.;
   inputs[2].TanBeta = 40.;
   inputs[2].Azero = 0.;
   inputs[2].MuInput = 700.;

   return inputs;
}

std::vector<CMSSM_input_parameters> initialize_two_scale_inputs()
{
   std::vector<CMSSM_input_parameters> inputs{2};

   inputs[0].m0 = 2800.;
   inputs[0].m12 = 660.;
   inputs[0].TanBeta = 40.;
   inputs[0].SignMu = 1;
   inputs[0].Azero = 0.;

   inputs[1].m0 = 1000.;
   inputs[1].m12 = 660.;
   inputs[1].TanBeta = 40.;
   inputs[1].SignMu = 1;
   inputs[1].Azero = 0.;

   return inputs;
}

BOOST_AUTO_TEST_CASE( test_semi_analytic_to_two_scale )
{
   const double precision = 1.0e-4;

   const std::vector<CMSSMSemiAnalytic_input_parameters> semi_analytic_inputs(
      initialize_semi_analytic_inputs());

   for (const auto& semi_analytic_input: semi_analytic_inputs) {
      softsusy::QedQcd qedqcd;

      Spectrum_generator_settings settings;
      settings.set(Spectrum_generator_settings::precision, precision);

      CMSSMSemiAnalytic_spectrum_generator<Semi_analytic> spectrum_generator;
      spectrum_generator.set_settings(settings);
      spectrum_generator.run(qedqcd, semi_analytic_input);

      const auto& semi_analytic_problems = spectrum_generator.get_problems();

      BOOST_CHECK_EQUAL(semi_analytic_problems.have_problem(), false);

      const double high_scale = spectrum_generator.get_high_scale();
      const double susy_scale = spectrum_generator.get_susy_scale();
      const double low_scale = spectrum_generator.get_low_scale();

      const CMSSMSemiAnalytic<Semi_analytic> semi_analytic_model = std::get<0>(spectrum_generator.get_models());

      CMSSM_scales scales;
      scales.HighScale = high_scale;
      scales.SUSYScale = susy_scale;
      scales.LowScale = low_scale;

      const CMSSM<Two_scale> two_scale_model(
         initialize_two_scale_model(semi_analytic_model, semi_analytic_input));

      const CMSSM<Two_scale> single_iteration_model =
         run_single_two_scale_iteration(two_scale_model, scales);

      // check that the originally found solution is an (approximate)
      // fixed point of the two-scale iteration as well
      const double susy_pars_rel_diff = max_susy_parameters_rel_diff(
         two_scale_model, single_iteration_model);
      const double soft_pars_rel_diff = max_soft_parameters_rel_diff(
         two_scale_model, single_iteration_model);
      const double mass_rel_diff = max_mass_rel_diff(
         two_scale_model, single_iteration_model);

      const double test_precision = 1.0e-2;
      BOOST_CHECK_LT(susy_pars_rel_diff, test_precision);
      BOOST_CHECK_LT(soft_pars_rel_diff, test_precision);
      BOOST_CHECK_LT(mass_rel_diff, test_precision);
   }
}

BOOST_AUTO_TEST_CASE( test_two_scale_to_semi_analytic )
{
   const double precision = 1.0e-4;

   const std::vector<CMSSM_input_parameters> two_scale_inputs(
      initialize_two_scale_inputs());

   for (const auto& two_scale_input: two_scale_inputs) {
      softsusy::QedQcd qedqcd;

      Spectrum_generator_settings settings;
      settings.set(Spectrum_generator_settings::precision, precision);

      CMSSM_spectrum_generator<Two_scale> spectrum_generator;
      spectrum_generator.set_settings(settings);
      spectrum_generator.run(qedqcd, two_scale_input);

      const auto& two_scale_problems = spectrum_generator.get_problems();

      BOOST_CHECK_EQUAL(two_scale_problems.have_problem(), false);

      const double high_scale = spectrum_generator.get_high_scale();
      const double susy_scale = spectrum_generator.get_susy_scale();
      const double low_scale = spectrum_generator.get_low_scale();

      const CMSSM<Two_scale> two_scale_model = std::get<0>(spectrum_generator.get_models());

      CMSSMSemiAnalytic_scales scales;
      scales.HighScale = high_scale;
      scales.SUSYScale = susy_scale;
      scales.LowScale = low_scale;

      CMSSM<Two_scale> high_scale_model(two_scale_model);
      high_scale_model.run_to(high_scale);

      const CMSSMSemiAnalytic<Semi_analytic> semi_analytic_model(
         initialize_semi_analytic_model(two_scale_model, two_scale_input,
                                        high_scale_model.get_Mu(),
                                        high_scale_model.get_BMu(), high_scale));

      const CMSSMSemiAnalytic<Semi_analytic> single_iteration_model =
         run_single_semi_analytic_iteration(semi_analytic_model, scales, precision);

      // check that the originally found solution is an (approximate)
      // fixed point of the semi-analytic iteration as well
      const double susy_pars_rel_diff = max_susy_parameters_rel_diff(
         semi_analytic_model, single_iteration_model);
      const double soft_pars_rel_diff = max_soft_parameters_rel_diff(
         semi_analytic_model, single_iteration_model);
      const double mass_rel_diff = max_mass_rel_diff(
         semi_analytic_model, single_iteration_model);

      const double test_precision = 1.0e-2;
      BOOST_CHECK_LT(susy_pars_rel_diff, test_precision);
      BOOST_CHECK_LT(soft_pars_rel_diff, test_precision);
      BOOST_CHECK_LT(mass_rel_diff, test_precision);
   }
}
