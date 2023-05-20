#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SMSemiAnalytic_consistent_solutions

#include <boost/test/unit_test.hpp>

#define private public

#include "SMSemiAnalytic_input_parameters.hpp"
#include "SMSemiAnalytic_slha_io.hpp"
#include "SMSemiAnalytic_semi_analytic_ewsb_solver.hpp"
#include "SMSemiAnalytic_semi_analytic_high_scale_constraint.hpp"
#include "SMSemiAnalytic_semi_analytic_low_scale_constraint.hpp"
#include "SMSemiAnalytic_semi_analytic_soft_parameters_constraint.hpp"
#include "SMSemiAnalytic_semi_analytic_spectrum_generator.hpp"
#include "SMSemiAnalytic_semi_analytic_susy_convergence_tester.hpp"
#include "SMSemiAnalytic_semi_analytic_susy_scale_constraint.hpp"

#include "SM_input_parameters.hpp"
#include "SM_slha_io.hpp"
#include "SM_two_scale_high_scale_constraint.hpp"
#include "SM_two_scale_low_scale_constraint.hpp"
#include "SM_two_scale_spectrum_generator.hpp"
#include "SM_two_scale_susy_scale_constraint.hpp"

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

   new_model.set_Lambdax(old_model.get_Lambdax());

   new_model.set_v(old_model.get_v());

   new_model.set_mu2(old_model.get_mu2());

   return new_model;
}

template <class NewInputs, class OldInputs>
NewInputs copy_input_parameters(
   const OldInputs& old_input)
{
   NewInputs input;

   input.LambdaIN = old_input.LambdaIN;
   input.Qin = old_input.Qin;
   input.QEWSB = old_input.QEWSB;

   return input;
}

SM<Two_scale> initialize_two_scale_model(
   const SMSemiAnalytic<Semi_analytic>& semi_analytic_model,
   const SMSemiAnalytic_input_parameters& semi_analytic_input)
{
   SM<Two_scale> two_scale_model
      = copy_parameters_from_model<SM<Two_scale> >(semi_analytic_model);

   two_scale_model.set_input_parameters(
      copy_input_parameters<SM_input_parameters>(semi_analytic_input));

   two_scale_model.calculate_DRbar_masses();

   return two_scale_model;
}

SMSemiAnalytic<Semi_analytic> initialize_semi_analytic_model(
   const SM<Two_scale>& two_scale_model,
   const SM_input_parameters& two_scale_input,
   double mu2, double high_scale)
{
   SMSemiAnalytic<Semi_analytic> semi_analytic_model
      = copy_parameters_from_model<SMSemiAnalytic<Semi_analytic> >(two_scale_model);

   semi_analytic_model.set_input_parameters(
      copy_input_parameters<SMSemiAnalytic_input_parameters>(two_scale_input));

   semi_analytic_model.set_mu20(mu2);

   semi_analytic_model.calculate_semi_analytic_solutions(high_scale);

   SMSemiAnalytic_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_semi_analytic_solutions(
      &(semi_analytic_model.get_semi_analytic_solutions()));
   semi_analytic_model.set_ewsb_solver(
      std::make_shared<SMSemiAnalytic_ewsb_solver<Semi_analytic> >(ewsb_solver));

   semi_analytic_model.calculate_DRbar_masses();

   return semi_analytic_model;
}

template <class Model>
double max_susy_parameters_rel_diff(const Model& old_model, const Model& new_model)
{
   std::array<double, 31> diff{};

   diff[0] = MaxRelDiff(old_model.get_g1(), new_model.get_g1());
   diff[1] = MaxRelDiff(old_model.get_g2(), new_model.get_g2());
   diff[2] = MaxRelDiff(old_model.get_g3(), new_model.get_g3());
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 3] = MaxRelDiff(old_model.get_Yd(i,j), new_model.get_Yd(i,j));
      }
   }
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 12] = MaxRelDiff(old_model.get_Ye(i,j), new_model.get_Ye(i,j));
      }
   }
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 21] = MaxRelDiff(old_model.get_Yu(i,j), new_model.get_Yu(i,j));
      }
   }
   diff[30] = MaxRelDiff(old_model.get_Lambdax(), new_model.get_Lambdax());

   return *std::max_element(diff.cbegin(), diff.cend());
}

template <class Model>
double max_soft_parameters_rel_diff(const Model& old_model, const Model& new_model)
{
   std::array<double,2> diff{};

   diff[0] = MaxRelDiff(old_model.get_mu2(), new_model.get_mu2());
   diff[1] = MaxRelDiff(old_model.get_v(), new_model.get_v());

   return *std::max_element(diff.cbegin(), diff.cend());
}

template <class Model>
double max_mass_rel_diff(const Model& old_model, const Model& new_model)
{
   std::array<double, 12> diff{};

   diff[0] = MaxRelDiff(old_model.get_Mhh(), new_model.get_Mhh());
   diff[1] = MaxRelDiff(old_model.get_MVZ(), new_model.get_MVZ());
   for (int i = 0; i < 3; ++i) {
      diff[i + 2] = MaxRelDiff(old_model.get_MFd(i), new_model.get_MFd(i));
   }
   for (int i = 0; i < 3; ++i) {
      diff[i + 5] = MaxRelDiff(old_model.get_MFu(i), new_model.get_MFu(i));
   }
   for (int i = 0; i < 3; ++i) {
      diff[i + 8] = MaxRelDiff(old_model.get_MFe(i), new_model.get_MFe(i));
   }
   diff[11] = MaxRelDiff(old_model.get_MVWp(), new_model.get_MVWp());

   return *std::max_element(diff.cbegin(), diff.cend());
}

SM<Two_scale> run_single_two_scale_iteration(const SM<Two_scale>& model, const SM_scales& scales)
{
   SM<Two_scale> next_model(model);

   SM_high_scale_constraint<Two_scale> high_scale_constraint;
   SM_susy_scale_constraint<Two_scale> susy_scale_constraint;
   SM_low_scale_constraint<Two_scale> low_scale_constraint;

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

SMSemiAnalytic<Semi_analytic> run_single_semi_analytic_iteration(
   const SMSemiAnalytic<Semi_analytic>& model,
   const SMSemiAnalytic_scales& scales, double precision)
{
   SMSemiAnalytic<Semi_analytic> next_model(model);

   SMSemiAnalytic_semi_analytic_solutions& solutions(
      next_model.get_semi_analytic_solutions());

   SMSemiAnalytic_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_semi_analytic_solutions(&solutions);
   next_model.set_ewsb_solver(
      std::make_shared<SMSemiAnalytic_ewsb_solver<Semi_analytic> >(ewsb_solver));

   SMSemiAnalytic_high_scale_constraint<Semi_analytic> high_scale_constraint;
   SMSemiAnalytic_susy_scale_constraint<Semi_analytic> susy_scale_constraint;
   SMSemiAnalytic_low_scale_constraint<Semi_analytic> low_scale_constraint;
   SMSemiAnalytic_soft_parameters_constraint<Semi_analytic> soft_constraint;

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

   SMSemiAnalytic_susy_convergence_tester<Semi_analytic> convergence_tester(
      &next_model, precision);

   Two_scale_increasing_precision running_precision(10., precision);

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

template <class Inputs>
std::vector<Inputs> initialize_inputs()
{
   std::vector<Inputs> inputs{1};

   inputs[0].LambdaIN = 0.192;
   inputs[0].Qin = 1000.;
   inputs[0].QEWSB = 173.3;

   return inputs;
}

BOOST_AUTO_TEST_CASE( test_semi_analytic_to_two_scale )
{
   const double precision = 1.e-4;

   const std::vector<SMSemiAnalytic_input_parameters> semi_analytic_inputs(
      initialize_inputs<SMSemiAnalytic_input_parameters>());

   for (const auto& semi_analytic_input: semi_analytic_inputs) {
      softsusy::QedQcd qedqcd;

      Spectrum_generator_settings settings;
      settings.set(Spectrum_generator_settings::precision, precision);
      settings.set(Spectrum_generator_settings::calculate_sm_masses, 1);

      SMSemiAnalytic_spectrum_generator<Semi_analytic> spectrum_generator;
      spectrum_generator.set_settings(settings);
      spectrum_generator.run(qedqcd, semi_analytic_input);

      const auto& semi_analytic_problems = spectrum_generator.get_problems();

      BOOST_CHECK_EQUAL(semi_analytic_problems.have_problem(), false);

      const double high_scale = spectrum_generator.get_high_scale();
      const double susy_scale = spectrum_generator.get_susy_scale();
      const double low_scale = spectrum_generator.get_low_scale();

      const SMSemiAnalytic<Semi_analytic> semi_analytic_model =
         std::get<0>(spectrum_generator.get_models());

      SM_scales scales;
      scales.HighScale = high_scale;
      scales.SUSYScale = susy_scale;
      scales.LowScale = low_scale;

      const SM<Two_scale> two_scale_model(
         initialize_two_scale_model(semi_analytic_model, semi_analytic_input));

      const SM<Two_scale> single_iteration_model =
         run_single_two_scale_iteration(two_scale_model, scales);

      // check that the originally found solution is an (approximate)
      // fixed point of the two-scale iteration as well
      const double susy_pars_rel_diff = max_susy_parameters_rel_diff(
         two_scale_model, single_iteration_model);
      const double soft_pars_rel_diff = max_soft_parameters_rel_diff(
         two_scale_model, single_iteration_model);
      const double mass_rel_diff = max_mass_rel_diff(
         two_scale_model, single_iteration_model);

      const double test_precision = 1.0e-3;
      BOOST_CHECK_LT(susy_pars_rel_diff, test_precision);
      BOOST_CHECK_LT(soft_pars_rel_diff, test_precision);
      BOOST_CHECK_LT(mass_rel_diff, test_precision);
   }
}

BOOST_AUTO_TEST_CASE( test_two_scale_to_semi_analytic )
{
   const double precision = 1.0e-4;

   const std::vector<SM_input_parameters> two_scale_inputs(
      initialize_inputs<SM_input_parameters>());

   for (const auto& two_scale_input: two_scale_inputs) {
      softsusy::QedQcd qedqcd;

      Spectrum_generator_settings settings;
      settings.set(Spectrum_generator_settings::precision, precision);
      settings.set(Spectrum_generator_settings::calculate_sm_masses, 1);

      SM_spectrum_generator<Two_scale> spectrum_generator;
      spectrum_generator.set_settings(settings);
      spectrum_generator.run(qedqcd, two_scale_input);

      const auto& two_scale_problems = spectrum_generator.get_problems();

      BOOST_CHECK_EQUAL(two_scale_problems.have_problem(), false);

      const double high_scale = spectrum_generator.get_high_scale();
      const double susy_scale = spectrum_generator.get_susy_scale();
      const double low_scale = spectrum_generator.get_low_scale();

      const SM<Two_scale> two_scale_model =
         std::get<0>(spectrum_generator.get_models());

      SMSemiAnalytic_scales scales;
      scales.HighScale = high_scale;
      scales.SUSYScale = susy_scale;
      scales.LowScale = low_scale;

      SM<Two_scale> high_scale_model(two_scale_model);
      high_scale_model.run_to(high_scale);

      const SMSemiAnalytic<Semi_analytic> semi_analytic_model(
         initialize_semi_analytic_model(two_scale_model, two_scale_input,
                                        high_scale_model.get_mu2(),
                                        high_scale));

      const SMSemiAnalytic<Semi_analytic> single_iteration_model =
         run_single_semi_analytic_iteration(semi_analytic_model, scales, precision);

      // check that the originally found solution is an (approximate)
      // fixed point of the semi-analytic iteration as well
      const double susy_pars_rel_diff = max_susy_parameters_rel_diff(
         semi_analytic_model, single_iteration_model);
      const double soft_pars_rel_diff = max_soft_parameters_rel_diff(
         semi_analytic_model, single_iteration_model);
      const double mass_rel_diff = max_mass_rel_diff(
         semi_analytic_model, single_iteration_model);

      const double test_precision = 1.0e-3;
      BOOST_CHECK_LT(susy_pars_rel_diff, test_precision);
      BOOST_CHECK_LT(soft_pars_rel_diff, test_precision);
      BOOST_CHECK_LT(mass_rel_diff, test_precision);
   }
}
