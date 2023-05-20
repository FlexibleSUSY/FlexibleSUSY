#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SSMSemiAnalytic_consistent_solutions

#include <boost/test/unit_test.hpp>

#define private public

#include "SSMSemiAnalytic_input_parameters.hpp"
#include "SSMSemiAnalytic_slha_io.hpp"
#include "SSMSemiAnalytic_semi_analytic_ewsb_solver.hpp"
#include "SSMSemiAnalytic_semi_analytic_high_scale_constraint.hpp"
#include "SSMSemiAnalytic_semi_analytic_low_scale_constraint.hpp"
#include "SSMSemiAnalytic_semi_analytic_soft_parameters_constraint.hpp"
#include "SSMSemiAnalytic_semi_analytic_spectrum_generator.hpp"
#include "SSMSemiAnalytic_semi_analytic_susy_convergence_tester.hpp"
#include "SSMSemiAnalytic_semi_analytic_susy_scale_constraint.hpp"

#include "SSM_input_parameters.hpp"
#include "SSM_slha_io.hpp"
#include "SSM_two_scale_high_scale_constraint.hpp"
#include "SSM_two_scale_low_scale_constraint.hpp"
#include "SSM_two_scale_spectrum_generator.hpp"
#include "SSM_two_scale_susy_scale_constraint.hpp"

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

   new_model.set_LambdaS(old_model.get_LambdaS());
   new_model.set_K2(old_model.get_K2());
   new_model.set_Lambdax(old_model.get_Lambdax());

   new_model.set_Kappa(old_model.get_Kappa());
   new_model.set_K1(old_model.get_K1());
   new_model.set_MS(old_model.get_MS());
   new_model.set_mu2(old_model.get_mu2());

   new_model.set_v(old_model.get_v());
   new_model.set_vS(old_model.get_vS());

   return new_model;
}

template <class NewInputs, class OldInputs>
NewInputs copy_input_parameters(
   const OldInputs& old_input)
{
   NewInputs input;

   input.Qin = old_input.Qin;
   input.QEWSB = old_input.QEWSB;
   input.Lambdainput = old_input.Lambdainput;
   input.LambdaSinput = old_input.LambdaSinput;
   input.Kappainput = old_input.Kappainput;
   input.K1input = old_input.K1input;
   input.K2input = old_input.K2input;
   input.vSInput = old_input.vSInput;

   return input;
}

SSM<Two_scale> initialize_two_scale_model(
   const SSMSemiAnalytic<Semi_analytic>& semi_analytic_model,
   const SSMSemiAnalytic_input_parameters& semi_analytic_input)
{
   SSM<Two_scale> two_scale_model
      = copy_parameters_from_model<SSM<Two_scale> >(semi_analytic_model);

   two_scale_model.set_input_parameters(
      copy_input_parameters<SSM_input_parameters>(semi_analytic_input));

   two_scale_model.calculate_DRbar_masses();

   return two_scale_model;
}

SSMSemiAnalytic<Semi_analytic> initialize_semi_analytic_model(
   const SSM<Two_scale>& two_scale_model,
   const SSM_input_parameters& two_scale_input,
   double MS, double mu2, double high_scale)
{
   SSMSemiAnalytic<Semi_analytic> semi_analytic_model
      = copy_parameters_from_model<SSMSemiAnalytic<Semi_analytic> >(two_scale_model);

   semi_analytic_model.set_input_parameters(
      copy_input_parameters<SSMSemiAnalytic_input_parameters>(two_scale_input));

   semi_analytic_model.set_MS0(MS);
   semi_analytic_model.set_mu20(mu2);

   semi_analytic_model.calculate_semi_analytic_solutions(high_scale);

   SSMSemiAnalytic_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_semi_analytic_solutions(
      &(semi_analytic_model.get_semi_analytic_solutions()));
   semi_analytic_model.set_ewsb_solver(
      std::make_shared<SSMSemiAnalytic_ewsb_solver<Semi_analytic> >(ewsb_solver));

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
   diff[31] = MaxRelDiff(old_model.get_LambdaS(), new_model.get_LambdaS());
   diff[32] = MaxRelDiff(old_model.get_K2(), new_model.get_K2());

   return *std::max_element(diff.cbegin(), diff.cend());
}

template <class Model>
double max_soft_parameters_rel_diff(const Model& old_model, const Model& new_model)
{
   std::array<double,6> diff{};

   diff[0] = MaxRelDiff(old_model.get_mu2(), new_model.get_mu2());
   diff[1] = MaxRelDiff(old_model.get_v(), new_model.get_v());
   diff[2] = MaxRelDiff(old_model.get_vS(), new_model.get_vS());
   diff[3] = MaxRelDiff(old_model.get_Kappa(), new_model.get_Kappa());
   diff[4] = MaxRelDiff(old_model.get_MS(), new_model.get_MS());
   diff[5] = MaxRelDiff(old_model.get_K1(), new_model.get_K1());

   return *std::max_element(diff.cbegin(), diff.cend());
}

template <class Model>
double max_mass_rel_diff(const Model& old_model, const Model& new_model)
{
   std::array<double, 13> diff{};

   diff[0] = MaxRelDiff(old_model.get_MVZ(), new_model.get_MVZ());
   for (int i = 0; i < 2; ++i) {
      diff[i + 1] = MaxRelDiff(old_model.get_Mhh(i), new_model.get_Mhh(i));
   }
   for (int i = 0; i < 3; ++i) {
      diff[i + 3] = MaxRelDiff(old_model.get_MFd(i), new_model.get_MFd(i));
   }
   for (int i = 0; i < 3; ++i) {
      diff[i + 6] = MaxRelDiff(old_model.get_MFu(i), new_model.get_MFu(i));
   }
   for (int i = 0; i < 3; ++i) {
      diff[i + 9] = MaxRelDiff(old_model.get_MFe(i), new_model.get_MFe(i));
   }
   diff[12] = MaxRelDiff(old_model.get_MVWp(), new_model.get_MVWp());

   return *std::max_element(diff.cbegin(), diff.cend());
}

SSM<Two_scale> run_single_two_scale_iteration(const SSM<Two_scale>& model, const SSM_scales& scales)
{
   SSM<Two_scale> next_model(model);

   SSM_high_scale_constraint<Two_scale> high_scale_constraint;
   SSM_susy_scale_constraint<Two_scale> susy_scale_constraint;
   SSM_low_scale_constraint<Two_scale> low_scale_constraint;

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

SSMSemiAnalytic<Semi_analytic> run_single_semi_analytic_iteration(
   const SSMSemiAnalytic<Semi_analytic>& model,
   const SSMSemiAnalytic_scales& scales, double precision)
{
   SSMSemiAnalytic<Semi_analytic> next_model(model);

   SSMSemiAnalytic_semi_analytic_solutions& solutions(
      next_model.get_semi_analytic_solutions());

   SSMSemiAnalytic_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_semi_analytic_solutions(&solutions);
   next_model.set_ewsb_solver(
      std::make_shared<SSMSemiAnalytic_ewsb_solver<Semi_analytic> >(ewsb_solver));

   SSMSemiAnalytic_high_scale_constraint<Semi_analytic> high_scale_constraint;
   SSMSemiAnalytic_susy_scale_constraint<Semi_analytic> susy_scale_constraint;
   SSMSemiAnalytic_low_scale_constraint<Semi_analytic> low_scale_constraint;
   SSMSemiAnalytic_soft_parameters_constraint<Semi_analytic> soft_constraint;

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

   SSMSemiAnalytic_susy_convergence_tester<Semi_analytic> convergence_tester(
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

   inputs[0].Lambdainput = 0.19;
   inputs[0].Qin = 1000.;
   inputs[0].QEWSB = 173.3;
   inputs[0].LambdaSinput = 0.1;
   inputs[0].Kappainput = 100.;
   inputs[0].K1input = -100.;
   inputs[0].K2input = 0.1;
   inputs[0].vSInput = 3.;

   return inputs;
}

BOOST_AUTO_TEST_CASE( test_semi_analytic_to_two_scale )
{
   const double precision = 1.e-4;

   const std::vector<SSMSemiAnalytic_input_parameters> semi_analytic_inputs(
      initialize_inputs<SSMSemiAnalytic_input_parameters>());

   for (const auto& semi_analytic_input: semi_analytic_inputs) {
      softsusy::QedQcd qedqcd;

      Spectrum_generator_settings settings;
      settings.set(Spectrum_generator_settings::precision, precision);
      settings.set(Spectrum_generator_settings::calculate_sm_masses, 1);

      SSMSemiAnalytic_spectrum_generator<Semi_analytic> spectrum_generator;
      spectrum_generator.set_settings(settings);
      spectrum_generator.run(qedqcd, semi_analytic_input);

      const auto& semi_analytic_problems = spectrum_generator.get_problems();

      BOOST_CHECK_EQUAL(semi_analytic_problems.have_problem(), false);

      const double high_scale = spectrum_generator.get_high_scale();
      const double susy_scale = spectrum_generator.get_susy_scale();
      const double low_scale = spectrum_generator.get_low_scale();

      const SSMSemiAnalytic<Semi_analytic> semi_analytic_model =
         std::get<0>(spectrum_generator.get_models());

      SSM_scales scales;
      scales.HighScale = high_scale;
      scales.SUSYScale = susy_scale;
      scales.LowScale = low_scale;

      const SSM<Two_scale> two_scale_model(
         initialize_two_scale_model(semi_analytic_model, semi_analytic_input));

      const SSM<Two_scale> single_iteration_model =
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

   const std::vector<SSM_input_parameters> two_scale_inputs(
      initialize_inputs<SSM_input_parameters>());

   for (const auto& two_scale_input: two_scale_inputs) {
      softsusy::QedQcd qedqcd;

      Spectrum_generator_settings settings;
      settings.set(Spectrum_generator_settings::precision, precision);
      settings.set(Spectrum_generator_settings::calculate_sm_masses, 1);

      SSM_spectrum_generator<Two_scale> spectrum_generator;
      spectrum_generator.set_settings(settings);
      spectrum_generator.run(qedqcd, two_scale_input);

      const auto& two_scale_problems = spectrum_generator.get_problems();

      BOOST_CHECK_EQUAL(two_scale_problems.have_problem(), false);

      const double high_scale = spectrum_generator.get_high_scale();
      const double susy_scale = spectrum_generator.get_susy_scale();
      const double low_scale = spectrum_generator.get_low_scale();

      const SSM<Two_scale> two_scale_model =
         std::get<0>(spectrum_generator.get_models());

      SSMSemiAnalytic_scales scales;
      scales.HighScale = high_scale;
      scales.SUSYScale = susy_scale;
      scales.LowScale = low_scale;

      SSM<Two_scale> high_scale_model(two_scale_model);
      high_scale_model.run_to(high_scale);

      const SSMSemiAnalytic<Semi_analytic> semi_analytic_model(
         initialize_semi_analytic_model(two_scale_model, two_scale_input,
                                        high_scale_model.get_MS(),
                                        high_scale_model.get_mu2(),
                                        high_scale));

      const SSMSemiAnalytic<Semi_analytic> single_iteration_model =
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
