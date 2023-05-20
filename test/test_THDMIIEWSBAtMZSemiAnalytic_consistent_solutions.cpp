#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_THDMIIEWSBAtMZSemiAnalytic_consistent_solutions

#include <boost/test/unit_test.hpp>

#define private public

#include "THDMIIEWSBAtMZSemiAnalytic_input_parameters.hpp"
#include "THDMIIEWSBAtMZSemiAnalytic_slha_io.hpp"
#include "THDMIIEWSBAtMZSemiAnalytic_semi_analytic_ewsb_solver.hpp"
#include "THDMIIEWSBAtMZSemiAnalytic_semi_analytic_low_scale_constraint.hpp"
#include "THDMIIEWSBAtMZSemiAnalytic_semi_analytic_soft_parameters_constraint.hpp"
#include "THDMIIEWSBAtMZSemiAnalytic_semi_analytic_spectrum_generator.hpp"
#include "THDMIIEWSBAtMZSemiAnalytic_semi_analytic_susy_convergence_tester.hpp"
#include "THDMIIEWSBAtMZSemiAnalytic_semi_analytic_susy_scale_constraint.hpp"

#include "THDMII_input_parameters.hpp"
#include "THDMII_slha_io.hpp"
#include "THDMII_two_scale_low_scale_constraint.hpp"
#include "THDMII_two_scale_spectrum_generator.hpp"
#include "THDMII_two_scale_susy_scale_constraint.hpp"

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

   new_model.set_Lambda6(old_model.get_Lambda6());
   new_model.set_Lambda5(old_model.get_Lambda5());
   new_model.set_Lambda7(old_model.get_Lambda7());
   new_model.set_Lambda1(old_model.get_Lambda1());
   new_model.set_Lambda4(old_model.get_Lambda4());
   new_model.set_Lambda3(old_model.get_Lambda3());
   new_model.set_Lambda2(old_model.get_Lambda2());

   new_model.set_M122(old_model.get_M122());
   new_model.set_M112(old_model.get_M112());
   new_model.set_M222(old_model.get_M222());

   new_model.set_v1(old_model.get_v1());
   new_model.set_v2(old_model.get_v2());

   return new_model;
}

template <class NewInputs, class OldInputs>
NewInputs copy_input_parameters(
   const OldInputs& old_input)
{
   NewInputs input;

   input.Qin = old_input.Qin;
   input.Lambda1IN = old_input.Lambda1IN;
   input.Lambda2IN = old_input.Lambda2IN;
   input.Lambda3IN = old_input.Lambda3IN;
   input.Lambda4IN = old_input.Lambda4IN;
   input.Lambda5IN = old_input.Lambda5IN;
   input.Lambda6IN = old_input.Lambda6IN;
   input.Lambda7IN = old_input.Lambda7IN;
   input.M122IN = old_input.M122IN;
   input.TanBeta = old_input.TanBeta;

   return input;
}

THDMII<Two_scale> initialize_two_scale_model(
   const THDMIIEWSBAtMZSemiAnalytic<Semi_analytic>& semi_analytic_model,
   const THDMIIEWSBAtMZSemiAnalytic_input_parameters& semi_analytic_input)
{
   THDMII<Two_scale> two_scale_model
      = copy_parameters_from_model<THDMII<Two_scale> >(semi_analytic_model);

   two_scale_model.set_input_parameters(
      copy_input_parameters<THDMII_input_parameters>(semi_analytic_input));

   two_scale_model.calculate_DRbar_masses();

   return two_scale_model;
}

THDMIIEWSBAtMZSemiAnalytic<Semi_analytic> initialize_semi_analytic_model(
   const THDMII<Two_scale>& two_scale_model,
   const THDMII_input_parameters& two_scale_input,
   double M112, double M222, double susy_scale)
{
   THDMIIEWSBAtMZSemiAnalytic<Semi_analytic> semi_analytic_model
      = copy_parameters_from_model<THDMIIEWSBAtMZSemiAnalytic<Semi_analytic> >(two_scale_model);

   semi_analytic_model.set_input_parameters(
      copy_input_parameters<THDMIIEWSBAtMZSemiAnalytic_input_parameters>(two_scale_input));

   semi_analytic_model.set_M112IN(M112);
   semi_analytic_model.set_M222IN(M222);

   semi_analytic_model.calculate_semi_analytic_solutions(susy_scale);

   THDMIIEWSBAtMZSemiAnalytic_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_semi_analytic_solutions(
      &(semi_analytic_model.get_semi_analytic_solutions()));
   semi_analytic_model.set_ewsb_solver(
      std::make_shared<THDMIIEWSBAtMZSemiAnalytic_ewsb_solver<Semi_analytic> >(ewsb_solver));

   semi_analytic_model.solve_ewsb();
   semi_analytic_model.calculate_DRbar_masses();

   return semi_analytic_model;
}

template <class Model>
double max_susy_parameters_rel_diff(const Model& old_model, const Model& new_model)
{
   std::array<double, 37> diff{};

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
   diff[30] = MaxRelDiff(old_model.get_Lambda6(), new_model.get_Lambda6());
   diff[31] = MaxRelDiff(old_model.get_Lambda5(), new_model.get_Lambda5());
   diff[32] = MaxRelDiff(old_model.get_Lambda7(), new_model.get_Lambda7());
   diff[33] = MaxRelDiff(old_model.get_Lambda1(), new_model.get_Lambda1());
   diff[34] = MaxRelDiff(old_model.get_Lambda4(), new_model.get_Lambda4());
   diff[35] = MaxRelDiff(old_model.get_Lambda3(), new_model.get_Lambda3());
   diff[36] = MaxRelDiff(old_model.get_Lambda2(), new_model.get_Lambda2());

   return *std::max_element(diff.cbegin(), diff.cend());
}

template <class Model>
double max_soft_parameters_rel_diff(const Model& old_model, const Model& new_model)
{
   std::array<double,5> diff{};

   diff[0] = MaxRelDiff(old_model.get_M122(), new_model.get_M122());
   diff[1] = MaxRelDiff(old_model.get_M112(), new_model.get_M112());
   diff[2] = MaxRelDiff(old_model.get_M222(), new_model.get_M222());
   diff[3] = MaxRelDiff(old_model.get_v1(), new_model.get_v1());
   diff[4] = MaxRelDiff(old_model.get_v2(), new_model.get_v2());

   return *std::max_element(diff.cbegin(), diff.cend());
}

template <class Model>
double max_mass_rel_diff(const Model& old_model, const Model& new_model)
{
   std::array<double, 17> diff{};

   diff[0] = MaxRelDiff(old_model.get_MVZ(), new_model.get_MVZ());
   for (int i = 0; i < 2; ++i) {
      diff[i + 1] = MaxRelDiff(old_model.get_Mhh(i), new_model.get_Mhh(i));
   }
   for (int i = 1; i < 2; ++i) {
      diff[i + 3] = MaxRelDiff(old_model.get_MAh(i), new_model.get_MAh(i));
   }
   for (int i = 1; i < 2; ++i) {
      diff[i + 5] = MaxRelDiff(old_model.get_MHm(i), new_model.get_MHm(i));
   }
   for (int i = 0; i < 3; ++i) {
      diff[i + 7] = MaxRelDiff(old_model.get_MFd(i), new_model.get_MFd(i));
   }
   for (int i = 0; i < 3; ++i) {
      diff[i + 10] = MaxRelDiff(old_model.get_MFu(i), new_model.get_MFu(i));
   }
   for (int i = 0; i < 3; ++i) {
      diff[i + 13] = MaxRelDiff(old_model.get_MFe(i), new_model.get_MFe(i));
   }
   diff[16] = MaxRelDiff(old_model.get_MVWm(), new_model.get_MVWm());

   return *std::max_element(diff.cbegin(), diff.cend());
}

THDMII<Two_scale> run_single_two_scale_iteration(
   const THDMII<Two_scale>& model, const THDMII_scales& scales)
{
   THDMII<Two_scale> next_model(model);

   THDMII_susy_scale_constraint<Two_scale> susy_scale_constraint;
   THDMII_low_scale_constraint<Two_scale> low_scale_constraint;

   susy_scale_constraint.set_model(&next_model);
   low_scale_constraint.set_model(&next_model);

   softsusy::QedQcd qedqcd;
   qedqcd.to(qedqcd.displayPoleMZ());

   susy_scale_constraint.set_sm_parameters(qedqcd);
   low_scale_constraint.set_sm_parameters(qedqcd);

   susy_scale_constraint.initialize();
   low_scale_constraint.initialize();

   susy_scale_constraint.scale = scales.SUSYScale;
   low_scale_constraint.scale = scales.LowScale;

   // apply constraints once
   next_model.run_to(scales.LowScale);
   low_scale_constraint.apply();

   next_model.run_to(scales.SUSYScale);
   susy_scale_constraint.apply();

   return next_model;
}

THDMIIEWSBAtMZSemiAnalytic<Semi_analytic> run_single_semi_analytic_iteration(
   const THDMIIEWSBAtMZSemiAnalytic<Semi_analytic>& model,
   const THDMIIEWSBAtMZSemiAnalytic_scales& scales, double precision)
{
   THDMIIEWSBAtMZSemiAnalytic<Semi_analytic> next_model(model);

   THDMIIEWSBAtMZSemiAnalytic_semi_analytic_solutions& solutions(
      next_model.get_semi_analytic_solutions());

   THDMIIEWSBAtMZSemiAnalytic_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_semi_analytic_solutions(&solutions);
   next_model.set_ewsb_solver(
      std::make_shared<THDMIIEWSBAtMZSemiAnalytic_ewsb_solver<Semi_analytic> >(ewsb_solver));

   THDMIIEWSBAtMZSemiAnalytic_susy_scale_constraint<Semi_analytic> susy_scale_constraint;
   THDMIIEWSBAtMZSemiAnalytic_low_scale_constraint<Semi_analytic> low_scale_constraint;
   THDMIIEWSBAtMZSemiAnalytic_soft_parameters_constraint<Semi_analytic> soft_constraint;

   susy_scale_constraint.set_model(&next_model);
   low_scale_constraint.set_model(&next_model);
   soft_constraint.set_model(&next_model);

   softsusy::QedQcd qedqcd;
   qedqcd.to(qedqcd.displayPoleMZ());

   susy_scale_constraint.set_sm_parameters(qedqcd);
   low_scale_constraint.set_sm_parameters(qedqcd);
   soft_constraint.set_sm_parameters(qedqcd);

   soft_constraint.set_boundary_scale(
      [&susy_scale_constraint] () {
         return susy_scale_constraint.get_scale(); });
   low_scale_constraint.set_scale(
      [&soft_constraint] () {
         return soft_constraint.get_scale(); });

   susy_scale_constraint.initialize();
   low_scale_constraint.initialize();
   soft_constraint.initialize();

   susy_scale_constraint.scale = scales.SUSYScale;
   low_scale_constraint.scale = scales.LowScale;
   soft_constraint.scale = scales.LowScale;

   THDMIIEWSBAtMZSemiAnalytic_susy_convergence_tester<Semi_analytic> convergence_tester(
      &next_model, precision);

   Two_scale_increasing_precision running_precision(10., precision);

   // apply constraints once
   RGFlow<Two_scale> inner_solver;
   inner_solver.reset();
   inner_solver.set_convergence_tester(&convergence_tester);
   inner_solver.set_running_precision(&running_precision);

   inner_solver.add(&low_scale_constraint, &next_model);
   inner_solver.add(&susy_scale_constraint, &next_model);

   inner_solver.solve();

   next_model.run_to(scales.LowScale);
   soft_constraint.apply();

   return next_model;
}

template <class Inputs>
std::vector<Inputs> initialize_inputs()
{
   std::vector<Inputs> inputs{1};

   inputs[0].Lambda1IN = 2.;
   inputs[0].Lambda2IN = 0.1;
   inputs[0].Lambda3IN = 0.5;
   inputs[0].Lambda4IN = 0.8;
   inputs[0].Lambda5IN = -2.;
   inputs[0].Lambda6IN = 0.;
   inputs[0].Lambda7IN = 0.;
   inputs[0].M122IN = 1000.;
   inputs[0].TanBeta = 10.;
   inputs[0].Qin = 126.;

   return inputs;
}

BOOST_AUTO_TEST_CASE( test_semi_analytic_to_two_scale )
{
   const double precision = 1.e-4;

   const std::vector<THDMIIEWSBAtMZSemiAnalytic_input_parameters> semi_analytic_inputs(
      initialize_inputs<THDMIIEWSBAtMZSemiAnalytic_input_parameters>());

   for (const auto& semi_analytic_input: semi_analytic_inputs) {
      softsusy::QedQcd qedqcd;

      Spectrum_generator_settings settings;
      settings.set(Spectrum_generator_settings::precision, precision);
      settings.set(Spectrum_generator_settings::calculate_sm_masses, 1);

      THDMIIEWSBAtMZSemiAnalytic_spectrum_generator<Semi_analytic> spectrum_generator;
      spectrum_generator.set_settings(settings);
      spectrum_generator.run(qedqcd, semi_analytic_input);

      const auto& semi_analytic_problems = spectrum_generator.get_problems();

      BOOST_CHECK_EQUAL(semi_analytic_problems.have_problem(), false);

      const double susy_scale = spectrum_generator.get_susy_scale();
      const double low_scale = spectrum_generator.get_low_scale();

      const THDMIIEWSBAtMZSemiAnalytic<Semi_analytic> semi_analytic_model =
         std::get<0>(spectrum_generator.get_models());

      THDMII_scales scales;
      scales.SUSYScale = susy_scale;
      scales.LowScale = low_scale;

      const THDMII<Two_scale> two_scale_model(
         initialize_two_scale_model(semi_analytic_model, semi_analytic_input));

      const THDMII<Two_scale> single_iteration_model =
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

   const std::vector<THDMII_input_parameters> two_scale_inputs(
      initialize_inputs<THDMII_input_parameters>());

   for (const auto& two_scale_input: two_scale_inputs) {
      softsusy::QedQcd qedqcd;

      Spectrum_generator_settings settings;
      settings.set(Spectrum_generator_settings::precision, precision);
      settings.set(Spectrum_generator_settings::calculate_sm_masses, 1);

      THDMII_spectrum_generator<Two_scale> spectrum_generator;
      spectrum_generator.set_settings(settings);
      spectrum_generator.run(qedqcd, two_scale_input);

      const auto& two_scale_problems = spectrum_generator.get_problems();

      BOOST_CHECK_EQUAL(two_scale_problems.have_problem(), false);

      const double susy_scale = spectrum_generator.get_susy_scale();
      const double low_scale = spectrum_generator.get_low_scale();

      const THDMII<Two_scale> two_scale_model =
         std::get<0>(spectrum_generator.get_models());

      THDMIIEWSBAtMZSemiAnalytic_scales scales;
      scales.SUSYScale = susy_scale;
      scales.LowScale = low_scale;

      THDMII<Two_scale> low_scale_model(two_scale_model);
      low_scale_model.run_to(low_scale);

      const THDMIIEWSBAtMZSemiAnalytic<Semi_analytic> semi_analytic_model(
         initialize_semi_analytic_model(low_scale_model, two_scale_input,
                                        two_scale_model.get_M112(),
                                        two_scale_model.get_M222(),
                                        susy_scale));

      const THDMIIEWSBAtMZSemiAnalytic<Semi_analytic> single_iteration_model =
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
