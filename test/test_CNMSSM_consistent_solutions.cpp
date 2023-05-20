#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CNMSSM_consistent_solutions

#include <boost/test/unit_test.hpp>

#define private public

#include "CNMSSM_input_parameters.hpp"
#include "CNMSSM_slha_io.hpp"
#include "CNMSSM_semi_analytic_ewsb_solver.hpp"
#include "CNMSSM_semi_analytic_high_scale_constraint.hpp"
#include "CNMSSM_semi_analytic_low_scale_constraint.hpp"
#include "CNMSSM_semi_analytic_soft_parameters_constraint.hpp"
#include "CNMSSM_semi_analytic_spectrum_generator.hpp"
#include "CNMSSM_semi_analytic_susy_convergence_tester.hpp"
#include "CNMSSM_semi_analytic_susy_scale_constraint.hpp"

#include "NMSSM_input_parameters.hpp"
#include "NMSSM_slha_io.hpp"
#include "NMSSM_two_scale_high_scale_constraint.hpp"
#include "NMSSM_two_scale_low_scale_constraint.hpp"
#include "NMSSM_two_scale_spectrum_generator.hpp"
#include "NMSSM_two_scale_susy_scale_constraint.hpp"

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

   new_model.set_Lambdax(old_model.get_Lambdax());
   new_model.set_Kappa(old_model.get_Kappa());

   new_model.set_g1(old_model.get_g1());
   new_model.set_g2(old_model.get_g2());
   new_model.set_g3(old_model.get_g3());

   new_model.set_vd(old_model.get_vd());
   new_model.set_vu(old_model.get_vu());
   new_model.set_vS(old_model.get_vS());

   new_model.set_TYu(old_model.get_TYu());
   new_model.set_TYd(old_model.get_TYd());
   new_model.set_TYe(old_model.get_TYe());

   new_model.set_TLambdax(old_model.get_TLambdax());
   new_model.set_TKappa(old_model.get_TKappa());

   new_model.set_mHd2(old_model.get_mHd2());
   new_model.set_mHu2(old_model.get_mHu2());
   new_model.set_ms2(old_model.get_ms2());

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

NMSSM_input_parameters initialize_two_scale_input(
   double m0Sq, const CNMSSM_input_parameters& semi_analytic_input)
{
   NMSSM_input_parameters input;

   input.LambdaInput = semi_analytic_input.LambdaInput;
   input.Azero = semi_analytic_input.Azero;
   input.SignvS = semi_analytic_input.SignvS;
   input.TanBeta = semi_analytic_input.TanBeta;
   input.m12 = semi_analytic_input.m12;
   input.m0 = Sqrt(m0Sq);

   return input;
}

NMSSM<Two_scale> initialize_two_scale_model(
   const CNMSSM<Semi_analytic>& semi_analytic_model,
   const CNMSSM_input_parameters& semi_analytic_input)
{
   NMSSM<Two_scale> two_scale_model
      = copy_parameters_from_model<NMSSM<Two_scale> >(semi_analytic_model);

   two_scale_model.set_input_parameters(
      initialize_two_scale_input(semi_analytic_model.get_m0Sq(),
                                 semi_analytic_input));

   two_scale_model.calculate_DRbar_masses();

   return two_scale_model;
}

CNMSSM_input_parameters initialize_semi_analytic_input(
   const NMSSM_input_parameters& two_scale_input)
{
   CNMSSM_input_parameters input;

   input.Azero = two_scale_input.Azero;
   input.SignvS = two_scale_input.SignvS;
   input.TanBeta = two_scale_input.TanBeta;
   input.m12 = two_scale_input.m12;
   input.LambdaInput = two_scale_input.LambdaInput;

   return input;
}

CNMSSM<Semi_analytic> initialize_semi_analytic_model(
   const NMSSM<Two_scale>& two_scale_model,
   const NMSSM_input_parameters& two_scale_input,
   double high_scale)
{
   CNMSSM<Semi_analytic> semi_analytic_model
      = copy_parameters_from_model<CNMSSM<Semi_analytic> >(two_scale_model);

   semi_analytic_model.set_input_parameters(
      initialize_semi_analytic_input(two_scale_input));

   semi_analytic_model.set_m0Sq(Sqr(two_scale_input.m0));

   semi_analytic_model.calculate_semi_analytic_solutions(high_scale);
   semi_analytic_model.calculate_DRbar_masses();

   return semi_analytic_model;
}

template <class Model>
double max_susy_parameters_rel_diff(const Model& old_model, const Model& new_model)
{
   std::array<double, 35> diff{};

   diff[0] = MaxRelDiff(old_model.get_g1(),new_model.get_g1());
   diff[1] = MaxRelDiff(old_model.get_g2(),new_model.get_g2());
   diff[2] = MaxRelDiff(old_model.get_g3(),new_model.get_g3());
   diff[3] = MaxRelDiff(old_model.get_vd(),new_model.get_vd());
   diff[4] = MaxRelDiff(old_model.get_vS(),new_model.get_vS());
   diff[5] = MaxRelDiff(old_model.get_vu(),new_model.get_vu());
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 6] = MaxRelDiff(old_model.get_Yd(i,j),new_model.get_Yd(i,j));
      }
   }
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 15] = MaxRelDiff(old_model.get_Ye(i,j),new_model.get_Ye(i,j));
      }
   }
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 24] = MaxRelDiff(old_model.get_Yu(i,j),new_model.get_Yu(i,j));
      }
   }
   diff[33] = MaxRelDiff(old_model.get_Kappa(),new_model.get_Kappa());
   diff[34] = MaxRelDiff(old_model.get_Lambdax(),new_model.get_Lambdax());

   return *std::max_element(diff.cbegin(), diff.cend());
}

template <class Model>
double max_soft_parameters_rel_diff(const Model& old_model, const Model& new_model)
{
   std::array<double,80> diff{};

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
   diff[27] = MaxRelDiff(old_model.get_TLambdax(), new_model.get_TLambdax());
   diff[28] = MaxRelDiff(old_model.get_TKappa(), new_model.get_TKappa());
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 29] = MaxRelDiff(old_model.get_mq2(i,j), new_model.get_mq2(i,j));
      }
   }
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 38] = MaxRelDiff(old_model.get_ml2(i,j), new_model.get_ml2(i,j));
      }
   }
   diff[47] = MaxRelDiff(old_model.get_mHd2(), new_model.get_mHd2());
   diff[48] = MaxRelDiff(old_model.get_mHu2(), new_model.get_mHu2());
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 49] = MaxRelDiff(old_model.get_md2(i,j), new_model.get_md2(i,j));
      }
   }
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 58] = MaxRelDiff(old_model.get_mu2(i,j), new_model.get_mu2(i,j));
      }
   }
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 67] = MaxRelDiff(old_model.get_me2(i,j), new_model.get_me2(i,j));
      }
   }
   diff[76] = MaxRelDiff(old_model.get_ms2(), new_model.get_ms2());
   diff[77] = MaxRelDiff(old_model.get_MassB(), new_model.get_MassB());
   diff[78] = MaxRelDiff(old_model.get_MassWB(), new_model.get_MassWB());
   diff[79] = MaxRelDiff(old_model.get_MassG(), new_model.get_MassG());

   return *std::max_element(diff.cbegin(), diff.cend());
}

template <class Model>
double max_mass_rel_diff(const Model& old_model, const Model& new_model)
{
   std::array<double, 37> diff{};

   diff[0] = MaxRelDiff(old_model.get_MGlu(),new_model.get_MGlu());
   for (int i = 0; i < 6; ++i) {
      diff[i + 1] = MaxRelDiff(old_model.get_MSd(i),new_model.get_MSd(i));
   }
   for (int i = 0; i < 3; ++i) {
      diff[i + 7] = MaxRelDiff(old_model.get_MSv(i),new_model.get_MSv(i));
   }
   for (int i = 0; i < 6; ++i) {
      diff[i + 10] = MaxRelDiff(old_model.get_MSu(i),new_model.get_MSu(i));
   }
   for (int i = 0; i < 6; ++i) {
      diff[i + 16] = MaxRelDiff(old_model.get_MSe(i),new_model.get_MSe(i));
   }
   for (int i = 0; i < 3; ++i) {
      diff[i + 22] = MaxRelDiff(old_model.get_Mhh(i),new_model.get_Mhh(i));
   }
   for (int i = 1; i < 3; ++i) {
      diff[i + 25] = MaxRelDiff(old_model.get_MAh(i),new_model.get_MAh(i));
   }
   for (int i = 1; i < 2; ++i) {
      diff[i + 28] = MaxRelDiff(old_model.get_MHpm(i),new_model.get_MHpm(i));
   }
   for (int i = 0; i < 5; ++i) {
      diff[i + 30] = MaxRelDiff(old_model.get_MChi(i),new_model.get_MChi(i));
   }
   for (int i = 0; i < 2; ++i) {
      diff[i + 35] = MaxRelDiff(old_model.get_MCha(i),new_model.get_MCha(i));
   }

   return *std::max_element(diff.cbegin(), diff.cend());
}

NMSSM<Two_scale> run_single_two_scale_iteration(const NMSSM<Two_scale>& model, const NMSSM_scales& scales)
{
   NMSSM<Two_scale> next_model(model);

   NMSSM_high_scale_constraint<Two_scale> high_scale_constraint;
   NMSSM_susy_scale_constraint<Two_scale> susy_scale_constraint;
   NMSSM_low_scale_constraint<Two_scale> low_scale_constraint;

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

CNMSSM<Semi_analytic> run_single_semi_analytic_iteration(
   const CNMSSM<Semi_analytic>& model,
   const CNMSSM_scales& scales,
   double precision)
{
   CNMSSM<Semi_analytic> next_model(model);

   CNMSSM_semi_analytic_solutions& solutions(
      next_model.get_semi_analytic_solutions());

   CNMSSM_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_semi_analytic_solutions(&solutions);
   next_model.set_ewsb_solver(
      std::make_shared<CNMSSM_ewsb_solver<Semi_analytic> >(ewsb_solver));

   CNMSSM_high_scale_constraint<Semi_analytic> high_scale_constraint;
   CNMSSM_susy_scale_constraint<Semi_analytic> susy_scale_constraint;
   CNMSSM_low_scale_constraint<Semi_analytic> low_scale_constraint;
   CNMSSM_soft_parameters_constraint<Semi_analytic> soft_constraint;

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

   CNMSSM_susy_convergence_tester<Semi_analytic> convergence_tester(
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

std::vector<CNMSSM_input_parameters> initialize_semi_analytic_inputs()
{
   std::vector<CNMSSM_input_parameters> inputs{1};

   inputs[0].Azero = -300.;
   inputs[0].SignvS = 1;
   inputs[0].TanBeta = 10.;
   inputs[0].m12 = 133.333;
   inputs[0].LambdaInput = -0.05;

   return inputs;
}

std::vector<NMSSM_input_parameters> initialize_two_scale_inputs()
{
   std::vector<NMSSM_input_parameters> inputs{1};

   inputs[0].LambdaInput = -0.05;
   inputs[0].Azero = -300.;
   inputs[0].SignvS = 1;
   inputs[0].TanBeta = 10.;
   inputs[0].m12 = 133.333;
   inputs[0].m0 = 36.31420167217364;

   return inputs;
}

BOOST_AUTO_TEST_CASE( test_semi_analytic_to_two_scale )
{
   const double precision = 1.0e-4;

   const std::vector<CNMSSM_input_parameters> semi_analytic_inputs(
      initialize_semi_analytic_inputs());

   for (const auto& semi_analytic_input: semi_analytic_inputs) {
      softsusy::QedQcd qedqcd;

      Spectrum_generator_settings settings;
      settings.set(Spectrum_generator_settings::precision, precision);

      CNMSSM_spectrum_generator<Semi_analytic> spectrum_generator;
      spectrum_generator.set_settings(settings);
      spectrum_generator.run(qedqcd, semi_analytic_input);

      const auto& semi_analytic_problems = spectrum_generator.get_problems();

      BOOST_CHECK_EQUAL(semi_analytic_problems.have_problem(), false);

      const double high_scale = spectrum_generator.get_high_scale();
      const double susy_scale = spectrum_generator.get_susy_scale();
      const double low_scale = spectrum_generator.get_low_scale();

      const CNMSSM<Semi_analytic> semi_analytic_model
         = std::get<0>(spectrum_generator.get_models());

      NMSSM_scales scales;
      scales.HighScale = high_scale;
      scales.SUSYScale = susy_scale;
      scales.LowScale = low_scale;

      const NMSSM<Two_scale> two_scale_model(
         initialize_two_scale_model(semi_analytic_model, semi_analytic_input));

      const NMSSM<Two_scale> single_iteration_model =
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
      BOOST_CHECK_LT(soft_pars_rel_diff, 5.0 * test_precision);
      BOOST_CHECK_LT(mass_rel_diff, test_precision);
   }
}

BOOST_AUTO_TEST_CASE( test_two_scale_to_semi_analytic )
{
   const double precision = 1.0e-6;

   const std::vector<NMSSM_input_parameters> two_scale_inputs(
      initialize_two_scale_inputs());

   for (const auto& two_scale_input: two_scale_inputs) {
      softsusy::QedQcd qedqcd;

      Spectrum_generator_settings settings;
      settings.set(Spectrum_generator_settings::precision, precision);

      NMSSM_spectrum_generator<Two_scale> spectrum_generator;
      spectrum_generator.set_settings(settings);
      spectrum_generator.run(qedqcd, two_scale_input);

      const auto& two_scale_problems = spectrum_generator.get_problems();

      BOOST_CHECK_EQUAL(two_scale_problems.have_problem(), false);

      const double high_scale = spectrum_generator.get_high_scale();
      const double susy_scale = spectrum_generator.get_susy_scale();
      const double low_scale = spectrum_generator.get_low_scale();

      const NMSSM<Two_scale> two_scale_model
         = std::get<0>(spectrum_generator.get_models());

      CNMSSM_scales scales;
      scales.HighScale = high_scale;
      scales.SUSYScale = susy_scale;
      scales.LowScale = low_scale;

      NMSSM<Two_scale> high_scale_model(two_scale_model);
      high_scale_model.run_to(high_scale);

      const CNMSSM<Semi_analytic> semi_analytic_model(
         initialize_semi_analytic_model(two_scale_model, two_scale_input,
                                        high_scale));

      const CNMSSM<Semi_analytic> single_iteration_model =
         run_single_semi_analytic_iteration(semi_analytic_model, scales, precision);

      // check that the originally found solution is an (approximate)
      // fixed point of the semi-analytic iteration as well
      const double susy_pars_rel_diff = max_susy_parameters_rel_diff(
         semi_analytic_model, single_iteration_model);
      const double soft_pars_rel_diff = max_soft_parameters_rel_diff(
         semi_analytic_model, single_iteration_model);
      const double mass_rel_diff = max_mass_rel_diff(
         semi_analytic_model, single_iteration_model);

      const double test_precision = 2.0e-2;
      BOOST_CHECK_LT(susy_pars_rel_diff, test_precision);
      BOOST_CHECK_LT(soft_pars_rel_diff, test_precision);
      BOOST_CHECK_LT(mass_rel_diff, test_precision);
   }
}
