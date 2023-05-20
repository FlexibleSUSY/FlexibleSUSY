#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CE6SSM_consistent_solutions

#include <boost/test/unit_test.hpp>

#define private public

#include "CE6SSM_input_parameters.hpp"
#include "CE6SSM_slha_io.hpp"
#include "CE6SSM_semi_analytic_ewsb_solver.hpp"
#include "CE6SSM_semi_analytic_high_scale_constraint.hpp"
#include "CE6SSM_semi_analytic_low_scale_constraint.hpp"
#include "CE6SSM_semi_analytic_soft_parameters_constraint.hpp"
#include "CE6SSM_semi_analytic_spectrum_generator.hpp"
#include "CE6SSM_semi_analytic_susy_convergence_tester.hpp"
#include "CE6SSM_semi_analytic_susy_scale_constraint.hpp"

#include "E6SSM_input_parameters.hpp"
#include "E6SSM_slha_io.hpp"
#include "E6SSM_two_scale_high_scale_constraint.hpp"
#include "E6SSM_two_scale_low_scale_constraint.hpp"
#include "E6SSM_two_scale_spectrum_generator.hpp"
#include "E6SSM_two_scale_susy_scale_constraint.hpp"

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

   new_model.set_Kappa(old_model.get_Kappa());
   new_model.set_Lambda12(old_model.get_Lambda12());
   new_model.set_Lambdax(old_model.get_Lambdax());

   new_model.set_g1(old_model.get_g1());
   new_model.set_g2(old_model.get_g2());
   new_model.set_g3(old_model.get_g3());
   new_model.set_gN(old_model.get_gN());

   new_model.set_MuPr(old_model.get_MuPr());

   new_model.set_vd(old_model.get_vd());
   new_model.set_vu(old_model.get_vu());
   new_model.set_vs(old_model.get_vs());

   new_model.set_TYu(old_model.get_TYu());
   new_model.set_TYd(old_model.get_TYd());
   new_model.set_TYe(old_model.get_TYe());

   new_model.set_TKappa(old_model.get_TKappa());
   new_model.set_TLambda12(old_model.get_TLambda12());
   new_model.set_TLambdax(old_model.get_TLambdax());

   new_model.set_mHd2(old_model.get_mHd2());
   new_model.set_mHu2(old_model.get_mHu2());
   new_model.set_ms2(old_model.get_ms2());
   new_model.set_mHp2(old_model.get_mHp2());
   new_model.set_mHpbar2(old_model.get_mHpbar2());
   new_model.set_BMuPr(old_model.get_BMuPr());

   new_model.set_mq2(old_model.get_mq2());
   new_model.set_mu2(old_model.get_mu2());
   new_model.set_md2(old_model.get_md2());
   new_model.set_ml2(old_model.get_ml2());
   new_model.set_me2(old_model.get_me2());
   new_model.set_mDx2(old_model.get_mDx2());
   new_model.set_mDxbar2(old_model.get_mDxbar2());
   new_model.set_mH1I2(old_model.get_mH1I2());
   new_model.set_mH2I2(old_model.get_mH2I2());
   new_model.set_msI2(old_model.get_msI2());

   new_model.set_MassB(old_model.get_MassB());
   new_model.set_MassBp(old_model.get_MassBp());
   new_model.set_MassWB(old_model.get_MassWB());
   new_model.set_MassG(old_model.get_MassG());

   return new_model;
}

E6SSM_input_parameters initialize_two_scale_input(
   double Azero, double m12, double m0Sq,
   double MuPr, double BMuPr,
   const CE6SSM_input_parameters& semi_analytic_input)
{
   E6SSM_input_parameters input;

   input.Lambda12Input = semi_analytic_input.Lambda12Input;
   input.vSInput = semi_analytic_input.vsInput;
   input.BmuPrimeInput = BMuPr;
   input.muPrimeInput = MuPr;
   input.KappaInput = semi_analytic_input.KappaInput;
   input.LambdaInput = semi_analytic_input.LambdaInput;
   input.Azero = Azero;
   input.TanBeta = semi_analytic_input.TanBeta;
   input.m12 = m12;
   input.m0 = Sqrt(m0Sq);

   return input;
}

E6SSM<Two_scale> initialize_two_scale_model(
   const CE6SSM<Semi_analytic>& semi_analytic_model,
   const CE6SSM_input_parameters& semi_analytic_input)
{
   E6SSM<Two_scale> two_scale_model
      = copy_parameters_from_model<E6SSM<Two_scale> >(semi_analytic_model);

   two_scale_model.set_input_parameters(
      initialize_two_scale_input(semi_analytic_model.get_Azero(),
                                 semi_analytic_model.get_m12(),
                                 semi_analytic_model.get_m0Sq(),
                                 semi_analytic_model.get_MuPr(),
                                 semi_analytic_model.get_BMuPr(),
                                 semi_analytic_input));

   two_scale_model.calculate_DRbar_masses();

   return two_scale_model;
}

CE6SSM_input_parameters initialize_semi_analytic_input(
   double MuPr, double BMuPr, const E6SSM_input_parameters& two_scale_input)
{
   CE6SSM_input_parameters input;

   input.TanBeta = two_scale_input.TanBeta;
   input.LambdaInput = two_scale_input.LambdaInput;
   input.KappaInput = two_scale_input.KappaInput;
   input.MuPrimeInput = MuPr;
   input.BMuPrimeInput = BMuPr;
   input.vsInput = two_scale_input.vSInput;
   input.Lambda12Input = two_scale_input.Lambda12Input;
   input.m0SqGuess = Sqr(two_scale_input.vSInput);
   input.m12Guess = two_scale_input.vSInput / two_scale_input.TanBeta;
   input.AzeroGuess = two_scale_input.vSInput / two_scale_input.TanBeta;

   return input;
}

CE6SSM<Semi_analytic> initialize_semi_analytic_model(
   const E6SSM<Two_scale>& two_scale_model,
   const E6SSM_input_parameters& two_scale_input,
   double MuPr, double BMuPr, double high_scale)
{
   CE6SSM<Semi_analytic> semi_analytic_model
      = copy_parameters_from_model<CE6SSM<Semi_analytic> >(two_scale_model);

   semi_analytic_model.set_input_parameters(
      initialize_semi_analytic_input(MuPr, BMuPr, two_scale_input));

   semi_analytic_model.set_Azero(two_scale_input.Azero);
   semi_analytic_model.set_m12(two_scale_input.m12);
   semi_analytic_model.set_m0Sq(Sqr(two_scale_input.m0));
   semi_analytic_model.set_MuPrBV(MuPr);

   semi_analytic_model.calculate_semi_analytic_solutions(high_scale);
   semi_analytic_model.calculate_DRbar_masses();

   return semi_analytic_model;
}

template <class Model>
double max_susy_parameters_rel_diff(const Model& old_model, const Model& new_model)
{
   std::array<double, 49> diff{};

   diff[0] = MaxRelDiff(old_model.get_g1(), new_model.get_g1());
   diff[1] = MaxRelDiff(old_model.get_g2(), new_model.get_g2());
   diff[2] = MaxRelDiff(old_model.get_g3(), new_model.get_g3());
   diff[3] = MaxRelDiff(old_model.get_gN(), new_model.get_gN());
   diff[4] = MaxRelDiff(old_model.get_vd(), new_model.get_vd());
   diff[5] = MaxRelDiff(old_model.get_vs(), new_model.get_vs());
   diff[6] = MaxRelDiff(old_model.get_vu(), new_model.get_vu());
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 7] = MaxRelDiff(old_model.get_Yd(i,j), new_model.get_Yd(i,j));
      }
   }
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 16] = MaxRelDiff(old_model.get_Ye(i,j), new_model.get_Ye(i,j));
      }
   }
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 25] = MaxRelDiff(old_model.get_Yu(i,j), new_model.get_Yu(i,j));
      }
   }
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 34] = MaxRelDiff(old_model.get_Kappa(i,j), new_model.get_Kappa(i,j));
      }
   }
   diff[43] = MaxRelDiff(old_model.get_Lambdax(), new_model.get_Lambdax());
   for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
         diff[i + 2*j + 44] = MaxRelDiff(old_model.get_Lambda12(i,j), new_model.get_Lambda12(i,j));
      }
   }
   diff[48] = MaxRelDiff(old_model.get_MuPr(), new_model.get_MuPr());

   return *std::max_element(diff.cbegin(), diff.cend());
}

template <class Model>
double max_soft_parameters_rel_diff(const Model& old_model, const Model& new_model)
{
   std::array<double,126> diff{};

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
         diff[i + 3*j + 18] = MaxRelDiff(old_model.get_TKappa(i,j), new_model.get_TKappa(i,j));
      }
   }
   for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
         diff[i + 2*j + 27] = MaxRelDiff(old_model.get_TLambda12(i,j), new_model.get_TLambda12(i,j));
      }
   }
   diff[31] = MaxRelDiff(old_model.get_TLambdax(),new_model.get_TLambdax());
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 32] = MaxRelDiff(old_model.get_TYu(i,j), new_model.get_TYu(i,j));
      }
   }
   diff[41] = MaxRelDiff(old_model.get_BMuPr(), new_model.get_BMuPr());
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 42] = MaxRelDiff(old_model.get_mq2(i,j), new_model.get_mq2(i,j));
      }
   }
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 51] = MaxRelDiff(old_model.get_ml2(i,j), new_model.get_ml2(i,j));
      }
   }
   diff[60] = MaxRelDiff(old_model.get_mHd2(), new_model.get_mHd2());
   diff[61] = MaxRelDiff(old_model.get_mHu2(), new_model.get_mHu2());
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 62] = MaxRelDiff(old_model.get_md2(i,j), new_model.get_md2(i,j));
      }
   }
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 71] = MaxRelDiff(old_model.get_mu2(i,j), new_model.get_mu2(i,j));
      }
   }
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 80] = MaxRelDiff(old_model.get_me2(i,j), new_model.get_me2(i,j));
      }
   }
   diff[89] = MaxRelDiff(old_model.get_ms2(),new_model.get_ms2());
   for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
         diff[i + 2*j + 90] = MaxRelDiff(old_model.get_mH1I2(i,j), new_model.get_mH1I2(i,j));
      }
   }
   for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
         diff[i + 2*j + 94] = MaxRelDiff(old_model.get_mH2I2(i,j), new_model.get_mH2I2(i,j));
      }
   }
   for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
         diff[i + 2*j + 98] = MaxRelDiff(old_model.get_msI2(i,j), new_model.get_msI2(i,j));
      }
   }
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 102] = MaxRelDiff(old_model.get_mDx2(i,j), new_model.get_mDx2(i,j));
      }
   }
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 111] = MaxRelDiff(old_model.get_mDxbar2(i,j), new_model.get_mDxbar2(i,j));
      }
   }
   diff[120] = MaxRelDiff(old_model.get_mHp2(), new_model.get_mHp2());
   diff[121] = MaxRelDiff(old_model.get_mHpbar2(), new_model.get_mHpbar2());
   diff[122] = MaxRelDiff(old_model.get_MassB(), new_model.get_MassB());
   diff[123] = MaxRelDiff(old_model.get_MassWB(), new_model.get_MassWB());
   diff[124] = MaxRelDiff(old_model.get_MassG(), new_model.get_MassG());
   diff[125] = MaxRelDiff(old_model.get_MassBp(), new_model.get_MassBp());

   return *std::max_element(diff.cbegin(), diff.cend());
}

template <class Model>
double max_mass_rel_diff(const Model& old_model, const Model& new_model)
{
   std::array<double, 73> diff{};

   diff[0] = MaxRelDiff(old_model.get_MGlu(), new_model.get_MGlu());
   diff[1] = MaxRelDiff(old_model.get_MChaP(), new_model.get_MChaP());
   diff[2] = MaxRelDiff(old_model.get_MVZp(), new_model.get_MVZp());
   for (int i = 0; i < 6; ++i) {
      diff[i + 3] = MaxRelDiff(old_model.get_MSd(i), new_model.get_MSd(i));
   }
   for (int i = 0; i < 3; ++i) {
      diff[i + 9] = MaxRelDiff(old_model.get_MSv(i), new_model.get_MSv(i));
   }
   for (int i = 0; i < 6; ++i) {
      diff[i + 12] = MaxRelDiff(old_model.get_MSu(i), new_model.get_MSu(i));
   }
   for (int i = 0; i < 6; ++i) {
      diff[i + 18] = MaxRelDiff(old_model.get_MSe(i), new_model.get_MSe(i));
   }
   for (int i = 0; i < 6; ++i) {
      diff[i + 24] = MaxRelDiff(old_model.get_MSDX(i), new_model.get_MSDX(i));
   }
   for (int i = 0; i < 3; ++i) {
      diff[i + 30] = MaxRelDiff(old_model.get_Mhh(i), new_model.get_Mhh(i));
   }
   for (int i = 2; i < 3; ++i) {
      diff[i + 33] = MaxRelDiff(old_model.get_MAh(i), new_model.get_MAh(i));
   }
   for (int i = 1; i < 2; ++i) {
      diff[i + 36] = MaxRelDiff(old_model.get_MHpm(i), new_model.get_MHpm(i));
   }
   for (int i = 0; i < 6; ++i) {
      diff[i + 38] = MaxRelDiff(old_model.get_MChi(i), new_model.get_MChi(i));
   }
   for (int i = 0; i < 2; ++i) {
      diff[i + 44] = MaxRelDiff(old_model.get_MCha(i), new_model.get_MCha(i));
   }
   for (int i = 0; i < 3; ++i) {
      diff[i + 46] = MaxRelDiff(old_model.get_MFDX(i), new_model.get_MFDX(i));
   }
   for (int i = 0; i < 4; ++i) {
      diff[i + 49] = MaxRelDiff(old_model.get_MSHI0(i), new_model.get_MSHI0(i));
   }
   for (int i = 0; i < 4; ++i) {
      diff[i + 53] = MaxRelDiff(old_model.get_MSHIp(i), new_model.get_MSHIp(i));
   }
   for (int i = 0; i < 2; ++i) {
      diff[i + 57] = MaxRelDiff(old_model.get_MChaI(i), new_model.get_MChaI(i));
   }
   for (int i = 0; i < 4; ++i) {
      diff[i + 59] = MaxRelDiff(old_model.get_MChiI(i), new_model.get_MChiI(i));
   }
   for (int i = 0; i < 2; ++i) {
      diff[i + 63] = MaxRelDiff(old_model.get_MSSI0(i), new_model.get_MSSI0(i));
   }
   for (int i = 0; i < 2; ++i) {
      diff[i + 65] = MaxRelDiff(old_model.get_MFSI(i), new_model.get_MFSI(i));
   }
   for (int i = 0; i < 2; ++i) {
      diff[i + 67] = MaxRelDiff(old_model.get_MSHp0(i), new_model.get_MSHp0(i));
   }
   for (int i = 0; i < 2; ++i) {
      diff[i + 69] = MaxRelDiff(old_model.get_MSHpp(i), new_model.get_MSHpp(i));
   }
   for (int i = 0; i < 2; ++i) {
      diff[i + 71] = MaxRelDiff(old_model.get_MChiP(i), new_model.get_MChiP(i));
   }

   return *std::max_element(diff.cbegin(), diff.cend());
}

E6SSM<Two_scale> run_single_two_scale_iteration(const E6SSM<Two_scale>& model, const E6SSM_scales& scales)
{
   E6SSM<Two_scale> next_model(model);

   E6SSM_high_scale_constraint<Two_scale> high_scale_constraint;
   E6SSM_susy_scale_constraint<Two_scale> susy_scale_constraint;
   E6SSM_low_scale_constraint<Two_scale> low_scale_constraint;

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

CE6SSM<Semi_analytic> run_single_semi_analytic_iteration(
   const CE6SSM<Semi_analytic>& model,
   const CE6SSM_scales& scales,
   double precision)
{
   CE6SSM<Semi_analytic> next_model(model);

   CE6SSM_semi_analytic_solutions& solutions(
      next_model.get_semi_analytic_solutions());

   CE6SSM_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_semi_analytic_solutions(&solutions);
   next_model.set_ewsb_solver(
      std::make_shared<CE6SSM_ewsb_solver<Semi_analytic> >(ewsb_solver));

   CE6SSM_high_scale_constraint<Semi_analytic> high_scale_constraint;
   CE6SSM_susy_scale_constraint<Semi_analytic> susy_scale_constraint;
   CE6SSM_low_scale_constraint<Semi_analytic> low_scale_constraint;
   CE6SSM_soft_parameters_constraint<Semi_analytic> soft_constraint;

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

   CE6SSM_susy_convergence_tester<Semi_analytic> convergence_tester(
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

std::vector<CE6SSM_input_parameters> initialize_semi_analytic_inputs()
{
   std::vector<CE6SSM_input_parameters> inputs{1};

   inputs[0].TanBeta = 10.;
   inputs[0].LambdaInput = 0.2;
   inputs[0].KappaInput = 0.15;
   inputs[0].MuPrimeInput = 10000.;
   inputs[0].BMuPrimeInput = 10000.;
   inputs[0].vsInput = 6000.;
   inputs[0].Lambda12Input = 0.2;
   inputs[0].m0SqGuess = Sqr(6000.);
   inputs[0].m12Guess = 600.;
   inputs[0].AzeroGuess = 600.;

   return inputs;
}

std::vector<E6SSM_input_parameters> initialize_two_scale_inputs()
{
   std::vector<E6SSM_input_parameters> inputs{1};

   inputs[0].Lambda12Input = 0.2;
   inputs[0].vSInput = 6000.;
   inputs[0].BmuPrimeInput = -1.42595046e7;
   inputs[0].muPrimeInput = 1.56425576e4;
   inputs[0].KappaInput = 0.15;
   inputs[0].LambdaInput = 0.2;
   inputs[0].Azero = 3.60559408e3;
   inputs[0].TanBeta = 10.;
   inputs[0].m12 = 1.92911828e3;
   inputs[0].m0 = Sqrt(1.85970254e7);

   return inputs;
}

BOOST_AUTO_TEST_CASE( test_semi_analytic_to_two_scale )
{
   const double precision = 1.0e-4;

   const std::vector<CE6SSM_input_parameters> semi_analytic_inputs(
      initialize_semi_analytic_inputs());

   for (const auto& semi_analytic_input: semi_analytic_inputs) {
      softsusy::QedQcd qedqcd;

      Spectrum_generator_settings settings;
      settings.set(Spectrum_generator_settings::precision, precision);

      CE6SSM_spectrum_generator<Semi_analytic> spectrum_generator;
      spectrum_generator.set_settings(settings);
      spectrum_generator.run(qedqcd, semi_analytic_input);

      const auto& semi_analytic_problems = spectrum_generator.get_problems();

      BOOST_CHECK_EQUAL(semi_analytic_problems.have_problem(), false);

      const double high_scale = spectrum_generator.get_high_scale();
      const double susy_scale = spectrum_generator.get_susy_scale();
      const double low_scale = spectrum_generator.get_low_scale();

      const CE6SSM<Semi_analytic> semi_analytic_model = std::get<0>(spectrum_generator.get_models());

      E6SSM_scales scales;
      scales.HighScale = high_scale;
      scales.SUSYScale = susy_scale;
      scales.LowScale = low_scale;

      const E6SSM<Two_scale> two_scale_model(
         initialize_two_scale_model(semi_analytic_model, semi_analytic_input));

      const E6SSM<Two_scale> single_iteration_model =
         run_single_two_scale_iteration(two_scale_model, scales);

      // check that the originally found solution is an (approximate)
      // fixed point of the two-scale iteration as well
      const double susy_pars_rel_diff = max_susy_parameters_rel_diff(
         two_scale_model, single_iteration_model);
      const double soft_pars_rel_diff = max_soft_parameters_rel_diff(
         two_scale_model, single_iteration_model);
      const double mass_rel_diff = max_mass_rel_diff(
         two_scale_model, single_iteration_model);

      const double test_precision = 2.0e-2;
      BOOST_CHECK_LT(susy_pars_rel_diff, test_precision);
      BOOST_CHECK_LT(soft_pars_rel_diff, test_precision);
      BOOST_CHECK_LT(mass_rel_diff, test_precision);
   }
}

BOOST_AUTO_TEST_CASE( test_two_scale_to_semi_analytic )
{
   const double precision = 1.0e-4;

   const std::vector<E6SSM_input_parameters> two_scale_inputs(
      initialize_two_scale_inputs());

   for (const auto& two_scale_input: two_scale_inputs) {
      softsusy::QedQcd qedqcd;

      Spectrum_generator_settings settings;
      settings.set(Spectrum_generator_settings::precision, precision);

      E6SSM_spectrum_generator<Two_scale> spectrum_generator;
      spectrum_generator.set_settings(settings);
      spectrum_generator.run(qedqcd, two_scale_input);

      const auto& two_scale_problems = spectrum_generator.get_problems();

      BOOST_CHECK_EQUAL(two_scale_problems.have_problem(), false);

      const double high_scale = spectrum_generator.get_high_scale();
      const double susy_scale = spectrum_generator.get_susy_scale();
      const double low_scale = spectrum_generator.get_low_scale();

      const E6SSM<Two_scale> two_scale_model = std::get<0>(spectrum_generator.get_models());

      CE6SSM_scales scales;
      scales.HighScale = high_scale;
      scales.SUSYScale = susy_scale;
      scales.LowScale = low_scale;

      E6SSM<Two_scale> high_scale_model(two_scale_model);
      high_scale_model.run_to(high_scale);

      const CE6SSM<Semi_analytic> semi_analytic_model(
         initialize_semi_analytic_model(two_scale_model, two_scale_input,
                                        high_scale_model.get_MuPr(),
                                        high_scale_model.get_BMuPr(),
                                        high_scale));

      const CE6SSM<Semi_analytic> single_iteration_model =
         run_single_semi_analytic_iteration(semi_analytic_model, scales, precision);

      // check that the originally found solution is an (approximate)
      // fixed point of the semi-analytic iteration as well
      const double susy_pars_rel_diff = max_susy_parameters_rel_diff(
         semi_analytic_model, single_iteration_model);
      const double soft_pars_rel_diff = max_soft_parameters_rel_diff(
         semi_analytic_model, single_iteration_model);
      const double mass_rel_diff = max_mass_rel_diff(
         semi_analytic_model, single_iteration_model);

      // note: the attainable precision is limited by how precisely
      // the inputs can be chosen so that the universality constraints end
      // up being satisfied at the high-scale
      const double test_precision = 0.5;
      BOOST_CHECK_LT(susy_pars_rel_diff, test_precision);
      BOOST_CHECK_LT(soft_pars_rel_diff, test_precision);
      BOOST_CHECK_LT(mass_rel_diff, test_precision);
   }
}
