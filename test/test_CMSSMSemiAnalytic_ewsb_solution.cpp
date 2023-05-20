#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSMSemiAnalytic_ewsb_solution

#include <boost/test/unit_test.hpp>

#include "test_CMSSMSemiAnalytic.hpp"
#include "CMSSMSemiAnalytic_semi_analytic_ewsb_solver.hpp"
#include "CMSSMSemiAnalytic_semi_analytic_solutions.hpp"

#include "CMSSM_mass_eigenstates.hpp"

using namespace flexiblesusy;

CMSSM_mass_eigenstates match_at_high_scale(
   const CMSSMSemiAnalytic_mass_eigenstates& model,
   const Boundary_values& values)
{
   CMSSM_mass_eigenstates matched;

   matched.set_scale(model.get_scale());
   matched.set_loops(model.get_loops());

   const auto Yu = model.get_Yu();
   const auto Yd = model.get_Yd();
   const auto Ye = model.get_Ye();

   matched.set_Yu(Yu);
   matched.set_Yd(Yd);
   matched.set_Ye(Ye);

   matched.set_g1(model.get_g1());
   matched.set_g2(model.get_g2());
   matched.set_g3(model.get_g3());

   matched.set_Mu(model.get_Mu());

   matched.set_vd(model.get_vd());
   matched.set_vu(model.get_vu());

   matched.set_TYu(Yu * values.Azero);
   matched.set_TYd(Yd * values.Azero);
   matched.set_TYe(Ye * values.Azero);

   matched.set_BMu(values.BMu0);

   matched.set_mq2(values.m0Sq * UNITMATRIX(3));
   matched.set_mu2(values.m0Sq * UNITMATRIX(3));
   matched.set_md2(values.m0Sq * UNITMATRIX(3));
   matched.set_ml2(values.m0Sq * UNITMATRIX(3));
   matched.set_me2(values.m0Sq * UNITMATRIX(3));
   matched.set_mHd2(values.m0Sq);
   matched.set_mHu2(values.m0Sq);

   matched.set_MassB(values.m12);
   matched.set_MassWB(values.m12);
   matched.set_MassG(values.m12);

   return matched;
}

BOOST_AUTO_TEST_CASE( test_tree_level_ewsb_solutions )
{
   const double precision = 1.0e-5;
   CMSSMSemiAnalytic_input_parameters sa_input;
   CMSSMSemiAnalytic_mass_eigenstates sa_model(sa_input);
   setup_CMSSMSemiAnalytic(sa_model, sa_input);

   const double high_scale = 2.e16;
   sa_model.run_to(high_scale);

   Boundary_values bv;
   setup_high_scale_CMSSMSemiAnalytic(sa_model, bv);

   CMSSMSemiAnalytic_mass_eigenstates high_scale_model(sa_model);

   CMSSMSemiAnalytic_semi_analytic_solutions solns;
   solns.set_input_scale(high_scale);
   solns.set_output_scale(Electroweak_constants::MZ);

   solns.calculate_coefficients(sa_model);

   sa_model.run_to(Electroweak_constants::MZ);

   CMSSMSemiAnalytic_ewsb_solver<Semi_analytic> sa_ewsb;
   sa_ewsb.set_loop_order(0);
   sa_ewsb.set_semi_analytic_solutions(&solns);
   sa_model.set_ewsb_solver(
      std::make_shared<CMSSMSemiAnalytic_ewsb_solver<Semi_analytic> >(sa_ewsb));

   sa_model.set_ewsb_loop_order(0);
   sa_model.set_ewsb_iteration_precision(precision);
   const int sa_error = sa_model.solve_ewsb_tree_level();
   BOOST_CHECK_EQUAL(sa_error, 0);

   BOOST_CHECK_SMALL(sa_model.get_ewsb_eq_hh_1(), precision);
   BOOST_CHECK_SMALL(sa_model.get_ewsb_eq_hh_2(), precision);

   const double m0Sq_sol = sa_model.get_m0Sq();
   const double BMu0_sol = sa_model.get_BMu0();

   bv.m0Sq = m0Sq_sol;
   bv.BMu0 = BMu0_sol;
   CMSSM_mass_eigenstates ts_model = match_at_high_scale(high_scale_model, bv);

   ts_model.run_to(Electroweak_constants::MZ);

   ts_model.set_ewsb_loop_order(0);
   ts_model.set_ewsb_iteration_precision(precision);
   const int ts_error = ts_model.solve_ewsb_tree_level();
   BOOST_CHECK_EQUAL(ts_error, 0);

   BOOST_CHECK_SMALL(ts_model.get_ewsb_eq_hh_1(), precision);
   BOOST_CHECK_SMALL(ts_model.get_ewsb_eq_hh_2(), precision);

   const double Mu_sol = ts_model.get_Mu();
   const double BMu_sol = ts_model.get_BMu();

   BOOST_CHECK_CLOSE(Mu_sol, sa_model.get_Mu(), 5.0e-1);
   BOOST_CHECK_CLOSE(BMu_sol, sa_model.get_BMu(), 5.0e-1);
   BOOST_CHECK_CLOSE(ts_model.get_mHd2(), sa_model.get_mHd2(), 5.0e-1);
   BOOST_CHECK_CLOSE(ts_model.get_mHu2(), sa_model.get_mHu2(), 5.0e-1);
}

BOOST_AUTO_TEST_CASE( test_one_loop_ewsb_solutions )
{
   const double precision = 1.0e-5;
   CMSSMSemiAnalytic_input_parameters sa_input;
   CMSSMSemiAnalytic_mass_eigenstates sa_model(sa_input);
   setup_CMSSMSemiAnalytic(sa_model, sa_input);

   const double high_scale = 2.e16;
   sa_model.run_to(high_scale);

   Boundary_values bv;
   setup_high_scale_CMSSMSemiAnalytic(sa_model, bv);

   CMSSMSemiAnalytic_mass_eigenstates high_scale_model(sa_model);

   CMSSMSemiAnalytic_semi_analytic_solutions solns;
   solns.set_input_scale(high_scale);
   solns.set_output_scale(Electroweak_constants::MZ);

   solns.calculate_coefficients(sa_model);

   sa_model.run_to(Electroweak_constants::MZ);
   sa_model.calculate_DRbar_masses();

   CMSSMSemiAnalytic_ewsb_solver<Semi_analytic> sa_ewsb;
   sa_ewsb.set_loop_order(1);
   sa_ewsb.set_semi_analytic_solutions(&solns);
   sa_model.set_ewsb_solver(
      std::make_shared<CMSSMSemiAnalytic_ewsb_solver<Semi_analytic> >(sa_ewsb));

   sa_model.set_ewsb_loop_order(1);
   sa_model.set_ewsb_iteration_precision(precision);
   const int sa_error = sa_model.solve_ewsb_one_loop();
   BOOST_CHECK_EQUAL(sa_error, 0);

   const std::complex<double> sa_tadpole_hh_1(sa_model.tadpole_hh_1loop(0));
   const std::complex<double> sa_tadpole_hh_2(sa_model.tadpole_hh_1loop(1));

   BOOST_CHECK_SMALL(Im(sa_tadpole_hh_1), 1.0e-12);
   BOOST_CHECK_SMALL(Im(sa_tadpole_hh_2), 1.0e-12);

   BOOST_CHECK_SMALL(sa_model.get_ewsb_eq_hh_1() - Re(sa_tadpole_hh_1), 5.);
   BOOST_CHECK_SMALL(sa_model.get_ewsb_eq_hh_2() - Re(sa_tadpole_hh_2), 5.);

   const double m0Sq_sol = sa_model.get_m0Sq();
   const double BMu0_sol = sa_model.get_BMu0();

   bv.m0Sq = m0Sq_sol;
   bv.BMu0 = BMu0_sol;
   CMSSM_mass_eigenstates ts_model = match_at_high_scale(high_scale_model, bv);

   ts_model.run_to(Electroweak_constants::MZ);
   ts_model.calculate_DRbar_masses();

   ts_model.set_ewsb_loop_order(1);
   ts_model.set_ewsb_iteration_precision(precision);
   const int ts_error = ts_model.solve_ewsb_one_loop();
   BOOST_CHECK_EQUAL(ts_error, 0);

   const std::complex<double> ts_tadpole_hh_1(ts_model.tadpole_hh_1loop(0));
   const std::complex<double> ts_tadpole_hh_2(ts_model.tadpole_hh_1loop(1));

   BOOST_CHECK_SMALL(Im(ts_tadpole_hh_1), 1.0e-12);
   BOOST_CHECK_SMALL(Im(ts_tadpole_hh_2), 1.0e-12);

   BOOST_CHECK_SMALL(ts_model.get_ewsb_eq_hh_1() - Re(ts_tadpole_hh_1), 40.);
   BOOST_CHECK_SMALL(ts_model.get_ewsb_eq_hh_2() - Re(ts_tadpole_hh_2), 20.);

   const double Mu_sol = ts_model.get_Mu();
   const double BMu_sol = ts_model.get_BMu();

   BOOST_CHECK_CLOSE(Mu_sol, sa_model.get_Mu(), 5.0e-1);
   BOOST_CHECK_CLOSE(BMu_sol, sa_model.get_BMu(), 5.0e-1);
   BOOST_CHECK_CLOSE(ts_model.get_mHd2(), sa_model.get_mHd2(), 5.0e-1);
   BOOST_CHECK_CLOSE(ts_model.get_mHu2(), sa_model.get_mHu2(), 5.0e-1);
}

BOOST_AUTO_TEST_CASE( test_two_loop_ewsb_solutions )
{
   const double precision = 1.0e-7;
   CMSSMSemiAnalytic_input_parameters sa_input;
   CMSSMSemiAnalytic_mass_eigenstates sa_model(sa_input);
   setup_CMSSMSemiAnalytic(sa_model, sa_input);

   const double high_scale = 2.e16;
   sa_model.run_to(high_scale);

   Boundary_values bv;
   setup_high_scale_CMSSMSemiAnalytic(sa_model, bv);

   CMSSMSemiAnalytic_mass_eigenstates high_scale_model(sa_model);

   CMSSMSemiAnalytic_semi_analytic_solutions solns;
   solns.set_input_scale(high_scale);
   solns.set_output_scale(Electroweak_constants::MZ);

   solns.calculate_coefficients(sa_model);

   sa_model.run_to(Electroweak_constants::MZ);
   sa_model.calculate_DRbar_masses();

   CMSSMSemiAnalytic_ewsb_solver<Semi_analytic> sa_ewsb;
   sa_ewsb.set_loop_order(2);
   sa_ewsb.set_semi_analytic_solutions(&solns);
   sa_model.set_ewsb_solver(
      std::make_shared<CMSSMSemiAnalytic_ewsb_solver<Semi_analytic> >(sa_ewsb));

   sa_model.set_ewsb_loop_order(2);
   sa_model.set_ewsb_iteration_precision(precision);
   const int sa_error = sa_model.solve_ewsb();
   BOOST_CHECK_EQUAL(sa_error, 0);

   const Eigen::Matrix<double,2,1> sa_tadpole_2l(sa_model.tadpole_hh_2loop());

   const std::complex<double> sa_tadpole_hh_1(
      sa_model.tadpole_hh_1loop(0) + sa_tadpole_2l(0));
   const std::complex<double> sa_tadpole_hh_2(
      sa_model.tadpole_hh_1loop(1) + sa_tadpole_2l(1));

   BOOST_CHECK_SMALL(Im(sa_tadpole_hh_1), 1.0e-12);
   BOOST_CHECK_SMALL(Im(sa_tadpole_hh_2), 1.0e-12);

   BOOST_CHECK_SMALL(sa_model.get_ewsb_eq_hh_1() - Re(sa_tadpole_hh_1), 5.);
   BOOST_CHECK_SMALL(sa_model.get_ewsb_eq_hh_2() - Re(sa_tadpole_hh_2), 5.);

   const double m0Sq_sol = sa_model.get_m0Sq();
   const double BMu0_sol = sa_model.get_BMu0();

   bv.m0Sq = m0Sq_sol;
   bv.BMu0 = BMu0_sol;
   CMSSM_mass_eigenstates ts_model = match_at_high_scale(high_scale_model, bv);

   ts_model.run_to(Electroweak_constants::MZ);
   ts_model.calculate_DRbar_masses();

   ts_model.set_ewsb_loop_order(2);
   ts_model.set_ewsb_iteration_precision(precision);
   const int ts_error = ts_model.solve_ewsb();
   BOOST_CHECK_EQUAL(ts_error, 0);

   const Eigen::Matrix<double,2,1> ts_tadpole_2l(ts_model.tadpole_hh_2loop());

   const std::complex<double> ts_tadpole_hh_1(
      ts_model.tadpole_hh_1loop(0) + ts_tadpole_2l(0));
   const std::complex<double> ts_tadpole_hh_2(
      ts_model.tadpole_hh_1loop(1) + ts_tadpole_2l(1));

   BOOST_CHECK_SMALL(Im(ts_tadpole_hh_1), 1.0e-12);
   BOOST_CHECK_SMALL(Im(ts_tadpole_hh_2), 1.0e-12);

   BOOST_CHECK_SMALL(ts_model.get_ewsb_eq_hh_1() - Re(ts_tadpole_hh_1), 40.);
   BOOST_CHECK_SMALL(ts_model.get_ewsb_eq_hh_2() - Re(ts_tadpole_hh_2), 20.);

   const double Mu_sol = ts_model.get_Mu();
   const double BMu_sol = ts_model.get_BMu();

   BOOST_CHECK_CLOSE(Mu_sol, sa_model.get_Mu(), 5.0e-1);
   BOOST_CHECK_CLOSE(BMu_sol, sa_model.get_BMu(), 5.0e-1);
   BOOST_CHECK_CLOSE(ts_model.get_mHd2(), sa_model.get_mHd2(), 5.0e-1);
   BOOST_CHECK_CLOSE(ts_model.get_mHu2(), sa_model.get_mHu2(), 5.0e-1);
}
