#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_VCMSSM_ewsb

#include <boost/test/unit_test.hpp>

#include "test_VCMSSM.hpp"
#include "CMSSM_two_scale_ewsb_solver.hpp"
#include "VCMSSM_two_scale_ewsb_solver.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_VCMSSM_ewsb_tree_level )
{
   VCMSSM_input_parameters input;
   VCMSSM_mass_eigenstates vcmssm(input);
   const double precision = 1.0e-5;
   setup_VCMSSM(vcmssm, input);

   // initial guess
   const double vev = Sqrt(Sqr(vcmssm.get_vd()) + Sqr(vcmssm.get_vu()));
   vcmssm.set_Mu(input.SignMu * 500.);
   vcmssm.set_vu(vev * Sin(ArcTan(input.TBGuess)));
   vcmssm.set_vd(vev * Cos(ArcTan(input.TBGuess)));

   vcmssm.set_ewsb_iteration_precision(precision);
   const int vcmssm_error = vcmssm.solve_ewsb_tree_level();

   BOOST_CHECK_EQUAL(vcmssm_error, 0);

   BOOST_CHECK_SMALL(vcmssm.get_ewsb_eq_hh_1(), precision);
   BOOST_CHECK_SMALL(vcmssm.get_ewsb_eq_hh_2(), precision);

   const double vcmssm_Mu_soln = vcmssm.get_Mu();
   const double vcmssm_BMu = vcmssm.get_BMu();

   // check that the EWSB solution respects the chosen sign of Mu
   BOOST_CHECK_EQUAL(input.SignMu, Sign(vcmssm_Mu_soln));

   CMSSM_mass_eigenstates cmssm;
   match_CMSSM_to_VCMSSM(cmssm, vcmssm);

   // initial guess: VCMSSM solution with small perturbation
   const double shift = 5.;
   cmssm.set_Mu(vcmssm_Mu_soln + shift);
   cmssm.set_BMu(vcmssm_BMu + shift);

   cmssm.set_ewsb_iteration_precision(precision);
   const int cmssm_error = cmssm.solve_ewsb_tree_level();

   BOOST_CHECK_EQUAL(cmssm_error, 0);

   BOOST_CHECK_SMALL(cmssm.get_ewsb_eq_hh_1(), precision);
   BOOST_CHECK_SMALL(cmssm.get_ewsb_eq_hh_2(), precision);

   const double cmssm_Mu_soln = cmssm.get_Mu();
   const double cmssm_BMu_soln = cmssm.get_BMu();

   BOOST_CHECK_CLOSE_FRACTION(vcmssm_Mu_soln, cmssm_Mu_soln, precision);
   BOOST_CHECK_CLOSE_FRACTION(vcmssm_BMu, cmssm_BMu_soln, precision);
}

BOOST_AUTO_TEST_CASE( test_VCMSSM_ewsb_tree_level_negative_Mu )
{
   VCMSSM_input_parameters input;
   VCMSSM_mass_eigenstates vcmssm(input);
   const double precision = 1.0e-5;
   setup_VCMSSM(vcmssm, input);

   input.SignMu = -1;
   vcmssm.set_input_parameters(input);

   // initial guess
   const double vev = Sqrt(Sqr(vcmssm.get_vd()) + Sqr(vcmssm.get_vu()));
   vcmssm.set_Mu(input.SignMu * 500.);
   vcmssm.set_vu(vev * Sin(ArcTan(input.TBGuess)));
   vcmssm.set_vd(vev * Cos(ArcTan(input.TBGuess)));

   vcmssm.set_ewsb_iteration_precision(precision);
   const int vcmssm_error = vcmssm.solve_ewsb_tree_level();

   BOOST_CHECK_EQUAL(vcmssm_error, 0);

   BOOST_CHECK_SMALL(vcmssm.get_ewsb_eq_hh_1(), precision);
   BOOST_CHECK_SMALL(vcmssm.get_ewsb_eq_hh_2(), precision);

   const double vcmssm_Mu_soln = vcmssm.get_Mu();
   const double vcmssm_BMu = vcmssm.get_BMu();

   // check that the EWSB solution respects the chosen sign of Mu
   BOOST_CHECK_EQUAL(input.SignMu, Sign(vcmssm_Mu_soln));

   CMSSM_mass_eigenstates cmssm;
   match_CMSSM_to_VCMSSM(cmssm, vcmssm);

   // initial guess: VCMSSM solution with small perturbation
   const double shift = 5.;
   cmssm.set_Mu(vcmssm_Mu_soln + shift);
   cmssm.set_BMu(vcmssm_BMu + shift);

   cmssm.set_ewsb_iteration_precision(precision);
   const int cmssm_error = cmssm.solve_ewsb_tree_level();

   BOOST_CHECK_EQUAL(cmssm_error, 0);

   BOOST_CHECK_SMALL(cmssm.get_ewsb_eq_hh_1(), precision);
   BOOST_CHECK_SMALL(cmssm.get_ewsb_eq_hh_2(), precision);

   const double cmssm_Mu_soln = cmssm.get_Mu();
   const double cmssm_BMu_soln = cmssm.get_BMu();

   BOOST_CHECK_CLOSE_FRACTION(vcmssm_Mu_soln, cmssm_Mu_soln, precision);
   BOOST_CHECK_CLOSE_FRACTION(vcmssm_BMu, cmssm_BMu_soln, precision);
}

BOOST_AUTO_TEST_CASE( test_VCMSSM_ewsb_one_loop )
{
   VCMSSM_input_parameters vcmssm_input;
   VCMSSM_mass_eigenstates vcmssm(vcmssm_input);
   const double precision = 1.0e-5;
   setup_VCMSSM(vcmssm, vcmssm_input);

   // initial guess
   const double vev = Sqrt(Sqr(vcmssm.get_vd()) + Sqr(vcmssm.get_vu()));
   vcmssm.set_Mu(vcmssm_input.SignMu * 100.);
   vcmssm.set_vu(vev * Sin(ArcTan(vcmssm_input.TBGuess)));
   vcmssm.set_vd(vev * Cos(ArcTan(vcmssm_input.TBGuess)));

   vcmssm.calculate_DRbar_masses();

   vcmssm.set_ewsb_iteration_precision(precision);
   const int vcmssm_error = vcmssm.solve_ewsb_one_loop();

   BOOST_CHECK_EQUAL(vcmssm_error, 0);

   const std::complex<double> vcmssm_tadpole_hh_1(vcmssm.tadpole_hh_1loop(0));
   const std::complex<double> vcmssm_tadpole_hh_2(vcmssm.tadpole_hh_1loop(1));

   BOOST_CHECK_SMALL(Im(vcmssm_tadpole_hh_1), 1.0e-12);
   BOOST_CHECK_SMALL(Im(vcmssm_tadpole_hh_2), 1.0e-12);

   BOOST_CHECK_SMALL(vcmssm.get_ewsb_eq_hh_1() - Re(vcmssm_tadpole_hh_1), 5.0);
   BOOST_CHECK_SMALL(vcmssm.get_ewsb_eq_hh_2() - Re(vcmssm_tadpole_hh_2), 5.0);

   const double vcmssm_Mu_soln = vcmssm.get_Mu();
   const double vcmssm_BMu = vcmssm.get_BMu();

   // check that the EWSB solution respects the chosen sign of Mu
   BOOST_CHECK_EQUAL(vcmssm_input.SignMu, Sign(vcmssm_Mu_soln));

   CMSSM_input_parameters cmssm_input;
   cmssm_input.m12 = vcmssm_input.m12;
   cmssm_input.m0 = vcmssm_input.m0;
   cmssm_input.Azero = vcmssm_input.Azero;
   cmssm_input.TanBeta = vcmssm.get_vu() / vcmssm.get_vd();
   cmssm_input.SignMu = vcmssm_input.SignMu;
   CMSSM_mass_eigenstates cmssm(cmssm_input);
   match_CMSSM_to_VCMSSM(cmssm, vcmssm);

   cmssm.calculate_DRbar_masses();

   // initial guess: VCMSSM solution with small perturbation
   const double shift = 5.;
   cmssm.set_Mu(vcmssm_Mu_soln + shift);
   cmssm.set_BMu(vcmssm_BMu + shift);

   cmssm.set_ewsb_iteration_precision(precision);
   const int cmssm_error = cmssm.solve_ewsb_one_loop();

   BOOST_CHECK_EQUAL(cmssm_error, 0);

   const std::complex<double> cmssm_tadpole_hh_1(cmssm.tadpole_hh_1loop(0));
   const std::complex<double> cmssm_tadpole_hh_2(cmssm.tadpole_hh_1loop(1));

   BOOST_CHECK_SMALL(Im(cmssm_tadpole_hh_1), 1.0e-12);
   BOOST_CHECK_SMALL(Im(cmssm_tadpole_hh_2), 1.0e-12);

   BOOST_CHECK_SMALL(cmssm.get_ewsb_eq_hh_1() - Re(cmssm_tadpole_hh_1), 2.0);
   BOOST_CHECK_SMALL(cmssm.get_ewsb_eq_hh_2() - Re(cmssm_tadpole_hh_2), 2.0);

   const double cmssm_Mu_soln = cmssm.get_Mu();
   const double cmssm_BMu_soln = cmssm.get_BMu();

   BOOST_CHECK_CLOSE_FRACTION(vcmssm_Mu_soln, cmssm_Mu_soln, precision);
   BOOST_CHECK_CLOSE_FRACTION(vcmssm_BMu, cmssm_BMu_soln, precision);
}

BOOST_AUTO_TEST_CASE( test_VCMSSM_ewsb_one_loop_negative_Mu )
{
   VCMSSM_input_parameters vcmssm_input;
   VCMSSM_mass_eigenstates vcmssm(vcmssm_input);
   const double precision = 1.0e-5;
   setup_VCMSSM(vcmssm, vcmssm_input);

   vcmssm_input.SignMu = -1;
   vcmssm.set_input_parameters(vcmssm_input);

   // initial guess
   const double vev = Sqrt(Sqr(vcmssm.get_vd()) + Sqr(vcmssm.get_vu()));
   vcmssm.set_Mu(vcmssm_input.SignMu * 100.);
   vcmssm.set_vu(vev * Sin(ArcTan(vcmssm_input.TBGuess)));
   vcmssm.set_vd(vev * Cos(ArcTan(vcmssm_input.TBGuess)));

   vcmssm.calculate_DRbar_masses();

   vcmssm.set_ewsb_iteration_precision(precision);
   const int vcmssm_error = vcmssm.solve_ewsb_one_loop();

   BOOST_CHECK_EQUAL(vcmssm_error, 0);

   const std::complex<double> vcmssm_tadpole_hh_1(vcmssm.tadpole_hh_1loop(0));
   const std::complex<double> vcmssm_tadpole_hh_2(vcmssm.tadpole_hh_1loop(1));

   BOOST_CHECK_SMALL(Im(vcmssm_tadpole_hh_1), 1.0e-12);
   BOOST_CHECK_SMALL(Im(vcmssm_tadpole_hh_2), 1.0e-12);

   BOOST_CHECK_CLOSE_FRACTION(vcmssm.get_ewsb_eq_hh_1(), Re(vcmssm_tadpole_hh_1), 5*precision);
   BOOST_CHECK_CLOSE_FRACTION(vcmssm.get_ewsb_eq_hh_2(), Re(vcmssm_tadpole_hh_2), 5*precision);

   const double vcmssm_Mu_soln = vcmssm.get_Mu();
   const double vcmssm_BMu = vcmssm.get_BMu();

   // check that the EWSB solution respects the chosen sign of Mu
   BOOST_CHECK_EQUAL(vcmssm_input.SignMu, Sign(vcmssm_Mu_soln));

   CMSSM_input_parameters cmssm_input;
   cmssm_input.m12 = vcmssm_input.m12;
   cmssm_input.m0 = vcmssm_input.m0;
   cmssm_input.Azero = vcmssm_input.Azero;
   cmssm_input.TanBeta = vcmssm.get_vu() / vcmssm.get_vd();
   cmssm_input.SignMu = vcmssm_input.SignMu;
   CMSSM_mass_eigenstates cmssm(cmssm_input);
   match_CMSSM_to_VCMSSM(cmssm, vcmssm);

   cmssm.calculate_DRbar_masses();

   // initial guess: VCMSSM solution with small perturbation
   const double shift = 5.;
   cmssm.set_Mu(vcmssm_Mu_soln + shift);
   cmssm.set_BMu(vcmssm_BMu + shift);

   cmssm.set_ewsb_iteration_precision(precision);
   const int cmssm_error = cmssm.solve_ewsb_one_loop();

   BOOST_CHECK_EQUAL(cmssm_error, 0);

   const std::complex<double> cmssm_tadpole_hh_1(cmssm.tadpole_hh_1loop(0));
   const std::complex<double> cmssm_tadpole_hh_2(cmssm.tadpole_hh_1loop(1));

   BOOST_CHECK_SMALL(Im(cmssm_tadpole_hh_1), 1.0e-12);
   BOOST_CHECK_SMALL(Im(cmssm_tadpole_hh_2), 1.0e-12);

   BOOST_CHECK_SMALL(cmssm.get_ewsb_eq_hh_1() - Re(cmssm_tadpole_hh_1), 2.0);
   BOOST_CHECK_SMALL(cmssm.get_ewsb_eq_hh_2() - Re(cmssm_tadpole_hh_2), 2.0);

   const double cmssm_Mu_soln = cmssm.get_Mu();
   const double cmssm_BMu_soln = cmssm.get_BMu();

   BOOST_CHECK_CLOSE_FRACTION(vcmssm_Mu_soln, cmssm_Mu_soln, precision);
   BOOST_CHECK_CLOSE_FRACTION(vcmssm_BMu, cmssm_BMu_soln, precision);
}
